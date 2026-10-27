"""Separate whole-field and contour fits of empirical Cai coordinate maps.

All coordinates are in exit diameters.  This module keeps the physical-point
mask, baseline density bands, and volume estimators fixed while coefficients
change.  There is no density amplitude parameter and no Jacobian multiplier.
The defaults describe the available cases; a whole-Kn holdout is opt-in.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import math
import time
from typing import Callable, Optional, Sequence

import numpy as np
from scipy.optimize import least_squares

from .scaling import CaiPredictor, CoordinateMap


CONTOUR_LEVELS = (0.1, 0.01, 0.001)
DEFAULT_ANGLE_EDGES = np.deg2rad(np.arange(0.0, 92.5, 2.5))
DENSITY_BAND_EDGES = (0.001, 0.01, 0.1)


@dataclass
class ScalingCase:
    """Raw matched samples; ``occupied`` means a meaningful DSMC estimator.

    Density and optional scalar translational temperature are normalized by
    exit values. Volumes must be real cell volumes, not an inferred mesh count.
    Missing occupancy defaults conservatively to positive density. A resolved
    zero can be retained by explicitly supplying ``occupied=True``. Optional
    ``metadata['domain_valid']`` excludes known solid/outside-domain cells.
    ``density_blocks`` in metadata contains independent time-block estimates;
    their independence must be established by the caller, not by occupancy.
    """

    name: str
    kn: float
    points: np.ndarray
    density: np.ndarray
    volumes: np.ndarray
    speed_ratio: float = 2.0
    temperature: Optional[np.ndarray] = None
    occupied: Optional[np.ndarray] = None
    suitable: bool = False
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        self.points = np.asarray(self.points, dtype=float).reshape(-1, 2)
        self.density = np.asarray(self.density, dtype=float).reshape(-1)
        self.volumes = np.asarray(self.volumes, dtype=float).reshape(-1)
        n = self.density.size
        if self.points.shape[0] != n or self.volumes.size != n:
            raise ValueError("points, density and actual cell volumes must match")
        if not float(self.kn) > 0 or not np.isfinite(self.kn):
            raise ValueError("a case needs a finite positive exit Kn")
        if self.occupied is None:
            self.occupied = self.density > 0
        self.occupied = np.asarray(self.occupied, dtype=bool).reshape(-1)
        if self.occupied.size != n:
            raise ValueError("occupied must match density")
        if self.temperature is not None:
            self.temperature = np.asarray(self.temperature, dtype=float).reshape(-1)
            if self.temperature.size != n:
                raise ValueError("temperature must match density")


@dataclass
class RayEstimator:
    """Fixed radial-shell estimator within one physical angular bin."""

    angle_bin: int
    angle_deg: float
    indices: list
    radii: np.ndarray
    gaps_after: np.ndarray

    def average(self, values, volumes):
        return np.asarray([np.average(values[i], weights=volumes[i])
                           for i in self.indices])


@dataclass
class PreparedCase:
    case: ScalingCase
    mask: np.ndarray
    points: np.ndarray
    density: np.ndarray
    volumes: np.ndarray
    baseline: np.ndarray
    weights: np.ndarray
    angle_bins: np.ndarray
    density_bands: np.ndarray
    coverage: dict
    rays: list
    angle_edges: np.ndarray
    floor: float


def _predictor(mapping, case, options=None, factory=None):
    return (factory or CaiPredictor)(mapping, speed_ratio=case.speed_ratio,
                                    **(options or {}))


def _fixed_weights(volumes, angles, bands):
    """Equal occupied angle-band strata, volume weights inside each stratum."""
    weights = np.zeros(volumes.size)
    strata = sorted(set(zip(angles.tolist(), bands.tolist())))
    for a, b in strata:
        m = (angles == a) & (bands == b)
        weights[m] = volumes[m] / volumes[m].sum() / len(strata)
    return weights, strata


def _prepare_rays(points, volumes, angles, angle_edges, radial_bins):
    """Volume averages in fixed angle and radius bins, with explicit gaps.

    Logarithmic shell edges cover the full observed positive-radius range.
    Empty shells interrupt a profile: crossings are never interpolated across
    an unsampled radial gap.  This is a matched cell estimator, not a claim
    that a cell centre is an exact sample on the bin's central ray.
    """
    radius = np.hypot(points[:, 0], points[:, 1])
    if not radius.size:
        return []
    lo, hi = radius.min(), radius.max()
    if hi <= lo:
        return []
    edges = np.geomspace(lo, np.nextafter(hi, np.inf), int(radial_bins) + 1)
    shells = np.clip(np.searchsorted(edges, radius, side="right") - 1,
                     0, radial_bins - 1)
    rays = []
    for a in np.unique(angles):
        selected = angles == a
        shell_ids = np.unique(shells[selected])
        indices = [np.flatnonzero(selected & (shells == b)) for b in shell_ids]
        radii = np.asarray([np.average(radius[i], weights=volumes[i]) for i in indices])
        theta = np.average(np.arctan2(points[selected, 1], points[selected, 0]),
                           weights=volumes[selected])
        rays.append(RayEstimator(int(a), float(np.rad2deg(theta)), indices,
                                 radii, np.diff(shell_ids) > 1))
    return rays


def prepare_case(case: ScalingCase, *, angle_edges=None, floor=1e-10,
                 radial_bins=80, evaluator_options=None,
                 evaluator_factory=None) -> PreparedCase:
    """Prepare fixed full-domain masks and fixed baseline-defined strata.

    There is no lip exclusion and no analytical-density cutoff.  The numerical
    log floor does not change sample selection. Exclusion counts are disjoint
    and add to the raw count, including undefined analytical exit-plane points.
    """
    if floor <= 0 or not np.isfinite(floor):
        raise ValueError("floor must be finite and positive")
    angle_edges = np.asarray(DEFAULT_ANGLE_EDGES if angle_edges is None else angle_edges,
                             dtype=float)
    if (angle_edges.ndim != 1 or angle_edges.size < 2 or
            np.any(np.diff(angle_edges) <= 0) or angle_edges[0] > 0 or
            angle_edges[-1] < math.pi / 2):
        raise ValueError("angle_edges must increase and cover [0, pi/2]")
    if int(radial_bins) < 2:
        raise ValueError("radial_bins must be at least two")
    n = case.density.size
    mask = np.ones(n, dtype=bool)
    exclusions = {}

    def exclude(label, bad):
        bad = np.asarray(bad, dtype=bool)
        exclusions[label] = int(np.count_nonzero(mask & bad))
        mask[bad] = False

    exclude("invalid_file_values", (~np.all(np.isfinite(case.points), axis=1) |
            ~np.isfinite(case.density) | (case.density < 0) |
            ~np.isfinite(case.volumes) | (case.volumes <= 0)))
    domain = np.asarray(case.metadata.get("domain_valid", np.ones(n)), dtype=bool)
    if domain.size != n:
        raise ValueError("domain_valid must match raw samples")
    exclude("solid_or_outside_domain", ~domain | (case.points[:, 0] < 0) |
            (case.points[:, 1] < 0))
    exclude("undefined_analytical_boundary", case.points[:, 0] == 0)
    exclude("no_meaningful_dsmc_estimator", ~case.occupied)
    baseline = np.full(n, np.nan)
    if np.any(mask):
        baseline[mask] = _predictor(CoordinateMap(), case, evaluator_options,
                                    evaluator_factory).density(case.points[mask], case.kn)
    exclude("undefined_baseline", ~np.isfinite(baseline) | (baseline < 0))
    points, density, volumes = case.points[mask], case.density[mask], case.volumes[mask]
    base = baseline[mask]
    theta = np.arctan2(points[:, 1], points[:, 0])
    angles = np.clip(np.searchsorted(angle_edges, theta, side="right") - 1,
                     0, angle_edges.size - 2)
    bands = np.searchsorted(DENSITY_BAND_EDGES, base, side="right")
    weights, strata = _fixed_weights(volumes, angles, bands)
    coverage = {
        "raw_samples": n, "valid_samples": int(mask.sum()), "excluded": exclusions,
        "coverage_fraction": float(mask.mean()) if n else 0.0,
        "occupied_strata": len(strata),
        "possible_strata": int((angle_edges.size - 1) * 4),
        "empty_strata": int((angle_edges.size - 1) * 4 - len(strata)),
        "empty_strata_policy": "omitted; remaining occupied strata equally weighted",
        "strata": [{"angle_bin": a, "density_band": b,
                    "samples": int(np.count_nonzero((angles == a) & (bands == b))),
                    "weight": float(weights[(angles == a) & (bands == b)].sum())}
                   for a, b in strata],
        "dsmc_floor_active": int(np.count_nonzero(density < floor)),
        "baseline_floor_active": int(np.count_nonzero(base < floor)),
        "resolved_zero_samples": int(np.count_nonzero(density == 0)),
        "lip_excluded": 0,
        "near_exit_samples": int(np.count_nonzero(points[:, 0] <= 1.0)),
        "near_lip_samples": int(np.count_nonzero(np.hypot(points[:, 0], points[:, 1] - .5) <= .1)),
        "floor": float(floor), "band_reference": "fixed direct Cai baseline",
    }
    return PreparedCase(case, mask, points, density, volumes, base, weights,
                        angles, bands, coverage,
                        _prepare_rays(points, volumes, angles, angle_edges, int(radial_bins)),
                        angle_edges, float(floor))


def extract_crossing(radii, density, level, *, gaps_after=None):
    """A unique decreasing crossing, interpolated in log radius/log density.

    No extrapolation, forced crossing, or automatic choice between branches.
    Multiple crossings are ``multiple``; a plateau exactly at the level is
    ``ambiguous``. A contour at an outer endpoint is domain-truncated rather
    than considered known. Missing shells invalidate a crossing across a gap.
    """
    r, d = np.asarray(radii, float), np.asarray(density, float)
    if r.size != d.size or r.ndim != 1:
        raise ValueError("one-dimensional radii and density must match")
    if level <= 0:
        raise ValueError("positive contour level required")
    if (r.size < 2 or np.any(~np.isfinite(r)) or np.any(r <= 0) or
            np.any(np.diff(r) <= 0) or np.any(~np.isfinite(d)) or np.any(d < 0)):
        return {"status": "missing", "radius": None, "crossings": 0}
    gaps = np.zeros(r.size - 1, bool) if gaps_after is None else np.asarray(gaps_after, bool)
    if gaps.size != r.size - 1:
        raise ValueError("gaps_after must have len(radii)-1 entries")
    exact = d == level
    if np.any(exact[:-1] & exact[1:]):
        return {"status": "ambiguous", "radius": None, "crossings": None}
    positions = []
    interrupted = False
    # Exact interior vertices count once. Endpoints are never accepted.
    for i in np.flatnonzero(exact[1:-1]) + 1:
        if gaps[i - 1] or gaps[i]:
            interrupted = True
        else:
            positions.append((float(r[i]), bool(d[i - 1] > d[i + 1])))
    for i in range(r.size - 1):
        if (d[i] - level) * (d[i + 1] - level) >= 0:
            continue
        if gaps[i]:
            interrupted = True
            continue
        if d[i] <= 0 or d[i + 1] <= 0:
            interrupted = True
            continue
        fraction = (math.log(level) - math.log(d[i])) / (math.log(d[i + 1]) - math.log(d[i]))
        positions.append((float(np.exp(np.log(r[i]) + fraction * np.log(r[i + 1] / r[i]))),
                          bool(d[i] > d[i + 1])))
    if len(positions) > 1:
        return {"status": "multiple", "radius": None, "crossings": len(positions)}
    if interrupted:
        return {"status": "missing_gap", "radius": None, "crossings": len(positions)}
    if len(positions) == 1:
        if not positions[0][1]:
            return {"status": "ambiguous_increasing", "radius": None, "crossings": 1}
        return {"status": "ok", "radius": positions[0][0], "crossings": 1}
    if d[-1] >= level:
        status = "domain_truncated_outer"
    elif d[0] <= level:
        status = "domain_truncated_inner"
    else:
        status = "missing"
    return {"status": status, "radius": None, "crossings": 0}


def contour_diagnostics(prepared: PreparedCase, predicted, *, levels=CONTOUR_LEVELS):
    """All requested levels/angular bins, including every missing branch flag."""
    predicted = np.asarray(predicted, float)
    if predicted.shape != prepared.density.shape:
        raise ValueError("predictions must use the fixed prepared mask")
    rows = []
    by_bin = {ray.angle_bin: ray for ray in prepared.rays}
    for a in range(prepared.angle_edges.size - 1):
        ray = by_bin.get(a)
        for level in levels:
            row = {"case": prepared.case.name, "kn": float(prepared.case.kn),
                   "angle_bin": a, "level": float(level),
                   "angle_deg": (float(np.rad2deg((prepared.angle_edges[a] + prepared.angle_edges[a + 1]) / 2))
                                 if ray is None else ray.angle_deg)}
            crossings = {}
            for label, field_values in (("dsmc", prepared.density),
                                        ("baseline", prepared.baseline), ("model", predicted)):
                cross = ({"status": "missing", "radius": None, "crossings": 0} if ray is None else
                         extract_crossing(ray.radii, ray.average(field_values, prepared.volumes),
                                          level, gaps_after=ray.gaps_after))
                row[f"{label}_radius"] = cross["radius"]
                row[f"{label}_status"] = cross["status"]
                crossings[label] = cross
            rd, rb, rm = row["dsmc_radius"], row["baseline_radius"], row["model_radius"]
            row["status"] = ("ok" if rd is not None and rm is not None else
                             f"dsmc:{row['dsmc_status']};model:{row['model_status']}")
            row["ratio_dsmc_to_baseline"] = rd / rb if rd is not None and rb is not None else None
            row["ratio_dsmc_to_model"] = rd / rm if rd is not None and rm is not None else None
            row["signed_displacement_D"] = rm - rd if rd is not None and rm is not None else None
            row["baseline_signed_displacement_D"] = rb - rd if rd is not None and rb is not None else None
            row["log_radius_residual"] = math.log(rm / rd) if rd is not None and rm is not None else None
            rows.append(row)
    return rows


def _linear_metrics(reference, prediction, volumes):
    if not len(reference):
        return {"count": 0, "normalized_L2": None, "bias": None, "mean_absolute": None,
                "rmse": None, "max_absolute": None}
    if np.any(~np.isfinite(prediction)):
        return {"count": len(reference), "invalid_model_samples": int(np.count_nonzero(~np.isfinite(prediction))),
                "normalized_L2": None, "bias": None, "mean_absolute": None,
                "rmse": None, "max_absolute": None}
    error = prediction - reference
    denom = float(np.dot(volumes, reference ** 2))
    return {"count": len(reference),
            "normalized_L2": math.sqrt(float(np.dot(volumes, error ** 2)) / denom) if denom > 0 else None,
            "bias": float(np.average(error, weights=volumes)),
            "mean_absolute": float(np.average(np.abs(error), weights=volumes)),
            "rmse": math.sqrt(float(np.average(error ** 2, weights=volumes))),
            "max_absolute": float(np.max(np.abs(error)))}


def _regions(prepared):
    x, p = prepared.points.T
    radius, theta = np.hypot(x, p), np.arctan2(p, x)
    return {"full_domain": np.ones(x.size, bool), "near_exit_X_le_1": x <= 1,
            "near_lip_distance_le_0p1": np.hypot(x, p - .5) <= .1,
            "radius_0_to_1": radius < 1, "radius_1_to_5": (radius >= 1) & (radius < 5),
            "radius_ge_5": radius >= 5, "centerline_theta_le_5deg": theta <= np.deg2rad(5),
            "off_axis_theta_gt_5deg": theta > np.deg2rad(5),
            "angle_5_to_30deg": (theta > np.deg2rad(5)) & (theta <= np.deg2rad(30)),
            "angle_30_to_60deg": (theta > np.deg2rad(30)) & (theta <= np.deg2rad(60)),
            "angle_60_to_90deg": theta > np.deg2rad(60)}


def independent_physics_diagnostics(prepared, predicted_density, *, evaluator_options=None):
    """Compare unchanged Cai velocity/flux when corresponding DSMC data exist.

    Cartesian transverse velocity needs the original cell azimuth.  If the
    portable data provide only meridional coordinates, only the known axial
    component is used unless cylindrical components are explicitly supplied.
    All fluxes use the same physical sample points. The density-map flux is
    explicitly hybrid and is never used as a fit target.
    """
    metadata = prepared.case.metadata
    velocity = metadata.get("velocity_sqrt_beta0")
    cylindrical = metadata.get("cylindrical_velocity_sqrt_beta0")
    if velocity is None and cylindrical is None:
        return {"available": False, "velocity_corrected": False,
                "reason": "normalized DSMC velocity not available"}
    from plumetools.cai2012 import analytical
    from .comparison import mass_flux_per_steradian, projections

    options = evaluator_options or {}
    chunk = int(options.get("chunk_size", 128))
    cai_velocity = np.empty((len(prepared.points), 3))
    for start in range(0, len(prepared.points), chunk):
        stop = start + chunk
        x, rho = prepared.points[start:stop].T
        cai_velocity[start:stop] = analytical.moments(
            x, rho, 0., .5, prepared.case.speed_ratio,
            n_radial=options.get("n_radial", 48),
            n_azimuthal=options.get("n_azimuthal", 96))["U_sqrt_beta0"]
    ds_ur = None
    if cylindrical is not None:
        data = np.asarray(cylindrical, float)
        if data.shape != (prepared.mask.size, 2):
            raise ValueError("cylindrical_velocity_sqrt_beta0 must be raw Nx2 ux,ur")
        data = data[prepared.mask]
        ds_ux, ds_ur = data.T
    else:
        data = np.asarray(velocity, float)
        if data.shape != (prepared.mask.size, 3):
            raise ValueError("velocity_sqrt_beta0 must be raw Nx3 Cartesian components")
        data = data[prepared.mask]
        ds_ux = data[:, 0]
        if "cartesian_centres_m" in metadata:
            centres = np.asarray(metadata["cartesian_centres_m"], float)[prepared.mask]
            ds_ur = projections(centres, data)["ur"]
    valid = np.isfinite(ds_ux) & np.isfinite(cai_velocity[:, 0])
    report = {"available": True, "velocity_corrected": False,
              "velocity_units": "u sqrt(beta_exit)",
              "axial_velocity": _linear_metrics(ds_ux[valid], cai_velocity[valid, 0], prepared.volumes[valid]),
              "angular_flux": {"available": False,
                               "reason": "transverse azimuth/cylindrical velocity unavailable"}}
    if ds_ur is None:
        return report
    valid &= np.isfinite(ds_ur) & np.isfinite(cai_velocity[:, 1])
    report["cylindrical_radial_velocity"] = _linear_metrics(
        ds_ur[valid], cai_velocity[valid, 1], prepared.volumes[valid])
    radius = np.hypot(prepared.points[:, 0], prepared.points[:, 1])
    theta = np.arctan2(prepared.points[:, 1], prepared.points[:, 0])
    ds_uR = np.cos(theta)*ds_ux + np.sin(theta)*ds_ur
    cai_uR = np.cos(theta)*cai_velocity[:, 0] + np.sin(theta)*cai_velocity[:, 1]
    args = (radius, .5, prepared.case.speed_ratio)
    ds_flux = mass_flux_per_steradian(prepared.density, ds_uR, *args)
    cai_flux = mass_flux_per_steradian(prepared.baseline, cai_uR, *args)
    hybrid_flux = mass_flux_per_steradian(predicted_density, cai_uR, *args)
    shell = valid & (np.abs(radius - 10.) <= .5)
    rows = []
    integrals = {"dsmc": 0., "cai_baseline": 0., "hybrid": 0.}
    solid_angle = 0.
    for a in range(prepared.angle_edges.size - 1):
        mask = shell & (prepared.angle_bins == a)
        if not mask.any():
            rows.append({"angle_bin": a, "status": "missing"})
            continue
        omega = float(2*np.pi*(np.cos(prepared.angle_edges[a]) - np.cos(prepared.angle_edges[a + 1])))
        row = {"angle_bin": a, "status": "ok", "samples": int(mask.sum()),
               "angle_deg": float(np.rad2deg(np.average(theta[mask], weights=prepared.volumes[mask])))}
        for label, flux in (("dsmc", ds_flux), ("cai_baseline", cai_flux), ("hybrid", hybrid_flux)):
            row[label] = float(np.average(flux[mask], weights=prepared.volumes[mask]))
            integrals[label] += omega * row[label]
        solid_angle += omega
        rows.append(row)
    report["angular_flux"] = {
        "available": bool(shell.any()), "radius_over_D": 10., "shell_half_width_over_D": .5,
        "profiles": rows, "partial_matched_integrals": integrals,
        "covered_solid_angle_sr": solid_angle,
        "hybrid_definition": "corrected density times unchanged Cai velocity at physical points",
        "hybrid_is_conservative_prediction": False,
        "interpretation": "independent diagnostic; partial cell-shell integral is not a conservation proof"}
    return report


def _contour_loss(rows, *, penalty=1e6):
    """Equal cases, then occupied levels, then valid reference angular bins."""
    case_losses = []
    invalid = 0
    for name in sorted({row["case"] for row in rows}):
        levels = []
        for level in CONTOUR_LEVELS:
            selected = [row for row in rows if row["case"] == name and
                        row["level"] == level and row["dsmc_status"] == "ok"]
            if selected:
                invalid += sum(row["model_status"] != "ok" for row in selected)
                levels.append(np.mean([row["log_radius_residual"] ** 2
                                       if row["model_status"] == "ok" else penalty
                                       for row in selected]))
        if levels:
            case_losses.append(float(np.mean(levels)))
    return (float(np.mean(case_losses)) if case_losses else None), invalid


def evaluate(cases, mapping, *, floor=1e-10, angle_edges=None, radial_bins=80,
             evaluator_options=None, evaluator_factory=None, temperature=True,
             floor_multipliers=(0.1, 1.0, 10.0)):
    """Evaluate both diagnostics and fixed regional masks for any fitted track.

    Statistical uncertainty is reported unavailable unless the caller supplies
    demonstrably independent time blocks. Occupancy is never a confidence band.
    A model-invalid sample invalidates the entire corresponding objective; it
    cannot disappear from a metric. Predictions remain tied to physical cells.
    """
    started = time.perf_counter()
    prepared = [c if isinstance(c, PreparedCase) else prepare_case(
        c, floor=floor, angle_edges=angle_edges, radial_bins=radial_bins,
        evaluator_options=evaluator_options, evaluator_factory=evaluator_factory) for c in cases]
    preprocess_seconds = time.perf_counter() - started
    evaluation_started = time.perf_counter()
    reports, contour_rows, field_losses = [], [], []
    floor_scores = {str(float(floor * m)): [] for m in floor_multipliers}
    for data in prepared:
        predictor = _predictor(mapping, data.case, evaluator_options, evaluator_factory)
        failure = None
        try:
            pred = predictor.density(data.points, data.case.kn)
        except (ValueError, FloatingPointError, OverflowError) as exc:
            pred = np.full(data.density.shape, np.nan)
            failure = str(exc)
        model_valid = bool(np.all(np.isfinite(pred)) and np.all(pred >= 0))
        rows = contour_diagnostics(data, pred)
        contour_rows.extend(rows)
        lf = None
        if data.density.size and model_valid:
            log_residual = np.log(np.maximum(pred, floor)) - np.log(np.maximum(data.density, floor))
            lf = float(np.dot(data.weights, log_residual ** 2))
            field_losses.append(lf)
            for f in floor_scores:
                eps = float(f)
                floor_scores[f].append(float(np.dot(data.weights, (
                    np.log(np.maximum(pred, eps)) - np.log(np.maximum(data.density, eps))) ** 2)))
        elif data.density.size:
            field_losses.append(1e6)
        lc, invalid_contours = _contour_loss(rows)
        regional = {}
        for region, mask in _regions(data).items():
            regional[region] = _linear_metrics(data.density[mask], pred[mask], data.volumes[mask])
            regional[region]["stratified_log_MSE"] = (float(np.average(
                (np.log(np.maximum(pred[mask], floor)) -
                 np.log(np.maximum(data.density[mask], floor))) ** 2,
                weights=data.weights[mask])) if mask.any() and model_valid else None)
        temp_report = {"available": False, "coefficients_refitted": False}
        if temperature and data.case.temperature is not None:
            temp = data.case.temperature[data.mask]
            valid = np.isfinite(temp) & (temp >= 0)
            try:
                predicted_temp = predictor.temperature(data.points, data.case.kn)
            except (ValueError, FloatingPointError, OverflowError):
                predicted_temp = np.full(temp.shape, np.nan)
            if np.any(valid) and np.all(np.isfinite(predicted_temp[valid])):
                temp_report = {"available": True, "coefficients_refitted": False,
                               "meaning": "density-derived coordinate-map transfer test; T/Te",
                               "metrics": _linear_metrics(temp[valid], predicted_temp[valid], data.volumes[valid]),
                               "excluded_temperature_samples": int((~valid).sum()),
                               "regions": {name: _linear_metrics(temp[m & valid], predicted_temp[m & valid],
                                                                  data.volumes[m & valid])
                                           for name, m in _regions(data).items()}}
        try:
            jac = np.asarray(mapping.jacobian(data.points, data.case.kn))
            jacobian = {"minimum": float(np.min(jac)) if jac.size else None,
                        "maximum": float(np.max(jac)) if jac.size else None,
                        "positive": bool(np.all(jac > 0)) if jac.size else None}
        except (ValueError, FloatingPointError):
            jacobian = {"minimum": None, "maximum": None, "positive": False}
        reports.append({"name": data.case.name, "kn": float(data.case.kn),
                        "suitable": bool(data.case.suitable), "coverage": data.coverage,
                        "L_F": lf, "L_C": lc, "invalid_model_contours": invalid_contours,
                        "model_evaluation_valid": model_valid,
                        "model_evaluation_failure": failure,
                        "model_floor_active": int(np.count_nonzero(pred < floor)),
                        "linear": _linear_metrics(data.density, pred, data.volumes),
                        "regions": regional, "temperature_transfer": temp_report,
                        "jacobian": jacobian,
                        "independent_physics": independent_physics_diagnostics(
                            data, pred, evaluator_options=evaluator_options),
                        "uncertainty": {"statistical": "unavailable",
                                        "reason": "no independently established time-block refit supplied",
                                        "density_blocks_present": "density_blocks" in data.case.metadata,
                                        "occupancy_is_confidence_interval": False}})
    contour_loss, invalid = _contour_loss(contour_rows)
    return {"L_F": float(np.mean(field_losses)) if field_losses else None,
            "L_C": contour_loss, "invalid_model_contours": invalid,
            "cases": reports, "contours": contour_rows,
            "floor_sensitivity": {f: float(np.mean(v)) if v else None for f, v in floor_scores.items()},
            "timing": {"preprocessing_seconds": preprocess_seconds,
                       "evaluation_seconds": time.perf_counter() - evaluation_started},
            "interpretation": "descriptive coordinate-only fit; no conservative or causal claim"}


def select_cases(cases, *, suitable_only=False, holdout_kn=()):
    """Split whole-Kn holdouts BEFORE any DSMC-dependent preprocessing."""
    train, held, excluded = [], [], []
    for case in cases:
        if suitable_only and not case.suitable:
            excluded.append(case)
        elif any(math.isclose(float(case.kn), float(k), rel_tol=1e-9) for k in holdout_kn):
            held.append(case)
        else:
            train.append(case)
    return train, held, excluded


def _assert_compatible(cases):
    """A Kn-only map is conditional on fixed source/gas/geometry settings."""
    if len({case.name for case in cases}) != len(cases):
        raise ValueError("case names must be unique within a joint fit")
    if cases and any(not math.isclose(case.speed_ratio, cases[0].speed_ratio,
                                      rel_tol=1e-10, abs_tol=1e-12) for case in cases):
        raise ValueError("joint Kn-only fit requires one speed ratio; stratify cases first")
    strata = {str(case.metadata["stratum"]) for case in cases
              if case.metadata.get("stratum") is not None}
    if len(strata) > 1:
        raise ValueError("joint Kn-only fit mixes source/gas/geometry strata; stratify cases first")


@dataclass
class FitResult:
    mapping: CoordinateMap
    track: str
    success: bool
    message: str
    training_cases: list
    holdout_cases: list
    identifiable: bool
    timings: dict
    diagnostics: dict

    def to_dict(self):
        return {"mapping": self.mapping.to_dict(), "track": self.track,
                "success": self.success, "message": self.message,
                "training_cases": self.training_cases, "holdout_cases": self.holdout_cases,
                "identifiable": self.identifiable, "timing": self.timings,
                "diagnostics": self.diagnostics}


def fit(cases: Sequence[ScalingCase], mapping: CoordinateMap, *, track="F",
        holdout_kn=(), suitable_only=False, max_nfev=100, regularization=0.0,
        coefficient_bound=None, floor=1e-10, angle_edges=None, radial_bins=80,
        evaluator_options=None, evaluator_factory=None, evaluate_temperature=True):
    """Direct joint fit across Kn and coordinates for one separate objective.

    Invalid map construction/evaluation returns a finite penalty on the fixed
    residual vector, never a changed sample mask. Coefficient bounds and the
    map's global log-derivative bounds are both enforced. The regularizer is
    explicit ``regularization * sum(coefficients**2)``; the reported physical
    field/contour losses exclude it. No whole-Kn holdout is enabled by default.
    """
    if track not in ("F", "C"):
        raise ValueError("track must be F (field) or C (contour)")
    if max_nfev < 1 or regularization < 0:
        raise ValueError("max_nfev must be positive and regularization nonnegative")
    started = time.perf_counter()
    train, held, excluded = select_cases(cases, suitable_only=suitable_only,
                                        holdout_kn=holdout_kn)
    _assert_compatible(train)
    prepared = [prepare_case(c, floor=floor, angle_edges=angle_edges,
                            radial_bins=radial_bins, evaluator_options=evaluator_options,
                            evaluator_factory=evaluator_factory) for c in train]
    usable = [c for c in prepared if c.density.size]
    nparams = mapping.parameter_count
    fixed_contours = {c.case.name: [r for r in contour_diagnostics(c, c.baseline)
                                  if r["dsmc_status"] == "ok"] for c in usable}
    if track == "C":
        usable = [c for c in usable if fixed_contours[c.case.name]]
    preprocess_seconds = time.perf_counter() - started
    usable_ids = {id(c) for c in usable}
    omitted = [c.case.name for c in prepared if id(c) not in usable_ids]
    count = (sum(c.density.size for c in usable) if track == "F" else
             sum(len(fixed_contours[c.case.name]) for c in usable))
    fit_started = time.perf_counter()
    penalty = 1e3
    bound = float(mapping.log_bound if coefficient_bound is None else coefficient_bound)
    if bound <= 0 or not np.isfinite(bound):
        raise ValueError("coefficient_bound must be finite and positive")
    last = {"invalid_candidates": 0}

    def make_map(parameters):
        document = mapping.to_dict()
        document["coefficients"] = np.asarray(parameters).tolist()
        return CoordinateMap.from_dict(document)

    def residual(parameters):
        try:
            current = make_map(parameters)
            values = []
            for data in usable:
                prediction = _predictor(current, data.case, evaluator_options,
                                        evaluator_factory).density(data.points, data.case.kn)
                if np.any(~np.isfinite(prediction)) or np.any(prediction < 0):
                    raise ValueError("unavailable mapped Cai evaluation")
                if track == "F":
                    values.extend((np.sqrt(data.weights / len(usable)) * (
                        np.log(np.maximum(prediction, floor)) -
                        np.log(np.maximum(data.density, floor)))).tolist())
                else:
                    rows = contour_diagnostics(data, prediction)
                    indexed = {(r["angle_bin"], r["level"]): r for r in rows}
                    reference = fixed_contours[data.case.name]
                    levels = sorted({r["level"] for r in reference})
                    per_level = {q: sum(r["level"] == q for r in reference) for q in levels}
                    for row in reference:
                        current_row = indexed[(row["angle_bin"], row["level"])]
                        value = (current_row["log_radius_residual"]
                                 if current_row["model_status"] == "ok" else penalty)
                        values.append(value / math.sqrt(len(usable) * len(levels) * per_level[row["level"]]))
            if regularization:
                values.extend((math.sqrt(regularization) * np.asarray(parameters)).tolist())
            return np.asarray(values)
        except (ValueError, FloatingPointError, OverflowError):
            last["invalid_candidates"] += 1
            # Fixed dimension even when an entire mapped field is unavailable.
            return np.full(count + (nparams if regularization else 0),
                           penalty + float(np.dot(parameters, parameters)))

    identified, success, message, objective, nfev, rank = False, False, "", None, 0, 0
    fitted = mapping
    if not usable:
        message = "No computable cases for this track; no coefficients fitted."
    elif nparams == 0:
        identified, success = True, True
        message = "Cai baseline: no fitted parameters."
        objective = float(np.dot(residual([]), residual([])))
    elif count < nparams:
        message = f"Only {count} residuals for {nparams} coefficients; model not identified."
    else:
        solution = least_squares(residual, mapping.coefficients, bounds=(-bound, bound),
                                 max_nfev=int(max_nfev), xtol=1e-8, ftol=1e-8, gtol=1e-8)
        fitted = make_map(solution.x)
        nfev = int(solution.nfev)
        data_jac = solution.jac[:count]
        singular = np.linalg.svd(data_jac, compute_uv=False)
        threshold = singular[0] * max(data_jac.shape) * np.finfo(float).eps if singular.size else 0
        rank = int(np.count_nonzero(singular > max(threshold, 1e-8)))
        identified = rank == nparams
        success = bool(solution.success and identified and solution.cost < penalty ** 2 / 2)
        message = str(solution.message)
        if not identified:
            message += f" Data Jacobian rank {rank}/{nparams}; coefficients not identified."
        objective = float(np.dot(solution.fun, solution.fun))
    fit_seconds = time.perf_counter() - fit_started
    evaluation_started = time.perf_counter()
    metrics = evaluate(prepared, fitted, floor=floor, evaluator_options=evaluator_options,
                       evaluator_factory=evaluator_factory, temperature=evaluate_temperature)
    min_kn = min((c.kn for c in train), default=None)
    max_kn = max((c.kn for c in train), default=None)
    validation = [{"name": c.name, "kn": float(c.kn),
                   "regime": ("interpolation" if min_kn is not None and min_kn < c.kn < max_kn else "extrapolation")}
                  for c in held]
    diagnostic = {"objective": objective, "regularization": float(regularization),
                  "regularization_meaning": "explicit coefficient-size penalty; empirical, not physics",
                  "coefficient_bound": bound, "optimizer_nfev": nfev,
                  "data_jacobian_rank": rank, "parameter_count": nparams,
                  "residual_count": count, "omitted_no_computable_objective": omitted,
                  "excluded_unsuitable": [c.name for c in excluded],
                  "invalid_candidates": last["invalid_candidates"],
                  "bound_diagnostics": fitted.bound_diagnostics(),
                  "coefficient_box_active": bool(np.any(np.isclose(np.abs(fitted.coefficients), bound,
                                                                            atol=1e-5, rtol=0))),
                  "training_metrics": metrics, "holdout": validation,
                  "predictive_validation": bool(held),
                  "preprocessing_policy": "whole-Kn split before DSMC-dependent preparation",
                  "uncertainty": "statistical uncertainty unavailable; see numerical sensitivities"}
    return FitResult(fitted, track, success, message, [c.name for c in train],
                     [c.name for c in held], identified,
                     {"preprocessing_seconds": preprocess_seconds, "fitting_seconds": fit_seconds,
                      "evaluation_seconds": time.perf_counter() - evaluation_started}, diagnostic)


def sensitivity_report(cases, mapping, *, evaluator_options=None,
                       evaluator_factory=None, floor=1e-10,
                       angle_widths_deg=(2.5, 5.0), radial_bin_counts=(40, 80)):
    """Numerical diagnostic sensitivities with fixed coefficients (no tuning).

    These quantify diagnostic variation, not a DSMC statistical interval. All
    settings are predeclared, and callers must pass training cases only when a
    whole-Kn holdout is enabled. Bound/regularization coefficient sensitivity
    requires separate configured refits; this function does not imply it was run.
    """
    rows = []
    for width in angle_widths_deg:
        if width <= 0:
            raise ValueError("positive angular width required")
        edges = np.deg2rad(np.linspace(0, 90, int(math.ceil(90 / width)) + 1))
        for bins in radial_bin_counts:
            metrics = evaluate(cases, mapping, floor=floor, angle_edges=edges,
                               radial_bins=bins, evaluator_options=evaluator_options,
                               evaluator_factory=evaluator_factory, temperature=False)
            rows.append({"angular_width_deg": float(width), "radial_bins": int(bins),
                         "L_F": metrics["L_F"], "L_C": metrics["L_C"],
                         "floor_sensitivity": metrics["floor_sensitivity"]})
    return {"kind": "fixed-coefficient numerical diagnostic sensitivity",
            "statistical_confidence_interval": False, "settings": rows}


def block_refits(cases, mapping, *, track="F", n_resamples=8, seed=1729,
                 max_nfev=30, holdout_kn=(), suitable_only=False,
                 min_block_to_correlation_ratio=5., **fit_options):
    """Opt-in joint time-block resampling with an explicit correlation gate.

    Every training case must supply ``density_blocks`` shaped (blocks, cells),
    ``independent_blocks=True``, a nonempty ``temporal_correlation_method``,
    ``block_duration_s``, and positive ``integrated_autocorrelation_time_s``.
    A block must span at least five correlation times by default. These are
    provenance declarations, not quantities estimated from particle occupancy.
    One block is independently drawn for each case per joint refit. Reported
    parameter spread is a small descriptive resampling ensemble, not a claimed
    confidence interval. Held-out Kn values never enter blocks or tuning.
    """
    if n_resamples < 2 or max_nfev < 1 or min_block_to_correlation_ratio <= 0:
        raise ValueError("at least two resamples, positive iterations and correlation ratio required")
    training, held, excluded = select_cases(cases, suitable_only=suitable_only,
                                           holdout_kn=holdout_kn)
    reasons, blocks = [], {}
    for case in training:
        m = case.metadata
        array = np.asarray(m.get("density_blocks", []), float)
        if array.ndim != 2 or array.shape[0] < 2 or array.shape[1] != case.density.size:
            reasons.append(f"{case.name}: two or more matched density blocks unavailable")
            continue
        if not np.all(np.isfinite(array)) or np.any(array < 0):
            reasons.append(f"{case.name}: invalid density-block samples")
            continue
        duration, correlation = m.get("block_duration_s"), m.get("integrated_autocorrelation_time_s")
        try:
            ratio = float(duration)/float(correlation)
            correlated = (not np.isfinite(ratio) or float(correlation) <= 0 or
                          ratio < min_block_to_correlation_ratio)
        except (TypeError, ValueError, ZeroDivisionError):
            correlated = True
        if (m.get("independent_blocks") is not True or
                not str(m.get("temporal_correlation_method", "")).strip() or correlated):
            reasons.append(f"{case.name}: temporal independence/correlation-duration criterion not established")
            continue
        blocks[case.name] = array
    if not training:
        reasons.append("no eligible training cases")
    if reasons:
        return {"available": False, "reasons": reasons, "resamples_run": 0,
                "holdout_cases": [c.name for c in held],
                "occupancy_is_confidence_interval": False}
    rng, rows = np.random.default_rng(seed), []
    started = time.perf_counter()
    for draw in range(int(n_resamples)):
        sampled, indices = [], {}
        for case in training:
            index = int(rng.integers(len(blocks[case.name])))
            indices[case.name] = index
            sampled.append(replace(case, density=blocks[case.name][index],
                                   metadata={key: value for key, value in case.metadata.items()
                                             if key != "density_blocks"}))
        result = fit(sampled, mapping, track=track, max_nfev=max_nfev,
                     evaluate_temperature=False, **fit_options)
        rows.append({"draw": draw, "block_indices": indices, "success": result.success,
                     "identifiable": result.identifiable,
                     "coefficients": result.mapping.coefficients.tolist()})
    valid = [row["coefficients"] for row in rows if row["success"] and row["identifiable"]]
    coefficients = np.asarray(valid, float)
    return {"available": len(valid) >= 2, "resamples_run": len(rows),
            "successful_refits": len(valid), "refits": rows,
            "coefficient_mean": coefficients.mean(axis=0).tolist() if len(valid) else None,
            "coefficient_sample_std": coefficients.std(axis=0, ddof=1).tolist() if len(valid) >= 2 else None,
            "seed": int(seed), "max_nfev": int(max_nfev),
            "minimum_block_to_correlation_ratio": float(min_block_to_correlation_ratio),
            "timing_seconds": time.perf_counter() - started,
            "holdout_cases": [c.name for c in held],
            "interpretation": "descriptive time-block resampling spread; not a confidence interval",
            "occupancy_is_confidence_interval": False}
