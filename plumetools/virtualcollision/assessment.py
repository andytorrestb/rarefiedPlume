"""Read-only case inventory and explicitly DSMC-referenced assessment metrics.

Legacy ``masked_metrics`` deliberately retains its analytical denominator.
These metrics use a supplied, shared mask and a DSMC global norm instead.
Nothing in this module generates cases or invokes OpenFOAM.
"""
from __future__ import annotations

import csv
import json
import math
import platform
import re
import subprocess
from pathlib import Path

import numpy as np
import yaml

from plumetools.cai2012.gas import (
    VhsSpecies, knudsen_from_number_density, number_density_from_kn,
    speed_from_speed_ratio,
)
from plumetools.foamio.times import time_names

REQUESTED_KN = (100., 10., 1., .1, .01, .001, .0001)


def load_yaml(path):
    path = Path(path)
    if not path.is_file():
        return {}
    try:
        value = yaml.safe_load(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {"unreadable": "not a mapping"}
    except (OSError, yaml.YAMLError) as exc:
        return {"unreadable": str(exc)}


def json_safe(value):
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path, document):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_safe(document), indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def write_rows(path, rows, *, columns=()):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys([*columns, *(key for row in rows for key in row)]))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (json.dumps(json_safe(v), sort_keys=True)
                                 if isinstance(v, (dict, list, tuple)) else v)
                             for k, v in json_safe(row).items()})


def source_revision(root):
    try:
        return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"],
                                       text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def dsmc_reference_metrics(reference, cai, corrected, *, weights, common_mask,
                           scale=None, baseline_floor=1e-12):
    """Score BOTH predictions on precisely the same requested reference mask.

    A failed corrected value invalidates the corrected score for that mask;
    it is never allowed to remove a difficult location from the baseline score.
    ``weights`` are physical volumes (or stated bin measures); a fixed ``scale``
    replaces a global DSMC RMS norm for components crossing zero.
    """
    ref, base, vc, w = [np.asarray(a, float).ravel()
                        for a in (reference, cai, corrected, weights)]
    mask = np.asarray(common_mask, bool).ravel()
    if not (ref.shape == base.shape == vc.shape == w.shape == mask.shape):
        raise ValueError("all metric arrays must have identical shapes")
    if np.any(w[mask] <= 0) or np.any(~np.isfinite(w[mask])):
        raise ValueError("metric weights must be finite and positive")
    use = mask & np.isfinite(ref)
    output = {"requested_points": int(mask.sum()), "reference_points": int(use.sum()),
              "common_mask": True, "weight_sum": float(w[use].sum()),
              "normalization": "DSMC weighted global RMS" if scale is None
                               else "fixed physical scale", "scale": scale}
    if not use.any():
        return {**output, "status": "empty-common-mask", "baseline_error": None,
                "corrected_error": None, "error_difference": None,
                "relative_improvement": None}
    denom = (float(np.sum(w[use] * ref[use] ** 2)) if scale is None
             else float(w[use].sum() * float(scale) ** 2))
    errors = []
    for name, values in (("baseline", base), ("corrected", vc)):
        failed = int(np.count_nonzero(~np.isfinite(values[use])))
        output[name + "_nonfinite"] = failed
        error = (None if failed or denom <= 0 else
                 float(np.sqrt(np.sum(w[use] * (values[use] - ref[use]) ** 2) / denom)))
        output[name + "_error"] = error
        output[name + "_signed_bias"] = (None if failed else
            float(np.sum(w[use] * (values[use] - ref[use])) / w[use].sum()))
        errors.append(error)
    b, c = errors
    output["error_difference"] = None if b is None or c is None else c - b
    output["relative_improvement"] = (None if b is None or c is None or b <= baseline_floor
                                       else (b - c) / b)
    output["status"] = ("failed-prediction-on-common-mask" if None in errors else "scored")
    return output


def time_block_standard_error(series, *, minimum_blocks=4):
    """Conservative blocking estimate for an INSTANTANEOUS stationary series.

    Blocks double until at most ``minimum_blocks`` remain; report the largest
    SE across those block sizes, preserving autocorrelation inflation. Running
    cumulative averages must not be supplied as independent samples.
    """
    values = np.asarray(series, float)
    if values.ndim != 1 or not np.all(np.isfinite(values)):
        raise ValueError("a finite one-dimensional instantaneous series is required")
    estimates = []
    size = 1
    while len(values) // size >= minimum_blocks:
        count = len(values) // size
        blocks = values[:count * size].reshape(count, size).mean(axis=1)
        estimates.append({"block_size": size, "blocks": count,
                          "standard_error": float(blocks.std(ddof=1) / np.sqrt(count))})
        size *= 2
    return {"status": "available" if estimates else "insufficient-instantaneous-samples",
            "method": "maximum SE over doubling non-overlapping block means",
            "standard_error": max((x["standard_error"] for x in estimates), default=None),
            "blocks": estimates}


def _uncomment(text):
    return re.sub(r"/\*.*?\*/|//[^\n]*", "", text, flags=re.S)


def _dictionary_value(text, key):
    match = re.search(r"\b" + re.escape(key) + r"\s+([^;{}]+);", text)
    return match.group(1).strip() if match else None


def _number(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def inventory(root, *, families=("cai2012", "cai2012-kn-extension")):
    """Discover ignored/generated case folders; templates are recorded separately.

    Configuration-derived Kn and diagnostic measured exit Kn are distinct. A
    dictionary corroborates the configured kernel, not its solver provenance.
    """
    root = Path(root).resolve()
    rows = []
    for family in families:
        family_root = root / "cases" / family
        if not family_root.is_dir():
            continue
        for path in sorted(family_root.rglob("case.yaml")):
            case = path.parent
            relative = case.relative_to(root).as_posix()
            cfg = load_yaml(path)
            ex, gas, nozzle = [cfg.get(k, {}) for k in ("exit", "gas", "nozzle")]
            geom, dsmc = [cfg.get(k, {}) for k in ("geometry", "dsmc")]
            interface = load_yaml(case / "interface.yaml")
            source_type = "interface" if interface.get("enabled", bool(interface)) else "direct-nozzle"
            template = any(part.lower().startswith("basecase") for part in case.parts)
            nominal = _number(ex.get("knudsen", cfg.get("meta", {}).get("knudsen")))
            result = case / "results"
            diag = load_yaml(result / "diagnostics_summary.yaml")
            station = load_yaml(result / "stationarity.yaml")
            compare = load_yaml(result / "comparison.yaml")
            convergence = load_yaml(result / "convergence.yaml")
            available = sorted(p.name for p in result.glob("*.csv")) if result.is_dir() else []
            frames = [case / name for name in time_names(case)]
            field_times = [p.name for p in frames if (p / "rhoNMean").is_file()]
            geometry_times = [p.name for p in frames if (p / "C").is_file() and (p / "V").is_file()]
            record = {"path": relative, "family": family, "source_type": source_type,
                      "nominal_kn": nominal, "original_six_point_inventory": nominal != .01,
                      "additional_cai_family_comparison": nominal == .01,
                      "kn_length_convention": ex.get("characteristic_length"),
                      "mean_free_path_convention": gas.get("mean_free_path_convention", "vhs"),
                      "source_distribution": "uniform drifting Maxwellian; outward one-way flux",
                      "gas": gas, "collision_law": dsmc.get("binary_collision_model"),
                      "collisions_enabled": dsmc.get("collisions_enabled"),
                      "exit": ex, "geometry": {**geom, **nozzle}, "available_csv": available,
                      "available_field_times": field_times, "available_geometry_times": geometry_times,
                      "exported_time": compare.get("case", {}).get("time", diag.get("time")),
                      "sampling": diag.get("sampling"),
                      "stationarity": station.get("assessment", {"verdict": "unavailable"}),
                      "occupancy": diag.get("particles"), "diagnostic_validity": diag.get("validity"),
                      "resolution": diag.get("resolution"), "measured_exit_diagnostics": diag.get("exit_plane"),
                      "production": compare.get("production", {}),
                      "non_production_marker": (case / "NON-PRODUCTION").is_file(),
                      "dsmc_uncertainty": {"status": "unavailable",
                          "reason": "averaged exports and running-average frame metrics are not independent instantaneous blocks"},
                      "dsmc_tensor": "unavailable unless tensor raw fields are explicitly supplied",
                      "dsmc_production_cost": load_yaml(result / "cost-model.yaml") or None}
            reasons = []
            record["solver_runtime_provenance"] = {"status": "unavailable",
                "reason": "configured law and dictionary headers do not identify the historic executable"}
            for log_name in ("log.dsmcFoam", "output.log"):
                log_path = case / log_name
                if log_path.is_file():
                    with log_path.open(encoding="utf-8", errors="replace") as stream:
                        header = "".join(stream.readline() for _ in range(100))
                    build = re.search(r"^Build\s*:\s*(.+)$", header, flags=re.M)
                    executable = re.search(r"^Exec\s*:\s*(.+)$", header, flags=re.M)
                    if build:
                        record["solver_runtime_provenance"] = {"status": "log-recorded",
                            "path": log_path.relative_to(root).as_posix(), "build": build.group(1),
                            "exec": executable.group(1) if executable else None,
                            "kernel_source_match": "not established from executable hash"}
                        break
            if cfg.get("model", "cai2012_circular_plume") != "cai2012_circular_plume":
                reasons.append("unsupported-source-model")
            if geom.get("plume_axis", [1, 0, 0]) != [1, 0, 0] or geom.get("symmetry_mode", "none") != "none":
                reasons.append("unsupported-geometry")
            if dsmc.get("collisions_enabled") is False:
                reasons.append("unsupported-reference-collision-state")
            try:
                for key in ("mass_kg", "diameter_m", "omega", "t_ref_K"):
                    if key not in gas:
                        raise ValueError("missing gas." + key)
                species = VhsSpecies(name=gas.get("species_name", "Ar"),
                                     mass_kg=float(gas.get("mass_kg", 6.63e-26)),
                                     diameter_m=float(gas.get("diameter_m", 4.17e-10)),
                                     omega=float(gas.get("omega", .74)),
                                     t_ref_K=float(gas.get("t_ref_K", 273)),
                                     internal_degrees_of_freedom=int(gas.get("internal_degrees_of_freedom", 0)))
                length = float(nozzle["diameter_m"])
                if ex.get("characteristic_length", "diameter") != "diameter":
                    raise ValueError("unsupported Kn length convention")
                density = number_density_from_kn(nominal, length, species, float(ex["T0_K"]),
                                                convention=record["mean_free_path_convention"])
                record["exit_density_m3"] = density
                record["exit_speed_m_s"] = speed_from_speed_ratio(float(ex["speed_ratio"]), species, float(ex["T0_K"]))
                record["verified_kn_from_configured_density"] = knudsen_from_number_density(
                    density, length, species, float(ex["T0_K"]), convention=record["mean_free_path_convention"])
                record["kn_verification"] = "configuration inverse check; independent injected density unavailable"
                boundary = case / "0" / ("boundaryNumberDensity_" + gas.get("species_name", "Ar"))
                record["verified_kn_from_injected_density"] = None
                record["injected_source_fields"] = {"status": "unavailable", "path": boundary.relative_to(root).as_posix()}
                if boundary.is_file() and source_type == "direct-nozzle":
                    text = _uncomment(boundary.read_text(encoding="utf-8", errors="replace"))
                    patch_name = cfg.get("mesh", {}).get("patch_names", {}).get("nozzle", "nozzle")
                    patch = re.search(r"\b"+re.escape(patch_name)+r"\s*\{([^{}]*)\}", text)
                    value = re.search(r"\bvalue\s+uniform\s+([^;]+);", patch.group(1)) if patch else None
                    injected = _number(value.group(1)) if value else None
                    if injected is None or injected <= 0:
                        reasons.append("unsupported-nonuniform-or-unreadable-exit-density")
                    else:
                        actual_kn = knudsen_from_number_density(injected, length, species, float(ex["T0_K"]),
                                                               convention=record["mean_free_path_convention"])
                        record["verified_kn_from_injected_density"] = actual_kn
                        record["injected_source_fields"] = {"status": "uniform-density-verified", "density_m3": injected,
                                                             "path": boundary.relative_to(root).as_posix()}
                        record["kn_verification"] = "independent injected nozzle-patch density checked using configured equilibrium MFP convention"
                        if not math.isclose(actual_kn, nominal, rel_tol=1e-5):
                            reasons.append("unsupported-injected-kn-mismatch")
            except (ValueError, TypeError, KeyError) as exc:
                reasons.append("unsupported-metadata: " + str(exc))
            dictionary = case / "constant" / "dsmcProperties"
            if dictionary.is_file():
                text = _uncomment(dictionary.read_text(encoding="utf-8", errors="replace"))
                verified = {k: _dictionary_value(text, k) for k in
                            ("BinaryCollisionModel", "Tref", "mass", "diameter", "omega", "internalDegreesOfFreedom", "InflowBoundaryModel")}
                type_ids = re.search(r"\btypeIdList\s*\(([^)]*)\)\s*;", text)
                verified["typeIdList"] = type_ids.group(1).split() if type_ids else None
                patches = re.search(r"\bpatches\s*\(([^)]*)\)\s*;", text)
                verified["inflow_patches"] = patches.group(1).split() if patches else None
                record["collision_dictionary"] = verified
                if verified["BinaryCollisionModel"] != dsmc.get("binary_collision_model"):
                    reasons.append("unsupported-collision-dictionary-mismatch")
                for key, expected in (("Tref", gas.get("t_ref_K")), ("mass", gas.get("mass_kg")),
                                      ("diameter", gas.get("diameter_m")), ("omega", gas.get("omega"))):
                    actual = _number(verified[key])
                    if actual is not None and expected is not None and not math.isclose(actual, float(expected), rel_tol=1e-6):
                        reasons.append("unsupported-collision-parameter-mismatch:" + key)
                if verified["typeIdList"] is not None and verified["typeIdList"] != [gas.get("species_name", "Ar")]:
                    reasons.append("unsupported-non-single-species-dictionary")
                if source_type == "direct-nozzle" and verified["InflowBoundaryModel"] is not None:
                    if verified["InflowBoundaryModel"] != "plumeFieldInflow" or verified["inflow_patches"] != ["nozzle"]:
                        reasons.append("unsupported-source-boundary-dictionary")
            else:
                record["collision_dictionary"] = {"status": "unavailable"}
            if dsmc.get("binary_collision_model") != "VariableHardSphere" or gas.get("internal_degrees_of_freedom", 0) != 0:
                reasons.append("unsupported-collision-law")
            if source_type == "interface":
                reasons.append("interface-source-changes-boundary-value-problem")
            if template:
                reasons.append("template-not-a-run")
            ready = compare.get("production", {}).get("production_ready")
            if nominal in (.001, .0001) and source_type == "direct-nozzle":
                reasons.append("known-inadequate-direct-exit-resolution-and-timestep")
            if ready is False or record["non_production_marker"]:
                reasons.append("non-production-reference")
            verdict = station.get("assessment", {}).get("verdict", "unavailable")
            if verdict != "stable-window-demonstrated":
                reasons.append("stationarity-" + verdict)
            if available and not (field_times and geometry_times):
                if (case / 'system' / 'blockMeshDict').is_file() and cfg.get('mesh', {}).get('type') == 'graded_cartesian':
                    record['compact_geometry'] = 'candidate for checked reconstruction; use compare_exports.py'
                    reasons.append('compact-exports-require-verified-mesh-reconstruction')
                else:
                    reasons.append("compact-exports-missing-cell-volumes-and-unaggregated-moments")
            if not available and not field_times:
                reasons.append("missing-dsmc-results")
            inadequate = any("inadequate" in s or "non-production" in s for s in reasons)
            unsupported = any(s.startswith("unsupported") for s in reasons) or source_type == "interface"
            quality = ("unsupported" if unsupported else "unavailable" if not available and not field_times
                       else "exploratory" if inadequate or verdict != "stable-window-demonstrated"
                       else "reference-quality")
            record.update({"quality": quality, "primary_included": quality == "reference-quality"
                           and bool(field_times and geometry_times), "reasons": reasons,
                           "kind": "template" if template else "case"})
            rows.append(record)
    for kn in REQUESTED_KN:
        if not any(r["nominal_kn"] == kn and r["kind"] == "case"
                   and r["source_type"] == "direct-nozzle" for r in rows):
            rows.append({"path": None, "family": "cai2012" if kn in (100, .1, .01) else "cai2012-kn-extension",
                         "source_type": "direct-nozzle", "nominal_kn": kn,
                         "original_six_point_inventory": kn != .01, "additional_cai_family_comparison": kn == .01,
                         "quality": "unavailable", "primary_included": False,
                         "reasons": ["requested-case-not-discovered"], "kind": "missing"})
    return {"source_revision": source_revision(root), "repository": str(root),
            "requested_kn": REQUESTED_KN, "cases": rows,
            "host": {"platform": platform.platform(), "processor": platform.processor()},
            "policy": "read-only aggregation across families; legacy generation/run gates unchanged"}


def aggregate_reference_bins(sampled, x_edges, r_edges, species_mass_kg):
    """Volume-integrate DSMC cell raw moments into cylindrical assessment bins.

    Cell centres choose a bin, so boundary cells introduce an explicitly measured
    finite-volume quadrature error. The scalar second raw moment is recoverable;
    its tensor is not inferred from scalar temperature. Cylindrical components
    are rotated before integration, matching the model's local ring basis.
    """
    from plumetools.cai2012.constants import BOLTZMANN_J_PER_K
    from plumetools.cai2012ext.comparison import projections
    x_edges, r_edges = np.asarray(x_edges, float), np.asarray(r_edges, float)
    if np.any(np.diff(x_edges) <= 0) or np.any(np.diff(r_edges) <= 0):
        raise ValueError("bin edges must increase strictly")
    points, volume = np.asarray(sampled.centres), np.asarray(sampled.volumes)
    radius = np.hypot(points[:, 1], points[:, 2])
    ix = np.searchsorted(x_edges, points[:, 0], side="right") - 1
    ir = np.searchsorted(r_edges, radius, side="right") - 1
    nx, nr = len(x_edges) - 1, len(r_edges) - 1
    inside = (ix >= 0) & (ix < nx) & (ir >= 0) & (ir < nr)
    if np.any(~np.isfinite(volume)) or np.any(volume <= 0):
        raise ValueError("physical cell volumes are required")
    ids = ix[inside] * nr + ir[inside]
    shape = nx * nr
    def total(value):
        return np.bincount(ids, weights=volume[inside] * np.asarray(value)[inside], minlength=shape)
    volumes = total(np.ones(len(points)))
    count = np.bincount(ids, minlength=shape)
    number = np.asarray(sampled.number_density)
    proj = projections(points, sampled.velocity)
    velocity = np.column_stack([proj[k] for k in ("ux", "ur", "uphi")])
    m0 = total(number)
    m1 = np.column_stack([total(number * velocity[:, i]) for i in range(3)])
    trace = total(number * (np.sum(np.asarray(sampled.velocity) ** 2, axis=1)
                           + 3 * BOLTZMANN_J_PER_K / species_mass_kg * sampled.temperature))
    with np.errstate(divide="ignore", invalid="ignore"):
        n = m0 / volumes
        u = m1 / m0[:, None]
        temperature = (trace / m0 - np.sum(u*u, axis=1)) * species_mass_kg / (3 * BOLTZMANN_J_PER_K)
    # A geometric sensitivity indicator, never a confidence interval: cells whose
    # half-diagonal reaches a bin boundary may straddle adjacent bins.
    half = .5 * np.sqrt(3) * np.cbrt(volume[inside])
    near = (np.minimum(points[inside, 0] - x_edges[ix[inside]],
                       x_edges[ix[inside] + 1] - points[inside, 0]) < half)
    near |= (np.minimum(radius[inside] - r_edges[ir[inside]],
                        r_edges[ir[inside] + 1] - radius[inside]) < half)
    boundary_volume = np.bincount(ids, weights=volume[inside]*near, minlength=shape)
    return {"M0": n, "M1": m1 / np.where(volumes > 0, volumes, np.nan)[:, None],
            "M2_trace": trace / np.where(volumes > 0, volumes, np.nan),
            "n": n, "U": u, "T_trans": temperature, "volumes": volumes,
            "cell_count": count,
            "boundary_cell_volume_fraction": boundary_volume / np.where(volumes > 0, volumes, np.nan),
            "boundary_indicator_status": "heuristic cube-volume half-diagonal; not a bound for anisotropic cells",
            "tensor_status": "unavailable: only scalar second-moment trace reconstructable"}


def contour_location_discrepancy(x, r, reference, prediction, levels):
    """Symmetric nearest-contour RMS distance in supplied coordinate units.

    No overlapping level returns an explicit missing result. Contours use the
    same grid; this measures location only, not topology or enclosed mass.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.spatial import cKDTree
    result = []
    for level in levels:
        contours = []
        for values in (reference, prediction):
            finite = np.asarray(values)[np.isfinite(values)]
            if not finite.size or not finite.min() < level < finite.max():
                contours.append(None)
                continue
            fig, ax = plt.subplots()
            cs = ax.contour(x, r, values, levels=[level])
            segments = [seg for seg in cs.allsegs[0] if len(seg)]
            contours.append(np.vstack(segments) if segments else None)
            plt.close(fig)
        a, b = contours
        if a is None or b is None:
            result.append({"level": level, "status": "missing-or-nonoverlapping-contour", "rms_distance": None})
        else:
            distance = np.concatenate([cKDTree(a).query(b)[0], cKDTree(b).query(a)[0]])
            result.append({"level": level, "status": "available",
                           "rms_distance": float(np.sqrt(np.mean(distance**2)))})
    return result


def run_reference_assessment(cfg, inventory_doc, output, *, smoke=True, plots=True):
    """Assess compatible existing full fields; compact exports are never invented.

    Uses identical cylindrical bins for density and raw-moment velocity/scalar
    temperature. DSMC cells are a finite-volume quadrature of each bin; the
    boundary-volume fraction reports sensitivity to centre-assigned cells.
    Interface sources are categorically excluded from this nozzle comparison.
    """
    from .source import DiskSource
    from .kernel import VHSKernel
    from .model import NumericalOptions, VirtualCollisionModel
    from plumetools.cai2012 import post
    from plumetools.cai2012ext import comparison
    root = Path(cfg.get("_root", inventory_doc["repository"]))
    output = Path(output)
    options_cfg = cfg.get("smoke" if smoke else "full", {})
    options = NumericalOptions(**{k: v for k, v in options_cfg.items()
                                  if k in NumericalOptions.__dataclass_fields__})
    settings = cfg.get("assessment_bins", {})
    nx, nr = int(settings.get("x_count", 20)), int(settings.get("r_count", 20))
    skipped, processed, metrics = [], [], []
    for row in inventory_doc["cases"]:
        if row["kind"] != "case" or row["source_type"] != "direct-nozzle" or row["quality"] == "unsupported":
            continue
        if smoke and row["nominal_kn"] not in options_cfg.get("kn", [10.]):
            continue
        reason = None
        if not row.get("available_field_times"):
            reason = "missing-unaggregated-DSMC-fields"
        elif not row.get("available_geometry_times"):
            reason = "missing-cell-centres-or-physical-volumes"
        if reason:
            skipped.append({"path": row["path"], "reason": reason,
                            "available_exports": row.get("available_csv", []),
                            "plots": "unavailable for matched three-way estimator"})
            continue
        case = root / row["path"]
        case_cfg = load_yaml(case / "case.yaml")
        gas = row["gas"]
        species = VhsSpecies(name=gas.get("species_name", "Ar"), mass_kg=gas["mass_kg"],
                             diameter_m=gas["diameter_m"], omega=gas["omega"],
                             t_ref_K=gas["t_ref_K"], internal_degrees_of_freedom=0)
        source = DiskSource(row["geometry"]["diameter_m"] / 2, row["exit_density_m3"],
                            row["exit"]["T0_K"], row["exit"]["speed_ratio"], species)
        model = VirtualCollisionModel(source, VHSKernel(species), options)
        d = 2 * source.radius_m
        xedges = np.linspace(0, d*min(float(settings.get("x_max_over_D", 10)),
                                    float(case_cfg.get("geometry", {}).get("x_max_over_D", 10))), nx + 1)
        redges = np.linspace(0, d*float(settings.get("r_max_over_D", 10)), nr + 1)
        xx, rr = np.meshgrid(.5*(xedges[:-1]+xedges[1:]), .5*(redges[:-1]+redges[1:]), indexing="ij")
        xb = np.repeat(np.column_stack([xedges[:-1], xedges[1:]]), nr, axis=0)
        rb = np.tile(np.column_stack([redges[:-1], redges[1:]]), (nx, 1))
        try:
            sampled = post.read_case(case, species.mass_kg, time=row["available_field_times"][-1],
                                     mesh_time=row["available_geometry_times"][-1])
            measured = aggregate_reference_bins(sampled, xedges, redges, species.mass_kg)
            model.preprocess()
            corrected = model.evaluate_bins(xb, rb)
            cai = model.evaluate_bins(xb, rb, collisions=False)
        except (ValueError, OSError, RuntimeError) as exc:
            skipped.append({"path": row["path"], "reason": "failed-numerical-evaluation", "detail": str(exc)})
            continue
        target = output / row["path"].replace("/", "__")
        target.mkdir(parents=True, exist_ok=True)
        arrays = {"x_bounds": xb, "r_bounds": rb, "DSMC_n": measured["n"],
                  "DSMC_U": measured["U"], "DSMC_T_trans": measured["T_trans"],
                  "DSMC_M2_trace": measured["M2_trace"], "DSMC_bin_volumes": measured["volumes"],
                  "DSMC_boundary_cell_volume_fraction": measured["boundary_cell_volume_fraction"]}
        for name, fields in (("Cai", cai), ("VC", corrected)):
            arrays.update({name+"_"+key: fields[key] for key in
                           ("M0", "M1", "M2", "n", "U", "T_tensor", "T_trans",
                            "P_tensor", "valid", "tensor_eigenvalues_K", "tensor_symmetric")})
            for key, values in fields.get("standard_error", {}).items():
                arrays[name+"_standard_error_"+key] = values
        ring_volume = np.pi * (rb[:, 1]**2-rb[:, 0]**2) * (xb[:, 1]-xb[:, 0])
        arrays["continuum_ring_volume"] = ring_volume
        arrays["DSMC_volume_over_continuum_ring_volume"] = measured["volumes"]/ring_volume
        points = np.column_stack([xx.ravel(), rr.ravel(), np.zeros(xx.size)])
        common = (cai["n"] / source.number_density_m3 >= float(settings.get("density_floor_over_n0", .001)))
        common &= measured["volumes"] > 0
        common &= ~comparison.lip_zone(points, source.radius_m, d)
        arrays["common_density_mask"] = common
        arrays["common_moment_mask"] = common & (measured["n"] > 0)
        np.savez_compressed(target / "raw.npz", **arrays)
        case_metrics = []
        for observable, ref, baseline, vc, scale in (
            ("density", measured["n"], cai["n"], corrected["n"], None),
            ("axial_velocity", measured["U"][:, 0], cai["U"][:, 0], corrected["U"][:, 0], source.thermal_speed*source.speed_ratio),
            ("cylindrical_radial_velocity", measured["U"][:, 1], cai["U"][:, 1], corrected["U"][:, 1], source.thermal_speed),
            ("scalar_temperature", measured["T_trans"], cai["T_trans"], corrected["T_trans"], source.temperature_K)):
            mask = common if observable == "density" else common & (measured["n"] > 0)
            score = dsmc_reference_metrics(ref, baseline, vc, weights=measured["volumes"], common_mask=mask, scale=scale)
            score.update({"path": row["path"], "kn": row["nominal_kn"], "observable": observable,
                          "quality": row["quality"], "headline_eligible": False,
                          "estimator": "same intended cylindrical bins; DSMC cell-centre volume quadrature versus continuum ring transport; convergence not established",
                          "estimator_limitation": "boundary-cell volume fractions and geometric coverage recorded; no headline claims until bin/cell convergence is established",
                          "corrected_invalid_bins_on_common_mask": int(np.sum(mask & ~corrected["valid"])),
                          "corrected_validity": "moment realizability only; invalid bins retained in raw arrays and errors",
                          "dsmc_uncertainty": "unavailable: instantaneous blocks not supplied"})
            uncertainty_key, component = {"density": ("n", None), "axial_velocity": ("U", 0),
                                          "cylindrical_radial_velocity": ("U", 1),
                                          "scalar_temperature": ("T_trans", None)}[observable]
            se = corrected.get("standard_error", {}).get(uncertainty_key)
            if se is not None:
                se = np.asarray(se)[:, component] if component is not None else np.asarray(se)
                score["model_batch_standard_error_rms_physical_units"] = (float(np.sqrt(np.average(se[mask]**2,
                    weights=measured["volumes"][mask]))) if mask.any() and np.isfinite(se[mask]).all() else None)
            score["model_uncertainty_scope"] = "scrambled-batch sampling only; domain, quadrature, estimator bias separate"
            case_metrics.append(score)
        write_json(target / "metrics.json", case_metrics)
        write_json(target / "diagnostics.json", corrected.get("diagnostics", {}))
        write_json(target / "estimator.json", {
            "status": "exploratory-until-spatial-estimator-convergence",
            "DSMC_volume_over_ring_volume": measured["volumes"]/ring_volume,
            "boundary_cell_volume_fraction": measured["boundary_cell_volume_fraction"],
            "boundary_indicator": measured["boundary_indicator_status"],
            "DSMC_tensor": measured["tensor_status"],
            "model_uncertainty": "batch standard errors in raw.npz, separate from unavailable DSMC uncertainty"})
        angular = angular_ring_profiles(points, ring_volume, measured, cai, corrected, source)
        write_rows(target / "angular.csv", angular)
        a = {key: np.array([value[key] for value in angular]) for key in angular[0]}
        common_angular = (a["contributing_ring_bins"] > 0)
        common_angular &= a["Cai_n_over_n0"] >= comparison.ANGULAR_FLOOR_FRACTION * np.nanmax(a["Cai_n_over_n0"], initial=0)
        for observable, suffix, scale in (("angular_density", "n_over_n0", None),
                ("angular_spherical_velocity", "uR_sqrt_beta0", source.speed_ratio),
                ("angular_particle_flux", "F", None),
                ("angular_scalar_temperature", "T_over_T0", 1.)):
            score = dsmc_reference_metrics(a["DSMC_"+suffix], a["Cai_"+suffix],
                    a["Virtual-collision model_"+suffix], weights=a["solid_angle_sr"],
                    common_mask=common_angular, scale=scale)
            score.update({"path": row["path"], "kn": row["nominal_kn"], "observable": observable,
                          "quality": row["quality"], "headline_eligible": False,
                          "estimator": "common coarse ring quadrature; existing 10D shell and 2.5deg footprints",
                          "dsmc_uncertainty": "unavailable: instantaneous blocks not supplied"})
            case_metrics.append(score)
        write_json(target / "metrics.json", case_metrics)
        metrics.extend(case_metrics)
        if plots:
            _plot_reference_bins(target, xx/d, rr/d, measured, cai, corrected, source, row)
            _plot_angular_ring_profiles(target, angular, row)
        processed.append({"path": row["path"], "quality": row["quality"], "headline_eligible": False,
                          "output": str(target), "bins": nx*nr, "tensor_validation": "unavailable"})
        if smoke:
            break
    write_json(output / "reference_skips.json", skipped)
    write_json(output / "reference_assessment.json", {"processed": processed, "skipped": skipped})
    write_rows(output / "assessment.csv", metrics, columns=("path", "kn", "observable", "quality",
        "headline_eligible", "status", "baseline_error", "corrected_error", "error_difference",
        "relative_improvement", "baseline_signed_bias", "corrected_signed_bias"))
    write_rows(output / "assessment_skips.csv", skipped,
               columns=("path", "reason", "available_exports", "plots"))
    if plots and metrics:
        _plot_error_kn(output, metrics)
    return {"processed_cases": len(processed), "primary_cases": sum(r["headline_eligible"] for r in processed),
            "processed": processed, "skipped": skipped,
            "uncertainty": "DSMC unavailable without instantaneous stationary block series",
            "compact_exports": "aggregated tables cannot reconstruct missing spatial weights or tensors"}


def _plot_reference_bins(output, x, r, measured, cai, corrected, source, row):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    shape = x.shape
    fields = [measured["n"], cai["n"], corrected["n"]]
    labels = ["DSMC", "Cai collisionless", "Virtual-collision model"]
    levels = np.array([.001, .002, .005, .01, .02, .05, .1, .2, .5, 1.])
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharex=True, sharey=True)
    for ax, values, label in zip(axes, fields, labels):
        raw = np.asarray(values).reshape(shape)/source.number_density_m3
        image = ax.contourf(x, r, np.ma.masked_invalid(raw), levels=levels, extend="both")
        ax.set_title(label)
        ax.set_xlabel("X/D")
    axes[0].set_ylabel("r/D")
    fig.suptitle(row["path"] + " | " + row["quality"].upper())
    fig.colorbar(image, ax=list(axes), label="n/n0")
    fig.savefig(output / "density_three_way.png", dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharex=True, sharey=True)
    differences = [(values - measured["n"]).reshape(shape)/source.number_density_m3
                   for values in fields[1:]]
    finite = np.concatenate([value[np.isfinite(value)] for value in differences])
    extent = float(np.max(np.abs(finite))) if finite.size else 1.
    extent = max(extent, 1e-15)
    for ax, difference, label in zip(axes, differences, labels[1:]):
        image = ax.pcolormesh(x, r, np.ma.masked_invalid(difference), cmap="coolwarm",
                              shading="nearest", vmin=-extent, vmax=extent)
        ax.set_title(label + " minus DSMC")
        ax.set_xlabel("X/D")
    axes[0].set_ylabel("r/D")
    fig.colorbar(image, ax=list(axes), label="signed density error / n0")
    fig.suptitle(row["quality"].upper() + " | spatial estimator convergence required")
    fig.savefig(output / "density_errors.png", dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, key, title, scale in zip(axes, ("n", "U", "T_trans"),
                                    ("Density n/n0", "Axial velocity U sqrt(beta0)", "Scalar T/T0"),
                                    (source.number_density_m3, source.thermal_speed, source.temperature_K)):
        for data, label in zip((measured, cai, corrected), labels):
            values = data[key][:, 0] if key == "U" else data[key]
            ax.plot(x[:, 0], np.asarray(values).reshape(shape)[:, 0]/scale, label=label)
        ax.set_xlabel("X/D")
        ax.set_title(title)
    axes[0].legend()
    fig.suptitle("Same finite-radius central bin | " + row["quality"].upper())
    fig.tight_layout()
    fig.savefig(output / "centerline_three_way.png", dpi=150)
    plt.close(fig)
    discrepancy = {label: contour_location_discrepancy(x, r,
                    measured["n"].reshape(shape)/source.number_density_m3,
                    values.reshape(shape)/source.number_density_m3, levels[:-1])
                   for label, values in zip(labels[1:], fields[1:])}
    write_json(output / "contour_location.json", {"units": "D", "metrics": discrepancy})


def angular_ring_profiles(points, bin_volumes, measured, cai, corrected, source):
    """Existing 10-D shell/2.5-degree conventions applied to common ring bins.

    This is a coarse ring quadrature of the established shell estimator, with
    the same footprint fractions and weights for every prediction. It does not
    claim equivalence to legacy angular CSVs, whose original cells are absent.
    """
    from plumetools.cai2012.constants import BOLTZMANN_J_PER_K
    from plumetools.cai2012ext import comparison
    points = np.asarray(points)
    d = 2*source.radius_m
    sphere = np.linalg.norm(points, axis=1)
    theta = np.arctan2(points[:, 1], points[:, 0])
    shell = np.abs(sphere - 10*d) <= .5*d + .5*np.cbrt(bin_volumes)
    edges = np.radians(np.arange(0, 90.01, 2.5))
    fractions = comparison.footprint_fractions(theta[shell], sphere[shell],
                                              np.asarray(bin_volumes)[shell], edges)
    shared = fractions * np.asarray(bin_volumes)[shell, None]
    volume = shared.sum(axis=0)
    def avg(values):
        values = np.asarray(values)[shell]
        # A missing contributing bin makes that angular estimate unavailable.
        missing = np.any((shared > 0) & ~np.isfinite(values[:, None]), axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            out = np.sum(shared * np.where(np.isfinite(values), values, 0)[:, None], axis=0)/volume
        return np.where(missing, np.nan, out)
    out = {"theta_deg": np.degrees(.5*(edges[1:]+edges[:-1])),
           "solid_angle_sr": 2*np.pi*(np.cos(edges[:-1])-np.cos(edges[1:])),
           "contributing_ring_bins": fractions.sum(axis=0)}
    ct, st = np.cos(theta), np.sin(theta)
    for label, data in (("DSMC", measured), ("Cai", cai), ("Virtual-collision model", corrected)):
        n = data["n"]
        m1 = data["M1"]
        radial_m1 = ct*m1[:, 0]+st*m1[:, 1]
        polar_m1 = -st*m1[:, 0]+ct*m1[:, 1]
        trace = data["M2_trace"] if "M2_trace" in data else np.trace(data["M2"], axis1=1, axis2=2)
        density = avg(n)
        with np.errstate(invalid="ignore", divide="ignore"):
            ur, ut, up = avg(radial_m1)/density, avg(polar_m1)/density, avg(m1[:, 2])/density
            temperature = (avg(trace)/density-ur**2-ut**2-up**2)*source.species.mass_kg/(3*BOLTZMANN_J_PER_K)
        flux = comparison.mass_flux_per_steradian(n/source.number_density_m3,
                    np.divide(radial_m1, n, out=np.full_like(n, np.nan), where=n != 0)/source.thermal_speed,
                    sphere, source.radius_m, source.speed_ratio)
        out.update({label+"_n_over_n0": density/source.number_density_m3,
                    label+"_uR_sqrt_beta0": ur/source.thermal_speed,
                    label+"_utheta_sqrt_beta0": ut/source.thermal_speed,
                    label+"_T_over_T0": temperature/source.temperature_K,
                    label+"_F": avg(flux)})
    return [{key: np.asarray(values)[i].item() for key, values in out.items()}
            for i in range(len(edges)-1)]


def _plot_angular_ring_profiles(output, rows, case):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, suffix, title in zip(axes, ("n_over_n0", "uR_sqrt_beta0", "F"),
                                  ("Angular density n/n0", "Spherical outward velocity", "Normalized particle flux / sr")):
        for label in ("DSMC", "Cai", "Virtual-collision model"):
            ax.plot([row["theta_deg"] for row in rows], [row[label+"_"+suffix] for row in rows],
                    label="Cai collisionless" if label == "Cai" else label)
        ax.set_xlabel("theta [degrees]")
        ax.set_title(title)
    axes[0].legend()
    fig.suptitle(case["quality"].upper()+" | common coarse rings in 10-D shell; estimator convergence required")
    fig.tight_layout()
    fig.savefig(output / "angular_three_way.png", dpi=150)
    plt.close(fig)


def _plot_error_kn(output, metrics):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, observable in zip(axes, ("density", "axial_velocity", "scalar_temperature")):
        rows = sorted([row for row in metrics if row["observable"] == observable], key=lambda row: row["kn"])
        for key, label in (("baseline_error", "Cai collisionless"), ("corrected_error", "Virtual-collision model")):
            valid = [row for row in rows if row[key] is not None]
            ax.scatter([row["kn"] for row in valid], [row[key] for row in valid], label=label)
        ax.set_xscale("log")
        ax.set_xlabel("Nozzle-diameter Kn (stratify by source/gas/geometry)")
        ax.set_title(observable.replace("_", " "))
    axes[0].legend()
    fig.suptitle("Exploratory numerical assessment | DSMC confidence intervals unavailable")
    fig.tight_layout()
    fig.savefig(output / "error_versus_kn.png", dpi=150)
    plt.close(fig)
