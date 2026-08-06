r"""DSMC mesh and time-step quality checks.

Most of these **warn**. A cell coarser than the local mean free path degrades the
collision statistics, but it is a judgement about accuracy, not a broken case, and
silently redesigning someone's mesh is worse than telling them about it.

Five conditions are **errors**, because a case with any of them cannot produce a
meaningful answer at all:

* non-positive or non-finite density or temperature;
* a NaN anywhere in the generated inflow fields;
* the inflow surface intersecting the cylinder;
* a time step that lets a particle cross more than one of the smallest cells;
* a required patch missing from the mesh.

Every result is machine-readable, so a case summary can carry the whole set
rather than a boolean.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from plumetools.markelov1999 import sourceflow as sf
from plumetools.markelov1999.geometry import MarkelovGeometry

#: Patch roles every generated mesh must have. Checked by name, from the config,
#: so renaming a patch in ``case.yaml`` moves the requirement with it.
REQUIRED_PATCH_ROLES = ("inflow", "cylinder", "plate", "outer", "symmetry",
                        "upstream_vacuum")


@dataclass(frozen=True)
class CheckResult:
    """One check outcome.

    Attributes:
        name: stable identifier, used as the key in a case summary.
        status: ``"pass"``, ``"warn"`` or ``"fail"``.
        message: one line for a human.
        value: the measured number, or ``None`` for a structural check.
        threshold: what it was compared against, or ``None``.
    """

    name: str
    status: str
    message: str
    value: float | None = None
    threshold: float | None = None

    @property
    def failed(self) -> bool:
        return self.status == "fail"

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "status": self.status,
            "message": self.message,
            "value": self.value,
            "threshold": self.threshold,
        }


class DsmcCheckError(ValueError):
    """One or more checks failed outright."""


@dataclass
class CheckReport:
    """The full set of check results."""

    results: list = field(default_factory=list)

    def add(self, result: CheckResult) -> None:
        self.results.append(result)

    @property
    def failures(self) -> list:
        return [r for r in self.results if r.status == "fail"]

    @property
    def warnings(self) -> list:
        return [r for r in self.results if r.status == "warn"]

    def raise_if_failed(self) -> None:
        """Raise :class:`DsmcCheckError` listing every failure, not just the first.

        All of them, because they usually share a cause -- a time step that is too
        large is often the same mistake as a refinement level that is too high --
        and fixing them one run at a time is slow.
        """
        if self.failures:
            lines = [f"{len(self.failures)} DSMC check(s) failed:"]
            lines += [f"  FAIL {r.name}: {r.message}" for r in self.failures]
            raise DsmcCheckError("\n".join(lines))

    def report(self) -> list[str]:
        symbols = {"pass": "  OK  ", "warn": "  WARN", "fail": "  FAIL"}
        lines = ["DSMC quality checks"]
        lines += [f"{symbols[r.status]} {r.name}: {r.message}" for r in self.results]
        if self.warnings:
            lines.append(f"  {len(self.warnings)} warning(s); "
                         f"see docs/markelov1999-case.md for what each one costs.")
        return lines

    def as_dict(self) -> dict:
        return {
            "passed": not self.failures,
            "n_failures": len(self.failures),
            "n_warnings": len(self.warnings),
            "results": [r.as_dict() for r in self.results],
        }


def run_checks(cfg, geom: MarkelovGeometry, *, inflow=None, patches=None,
               estimate=None, cell_sizes=None) -> CheckReport:
    """Run every configured check.

    Args:
        cfg: the case config.
        geom: the geometry.
        inflow: an ``InflowResult``, if the mesh exists. Skipped if ``None``.
        patches: the mesh's patches. Skipped if ``None``.
        estimate: a :class:`plumetools.markelov1999.resolution.ResolutionEstimate`.
        cell_sizes: cell size at each refined surface [m].

    Returns:
        A :class:`CheckReport`. Nothing is raised here -- call
        :meth:`CheckReport.raise_if_failed` to stop on an error, so a caller can
        record the full report first.
    """
    report = CheckReport()

    if not cfg.checks.enabled:
        report.add(CheckResult(
            "checks_enabled", "warn",
            "checks.enabled is false, so no DSMC quality check was run"))
        return report

    _check_positive_inputs(cfg, report)
    _check_inflow_clearance(geom, report)
    if patches is not None:
        _check_patches(cfg, patches, report)
    if cell_sizes is not None:
        _check_time_step(cfg, cell_sizes, report)
    if estimate is not None:
        _check_mean_free_path(cfg, estimate, report)
        _check_occupancy(cfg, estimate, report)
    if inflow is not None:
        _check_inflow_fields(cfg, inflow, report)

    return report


def _check_positive_inputs(cfg, report: CheckReport) -> None:
    """Density and temperature inputs must be finite and positive. HARD ERROR."""
    checks = {
        "stagnation.p0_pa": cfg.stagnation.p0_pa,
        "stagnation.T0_K": cfg.stagnation.T0_K,
        "stagnation.throat_radius_m": cfg.stagnation.throat_radius_m,
        "dsmc.wall_temperature_K": cfg.dsmc.wall_temperature_K,
        "dsmc.delta_t_s": cfg.dsmc.delta_t_s,
        "dsmc.end_time_s": cfg.dsmc.end_time_s,
    }
    bad = [name for name, value in checks.items()
           if value is None or not math.isfinite(value) or value <= 0.0]
    if bad:
        report.add(CheckResult(
            "positive_inputs", "fail",
            f"{bad} must be finite and positive; a non-positive pressure or "
            f"temperature has no physical meaning and produces a zero or negative "
            f"density"))
    else:
        report.add(CheckResult(
            "positive_inputs", "pass",
            f"pressure {cfg.stagnation.p0_pa:.6g} Pa, temperature "
            f"{cfg.stagnation.T0_K:g} K, time step {cfg.dsmc.delta_t_s:g} s "
            f"all finite and positive"))


def _check_inflow_clearance(geom: MarkelovGeometry, report: CheckReport) -> None:
    """The inflow surface must not reach the cylinder. HARD ERROR.

    Also reports the minimum distance, which requirement 5 asks for explicitly.
    """
    clearance = geom.inflow_to_cylinder_clearance_m()
    if clearance <= 0.0:
        report.add(CheckResult(
            "inflow_cylinder_clearance", "fail",
            f"the inflow hemisphere (R = {geom.inflow_radius_m:.6g} m) intersects "
            f"the cylinder: clearance {clearance:.6g} m. The source-flow solution "
            f"would be imposed inside a solid body.",
            value=clearance, threshold=0.0))
    else:
        report.add(CheckResult(
            "inflow_cylinder_clearance", "pass",
            f"minimum distance from the inflow surface to the cylinder is "
            f"{clearance:.6f} m ({clearance / 0.0254:.2f} in)",
            value=clearance, threshold=0.0))


def _check_patches(cfg, patches, report: CheckReport) -> None:
    """Required patches present, no duplicates, no empty physical patches."""
    names = cfg.mesh.patch_names
    expected = {role: names[role] for role in REQUIRED_PATCH_ROLES if role in names}

    missing = sorted(name for name in expected.values() if name not in patches)
    if missing:
        report.add(CheckResult(
            "required_patches", "fail",
            f"the mesh is missing patch(es) {missing}; it has {sorted(patches)}"))
    else:
        report.add(CheckResult(
            "required_patches", "pass",
            f"all {len(expected)} required patches present: "
            f"{sorted(expected.values())}"))

    # Duplicate NAMES cannot occur -- read_boundary returns a dict -- but a
    # duplicated role mapping can, and it silently merges two boundaries with
    # different physics.
    used = list(expected.values())
    duplicated = sorted({n for n in used if used.count(n) > 1})
    if duplicated:
        report.add(CheckResult(
            "duplicate_patches", "fail",
            f"patch name(s) {duplicated} serve more than one role"))
    else:
        report.add(CheckResult(
            "duplicate_patches", "pass",
            "no patch name serves two roles"))

    empty = sorted(name for name, info in patches.items()
                   if info.n_faces == 0 and name in expected.values())
    if empty:
        report.add(CheckResult(
            "empty_patches", "fail",
            f"physical patch(es) {empty} have zero faces. snappyHexMesh creates a "
            f"patch even when it carved nothing, so an empty body patch means the "
            f"body is not in the mesh at all -- the run would proceed with no "
            f"obstruction where the geometry says there is one."))
    else:
        report.add(CheckResult(
            "empty_patches", "pass",
            "every required patch has faces: "
            + ", ".join(f"{n}({patches[n].n_faces})" for n in sorted(expected.values()))))


def _check_time_step(cfg, cell_sizes: dict, report: CheckReport) -> None:
    """A particle must not cross more than one of the smallest cells. HARD ERROR.

    The characteristic speed is the limiting speed plus three thermal standard
    deviations, which bounds essentially the whole distribution rather than its
    mean -- the mean would let the fast tail cross several cells unnoticed.
    """
    from plumetools.markelov1999.inflow import molecular_mass_from_config

    mass = molecular_mass_from_config(cfg)
    v_limit = sf.limiting_velocity(cfg.gas.gamma, cfg.stagnation.T0_K, mass)
    thermal_sigma = math.sqrt(
        1.380649e-23 * cfg.stagnation.T0_K / mass)
    speed = v_limit + 3.0 * thermal_sigma

    smallest = min(cell_sizes.values())
    courant = speed * cfg.dsmc.delta_t_s / smallest

    if courant > cfg.checks.max_courant:
        report.add(CheckResult(
            "time_step_courant", "fail",
            f"a particle at {speed:.0f} m/s crosses {courant:.2f} of the smallest "
            f"cell ({smallest:.6g} m) per step of {cfg.dsmc.delta_t_s:g} s. Above "
            f"{cfg.checks.max_courant:g} the collision sampling is invalid: "
            f"particles skip cells without being offered a collision partner in "
            f"them. Reduce dsmc.delta_t_s below "
            f"{cfg.checks.max_courant * smallest / speed:.3g} s.",
            value=courant, threshold=cfg.checks.max_courant))
    elif courant > cfg.checks.warn_courant:
        report.add(CheckResult(
            "time_step_courant", "warn",
            f"a particle crosses {courant:.2f} of the smallest cell per step; "
            f"below {cfg.checks.warn_courant:g} is preferable",
            value=courant, threshold=cfg.checks.warn_courant))
    else:
        report.add(CheckResult(
            "time_step_courant", "pass",
            f"a particle at {speed:.0f} m/s crosses {courant:.3f} of the smallest "
            f"cell ({smallest:.6g} m) per step",
            value=courant, threshold=cfg.checks.warn_courant))


def _check_mean_free_path(cfg, estimate, report: CheckReport) -> None:
    """Cell size against the local mean free path. WARNING.

    DSMC theory wants a collision cell well inside one mean free path so that
    collision partners are genuinely local. Exceeding it biases the collision rate
    rather than breaking the run, so this warns and reports by how much.
    """
    worst = max(estimate.regions, key=lambda r: r.cell_over_mfp)
    ratio = worst.cell_over_mfp

    if ratio > cfg.checks.max_cell_over_mfp:
        report.add(CheckResult(
            "cell_over_mean_free_path", "warn",
            f"the '{worst.region.name}' cell is {ratio:.2f} local mean free paths "
            f"({worst.region.cell_size_m:.6g} m against "
            f"{worst.mean_free_path_m:.6g} m). Above "
            f"{cfg.checks.max_cell_over_mfp:g} the collision rate is biased; one "
            f"more refinement level there would fix the ratio and multiply the "
            f"particle count by eight.",
            value=ratio, threshold=cfg.checks.max_cell_over_mfp))
    else:
        report.add(CheckResult(
            "cell_over_mean_free_path", "pass",
            f"the worst cell-to-mean-free-path ratio is {ratio:.3f}, at the "
            f"'{worst.region.name}' region",
            value=ratio, threshold=cfg.checks.max_cell_over_mfp))


def _check_occupancy(cfg, estimate, report: CheckReport) -> None:
    """Estimated particles per cell in the regions of interest. WARNING.

    An estimate, so never an error: only the post-run audit can say what the
    occupancy was.
    """
    worst = min(estimate.regions, key=lambda r: r.particles_per_cell)
    value = worst.particles_per_cell

    if value < cfg.checks.min_particles_per_cell:
        report.add(CheckResult(
            "particles_per_cell", "warn",
            f"the '{worst.region.name}' region is estimated at {value:.1f} "
            f"particles/cell, below the {cfg.checks.min_particles_per_cell:g} "
            f"floor. Standard dsmcFoam has one weight for the whole domain, so "
            f"meeting the target here means raising it everywhere; set "
            f"resolution.sizing_region to '{worst.region.name}' if that is the "
            f"trade you want.",
            value=value, threshold=cfg.checks.min_particles_per_cell))
    elif value < estimate.target_particles_per_cell:
        report.add(CheckResult(
            "particles_per_cell", "warn",
            f"the '{worst.region.name}' region is estimated at {value:.1f} "
            f"particles/cell, under the "
            f"{estimate.target_particles_per_cell:g} target but above the floor",
            value=value, threshold=estimate.target_particles_per_cell))
    else:
        report.add(CheckResult(
            "particles_per_cell", "pass",
            f"every region of interest is estimated at or above the "
            f"{estimate.target_particles_per_cell:g}-particle target "
            f"(worst {value:.1f}, at '{worst.region.name}')",
            value=value, threshold=estimate.target_particles_per_cell))


def _check_inflow_fields(cfg, inflow, report: CheckReport) -> None:
    """Generated inflow values finite and positive, and inside the plume cone."""
    arrays = {
        "rhoN": inflow.rhoN,
        "T": inflow.T,
        "U": inflow.U,
        "areas": inflow.areas,
    }
    non_finite = {name: int(np.count_nonzero(~np.isfinite(a)))
                  for name, a in arrays.items()}
    if any(non_finite.values()):
        report.add(CheckResult(
            "inflow_finite", "fail",
            f"non-finite values in the generated inflow fields: "
            f"{ {k: v for k, v in non_finite.items() if v} }"))
    else:
        report.add(CheckResult(
            "inflow_finite", "pass",
            f"all {inflow.n_faces} inflow faces carry finite values"))

    non_positive = int(np.count_nonzero(inflow.rhoN <= 0.0))
    if non_positive:
        report.add(CheckResult(
            "inflow_density_positive", "fail",
            f"{non_positive} of {inflow.n_faces} inflow faces have a non-positive "
            f"number density. On a hemispherical source with a diatomic gas the "
            f"whole surface is inside the {math.degrees(sf.limiting_angle(cfg.gas.gamma)):.1f} "
            f"deg plume cone, so this means the inflow patch is not the hemisphere "
            f"the configuration describes.",
            value=float(non_positive)))
    else:
        report.add(CheckResult(
            "inflow_density_positive", "pass",
            f"number density {inflow.rhoN.min():.4e} .. {inflow.rhoN.max():.4e} 1/m^3, "
            f"all positive"))

    if np.any(inflow.T <= 0.0):
        report.add(CheckResult(
            "inflow_temperature_positive", "fail",
            "non-positive inflow temperature; the solver aborts on this with "
            "'Zero boundary temperature detected'"))
    else:
        report.add(CheckResult(
            "inflow_temperature_positive", "pass",
            f"inflow temperature {inflow.T.min():g} .. {inflow.T.max():g} K"))

    theta_l = sf.limiting_angle(cfg.gas.gamma)
    outside = int(np.count_nonzero(inflow.theta >= theta_l))
    if outside:
        report.add(CheckResult(
            "inflow_inside_plume_cone", "fail",
            f"{outside} of {inflow.n_faces} inflow faces lie at or beyond the "
            f"limiting angle {math.degrees(theta_l):.2f} deg, where the density is "
            f"exactly zero. They would inject nothing.",
            value=float(outside)))
    else:
        report.add(CheckResult(
            "inflow_inside_plume_cone", "pass",
            f"every inflow face is inside the {math.degrees(theta_l):.2f} deg plume "
            f"cone (worst {math.degrees(inflow.theta.max()):.2f} deg)",
            value=float(math.degrees(inflow.theta.max())),
            threshold=float(math.degrees(theta_l))))


def write_summary(case_dir: Path, summary: dict) -> Path:
    """Write the machine-readable case summary.

    Args:
        case_dir: the case directory.
        summary: anything JSON-serialisable.

    Returns:
        The path written, ``<case>/case-summary.json``.

    JSON rather than YAML: it is what a downstream comparison script will read,
    and there is no need for comments in a generated record.
    """
    import json

    path = Path(case_dir) / "case-summary.json"
    path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
    return path
