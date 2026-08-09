r"""DSMC mesh, time-step and boundary-type quality checks.

Most of these **warn**. A cell coarser than the local mean free path biases the
collision rate, but that is a judgement about accuracy, not a broken case, and
silently redesigning someone's mesh is worse than telling them about it.

These are **errors**, because a case with any of them cannot produce a
meaningful answer at all:

* a non-positive or non-finite density, temperature, speed or time step;
* a time step that lets a molecule cross more than one of the smallest cells;
* a missing required patch, or a nozzle patch with no faces -- nothing would
  ever be injected;
* a nozzle patch whose area misses :math:`\pi R_0^2` by more than
  ``checks.max_nozzle_area_error`` -- the injected mass flow is proportional to
  it, so the whole solution scales with the error;
* a vacuum boundary declared ``wall`` -- the plume would expand into a
  reflecting box rather than into vacuum, which is a different problem;
* the nozzle declared ``wall`` -- ``plumeFieldInflow`` refuses to inject through
  one, so this fails anyway, but much later and much less clearly.

Every result is machine-readable, so a case summary carries the whole set rather
than a boolean.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path


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
        return {"name": self.name, "status": self.status, "message": self.message,
                "value": self.value, "threshold": self.threshold}


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
        """Raise :class:`DsmcCheckError` listing **every** failure.

        All of them, not the first: they usually share a cause -- an over-large
        time step and an over-refined mesh are often the same mistake -- and
        fixing them one run at a time is slow.
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
            lines.append(f"  {len(self.warnings)} warning(s); see "
                         f"docs/cai2012-case.md for what each one costs.")
        return lines

    def as_dict(self) -> dict:
        return {
            "passed": not self.failures,
            "n_failures": len(self.failures),
            "n_warnings": len(self.warnings),
            "results": [r.as_dict() for r in self.results],
        }


#: Patch roles every generated mesh must have.
REQUIRED_PATCH_ROLES = ("nozzle", "upstream_vacuum", "outer")


def run_checks(cfg, geom, exit_state, plan=None, run=None, *,
               patches=None, nozzle=None, particles=None) -> CheckReport:
    """Run every configured check.

    Args:
        cfg: the case config.
        geom: the :class:`~plumetools.cai2012.geometry.CaiGeometry`.
        exit_state: the :class:`~plumetools.cai2012.inflow.ExitState`.
        plan: the :class:`~plumetools.cai2012.mesh.MeshPlan`, if there is one.
        run: the :class:`~plumetools.cai2012.inflow.RunSettings`, if there is one.
        patches: the mesh's patches, once it exists.
        nozzle: a :class:`~plumetools.cai2012.mesh.NozzlePrediction`, or the
            measured equivalent.
        particles: the dict :func:`plumetools.cai2012.mesh.estimate_particles`
            returns.

    Returns:
        A :class:`CheckReport`. Nothing is raised here -- call
        :meth:`CheckReport.raise_if_failed` so a caller can record the full
        report first.
    """
    report = CheckReport()

    if not cfg.checks.enabled:
        report.add(CheckResult(
            "checks_enabled", "warn",
            "checks.enabled is false, so no DSMC quality check was run"))
        return report

    _check_exit_state(exit_state, report)
    if plan is not None:
        _check_cell_over_mfp(cfg, plan, report)
        _check_cell_budget(cfg, plan, report)
        _check_grading(plan, report)
    if plan is not None and run is not None:
        _check_time_step(cfg, plan, run, exit_state, report)
        _check_occupancy(cfg, run, report)
    if nozzle is not None:
        _check_nozzle(cfg, nozzle, report)
    if patches is not None:
        _check_patches(cfg, patches, report)
    if particles is not None:
        _check_particle_count(particles, report)
    return report


def _check_exit_state(exit_state, report: CheckReport) -> None:
    """Every derived exit quantity finite and positive. HARD ERROR."""
    values = {
        "n0": exit_state.number_density_per_m3,
        "T0": exit_state.T0_K,
        "lambda0": exit_state.mean_free_path_m,
        "Kn": exit_state.knudsen,
        "number flux": exit_state.number_flux_per_m2_s,
    }
    bad = sorted(name for name, v in values.items()
                 if v is None or not math.isfinite(float(v)) or float(v) <= 0.0)
    if bad:
        report.add(CheckResult(
            "exit_state_positive", "fail",
            f"{bad} must be finite and positive. The exit state is derived from "
            f"Kn, S0 and T0, so a bad value here means one of those three is "
            f"wrong in case.yaml."))
    else:
        report.add(CheckResult(
            "exit_state_positive", "pass",
            f"Kn = {exit_state.knudsen:g} -> lambda0 = "
            f"{exit_state.mean_free_path_m:.4e} m -> n0 = "
            f"{exit_state.number_density_per_m3:.4e} 1/m^3, "
            f"U0 = {exit_state.velocity_m_per_s:.1f} m/s, T0 = {exit_state.T0_K:g} K"))


def _check_cell_over_mfp(cfg, plan, report: CheckReport) -> None:
    """Core cell against the exit mean free path. WARNING.

    Cai's own criterion is ``dx / lambda0 = 1``. Exceeding it biases the
    collision rate rather than breaking the run, and at ``Kn = 100`` the mean
    free path is a hundred domain diameters so the ratio is meaningless -- it
    warns and reports by how much, and the documentation carries the number for
    every case.
    """
    ratio = plan.cell_over_mfp
    limit = float(cfg.checks.max_cell_over_mfp)
    if ratio > limit:
        report.add(CheckResult(
            "cell_over_mean_free_path", "warn",
            f"the core cell is {ratio:.2f} exit mean free paths "
            f"({plan.core_cell_size_m:.4e} m against "
            f"{plan.mean_free_path_m:.4e} m). Cai's criterion is 1. Halving the "
            f"cell fixes the ratio and multiplies the cell count by eight; the "
            f"budget is mesh.max_cells = {plan.max_cells:,}.",
            value=ratio, threshold=limit))
    else:
        report.add(CheckResult(
            "cell_over_mean_free_path", "pass",
            f"the core cell is {ratio:.3g} exit mean free paths, within Cai's "
            f"criterion of {limit:g}",
            value=ratio, threshold=limit))


def _check_cell_budget(cfg, plan, report: CheckReport) -> None:
    """Whether the mesh had to be coarsened to fit the budget. WARNING."""
    if plan.coarsened:
        report.add(CheckResult(
            "cell_budget", "warn",
            f"the DSMC criterion asked for a {plan.requested_cell_size_m:.4e} m "
            f"core cell; that did not fit in mesh.max_cells = "
            f"{plan.max_cells:,}, so {plan.core_cell_size_m:.4e} m was used "
            f"({plan.n_cells:,} cells). This is a resolution compromise and is "
            f"recorded in the manifest.",
            value=float(plan.n_cells), threshold=float(plan.max_cells)))
    else:
        report.add(CheckResult(
            "cell_budget", "pass",
            f"{plan.n_cells:,} cells at the requested "
            f"{plan.core_cell_size_m:.4e} m core cell, within the "
            f"{plan.max_cells:,} budget",
            value=float(plan.n_cells), threshold=float(plan.max_cells)))


def _check_grading(plan, report: CheckReport) -> None:
    """Continuity of cell size across the core boundary. WARNING."""
    jump = plan.max_grading_jump
    if jump > 1.001:
        report.add(CheckResult(
            "grading_continuity", "warn",
            f"the first cell of an expanding segment is {jump:.2f}x the core "
            f"cell, because mesh.outer_expansion capped the ratio. A step in "
            f"cell size shows up in the sampled density as a feature that looks "
            f"physical. Raise mesh.outer_expansion or add outer cells.",
            value=jump, threshold=1.0))
    else:
        report.add(CheckResult(
            "grading_continuity", "pass",
            f"cell size is continuous across the core boundary "
            f"(largest jump {jump:.4f}x); the expansion ratios were solved for, "
            f"not assumed",
            value=jump, threshold=1.0))


def _check_time_step(cfg, plan, run, exit_state, report: CheckReport) -> None:
    """A molecule must not cross more than one of the smallest cells. HARD ERROR.

    The characteristic speed is ``U0 + 3 sigma``, which bounds essentially the
    whole distribution. The mean would let the fast tail skip cells unnoticed,
    and it is the tail that breaks the collision sampling.
    """
    courant = run.courant
    speed = exit_state.characteristic_speed_m_per_s
    smallest = plan.min_cell_size_m

    if courant > float(cfg.checks.max_courant):
        report.add(CheckResult(
            "time_step_courant", "fail",
            f"a molecule at {speed:.0f} m/s crosses {courant:.2f} of the "
            f"smallest cell ({smallest:.4e} m) per step of {run.delta_t_s:g} s. "
            f"Above {cfg.checks.max_courant:g} the collision sampling is "
            f"invalid: molecules skip cells without being offered a partner in "
            f"them. Reduce dsmc.courant_target, or set dsmc.delta_t_s below "
            f"{cfg.checks.max_courant * smallest / speed:.3g} s.",
            value=courant, threshold=float(cfg.checks.max_courant)))
    elif courant > float(cfg.checks.warn_courant):
        report.add(CheckResult(
            "time_step_courant", "warn",
            f"a molecule crosses {courant:.2f} of the smallest cell per step; "
            f"below {cfg.checks.warn_courant:g} is preferable",
            value=courant, threshold=float(cfg.checks.warn_courant)))
    else:
        report.add(CheckResult(
            "time_step_courant", "pass",
            f"a molecule at {speed:.0f} m/s crosses {courant:.3f} of the "
            f"smallest cell ({smallest:.4e} m) per step",
            value=courant, threshold=float(cfg.checks.warn_courant)))

    # Cai's own step, reported whether or not it was used.
    if run.cai_courant > float(cfg.checks.max_courant):
        report.add(CheckResult(
            "cai_time_step", "warn",
            f"Cai's dt/t0 = 1 would be {run.cai_delta_t_s:.4e} s, i.e. Courant "
            f"{run.cai_courant:.2f} on this mesh -- above the "
            f"{cfg.checks.max_courant:g} limit. Cai does not define t0 precisely "
            f"enough to reconstruct his step, so this case derives one from "
            f"dsmc.courant_target instead and reports the difference rather than "
            f"claiming equivalence.",
            value=run.cai_courant, threshold=float(cfg.checks.max_courant)))
    else:
        report.add(CheckResult(
            "cai_time_step", "pass",
            f"Cai's dt/t0 = 1 ({run.cai_delta_t_s:.4e} s) would give Courant "
            f"{run.cai_courant:.3f} here, within the limit",
            value=run.cai_courant, threshold=float(cfg.checks.max_courant)))


def _check_occupancy(cfg, run, report: CheckReport) -> None:
    """Estimated particles in the exit cell. WARNING.

    An estimate, so never an error: only a post-run audit of the sampled
    ``dsmcRhoN`` field can say what the occupancy actually was.
    """
    value = run.exit_particles_per_cell
    floor = float(cfg.checks.min_particles_per_cell)
    target = float(cfg.resolution.target_particles_per_cell)
    if value < floor:
        report.add(CheckResult(
            "particles_per_cell", "warn",
            f"the exit cell is estimated at {value:.1f} particles/cell, below "
            f"the {floor:g} floor. This is the DENSEST cell in the domain, so "
            f"everywhere else is worse. dsmc.n_equivalent_particles is pinned in "
            f"case.yaml; setting it to null derives a weight that meets the "
            f"{target:g} target.",
            value=value, threshold=floor))
    elif value < target:
        report.add(CheckResult(
            "particles_per_cell", "warn",
            f"the exit cell is estimated at {value:.1f} particles/cell, under "
            f"the {target:g} target but above the floor",
            value=value, threshold=target))
    else:
        report.add(CheckResult(
            "particles_per_cell", "pass",
            f"the exit cell is estimated at {value:.1f} particles/cell "
            f"(target {target:g}). Standard dsmcFoam has one weight for the "
            f"whole domain, so the target holds here and falls off downstream "
            f"with the density.",
            value=value, threshold=target))


def _check_nozzle(cfg, nozzle, report: CheckReport) -> None:
    """Face count and area of the staircase nozzle. HARD ERROR on both."""
    if nozzle.n_faces < int(cfg.checks.min_nozzle_faces):
        report.add(CheckResult(
            "nozzle_faces", "fail",
            f"the nozzle patch has {nozzle.n_faces} face(s), below the "
            f"{cfg.checks.min_nozzle_faces} minimum. A disk resolved by a "
            f"handful of squares is not a circular exit -- the injected profile "
            f"would be square. Reduce mesh.core_cell_size_m, or raise "
            f"mesh.max_cells so the derived cell can shrink.",
            value=float(nozzle.n_faces),
            threshold=float(cfg.checks.min_nozzle_faces)))
    else:
        report.add(CheckResult(
            "nozzle_faces", "pass",
            f"the nozzle patch has {nozzle.n_faces} faces, "
            f"{nozzle.cells_across_diameter:.1f} cells across the diameter",
            value=float(nozzle.n_faces),
            threshold=float(cfg.checks.min_nozzle_faces)))

    error = abs(nozzle.area_error)
    limit = float(cfg.checks.max_nozzle_area_error)
    if error > limit:
        report.add(CheckResult(
            "nozzle_area", "fail",
            f"the meshed nozzle area is {nozzle.area_m2:.6e} m^2 against "
            f"pi R0^2 = {nozzle.exact_area_m2:.6e} m^2, an error of "
            f"{100.0 * nozzle.area_error:+.2f}% -- past the "
            f"{100.0 * limit:.1f}% limit. The injected molecule flow is "
            f"proportional to the area, so the whole solution scales with this. "
            f"A Cartesian grid cannot represent a circle exactly; the error "
            f"falls as the core cell shrinks.",
            value=error, threshold=limit))
    else:
        report.add(CheckResult(
            "nozzle_area", "pass",
            f"the meshed nozzle area is {nozzle.area_m2:.6e} m^2, "
            f"{100.0 * nozzle.area_error:+.2f}% from pi R0^2 -- within the "
            f"{100.0 * limit:.1f}% limit for a staircase disk",
            value=error, threshold=limit))


def _check_patches(cfg, patches, report: CheckReport) -> None:
    """Required patches present, present with faces, and of the right type."""
    names = cfg.mesh.patch_names
    expected = {role: names[role] for role in REQUIRED_PATCH_ROLES if role in names}

    missing = sorted(name for name in expected.values() if name not in patches)
    if missing:
        report.add(CheckResult(
            "required_patches", "fail",
            f"the mesh is missing patch(es) {missing}; it has {sorted(patches)}. "
            f"'{names.get('nozzle')}' is created by topoSet + createPatch, not "
            f"by blockMesh."))
        return
    report.add(CheckResult(
        "required_patches", "pass",
        f"all {len(expected)} required patches present: "
        + ", ".join(f"{n}({patches[n].n_faces})"
                    for n in sorted(expected.values()))))

    nozzle = names["nozzle"]
    if patches[nozzle].n_faces == 0:
        report.add(CheckResult(
            "nozzle_not_empty", "fail",
            f"the nozzle patch {nozzle!r} has zero faces, so nothing would ever "
            f"be injected and the domain would stay a vacuum."))
    else:
        report.add(CheckResult(
            "nozzle_not_empty", "pass",
            f"the nozzle patch has {patches[nozzle].n_faces} faces"))

    # Geometric types. These decide the physics, not the labelling.
    walls = sorted(name for name, info in patches.items() if info.type == "wall")
    if walls:
        report.add(CheckResult(
            "no_wall_boundaries", "fail",
            f"patch(es) {walls} have geometric type 'wall'. Every boundary here "
            f"is either the nozzle or open vacuum: a wall REFLECTS molecules "
            f"(DSMCParcel::hitWallPatch runs the WallInteractionModel), so the "
            f"plume would expand into a closed box instead of into vacuum."))
    else:
        report.add(CheckResult(
            "no_wall_boundaries", "pass",
            "no boundary is a 'wall': every non-nozzle patch is a plain 'patch', "
            "so particle::hitBoundaryFace deletes molecules that reach it"))

    if patches[nozzle].type != "patch":
        report.add(CheckResult(
            "nozzle_patch_type", "fail",
            f"the nozzle is geometric type {patches[nozzle].type!r}, not "
            f"'patch'. plumeFieldInflow refuses to inject across a wall, and a "
            f"constraint type would not be an inflow at all."))
    else:
        report.add(CheckResult(
            "nozzle_patch_type", "pass",
            "the nozzle is geometric type 'patch', which plumeFieldInflow "
            "accepts as an injection surface"))


def _check_particle_count(particles: dict, report: CheckReport) -> None:
    """Order of magnitude of the parcel count. WARNING above 1e8.

    A DSMC parcel costs roughly 10^2 bytes in OpenFOAM, so 10^8 parcels is tens
    of gigabytes before the mesh is counted. Worth knowing before a run is
    submitted rather than when the machine starts swapping.
    """
    total = float(particles.get("total_particles", 0.0))
    if total > 1.0e8:
        report.add(CheckResult(
            "particle_count", "warn",
            f"an estimated {total:.3e} parcels. At roughly 10^2 bytes each that "
            f"is tens of gigabytes of parcel storage alone. Raise the particle "
            f"weight (lower resolution.target_particles_per_cell) or coarsen the "
            f"mesh before submitting this.",
            value=total, threshold=1.0e8))
    else:
        report.add(CheckResult(
            "particle_count", "pass",
            f"an estimated {total:.3e} parcels at steady state "
            f"({particles.get('mean_particles_per_cell', 0.0):.2f} per cell on "
            f"average, most of them in the plume core)",
            value=total, threshold=1.0e8))


def write_summary(case_dir: Path, summary: dict) -> Path:
    """Write ``<case>/case-summary.json``.

    JSON rather than YAML: it is what the study table reads, and a generated
    record needs no comments.
    """
    import json

    path = Path(case_dir) / "case-summary.json"
    path.write_text(json.dumps(summary, indent=2, sort_keys=True, default=float)
                    + "\n", encoding="utf-8")
    return path
