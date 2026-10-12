r"""DSMC mesh, time-step and boundary-type quality checks.

Most of these **warn**, and for a real nozzle exit several of them always will.
An exit at a Knudsen number of :math:`10^{-4}` has a mean free path thousands of
times smaller than its diameter; a mesh that resolved it would have
:math:`10^{10}` cells in the exit plane alone. So the cell-size and
collision-time criteria are stated, measured and reported -- with how far off
they are and what that costs -- rather than enforced. A check that could only
ever fail teaches nobody anything, and a run that hides how under-resolved it is
teaches something false.

These are **errors**, because a case with any of them cannot produce a
meaningful answer at all:

* a non-positive or non-finite density, temperature, speed or time step;
* a time step that lets a molecule cross more than one of the smallest cells;
* a missing required patch, or a nozzle patch with no faces;
* a nozzle patch that is not the one the plan predicted -- the mesh is then not
  the mesh that was checked;
* a nozzle patch whose area misses the exit area by more than
  ``checks.max_nozzle_area_error`` -- the injected mass flow is proportional to
  it, so the whole solution scales with the error;
* an open boundary declared ``wall`` -- the plume would expand into a reflecting
  box rather than into vacuum;
* the plate declared anything but ``wall`` -- it would absorb the plume and feel
  no load;
* the nozzle declared ``wall`` -- ``plumeFieldInflow`` refuses to inject through
  one, so this fails anyway, but much later and much less clearly.

Every result is machine-readable, so a case summary carries the whole set
rather than a boolean.
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

    def add(self, name: str, status: str, message: str, value=None,
            threshold=None) -> None:
        self.results.append(CheckResult(name, status, message, value, threshold))

    @property
    def failures(self) -> list:
        return [r for r in self.results if r.status == "fail"]

    @property
    def warnings(self) -> list:
        return [r for r in self.results if r.status == "warn"]

    def raise_if_failed(self) -> None:
        """Raise :class:`DsmcCheckError` listing **every** failure.

        All of them, not the first: they usually share a cause, and fixing them
        one run at a time is slow.
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
            lines.append(
                f"  {len(self.warnings)} warning(s). A warning is a statement about "
                f"accuracy, not a broken case: the run proceeds.")
        return lines

    def as_dict(self) -> dict:
        return {
            "passed": not self.failures,
            "n_failures": len(self.failures),
            "n_warnings": len(self.warnings),
            "results": [r.as_dict() for r in self.results],
        }


def run_checks(cfg, geom, exit_state, plan=None, run=None, *, nozzle=None,
               patches=None) -> CheckReport:
    """Run every configured check.

    Args:
        cfg: the case config.
        geom: the :class:`~plumetools.nozzleplume.geometry.PlumeGeometry`.
        exit_state: the :class:`~plumetools.nozzleplume.exitstate.ExitState`.
        plan: the :class:`~plumetools.nozzleplume.mesh.MeshPlan`, if there is one.
        run: the :class:`~plumetools.nozzleplume.run.RunSettings`, if there is one.
        nozzle: the :class:`~plumetools.nozzleplume.mesh.NozzlePrediction`.
        patches: the mesh's patches, once it exists.

    Returns:
        A :class:`CheckReport`. Nothing is raised here -- call
        :meth:`CheckReport.raise_if_failed` so a caller can record the full
        report first.
    """
    report = CheckReport()

    if not cfg.checks.enabled:
        report.add("checks_enabled", "warn",
                   "checks.enabled is false, so no DSMC quality check was run")
        return report

    _check_exit_state(exit_state, report)
    _check_gas_model(cfg, report)
    if plan is not None:
        _check_cell_over_mfp(cfg, plan, report)
    if plan is not None and run is not None:
        _check_time_step(cfg, plan, run, report)
        _check_occupancy(cfg, run, report)
        _check_parcel_count(cfg, run, report)
        _check_transient(geom, run, report)
    if nozzle is not None:
        _check_nozzle(cfg, nozzle, report)
    if patches is not None:
        _check_patches(cfg, geom, patches, nozzle, report)
    return report


def _check_exit_state(exit_state, report: CheckReport) -> None:
    """Every exit quantity finite and positive. HARD ERROR."""
    values = {
        "number density": exit_state.number_density_per_m3,
        "temperature": exit_state.temperature_K,
        "speed": exit_state.speed_m_per_s,
        "mean free path": exit_state.mean_free_path_m,
        "number flux": exit_state.number_flux_per_m2_s,
    }
    bad = sorted(name for name, v in values.items()
                 if v is None or not math.isfinite(float(v)) or float(v) <= 0.0)
    if bad:
        report.add(
            "exit_state_positive", "fail",
            f"{bad} must be finite and positive; check the exit and gas "
            f"sections of case.yaml")
    else:
        report.add(
            "exit_state_positive", "pass",
            f"n = {exit_state.number_density_per_m3:.4e} 1/m^3, "
            f"U = {exit_state.speed_m_per_s:.1f} m/s, "
            f"T = {exit_state.temperature_K:g} K -> S = "
            f"{exit_state.speed_ratio:.2f}, Kn = {exit_state.knudsen:.3e}")


def _check_gas_model(cfg, report: CheckReport) -> None:
    """Internal modes against the collision model. WARNING.

    Elastic VHS never exchanges energy with a molecule's internal modes. A gas
    given internal degrees of freedom under it keeps whatever internal energy it
    was injected with for ever -- "frozen" -- which is a model, but rarely the
    one intended.
    """
    dof = int(cfg.gas.internal_degrees_of_freedom)
    model = cfg.dsmc.collision_model
    if dof > 0 and model != "LarsenBorgnakkeVariableHardSphere":
        report.add(
            "internal_energy_exchange", "warn",
            f"gas.internal_degrees_of_freedom is {dof} but dsmc.collision_model "
            f"is {model}, which never exchanges energy with internal modes: each "
            f"molecule keeps the internal energy it was injected with. Use "
            f"LarsenBorgnakkeVariableHardSphere for relaxing internal modes.",
            value=float(dof))
    else:
        report.add(
            "internal_energy_exchange", "pass",
            f"{dof} internal degree(s) of freedom with {model}"
            + (" (monatomic-like: gamma = 5/3 in the simulated gas)"
               if dof == 0 else ""),
            value=float(dof))


def _check_cell_over_mfp(cfg, plan, report: CheckReport) -> None:
    """Exit cell against the exit mean free path. WARNING.

    The exit is the densest place in the domain, so it is the worst case: the
    ratio falls downstream as the plume thins.
    """
    ratio = plan.exit_cell_over_mfp
    limit = float(cfg.checks.max_cell_over_mfp)
    if ratio > limit:
        levels = math.ceil(math.log2(ratio / limit))
        report.add(
            "cell_over_mean_free_path", "warn",
            f"the exit cell is {ratio:.3g} exit mean free paths "
            f"({plan.exit_cell_size_m:.4g} m against "
            f"{plan.mean_free_path_m:.3e} m); the DSMC criterion is {limit:g}. "
            f"Collision partners are drawn from across a cell, so where it is "
            f"this much larger than a mean free path the gas behaves as more "
            f"viscous than it is and gradients finer than a cell are smeared. "
            f"Meeting the criterion at the exit would take {levels} more "
            f"refinement level(s), x{8 ** levels:.3g} cells there. The ratio "
            f"falls with the density downstream.",
            value=ratio, threshold=limit)
    else:
        report.add(
            "cell_over_mean_free_path", "pass",
            f"the exit cell is {ratio:.3g} exit mean free paths, within the "
            f"criterion of {limit:g}",
            value=ratio, threshold=limit)


def _check_time_step(cfg, plan, run, report: CheckReport) -> None:
    """Courant number (HARD ERROR above 1) and collision time (WARNING).

    The characteristic speed is ``U + 3 sigma``, which bounds essentially the
    whole distribution. The mean would let the fast tail skip cells unnoticed.
    """
    courant = run.courant
    speed = run.characteristic_speed_m_per_s
    smallest = plan.min_cell_size_m

    if courant > float(cfg.checks.max_courant):
        report.add(
            "time_step_courant", "fail",
            f"a molecule at {speed:.0f} m/s crosses {courant:.2f} of the "
            f"smallest cell ({smallest:.4g} m) per step of {run.delta_t_s:g} s. "
            f"Above {cfg.checks.max_courant:g} molecules skip cells without "
            f"being offered a collision partner in them. Reduce "
            f"dsmc.courant_target, or set dsmc.delta_t_s below "
            f"{cfg.checks.max_courant * smallest / speed:.3g} s.",
            value=courant, threshold=float(cfg.checks.max_courant))
    elif courant > float(cfg.checks.warn_courant):
        report.add(
            "time_step_courant", "warn",
            f"a molecule crosses {courant:.2f} of the smallest cell per step; "
            f"below {cfg.checks.warn_courant:g} is preferable",
            value=courant, threshold=float(cfg.checks.warn_courant))
    else:
        report.add(
            "time_step_courant", "pass",
            f"a molecule at {speed:.0f} m/s crosses {courant:.3f} of the "
            f"smallest cell ({smallest:.4g} m) per step",
            value=courant, threshold=float(cfg.checks.warn_courant))

    ratio = run.time_step_over_collision_time
    limit = float(cfg.checks.max_time_step_over_collision_time)
    if ratio > limit:
        report.add(
            "time_step_collision_time", "warn",
            f"deltaT is {ratio:.3g} exit mean collision times "
            f"({run.delta_t_s:.3e} s against {run.collision_time_s:.3e} s); the "
            f"DSMC criterion is {limit:g}. A step this long couples moving and "
            f"colliding at the exit, where a molecule should have collided "
            f"several times within it. Like the cell size, this eases "
            f"downstream as the collision time grows.",
            value=ratio, threshold=limit)
    else:
        report.add(
            "time_step_collision_time", "pass",
            f"deltaT is {ratio:.3g} exit mean collision times",
            value=ratio, threshold=limit)


def _check_occupancy(cfg, run, report: CheckReport) -> None:
    """Parcels in the exit cell. WARNING.

    Computed, not sampled, so never an error: only a post-run look at the
    sampled ``dsmcRhoNMean`` field can say what the occupancy actually was.
    """
    value = run.exit_particles_per_cell
    floor = float(cfg.checks.min_particles_per_cell)
    if value < floor:
        report.add(
            "particles_per_cell", "warn",
            f"an exit cell holds {value:.3g} parcels, below the {floor:g} floor. "
            f"Collisions are sampled between parcels in the same cell, so with "
            f"fewer than a handful the collision rate itself becomes noisy. "
            f"nEquivalentParticles = "
            f"{run.n_equivalent_particles * value / floor:.3e} would put "
            f"{floor:g} there -- at {floor / value:.3g}x the parcels everywhere.",
            value=value, threshold=floor)
    else:
        report.add(
            "particles_per_cell", "pass",
            f"an exit cell holds {value:.3g} parcels (floor {floor:g}). One "
            f"weight serves the whole domain, so occupancy elsewhere follows the "
            f"local density and cell volume.",
            value=value, threshold=floor)


def _check_parcel_count(cfg, run, report: CheckReport) -> None:
    """Order of magnitude of the parcel count. WARNING.

    A DSMC parcel costs roughly 10^2 bytes in OpenFOAM, so 10^8 parcels is tens
    of gigabytes before the mesh is counted. Worth knowing before a run is
    started rather than when the machine starts swapping.
    """
    total = float(run.estimated_parcels)
    limit = float(cfg.checks.warn_parcels)
    if total > limit:
        report.add(
            "parcel_count", "warn",
            f"roughly {total:.2e} parcels at steady state, over the "
            f"{limit:.2g} threshold -- about {total * 100.0 / 1.0e9:.1f} GB of "
            f"parcel storage at ~100 bytes each, and a step time to match. "
            f"Raise dsmc.n_equivalent_particles to thin them out.",
            value=total, threshold=limit)
    else:
        report.add(
            "parcel_count", "pass",
            f"roughly {total:.2e} parcels at steady state (an order-of-magnitude "
            f"estimate; the solver log is the measurement)",
            value=total, threshold=limit)


def _check_transient(geom, run, report: CheckReport) -> None:
    """Whether averaging starts after the flow has had time to settle. WARNING."""
    settle = run.settling_time_s
    start = run.average_start_s
    if start < settle:
        report.add(
            "transient_length", "warn",
            f"averaging starts at {start:.3e} s, before the roughly "
            f"{settle:.3e} s the flow needs to settle"
            + (" (the jet's flight to the plate plus the time gas is held up "
               "there)" if geom.plate is not None else
               " (one flight of the jet across the domain)")
            + ". The mean fields would include the start-up. Raise "
              "dsmc.transient_domain_transits or dsmc.average_start_s.",
            value=start, threshold=settle)
    else:
        report.add(
            "transient_length", "pass",
            f"averaging starts at {start:.3e} s, after the roughly "
            f"{settle:.3e} s the flow needs to settle",
            value=start, threshold=settle)


def _check_nozzle(cfg, nozzle, report: CheckReport) -> None:
    """Face count and area of the staircase nozzle. HARD ERROR on both."""
    if nozzle.n_faces < int(cfg.checks.min_nozzle_faces):
        report.add(
            "nozzle_faces", "fail",
            f"the nozzle patch would have {nozzle.n_faces} face(s), below the "
            f"{cfg.checks.min_nozzle_faces} minimum. A disk resolved by a "
            f"handful of squares is not a circular exit. Put a refinement "
            f"region over the exit, or raise its level.",
            value=float(nozzle.n_faces),
            threshold=float(cfg.checks.min_nozzle_faces))
    else:
        report.add(
            "nozzle_faces", "pass",
            f"the nozzle patch has {nozzle.n_faces} faces, "
            f"{nozzle.cells_across_diameter:.1f} cells across the diameter",
            value=float(nozzle.n_faces),
            threshold=float(cfg.checks.min_nozzle_faces))

    error = abs(nozzle.area_error)
    limit = float(cfg.checks.max_nozzle_area_error)
    if error > limit:
        report.add(
            "nozzle_area", "fail",
            f"the meshed nozzle area is {nozzle.area_m2:.6e} m^2 against the "
            f"exact {nozzle.exact_area_m2:.6e} m^2, an error of "
            f"{100.0 * nozzle.area_error:+.2f}% -- past the "
            f"{100.0 * limit:.1f}% limit. The injected mass flow is "
            f"proportional to the area, so the whole solution scales with this. "
            f"The error wanders as the exit cell shrinks -- it is not monotone "
            f"-- but its bound falls with the cell, so try the next level.",
            value=error, threshold=limit)
    else:
        report.add(
            "nozzle_area", "pass",
            f"the meshed nozzle area is {nozzle.area_m2:.6e} m^2, "
            f"{100.0 * nozzle.area_error:+.2f}% from the exact area -- within "
            f"the {100.0 * limit:.1f}% limit for a staircase disk",
            value=error, threshold=limit)


def expected_patch_types(cfg, geom) -> dict:
    """Patch name to the geometric type the built mesh must give it."""
    names = cfg.mesh.patch_names
    expected = {
        names["nozzle"]: "patch",
        names["upstream_vacuum"]: "patch",
        names["outer"]: "patch",
    }
    if geom.plate is not None:
        expected[names["plate"]] = "wall"
    if geom.symmetry_planes:
        expected[names["symmetry"]] = "symmetry"
    return expected


def _check_patches(cfg, geom, patches, nozzle, report: CheckReport) -> None:
    """Required patches present, non-empty, and of the right geometric type."""
    names = cfg.mesh.patch_names
    expected = expected_patch_types(cfg, geom)

    missing = sorted(name for name in expected if name not in patches)
    if missing:
        report.add(
            "required_patches", "fail",
            f"the mesh is missing patch(es) {missing}; it has {sorted(patches)}. "
            f"'{names['nozzle']}' is created by topoSet + createPatch and "
            f"'{names['plate']}' by snappyHexMesh, not by blockMesh -- re-run "
            f"./Allmesh.")
        return
    report.add(
        "required_patches", "pass",
        f"all {len(expected)} required patches present: "
        + ", ".join(f"{n}({patches[n].n_faces})" for n in expected))

    # Geometric types. These decide the physics, not the labelling.
    reasons = {
        "patch": "an open boundary or an inlet must be 'patch': a molecule "
                 "reaching a plain patch is deleted, and plumeFieldInflow will "
                 "not inject across a wall",
        "wall": "a body must be 'wall', or hitWallPatch never runs: the plate "
                "would absorb the plume and feel no load",
        "symmetry": "a symmetry plane must be 'symmetry', or molecules that "
                    "reach it are deleted instead of reflected",
    }
    wrong = [(name, kind, patches[name].type) for name, kind in expected.items()
             if patches[name].type != kind]
    if wrong:
        report.add(
            "patch_types", "fail",
            "; ".join(f"{name!r} is geometric type {actual!r}, not {kind!r} -- "
                      f"{reasons[kind]}" for name, kind, actual in wrong))
    else:
        report.add(
            "patch_types", "pass",
            ", ".join(f"{name} is '{kind}'" for name, kind in expected.items()))

    stray = sorted(name for name, info in patches.items()
                   if info.type == "wall" and expected.get(name) != "wall")
    if stray:
        report.add(
            "no_stray_walls", "fail",
            f"patch(es) {stray} are geometric type 'wall' but are not the plate. "
            f"A wall REFLECTS molecules, so the plume would expand into a "
            f"partly closed box instead of into vacuum.")

    empty = sorted(name for name in expected if patches[name].n_faces == 0)
    if empty:
        report.add(
            "patches_not_empty", "fail",
            f"patch(es) {empty} have no faces. An empty nozzle injects nothing; "
            f"an empty plate means snappyHexMesh carved nothing, so there is no "
            f"plate in the mesh.")

    if nozzle is not None:
        meshed = patches[names["nozzle"]].n_faces
        if meshed != nozzle.n_faces:
            report.add(
                "nozzle_matches_plan", "fail",
                f"the nozzle patch has {meshed} faces where the plan predicted "
                f"{nozzle.n_faces}. The Python rule (face centre inside "
                f"r <= R0, on the exit-level lattice) and topoSet's "
                f"cylinderToFace are the same rule, so the mesh is not the mesh "
                f"that was planned -- most likely the exit is not at one "
                f"refinement level. Re-run ./Allmesh.",
                value=float(meshed), threshold=float(nozzle.n_faces))
        else:
            report.add(
                "nozzle_matches_plan", "pass",
                f"the nozzle patch has the {meshed} faces the plan predicted",
                value=float(meshed), threshold=float(nozzle.n_faces))


def write_summary(case_dir: Path, summary: dict) -> Path:
    """Write ``<case>/case-summary.json``.

    JSON rather than YAML: a generated record needs no comments, and anything
    downstream can read it.
    """
    import json

    path = Path(case_dir) / "case-summary.json"
    path.write_text(json.dumps(summary, indent=2, sort_keys=True, default=float)
                    + "\n", encoding="utf-8")
    return path
