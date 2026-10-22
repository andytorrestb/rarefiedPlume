r"""Is a case fit for production, and what will it cost -- before it is meshed.

Every number here comes from the existing planners: :func:`plumetools.cai2012.mesh.plan`,
:func:`plumetools.cai2012.inflow.derive_run_settings`,
:func:`plumetools.cai2012.mesh.estimate_particles`,
:func:`plumetools.cai2012.checks.run_checks` and the cost model the health
family fitted to its measured runs (:mod:`plumetools.cai2012.health`). This
module only reads them against a list of criteria and writes down which ones a
case misses.

The verdict
-----------
``production-candidate``
    every pre-run criterion below is met. That is **not** a validated result:
    convergence and statistical health still have to be demonstrated from the
    sampled output (``Allaudit.py``, ``Allconverge.py``).
``unsuitable-for-production``
    at least one criterion is missed. The configuration is generated, recorded
    and can be meshed, but ``./AllrunCases`` will not run it without
    ``--allow-unsuitable``, and anything it produces is labelled
    non-production. **A finer parcel population does not fix any of the
    numerical criteria**: occupancy is statistics, the cell and the step are
    discretisation.

The criteria: the nozzle-plume family's pre-run checks
------------------------------------------------------
The list is :mod:`plumetools.nozzleplume.checks`' -- the same names, the same
definitions and the same default thresholds -- so a Cai-extension case and an
R-4D case are judged by one set of resolution rules. Two rows are this family's
own (``dsmc_checks``, ``mesh_budget``).

============================  ==========  =======================================
``dsmc_checks``               numerical   no hard failure from
                                          ``plumetools.cai2012.checks`` (Courant
                                          above 1, nozzle faces and area, ...)
``mesh_budget``               numerical   the mesh was NOT coarsened to fit
                                          ``mesh.max_cells``
``cell_over_mean_free_path``  numerical   exit cell <= ``checks.max_cell_over_mfp``
                                          exit mean free paths (Cai's dx/lambda0)
``time_step_collision_time``  numerical   deltaT <= ``max_time_step_over_
                                          collision_time`` exit collision times
``time_step_courant``         numerical   a molecule at U0 + 3 sigma crosses at
                                          most ``checks.warn_courant`` of the
                                          smallest cell
``transient_length``          numerical   averaging starts after the flow has
                                          settled: one flight of the jet across
                                          the domain, there being no plate
``internal_energy_exchange``  numerical   no internal modes left frozen by the
                                          collision model
``particles_per_cell``        statistics  exit-cell parcels >=
                                          ``checks.min_particles_per_cell``
``parcel_count``              resource    estimated parcels <= ``warn_parcels``
============================  ==========  =======================================

Thresholds the cai2012 ``checks`` section has are read from the case; the two it
has no key for -- ``max_time_step_over_collision_time`` and ``warn_parcels`` --
are the nozzle-plume family's defaults, :data:`NOZZLE_CHECKS`, imported rather
than retyped. Neither is a machine limit.

Cost
----
Wall clock is the health family's fitted two-term model,
``(a * parcels + b * cells) * steps``, measured on that family's machine and
rank count. It was fitted to Kn = 100 runs, where the solver reports "No
collisions" on most steps, so for the collisional cases here it is a **lower
bound**. Disk is the health family's measured bytes per cell and per parcel.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from pathlib import Path

from plumetools.cai2012 import health, mesh
from plumetools.nozzleplume.config import ChecksConfig as NozzleChecksConfig

#: The nozzle-plume family's default check thresholds.
NOZZLE_CHECKS = NozzleChecksConfig()

#: deltaT may not exceed this many exit mean collision times.
#:
#: ``plumetools.nozzleplume``'s ``checks.max_time_step_over_collision_time``
#: default: the standard DSMC requirement that collisions be uncoupled from
#: motion only over a step short against the time between them. Against the
#: case's OWN exit collision time ``lambda0 / c_bar`` -- the densest state in
#: the domain, so everywhere else is easier. Cai's ``dt/t0 = 1`` uses the
#: Kn = 0.01 reference instead and is reported alongside, not substituted.
MAX_DT_OVER_COLLISION_TIME = float(NOZZLE_CHECKS.max_time_step_over_collision_time)

#: Estimated parcels above which a case is refused as a resource risk:
#: ``plumetools.nozzleplume``'s ``checks.warn_parcels`` default.
PARCEL_LIMIT = float(NOZZLE_CHECKS.warn_parcels)

#: Bytes per parcel in memory, the order of magnitude ``checks.py`` uses.
PARCEL_MEMORY_BYTES = 1.0e2

#: Reconstructed bytes per written field set, per cell, and per parcel for the
#: cloud. MEASURED on ``cases/cai2012/Cases/Kn100`` and used by
#: ``cases/cai2012-health/generate_cases.py``: 255 MB of fields for 1.31e6
#: cells, 200 MB of cloud for 1.2e6 parcels.
FIELD_BYTES_PER_CELL = 255.0e6 / 1.31e6
CLOUD_BYTES_PER_PARCEL = 200.0e6 / 1.2e6

STATUS_CANDIDATE = "production-candidate"
STATUS_UNSUITABLE = "unsuitable-for-production"


@dataclass(frozen=True)
class Criterion:
    """One pre-run criterion and how the case fared against it."""

    name: str
    kind: str
    value: float | None
    limit: float | None
    passed: bool
    message: str

    def as_dict(self) -> dict:
        return {"name": self.name, "kind": self.kind,
                "value": None if self.value is None else float(self.value),
                "limit": None if self.limit is None else float(self.limit),
                "passed": bool(self.passed), "message": self.message}


def criterion_mesh(cfg, geom, exit_state) -> mesh.MeshPlan:
    """The mesh the collision criterion asks for, with no cell budget.

    Planned and never written: it exists to say how far the budgeted mesh is
    from it. Cheap -- :func:`plumetools.cai2012.mesh.plan` builds segment
    counts, not cell arrays.
    """
    unbounded = dataclasses.replace(
        cfg, mesh=dataclasses.replace(cfg.mesh, max_cells=10 ** 18))
    return mesh.plan(unbounded, geom, exit_state)


def load_cost_model(family_root: Path, health_root: Path):
    """The cost model to estimate with, and where it came from.

    This family's own measurements once ``./AllrunCases`` has recorded any;
    until then the health family's, which were measured on the same solver,
    case template and machine.

    Returns:
        ``(model, origin)``.
    """
    own = Path(family_root) / "results" / "cost-model.yaml"
    model = health.load_cost_model(own)
    if model.n_measurements > 0:
        return model, str(own)
    inherited = Path(health_root) / "results" / "cost-model.yaml"
    return health.load_cost_model(inherited), str(inherited)


def record_measurement(path: Path, entry: dict) -> Path:
    """Add one measured run to this family's ``results/cost-model.yaml``.

    The health family's bookkeeping -- a resumed case replaces its entry, and
    the model is refitted with :func:`plumetools.cai2012.health.fit_cost_model`
    -- under this family's header. :func:`load_cost_model` prefers this file as
    soon as it holds a measurement, so the estimates move from the
    collisionless health runs to this family's own.

    Not locked: jobs that can finish together call it through
    :func:`plumetools.cai2012ext.hpc.record_entry`, which is.
    """
    import yaml

    path = Path(path)
    document = {}
    if path.is_file():
        try:
            document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        except yaml.YAMLError:
            document = {}
    measurements = [m for m in (document.get("measurements") or [])
                    if m.get("case") != entry.get("case")]
    measurements.append(dict(entry))
    measurements.sort(key=lambda m: str(m.get("case", "")))
    model = health.fit_cost_model(measurements)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# GENERATED by cases/cai2012-kn-extension/AllrunCases -- do not edit.\n"
        "# Measured wall clock per case, and the cost model fitted to it.\n"
        "# Read by ./generate_cases.py in place of cases/cai2012-health's.\n\n"
        + yaml.safe_dump({"measurements": measurements,
                          "model": model.as_dict()},
                         sort_keys=False, default_flow_style=False),
        encoding="utf-8", newline="\n")
    return path


def assess(cfg, geom, exit_state, plan, run, nozzle, particles, report, *,
           cost_model=None, cost_origin: str | None = None) -> dict:
    """Read one fully derived case against the production criteria.

    Args:
        cfg, geom, exit_state, plan, run, nozzle, particles, report: exactly
            what ``generate_cases.derive`` returns for the case.
        cost_model: a :class:`plumetools.cai2012.health.CostModel`, or ``None``
            to skip the wall-clock estimate.
        cost_origin: where the model's measurements live, for the record.

    Returns:
        A dict for the manifest: ``status``, ``production_ready``,
        ``unresolved`` (one line per missed criterion), ``criteria``, and the
        ``resolution``, ``time_step``, ``particles`` and ``cost`` figures the
        verdict was reached from.
    """
    criteria = []

    def add(name, kind, value, limit, passed, message):
        criteria.append(Criterion(name, kind, value, limit, bool(passed), message))

    failures = [r.name for r in report.failures]
    add("dsmc_checks", "numerical", float(len(failures)), 0.0, not failures,
        "plumetools.cai2012.checks reports no hard failure" if not failures
        else f"plumetools.cai2012.checks FAILS {failures}")

    ideal = criterion_mesh(cfg, geom, exit_state)
    add("mesh_budget", "numerical", float(ideal.n_cells), float(plan.max_cells),
        not plan.coarsened,
        f"the criterion needs {ideal.n_cells:.3e} cells at a "
        f"{ideal.core_cell_size_m:.4e} m core; mesh.max_cells = "
        f"{plan.max_cells:,} " + (
            f"forced a {plan.core_cell_size_m:.4e} m core instead "
            f"({plan.core_cell_size_m / ideal.core_cell_size_m:.3g}x coarser)"
            if plan.coarsened else "accommodates it"))

    # --- the nozzle-plume family's checks, in its order -----------------------
    ratio = float(plan.cell_over_mfp)
    ratio_limit = float(cfg.checks.max_cell_over_mfp)
    add("cell_over_mean_free_path", "numerical", ratio, ratio_limit,
        ratio <= ratio_limit,
        f"the exit cell {plan.core_cell_size_m:.4e} m is {ratio:.3g} exit mean "
        f"free paths (lambda0 = {plan.mean_free_path_m:.4e} m); the criterion "
        f"is {ratio_limit:g}")

    exit_time = float(exit_state.collision_time_s())
    dt_ratio = float(run.delta_t_s) / exit_time
    add("time_step_collision_time", "numerical", dt_ratio,
        MAX_DT_OVER_COLLISION_TIME, dt_ratio <= MAX_DT_OVER_COLLISION_TIME,
        f"deltaT {run.delta_t_s:.4e} s is {dt_ratio:.3g} exit mean collision "
        f"times (tau0 = lambda0/c_bar = {exit_time:.4e} s); the criterion is "
        f"{MAX_DT_OVER_COLLISION_TIME:g}")

    courant_limit = float(cfg.checks.warn_courant)
    add("time_step_courant", "numerical", float(run.courant), courant_limit,
        float(run.courant) <= courant_limit,
        f"a molecule at U0 + 3 sigma crosses {run.courant:.3f} of the smallest "
        f"cell ({plan.min_cell_size_m:.4e} m) per step")

    # No plate: the flow is settled once the jet has flown the domain once --
    # plumetools.nozzleplume.run's settling time, flow_length / U.
    settling = float(run.domain_transit_s)
    add("transient_length", "numerical", float(run.average_start_s), settling,
        float(run.average_start_s) >= settling,
        f"averaging starts at {run.average_start_s:.4e} s, "
        f"{run.average_start_s / settling:.3g} times the {settling:.4e} s one "
        f"flight of the jet across the domain takes")

    dof = int(cfg.gas.internal_degrees_of_freedom)
    model = cfg.dsmc.binary_collision_model
    frozen = dof > 0 and model != "LarsenBorgnakkeVariableHardSphere"
    add("internal_energy_exchange", "numerical", float(dof), 0.0, not frozen,
        f"{dof} internal degree(s) of freedom with {model}"
        + (" -- internal energy would stay frozen" if frozen else ""))

    occupancy = float(run.exit_particles_per_cell)
    floor = float(cfg.checks.min_particles_per_cell)
    add("particles_per_cell", "statistics", occupancy, floor, occupancy >= floor,
        f"{occupancy:.1f} parcels in the exit cell by construction (floor "
        f"{floor:g}); downstream occupancy falls with the density and is "
        f"measured after the run (Allaudit.py, Alldiagnose.py)")

    parcels = float(particles["total_particles"])
    add("parcel_count", "resource", parcels, PARCEL_LIMIT, parcels <= PARCEL_LIMIT,
        f"an estimated {parcels:.3e} parcels at steady state (threshold "
        f"{PARCEL_LIMIT:.3g})")

    unresolved = [f"{c.name}: {c.message}" for c in criteria if not c.passed]
    ready = not unresolved

    steps = int(run.n_steps)
    fields_bytes = FIELD_BYTES_PER_CELL * float(plan.n_cells)
    cloud_bytes = CLOUD_BYTES_PER_PARCEL * parcels
    writes = int(run.n_writes)
    cost = {
        "n_cells": int(plan.n_cells),
        "estimated_parcels": parcels,
        "n_steps": steps,
        "parcel_steps": parcels * steps,
        "cell_steps": float(plan.n_cells) * steps,
        "n_writes": writes,
        "n_sampled_writes": int(run.n_sampled_writes),
        "n_transient_writes": writes - int(run.n_sampled_writes),
        "fields_gb_per_write": fields_bytes / 1024.0 ** 3,
        "cloud_gb_per_write": cloud_bytes / 1024.0 ** 3,
        # processor*/ keeps fields and cloud for every write; reconstructPar
        # duplicates the fields, and AllrunCases prunes the reconstructed
        # clouds of all but the newest write.
        "disk_gb_estimate": writes * (2.0 * fields_bytes + cloud_bytes)
        / 1024.0 ** 3,
        "parcel_memory_gb_order": parcels * PARCEL_MEMORY_BYTES / 1024.0 ** 3,
        "disk_basis": "cases/cai2012-health measured bytes per cell and per "
                      "parcel; every write kept (purgeWrite 0)",
    }
    if cost_model is not None:
        seconds = cost_model.seconds(parcels=parcels, cells=plan.n_cells,
                                     steps=steps)
        cost.update({
            "predicted_wall_hours": seconds / 3600.0,
            "cost_model": cost_model.as_dict(),
            "cost_model_origin": cost_origin,
            "cost_caveat": "fitted to collisionless (Kn = 100) runs on the "
                           "health family's machine and rank count; collision "
                           "work is not in it, so collisional cases will take "
                           "longer than this",
        })

    return {
        "status": STATUS_CANDIDATE if ready else STATUS_UNSUITABLE,
        "production_ready": ready,
        "unresolved": unresolved,
        "criteria": [c.as_dict() for c in criteria],
        "resolution": {
            "n_cells": int(plan.n_cells),
            "core_cell_size_m": float(plan.core_cell_size_m),
            "min_cell_size_m": float(plan.min_cell_size_m),
            "requested_cell_size_m": float(plan.requested_cell_size_m),
            "coarsened_for_budget": bool(plan.coarsened),
            "max_cells_budget": int(plan.max_cells),
            "cell_over_mean_free_path": ratio,
            "cells_to_meet_criterion": float(ideal.n_cells),
            "nozzle_faces": int(nozzle.n_faces),
            "nozzle_area_error": float(nozzle.area_error),
        },
        "time_step": {
            "delta_t_s": float(run.delta_t_s),
            "courant": float(run.courant),
            "exit_collision_time_s": exit_time,
            "dt_over_exit_collision_time": dt_ratio,
            "reference_collision_time_s": float(run.reference_collision_time_s),
            "dt_over_reference_collision_time": float(run.delta_t_s)
            / float(run.reference_collision_time_s),
            "cai_delta_t_s": float(run.cai_delta_t_s),
            "cai_courant": float(run.cai_courant),
            "settling_time_s": settling,
            "average_start_over_settling_time": float(run.average_start_s)
            / settling,
        },
        "particles": {
            "estimated_parcels": parcels,
            "exit_particles_per_cell": occupancy,
            "mean_particles_per_cell": float(particles["mean_particles_per_cell"]),
            "estimate_basis": "collisionless density integrated over the box "
                              "(plumetools.cai2012.mesh.estimate_particles)",
        },
        "cost": cost,
        "note": ("production-candidate means every pre-run criterion is met. It "
                 "is not a validated result: convergence and occupancy must be "
                 "shown from the sampled output."),
    }


def describe(assessment: dict) -> list[str]:
    """Report lines for one assessed case."""
    resolution = assessment["resolution"]
    step = assessment["time_step"]
    cost = assessment["cost"]
    lines = [
        f"          verdict  {assessment['status'].upper()}",
        f"          mesh     {resolution['n_cells']:,} cells, core "
        f"{resolution['core_cell_size_m']:.4e} m, cell/lambda0 "
        f"{resolution['cell_over_mean_free_path']:.3g}"
        + ("  (COARSENED to the budget; the criterion needs "
           f"{resolution['cells_to_meet_criterion']:.2e} cells)"
           if resolution["coarsened_for_budget"] else ""),
        f"          step     deltaT {step['delta_t_s']:.4e} s, Courant "
        f"{step['courant']:.3f}, dt/tau0 {step['dt_over_exit_collision_time']:.3g}",
        f"          cost     ~{cost['estimated_parcels']:.2e} parcels x "
        f"{cost['n_steps']:,} steps, {cost['n_writes']} writes "
        f"({cost['n_sampled_writes']} sampled), ~{cost['disk_gb_estimate']:.0f} GB"
        + (f", >= {cost['predicted_wall_hours']:.1f} h"
           if "predicted_wall_hours" in cost else ""),
    ]
    for reason in assessment["unresolved"]:
        lines.append(f"          UNRESOLVED {reason}")
    return lines
