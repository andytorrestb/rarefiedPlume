"""The pre-run verdict of an interface case: criteria, resources, cost.

Reads a :class:`~plumetools.cai2012ext.interface_case.Derived` and writes
nothing. Numerical suitability and coupling validation are judged separately
(:mod:`plumetools.cai2012ext.interface_evidence`) and neither implies the
other; the verdicts are described in :mod:`plumetools.cai2012ext.interface_case`.
Split out of that module.
"""

from __future__ import annotations

import numpy as np

from plumetools.cai2012ext import feasibility
from plumetools.cai2012ext.interface_evidence import coupling_validation

STATUS_INFEASIBLE = "interface-infeasible"
STATUS_UNSUITABLE = feasibility.STATUS_UNSUITABLE
STATUS_PENDING = "pending-coupling-validation"
STATUS_CANDIDATE = feasibility.STATUS_CANDIDATE
STATUS_VALIDATION = "interface-validation"
STATUS_SMOKE = "development-smoke"


def assess(derived, *, cost_model=None, cost_origin=None) -> dict:
    """The verdict, the criteria behind it, and what the case would cost."""
    settings = derived.settings
    cfg = derived.cfg
    selection = derived.selection
    criteria = []

    def add(name, kind, value, limit, passed, message, blocking=True):
        criteria.append({"name": name, "kind": kind,
                         "value": None if value is None else float(value),
                         "limit": None if limit is None else float(limit),
                         "passed": bool(passed), "blocking": bool(blocking),
                         "message": message})

    missing = settings.missing_decisions() if settings.scientific else []
    add("configuration_complete", "configuration", float(len(missing)), 0.0,
        not missing, "every scientific decision is stated" if not missing
        else " ".join(missing))
    add("interface_selection", "interface", selection.max_core_bird_at_radius
        if selection.admissible else selection.max_core_bird_at_r_min,
        selection.bird_threshold, selection.admissible,
        f"{selection.status}: {selection.reason}")
    if selection.admissible and derived.mesh_error:
        add("interface_mesh", "numerical", None, None, False, derived.mesh_error)

    planned = derived.planned
    if planned is not None:
        plan, run, ref = planned.mesh_plan, planned.run, planned.reference
        add("interface_radius_set_by", "interface", None, None, True,
            f"R_interface is set by the {selection.limited_by}"
            + ("" if selection.limited_by == "bird" else
               " -- the B bound would allow a larger sphere; enlarge "
               "geometry.* (and diagnostic_radius_over_D) explicitly if more "
               "of the far plume should be continuum-injected"), blocking=False)
        add("mesh_budget", "numerical", float(plan.retained_cells),
            float(plan.max_cells), not plan.coarsened,
            f"{plan.retained_cells:,} retained cells at a "
            f"{plan.core_cell_size_m:.4e} m core"
            + (f"; lambda at the interface asked for "
               f"{plan.requested_cell_size_m:.4e} m, which needs "
               + (f"{planned.resolution['cells_to_meet_criterion']:.3e}"
                  if planned.resolution.get("cells_to_meet_criterion") else "more than")
               + f" retained cells against mesh.max_cells = {plan.max_cells:,}"
               if plan.coarsened else ""))
        limit = float(cfg.checks.max_cell_over_mfp)
        add("cell_over_mean_free_path", "numerical", plan.cell_over_mfp, limit,
            plan.cell_over_mfp <= limit,
            f"the interface core cell {plan.core_cell_size_m:.4e} m is "
            f"{plan.cell_over_mfp:.3g} closure mean free paths at (R, 0)")
        worst = max(planned.resolution["axis"]["max_cell_over_lambda"],
                    planned.resolution["core_edge"]["max_cell_over_lambda"])
        add("downstream_resolution_estimate", "estimate", worst, limit,
            worst <= limit,
            f"largest planned cell/lambda along the axis and core edge is "
            f"{worst:.3g} (continuum estimate; measured by Alldiagnose.py)",
            blocking=False)
        tau = float(ref["collision_time_s"])
        ratio = float(run.delta_t_s) / tau
        add("time_step_collision_time", "numerical", ratio,
            feasibility.MAX_DT_OVER_COLLISION_TIME,
            ratio <= feasibility.MAX_DT_OVER_COLLISION_TIME * (1.0 + 1e-12),
            f"deltaT {run.delta_t_s:.4e} s is {ratio:.3g} reference-cell "
            f"collision times (tau = {tau:.4e} s; source: {run.delta_t_source})")
        courant_limit = float(cfg.checks.warn_courant)
        add("time_step_courant", "numerical", float(run.courant), courant_limit,
            float(run.courant) <= courant_limit,
            f"a molecule at |U| + 3 sqrt(RT) of the interface state "
            f"({ref['characteristic_speed_m_per_s']:.1f} m/s) crosses "
            f"{run.courant:.3f} of the smallest cell per step")
        add("transient_length", "numerical", float(run.average_start_s),
            float(run.domain_transit_s),
            float(run.average_start_s) >= float(run.domain_transit_s),
            f"averaging starts at {run.average_start_s:.4e} s, "
            f"{run.average_start_s / run.domain_transit_s:.3g} transits "
            f"(x_max/U_e; transient basis '{run.transient_basis}'"
            + ("; Cai's 10^4 reference collision times, unchanged by the "
               "interface)" if run.transient_basis == "cai" else ")"))
        dof = int(cfg.gas.internal_degrees_of_freedom)
        add("internal_energy_exchange", "numerical", float(dof), 0.0, dof == 0,
            f"{dof} internal degree(s) of freedom; the closure is monatomic")
        floor = float(cfg.checks.min_particles_per_cell)
        add("particles_per_cell", "statistics", float(run.exit_particles_per_cell),
            floor, float(run.exit_particles_per_cell) >= floor,
            f"{run.exit_particles_per_cell:.1f} parcels in the declared reference "
            f"cell ({ref['identity']}); occupancy elsewhere is measured after "
            f"the run")
        add("parcel_count", "resource", planned.parcels, feasibility.PARCEL_LIMIT,
            planned.parcels <= feasibility.PARCEL_LIMIT,
            f"an estimated {planned.parcels:.3e} parcels at steady state")
        balance = planned.ideal_flux
        residual = max(abs(balance["net_mass_residual"]),
                       abs(balance["gross_mass_residual"]),
                       abs(balance["gross_energy_residual"]))
        add("ideal_sphere_flux_balance", "coupling", residual,
            settings.meshed_flux_tolerance,
            residual <= settings.meshed_flux_tolerance,
            f"on the exact sphere the core injects gross mass "
            f"{balance['gross_mass_residual']:+.2e} and energy "
            f"{balance['gross_energy_residual']:+.2e} relative to the core "
            f"targets (net mass {balance['net_mass_residual']:+.2e}); the meshed "
            f"patch is checked again by runInflow.py")
        # Not a gate -- the collisionless solution is not the truth at low Kn,
        # the overlap study is -- but the size of the angular-law mismatch is
        # stated before any run, never discovered afterwards.
        from plumetools.cai2012 import analytical

        radius = float(derived.selection.radius_m)
        closure = float(derived.source.state(np.array([radius]),
                                             np.array([0.0]))["n"][0]) / derived.source.n_e
        cai = float(analytical.centerline_density_ratio(
            np.array([radius]), 0.5 * derived.source.diameter_m,
            float(cfg.exit.speed_ratio))[0])
        add("angular_law_axis_density", "coupling", closure / cai, 1.0, True,
            f"at R the closure's axis density is {closure / cai:.3f} x Cai's "
            f"collisionless value (Markelov's sonic-orifice f(theta) is broader "
            f"than this supersonic exit's plume; CP-3 is not tuned -- the "
            f"overlap study measures the consequence)", blocking=False)
        memory = _memory(derived)
        add("memory_ceiling", "resource", memory.get("estimate_bytes"),
            memory.get("ceiling_bytes"), memory["within"] is not False,
            memory["message"], blocking=memory["within"] is False)

    numerical_unresolved = [f"{c['name']}: {c['message']}" for c in criteria
                            if c["blocking"] and not c["passed"]]
    if not derived.meshable:
        numeric_status = "not-assessed"
    else:
        numeric_status = STATUS_CANDIDATE if not numerical_unresolved else STATUS_UNSUITABLE
    numerical = {"status": numeric_status, "unresolved": numerical_unresolved}
    coupling = coupling_validation(derived)

    if not selection.admissible or missing or not derived.meshable:
        status = STATUS_INFEASIBLE
    elif settings.purpose == "development-smoke":
        status = STATUS_SMOKE
    elif settings.purpose == "interface-validation":
        status = STATUS_VALIDATION
    elif numerical_unresolved:
        status = STATUS_UNSUITABLE
    elif coupling["status"] != "passed":
        status = STATUS_PENDING
    else:
        status = STATUS_CANDIDATE

    unresolved = list(numerical_unresolved)
    if coupling["required"] and coupling["status"] != "passed":
        unresolved.append(f"coupling_validation: {coupling['status']} -- "
                          f"{coupling.get('note', coupling.get('evidence', ''))}")
    if settings.purpose != "production":
        unresolved.append(f"purpose: {settings.purpose} -- never a production "
                          f"result")
    if settings.run_mode == "exploratory":
        unresolved.append("run_mode: exploratory -- runs with unresolved "
                          "criteria, labelled NON-PRODUCTION")

    runnable, refusal = _runnable(derived, status, numerical, coupling)
    cost = _cost(derived, cost_model, cost_origin) if derived.meshable else None
    return {
        "status": status,
        "production_ready": status == STATUS_CANDIDATE,
        "purpose": settings.purpose,
        "run_mode": settings.run_mode,
        "runnable": runnable,
        "refusal": refusal,
        "unresolved": unresolved,
        "numerical": numerical,
        "coupling_validation": coupling,
        "criteria": criteria,
        "cost": cost,
        "note": "numerical suitability and coupling validation are separate; "
                "a pending HPC validation stays pending",
    }


def _runnable(derived, status, numerical, coupling) -> tuple[bool, str | None]:
    settings = derived.settings
    if status == STATUS_INFEASIBLE:
        return False, "the interface is infeasible as configured"
    if settings.run_mode == "exploratory":
        return True, None
    if settings.purpose == "production" and status != STATUS_CANDIDATE:
        return False, (f"{status}: a production run needs every numerical "
                       f"criterion met and the coupling validated "
                       f"(coupling: {coupling['status']}). Set run_mode: "
                       f"exploratory to run it labelled NON-PRODUCTION.")
    if settings.purpose == "interface-validation" and numerical["unresolved"]:
        return False, ("an overlap validation on an unresolved mesh measures "
                       "the mesh, not the coupling. Fix the criteria, or set "
                       "run_mode: exploratory.")
    return True, None


def _memory(derived) -> dict:
    resources = derived.raw.get("resources") or {}
    per_cell = resources.get("bytes_per_cell")
    per_parcel = resources.get("bytes_per_parcel")
    ceiling = resources.get("memory_ceiling_bytes")
    if per_cell is None or per_parcel is None:
        return {"within": None, "ceiling_bytes": ceiling,
                "message": "UNCALIBRATED: no measured bytes per cell / parcel "
                           "for this solver at these densities "
                           "(resources.bytes_per_cell, bytes_per_parcel); collect "
                           "them from the smoke run's log and set them"}
    planned = derived.planned
    estimate = (float(per_cell) * planned.mesh_plan.retained_cells
                + float(per_parcel) * planned.parcels)
    within = None if ceiling is None else estimate <= float(ceiling)
    return {"within": within, "estimate_bytes": estimate, "ceiling_bytes": ceiling,
            "message": f"~{estimate / 1024 ** 3:.1f} GiB from the user-supplied "
                       f"calibration" + ("" if ceiling is None else
                                         f" against {float(ceiling) / 1024 ** 3:.0f} GiB")}


def _cost(derived, cost_model, cost_origin) -> dict:
    planned = derived.planned
    plan, run, parcels = planned.mesh_plan, planned.run, planned.parcels
    fields = feasibility.FIELD_BYTES_PER_CELL * plan.retained_cells
    cloud = feasibility.CLOUD_BYTES_PER_PARCEL * parcels
    cost = {
        "n_cells": int(plan.retained_cells),
        "lattice_cells": int(plan.lattice.n_cells),
        "estimated_parcels": float(parcels),
        "n_steps": int(run.n_steps),
        "n_writes": int(run.n_writes),
        "n_sampled_writes": int(run.n_sampled_writes),
        "disk_gb_estimate": run.n_writes * 2.0 * (fields + cloud) / 1024 ** 3,
        "disk_basis": "cases/cai2012-health bytes per cell and per parcel "
                      "(collisionless Kn = 100); every write kept, decomposed "
                      "and reconstructed",
        "memory": _memory(derived),
    }
    if cost_model is not None:
        seconds = cost_model.seconds(parcels=parcels, cells=plan.retained_cells,
                                     steps=run.n_steps)
        cost.update({"predicted_wall_hours": seconds / 3600.0,
                     "cost_model_origin": cost_origin,
                     "cost_caveat": "fitted to collisionless runs; no collision "
                                    "work in it -- a lower bound, uncalibrated "
                                    "for interface cases"})
    return cost
