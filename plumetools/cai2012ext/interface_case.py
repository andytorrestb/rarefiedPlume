r"""One case in interface mode: derive it, assess it, mesh it, write its inflow.

The mode is chosen in exactly one place, :func:`plumetools.cai2012ext.interface.load_interface`
-- ``interface.yaml`` absent or ``enabled: false`` is the direct exit, unchanged
-- and an interface case is derived in exactly one place, :func:`derive`.
Generation, the manifest, ``Allmesh``, ``runInflow.py``, the PBS guard and the
NON-PRODUCTION label all call it, so the sphere they mesh, the state they
inject, the weight they write and the verdict they print are the same numbers.

.. code-block:: text

    python -m plumetools.cai2012ext.interface_case report       <case>
    python -m plumetools.cai2012ext.interface_case mesh-dicts   <case>
    python -m plumetools.cai2012ext.interface_case verify-mesh  <case>
    python -m plumetools.cai2012ext.interface_case inflow       <case>
    python -m plumetools.cai2012ext.interface_case check-decomposed <case>
    python -m plumetools.cai2012ext.interface_case collect      <case>

The verdicts
------------
==============================  ==============================================
``interface-infeasible``        no admissible sphere, or the configuration is
                                missing a decision; nothing can be meshed
``unsuitable-for-production``   a production case misses a numerical criterion
``pending-coupling-validation`` numerically sound, but no overlap evidence at
                                this threshold has passed yet
``production-candidate``        every criterion met AND the coupling validated
``interface-validation``        a truncated Kn = 0.01 run for the overlap study
``development-smoke``           a mechanics check; never evidence
==============================  ==============================================

Numerical suitability (``numerical.status``) is reported separately from the
coupling validation (``coupling_validation.status``), and neither implies the
other.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
import sys
from collections import deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

from plumetools import geometry as facegeom
from plumetools.cai2012 import dictionaries, foamfields
from plumetools.cai2012 import inflow as cai_inflow
from plumetools.cai2012.config import CaiCaseConfig, load_case_config
from plumetools.cai2012.geometry import CaiGeometry
from plumetools.cai2012.geometry import from_config as cai_geometry
from plumetools.cai2012.inflow import ExitState, RunSettings
from plumetools.cai2012ext import feasibility, interface, interface_mesh, policy
from plumetools.foamio import volfield
from plumetools.foamio.patchfields import nonuniform_value_entry
from plumetools.inflow import read_patch_geometry
from plumetools.mesh.boundary import read_boundary

STATUS_INFEASIBLE = "interface-infeasible"
STATUS_UNSUITABLE = feasibility.STATUS_UNSUITABLE
STATUS_PENDING = "pending-coupling-validation"
STATUS_CANDIDATE = feasibility.STATUS_CANDIDATE
STATUS_VALIDATION = "interface-validation"
STATUS_SMOKE = "development-smoke"

#: Written beside ``case.yaml`` by ``Allmesh``; checked before every run.
IDENTITY_FILE = "interface-identity.yaml"

#: The per-case report, written by ``Allmesh`` and ``runInflow.py``.
REPORT_FILE = Path("results") / "interface-report.yaml"

#: Excluded Kn an interface-VALIDATION case may sit at: Cai's 0.01, where the
#: full-domain benchmark exists. 0.1 and 100 stay refused for every purpose.
VALIDATION_KNUDSEN = (0.01,)

#: The parallel tree, beside ``Cases/``.
TREE = "Interface"
MANIFEST = "interface-manifest.yaml"

#: Fraction of R a snapped vertex may sit off the sphere.
SNAP_TOLERANCE = 0.02

#: Bytes per cell and per parcel: NOT calibrated for this solver at these
#: densities. Only ``resources`` in ``interface.yaml`` sets them.
RESOURCE_KEYS = ("bytes_per_cell", "bytes_per_parcel", "memory_ceiling_bytes",
                 "memory_ceiling_pbs")


# --------------------------------------------------------------------------- #
# derivation
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class _RetainedState:
    """What :func:`plumetools.cai2012.inflow.derive_run_settings` reads of a state.

    ``velocity_m_per_s`` stays the EXIT speed: one domain transit is
    ``x_max / U_e`` for both modes, so the averaging window means the same
    thing as in the direct cases. The density and the characteristic speed are
    the declared reference cell's.
    """

    number_density_per_m3: float
    characteristic_speed_m_per_s: float
    velocity_m_per_s: float


@dataclass(frozen=True)
class InterfacePlan:
    """What a meshable interface case runs at -- all of it, or none of it.

    Attributes:
        mesh_plan: the lattice and snappyHexMesh plan around the sphere.
        run: time step, averaging window and particle weight, set at the
            declared reference cell rather than the exit cell.
        reference: the reference cell's state: ``identity``,
            ``number_density_per_m3``, ``characteristic_speed_m_per_s``,
            ``collision_time_s``.
        parcels: estimated parcels in the retained box at steady state.
        ideal_flux: the flux balance on the exact (unmeshed) sphere.
        resolution: planned cell / mean free path along the axis and the core
            edge, and the cells the criterion would need.
    """

    mesh_plan: interface_mesh.InterfaceMeshPlan
    run: RunSettings
    reference: dict
    parcels: float
    ideal_flux: dict
    resolution: dict


@dataclass(frozen=True)
class Derived:
    """Everything an interface case runs at. Built only by :func:`derive`.

    ``planned`` is set exactly when an admissible sphere was found and its
    mesh planned, and then holds every planned quantity. When it is ``None``,
    ``selection`` (no admissible sphere) or ``mesh_error`` (planning failed)
    says why; nothing downstream has to infer which from the fields that are
    missing.

    ``mesh_plan``, ``run``, ``reference``, ``parcels``, ``ideal_flux`` and
    ``resolution`` read through to ``planned`` -- ``None`` when it is -- for
    code written against the earlier flat layout.
    """

    case_dir: Path
    cfg: CaiCaseConfig
    geom: CaiGeometry
    exit_state: ExitState
    settings: interface.InterfaceConfig
    raw: dict
    source: interface.EquivalentSource
    selection: interface.SphereSelection
    planned: InterfacePlan | None = None
    mesh_error: str | None = None

    @property
    def meshable(self) -> bool:
        return self.planned is not None

    @property
    def mesh_plan(self) -> interface_mesh.InterfaceMeshPlan | None:
        return None if self.planned is None else self.planned.mesh_plan

    @property
    def run(self) -> RunSettings | None:
        return None if self.planned is None else self.planned.run

    @property
    def reference(self) -> dict | None:
        return None if self.planned is None else self.planned.reference

    @property
    def parcels(self) -> float | None:
        return None if self.planned is None else self.planned.parcels

    @property
    def ideal_flux(self) -> dict | None:
        return None if self.planned is None else self.planned.ideal_flux

    @property
    def resolution(self) -> dict | None:
        return None if self.planned is None else self.planned.resolution


def estimate_parcels(source, settings, geom, radius_m: float, weight: float,
                     n_theta: int = 90, n_phi: int = 72) -> float:
    """Parcels in the retained box at steady state, from the closure.

    The number between ``R`` and the box along a ray is the ray's number flow
    times its transit time ``int dr / U(r)``, so the total needs no volume
    grid. Core rays only (the tail is not injected); collisionless spreading
    out of the cone moves parcels but does not create them.
    """
    theta_c = settings.core_half_angle_rad
    if settings.tails == "inject_exploratory":
        theta_c = source.theta_L
    t_edges = np.linspace(0.0, theta_c, n_theta + 1)
    t = 0.5 * (t_edges[:-1] + t_edges[1:])
    p = (np.arange(n_phi) + 0.5) * 2.0 * math.pi / n_phi
    tt, pp = np.meshgrid(t, p, indexing="ij")
    direction = np.stack([np.cos(tt), np.sin(tt) * np.cos(pp),
                          np.sin(tt) * np.sin(pp)], axis=-1)
    with np.errstate(divide="ignore"):
        exits = np.stack([
            np.where(direction[..., 0] > 0, geom.x_max_m / direction[..., 0], np.inf),
            np.where(direction[..., 1] > 0, geom.y_max_m / direction[..., 1],
                     np.where(direction[..., 1] < 0, geom.y_min_m / direction[..., 1],
                              np.inf)),
            np.where(direction[..., 2] > 0, geom.z_max_m / direction[..., 2],
                     np.where(direction[..., 2] < 0, geom.z_min_m / direction[..., 2],
                              np.inf)),
        ]).min(axis=0)
    solid = (np.cos(t_edges[:-1]) - np.cos(t_edges[1:]))[:, None] * (2.0 * math.pi / n_phi)
    flow = source.j(tt) * solid / source.mass_kg            # molecules/s per ray
    r_far = float(np.max(exits))
    grid = np.geomspace(radius_m, max(r_far, radius_m * 1.0001), 400)
    inverse_u = 1.0 / source.state(grid, np.zeros_like(grid))["U"]
    cumulative = np.concatenate([[0.0], np.cumsum(
        0.5 * (inverse_u[1:] + inverse_u[:-1]) * np.diff(grid))])
    transit = np.interp(np.clip(exits, radius_m, r_far), grid, cumulative)
    return float(np.sum(flow * transit) / weight)


def derive(case_dir: Path) -> Derived:
    """Derive an interface case from its own ``case.yaml`` and ``interface.yaml``.

    Raises:
        policy.ExtensionError: if the case is not in interface mode.
        interface.InterfaceConfigError: for a malformed ``interface.yaml``.
    """
    case_dir = Path(case_dir)
    raw = interface.read_settings(case_dir)
    if raw is None:
        raise policy.ExtensionError(f"{case_dir}: no {interface.SETTINGS_FILE}")
    settings = interface.parse_interface(raw["interface"],
                                         str(case_dir / interface.SETTINGS_FILE))
    if not settings.enabled:
        raise policy.ExtensionError(
            f"{case_dir}: {interface.SETTINGS_FILE} says enabled: false -- this "
            f"is a direct-exit case")
    _check_resources(raw.get("resources"), case_dir)
    cfg = load_case_config(case_dir)
    geom = cai_geometry(cfg)
    exit_state = cai_inflow.from_config(cfg)
    source = interface.equivalent_source(
        cfg, exponent_offset=settings.angular_exponent_offset)
    r_max, why = interface.domain_radius_limit(cfg, settings)
    selection = interface.select_sphere(source, settings, r_max_m=r_max,
                                        r_max_reason=why)
    base = dict(case_dir=case_dir, cfg=cfg, geom=geom, exit_state=exit_state,
                settings=settings, raw=raw, source=source, selection=selection)
    if not selection.admissible:
        return Derived(**base)

    radius = float(selection.radius_m)
    at_interface = source.state(np.array([radius]), np.array([0.0]))
    try:
        mesh_plan = interface_mesh.plan(
            cfg, geom, settings, radius_m=radius,
            mean_free_path_m=float(at_interface["lambda"][0]))
    except ValueError as exc:
        return Derived(**base, mesh_error=str(exc))

    reference = interface.reference_state(source, radius, mesh_plan.core_cell_size_m)
    retained = _RetainedState(
        number_density_per_m3=reference["number_density_per_m3"],
        characteristic_speed_m_per_s=reference["characteristic_speed_m_per_s"],
        velocity_m_per_s=float(exit_state.velocity_m_per_s))
    run = cai_inflow.derive_run_settings(
        cfg, retained, geom, min_cell_size_m=mesh_plan.min_cell_size_m,
        exit_cell_volume_m3=mesh_plan.core_cell_volume_m3)
    bound = feasibility.MAX_DT_OVER_COLLISION_TIME * reference["collision_time_s"]
    if cfg.dsmc.delta_t_s is None and run.delta_t_s > bound:
        pinned = dataclasses.replace(cfg, dsmc=dataclasses.replace(
            cfg.dsmc, delta_t_s=bound))
        run = dataclasses.replace(
            cai_inflow.derive_run_settings(
                pinned, retained, geom, min_cell_size_m=mesh_plan.min_cell_size_m,
                exit_cell_volume_m3=mesh_plan.core_cell_volume_m3),
            delta_t_source="interface reference-cell collision time")

    parcels = estimate_parcels(source, settings, geom, radius,
                               run.n_equivalent_particles)
    centres, normals, areas = interface.ideal_sphere_faces(radius, 0.5 * math.pi)
    ideal = interface.flux_balance(
        source, interface.face_inflow(source, settings, centres=centres,
                                      outward_normals=normals, areas=areas),
        selection.core_fraction if settings.tails == "exclude" else 1.0)
    resolution = interface_mesh.resolution_estimate(
        source, mesh_plan, geom, settings.core_half_angle_rad)
    if mesh_plan.coarsened:
        # What the criterion would need with no budget -- planned, never written.
        unbudgeted = dataclasses.replace(cfg, mesh=dataclasses.replace(
            cfg.mesh, max_cells=2 ** 62))
        try:
            resolution["cells_to_meet_criterion"] = interface_mesh.plan(
                unbudgeted, geom, settings, radius_m=radius,
                mean_free_path_m=float(at_interface["lambda"][0])).retained_cells
        except ValueError:
            resolution["cells_to_meet_criterion"] = None
    return Derived(**base, planned=InterfacePlan(
        mesh_plan=mesh_plan, run=run, reference=reference, parcels=parcels,
        ideal_flux=ideal, resolution=resolution))


def _check_resources(raw, case_dir) -> None:
    if raw is None:
        return
    if not isinstance(raw, dict):
        raise interface.InterfaceConfigError(
            f"{case_dir}: 'resources' must be a mapping")
    unknown = sorted(set(raw) - set(RESOURCE_KEYS))
    if unknown:
        raise interface.InterfaceConfigError(
            f"{case_dir}: unknown resources key(s) {unknown}; valid: "
            f"{list(RESOURCE_KEYS)}")
    for key in ("bytes_per_cell", "bytes_per_parcel", "memory_ceiling_bytes"):
        value = raw.get(key)
        if value is not None and (isinstance(value, bool)
                                  or not isinstance(value, (int, float))
                                  or not value > 0):
            raise interface.InterfaceConfigError(
                f"{case_dir}: resources.{key} must be a positive number or null")


# --------------------------------------------------------------------------- #
# assessment
# --------------------------------------------------------------------------- #

def _coupling_validation(derived: Derived) -> dict:
    """Whether overlap evidence at THIS threshold has passed."""
    settings = derived.settings
    if settings.purpose != "production":
        return {"required": False, "status": "not-applicable",
                "note": f"a {settings.purpose} run does not need it"}
    path = settings.validation_evidence
    if path is None:
        return {"required": True, "status": "pending",
                "note": "no overlap evidence named (interface.validation_evidence)"}
    evidence = (derived.case_dir / path).resolve()
    if not evidence.is_file():
        return {"required": True, "status": "missing", "evidence": str(evidence),
                "note": "the named overlap result is not on this machine"}
    try:
        doc = yaml.safe_load(evidence.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        return {"required": True, "status": "unreadable", "evidence": str(evidence),
                "note": str(exc)}
    verdict = str(doc.get("verdict", "missing"))
    threshold = (doc.get("interface") or {}).get("bird_threshold")
    same = threshold is not None and math.isclose(
        float(threshold), float(settings.bird_threshold), rel_tol=1e-12)
    status = "passed" if verdict == "pass" and same else (
        "threshold-mismatch" if verdict == "pass" else verdict)
    return {"required": True, "status": status, "evidence": str(evidence),
            "evidence_threshold": threshold, "verdict": verdict}


def assess(derived: Derived, *, cost_model=None, cost_origin=None) -> dict:
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
    coupling = _coupling_validation(derived)

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


def _memory(derived: Derived) -> dict:
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


def _cost(derived: Derived, cost_model, cost_origin) -> dict:
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


# --------------------------------------------------------------------------- #
# the report and the manifest row
# --------------------------------------------------------------------------- #

def report(derived: Derived, assessment: dict | None = None, *,
           decision_table: bool = False) -> dict:
    """Everything about the case, for ``results/interface-report.yaml``."""
    assessment = assessment or assess(derived)
    source, selection, settings = derived.source, derived.selection, derived.settings
    d = source.diameter_m
    out = {
        "model": interface.MODEL_NAME,
        "purpose": settings.purpose,
        "non_production": assessment["status"] != STATUS_CANDIDATE,
        "settings": settings.as_dict(),
        "source": source.as_dict(),
        "selection": selection.as_dict(d, source.r_a),
        "diagnostic_radius_over_D": settings.diagnostic_radius_over_D,
        "diagnostic_radius_reproduces_cai": settings.diagnostic_radius_over_D == 10.0,
        "assessment": assessment,
        "terminal_model_comparison": interface.compare_with_terminal_model(
            source, sorted({1.0, 2.0, 5.0, 10.0} | (
                {selection.radius_m / d} if selection.admissible else set()))),
    }
    if decision_table:
        out["angular_shape_comparison"] = interface.angular_shape_comparison(source)
    if derived.meshable:
        out.update({
            "mesh": derived.planned.mesh_plan.as_dict(),
            "run": derived.planned.run.as_dict(),
            "reference_cell": derived.planned.reference,
            "estimated_parcels": derived.planned.parcels,
            "ideal_sphere_flux_balance": derived.planned.ideal_flux,
            "resolution_estimate": derived.planned.resolution,
        })
    if decision_table:
        r_max, _ = interface.domain_radius_limit(derived.cfg, settings)
        r_valid = settings.source_validity_min_radius_over_D or (source.r_a / d)
        out["decision_table"] = {
            "note": "what each (core angle, threshold) WOULD give at "
                    f"source_validity_min_radius_over_D = {r_valid:.4g}; for "
                    "choosing them before any DSMC result, not against one",
            "rows": interface.radius_table(
                source, settings, angles_deg=(10.0, 20.0, 30.0, 45.0, 60.0),
                thresholds=sorted({0.02, 0.03, 0.05, settings.bird_threshold}),
                r_min_over_D=r_valid, r_max_m=r_max),
        }
    return out


def manifest_entry(root: Path, name: str, derived: Derived, assessment: dict,
                   *, group: str) -> dict:
    """One ``interface-manifest.yaml`` row."""
    selection, source, settings = derived.selection, derived.source, derived.settings
    d = source.diameter_m
    entry = {
        "name": name,
        "path": derived.case_dir.resolve().relative_to(Path(root).resolve()).as_posix(),
        "group": group,
        "mode": "interface",
        "model": interface.MODEL_NAME,
        "Kn": float(derived.cfg.exit.knudsen),
        "purpose": settings.purpose,
        "run_mode": settings.run_mode,
        "bird_threshold": float(settings.bird_threshold),
        "bird_definition": settings.bird_definition,
        "core_half_angle_deg": settings.core_half_angle_deg,
        "selection_status": selection.status,
        "R_interface_over_D": None if selection.radius_m is None
        else selection.radius_m / d,
        "r_a_over_D": source.r_a / d,
        "max_core_bird_at_r_min": selection.max_core_bird_at_r_min,
        "axis_bird_at_r_min": selection.axis_bird_at_r_min,
        "core_mass_fraction": selection.core_fraction,
        "status": assessment["status"],
        "numerical_status": assessment["numerical"]["status"],
        "coupling_validation": assessment["coupling_validation"]["status"],
        "runnable": assessment["runnable"],
        "unresolved": assessment["unresolved"],
        "exit": {"D_m": d, "S0": float(derived.cfg.exit.speed_ratio),
                 "T0_K": float(derived.cfg.exit.T0_K), "n0_per_m3": source.n_e,
                 "U0_m_per_s": source.U_e,
                 "mdot_target_kg_per_s": source.targets.mdot_kg_per_s,
                 "h_target_J_per_kg": source.targets.h_J_per_kg},
    }
    if derived.meshable:
        entry.update({
            "n_cells": int(derived.planned.mesh_plan.retained_cells),
            "core_cell_size_m": derived.planned.mesh_plan.core_cell_size_m,
            "cell_over_interface_mfp": derived.planned.mesh_plan.cell_over_mfp,
            "deltaT_s": float(derived.planned.run.delta_t_s),
            "n_equivalent_particles": float(derived.planned.run.n_equivalent_particles),
            "average_start_s": float(derived.planned.run.average_start_s),
            "end_time_s": float(derived.planned.run.end_time_s),
            "n_steps": int(derived.planned.run.n_steps),
            "estimated_parcels": float(derived.planned.parcels),
        })
    return entry


def describe(derived: Derived, assessment: dict) -> list[str]:
    """Report lines for the console."""
    s, src = derived.selection, derived.source
    d = src.diameter_m
    lines = [
        f"          interface {s.status}: {s.reason}",
        f"          source    r_a = {src.r_a / d:.4f} D, M_a = {src.M_a:.4f}, "
        f"T_t = {src.T_t:.2f} K, V_inf = {src.V_inf:.2f} m/s, core fraction "
        f"{s.core_fraction if s.core_fraction is None else round(s.core_fraction, 4)}",
    ]
    if derived.meshable:
        p, r = derived.planned.mesh_plan, derived.planned.run
        lines += [
            f"          mesh      R = {p.radius_m / d:.4f} D, {p.retained_cells:,} "
            f"cells (lattice {p.lattice.n_cells:,}), core {p.core_cell_size_m:.4e} m, "
            f"cell/lambda {p.cell_over_mfp:.3g}",
            f"          step      deltaT {r.delta_t_s:.4e} s ({r.delta_t_source}), "
            f"Courant {r.courant:.3f}, weight {r.n_equivalent_particles:.4e}, "
            f"~{derived.planned.parcels:.2e} parcels",
        ]
    lines.append(f"          verdict   {assessment['status'].upper()}"
                 f" (numerical: {assessment['numerical']['status']}, coupling: "
                 f"{assessment['coupling_validation']['status']})")
    for reason in assessment["unresolved"]:
        lines.append(f"          UNRESOLVED {reason}")
    return lines


def write_report(case_dir: Path, document: dict) -> Path:
    path = Path(case_dir) / REPORT_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("# GENERATED by plumetools.cai2012ext.interface_case -- do "
                    "not edit.\n" + yaml.safe_dump(_plain(document), sort_keys=False,
                                                   default_flow_style=False),
                    encoding="utf-8", newline="\n")
    return path


def _plain(value):
    """numpy scalars and tuples to plain YAML types."""
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


# --------------------------------------------------------------------------- #
# identity: never resume old parcels on a different mesh or source
# --------------------------------------------------------------------------- #

def identity(derived: Derived) -> dict:
    """What a mesh and its results were built for."""
    case = derived.case_dir
    physics = dict(derived.raw.get("interface") or {})
    physics.pop("validation_evidence", None)
    physics.pop("label", None)
    canonical = json.dumps(physics, sort_keys=True, default=str)
    return {
        "mode": "interface",
        "model": interface.MODEL_NAME,
        "case_yaml_sha256": hashlib.sha256(
            (case / "case.yaml").read_bytes()).hexdigest(),
        "interface_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
        "R_interface_m": float(derived.selection.radius_m),
        "core_cell_size_m": float(derived.planned.mesh_plan.core_cell_size_m),
        "divisions": [int(v) for v in derived.planned.mesh_plan.lattice.divisions],
        "n_equivalent_particles": float(derived.planned.run.n_equivalent_particles),
    }


def write_identity(derived: Derived) -> Path:
    path = derived.case_dir / IDENTITY_FILE
    path.write_text("# GENERATED by Allmesh (interface mode): what this mesh was "
                    "built for.\n# A run whose results disagree with it is "
                    "refused.\n" + yaml.safe_dump(identity(derived), sort_keys=True),
                    encoding="utf-8", newline="\n")
    return path


def identity_problem(derived: Derived) -> str | None:
    """Why the case's mesh/results are not this configuration's, or ``None``."""
    path = derived.case_dir / IDENTITY_FILE
    if not path.is_file():
        return (f"no {IDENTITY_FILE}: the mesh here was not built by the "
                f"interface Allmesh (a direct-exit mesh, or none). Old parcels are "
                f"never resumed on a different mesh: ./Allclean --mesh && ./Allmesh.")
    stored = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    want = identity(derived)
    differ = sorted(k for k in want if stored.get(k) != want[k])
    if differ:
        return (f"{IDENTITY_FILE} disagrees with the current configuration in "
                f"{differ}: the mesh (and any parcels on it) belong to another "
                f"interface, source or weight. Start over: ./Allclean --mesh && "
                f"./Allmesh. Old parcels are never resumed on a different mesh.")
    return None


def knudsen_problem(cfg, settings, where: str) -> str | None:
    """The Kn rule for interface cases: Cai's 0.01 for validation only."""
    kn = cfg.exit.knudsen
    if not policy.is_excluded(kn):
        return None
    if settings is not None and settings.purpose == "interface-validation" and any(
            policy.same_knudsen(kn, k) for k in VALIDATION_KNUDSEN):
        return None
    return (policy._excluded_message(kn, where)
            + " Only an interface-validation case may sit at Kn = 0.01, against "
              "the full-domain benchmark; nothing may duplicate a published "
              "case as production.")


def guard(case_dir: Path) -> str | None:
    """Why this interface case must not be meshed or run, or ``None``.

    Shared by ``runMesh.sh`` and ``runCase.sh``: the Kn rule, an infeasible
    interface, and results that belong to another mesh. Whether a
    numerically sound case may run as PRODUCTION is decided by :func:`label`,
    which only ``runCase.sh`` calls.
    """
    from plumetools.cai2012ext.generation import has_solver_output

    case_dir = Path(case_dir)
    derived = derive(case_dir)
    problem = knudsen_problem(derived.cfg, derived.settings,
                              f"{case_dir.resolve().name}/case.yaml")
    if problem:
        return problem
    if not derived.meshable:
        return (f"{derived.selection.status}: "
                f"{derived.mesh_error or derived.selection.reason}")
    missing = derived.settings.missing_decisions() if derived.settings.scientific else []
    if missing:
        return f"{STATUS_INFEASIBLE}: " + " ".join(missing)
    if has_solver_output(case_dir):
        return identity_problem(derived)
    return None


def label_text(assessment: dict) -> str:
    return ("This interface case was run by its runCase.sh job.\n"
            f"Status: {assessment['status']} (purpose {assessment['purpose']}, "
            f"run_mode {assessment['run_mode']}). Nothing from it is a "
            "validated production result.\n\nUnresolved:\n"
            + "".join(f"  - {line}\n" for line in assessment["unresolved"]))


def label(case_dir: Path, marker_name: str) -> tuple[str, bool]:
    """Write the persistent NON-PRODUCTION marker; say whether it may run.

    Returns ``(text, runnable)``. A production case that is not a candidate
    is refused here unless ``run_mode: exploratory``; validation and smoke
    runs are always marked.
    """
    case_dir = Path(case_dir)
    assessment = assess(derive(case_dir))
    name = case_dir.resolve().name
    if assessment["status"] == STATUS_CANDIDATE:
        return (f"{name}: {STATUS_CANDIDATE} (interface) -- every criterion met "
                f"and the coupling validated; not yet a converged result\n", True)
    text = label_text(assessment)
    if not assessment["runnable"]:
        return (text + f"\nREFUSED: {assessment['refusal']}\n", False)
    marker = case_dir / marker_name
    marker.write_text(text, encoding="utf-8", newline="\n")
    return text + f"\nmarked {marker_name}: {marker.resolve()}\n", True


# --------------------------------------------------------------------------- #
# Allmesh
# --------------------------------------------------------------------------- #

def write_mesh_dictionaries(derived: Derived) -> list[Path]:
    """Every dictionary ``Allmesh`` writes, from the one derivation."""
    if not derived.meshable:
        raise policy.ExtensionError(
            f"nothing to mesh: {derived.selection.status} -- "
            f"{derived.mesh_error or derived.selection.reason}")
    case = derived.case_dir
    paths = interface_mesh.write_mesh_setup(case, derived.cfg, derived.geom,
                                            derived.planned.mesh_plan)
    paths += dictionaries.write_case_dictionaries(case, derived.cfg,
                                                  derived.exit_state, derived.planned.run)
    properties = case / "constant" / "dsmcProperties"
    text = properties.read_text(encoding="utf-8")
    note = ("// INTERFACE MODE: 'nozzle' is the downstream continuum sphere "
            f"R = {derived.selection.radius_m:.6g} m,\n"
            "// not the exit disk. nEquivalentParticles is set at the declared "
            "reference cell\n// (results/interface-report.yaml), not the exit "
            "cell the comments below describe.\n\n")
    marker = "// General Properties"
    if marker in text:
        properties.write_text(text.replace(marker, note + marker, 1),
                              encoding="utf-8", newline="\n")
    return paths


def openfoam_normals(faces, points) -> np.ndarray:
    """Unit area vectors by OpenFOAM's right-hand rule: out of the owner cell.

    For a boundary face the owner is the fluid, so these point OUT of the DSMC
    domain -- toward the source, on the interface cavity. Deliberately not
    :func:`plumetools.geometry.normals`, whose legacy ``(p0-p2) x (p0-p1)``
    orientation is the opposite one.
    """
    out = np.empty((len(faces), 3), dtype=np.float64)
    for i, face in enumerate(faces):
        v = points[np.asarray(face)]
        fan = np.cross(v[1:-1] - v[0], v[2:] - v[0]).sum(axis=0)
        out[i] = fan / np.linalg.norm(fan)
    return out


def verify_mesh(derived: Derived) -> dict:
    """The meshed inlet patch, measured; raises if it is not the sphere."""
    case = derived.case_dir
    name = derived.cfg.mesh.patch_names["nozzle"]
    patches = read_boundary(case)
    if name not in patches or patches[name].n_faces == 0:
        raise policy.ExtensionError(f"no '{name}' patch with faces in the mesh; "
                                    f"snappyHexMesh did not carve the sphere")
    if patches[name].type != "patch":
        raise policy.ExtensionError(
            f"'{name}' is a {patches[name].type!r}, not a 'patch': DSMC particles "
            f"crossing back through the interface must be deleted, not reflected")
    others = [n for n, info in patches.items() if n != name
              and info.type in ("wall", "symmetry", "symmetryPlane")]
    if others:
        raise policy.ExtensionError(
            f"patches {others} are walls or symmetry planes; an open free plume "
            f"has neither")
    labels, faces, points = read_patch_geometry(case, name)
    centres = facegeom.centroids(faces, points)
    patch_points = points[np.unique(np.concatenate([np.asarray(f) for f in faces]))]
    measured = interface_mesh.patch_report(
        centres=centres, outward_normals=openfoam_normals(faces, points),
        areas=facegeom.areas(faces, points), points=patch_points,
        radius_m=derived.selection.radius_m, diameter_m=derived.source.diameter_m,
        theta_core=derived.settings.core_half_angle_rad,
        theta_limit=derived.source.theta_L)
    measured["patches"] = {n: {"type": i.type, "faces": i.n_faces}
                           for n, i in patches.items()}
    measured["problems"] = interface_mesh.patch_problems(
        measured, snap_tolerance=SNAP_TOLERANCE,
        min_core_faces=int(derived.cfg.checks.min_nozzle_faces),
        area_tolerance=float(derived.cfg.checks.max_nozzle_area_error))
    return measured


# --------------------------------------------------------------------------- #
# runInflow.py
# --------------------------------------------------------------------------- #

def _nonuniform(patch: str, list_type: str, rows) -> list[str]:
    return nonuniform_value_entry(patch, list_type, rows)


def write_inflow(derived: Derived) -> tuple[list[Path], dict]:
    """Write ``0/`` from the meshed patch, check it, and return the balance.

    Raises:
        policy.ExtensionError: wrong mesh, a patch that is not the sphere, a
            written field that does not read back in patch order, or a
            scientific case whose meshed flux misses the core target by more
            than ``meshed_flux_tolerance``.
    """
    problem = identity_problem(derived)
    if problem:
        raise policy.ExtensionError(problem)
    measured = verify_mesh(derived)
    if measured["problems"]:
        raise policy.ExtensionError("the meshed interface patch is not the "
                                    "sphere: " + "; ".join(measured["problems"]))
    case, cfg = derived.case_dir, derived.cfg
    name = cfg.mesh.patch_names["nozzle"]
    labels, faces, points = read_patch_geometry(case, name)
    faces_inflow = interface.face_inflow(
        derived.source, derived.settings, centres=facegeom.centroids(faces, points),
        outward_normals=openfoam_normals(faces, points),
        areas=facegeom.areas(faces, points))
    fraction = (derived.selection.core_fraction
                if derived.settings.tails == "exclude" else 1.0)
    balance = interface.flux_balance(derived.source, faces_inflow, fraction)
    worst = max(abs(balance["gross_mass_residual"]),
                abs(balance["gross_energy_residual"]))
    balance["within_tolerance"] = worst <= derived.settings.meshed_flux_tolerance
    if not balance["within_tolerance"] and derived.settings.scientific:
        raise policy.ExtensionError(
            f"the meshed patch injects {balance['gross_mass_residual']:+.3%} mass "
            f"and {balance['gross_energy_residual']:+.3%} energy against the core "
            f"targets, beyond meshed_flux_tolerance "
            f"{derived.settings.meshed_flux_tolerance:g}. A failed coupling, not "
            f"something to rescale: refine the mesh to separate faceting from "
            f"the model.")

    patches = read_boundary(case)
    others = {n: info for n, info in patches.items() if n != name}
    out = case / "0"
    species = cfg.gas.species_name
    u_rows = [f"            ({u:.10g} {v:.10g} {w:.10g})"
              for u, v, w in faces_inflow.velocity]
    t_rows = [f"            {t:.10g}" for t in faces_inflow.temperature]
    n_rows = [f"            {n:.10g}" for n in faces_inflow.number_density]
    note = (f"INTERFACE MODE: the closure state on the sphere R = "
            f"{derived.selection.radius_m:.6g} m (patch '{name}'), per face in "
            f"patch order. Zero density beyond the {derived.settings.core_half_angle_deg:g} "
            f"deg core; T stays positive there for plumeFieldInflow.")
    written = [
        foamfields._write_field(
            out / "boundaryU", "volVectorField", "boundaryU", "[0 1 -1 0 0 0 0]",
            "(0 0 0)", _nonuniform(name, "vector", u_rows)
            + foamfields._boundary_entries(others, "", None), note),
        foamfields._write_field(
            out / "boundaryT", "volScalarField", "boundaryT", "[0 0 0 1 0 0 0]",
            "0", _nonuniform(name, "scalar", t_rows)
            + foamfields._boundary_entries(others, "", None), note),
        foamfields._write_field(
            out / f"boundaryNumberDensity_{species}", "volScalarField",
            f"boundaryNumberDensity_{species}", "[0 -3 0 0 0 0 0]", "0",
            _nonuniform(name, "scalar", n_rows)
            + foamfields._boundary_entries(others, "", None), note),
    ]
    written += foamfields.write_measurement_fields(case, patches)

    # Read back in patch order: the ordering contract with the solver.
    back_n = volfield.read_patch_field(out / f"boundaryNumberDensity_{species}", name)
    back_u = volfield.read_patch_field(out / "boundaryU", name)
    back_t = volfield.read_patch_field(out / "boundaryT", name)
    for label, got, want in (("n", back_n, faces_inflow.number_density),
                             ("U", back_u, faces_inflow.velocity),
                             ("T", back_t, faces_inflow.temperature)):
        if got.shape != np.shape(want) or not np.allclose(got, want, rtol=1e-9,
                                                          atol=0.0):
            raise policy.ExtensionError(
                f"0/ {label} on '{name}' does not read back as written, in patch "
                f"order ({got.shape} vs {np.shape(want)})")
    if np.any(back_t <= 0.0):
        raise policy.ExtensionError("a non-positive boundaryT on the interface")
    reread = interface.flux_balance(derived.source, dataclasses.replace(
        faces_inflow, number_density=back_n, velocity=back_u, temperature=back_t,
        u_normal=np.einsum("ij,ij->i", back_u, faces_inflow.inward_normals)),
        fraction)
    balance["from_written_fields_gross_mass_kg_per_s"] = \
        reread["gross_injected_mass_kg_per_s"]
    return written, {"patch": measured, "flux_balance": balance}


def check_decomposed(derived: Derived) -> dict:
    """The decomposed ``0/`` fields carry the same injection as the whole.

    Sums ``A m n J1`` over every ``processor*/`` share of the interface patch
    and compares it with the reconstructed patch's total.
    """
    case, name = derived.case_dir, derived.cfg.mesh.patch_names["nozzle"]
    species = derived.cfg.gas.species_name

    def total(directory: Path) -> tuple[float, int]:
        patches = read_boundary(directory)
        if name not in patches or patches[name].n_faces == 0:
            return 0.0, 0
        _, faces, points = read_patch_geometry(directory, name)
        n = volfield.read_patch_field(
            directory / "0" / f"boundaryNumberDensity_{species}", name)
        u = volfield.read_patch_field(directory / "0" / "boundaryU", name)
        t = volfield.read_patch_field(directory / "0" / "boundaryT", name)
        normals = openfoam_normals(faces, points)
        areas = facegeom.areas(faces, points)
        j1, _ = interface.one_way_moments(np.einsum("ij,ij->i", u, -normals), t,
                                          derived.source.gas_constant)
        return float(np.sum(areas * derived.source.mass_kg * n * j1)), len(areas)

    whole, faces = total(case)
    parts = sorted(case.glob("processor[0-9]*"))
    if not parts:
        return {"decomposed": False}
    pieces = [total(p) for p in parts]
    summed = sum(m for m, _ in pieces)
    return {"decomposed": True, "processors": len(parts),
            "faces_reconstructed": faces, "faces_decomposed": sum(f for _, f in pieces),
            "gross_mass_reconstructed": whole, "gross_mass_decomposed": summed,
            "relative_difference": summed / whole - 1.0 if whole else None}


# --------------------------------------------------------------------------- #
# the per-case post-processing, in either mode
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class PostContext:
    """What ``Allcompare.py`` / ``Allconverge.py`` need of a case, either mode.

    ``refusal`` is set when the case must not be post-processed (an excluded
    Kn outside the interface-validation exception). ``interface`` is ``None``
    for a direct-exit case and otherwise the mask and provenance every
    comparison carries.
    """

    mode: str
    cfg: object
    geom: object
    exit_state: object
    run: object
    verdict: dict
    nozzle_area_error: float | None
    interface: dict | None
    refusal: str | None


def post_context(case_dir: Path) -> PostContext:
    """The direct-exit derivation unchanged, or the interface one."""
    from plumetools.cai2012ext import direct_case

    case_dir = Path(case_dir)
    if interface.load_interface(case_dir) is None:
        case = direct_case.derive(case_dir)
        verdict = feasibility.assess_case(case)
        refusal = policy.direct_exit_problem(case.cfg, str(case_dir))
        return PostContext("direct", case.cfg, case.geom, case.exit_state, case.run,
                           verdict, float(case.nozzle.area_error), None, refusal)
    derived = derive(case_dir)
    assessment = assess(derived)
    refusal = knudsen_problem(derived.cfg, derived.settings, case_dir.name)
    if not derived.meshable:
        refusal = refusal or f"{assessment['status']}: nothing was meshed"
    d = derived.source.diameter_m
    radius = derived.selection.radius_m
    provenance = {
        "model": interface.MODEL_NAME,
        "purpose": derived.settings.purpose,
        "bird_threshold": derived.settings.bird_threshold,
        "bird_definition": derived.settings.bird_definition,
        "core_half_angle_deg": derived.settings.core_half_angle_deg,
        "R_interface_over_D": None if radius is None else radius / d,
        "source_point_m": list(derived.source.source_point),
        "exit_origin_m": list(derived.source.exit_origin),
        "core_mass_fraction": derived.selection.core_fraction,
        "mask": ("R < R_interface (exit-centred) is not in this mesh: no cell "
                 "there is compared, and none is scored as zero. Cells beyond "
                 "the core angle lack the excluded tail by construction."),
        "diagnostic_radius_over_D": derived.settings.diagnostic_radius_over_D,
        "reference_note": "Cai's solution is still launched from the PHYSICAL "
                          "exit; the interface changes what is simulated, not "
                          "what it is compared with",
    }
    verdict = {"status": assessment["status"],
               "production_ready": assessment["production_ready"],
               "unresolved": assessment["unresolved"]}
    return PostContext("interface", derived.cfg, derived.geom, derived.exit_state,
                       derived.run, verdict, None, provenance, refusal)


# --------------------------------------------------------------------------- #
# collect: a small file to bring back from the cluster
# --------------------------------------------------------------------------- #

def collect(case_dir: Path) -> dict:
    """Concise evidence from a run: verdicts, metrics, resources, log tails.

    Reads only small files the case already wrote -- never a field -- so it
    can be run on the cluster and the output copied to the laptop.
    """
    case_dir = Path(case_dir)
    out = {"case": case_dir.resolve().name}
    for key, relative in (("interface_report", REPORT_FILE),
                          ("inflow", Path("results") / "interface-inflow.yaml"),
                          ("overlap", Path("results") / "overlap.yaml"),
                          ("comparison", Path("results") / "comparison.yaml"),
                          ("diagnostics", Path("results") / "diagnostics_summary.yaml"),
                          ("stationarity", Path("results") / "stationarity.yaml")):
        path = case_dir / relative
        if path.is_file():
            try:
                doc = yaml.safe_load(path.read_text(encoding="utf-8"))
            except yaml.YAMLError:
                doc = {"unreadable": str(path)}
            out[key] = _summarise(key, doc)
    marker = case_dir / "NON-PRODUCTION"
    out["non_production_marker"] = marker.read_text(encoding="utf-8") \
        if marker.is_file() else None
    for log in ("log.dsmcFoam", "log.snappyHexMesh", "log.checkMesh",
                "output.log", "error.log"):
        path = case_dir / log
        if path.is_file():
            out.setdefault("logs", {})[log] = log_summary(path)
    return out


def log_summary(path: Path) -> dict:
    """A log's line count, last 15 lines, last 5 memory lines, last ExecutionTime.

    One streamed pass holding only those few lines: a ``log.dsmcFoam`` from a
    long run is gigabytes, and reading it whole to keep twenty lines of it
    needed several times its size in memory. Lines are split exactly as
    ``str.splitlines`` splits the whole text.
    """
    count = 0
    tail, memory, execution = deque(maxlen=15), deque(maxlen=5), deque(maxlen=1)
    with open(path, encoding="utf-8", errors="replace", newline="") as handle:
        for physical in handle:
            for line in physical.splitlines():
                count += 1
                tail.append(line)
                if "MB" in line or "memory" in line.lower():
                    memory.append(line)
                if line.startswith("ExecutionTime"):
                    execution.append(line)
    return {"lines": count, "tail": list(tail), "memory_lines": list(memory),
            "execution_time": list(execution)}


def _summarise(key: str, doc):
    if not isinstance(doc, dict):
        return doc
    if key == "interface_report":
        keep = ("purpose", "non_production", "selection", "assessment",
                "reference_cell", "estimated_parcels", "ideal_sphere_flux_balance")
        summary = {k: doc[k] for k in keep if k in doc}
        if "assessment" in summary:
            summary["assessment"] = {k: summary["assessment"].get(k) for k in
                                     ("status", "runnable", "unresolved",
                                      "numerical", "coupling_validation")}
        if "selection" in summary:
            summary["selection"].pop("scan", None)
        return summary
    return doc


# --------------------------------------------------------------------------- #
# the command line
# --------------------------------------------------------------------------- #

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m plumetools.cai2012ext.interface_case",
        description="Derive, mesh and inject one interface-mode case.")
    steps = parser.add_subparsers(dest="step", required=True)
    for name in ("report", "mesh-dicts", "verify-mesh", "inflow",
                 "check-decomposed", "collect"):
        step = steps.add_parser(name)
        step.add_argument("case_dir", type=Path)
        if name == "report":
            step.add_argument("--decision-table", action="store_true")
        if name == "collect":
            step.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)
    case = args.case_dir

    if args.step == "collect":
        document = _plain(collect(case))
        text = yaml.safe_dump(document, sort_keys=False, default_flow_style=False)
        if args.output:
            args.output.write_text(text, encoding="utf-8", newline="\n")
            print(f"wrote {args.output}")
        else:
            print(text, end="")
        return 0

    try:
        derived = derive(case)
    except (policy.ExtensionError, interface.InterfaceConfigError, ValueError) as exc:
        print(f"interface_case {args.step}: {exc}", file=sys.stderr)
        return 1
    assessment = assess(derived)
    print("\n".join(describe(derived, assessment)))

    if args.step == "report":
        path = write_report(case, report(derived, assessment,
                                         decision_table=args.decision_table))
        print(f"    wrote {path}")
        return 0 if derived.meshable else 3
    if args.step == "mesh-dicts":
        try:
            paths = write_mesh_dictionaries(derived)
        except policy.ExtensionError as exc:
            print(f"interface_case mesh-dicts: REFUSED. {exc}", file=sys.stderr)
            write_report(case, report(derived, assessment, decision_table=True))
            return 3
        print("\n".join(derived.planned.mesh_plan.describe()))
        for path in paths + [write_identity(derived),
                             write_report(case, report(derived, assessment))]:
            print(f"    wrote {path}")
        return 0
    if args.step == "verify-mesh":
        measured = verify_mesh(derived)
        print(yaml.safe_dump(_plain(measured), sort_keys=False))
        return 1 if measured["problems"] else 0
    if args.step == "inflow":
        try:
            written, document = write_inflow(derived)
        except policy.ExtensionError as exc:
            print(f"interface_case inflow: REFUSED. {exc}", file=sys.stderr)
            return 1
        path = case / "results" / "interface-inflow.yaml"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump(_plain(document), sort_keys=False),
                        encoding="utf-8", newline="\n")
        balance = document["flux_balance"]
        print(f"  inlet faces              {balance['faces']} "
              f"({balance['injecting_faces']} injecting)")
        print(f"  gross mass vs core target {balance['gross_mass_residual']:+.3e}"
              f"   energy {balance['gross_energy_residual']:+.3e}")
        print(f"  net mass vs core target   {balance['net_mass_residual']:+.3e}"
              f"   (Maxwellian return half "
              f"{balance['maxwellian_return_mass_kg_per_s']:.3e} kg/s)")
        for p in written + [path]:
            print(f"    wrote {p}")
        return 0
    if args.step == "check-decomposed":
        result = check_decomposed(derived)
        print(yaml.safe_dump(_plain(result), sort_keys=False))
        diff = result.get("relative_difference")
        return 0 if diff is None or abs(diff) < 1e-9 else 1
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
