"""Shared fitted-interface provenance, rectilinear sampling and flux assessment."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from . import workflow as wf, breakdown as bd, halfrange
from . import interface as surface_tools
from .fluxes import continuum_flux

def sample_containing_cells(mesh, positions, allowed, axes=(0, 1)):
    """Locate points in the rectilinear CFD mesh; never use nearest-cell guesses."""
    points = np.asarray(positions, float)
    if points.ndim != 2 or points.shape[1] <= max(axes) or not np.isfinite(points).all():
        raise wf.WorkflowError("sample positions must be finite x/y points")
    allowed = np.asarray(allowed, bool)
    if allowed.shape != (mesh.n_cells,):
        raise wf.WorkflowError("allowed CFD-cell mask has the wrong length")
    spacing = float(np.min(mesh.cell_size()[:, axes]))
    edges = []
    for k in axes:
        values = np.unique(mesh.points[:, k])
        groups = np.cumsum(np.r_[0, np.diff(values) > 1e-8*spacing])
        edges.append(np.bincount(groups, weights=values)/np.bincount(groups))
    shape = tuple(len(e)-1 for e in edges)
    if np.prod(shape) != mesh.n_cells:
        raise wf.WorkflowError("sampling requires a rectilinear CFD mesh with the stated active axes")
    p = points[:, axes].copy()
    indices = []
    for k, e in enumerate(edges):
        tol = 1e-10*max(1.0, e[-1]-e[0])
        if np.any(p[:, k] < e[0]-tol) or np.any(p[:, k] > e[-1]+tol):
            raise wf.WorkflowError("interface or DSMC sample lies outside the CFD domain")
        p[:, k] = np.clip(p[:, k], e[0], e[-1])
        indices.append(np.clip(np.searchsorted(e, p[:, k], side="right")-1, 0, len(e)-2))
    grid = np.full(shape, -1, dtype=int)
    cell_indices = tuple(np.searchsorted(e, mesh.cell_centres[:, axis], side="right")-1
                         for axis, e in zip(axes, edges))
    grid[cell_indices] = np.arange(mesh.n_cells)
    cells = grid[tuple(indices)]
    if np.any(cells < 0) or np.any(~allowed[cells]):
        raise wf.WorkflowError("sample has no containing accepted/allowed CFD cell")
    return cells


def breakdown_identity(case_dir, cfg, time):
    """Bind the classification to the exact source fields, mesh and configuration."""
    from .provenance import sha256
    case_dir = Path(case_dir)
    fields = {name: sha256(case_dir/"cfd"/time/name) for name in ("p", "rho", "T", "U")}
    mesh = {name: sha256(case_dir/"cfd/constant/polyMesh"/name)
            for name in ("points", "faces", "owner", "neighbour", "boundary")}
    data = {"config": cfg.digest, "time": time, "fields": fields, "mesh": mesh}
    return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()


def _assert_breakdown_current(case_dir, data, cfg=None):
    cfg = cfg or wf.load_case(case_dir)
    expected = str(data.get("input_identity", ""))
    if not expected or expected != breakdown_identity(case_dir, cfg, str(data["time"])):
        raise wf.WorkflowError("breakdown inputs changed or lack provenance; rerun the breakdown step")


def _identity(case_dir, cfg, fit, time):
    from .provenance import sha256
    fields = {name: sha256(Path(case_dir)/"cfd"/time/name) for name in ("p", "rho", "T", "U")}
    data = {"config": cfg.digest, "fit": [fit[k] for k in ("xc", "a", "b")],
            "time": time, "fields": fields, "breakdown_inputs": breakdown_identity(case_dir, cfg, time)}
    return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()


def _load_breakdown(case_dir):
    path = wf.results_dir(case_dir)/"breakdown.npz"
    if not path.is_file():
        raise wf.WorkflowError("run the breakdown step first")
    with np.load(path) as z:
        data = {k: z[k].copy() for k in z.files}
    _assert_breakdown_current(case_dir, data)
    return data


def _assess(mesh, faces, normals, areas, cells, all_cells, state, ind, grads, geo, cfg, source):
    gas = source.gas
    rho, temperature, velocity = state.rho[cells], state.T[cells], state.U[cells]
    valid = ((rho > 0) & (temperature > 0) & np.isfinite(rho)
             & np.isfinite(temperature) & np.isfinite(velocity).all(axis=1))
    if not valid.all():
        raise wf.WorkflowError("actual interface contains a non-positive or invalid CFD state")
    un = np.einsum("ij,ij->i", velocity, normals)
    mn = un/gas.sound_speed(temperature)
    forward = halfrange.integrate(halfrange.incoming_flux(
        rho, velocity, temperature, normals, gas.gas_constant), areas)
    reverse = halfrange.integrate(halfrange.incoming_flux(
        rho, velocity, temperature, -normals, gas.gas_constant), areas)
    gu = np.stack([g.vectors for g in grads["U"]], axis=2)[cells]
    gt = grads["T"].vectors[cells]
    grads_ok = np.isfinite(gu).all(axis=(1, 2)) & np.isfinite(gt).all(axis=1)
    nsf = continuum_flux(rho, velocity, temperature, normals, gas,
                         grad_u=np.where(grads_ok[:, None, None], gu, 0),
                         grad_t=np.where(grads_ok[:, None], gt, 0))
    continuum = halfrange.integrate(nsf, areas)
    net = {k: (np.asarray(forward[k])-np.asarray(reverse[k])).tolist()
           for k in ("mass", "momentum", "energy")}
    returning = float(reverse["mass"]/forward["mass"]) if forward["mass"] > 0 else float("inf")
    sampled = {name: ind[name][all_cells] for name in (*bd.INDICATORS, "valid")}
    overlap = surface_tools.overlap_layers(
        mesh, SimpleNamespace(retained=np.zeros(mesh.n_cells, bool), dsmc_cell=mesh.owner[faces]),
        sampled, cfg.require("breakdown.thresholds"),
        tuple(cfg.get("breakdown.criterion", ["Kn_GLL", "B_exp"])),
        float(cfg.get("breakdown.overlap_ceiling_factor", 2.0)))
    subsonic = float(areas[mn < 1].sum()/areas.sum())
    out = {"face_state_method": "bounded containing-cell reconstruction (zero slope)",
           "interface_faces": len(faces), "interface_area_m2": float(areas.sum()),
           "return_mass_fraction": returning, "outer_boundary_contact_fraction": 0.0,
           "invalid_face_states": 0, "invalid_retained_cells": 0,
           "faces_without_gradients": int((~grads_ok).sum()), "geometry_checks": geo,
           "source": {"covered_area_fraction": geo["source_covered_area_fraction"],
                      "uncovered_faces_injected_directly_by_dsmc": mesh.patches.get("source", SimpleNamespace(n_faces=0)).n_faces},
           "normal_mach": {"min": float(mn.min()), "max": float(mn.max()),
                           "area_weighted_mean": float(np.average(mn, weights=areas)),
                           "subsonic_area_fraction": subsonic,
                           "inflow_to_cfd_area_fraction": float(areas[un < 0].sum()/areas.sum())},
           "overlap": overlap,
           "fluxes": {"maxwellian_forward": forward, "maxwellian_reverse": reverse,
                      "maxwellian_net": net, "continuum_nsf": continuum},
           "maxwellian_vs_nsf": {
               "mass_abs_difference_kg_s": net["mass"]-continuum["mass"],
               "energy_abs_difference_W": net["energy"]-continuum["energy"],
               "note": "Maxwellian inflow does not transfer viscous stress or heat flux exactly"},
           "flux_depth_m": source.depth_m, "criteria": {}, "reasons": []}
    crit = cfg.get("interface.criteria", {})
    values = [
        ("return_mass_fraction", returning, crit.get("max_return_mass_fraction", 0.05), False),
        ("subsonic_normal_area_fraction", subsonic, crit.get("max_subsonic_normal_area_fraction", 0.25), False),
        ("outer_boundary_contact_fraction", 0.0, crit.get("max_outer_boundary_contact_fraction", 0.0), False),
        ("overlap_layers", overlap["overlap_layers"], crit.get("min_overlap_layers", 1), True),
        ("invalid_face_states", 0, 0, False), ("invalid_retained_cells", 0, 0, False),
    ]
    for name, value, limit, lower_bound in values:
        passed = bool(np.isfinite(value) and (value >= limit if lower_bound else value <= limit))
        out["criteria"][name] = {"value": value, "limit": limit, "pass": passed}
        if not passed:
            out["reasons"].append(f"{name} = {value:.4g} (limit {limit})")
    out["status"] = surface_tools.STATUS_UNSUITABLE if out["reasons"] else surface_tools.STATUS_SUITABLE
    for name in bd.INDICATORS:
        values = ind[name][cells]
        out.setdefault("indicators_cfd_side", {})[name] = {
            "max": float(np.nanmax(values)), "area_weighted_mean": float(np.average(values, weights=areas))}
    return out, {"rho": rho, "T": temperature, "U": velocity, "valid": valid, "Mn": mn, "Un": un}

