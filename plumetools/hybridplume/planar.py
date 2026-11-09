"""Opt-in planar hybrid workflow with a separately meshed elliptical interface."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.spatial import cKDTree

from . import breakdown as bd, config, dsmc, halfrange, mapping, resources
from . import interface as surface_tools
from . import workflow as wf
from .ellipse import fit_ellipse, footprint_clearance
from .fluxes import continuum_flux
from .meshgeom import read_mesh
from .planar_mesh import render_dsmc_mesh


def _footprints(mesh):
    half = mesh.cell_size()[:, :2] / 2
    return mesh.cell_centres[:, :2] - half, mesh.cell_centres[:, :2] + half


def sample_containing_cells(mesh, positions, allowed):
    """Locate points in the rectilinear CFD mesh; never use nearest-cell guesses."""
    points = np.asarray(positions, float)
    if points.ndim != 2 or points.shape[1] not in (2, 3) or not np.isfinite(points).all():
        raise wf.WorkflowError("sample positions must be finite x/y points")
    allowed = np.asarray(allowed, bool)
    if allowed.shape != (mesh.n_cells,):
        raise wf.WorkflowError("allowed CFD-cell mask has the wrong length")
    spacing = float(np.min(mesh.cell_size()[:, :2]))
    edges = []
    for k in (0, 1):
        values = np.unique(mesh.points[:, k])
        groups = np.cumsum(np.r_[0, np.diff(values) > 1e-8*spacing])
        edges.append(np.bincount(groups, weights=values)/np.bincount(groups))
    nx, ny = (len(e)-1 for e in edges)
    if nx*ny != mesh.n_cells:
        raise wf.WorkflowError("planar sampling requires a one-layer rectilinear CFD mesh")
    p = points[:, :2].copy()
    indices = []
    for k, e in enumerate(edges):
        tol = 1e-10*max(1.0, e[-1]-e[0])
        if np.any(p[:, k] < e[0]-tol) or np.any(p[:, k] > e[-1]+tol):
            raise wf.WorkflowError("interface or DSMC sample lies outside the CFD domain")
        p[:, k] = np.clip(p[:, k], e[0], e[-1])
        indices.append(np.clip(np.searchsorted(e, p[:, k], side="right")-1, 0, len(e)-2))
    grid = np.full((nx, ny), -1, dtype=int)
    ix = np.searchsorted(edges[0], mesh.cell_centres[:, 0], side="right")-1
    iy = np.searchsorted(edges[1], mesh.cell_centres[:, 1], side="right")-1
    grid[ix, iy] = np.arange(mesh.n_cells)
    cells = grid[indices[0], indices[1]]
    if np.any(cells < 0) or np.any(~allowed[cells]):
        raise wf.WorkflowError("sample has no containing accepted/allowed CFD cell")
    return cells


def _identity(case_dir, cfg, fit, time):
    from .provenance import sha256
    fields = {name: sha256(Path(case_dir)/"cfd"/time/name) for name in ("p", "rho", "T", "U")}
    data = {"config": cfg.digest, "fit": [fit[k] for k in ("xc", "a", "b")],
            "time": time, "fields": fields}
    return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()


def _load_breakdown(case_dir):
    path = wf.results_dir(case_dir)/"breakdown.npz"
    if not path.is_file():
        raise wf.WorkflowError("run the breakdown step first")
    with np.load(path) as z:
        return {k: z[k].copy() for k in z.files}


def interface(case_dir, time=None, allow_unsuitable=False):
    case_dir = Path(case_dir)
    cfg = wf.load_case(case_dir)
    source = config.build_source(cfg)
    res = wf.results_dir(case_dir)
    data = _load_breakdown(case_dir)
    saved_time = str(data["time"])
    if time is not None and str(time) != saved_time:
        raise wf.WorkflowError("requested CFD time differs from breakdown; rerun breakdown at that time")
    mesh = read_mesh(case_dir/"cfd")
    lower, upper = _footprints(mesh)
    fit = fit_ellipse(lower, upper, data["retained"], source_radius=source.radius_m,
                      clearance=float(cfg.get("interface.fit.clearance_cells", 0.25))
                      * float(cfg.require("mesh.core_cell_over_D"))*source.diameter_m)
    doc = {"mode": "fitted-ellipse", "geometry": "planar", "time": saved_time,
           "label": cfg.get("label", cfg.profile), "fit": fit,
           "allow_unsuitable": bool(allow_unsuitable or cfg.get("interface.allow_unsuitable", False))}
    if fit["status"] != "fitted":
        doc.update(status=surface_tools.STATUS_NONE, reasons=[fit["reason"]])
        wf.write_json(res/"ellipse.json", fit)
        wf.write_json(res/"interface.json", doc)
        return doc
    identity = _identity(case_dir, cfg, fit, saved_time)
    old = json.loads((res/"interface.json").read_text()) if (res/"interface.json").is_file() else {}
    if (case_dir/"dsmc/constant/polyMesh/owner").is_file() and old.get("planar_identity") != identity:
        raise wf.WorkflowError("the existing DSMC mesh belongs to another fit/CFD state; use a new case directory")
    fit["time"] = saved_time
    doc.update(status="mesh_pending", reasons=["actual curved mesh awaits final suitability assessment"],
               planar_identity=identity)
    system = case_dir/"dsmc/system"
    system.mkdir(parents=True, exist_ok=True)
    (system/"blockMeshDict").write_text(render_dsmc_mesh(cfg, source, fit))
    wf.write_json(res/"ellipse.json", fit)
    wf.write_json(res/"interface.json", doc)
    np.save(res/"retained.npy", (((mesh.cell_centres[:, 0]-fit["xc"])/fit["a"])**2
                               + (mesh.cell_centres[:, 1]/fit["b"])**2 <= 1))
    return doc


def dsmc_setup(case_dir, np_ranks=None):
    case_dir = Path(case_dir)
    cfg = wf.load_case(case_dir)
    source = config.build_source(cfg)
    res = wf.results_dir(case_dir)
    if not (res/"ellipse.json").is_file():
        raise wf.WorkflowError("run interface first to fit the planar ellipse")
    fit = json.loads((res/"ellipse.json").read_text())
    if fit.get("status") != "fitted":
        raise wf.WorkflowError("no feasible ellipse exists; a DSMC interface cannot be built")
    doc = json.loads((res/"interface.json").read_text())
    if doc.get("planar_identity") != _identity(case_dir, cfg, fit, doc["time"]):
        raise wf.WorkflowError("CFD inputs changed since the ellipse fit; use a fresh case")
    if doc["status"] == surface_tools.STATUS_UNSUITABLE and not doc.get("allow_unsuitable"):
        raise wf.WorkflowError("interface is unsuitable; explicitly allow a machinery smoke test")
    ranks = int(np_ranks or cfg.get("resources.default_ranks", 1))
    if not 1 <= ranks <= int(cfg.get("resources.max_ranks", 4)):
        raise wf.WorkflowError("requested rank count is outside the planar profile's resource limits")
    from .mesh import plan
    estimate = resources.assess_configuration(cfg, source, plan(cfg, source))["full_dsmc"]
    run = dsmc.plan_run(cfg, source, estimate, float(cfg.require("domain.x_max_over_D"))*source.diameter_m)
    patches = ["interface"]
    if fit["source_half_height_m"] < source.radius_m-1e-9*source.diameter_m:
        patches.append("source")
    dsmc.write_case(case_dir/"dsmc", cfg, source, run, patches, ranks)
    (case_dir/"dsmc/system/blockMeshDict").write_text(render_dsmc_mesh(cfg, source, fit))
    (case_dir/"dsmc/planar.foam").touch()
    out = {"label": cfg.get("label", cfg.profile), "geometry": "planar",
           "interface_status": doc["status"], "inflow_patches": patches, "ranks": ranks,
           "run": run.describe(), "pending_actual_geometry": True,
           "planar_identity": doc["planar_identity"]}
    wf.write_json(res/"dsmc_setup.json", out)
    return out


def _check_geometry(mesh, cfd_mesh, retained, fit, source):
    if np.any(mesh.cell_volumes <= 0):
        raise wf.WorkflowError("curved DSMC mesh contains non-positive cells")
    h = float(np.min(cfd_mesh.cell_size()[:, :2]))
    z = mesh.points[:, 2]
    lo, hi = z.min(), z.max()
    front = mesh.points[np.abs(z-lo) < 1e-8*h, :2]
    back = mesh.points[np.abs(z-hi) < 1e-8*h, :2]
    distance, rows = cKDTree(front).query(back)
    if (len(front) != len(back) or len(front)+len(back) != len(z)
            or np.unique(rows).size != len(rows) or np.any(distance > 1e-7*h)
            or not np.isclose(hi-lo, source.depth_m, rtol=1e-7)):
        raise wf.WorkflowError("curved mesh must have matching front/back points and one uniform depth layer")
    for name in ("front", "back"):
        if name not in mesh.patches or mesh.patches[name].type != "empty":
            raise wf.WorkflowError("planar depth patches must be empty")
    sl = mesh.patch_slice("interface")
    faces = np.arange(sl.start, sl.stop)
    if not len(faces):
        raise wf.WorkflowError("curved mesh has no interface faces")
    vertices = np.unique(np.concatenate([mesh.labels[mesh.offsets[f]:mesh.offsets[f+1]] for f in faces]))
    xy = mesh.points[vertices, :2]
    q = ((xy[:, 0]-fit["xc"])/fit["a"])**2 + (xy[:, 1]/fit["b"])**2
    if np.max(np.abs(q-1)) > 2e-6:
        raise wf.WorkflowError("actual interface vertices do not match the fitted ellipse")
    lower, upper = _footprints(cfd_mesh)
    bad_lo, bad_hi = lower[~retained], upper[~retained]
    if np.any(footprint_clearance([fit[k] for k in ("xc", "a", "b")], bad_lo, bad_hi) < 1-1e-7):
        raise wf.WorkflowError("fitted ellipse intersects an excluded CFD footprint")
    sample_containing_cells(cfd_mesh, xy, retained)
    areas = mesh.face_areas[faces]
    magnitude = np.linalg.norm(areas, axis=1)
    normals = -areas/magnitude[:, None]
    gradient = (mesh.face_centres[faces, :2]-[fit["xc"], 0])/[fit["a"]**2, fit["b"]**2]
    if np.any(np.sum(normals[:, :2]*gradient, axis=1) <= 0) or np.max(np.abs(normals[:, 2])) > 1e-7:
        raise wf.WorkflowError("actual DSMC interface normals do not point out of the continuum core")
    covered = 2*fit["source_half_height_m"]*source.depth_m
    uncovered = np.linalg.norm(mesh.face_areas[mesh.patch_slice("source")], axis=1).sum()
    if not np.isclose(covered+uncovered, source.area_m2, rtol=1e-6):
        raise wf.WorkflowError("covered and uncovered source faces do not preserve the complete slit")
    total = np.sum(normals*magnitude[:, None], axis=0)-[covered, 0, 0]
    closure = float(np.linalg.norm(total)/(magnitude.sum()+covered))
    if closure > 1e-6:
        raise wf.WorkflowError("curved continuum region does not geometrically close")
    return faces, normals, magnitude, {"matching_depth_planes": True,
        "positive_volumes": True, "closure_relative": closure,
        "source_covered_area_m2": covered, "source_uncovered_area_m2": float(uncovered),
        "source_covered_area_fraction": covered/source.area_m2}


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
                      "uncovered_faces_injected_directly_by_dsmc": mesh.patches["source"].n_faces},
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


def _tracking_width(mesh):
    """Conservative in-plane cell height V / side-face area, excluding empty faces."""
    area = np.linalg.norm(mesh.face_areas, axis=1)
    side = np.abs(mesh.face_areas[:, 2]) < 1e-8*area
    width = np.full(mesh.n_cells, np.inf)
    owner = mesh.owner[side]
    np.minimum.at(width, owner, mesh.cell_volumes[owner]/area[side])
    internal = side[:mesh.n_internal]
    neighbour = mesh.neighbour[internal]
    np.minimum.at(width, neighbour, mesh.cell_volumes[neighbour]/area[:mesh.n_internal][internal])
    if not np.isfinite(width).all() or np.any(width <= 0):
        raise wf.WorkflowError("cannot determine positive in-plane DSMC tracking widths")
    return width


def dsmc_fields(case_dir):
    case_dir = Path(case_dir)
    cfg = wf.load_case(case_dir)
    source = config.build_source(cfg)
    res = wf.results_dir(case_dir)
    doc = json.loads((res/"interface.json").read_text())
    fit = json.loads((res/"ellipse.json").read_text())
    if doc.get("planar_identity") != _identity(case_dir, cfg, fit, doc["time"]):
        raise wf.WorkflowError("CFD state changed since interface fitting; use a fresh case")
    data = _load_breakdown(case_dir)
    if str(data["time"]) != doc["time"]:
        raise wf.WorkflowError("breakdown and fitted interface refer to different CFD times")
    cmesh, mesh = read_mesh(case_dir/"cfd"), read_mesh(case_dir/"dsmc")
    faces, normals, areas, geo = _check_geometry(mesh, cmesh, data["retained"], fit, source)
    cells = sample_containing_cells(cmesh, mesh.face_centres[faces], data["retained"])
    state = bd.read_state(case_dir/"cfd", doc["time"], source.gas)
    state_valid = ((state.rho > 0) & (state.T > 0) & np.isfinite(state.rho)
                   & np.isfinite(state.T) & np.isfinite(state.U).all(axis=1))
    all_cells = sample_containing_cells(cmesh, mesh.cell_centres, state_valid)
    ind = {k: data[k] for k in (*bd.INDICATORS, "valid")}
    assessment, arrays = _assess(mesh, faces, normals, areas, cells, all_cells, state, ind,
                                bd.gradients(cmesh, state), geo, cfg, source)
    doc.update(assessment)
    doc["smoke_only"] = doc["status"] != surface_tools.STATUS_SUITABLE
    fields = {"numberDensity": arrays["rho"]/source.gas.mass_kg,
              "massDensity": arrays["rho"], "velocity": arrays["U"],
              "translationalTemperature": arrays["T"], "faceArea": areas,
              "faceNormal": normals, "faceCentre": mesh.face_centres[faces],
              "valid": arrays["valid"], "cfdFace": np.full(len(faces), -1, dtype=int),
              "cfdCell": cells, "interpolationWeight": np.ones(len(faces))}
    mapping.save(res/"interface_fields.npz", fields)
    surface_tools.write_vtk(res/"interface.vtk", mesh, faces,
                            {"Mn": arrays["Mn"], "normal": normals,
                             "Kn_GLL_cfd_side": ind["Kn_GLL"][cells]})
    doc["state_consistency"] = mapping.consistency(fields, source.gas, state.p[cells])
    doc["artifacts"] = {"surface_vtk": str(res/"interface.vtk"),
                        "canonical_fields": str(res/"interface_fields.npz")}
    wf.write_json(res/"interface.json", doc)
    if doc["smoke_only"] and not doc.get("allow_unsuitable"):
        raise wf.WorkflowError("actual curved interface is unsuitable: "
                               + "; ".join(doc["reasons"])
                               + ". Re-run post-cfd --allow-unsuitable only for a machinery smoke test")
    setup_doc = json.loads((res/"dsmc_setup.json").read_text())
    cost_cfg = config.CaseConfig(config.deep_merge(cfg.doc, {"resources": {
        "max_ranks": setup_doc["ranks"]}}), cfg.path, cfg.profile)
    width = _tracking_width(mesh)
    estimate = resources.dsmc_estimate(
        cost_cfg, source, number_density=state.rho[all_cells]/source.gas.mass_kg,
        temperature=state.T[all_cells], speed=state.speed[all_cells],
        volumes=mesh.cell_volumes, sizes=mesh.cell_size()[:, :2], tracking_sizes=width,
        length_scale_m=float(mesh.points[:, 0].max()),
        basis="actual curved planar DSMC cells, containing-cell CFD states, real depth")
    run = dsmc.plan_run(cfg, source, estimate, float(mesh.points[:, 0].max()))
    patches = dsmc.inflow_patches(True, mesh.patches["source"].n_faces)
    dsmc.write_case(case_dir/"dsmc", cfg, source, run, patches, setup_doc["ranks"])
    inflow = {"interface": {"n": fields["numberDensity"], "U": fields["velocity"],
                             "T": fields["translationalTemperature"]}}
    if mesh.patches["source"].n_faces:
        inflow["source"] = mapping.source_inflow(mesh.patches["source"].n_faces, source)
    dsmc.write_fields(case_dir/"dsmc", mesh.patches, source.gas.name, inflow, source)
    setup_doc.update(interface_status=doc["status"], smoke_only=doc["smoke_only"],
                     run=run.describe(), dsmc_cells=mesh.n_cells, inflow_patches=patches,
                     pending_actual_geometry=False, flux_depth_m=source.depth_m)
    wf.write_json(res/"dsmc_setup.json", setup_doc)
    result = {"geometry": "planar", "interface_status": doc["status"], "mapping": {
        "method": "containing accepted CFD cell; bounded zero-slope reconstruction",
        "faces": len(faces), "shared_cfd_faces": False, "positive_states": True,
        "source_coverage_verified": True}, "geometry_checks": geo,
        "min_tracking_width_m": float(width.min()), "state_consistency": doc["state_consistency"],
        "patches": {k: {"type": v.type, "n_faces": v.n_faces} for k, v in mesh.patches.items()}}
    wf.write_json(res/"dsmc_fields.json", result)
    (case_dir/"dsmc/planar.foam").touch()
    return result


def reference(case_dir):
    doc = {"status": "not_applicable",
           "reason": "the circular-disk collisionless reference does not describe a planar slit plume"}
    wf.write_json(wf.results_dir(case_dir)/"reference.json", doc)
    return doc


def dsmc_post(case_dir, time=None):
    from . import post, particles
    case_dir = Path(case_dir)
    cfg = wf.load_case(case_dir)
    source = config.build_source(cfg)
    res = wf.results_dir(case_dir)
    dsmc_dir = case_dir/"dsmc"
    mesh = read_mesh(dsmc_dir)
    time = time or wf._latest(dsmc_dir)
    setup_doc = json.loads((res/"dsmc_setup.json").read_text())
    try:
        fields = post.dsmc_fieldaverage(dsmc_dir/time, source.gas)
    except (OSError, ValueError) as exc:
        raise wf.WorkflowError(f"planar DSMC mean fields are unavailable at {time}: {exc}") from exc
    if len(fields["n"]) != mesh.n_cells:
        raise wf.WorkflowError("DSMC mean field length does not match its curved mesh")
    valid = ((fields["n"] > 0) & (fields["T"] > 0) & np.isfinite(fields["n"])
             & np.isfinite(fields["T"]) & np.isfinite(fields["U"]).all(axis=1))
    np.savez_compressed(res/"planar_fields.npz", centres=mesh.cell_centres,
                        volumes=mesh.cell_volumes, valid=valid, **fields)
    h = float(cfg.require("mesh.core_cell_over_D"))*source.diameter_m
    axis = np.abs(mesh.cell_centres[:, 1]) <= .51*h
    profile = post.line_profile(mesh.cell_centres, mesh.cell_volumes, axis, {
        "n_over_ne": fields["n"]/source.number_density,
        "T_over_Te": fields["T"]/source.temperature_K,
        "Ux_over_ae": fields["U"][:, 0]/source.thermal_speed})
    profile["x_over_D"] = profile.pop("x")/source.diameter_m
    post.write_csv(res/"centreline.csv", {"x_over_D": profile.pop("x_over_D"), **profile})
    bd.write_fields(dsmc_dir, time, mesh, {
        "planar_n_over_ne": np.where(valid, fields["n"]/source.number_density, np.nan),
        "planar_T_over_Te": np.where(valid, fields["T"]/source.temperature_K, np.nan),
        "planar_Ux_over_ae": np.where(valid, fields["U"][:, 0]/source.thermal_speed, np.nan)})
    doc = {"geometry": "planar", "time": time, "label": cfg.get("label", cfg.profile),
           "smoke_only": setup_doc.get("smoke_only", False), "cells": mesh.n_cells,
           "valid_cells": int(valid.sum()), "invalid_or_unoccupied_cells": int((~valid).sum()),
           "flux_depth_m": source.depth_m,
           "sampling": {"parcels_per_cell_mean": float(np.mean(fields["parcels_per_step"])),
                        "note": "fieldAverage moments on the actual DSMC mesh; planar slit, three molecular velocity components"},
           "analytical_reference": reference(case_dir),
           "artifacts": {"planar_fields": str(res/"planar_fields.npz"),
                         "dsmc_centreline": str(res/"centreline.csv"),
                         "paraview": str(dsmc_dir/"planar.foam")}}
    log = dsmc_dir/"log.dsmcFoam"
    count = wf._log_series(log, r"Number of molecules\s*=\s*([0-9.eE+]+)")
    if len(count) > 10 and count[-1] > 0:
        quarter = max(2, len(count)//4)
        doc["stationarity"] = {"parcels_last_quarter_rel_change":
                              float((count[-1]-count[-quarter])/count[-1]),
                              "parcels_final": float(count[-1])}
    inserted = particles.inserted_from_log(log) if log.is_file() else np.empty(0)
    run = setup_doc.get("run", {})
    dt, weight = run.get("dt_s"), run.get("weight")
    average_start = float(run.get("average_start_s", 0))
    if inserted.size and dt and weight and float(time) > average_start:
        window = int(round((float(time)-average_start)/dt))
        if 0 < window <= len(inserted):
            doc["injection"] = {"measured_gross_number_rate_per_s":
                                float(inserted[-window:].mean()*weight/dt),
                                "note": "gross injection, not net mass transport"}
    (dsmc_dir/"planar.foam").touch()
    wf.write_json(res/"dsmc_post.json", doc)
    return doc


def plot_interface(case_dir):
    """Raw CFD boundary / fitted ellipse, then fields on the actual DSMC mesh."""
    case_dir = Path(case_dir)
    cfg = wf.load_case(case_dir)
    source = config.build_source(cfg)
    res = wf.results_dir(case_dir)
    if not (res/"breakdown.npz").is_file() and not (res/"planar_fields.npz").is_file():
        return []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from matplotlib.collections import LineCollection
    except ImportError:
        wf.write_json(res/"plot_status.json", {"status": "unavailable",
                      "reason": "optional matplotlib dependency is not installed; inspect the .foam cases"})
        return []
    from .ellipse import ellipse_points
    d = source.diameter_m
    out = []
    fit = json.loads((res/"ellipse.json").read_text()) if (res/"ellipse.json").is_file() else {}
    if (res/"breakdown.npz").is_file():
        data = _load_breakdown(case_dir)
        mesh = read_mesh(case_dir/"cfd")
        fig, ax = plt.subplots(figsize=(8, 4.6), constrained_layout=True)
        values = data.get("Kn_GLL", np.full(mesh.n_cells, np.nan))
        image = ax.scatter(mesh.cell_centres[:, 0]/d, mesh.cell_centres[:, 1]/d,
                           c=np.maximum(values, 1e-8), s=14, marker="s",
                           cmap="viridis", norm=LogNorm(vmin=1e-4, vmax=1))
        fig.colorbar(image, ax=ax, label="Kn_GLL")
        retained = data["retained"]
        cuts = np.flatnonzero(retained[mesh.owner[:mesh.n_internal]] != retained[mesh.neighbour])
        segments = []
        for f in cuts:
            points = np.unique(np.round(mesh.face_vertices(f)[:, :2], 12), axis=0)
            if len(points) == 2:
                segments.append(points/d)
        if segments:
            ax.add_collection(LineCollection(segments, colors="#c44e52", linewidths=1,
                                             label="raw source-connected breakdown boundary"))
        if fit.get("status") == "fitted":
            arc = ellipse_points(fit, 401)/d
            ax.plot(arc[:, 0], arc[:, 1], color="#202020", lw=2, label="conservative fitted ellipse")
        ax.plot([0, 0], [-.5, .5], color="#2369b4", lw=3, label="source slit")
        ax.set(xlabel="x / D", ylabel="y / D", title="Planar prototype: raw breakdown and smooth interface")
        ax.set_aspect("equal")
        ax.legend(fontsize=8, loc="upper right")
        path = res/"breakdown_ellipse.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        out.append(path)
    if (res/"planar_fields.npz").is_file():
        with np.load(res/"planar_fields.npz") as z:
            centres, valid = z["centres"], z["valid"]
            values = (z["n"]/source.number_density, z["T"]/source.temperature_K,
                      z["U"][:, 0]/source.thermal_speed)
        fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
        for ax, value, label in zip(axes, values, ("n / n_e", "T / T_e", "Ux / a_e")):
            image = ax.scatter(centres[valid, 0]/d, centres[valid, 1]/d,
                               c=value[valid], s=5, cmap="viridis")
            fig.colorbar(image, ax=ax, label=label)
            ax.set(xlabel="x / D", ylabel="y / D")
            ax.set_aspect("equal")
        fig.suptitle("Planar DSMC moments on the curved mesh")
        path = res/"planar_dsmc.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        out.append(path)
    wf.write_json(res/"plot_status.json", {"status": "written", "artifacts": out})
    return out
