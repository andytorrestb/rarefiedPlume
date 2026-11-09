"""Fully 3D fitted-ellipsoid hybrid path; native meshing and actual-face transfer."""
from __future__ import annotations
import json
import math
import re
import shutil
from pathlib import Path
import numpy as np

from . import config, breakdown as bd, dsmc, mapping, resources, fitted
from . import workflow as wf
from .ellipsoid import fit_ellipsoid, verify_ellipsoid
from .ellipsoid_mesh import write_mesh_inputs
from .meshgeom import read_mesh
from . import interface as itf


def _footprints(mesh):
    half = mesh.cell_size()/2
    return mesh.cell_centres-half, mesh.cell_centres+half


def _mesh_identity(case_dir, expected=None):
    directory = Path(case_dir)/"dsmc/constant/polyMesh"
    actual = {name: wf.provenance.sha256(directory/name)
              for name in ("points", "faces", "owner", "neighbour", "boundary")}
    if expected is not None and actual != expected:
        raise wf.WorkflowError("DSMC mesh identity changed; use a fresh case instead of reusing the cloud")
    return actual


def certify_facets(vertices, fit, lower, upper, retained, tolerance):
    axes = np.array([fit["a"], fit["b"], fit["b"]])
    q = np.sum(((vertices-np.array([fit["xc"], 0, 0]))/axes)**2, axis=1)
    deviation = float(np.max(np.abs(np.sqrt(q)-1)))
    if deviation > tolerance:
        raise wf.WorkflowError("snapped facets exceed the fitted geometry tolerance")
    enclosure = dict(fit, clearance_m=0)
    growth = max(1.0, float(np.sqrt(q.max())))*(1+1e-10)
    enclosure.update(a=fit["a"]*growth, b=fit["b"]*growth)
    check = verify_ellipsoid(enclosure, lower, upper, retained)
    if not check["feasible"]:
        raise wf.WorkflowError("snapped-facet enclosure enters an excluded continuum voxel")
    if np.any(vertices < lower.min(axis=0)-1e-9) or np.any(vertices > upper.max(axis=0)+1e-9):
        raise wf.WorkflowError("snapped interface extends outside the CFD domain")
    return {"max_relative_radius_deviation": deviation, "enclosure_growth": growth,
            "voxel_verification": check}


def interface(case_dir, time=None, allow_unsuitable=False):
    case_dir = Path(case_dir)
    cfg = wf.load_case(case_dir)
    source = config.build_source(cfg)
    res = wf.results_dir(case_dir)
    data = fitted._load_breakdown(case_dir)
    saved_time = str(data["time"])
    if time is not None and time != saved_time:
        raise wf.WorkflowError("requested time differs from breakdown; rerun breakdown")
    mesh = read_mesh(case_dir/"cfd")
    lower, upper = _footprints(mesh)
    fit = fit_ellipsoid(lower, upper, data["retained"], source_radius=source.radius_m,
                        clearance=float(cfg.get("interface.fit.clearance_cells", .25))
                        *float(cfg.require("mesh.core_cell_over_D"))*source.diameter_m)
    fitted._assert_breakdown_current(case_dir, data, cfg)
    doc = {"mode": "fitted-ellipsoid", "geometry": "cartesian", "time": saved_time,
           "label": cfg.get("label", cfg.profile), "fit": fit,
           "allow_unsuitable": bool(allow_unsuitable or cfg.get("interface.allow_unsuitable", False))}
    if fit["status"] != "fitted":
        doc.update(status=itf.STATUS_NONE, reasons=[fit["reason"]])
    else:
        identity = fitted._identity(case_dir, cfg, fit, saved_time)
        old = json.loads((res/"interface.json").read_text()) if (res/"interface.json").is_file() else {}
        if (case_dir/"dsmc/constant/polyMesh/owner").is_file() and old.get("fitted_identity") != identity:
            raise wf.WorkflowError("existing DSMC mesh has a different fit/state; use a fresh case")
        fit["time"] = saved_time
        doc.update(status="mesh_pending", fitted_identity=identity,
                   reasons=["actual 3D mesh awaits suitability assessment"])
        write_mesh_inputs(case_dir, cfg, source, fit)
    wf.write_json(res/"ellipsoid.json", fit)
    wf.write_json(res/"interface.json", doc)
    return doc


def dsmc_setup(case_dir, np_ranks=None):
    from .mesh import plan
    case_dir = Path(case_dir)
    cfg = wf.load_case(case_dir)
    source = config.build_source(cfg)
    res = wf.results_dir(case_dir)
    if not (res/"ellipsoid.json").is_file():
        raise wf.WorkflowError("run interface first to fit the ellipsoid")
    fit = json.loads((res/"ellipsoid.json").read_text())
    if fit.get("status") != "fitted":
        raise wf.WorkflowError("no feasible ellipsoid; DSMC cannot be initialized")
    doc = json.loads((res/"interface.json").read_text())
    if doc.get("fitted_identity") != fitted._identity(case_dir, cfg, fit, doc["time"]):
        raise wf.WorkflowError("CFD inputs changed since fitting; use a fresh case")
    ranks = int(np_ranks or cfg.get("resources.default_ranks", 192))
    if not 1 <= ranks <= int(cfg.get("resources.max_ranks", 192)):
        raise wf.WorkflowError("requested ranks exceed the profile limit")
    old = json.loads((res/"dsmc_setup.json").read_text()) if (res/"dsmc_setup.json").is_file() else {}
    if old.get("fitted_identity") not in (None, doc["fitted_identity"]):
        raise wf.WorkflowError("DSMC checkpoints belong to another interface")
    mesh_identity = None
    if (case_dir/"dsmc/constant/polyMesh/owner").is_file():
        if old.get("mesh_identity") is None and (
                any((case_dir/"dsmc").glob("processor[0-9]*"))
                or any((case_dir/"dsmc").glob("*/lagrangian"))):
            raise wf.WorkflowError("saved cloud has no verified DSMC mesh identity; use a fresh case")
        mesh_identity = _mesh_identity(case_dir, old.get("mesh_identity"))
    else:
        build = (case_dir/"mesh-build").resolve()
        for generated in [build/"constant/polyMesh", *build.glob("processor[0-9]*")]:
            if generated.exists():
                if generated.is_symlink() or not generated.resolve().is_relative_to(build):
                    raise wf.WorkflowError("mesher scratch path escapes its dedicated workspace")
                shutil.rmtree(generated)
    estimate = resources.assess_configuration(cfg, source, plan(cfg, source))["full_dsmc"]
    run = dsmc.plan_run(cfg, source, estimate, float(cfg.require("domain.x_max_over_D"))*source.diameter_m)
    dsmc.write_case(case_dir/"dsmc", cfg, source, run, ["interface", "source"], ranks)
    write_mesh_inputs(case_dir, cfg, source, fit)
    (case_dir/"dsmc/spatial.foam").touch()
    result = {"geometry": "cartesian", "ranks": ranks, "run": run.describe(),
              "fitted_identity": doc["fitted_identity"], "pending_actual_geometry": True,
              "mesh_identity": mesh_identity, "inflow_patches": ["interface", "source"],
              "interface_status": doc["status"]}
    wf.write_json(res/"dsmc_setup.json", result)
    return result


def publish_mesh(case_dir):
    """Publish a complete reconstructed mesh, never mesher processor directories."""
    case_dir = Path(case_dir)
    build = case_dir/"mesh-build/constant/polyMesh"
    for name in ("points", "faces", "owner", "neighbour", "boundary"):
        if not (build/name).is_file():
            raise wf.WorkflowError("meshing did not produce a complete reconstructed mesh")
    boundary = build/"boundary"
    text = boundary.read_text()
    text = re.sub(r"(?s)(\bexitPlane\s*\{.*?\btype\s+)\w+;", r"\g<1>patch;", text)
    boundary.write_text(text)
    target = case_dir/"dsmc/constant/polyMesh"
    if target.exists():
        raise wf.WorkflowError("refusing to overwrite an existing DSMC mesh")
    pending = target.parent/".ellipsoidMesh.pending"
    if pending.exists():
        if pending.is_symlink() or pending.resolve().parent != target.parent.resolve():
            raise wf.WorkflowError("invalid pending mesh path")
        shutil.rmtree(pending)
    shutil.copytree(build, pending)
    pending.rename(target)


def _geometry(mesh, cmesh, data, fit, cfg, source):
    if np.any(mesh.cell_volumes <= 0) or any(p.type in ("empty", "wedge") for p in mesh.patches.values()):
        raise wf.WorkflowError("ellipsoid DSMC mesh must have positive full-3D cells")
    sl = mesh.patch_slice("interface")
    faces = np.arange(sl.start, sl.stop)
    if not len(faces):
        raise wf.WorkflowError("snapped mesh has no interface")
    vertices = np.unique(np.concatenate([mesh.labels[mesh.offsets[f]:mesh.offsets[f+1]] for f in faces]))
    lo, hi = _footprints(cmesh)
    checks = certify_facets(mesh.points[vertices], fit, lo, hi, data["retained"],
                             float(cfg.get("interface.fit.geometry_relative_tolerance", .02)))
    area_vectors = -mesh.face_areas[faces]
    areas = np.linalg.norm(area_vectors, axis=1)
    normals = area_vectors/areas[:, None]
    gradient = (mesh.face_centres[faces]-[fit["xc"], 0, 0])/[fit["a"]**2, fit["b"]**2, fit["b"]**2]
    if np.any(np.sum(normals*gradient, axis=1) <= 0):
        raise wf.WorkflowError("snapped interface has reversed normals")
    covered = float(area_vectors[:, 0].sum())
    source_info = mesh.patches.get("source")
    direct = float(np.linalg.norm(mesh.face_areas[mesh.patch_slice("source")], axis=1).sum()) if source_info else 0.0
    error = abs((covered+direct)/source.area_m2-1)
    tolerance = float(cfg.get("interface.fit.max_source_area_error", .01))
    if covered <= 0 or error > tolerance:
        raise wf.WorkflowError(f"source coverage area error {error:.4g} exceeds tolerance")
    cfd_area = float(np.linalg.norm(cmesh.face_areas[cmesh.patch_slice("source")], axis=1).sum())
    cfd_error = abs(cfd_area/source.area_m2-1)
    if not math.isfinite(cfd_error) or cfd_error > tolerance:
        raise wf.WorkflowError(f"CFD source area error {cfd_error:.4g} exceeds tolerance {tolerance:.4g}")
    closure = float(np.linalg.norm(area_vectors.sum(axis=0)-[covered, 0, 0])/(areas.sum()+covered))
    if closure > 1e-6:
        raise wf.WorkflowError("retained region does not close at the source plane")
    checks.update(positive_volumes=True, closure_relative=closure,
                  source_covered_area_m2=covered, source_uncovered_area_m2=direct,
                  source_covered_area_fraction=covered/source.area_m2, source_area_error=error,
                  cfd_source_area_m2=cfd_area, cfd_source_area_error=cfd_error)
    return faces, normals, areas, checks


def _widths(mesh):
    area = np.linalg.norm(mesh.face_areas, axis=1)
    width = np.full(mesh.n_cells, np.inf)
    np.minimum.at(width, mesh.owner, mesh.cell_volumes[mesh.owner]/area)
    np.minimum.at(width, mesh.neighbour, mesh.cell_volumes[mesh.neighbour]/area[:mesh.n_internal])
    return width


def dsmc_fields(case_dir):
    case_dir = Path(case_dir)
    cfg = wf.load_case(case_dir)
    source = config.build_source(cfg)
    res = wf.results_dir(case_dir)
    doc = json.loads((res/"interface.json").read_text())
    fit = json.loads((res/"ellipsoid.json").read_text())
    if doc.get("fitted_identity") != fitted._identity(case_dir, cfg, fit, doc["time"]):
        raise wf.WorkflowError("CFD inputs changed since fitting")
    data = fitted._load_breakdown(case_dir)
    setup = json.loads((res/"dsmc_setup.json").read_text())
    mesh_identity = _mesh_identity(case_dir, setup.get("mesh_identity"))
    cmesh, mesh = read_mesh(case_dir/"cfd"), read_mesh(case_dir/"dsmc")
    faces, normals, areas, geo = _geometry(mesh, cmesh, data, fit, cfg, source)
    cells = fitted.sample_containing_cells(cmesh, mesh.face_centres[faces], data["retained"], axes=(0, 1, 2))
    state = bd.read_state(case_dir/"cfd", doc["time"], source.gas)
    valid = (state.rho > 0) & (state.T > 0) & np.isfinite(state.rho) & np.isfinite(state.T) & np.isfinite(state.U).all(axis=1)
    all_cells = fitted.sample_containing_cells(cmesh, mesh.cell_centres, valid, axes=(0, 1, 2))
    ind = {k: data[k] for k in (*bd.INDICATORS, "valid")}
    assessment, arrays = fitted._assess(mesh, faces, normals, areas, cells, all_cells,
                                       state, ind, bd.gradients(cmesh, state), geo, cfg, source)
    fitted._assert_breakdown_current(case_dir, data, cfg)
    doc.update(assessment)
    doc.pop("flux_depth_m", None)
    doc["smoke_only"] = doc["status"] != itf.STATUS_SUITABLE
    fields = {"numberDensity": arrays["rho"]/source.gas.mass_kg, "massDensity": arrays["rho"],
              "velocity": arrays["U"], "translationalTemperature": arrays["T"],
              "faceArea": areas, "faceNormal": normals, "faceCentre": mesh.face_centres[faces],
              "valid": arrays["valid"], "cfdCell": cells, "cfdFace": np.full(len(cells), -1, int),
              "interpolationWeight": np.ones(len(cells))}
    mapping.save(res/"interface_fields.npz", fields)
    itf.write_vtk(res/"interface.vtk", mesh, faces, {"Mn": arrays["Mn"], "normal": normals})
    doc["state_consistency"] = mapping.consistency(fields, source.gas, state.p[cells])
    wf.write_json(res/"interface.json", doc)
    if doc["smoke_only"] and not doc.get("allow_unsuitable"):
        raise wf.WorkflowError("3D interface unsuitable: "+"; ".join(doc["reasons"])+"; explicit smoke-only override required")
    cost_cfg = config.CaseConfig(config.deep_merge(cfg.doc, {"resources": {"max_ranks": setup["ranks"]}}), cfg.path, cfg.profile)
    estimate = resources.dsmc_estimate(cost_cfg, source, number_density=state.rho[all_cells]/source.gas.mass_kg,
        temperature=state.T[all_cells], speed=state.speed[all_cells], volumes=mesh.cell_volumes,
        sizes=mesh.cell_size(), tracking_sizes=_widths(mesh), length_scale_m=float(mesh.points[:, 0].max()),
        basis="actual full-3D remeshed cells and containing CFD states")
    peak = resources.peak_memory_estimate(cfg, cmesh.n_cells, estimate, ranks=setup["ranks"])
    wf.write_json(res/"memory_estimate.json", peak)
    if not peak["fits_memory_budget"]:
        raise wf.WorkflowError("estimated peak memory exceeds configured node budget; DSMC not initialized")
    _mesh_identity(case_dir, mesh_identity)
    run = dsmc.plan_run(cfg, source, estimate, float(mesh.points[:, 0].max()))
    n_source = mesh.patches["source"].n_faces if "source" in mesh.patches else 0
    patches = dsmc.inflow_patches(True, n_source)
    dsmc.write_case(case_dir/"dsmc", cfg, source, run, patches, setup["ranks"])
    inflow = {"interface": {"n": fields["numberDensity"], "U": fields["velocity"], "T": fields["translationalTemperature"]}}
    if n_source:
        inflow["source"] = mapping.source_inflow(n_source, source)
    dsmc.write_fields(case_dir/"dsmc", mesh.patches, source.gas.name, inflow, source)
    setup.update(run=run.describe(), dsmc_cells=mesh.n_cells, inflow_patches=patches,
                 pending_actual_geometry=False, interface_status=doc["status"],
                 mesh_identity=mesh_identity, smoke_only=doc["smoke_only"], memory=peak)
    wf.write_json(res/"dsmc_setup.json", setup)
    result = {"geometry": "cartesian", "interface_faces": len(faces), "geometry_checks": geo,
              "mapping": "containing accepted 3D CFD cells", "memory": peak}
    wf.write_json(res/"dsmc_fields.json", result)
    return result


def dsmc_post(case_dir, time=None):
    """Reuse 3D moment/flux readers without scattering to CFD cell numbering."""
    return wf._dsmc_post_3d(Path(case_dir), time=time)


def plots(case_dir):
    """Thin meridional sections on each mesh; never project all 3D cells."""
    case_dir = Path(case_dir)
    cfg = wf.load_case(case_dir)
    source = config.build_source(cfg)
    res = wf.results_dir(case_dir)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from matplotlib.collections import LineCollection
    except ImportError:
        wf.write_json(res/"plot_status.json", {"status": "unavailable",
                      "reason": "optional matplotlib missing; inspect spatial.foam and interface.vtk"})
        return []
    from .ellipse import ellipse_points
    from . import report
    d = source.diameter_m
    label = cfg.get("label", cfg.profile)
    doc = json.loads((res/"interface.json").read_text()) if (res/"interface.json").is_file() else {}
    if doc.get("smoke_only"):
        label += " — MACHINERY SMOKE"
    fit = json.loads((res/"ellipsoid.json").read_text()) if (res/"ellipsoid.json").is_file() else {}
    arc = ellipse_points(fit, 401)/d if fit.get("status") == "fitted" else None
    out = []
    if (res/"breakdown.npz").is_file():
        data = fitted._load_breakdown(case_dir)
        mesh = read_mesh(case_dir/"cfd")
        c, widths = mesh.cell_centres, mesh.cell_size()
        select = np.abs(c[:, 2]) <= .51*widths[:, 2]
        fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
        values = data["Kn_GLL"][select]
        values = np.where(np.isfinite(values) & (values > 0), values, np.nan)
        image = ax.scatter(c[select, 0]/d, c[select, 1]/d, c=values, s=6, marker="s",
                           cmap="viridis", norm=LogNorm(vmin=1e-4, vmax=1))
        fig.colorbar(image, ax=ax, label="Kn_GLL")
        retained = data["retained"]
        cuts = np.flatnonzero((retained[mesh.owner[:mesh.n_internal]] != retained[mesh.neighbour])
                             & select[mesh.owner[:mesh.n_internal]])
        segments = []
        for face in cuts:
            v = mesh.face_vertices(face)
            if v[:, 2].min() > 1e-10 or v[:, 2].max() < -1e-10 or np.ptp(v[:, 2]) < 1e-10:
                continue
            cross = []
            for a, b in zip(v, np.roll(v, -1, axis=0)):
                if abs(a[2]) < 1e-10:
                    cross.append(a[:2])
                elif a[2]*b[2] < 0:
                    cross.append((a+(b-a)*(-a[2]/(b[2]-a[2])))[:2])
            points = np.unique(np.round(cross, 12), axis=0)
            if len(points) == 2:
                segments.append(points/d)
        if segments:
            ax.add_collection(LineCollection(segments, colors="#c44e52", linewidths=1,
                                             label="raw continuum boundary at z = 0"))
        if arc is not None:
            ax.plot(arc[:, 0], arc[:, 1], color="#202020", lw=2, label="fitted spheroid section")
        ax.plot([0, 0], [-.5, .5], color="#2369b4", lw=3, label="circular source section")
        ax.set(xlabel="x / D", ylabel="y / D", title=label+"\nCFD cells adjacent to z = 0")
        ax.set_aspect("equal")
        ax.legend(fontsize=8)
        path = res/"breakdown_ellipsoid.png"
        fig.savefig(path, dpi=130)
        plt.close(fig)
        out.append(path)
    if (res/"spatial_fields.npz").is_file():
        with np.load(res/"spatial_fields.npz") as fields:
            c = fields["centres"]
            select = fields["valid"] & fields["meridional"]
            values = (fields["n"]/source.number_density, fields["T"]/source.temperature_K,
                      fields["U"][:, 0]/source.thermal_speed)
        fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
        for ax, value, name in zip(axes, values, ("n / n_e", "T / T_e", "Ux / a_e")):
            image = ax.scatter(c[select, 0]/d, c[select, 1]/d, c=value[select], s=4, cmap="viridis")
            fig.colorbar(image, ax=ax, label=name)
            if arc is not None:
                ax.plot(arc[:, 0], arc[:, 1], color="#202020", lw=1)
            ax.set(xlabel="x / D", ylabel="y / D")
            ax.set_aspect("equal")
        fig.suptitle(label+"\nDSMC cells adjacent to z = 0; unoccupied/invalid cells omitted")
        path = res/"spatial_dsmc.png"
        fig.savefig(path, dpi=130)
        plt.close(fig)
        out.append(path)
    if (res/"centreline.csv").is_file():
        data = np.atleast_1d(np.genfromtxt(res/"centreline.csv", delimiter=",", names=True))
        fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
        for ax, quantity, suffix, ylabel in zip(axes, ("n", "T", "U"),
                ("over_ne", "over_Te", "over_ae_x"), ("n / n_e", "T / T_e", "Ux / a_e")):
            for key, name, style in (("fa", "DSMC fieldAverage", "-"),
                                     ("ref", "collisionless limit", ":")):
                column = f"{quantity}_{key}_{suffix}"
                if column in data.dtype.names:
                    ax.plot(data["x_over_D"], data[column], style, label=name, color=report.SERIES[f"dsmc_{key}" if key == "fa" else key])
            report._style(ax, "x / D (finite axis tube; DSMC portion)", ylabel)
        axes[0].legend(frameon=False, fontsize=8)
        fig.suptitle(label)
        path = res/"centreline.png"
        fig.savefig(path, dpi=130)
        plt.close(fig)
        out.append(path)
    wf.write_json(res/"plot_status.json", {"status": "written", "artifacts": out,
                  "section": "cells with |z centre| <= 0.51 of their z bounding width"})
    return out
