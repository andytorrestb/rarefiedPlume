"""The explicit steps of one hybrid case, each a function and a CLI command.

Layout of a case directory ``<case>``::

    <case>/hybrid.resolved.yaml   the merged configuration every step reads
    <case>/cfd/                   the rhoCentralFoam precursor (OpenFOAM case)
    <case>/cfd_bg<ratio>/         the same precursor at another background floor
    <case>/dsmc/                  the DSMC case: the CFD mesh minus the retained
                                  continuum region (or all of it: full DSMC)
    <case>/results/               JSON reports, CSV profiles, PNG plots, VTK

Steps (``python -m plumetools.hybridplume <step> <case>``):

==================  ==========================================================
``setup``           resolve the profile, assess resources, write the CFD case
``check-mesh``      built source patch == predicted staircase (count, area)
``cfd-report``      budgets of the converged precursor; background sensitivity
``breakdown``       3-D Kn_rho, Kn_T, Kn_U, Kn_GLL, B_exp on the CFD cells
``interface``       select (automatic/manual/none) and assess the surface
``dsmc-setup``      DSMC dictionaries + the cellSet that subsetMesh removes
``dsmc-fields``     map CFD states onto the DSMC inflow faces; write 0/
==================  ==========================================================

OpenFOAM commands are run by ``cases/cai2012-hybrid/Allrun``, not from here:
each step either writes inputs for them or reads what they wrote.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import yaml

from . import cfd as cfd_mod
from . import config as config_mod
from . import mesh as mesh_mod
from . import provenance

RESOLVED = "hybrid.resolved.yaml"


class WorkflowError(RuntimeError):
    """A step that cannot proceed; the message says what is missing."""


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, Path):
        return str(o)
    raise TypeError(f"not JSON-serialisable: {type(o)}")


def write_json(path: Path, doc) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, indent=2, default=_json_default, allow_nan=True) + "\n")
    return path


def load_case(case_dir: Path) -> config_mod.CaseConfig:
    path = Path(case_dir) / RESOLVED
    if not path.is_file():
        raise WorkflowError(f"{path} does not exist; run the 'setup' step first")
    doc = yaml.safe_load(path.read_text())
    return config_mod.CaseConfig(doc=doc["config"], path=path, profile=doc["profile"])


def results_dir(case_dir: Path) -> Path:
    out = Path(case_dir) / "results"
    out.mkdir(parents=True, exist_ok=True)
    return out


# --------------------------------------------------------------------------- #
# setup
# --------------------------------------------------------------------------- #

def setup(case_dir: Path, config_path: Path, profile: str | None, *,
          n_ranks: int | None = None, overrides: dict | None = None,
          background_ratio: float | None = None, subdir: str = "cfd") -> dict:
    """Resolve the configuration, assess resources, write the CFD case.

    ``background_ratio`` writes a sensitivity copy of the precursor at another
    vacuum floor into ``<case>/<subdir>``; everything else is identical.
    """
    from . import resources

    case_dir = Path(case_dir)
    case_dir.mkdir(parents=True, exist_ok=True)
    over = dict(overrides or {})
    if background_ratio is not None:
        over = config_mod.deep_merge(over, {"cfd": {"background_pressure_ratio":
                                                    float(background_ratio)}})
    cfg = config_mod.load(config_path, profile, over)
    resolved = case_dir / RESOLVED
    if subdir == "cfd":
        if resolved.is_file():
            old = yaml.safe_load(resolved.read_text())
            if old.get("digest") != cfg.digest and any(
                    (case_dir / d).exists() for d in ("dsmc",)):
                raise WorkflowError(
                    f"{case_dir} already has a DSMC stage built from a different "
                    "configuration; start a new case directory or clean it")
        resolved.write_text(yaml.safe_dump({
            "profile": cfg.profile, "digest": cfg.digest,
            "source_config": str(config_path), "config": cfg.doc}, sort_keys=False))

    source = config_mod.build_source(cfg)
    plan = mesh_mod.plan(cfg, source)
    pred = mesh_mod.predicted_source(plan, source)
    ranks = int(n_ranks or cfg.get("resources.max_ranks", 1))
    ranks = max(1, min(ranks, int(cfg.get("resources.max_ranks", ranks))))

    target = case_dir / subdir
    mesh_mod.write_mesh_setup(target, plan, source)
    written = cfd_mod.write_case(target, cfg, source, plan, n_ranks=ranks)
    assessment = resources.assess_configuration(cfg, source, plan)
    doc = {
        "label": cfg.get("label", cfg.profile),
        "profile": cfg.profile, "config_digest": cfg.digest,
        "source": source.describe(),
        "departures_from_paper": config_mod.departures_from_paper(cfg),
        "mesh": plan.describe(), "source_discretisation": pred.describe(),
        "cfd": {k: v for k, v in written.items() if k != "written"},
        "resources": assessment, "n_ranks": ranks,
        "provenance": provenance.collect(Path(__file__).resolve().parents[2]),
    }
    name = "setup.json" if subdir == "cfd" else f"setup_{subdir}.json"
    write_json(results_dir(case_dir) / name, doc)
    return doc


# --------------------------------------------------------------------------- #
# check-mesh
# --------------------------------------------------------------------------- #

def check_mesh(case_dir: Path, subdir: str = "cfd") -> dict:
    """The built ``source`` patch must be the predicted staircase."""
    from plumetools.mesh.patches import read_patches

    cfg = load_case(case_dir)
    source = config_mod.build_source(cfg)
    plan = mesh_mod.plan(cfg, source)
    pred = mesh_mod.predicted_source(plan, source)
    built = read_patches(Path(case_dir) / subdir, [mesh_mod.PATCHES["source"]])
    geo = built[mesh_mod.PATCHES["source"]]
    normals = geo.normals
    out = {
        "predicted": pred.describe(),
        "built": {"n_faces": geo.n_faces, "area_m2": geo.area_m2,
                  "relative_area_error": geo.area_m2 / source.area_m2 - 1.0,
                  "outward_normal_mean": normals.mean(axis=0).tolist(),
                  "centroid_m": (geo.centres * geo.areas[:, None]).sum(0) / geo.area_m2},
    }
    out["matches_prediction"] = bool(
        geo.n_faces == pred.n_faces
        and abs(geo.area_m2 - pred.area_m2) <= 1e-9 * pred.area_m2)
    out["outward_normal_is_minus_x"] = bool(np.allclose(normals, [-1.0, 0.0, 0.0], atol=1e-9))
    write_json(results_dir(case_dir) / f"mesh_check_{subdir}.json", out)
    if not out["matches_prediction"]:
        raise WorkflowError(f"built source patch ({geo.n_faces} faces, {geo.area_m2:.6e} m^2) "
                            f"differs from the prediction ({pred.n_faces}, {pred.area_m2:.6e})")
    if not out["outward_normal_is_minus_x"]:
        raise WorkflowError("source patch normals are not all -x")
    return out


def background_label(ratio: float) -> str:
    """``cfd_bg1e-03`` style directory name for a sensitivity copy."""
    return "cfd_bg" + f"{ratio:.0e}".replace("+", "")


# --------------------------------------------------------------------------- #
# cfd-report
# --------------------------------------------------------------------------- #

def _latest(case: Path) -> str:
    from plumetools.foamio.times import latest_time
    t = latest_time(case)
    if float(t) <= 0.0:
        raise WorkflowError(f"{case} has no solution time beyond 0")
    return t


def _latest_or_none(case: Path):
    from plumetools.foamio.times import latest_time
    try:
        return latest_time(case)
    except Exception:
        return None


def _postprocessing_series(case: Path, name: str) -> np.ndarray | None:
    """Concatenate a function object's .dat files over restarts (by time)."""
    base = Path(case) / "postProcessing" / name
    if not base.is_dir():
        return None
    rows = []
    for f in sorted(base.rglob("*.dat")):
        for line in f.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            rows.append([float(v) for v in line.replace("(", " ").replace(")", " ").split()])
    if not rows:
        return None
    width = min(len(r) for r in rows)
    arr = np.asarray([r[:width] for r in rows])
    _, keep = np.unique(arr[:, 0], return_index=True)
    return arr[np.sort(keep)]


def cfd_report(case_dir: Path, time: str | None = None, subdir: str = "cfd") -> dict:
    """Mass/momentum/energy budgets and stationarity of the converged precursor."""
    from .budgets import cfd_budgets
    from .meshgeom import read_mesh
    case_dir = Path(case_dir)
    cfg = load_case(case_dir)
    source = config_mod.build_source(cfg)
    cfd = case_dir / subdir
    time = time or _latest(cfd)
    mesh = read_mesh(cfd)
    series = {k: _postprocessing_series(cfd, k) for k in (
        "boundaryMassFlux", "exitPlaneMassFlux", "farFieldMassFlux", "domainMass")}
    doc = cfd_budgets(cfd, time, mesh, source, series, cfg)
    doc["time"] = time
    doc["label"] = cfg.get("label", cfg.profile)
    doc["run_record"] = _run_record(cfd, cfg)
    write_json(results_dir(case_dir) / f"cfd_report_{subdir}.json", doc)
    return doc


def _run_record(case: Path, cfg, extra: dict | None = None) -> dict:
    """Configuration digest, mesh identity, repository and toolchain state, now."""
    mesh_files = {name: provenance.sha256(Path(case) / "constant" / "polyMesh" / name)
                  for name in ("points", "faces", "owner", "neighbour", "boundary")}
    rec = {"config_digest": cfg.digest, "profile": cfg.profile,
           "mesh_sha256": mesh_files,
           "provenance": provenance.collect(Path(__file__).resolve().parents[2])}
    rec.update(extra or {})
    return rec


# --------------------------------------------------------------------------- #
# breakdown
# --------------------------------------------------------------------------- #

def breakdown(case_dir: Path, time: str | None = None) -> dict:
    from . import breakdown as bd
    from .meshgeom import read_mesh
    case_dir = Path(case_dir)
    cfg = load_case(case_dir)
    source = config_mod.build_source(cfg)
    gas = source.gas
    cfd = case_dir / "cfd"
    time = time or _latest(cfd)
    mesh = read_mesh(cfd)
    graph = mesh.cell_graph()
    state = bd.read_state(cfd, time, gas)
    grads = bd.gradients(mesh, state, graph)
    conv = cfg.get("breakdown.mean_free_path_convention", "vhs")
    ind = bd.indicators(state, grads, gas, conv)
    thresholds = cfg.require("breakdown.thresholds")
    criterion = tuple(cfg.get("breakdown.criterion", ["Kn_GLL", "B_exp"]))
    mask = bd.continuum_mask(ind, thresholds, criterion)
    reg = bd.source_connected(mesh, mask, mesh_mod.PATCHES["source"])
    lip = bd.lip_distance(mesh, source.radius_m)
    near_lip = lip < float(cfg.get("breakdown.lip_exclusion_over_D", 0.05)) * source.diameter_m
    doc = {
        "time": time, "label": cfg.get("label", cfg.profile),
        "definitions": {
            "lambda": f"Gas.mean_free_path, convention {conv!r}",
            "Kn_rho": "lambda |grad rho| / rho", "Kn_T": "lambda |grad T| / T",
            "Kn_U": "lambda |grad |U|| / max(|U|, c)",
            "Kn_GLL": "max(Kn_rho, Kn_T, Kn_U)",
            "B_exp": "-(lambda / c_bar) U . grad(rho) / rho  (signed; >0 expansion)",
            "gradient": "weighted least squares over face neighbours (OpenFOAM leastSquares); "
                        "one-sided at boundaries; unaveraged CFD cell values",
            "criterion": list(criterion), "thresholds": thresholds,
            "threshold_status": "sensitivity-study settings, not validity boundaries",
        },
        "summary": bd.summarise(mesh, ind, reg["retained"], near_lip),
        "regions": {k: v for k, v in reg.items() if k not in ("retained", "labels")},
        "threshold_sensitivity": bd.threshold_sensitivity(
            mesh, ind, thresholds, criterion,
            cfg.get("breakdown.sensitivity_factors", [0.5, 1.0, 2.0]),
            mesh_mod.PATCHES["source"]),
    }
    # Background-floor sensitivity of the retained core: same mesh, other floor.
    sens = []
    for other in sorted(case_dir.glob("cfd_bg*")):
        t_other = _latest_or_none(other)
        if t_other is None or float(t_other) <= 0:
            sens.append({"case": other.name, "available": False,
                         "reason": "no solution beyond t = 0"})
            continue
        st2 = bd.read_state(other, t_other, gas)
        r = reg["retained"]
        ind2 = bd.indicators(st2, bd.gradients(mesh, st2, graph), gas, conv)
        reg2 = bd.source_connected(mesh, bd.continuum_mask(ind2, thresholds, criterion),
                                   mesh_mod.PATCHES["source"])

        def rel(a, b):
            return float(np.nanmax(np.abs(a[r] / b[r] - 1.0))) if r.any() else None
        sens.append({
            "case": other.name, "time": t_other, "available": True,
            "max_rel_diff_in_retained": {"rho": rel(st2.rho, state.rho),
                                         "T": rel(st2.T, state.T),
                                         "speed": rel(st2.speed, state.speed)},
            "retained_cells": int(reg2["retained"].sum()),
            "retained_cells_reference": int(r.sum()),
            "retained_symmetric_difference_cells": int((reg2["retained"] ^ r).sum())})
    doc["background_sensitivity"] = sens
    fields = {k: ind[k] for k in bd.INDICATORS}
    fields["lambda"] = ind["lambda"]
    fields["continuumMask"] = mask.astype(float)
    fields["retainedRegion"] = reg["retained"].astype(float)
    bd.write_fields(cfd, time, mesh, fields)
    np.savez_compressed(results_dir(case_dir) / "breakdown.npz", retained=reg["retained"],
                        mask=mask, near_lip=near_lip, time=time,
                        valid=ind["valid"], **{k: ind[k] for k in bd.INDICATORS})
    write_json(results_dir(case_dir) / "breakdown.json", doc)
    return doc


# --------------------------------------------------------------------------- #
# interface
# --------------------------------------------------------------------------- #

def interface(case_dir: Path, time: str | None = None, allow_unsuitable: bool = False) -> dict:
    from . import breakdown as bd
    from . import interface as itf
    from . import mapping, resources
    from .cfd import background_state
    from .meshgeom import read_mesh
    case_dir = Path(case_dir)
    cfg = load_case(case_dir)
    source = config_mod.build_source(cfg)
    gas = source.gas
    mode = cfg.get("interface.mode", "automatic")
    res = results_dir(case_dir)
    if mode == "none":
        doc = {"mode": "none", "status": itf.STATUS_FULL_DSMC,
               "label": cfg.get("label", cfg.profile)}
        write_json(res / "interface.json", doc)
        return doc
    cfd = case_dir / "cfd"
    mesh = read_mesh(cfd)
    bdz = res / "breakdown.npz"
    if not bdz.is_file():
        raise WorkflowError("run the 'breakdown' step first")
    with np.load(bdz) as z:
        time = time or str(z["time"])
        ind = {k: z[k] for k in bd.INDICATORS}
        ind["valid"] = z["valid"]
        retained = z["retained"]
    state = bd.read_state(cfd, time, gas)
    grads = bd.gradients(mesh, state)
    if mode == "manual":
        region = itf.manual_region(mesh, cfg.require("interface.manual"), source.diameter_m)
        retained = itf.connected_to_source(mesh, region)
    surf = itf.extract(mesh, retained)
    doc = itf.assess(mesh, surf, state, grads, ind, cfg, source, mode=mode)
    arrays = doc.pop("_face_arrays", None)
    doc["time"] = time
    doc["label"] = cfg.get("label", cfg.profile)
    doc["allow_unsuitable"] = bool(allow_unsuitable or cfg.get("interface.allow_unsuitable"))
    if arrays is not None:
        fields = mapping.canonical(surf, arrays, gas)
        p_face = itf.face_values(mesh, surf, state.p, cfg.get("interface.face_state", "linear"))
        doc["state_consistency"] = mapping.consistency(fields, gas, p_face)
        mapping.save(res / "interface_fields.npz", fields)
        itf.write_vtk(res / "interface.vtk", mesh, surf.faces, {
            "Mn": arrays["Mn"], "Un": arrays["Un"], "rho": arrays["rho"], "T": arrays["T"],
            "Kn_GLL_cfd_side": ind["Kn_GLL"][surf.cfd_cell],
            "Kn_GLL_dsmc_side": ind["Kn_GLL"][surf.dsmc_cell],
            "B_exp_cfd_side": ind["B_exp"][surf.cfd_cell],
            "normal": surf.normal, "U": arrays["U"]})
        doc["artifacts"] = {"surface_vtk": str(res / "interface.vtk"),
                            "canonical_fields": str(res / "interface_fields.npz")}
        doc["geometry"] = {
            "centroid_m": (surf.centre * surf.area[:, None]).sum(0) / surf.area.sum(),
            "x_range_m": [float(surf.centre[:, 0].min()), float(surf.centre[:, 0].max())],
            "max_radius_m": float(np.hypot(surf.centre[:, 1], surf.centre[:, 2]).max()),
        }
    keep = ~retained
    floor = 3.0 * background_state(cfg, source)["density_kg_m3"]
    use = keep & (state.rho > floor) & np.isfinite(state.rho)
    if use.any():
        doc["dsmc_estimate"] = resources.dsmc_estimate(
            cfg, source, number_density=state.rho[use] / gas.mass_kg, temperature=state.T[use],
            speed=state.speed[use], volumes=mesh.cell_volumes[use],
            sizes=mesh.cell_size()[use], length_scale_m=float(mesh.points[:, 0].max()),
            basis="CFD state on the DSMC cells, excluding cells within 3x of the "
                  "background floor (the DSMC there is near-vacuum)")
        doc["dsmc_estimate"]["cells_total_in_dsmc_region"] = int(keep.sum())
    np.save(res / "retained.npy", retained)
    write_json(res / "interface.json", doc)
    return doc


# --------------------------------------------------------------------------- #
# dsmc-setup / dsmc-fields
# --------------------------------------------------------------------------- #

def _add_zero_face_patch(case: Path, name: str) -> None:
    """Append ``name`` (type patch, 0 faces) so subsetMesh does not make it ``empty``."""
    from plumetools.mesh.boundary import read_boundary
    patches = read_boundary(case)
    if name in patches:
        return
    path = Path(case) / "constant" / "polyMesh" / "boundary"
    n_faces = max(p.start_face + p.n_faces for p in patches.values())
    text = path.read_text()
    count = len(patches)
    marker = f"\n{count}\n("
    if marker not in text:
        raise WorkflowError(f"cannot find the patch count {count} in {path}")
    head, _, tail = text.rpartition(")")
    head = head.replace(marker, f"\n{count + 1}\n(", 1)
    entry = (f"    {name}\n    {{\n        type            patch;\n"
             f"        nFaces          0;\n        startFace       {n_faces};\n    }}\n")
    path.write_text(head + entry + ")" + tail)


def dsmc_setup(case_dir: Path, np_ranks: int | None = None, **kwargs) -> dict:
    import shutil
    from . import breakdown as bd
    from . import dsmc, resources
    from .cfd import background_state
    from .interface import STATUS_FULL_DSMC, STATUS_SUITABLE
    from .meshgeom import read_mesh
    np_ranks = kwargs.get("np", np_ranks)
    case_dir = Path(case_dir)
    cfg = load_case(case_dir)
    source = config_mod.build_source(cfg)
    res = results_dir(case_dir)
    itf_path = res / "interface.json"
    mode = cfg.get("interface.mode", "automatic")
    if itf_path.is_file():
        itf = json.loads(itf_path.read_text())
    elif mode == "none":
        itf = {"mode": "none", "status": STATUS_FULL_DSMC}
    else:
        raise WorkflowError("run the 'interface' step first")
    status = itf["status"]
    if status not in (STATUS_SUITABLE, STATUS_FULL_DSMC) and not itf.get("allow_unsuitable"):
        raise WorkflowError(f"interface status is {status!r} ({itf.get('reasons')}); "
                            "re-run 'interface --allow-unsuitable' to build it as a "
                            "labelled machinery smoke test, or choose full DSMC")
    cfd = case_dir / "cfd"
    dsmc_dir = case_dir / "dsmc"
    mesh = read_mesh(cfd)
    retained = (np.load(res / "retained.npy")
                if status != STATUS_FULL_DSMC and (res / "retained.npy").is_file()
                else np.zeros(mesh.n_cells, bool))
    keep = ~retained
    has_interface = bool(retained.any())
    src_sl = mesh.patch_slice(mesh_mod.PATCHES["source"])
    src_left = int((~retained[mesh.owner[src_sl]]).sum())

    t = _latest_or_none(cfd)
    if t is not None and float(t) > 0:
        state = bd.read_state(cfd, t, source.gas)
        floor = 3.0 * background_state(cfg, source)["density_kg_m3"]
        use = keep & (state.rho > floor)
        est = resources.dsmc_estimate(
            cfg, source, number_density=state.rho[use] / source.gas.mass_kg,
            temperature=state.T[use], speed=state.speed[use],
            volumes=mesh.cell_volumes[use], sizes=mesh.cell_size()[use],
            length_scale_m=float(mesh.points[:, 0].max()),
            basis="CFD state on the DSMC cells")
    else:
        est = resources.assess_configuration(cfg, source,
                                             mesh_mod.plan(cfg, source))["full_dsmc"]
    run = dsmc.plan_run(cfg, source, est, float(mesh.points[:, 0].max()))
    ranks = int(np_ranks or cfg.get("resources.max_ranks", 1))
    patches = dsmc.inflow_patches(has_interface, src_left)
    dsmc_dir.mkdir(parents=True, exist_ok=True)
    poly = dsmc_dir / "constant" / "polyMesh"
    if not (poly / "owner").is_file():
        shutil.copytree(cfd / "constant" / "polyMesh", poly, dirs_exist_ok=True)
        shutil.rmtree(poly / "sets", ignore_errors=True)
        if has_interface:
            _add_zero_face_patch(dsmc_dir, dsmc.INTERFACE_PATCH)
            dsmc.write_keep_set(dsmc_dir, np.nonzero(keep)[0])
    dsmc.write_case(dsmc_dir, cfg, source, run, patches, ranks)
    doc = {"label": cfg.get("label", cfg.profile), "interface_status": status,
           "smoke_only": status not in (STATUS_SUITABLE, STATUS_FULL_DSMC),
           "inflow_patches": patches, "dsmc_cells": int(keep.sum()),
           "source_faces_injected_by_dsmc": src_left, "run": run.describe(), "ranks": ranks}
    write_json(res / "dsmc_setup.json", doc)
    return doc


def dsmc_fields(case_dir: Path) -> dict:
    """After subsetMesh: write 0/ with the mapped interface and residual source."""
    from plumetools.mesh.boundary import read_boundary
    from . import dsmc, mapping
    from .meshgeom import read_mesh
    case_dir = Path(case_dir)
    cfg = load_case(case_dir)
    source = config_mod.build_source(cfg)
    res = results_dir(case_dir)
    dsmc_dir = case_dir / "dsmc"
    patches = read_boundary(dsmc_dir)
    setup_doc = json.loads((res / "dsmc_setup.json").read_text())
    inflow = {}
    doc = {"label": cfg.get("label", cfg.profile)}
    dmesh = read_mesh(dsmc_dir)
    cfd_mesh = read_mesh(case_dir / "cfd")
    cell_map = mapping.match_cells(dmesh.cell_centres, cfd_mesh.cell_centres,
                                   cfd_mesh.cell_size().max(axis=1))
    np.save(res / "dsmc_to_cfd_cells.npy", cell_map)
    doc["cell_map"] = {"dsmc_cells": int(cell_map.size), "verified": "centre coincidence"}
    if dsmc.INTERFACE_PATCH in setup_doc["inflow_patches"]:
        info = patches.get(dsmc.INTERFACE_PATCH)
        if info is None or info.type != "patch" or info.n_faces == 0:
            raise WorkflowError(f"DSMC mesh interface patch is {info}; expected a non-empty "
                                "'patch' (an 'empty' patch means it was not pre-created)")
        fields = mapping.load(res / "interface_fields.npz")
        match = mapping.match_interface(dmesh, fields, dsmc.INTERFACE_PATCH)
        inflow[dsmc.INTERFACE_PATCH] = mapping.inflow_from_fields(fields, match)
        doc["mapping"] = {k: v for k, v in match.items()
                          if k not in ("rows", "dsmc_faces", "dsmc_area_vectors")}
        doc["mapping_flux_change"] = mapping.flux_change(fields, match, source.gas)
        np.save(res / "interface_match_rows.npy", match["rows"])
    src = mesh_mod.PATCHES["source"]
    if src in setup_doc["inflow_patches"]:
        inflow[src] = mapping.source_inflow(patches[src].n_faces, source)
        doc["source_faces"] = patches[src].n_faces
    dsmc.write_fields(dsmc_dir, patches, source.gas.name, inflow, source)
    doc["patches"] = {k: {"type": v.type, "n_faces": v.n_faces} for k, v in patches.items()}
    write_json(res / "dsmc_fields.json", doc)
    return doc
