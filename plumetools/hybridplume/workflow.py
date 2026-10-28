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
``dsmc-post``       moments, measured fluxes, budgets, comparisons
``crossings``       exact exits over a per-step parcel window after the run
``overlap``         hybrid vs a full-DSMC case on the same mesh, cell by cell
``reference``       the collisionless reference at the report stations
``report``          collect everything into results/summary.json + plots
``hpc``             resumable batch scripts sized from the resource estimates
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


def _time_dirs(case: Path) -> list:
    from plumetools.foamio.times import time_names
    return sorted(time_names(case), key=float)


def _log_series(log: Path, pattern: str) -> np.ndarray:
    import re
    if not log.is_file():
        return np.empty(0)
    return np.array([float(v) for v in re.findall(pattern, log.read_text(errors="replace"))])


def dsmc_post(case_dir: Path, time: str | None = None) -> dict:
    """Sampled DSMC solution, measured injection/return, budgets, comparisons."""
    from . import collisionless, halfrange, mapping, particles, post, resources
    from .meshgeom import read_mesh
    case_dir = Path(case_dir)
    cfg = load_case(case_dir)
    source = config_mod.build_source(cfg)
    gas = source.gas
    res = results_dir(case_dir)
    dsmc_dir = case_dir / "dsmc"
    setup_doc = json.loads((res / "dsmc_setup.json").read_text())
    run = setup_doc["run"]
    weight, dt = float(run["weight"]), float(run["dt_s"])
    t_end, t_avg = float(run["end_time_s"]), float(run["average_start_s"])
    times = [t for t in _time_dirs(dsmc_dir) if float(t) <= t_end * (1 + 1e-9)]
    if not times or float(times[-1]) <= t_avg:
        raise WorkflowError(f"{dsmc_dir} has not reached its averaging window "
                            f"(start {t_avg:.4e} s)")
    t_final = time or times[-1]
    dmesh = read_mesh(dsmc_dir)
    sizes = dmesh.cell_size()
    doc = {"label": cfg.get("label", cfg.profile), "time": t_final,
           "smoke_only": setup_doc.get("smoke_only", False),
           "interface_status": setup_doc.get("interface_status")}

    # -- fieldAverage moments ------------------------------------------------
    fa = post.dsmc_fieldaverage(dsmc_dir / t_final, gas)
    occ = fa["parcels_per_step"]
    doc["sampling"] = {
        "average_start_s": t_avg, "end_s": float(t_final),
        "steps_averaged": int(round((float(t_final) - t_avg) / dt)),
        "parcels_per_cell_per_step": {"median": float(np.median(occ)),
                                      "volume_fraction_below_1": float(
                                          dmesh.cell_volumes[occ < 1].sum()
                                          / dmesh.cell_volumes.sum())},
        "weight_check": {"max_rel": float(np.nanmax(np.abs(
            fa["n"][occ > 0] * dmesh.cell_volumes[occ > 0] / occ[occ > 0] / weight - 1)))},
    }
    if "n_cv_per_step" in fa:
        cv = fa["n_cv_per_step"]
        doc["sampling"]["n_cv_per_step_median"] = float(np.nanmedian(cv[occ > 1]))

    # -- stationarity from the log ------------------------------------------
    log = dsmc_dir / "log.dsmcFoam"
    nmol = _log_series(log, r"Number of molecules\s*=\s*([0-9.eE+]+)")
    if nmol.size > 10:
        q = max(2, nmol.size // 4)
        doc["stationarity"] = {"parcels_last_quarter_rel_change":
                               float((nmol[-1] - nmol[-q]) / nmol[-1]),
                               "parcels_final": float(nmol[-1])}
    # -- measured vs expected injection ---------------------------------------
    inserted = particles.inserted_from_log(log)
    n_window = int(round((float(t_final) - t_avg) / dt))
    expected = 0.0
    if (res / "interface_fields.npz").is_file() and "interface" in setup_doc["inflow_patches"]:
        f = mapping.load(res / "interface_fields.npz")
        expected += float(np.dot(halfrange.incoming_flux(
            f["numberDensity"], f["velocity"], f["translationalTemperature"], f["faceNormal"],
            gas.gas_constant)["mass"], f["faceArea"]))
    if "source" in setup_doc["inflow_patches"]:
        sl = dmesh.patch_slice(mesh_mod.PATCHES["source"])
        a_src = float(np.linalg.norm(dmesh.face_areas[sl], axis=1).sum())
        expected += halfrange.gross_source_number_flux(
            source.number_density, source.speed_ratio, source.thermal_speed) * a_src
    if inserted.size >= n_window > 0:
        measured = float(inserted[-n_window:].mean()) * weight / dt
        doc["injection"] = {
            "expected_gross_number_rate_per_s": expected,
            "measured_gross_number_rate_per_s": measured,
            "relative_difference": measured / expected - 1.0 if expected > 0 else None,
            "expected_counting_stderr_rel": 1.0 / np.sqrt(inserted[-n_window:].sum()),
            "note": "gross half-range injection, NOT the net continuum mass transport"}

    # -- tensor moments from parcel snapshots, in independent blocks ----------
    snaps = [t for t in times if float(t) >= t_avg - 1e-12 * t_end]
    n_blocks = int(cfg.get("dsmc.uncertainty_blocks", 4))
    blocks = np.array_split(np.array(snaps), max(1, min(n_blocks, len(snaps))))
    total = particles.MomentAccumulator(dmesh.n_cells)
    per_block = []
    for blk in blocks:
        acc = particles.MomentAccumulator(dmesh.n_cells)
        for t in blk:
            cl = particles.read_cloud(dsmc_dir / t)
            acc.add(cl, weight)
            total.add(cl, weight)
        per_block.append(acc.fields(dmesh.cell_volumes, gas))
    lag = total.fields(dmesh.cell_volumes, gas)
    doc["snapshot_moments"] = {"snapshots": lag["snapshots"], "blocks": len(blocks)}

    # -- profiles, on the full mesh: DSMC cells scattered into CFD numbering ----
    # In a hybrid case the axis may lie inside the retained CFD region, so the
    # profile is taken on the WHOLE CFD mesh: "hybrid" = CFD where retained,
    # DSMC elsewhere; the DSMC-only and CFD-only columns are kept beside it.
    dsmc_fields_ = {"n_fa": fa["n"], "T_fa": fa["T"], "U_fa": fa["U"], "n_lag": lag["n"],
                    "T_lag": lag["T"], "Tpar_lag": lag["T_parallel"],
                    "Tperp_lag": lag["T_perpendicular"]}
    for k, name in (("n", "n_lag_se"), ("T", "T_lag_se")):
        dsmc_fields_[name] = particles.block_statistics([b[k] for b in per_block])["stderr"]
    cell_map = res / "dsmc_to_cfd_cells.npy"
    t_cfd = _latest_or_none(case_dir / "cfd") if (case_dir / "cfd").is_dir() else None
    hybrid = cell_map.is_file() and t_cfd is not None and float(t_cfd) > 0
    if hybrid:
        from . import breakdown as bd
        cmesh = read_mesh(case_dir / "cfd")
        cm = np.load(cell_map)
        centres, volumes, n_all = cmesh.cell_centres, cmesh.cell_volumes, cmesh.n_cells
        fields = {}
        for k, v in dsmc_fields_.items():
            full = np.full((n_all,) + np.shape(v)[1:], np.nan)
            full[cm] = v
            fields[k] = full
        st = bd.read_state(case_dir / "cfd", t_cfd, gas)
        fields["n_cfd"] = st.rho / gas.mass_kg
        fields["T_cfd"] = st.T
        fields["U_cfd"] = st.U
        retained_c = (np.load(res / "retained.npy") if (res / "retained.npy").is_file()
                      else np.zeros(n_all, bool))
        fields["n_hyb"] = np.where(retained_c, fields["n_cfd"], fields["n_fa"])
        fields["T_hyb"] = np.where(retained_c, fields["T_cfd"], fields["T_fa"])
        fields["U_hyb"] = np.where(retained_c[:, None], fields["U_cfd"], fields["U_fa"])
        from .breakdown import write_fields as _wf
        _wf(case_dir / "cfd", t_cfd, cmesh, {
            "hybrid_n_over_ne": fields["n_hyb"] / source.number_density,
            "hybrid_T_over_Te": fields["T_hyb"] / source.temperature_K,
            "hybrid_Ux_over_ae": fields["U_hyb"][:, 0] / source.thermal_speed,
            "dsmc_n_over_ne": fields["n_fa"] / source.number_density})
        h = float(np.min(cmesh.cell_size()[:, 1]))
    else:
        centres, volumes, n_all = dmesh.cell_centres, dmesh.cell_volumes, dmesh.n_cells
        fields = dict(dsmc_fields_)
        h = float(np.min(sizes[:, 1]))
    axis = post.axis_cells(centres, h * 0.51)
    ref = collisionless.reference(centres[axis], source)
    for k in ("n", "T", "U"):
        full = np.full((n_all,) + np.shape(ref[k])[1:], np.nan)
        full[axis] = ref[k]
        fields[f"{k}_ref"] = full
    prof = post.line_profile(centres, volumes, axis, fields)
    norm = {"x_over_D": prof["x"] / source.diameter_m}
    for k, v in prof.items():
        if k == "x":
            continue
        if k.startswith("n_"):
            norm[k + "_over_ne"] = v / source.number_density
        elif k.startswith(("T", "Tpar", "Tperp")):
            norm[k + "_over_Te"] = v / source.temperature_K
        elif k.startswith("U_"):
            norm[k + "_over_ae"] = v / source.thermal_speed
    post.write_csv(res / "centreline.csv", norm)
    cmp = {
        "dsmc_fieldaverage_vs_collisionless": {
            "n": post.compare(prof["n_fa"], prof["n_ref"]),
            "T": post.compare(prof["T_fa"], prof["T_ref"]),
            "Ux": post.compare(prof["U_fa"][:, 0], prof["U_ref"][:, 0])},
        "dsmc_snapshot_vs_fieldaverage": {
            "n": post.compare(prof["n_lag"], prof["n_fa"]),
            "T": post.compare(prof["T_lag"], prof["T_fa"])},
    }
    if "n_cfd" in prof:
        cmp["dsmc_vs_cfd_precursor"] = {
            "n": post.compare(prof["n_fa"], prof["n_cfd"]),
            "T": post.compare(prof["T_fa"], prof["T_cfd"]),
            "Ux": post.compare(prof["U_fa"][:, 0], prof["U_cfd"][:, 0])}
        cmp["hybrid_vs_collisionless"] = {
            "n": post.compare(prof["n_hyb"], prof["n_ref"]),
            "T": post.compare(prof["T_hyb"], prof["T_ref"])}
    cmp["note"] = ("the collisionless reference is exact only in the collisionless limit; a "
                   "collisional difference is physics, not numerical error")
    doc["centreline_comparison"] = cmp
    doc["centreline_snapshot_stderr_rel_median"] = {
        "n": float(np.nanmedian(prof["n_lag_se"] / prof["n_lag"])),
        "T": float(np.nanmedian(prof["T_lag_se"] / prof["T_lag"]))}

    # -- azimuthal structure (unaveraged spread reported) -----------------------
    d = source.diameter_m
    good = occ >= 5
    doc["azimuthal"] = {
        name: {k: v for k, v in post.azimuthal_statistics(
            dmesh.cell_centres, vals, 0.25 * d, 0.1 * d, mask=good).items()
               if k in ("max_rel_std", "median_rel_std")}
        for name, vals in (("n", fa["n"]), ("T", fa["T"]))}
    doc["azimuthal"]["note"] = ("relative std over azimuth within (x, r) bins of cells with "
                                ">= 5 parcels/step; includes sampling noise and the "
                                "Cartesian staircase")
    # -- off-axis rings and the spherical arc, azimuthally averaged for reporting -
    # Rings: cells whose radius is within half a cell of r0, grouped by x; the
    # unaveraged azimuthal spread of the same cells is in "azimuthal" above.
    c = centres
    r_cell = np.hypot(c[:, 1], c[:, 2])
    off = {}
    for r0 in cfg.get("post.off_axis_r_over_D", [0.25, 0.5, 0.75]):
        ring = np.abs(r_cell - float(r0) * d) <= 0.5 * h
        if ring.sum() < 4:
            continue
        ref_r = collisionless.reference(c[ring], source)
        rf = {"n": np.full(n_all, np.nan), "T": np.full(n_all, np.nan)}
        rf["n"][ring], rf["T"][ring] = ref_r["n"], ref_r["T"]
        pr = post.line_profile(c, volumes, ring, {
            "n": fields.get("n_hyb", fields["n_fa"]), "T": fields.get("T_hyb", fields["T_fa"]),
            "n_ref": rf["n"], "T_ref": rf["T"]})
        post.write_csv(res / f"offaxis_r{float(r0):g}D.csv", {
            "x_over_D": pr["x"] / d, "n_over_ne": pr["n"] / source.number_density,
            "T_over_Te": pr["T"] / source.temperature_K,
            "n_ref_over_ne": pr["n_ref"] / source.number_density,
            "T_ref_over_Te": pr["T_ref"] / source.temperature_K})
        off[f"r{float(r0):g}D"] = {"n_vs_collisionless": post.compare(pr["n"], pr["n_ref"]),
                                   "T_vs_collisionless": post.compare(pr["T"], pr["T_ref"])}
    doc["off_axis"] = off
    r_arc = float(cfg.get("post.arc_r_over_D", 10.0)) * d
    reach = min(float(c[:, 0].max()), float(r_cell.max()))
    if reach >= r_arc:
        shell = np.abs(np.linalg.norm(c, axis=1) - r_arc) <= 0.5 * h
        ang = np.degrees(np.arctan2(r_cell[shell], c[shell, 0]))
        doc["arc"] = {"r_over_D": r_arc / d, "cells": int(shell.sum()),
                      "angle_range_deg": [float(ang.min()), float(ang.max())]}
    else:
        doc["arc"] = {"r_over_D": r_arc / d, "available": False,
                      "reason": f"domain reaches {reach / d:.2f} D on the axis/radius; "
                                "the reference arc is in reference_arc_r10D.csv"}

    # transverse momentum symmetry: net y, z momentum of the sampled field
    mom = (fa["rho"][:, None] * fa["U"] * dmesh.cell_volumes[:, None])
    mom = np.where(np.isfinite(mom), mom, 0.0)
    doc["transverse_momentum_over_axial"] = (mom[:, 1:].sum(0) / mom[:, 0].sum()).tolist()

    # -- returns through the interface (estimate from snapshots) --------------
    if "interface" in setup_doc["inflow_patches"]:
        sl = dmesh.patch_slice("interface")
        faces = np.arange(sl.start, sl.stop)
        a = np.linalg.norm(dmesh.face_areas[faces], axis=1)
        n_out = dmesh.face_areas[faces] / a[:, None]
        rates = []
        for t in snaps:
            est = particles.cell_half_range_estimate(
                particles.read_cloud(dsmc_dir / t), dmesh.owner[faces], n_out, a,
                dmesh.cell_volumes, weight, gas.mass_kg)
            rates.append(est["mass_kg_s"])
        itf = json.loads((res / "interface.json").read_text())
        exp_ret = itf["fluxes"]["maxwellian_reverse"]["mass"]
        doc["return_flux"] = {
            "estimated_mass_kg_s": float(np.mean(rates)),
            "estimate_stderr_kg_s": float(np.std(rates, ddof=1) / np.sqrt(len(rates)))
            if len(rates) > 1 else None,
            "expected_from_cfd_maxwellian_kg_s": exp_ret,
            "forward_injection_expected_kg_s": itf["fluxes"]["maxwellian_forward"]["mass"],
            "method": "cell half-range estimator from snapshots (an ESTIMATE); run the "
                      "'crossings' stage for exact per-step crossings"}
    cw = res / "crossing_window.json"
    if cw.is_file():
        doc["crossing_window"] = json.loads(cw.read_text())

    # -- resources measured -------------------------------------------------
    doc["measured_performance"] = resources.measured_rates(
        log, float(particles.read_cloud(dsmc_dir / t_final).size), dmesh.n_cells,
        int(setup_doc["ranks"]))
    doc["rank_imbalance"] = resources.rank_imbalance(dsmc_dir, t_final)
    doc["run_record"] = _run_record(dsmc_dir, cfg, extra={
        "averaging_interval_s": [t_avg, float(t_final)], "particle_weight": weight,
        "random_seed": "not configurable in standard dsmcFoam: each rank's Random is "
                       "seeded deterministically, so a rerun on the same decomposition "
                       "repeats; independence comes from blocks, not seeds"})
    write_json(res / "dsmc_post.json", doc)
    return doc


def crossings(case_dir: Path, **_) -> dict:
    """Exact exits over the per-step window written by the Allrun 'crossings' stage."""
    from . import particles
    from .meshgeom import read_mesh
    case_dir = Path(case_dir)
    cfg = load_case(case_dir)
    source = config_mod.build_source(cfg)
    res = results_dir(case_dir)
    setup_doc = json.loads((res / "dsmc_setup.json").read_text())
    dt = float(setup_doc["run"]["dt_s"])
    t_end = float(setup_doc["run"]["end_time_s"])
    dsmc_dir = case_dir / "dsmc"
    times = [t for t in _time_dirs(dsmc_dir) if float(t) >= t_end * (1 - 1e-9)]
    # Measured 2026-10-28: every dsmcFoam restart with plumeFieldInflow injects a
    # burst (~+1300 parcels over ~10 steps on smoke-kn1em3, identical inflow
    # fields, with or without per-step writes). The first steps after the
    # restart are therefore discarded; the insertion series is reported.
    discard = int(cfg.get("dsmc.crossing_discard_steps", 12))
    times = times[discard:]
    if len(times) < 3:
        raise WorkflowError("too few per-step snapshots after the post-restart discard "
                            f"({discard}); run a longer 'crossings' window")
    snaps = [(float(t), particles.read_cloud(dsmc_dir / t)) for t in times]
    cfd_mesh = read_mesh(case_dir / "cfd")
    retained = (np.load(res / "retained.npy") if (res / "retained.npy").is_file()
                else np.zeros(cfd_mesh.n_cells, bool))
    plan = mesh_mod.plan(cfg, source)
    in_region = particles.lattice_lookup(cfd_mesh.cell_centres, retained,
                                         tuple(plan.edges(a) for a in range(3)))
    lo = cfd_mesh.points.min(axis=0)
    hi = cfd_mesh.points.max(axis=0)
    doc = particles.crossing_window(snaps, dt, float(setup_doc["run"]["weight"]), source.gas,
                                    (lo, hi), in_region, source.radius_m)
    inserted = particles.inserted_from_log(dsmc_dir / "log.dsmcFoam.crossings")
    doc["discarded_post_restart_steps"] = discard
    if inserted.size:
        w = float(setup_doc["run"]["weight"])
        doc["insertion_series_after_restart"] = inserted.tolist()
        doc["injected_from_log"] = {
            "steps": int(min(inserted.size, doc["steps"])),
            "mass_kg_s": float(inserted[-doc["steps"]:].sum() * w * source.gas.mass_kg
                               / doc["duration_s"])}
        itf_path = res / "interface.json"
        if itf_path.is_file():
            itf = json.loads(itf_path.read_text())
            if "fluxes" in itf:
                doc["expected_from_cfd_maxwellian_kg_s"] = {
                    "forward": itf["fluxes"]["maxwellian_forward"]["mass"],
                    "reverse": itf["fluxes"]["maxwellian_reverse"]["mass"]}
    write_json(res / "crossing_window.json", doc)
    return doc


def reference(case_dir: Path, **_) -> dict:
    """The collisionless reference at the report stations, with quadrature checks."""
    from . import collisionless, post
    from plumetools.cai2012 import analytical
    case_dir = Path(case_dir)
    cfg = load_case(case_dir)
    source = config_mod.build_source(cfg)
    pts = collisionless.centreline_and_angular_points(source)
    res = results_dir(case_dir)
    out = {"label": cfg.get("label", cfg.profile),
           "status": "exact only in the collisionless limit (Kn -> infinity)"}
    for name, p in (("centreline", pts["centreline"]), ("arc_r10D", pts["arc"])):
        ref = collisionless.converged(p, source)
        cf = analytical.moments(p[:, 0], p[:, 1], p[:, 2], source.radius_m, source.speed_ratio)
        cols = {"x": p[:, 0], "y": p[:, 1], "z": p[:, 2], "n_over_ne": ref["n_over_ne"],
                "U_over_ae": ref["U_over_ae"], "T_over_Te": ref["T_over_Te"],
                "Tdir_over_Te": ref["T_directional"] / source.temperature_K,
                "quadrature_rel_change": ref["quadrature_rel_change"]}
        if name == "arc_r10D":
            cols["angle_deg"] = pts["arc_angle_deg"]
        post.write_csv(res / f"reference_{name}.csv", cols)
        out[name] = {
            "points": int(p.shape[0]),
            "max_quadrature_rel_change": float(np.max(ref["quadrature_rel_change"])),
            "vs_cai2012_closed_form": {
                "n": post.compare(ref["n_over_ne"], cf["n_over_n0"]),
                "T": post.compare(ref["T_over_Te"], cf["T_over_T0"]),
                "Ux": post.compare(ref["U_over_ae"][:, 0], cf["U_sqrt_beta0"][:, 0])}}
    # rotational symmetry: same (x, r) at two azimuths, not on a quadrature axis
    p1 = np.array([[0.3 * source.diameter_m, 0.4 * source.diameter_m, 0.0]])
    ang = np.radians(37.0)
    p2 = np.array([[p1[0, 0], 0.4 * source.diameter_m * np.cos(ang),
                    0.4 * source.diameter_m * np.sin(ang)]])
    r1 = collisionless.reference(p1, source)
    r2 = collisionless.reference(p2, source)
    out["rotational_symmetry_rel_diff"] = {
        "n": float(abs(r1["n"][0] / r2["n"][0] - 1)),
        "T": float(abs(r1["T"][0] / r2["T"][0] - 1))}
    write_json(res / "reference.json", out)
    return out


def overlap(case_dir: Path, reference: str | None = None, **_) -> dict:
    """Hybrid vs a full-DSMC case on the SAME mesh, cell by cell.

    Three regions, kept apart: the retained CFD region (CFD vs full DSMC: is
    the continuum solution valid where it was kept?), the hybrid's DSMC region
    (hybrid DSMC vs full DSMC: what the frozen interface did), and the axis.
    The sampling-noise floor is estimated as cv_step / sqrt(steps) -- a lower
    bound, since successive DSMC steps are correlated.
    """
    from . import breakdown as bd
    from . import post
    from .meshgeom import read_mesh
    case_dir = Path(case_dir)
    if not reference:
        raise WorkflowError("overlap needs --reference <full-DSMC case directory>")
    ref_dir = Path(reference)
    cfg = load_case(case_dir)
    source = config_mod.build_source(cfg)
    gas = source.gas
    res = results_dir(case_dir)
    cmesh = read_mesh(case_dir / "cfd")
    fmesh = read_mesh(ref_dir / "dsmc")
    if fmesh.n_cells != cmesh.n_cells or not np.allclose(
            fmesh.cell_centres, cmesh.cell_centres, rtol=0, atol=1e-9):
        raise WorkflowError("the reference is not on the same mesh as this case's CFD")
    hyb_post = json.loads((res / "dsmc_post.json").read_text())
    ref_post = json.loads((ref_dir / "results" / "dsmc_post.json").read_text())
    full = post.dsmc_fieldaverage(ref_dir / "dsmc" / ref_post["time"], gas)
    hyb = post.dsmc_fieldaverage(case_dir / "dsmc" / hyb_post["time"], gas)
    cm = np.load(res / "dsmc_to_cfd_cells.npy")
    retained = np.load(res / "retained.npy")
    st = bd.read_state(case_dir / "cfd", _latest(case_dir / "cfd"), gas)
    hyb_n = np.full(cmesh.n_cells, np.nan)
    hyb_t = np.full(cmesh.n_cells, np.nan)
    hyb_n[cm], hyb_t[cm] = hyb["n"], hyb["T"]
    occ_ref = full["parcels_per_step"]
    good = occ_ref >= float(cfg.get("dsmc.overlap_min_parcels_per_step", 1.0))
    v = cmesh.cell_volumes

    def stats(a, b, mask):
        m = mask & np.isfinite(a) & np.isfinite(b) & (b > 0)
        if not m.any():
            return {"cells": 0}
        rel = a[m] / b[m] - 1.0
        w = v[m] / v[m].sum()
        return {"cells": int(m.sum()), "volume_weighted_mean_rel": float(np.sum(w * rel)),
                "volume_weighted_mean_abs_rel": float(np.sum(w * np.abs(rel))),
                "p95_abs_rel": float(np.percentile(np.abs(rel), 95))}

    noise = {k: (d["sampling"].get("n_cv_per_step_median") or float("nan"))
             / np.sqrt(max(d["sampling"]["steps_averaged"], 1))
             for k, d in (("hybrid", hyb_post), ("reference", ref_post))}
    doc = {
        "label": cfg.get("label", cfg.profile), "reference_case": str(ref_dir),
        "reference_label": ref_post.get("label"),
        "interface_status": hyb_post.get("interface_status"),
        "retained_cfd_vs_full_dsmc": {
            "n": stats(st.rho / gas.mass_kg, full["n"], retained & good),
            "T": stats(st.T, full["T"], retained & good)},
        "hybrid_dsmc_vs_full_dsmc_in_dsmc_region": {
            "n": stats(hyb_n, full["n"], ~retained & good),
            "T": stats(hyb_t, full["T"], ~retained & good)},
        "sampling_noise_floor_rel_lower_bound": noise,
        "note": "a smoke comparison on an under-resolved mesh measures the machinery; "
                "its differences are not a physical accuracy statement",
    }
    h = float(np.min(cmesh.cell_size()[:, 1]))
    axis = post.axis_cells(cmesh.cell_centres, 0.51 * h)
    hyb_full_n = np.where(retained, st.rho / gas.mass_kg, hyb_n)
    hyb_full_t = np.where(retained, st.T, hyb_t)
    prof = post.line_profile(cmesh.cell_centres, v, axis,
                             {"n_hybrid": hyb_full_n, "n_full": full["n"],
                              "T_hybrid": hyb_full_t, "T_full": full["T"]})
    doc["centreline_hybrid_vs_full"] = {"n": post.compare(prof["n_hybrid"], prof["n_full"]),
                                        "T": post.compare(prof["T_hybrid"], prof["T_full"])}
    post.write_csv(res / "overlap_centreline.csv", {
        "x_over_D": prof["x"] / source.diameter_m,
        "n_hybrid_over_ne": prof["n_hybrid"] / source.number_density,
        "n_full_over_ne": prof["n_full"] / source.number_density,
        "T_hybrid_over_Te": prof["T_hybrid"] / source.temperature_K,
        "T_full_over_Te": prof["T_full"] / source.temperature_K})
    write_json(res / "overlap.json", doc)
    return doc


def hpc(case_dir: Path, scheduler: str = "pbs", np: int | None = None, **_) -> dict:
    """Write resumable batch scripts sized from the latest resource estimates."""
    from . import hpc as hpc_mod
    case_dir = Path(case_dir)
    cfg = load_case(case_dir)
    res = results_dir(case_dir)
    est = {}
    setup_doc = json.loads((res / "setup.json").read_text()) if (res / "setup.json").is_file() else {}
    est["cfd"] = (setup_doc.get("resources", {}).get("cfd", {}) or {}).get(
        "runtime_estimate_s_at_max_ranks")
    for name in ("dsmc_setup.json", "interface.json"):
        p = res / name
        if p.is_file():
            d = json.loads(p.read_text())
            e = d.get("run", {}).get("estimate") or d.get("dsmc_estimate")
            if e:
                est["dsmc"] = e.get("runtime_estimate_s_at_max_ranks")
                break
    else:
        est["dsmc"] = setup_doc.get("resources", {}).get("full_dsmc", {}).get(
            "runtime_estimate_s_at_max_ranks")
    doc = hpc_mod.write_scripts(case_dir, Path(__file__).resolve().parents[2], cfg,
                                scheduler, np, est)
    write_json(res / "hpc_scripts.json", doc)
    return doc


def report(case_dir: Path, **_) -> dict:
    """Collect every step's JSON into results/summary.json and draw the plots."""
    from . import report as rep
    case_dir = Path(case_dir)
    doc = rep.summarise(case_dir)
    rep.plots(case_dir)
    write_json(results_dir(case_dir) / "summary.json", doc)
    return {k: doc[k] for k in ("label", "stages_present", "headline") if k in doc}


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
