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
