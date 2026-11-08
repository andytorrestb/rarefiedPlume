"""Immutable measured target plans, portable artifacts, and source verification."""
from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import tempfile
import time
from typing import Optional

import numpy as np

from .assess import usable_window
from .geometry import file_sha256, mesh_fingerprint
from .model import Geometry, RefinementPlan

SCHEMA_VERSION = 1
ACTIONS = {"criteria_satisfied", "needs_mesh", "needs_parcels", "unresolved"}


def _publish_stage(stage, destination):
    """Atomically publish; tolerate bounded transient Windows sharing locks."""
    for attempt in range(8):
        if destination.exists():
            raise FileExistsError(f'destination appeared during staging: {destination}')
        try:
            stage.rename(destination)
            return
        except PermissionError as exc:
            if getattr(exc,'winerror',None) not in (5,32) or attempt==7:
                raise
            time.sleep(min(.025*2**attempt,.4))


def json_safe(value):
    """Use JSON null for unavailable diagnostics, never non-standard NaN tokens."""
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def plan_refinement(geometry, measured, assessment, *, particle_weight: Optional[float] = None) -> RefinementPlan:
    """Select only trusted failures, with no executable mesh guess for data holes."""
    if (assessment.summary.get("mesh_sha256") != geometry.mesh_sha256
            or assessment.summary.get("source_sha256") != measured.source_sha256):
        raise ValueError("assessment belongs to different geometry or measurements")
    action = assessment.summary["action"]
    selected = np.flatnonzero(assessment.mesh_status == 0) if action == "needs_mesh" else np.empty(0, dtype=np.int64)
    weight = None
    p = measured.mean_parcels
    if particle_weight is not None:
        if not math.isfinite(float(particle_weight)) or particle_weight <= 0:
            raise ValueError("particle_weight must be finite and positive")
        usable = (usable_window(measured.window) and np.all(np.isfinite(p) & (p > 0))
                  and np.all(np.isfinite(measured.number_density) & (measured.number_density > 0))
                  and np.all(np.isfinite(measured.temperature_K) & (measured.temperature_K > 0)))
        if usable and action != "unresolved":
            weight = float(particle_weight) * min(1.0, float(np.min(p / assessment.targets.min_mean_parcels)))
    provenance = json_safe({**measured.provenance, "window": measured.window,
                            "assessment": assessment.summary, "targets": asdict(assessment.targets),
                            "density_estimate_kind": "parent-measured", "requires_flow_rerun": True})
    return RefinementPlan(action, geometry.mesh_sha256, measured.source_sha256, geometry,
                          assessment.target_edge_m.copy(), measured.number_density.copy(), selected,
                          assessment.targets.min_mean_parcels, weight, provenance)


def verify_sources(case: Path, plan: RefinementPlan) -> None:
    """Refuse a stale plan even if the new mesh/fields have the same cell count."""
    case = Path(case)
    if mesh_fingerprint(case) != plan.mesh_sha256:
        raise ValueError("mesh changed since planning; reassess the case")
    for relative, expected in plan.provenance.get("source_files", {}).items():
        path = (case / relative).resolve()
        try:
            path.relative_to(case.resolve())
        except ValueError:
            raise ValueError(f"invalid source reference {relative!r}") from None
        if not path.is_file() or file_sha256(path) != expected:
            raise ValueError(f"source {relative} changed or is missing; the plan is stale")


def _arrays(plan):
    g = plan.source_geometry
    return {"centres_m": g.centres_m, "volumes_m3": g.volumes_m3, "max_edge_m": g.max_edge_m,
            "target_edge_m": plan.target_edge_m, "density_estimate": plan.density_estimate,
            "selected_cells": plan.selected_cells}


def _metadata(plan):
    return {"kind": "plumetools.refinement.plan", "schema_version": SCHEMA_VERSION,
            "action": plan.action, "mesh_sha256": plan.mesh_sha256, "source_sha256": plan.source_sha256,
            "target_parcels": plan.target_parcels, "recommended_weight": plan.recommended_weight,
            "provenance": json_safe(plan.provenance)}


def plan_identity(plan):
    """Hash semantic metadata and array content rather than archive timestamps."""
    digest = hashlib.sha256(json.dumps(_metadata(plan), sort_keys=True, allow_nan=False).encode())
    for name, array in sorted(_arrays(plan).items()):
        a = np.asarray(array, dtype="<i8" if name == "selected_cells" else "<f8")
        digest.update(name.encode() + str(a.shape).encode() + a.tobytes())
    return digest.hexdigest()


def _validate(plan):
    n = len(plan.source_geometry.volumes_m3)
    if plan.action not in ACTIONS or plan.mesh_sha256 != plan.source_geometry.mesh_sha256:
        raise ValueError("invalid plan action or mesh identity")
    for name in ("target_edge_m", "density_estimate"):
        if np.asarray(getattr(plan, name)).shape != (n,):
            raise ValueError(f"plan {name} does not match its cell geometry")
    selected = np.asarray(plan.selected_cells)
    if selected.ndim != 1 or not np.issubdtype(selected.dtype, np.integer):
        raise ValueError("selected_cells must be an integer vector")
    if np.any(selected < 0) or np.any(selected >= n) or len(np.unique(selected)) != selected.size:
        raise ValueError("invalid selected cell identities")
    if plan.action != "needs_mesh" and selected.size:
        raise ValueError("non-mesh plan must not contain executable refinement cells")
    if plan.action == "needs_mesh" and not np.all(np.isfinite(plan.target_edge_m) & (plan.target_edge_m > 0)):
        raise ValueError("mesh plan contains unavailable target edges")
    if not math.isfinite(float(plan.target_parcels)) or plan.target_parcels <= 0:
        raise ValueError("plan target_parcels must be finite and positive")
    if plan.recommended_weight is not None and (not math.isfinite(plan.recommended_weight) or plan.recommended_weight <= 0):
        raise ValueError("invalid recommended particle weight")


def write_plan(plan, directory: Path) -> Path:
    """Write a fresh artifact atomically; preserve an identical existing artifact."""
    from .report import write_cells_csv
    _validate(plan)
    directory = Path(directory)
    identity = plan_identity(plan)
    if directory.exists() and any(directory.iterdir()):
        old = directory / "plan.json"
        if old.is_file() and plan_identity(read_plan(old)) == identity:
            return old
        raise ValueError(f"output already exists: {directory}; choose a fresh directory")
    directory.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".refinement-plan-", dir=directory.parent) as temporary:
        stage = Path(temporary) / "artifact"; stage.mkdir()
        np.savez_compressed(stage / "targets.npz", **_arrays(plan))
        document = {**_metadata(plan), "plan_id": identity, "arrays": "targets.npz",
                    "arrays_sha256": file_sha256(stage / "targets.npz")}
        (stage / "plan.json").write_text(json.dumps(document, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        write_cells_csv(stage / "cells.csv", plan.source_geometry, target=plan.target_edge_m,
                        density=plan.density_estimate, selected=plan.selected_cells)
        if directory.exists():
            directory.rmdir()  # Only the explicitly supplied empty destination.
        _publish_stage(stage,directory)
    return directory / "plan.json"


def read_plan(path: Path) -> RefinementPlan:
    """Load a versioned plan and verify archive integrity before using arrays."""
    path = Path(path)
    if path.is_dir():
        path = path / "plan.json"
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
        if document.get("kind") != "plumetools.refinement.plan" or document.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("unsupported refinement plan schema")
        if document.get("arrays") != "targets.npz":
            raise ValueError("invalid plan array reference")
        archive = path.parent / "targets.npz"
        if file_sha256(archive) != document["arrays_sha256"]:
            raise ValueError("plan arrays checksum does not match; artifact is corrupt")
        with np.load(archive, allow_pickle=False) as a:
            g = Geometry(a["centres_m"], a["volumes_m3"], a["max_edge_m"], document["mesh_sha256"])
            plan = RefinementPlan(document["action"], document["mesh_sha256"], document["source_sha256"], g,
                                  a["target_edge_m"], a["density_estimate"], a["selected_cells"],
                                  document["target_parcels"], document["recommended_weight"], document["provenance"])
        _validate(plan)
        if plan_identity(plan) != document["plan_id"]:
            raise ValueError("plan content hash does not match; artifact is corrupt")
        return plan
    except (OSError, KeyError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read refinement plan {path}: {exc}") from None
