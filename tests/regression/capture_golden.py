"""Capture golden inflow data from the UNMODIFIED legacy processInflowData.py.

This runs once, at commit `pre-refactor-baseline`, to freeze the exact numbers the
pre-refactor code produced for cases/3d-inflow. Everything the refactor does is
then checked against this snapshot.

The legacy script is imported without editing it. Two things are worked around:

  * it calls `plumeSourceFlowModel()` at module scope (line 556), so the trailing
    call is stripped before exec and invoked explicitly afterwards;
  * `readMeshStats()` shells out to `checkMesh`, which needs OpenFOAM. It is
    replaced by a stub whose `points`/`faces` values are READ FROM THE MESH, not
    hard-coded -- those are the only two keys the script actually consumes.

Face insertion order is preserved: printInflow.py writes values in dict-iteration
order, and dsmcFoam+ expects that order for the patch.

Usage (from the repository root):
    python tests/regression/capture_golden.py
"""

from __future__ import annotations

import contextlib
import io
import os
import shutil
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "cases" / "3d-inflow"
GOLDEN = Path(__file__).resolve().parent / "golden"

FIELD_FILES = ("boundaryU", "boundaryT", "boundaryNumberDensity_Ar")


def count_from_header(path: Path) -> int:
    """Read an OpenFOAM list count: the integer on the line before the opening '('."""
    lines = path.read_text().splitlines()
    open_idx = next(i for i, line in enumerate(lines) if line.strip() == "(")
    return int(lines[open_idx - 1].strip())


def build_mesh_stats() -> dict:
    """Stand in for readMeshStats() without invoking checkMesh.

    Only 'points' and 'faces' are read by the legacy script; both are taken from
    the mesh files themselves so this cannot drift from the committed mesh.
    """
    n_points = count_from_header(CASE / "constant/polyMesh/points")
    n_faces = count_from_header(CASE / "constant/polyMesh/faces")
    n_internal = count_from_header(CASE / "constant/polyMesh/neighbour")
    print(f"  mesh: {n_points} points, {n_faces} faces, {n_internal} internal")
    return {
        "points": n_points,
        "faces": n_faces,
        "internal_faces": n_internal,
        "cells": 0,
        "faces_per_cell": 0,
        "boundary_patches": 0,
        "point_zones": 0,
        "face_zones": 0,
        "cell_zones": 0,
    }


def load_legacy_module():
    """Exec processInflowData.py with its module-scope driver call removed."""
    # Retired to *.legacy in Step 5; still the authority for this snapshot.
    source_path = CASE / "processInflowData.py.legacy"
    if not source_path.is_file():
        source_path = CASE / "processInflowData.py"
    source = source_path.read_text()

    marker = "plumeSourceFlowModel()"
    idx = source.rindex(marker)
    if source[idx + len(marker):].strip():
        raise RuntimeError("unexpected code after the driver call; refusing to guess")
    source = source[:idx]

    namespace = {"__name__": "legacy_processInflowData", "__file__": str(source_path)}
    exec(compile(source, str(source_path), "exec"), namespace)
    return namespace


def main() -> int:
    GOLDEN.mkdir(parents=True, exist_ok=True)
    (CASE / "0").mkdir(exist_ok=True)

    print(f"capturing golden from {CASE}")
    original_cwd = Path.cwd()
    # The legacy script uses relative paths for both the mesh and its 0/ output.
    os.chdir(CASE)
    sys.path.insert(0, str(CASE))
    try:
        namespace = load_legacy_module()
        namespace["readMeshStats"] = build_mesh_stats

        captured = {}
        real_print_inflow = namespace["printInflowSurface"]

        def capturing_print_inflow(inflow):
            captured["rhoN"], captured["U"], captured["T"] = inflow
            return real_print_inflow(inflow)

        namespace["printInflowSurface"] = capturing_print_inflow

        # The legacy code prints several lines per face; discard that noise but keep
        # stderr so genuine failures still surface.
        print("  running legacy plumeSourceFlowModel() (quiet, ~1 min)...")
        with contextlib.redirect_stdout(io.StringIO()):
            namespace["plumeSourceFlowModel"]()
    finally:
        os.chdir(original_cwd)
        sys.path.remove(str(CASE))

    face_rhoN, face_U, face_T = captured["rhoN"], captured["U"], captured["T"]

    # Insertion order is the contract with dsmcFoam+; store it explicitly.
    labels = list(face_rhoN.keys())
    if list(face_U.keys()) != labels or list(face_T.keys()) != labels:
        raise RuntimeError("rhoN/U/T dictionaries disagree on face order")

    rhoN = np.array([float(face_rhoN[f]) for f in labels], dtype=np.float64)
    U = np.array([[float(c) for c in face_U[f]] for f in labels], dtype=np.float64)
    T = np.array([float(face_T[f]) for f in labels], dtype=np.float64)

    out = GOLDEN / "3d_inflow_v0.npz"
    np.savez_compressed(out, labels=np.array(labels), rhoN=rhoN, U=U, T=T)

    for name in FIELD_FILES:
        shutil.copyfile(CASE / "0" / name, GOLDEN / f"{name}.txt")

    print(f"\n  faces          : {len(labels)}")
    print(f"  rhoN [min,max] : {rhoN.min():.6e} {rhoN.max():.6e}")
    print(f"  |U|  [min,max] : {np.linalg.norm(U, axis=1).min():.6f} "
          f"{np.linalg.norm(U, axis=1).max():.6f}")
    print(f"  T    [min,max] : {T.min()} {T.max()}")
    print(f"  finite         : rhoN={np.isfinite(rhoN).all()} U={np.isfinite(U).all()}")
    print(f"\n  wrote {out}")
    for name in FIELD_FILES:
        print(f"  wrote {GOLDEN / (name + '.txt')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
