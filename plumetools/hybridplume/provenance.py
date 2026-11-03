"""What produced a result: repository SHAs, source hashes, solver binaries, toolchain.

Collected at every step that writes a report, so a results directory can be
traced back without trusting anyone's memory. Every probe is allowed to fail:
a missing tool is recorded as ``null`` with the reason, never invented.
"""

from __future__ import annotations

import hashlib
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

#: Source files whose content defines what this family computes.
HASHED_SOURCES = (
    "plumetools/hybridplume",
    "applications/dsmcBoundaryModels/plumeFieldInflow/plumeFieldInflow.C",
    "applications/thermophysicalModels/hybridPlumeThermo/hybridPlumeThermos.C",
    "cases/cai2012-hybrid/hybrid.yaml",
    "cases/cai2012-hybrid/Allrun",
)


def _run(cmd, cwd=None, timeout=20) -> str | None:
    try:
        out = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                             timeout=timeout, check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    text = (out.stdout or "") + (out.stderr or "")
    return text.strip() or None


def sha256(path: Path) -> str | None:
    path = Path(path)
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_state(repo: Path) -> dict:
    repo = Path(repo)
    if not (repo / ".git").exists():
        return {"path": str(repo), "available": False}
    sha = _run(["git", "rev-parse", "HEAD"], cwd=repo)
    branch = _run(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=repo)
    status = _run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=repo)
    return {"path": str(repo), "available": sha is not None, "sha": sha,
            "branch": branch, "dirty_tracked_files": bool(status)}


def source_hashes(repo: Path) -> dict:
    out = {}
    for rel in HASHED_SOURCES:
        path = Path(repo) / rel
        files = sorted(path.rglob("*.py")) if path.is_dir() else [path]
        for f in files:
            if "__pycache__" in f.parts:
                continue
            out[str(f.relative_to(repo)).replace("\\", "/")] = sha256(f)
    return out


def toolchain() -> dict:
    """OpenFOAM, compiler and MPI as the current environment sees them."""
    def which(name):
        return shutil.which(name)
    binaries = {b: which(b) for b in ("rhoCentralFoam", "dsmcFoam", "dsmcInitialise",
                                       "subsetMesh", "topoSet", "blockMesh",
                                       "decomposePar", "mpirun")}
    libs = {}
    for name in ("libplumeDsmcBoundaryModels.so", "libhybridPlumeThermo.so"):
        for var in ("FOAM_USER_LIBBIN", "FOAM_SITE_LIBBIN", "FOAM_LIBBIN"):
            d = os.environ.get(var)
            if d and (Path(d) / name).is_file():
                libs[name] = {"path": str(Path(d) / name), "sha256": sha256(Path(d) / name)}
                break
        else:
            libs[name] = None
    return {
        "openfoam": {k: os.environ.get(k) for k in (
            "WM_PROJECT", "WM_PROJECT_VERSION", "WM_OPTIONS", "WM_COMPILER",
            "WM_MPLIB", "FOAM_MPI", "WM_PROJECT_DIR")},
        "binaries": {k: ({"path": v, "sha256": sha256(Path(v))} if v else None)
                     for k, v in binaries.items()},
        "libraries": libs,
        "compiler": (_run(["g++", "--version"]) or "").splitlines()[:1] or None,
        "mpi": (_run(["mpirun", "--version"]) or "").splitlines()[:1] or None,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "hostname": platform.node(),
    }


def hystrath_state() -> dict:
    """The hyStrath checkout, if one is found. This family does not use it."""
    candidates = [os.environ.get("HYSTRATH_DIR")]
    candidates += [str(Path.home() / "Dev" / "hyStrath"), "/mnt/c/Users/andy9/Dev/hyStrath"]
    for c in candidates:
        if c and Path(c, ".git").exists():
            state = git_state(Path(c))
            state["used_by_this_family"] = False
            return state
    return {"available": False, "used_by_this_family": False}


def collect(repo: Path, *, with_toolchain: bool = True) -> dict:
    import numpy
    import scipy
    doc = {
        "rarefiedPlume": git_state(repo),
        "hyStrath": hystrath_state(),
        "source_sha256": source_hashes(repo),
        "python_packages": {"numpy": numpy.__version__, "scipy": scipy.__version__},
    }
    if with_toolchain:
        doc["toolchain"] = toolchain()
    return doc
