"""Shared fixtures for the OpenFOAM integration tests.

Every test here is marked ``needs_openfoam`` and is deselected by default::

    pytest -m "not needs_openfoam"    # the fast suite; no OpenFOAM required
    pytest -m needs_openfoam          # these

They run the real toolchain, so they are slow: one meshed case is built once per
session and shared, and the DSMC smoke test runs a deliberately short case.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BASE_CASE = REPO / "cases" / "markelov1999" / "baseCase"

#: Commands the mesh pipeline needs.
MESH_COMMANDS = ("blockMesh", "snappyHexMesh", "checkMesh")

#: Commands the solver needs.
SOLVER_COMMANDS = ("dsmcInitialise", "dsmcFoam")

#: The library that provides plumeFieldInflow.
PLUME_LIB = "libplumeDsmcBoundaryModels.so"


def have(command: str) -> bool:
    return shutil.which(command) is not None


def missing(commands) -> list:
    return [c for c in commands if not have(c)]


def plume_library_path() -> Path | None:
    """Where the custom inflow library is installed, or ``None``."""
    for variable in ("FOAM_USER_LIBBIN", "FOAM_SITE_LIBBIN", "FOAM_LIBBIN"):
        directory = os.environ.get(variable)
        if directory and (Path(directory) / PLUME_LIB).is_file():
            return Path(directory) / PLUME_LIB
    return None


def run(command, cwd: Path, timeout: int = 3600):
    """Run a command, returning the CompletedProcess with output captured."""
    return subprocess.run(
        command, cwd=str(cwd), capture_output=True, text=True,
        timeout=timeout, check=False)


def run_or_fail(command, cwd: Path, timeout: int = 3600):
    """Run a command and fail the test with its output if it does not succeed."""
    result = run(command, cwd, timeout)
    if result.returncode != 0:
        tail = "\n".join((result.stdout + result.stderr).splitlines()[-40:])
        pytest.fail(f"{' '.join(command)} failed with {result.returncode}:\n{tail}")
    return result


@pytest.fixture(scope="session")
def openfoam():
    """Skip the whole module unless the mesh toolchain is available."""
    absent = missing(MESH_COMMANDS)
    if absent:
        pytest.skip(f"OpenFOAM mesh utilities not on PATH: {absent}")
    return True


@pytest.fixture(scope="session")
def meshed_case(openfoam, tmp_path_factory):
    """A meshed p005psi case, built once and shared by every test that reads it.

    Built by invoking the case's own ``./Allmesh``, so the test exercises the
    script a user would run rather than a reimplementation of it.
    """
    work = tmp_path_factory.mktemp("markelov")
    case = work / "p005psi"
    shutil.copytree(BASE_CASE, case)

    for script in case.glob("All*"):
        script.chmod(0o755)

    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(REPO), env.get("PYTHONPATH", "")]).rstrip(os.pathsep)

    result = subprocess.run(
        ["bash", "./Allmesh"], cwd=str(case), capture_output=True, text=True,
        timeout=3600, env=env, check=False)
    if result.returncode != 0:
        tail = "\n".join((result.stdout + result.stderr).splitlines()[-40:])
        pytest.fail(f"./Allmesh failed:\n{tail}")

    return case


@pytest.fixture(scope="session")
def solved_case(meshed_case, tmp_path_factory):
    """A short DSMC run of the meshed case, for the smoke tests.

    Deliberately brief -- a few hundred steps rather than a converged run -- so
    the suite stays usable. That is enough to answer the questions these tests
    ask: does the custom model load and select the right patch, does it read a
    spatially varying density, and do particles leave through an open boundary?
    """
    absent = missing(SOLVER_COMMANDS)
    if absent:
        pytest.skip(f"dsmcFoam utilities not on PATH: {absent}")
    if plume_library_path() is None:
        pytest.skip(
            f"{PLUME_LIB} is not built; run applications/dsmcBoundaryModels/Allwmake")

    work = tmp_path_factory.mktemp("markelov-run")
    case = work / "p005psi"
    shutil.copytree(meshed_case, case)
    for script in case.glob("All*"):
        script.chmod(0o755)

    # Shorten the run. deltaT is untouched -- the Courant check depends on it --
    # so this only reduces the number of steps.
    for entry, value in (("endTime", "4e-05"), ("writeInterval", "2e-05"),
                         ("functions/fieldAverage1/timeStart", "0")):
        run_or_fail(["foamDictionary", "system/controlDict",
                     "-entry", entry, "-set", value], case)

    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(REPO), env.get("PYTHONPATH", "")]).rstrip(os.pathsep)

    result = subprocess.run(
        ["bash", "./Allrun", "--serial"], cwd=str(case), capture_output=True,
        text=True, timeout=7200, env=env, check=False)
    if result.returncode != 0:
        tail = "\n".join((result.stdout + result.stderr).splitlines()[-40:])
        pytest.fail(f"./Allrun failed:\n{tail}")

    return case
