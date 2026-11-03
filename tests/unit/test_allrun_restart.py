"""The restart decisions of a case's ``Allrun``, with stub OpenFOAM tools.

``Allrun`` takes what to resume from, and whether to decompose, from
``python -m plumetools.foamio.checkpoints plan``. The stubs below only log
their arguments, so these tests check which tools ``Allrun`` would run, in
which order, on which saved state -- not the tools themselves (the OpenFOAM
tier does that).

Defect 1, as it was: ``Allrun`` read the newest time of ``processor0`` alone.
A job killed while writing left processor0 one time ahead of another rank, and
the next run resumed from that time -- which aborts in parallel, or, with the
directory present but its parcel cloud missing, resumes that subdomain empty.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
ALLRUN = REPO / "cases" / "cai2012" / "baseCase" / "Allrun"
FOOTER = "\n// " + "*" * 73 + " //\n"

pytestmark = pytest.mark.skipif(
    sys.platform == "win32" or shutil.which("bash") is None,
    reason="bash on native Windows is WSL's and cannot open a Windows path; "
           "the suite runs under WSL")


def checkpoint(base: Path, name: str) -> None:
    (base / name / "uniform").mkdir(parents=True)
    (base / name / "uniform" / "time").write_text(f"value {name};" + FOOTER)
    (base / name / "rhoNMean").write_text("internalField uniform 0;" + FOOTER)
    (base / name / "lagrangian" / "dsmc").mkdir(parents=True)
    (base / name / "lagrangian" / "dsmc" / "positions").write_text("0()" + FOOTER)


@pytest.fixture
def case(tmp_path):
    """A case directory with a mesh, two ranks' worth of decomposition, and
    stub tools that record their calls in ``calls``."""
    case = tmp_path / "case"
    (case / "constant" / "polyMesh").mkdir(parents=True)
    (case / "system").mkdir()
    shutil.copy2(ALLRUN, case / "Allrun")
    (case / "case.yaml").write_text("{}\n")
    (case / "runInflow.py").write_text("raise SystemExit('not on a resume')\n")
    (case / "constant" / "dsmcProperties").write_text("InflowBoundaryModel FreeStream;\n")
    (case / "system" / "controlDict").write_text(
        "startFrom latestTime;\nendTime 1;\ndeltaT 0.01;\n")
    (case / "system" / "decomposeParDict").write_text(
        "numberOfSubdomains 2;\nmethod simple;\ncoeffs\n{\n    n (2 1 1);\n}\n")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    calls = tmp_path / "calls"
    for tool in ("dsmcFoam", "dsmcInitialise", "decomposePar", "reconstructPar",
                 "mpirun", "foamDictionary"):
        stub = bin_dir / tool
        stub.write_text(f'#!/bin/sh\necho "{tool} $*" >> "{calls.as_posix()}"\n')
        stub.chmod(0o755)
    python = bin_dir / "python"
    python.write_text(f'#!/bin/sh\nexec "{Path(sys.executable).as_posix()}" "$@"\n')
    python.chmod(0o755)
    env = dict(os.environ)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [str(REPO), env.get("PYTHONPATH")]))
    for key in ("FOAM_APPBIN", "FOAM_USER_LIBBIN", "FOAM_SITE_LIBBIN", "FOAM_LIBBIN"):
        env.pop(key, None)

    def run(*args):
        result = subprocess.run(["bash", "./Allrun", *args], cwd=case,
                                capture_output=True, text=True, env=env, timeout=300)
        logged = calls.read_text().splitlines() if calls.is_file() else []
        return result, logged

    return case, run


def test_a_rank_behind_processor0_stops_the_run_before_any_tool(case):
    directory, run = case
    for rank in (0, 1):
        checkpoint(directory / f"processor{rank}", "0.1")
    checkpoint(directory / "processor0", "0.2")        # processor1 was killed
    result, calls = run()
    assert result.returncode != 0
    assert "processor1 has no 0.2/" in result.stderr
    assert "cannot be continued as it is" in result.stderr
    assert calls == []
    assert (directory / "processor0" / "0.2").is_dir()  # nothing was changed


def test_complete_ranks_are_resumed_in_place(case):
    directory, run = case
    for rank in (0, 1):
        checkpoint(directory / f"processor{rank}", "0.1")
        checkpoint(directory / f"processor{rank}", "0.2")
    result, calls = run()
    assert result.returncode == 0, result.stderr
    assert "already hold t = 0.2; not re-decomposing" in result.stdout
    assert [c.split()[0] for c in calls] == ["mpirun", "reconstructPar"]


def test_a_new_decomposition_reconstructs_the_newer_state_first(case):
    """decomposePar -force deletes processor*/: their newer state is kept."""
    directory, run = case
    checkpoint(directory, "0.1")                       # reconstructed earlier
    for rank in (0, 1):
        checkpoint(directory / f"processor{rank}", "0.1")
        checkpoint(directory / f"processor{rank}", "0.2")
    result, calls = run("-np", "4")
    assert result.returncode == 0, result.stderr
    tools = [c.split()[0] for c in calls if not c.startswith("foamDictionary")]
    assert tools == ["reconstructPar", "decomposePar", "mpirun", "reconstructPar"]
    assert "-newTimes" in calls[[c.split()[0] for c in calls].index("reconstructPar")]
    assert any("decomposePar -force -latestTime" in c for c in calls)
    assert "reconstructing t = 0.2 from the old decomposition first" in result.stdout
