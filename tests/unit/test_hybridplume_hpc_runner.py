"""Exercise the PBS-to-Allrun lifecycle using stubs, without OpenFOAM or PBS."""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.skipif(sys.platform == "win32" or shutil.which("bash") is None,
                              reason="run shell lifecycle checks under WSL/Linux")

# Python remains real for YAML/heredoc work; only CFD/DSMC and checkpoint CLI
# endpoints are stubbed. The real Allrun controls their ordering and refusal.
PYTHON_STUB = r'''
import json, os, pathlib, sys
args = sys.argv[1:]
if args == ["-c", "import plumetools.hybridplume"]:
    sys.exit(0)
def record(text):
    with open(os.environ["CALLS"], "a") as stream:
        stream.write(text+"\n")
if args[:2] == ["-m", "plumetools.hybridplume"]:
    stage, case = args[2], pathlib.Path(args[3])
    record("HP "+" ".join(args[2:]))
    res = case/"results"
    res.mkdir(exist_ok=True)
    if stage == "setup":
        import yaml
        (case/"hybrid.resolved.yaml").write_text(yaml.safe_dump({"config": {
            "mesh": {"geometry": "cartesian"}, "interface": {"mode": "fitted-ellipsoid"},
            "resources": {"default_ranks": 192, "max_ranks": 192}}}))
        for kind in ("cfd", "dsmc"):
            (case/kind/"constant/polyMesh").mkdir(parents=True, exist_ok=True)
            (case/kind/"constant/polyMesh/owner").touch()
            (case/kind/"system").mkdir(exist_ok=True)
        (case/"dsmc/.mesh-checked").touch()
    if stage == "interface":
        (res/"interface.json").write_text(json.dumps({
            "allow_unsuitable": "--allow-unsuitable" in args}))
    if stage == "dsmc-fields" and os.environ.get("FAIL_ASSESSMENT") == "1":
        allowed = json.loads((res/"interface.json").read_text())["allow_unsuitable"]
        if not allowed:
            sys.exit(1)
        record("SMOKE_ONLY true")
    sys.exit(0)
if args[:2] == ["-m", "plumetools.foamio.checkpoints"]:
    if args[2] == "finished":
        sys.exit(1)
    print("DECOMPOSE=all\nPROC_TIME=0\nRECONSTRUCT_FIRST=0")
    sys.exit(0)
os.execv(os.environ["REAL_PYTHON"], [os.environ["REAL_PYTHON"], *args])
'''


@pytest.mark.parametrize("from_root,unsuitable,allow", [
    (True, False, False), (False, False, False),
    (True, True, False), (True, True, True)])
def test_pbs_launch_runs_full_chain_and_preserves_coupling_gate(tmp_path, from_root,
                                                              unsuitable, allow):
    repo = tmp_path/"checkout"
    scripts = repo/"cases/cai2012-hybrid"
    scripts.mkdir(parents=True)
    for name in ("Allrun", "run3DCase.sh", "hybrid.yaml"):
        shutil.copy2(ROOT/"cases/cai2012-hybrid"/name, scripts/name)
    bin_dir = tmp_path/"bin"
    lib = tmp_path/"lib"
    bin_dir.mkdir()
    lib.mkdir()
    for name in ("libhybridPlumeThermo.so", "libplumeDsmcBoundaryModels.so"):
        (lib/name).touch()
    stub = tmp_path/"python-stub.py"
    stub.write_text(PYTHON_STUB)
    for name in ("python", "python3"):
        path = bin_dir/name
        path.write_text(f'#!/bin/sh\nexec "{sys.executable}" "{stub}" "$@"\n')
        path.chmod(0o755)
    for name in ("foamDictionary", "mpirun", "decomposePar", "reconstructPar", "dsmcInitialise"):
        path = bin_dir/name
        path.write_text(f'#!/bin/sh\nprintf "%s\\n" "{name} $*" >> "$CALLS"\n')
        path.chmod(0o755)
    env_file = tmp_path/"env.sh"
    env_file.write_text(f'export FOAM_APPBIN="{bin_dir}"\nexport FOAM_USER_LIBBIN="{lib}"\n')
    calls = tmp_path/"calls"
    env = dict(os.environ, PBS_O_WORKDIR=str(repo if from_root else scripts),
               HYBRID_HPC_ENV=str(env_file), CALLS=str(calls), REAL_PYTHON=sys.executable,
               FAIL_ASSESSMENT=str(int(unsuitable)), ALLOW_UNSUITABLE=str(int(allow)),
               PATH=f"{bin_dir}:{os.environ.get('PATH', '')}")
    env.pop("CASE_DIR", None)
    env.pop("PROFILE", None)
    result = subprocess.run(["bash", str(scripts/"run3DCase.sh")], cwd=repo,
                            env=env, text=True, capture_output=True, timeout=30)
    events = calls.read_text().splitlines()
    if unsuitable and not allow:
        assert result.returncode != 0
        assert not any("dsmcInitialise" in event or "dsmcFoam" in event for event in events)
        assert not any(event.startswith("HP dsmc-post") for event in events)
    else:
        assert result.returncode == 0, result.stderr
        assert events.index("dsmcInitialise ") < events.index("mpirun -np 192 dsmcFoam -parallel")
        assert next(i for i, event in enumerate(events) if event.startswith("HP dsmc-post")) > events.index("mpirun -np 192 dsmcFoam -parallel")
        assert any(event.startswith("HP report") for event in events)
    setup = next(event for event in events if event.startswith("HP setup"))
    assert "--profile hpc-kn1em3-ellipsoid" in setup and "--np 192" in setup
    if allow:
        assert "SMOKE_ONLY true" in events
