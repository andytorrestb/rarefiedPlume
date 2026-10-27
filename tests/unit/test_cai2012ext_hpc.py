"""The extension's PBS job scripts, and what runCase.sh calls around ./Allrun.

One job per case: ``baseCase/run*.sh`` are copied into every case and submitted
from it; ``runMeshAll.sh`` and ``runAnalyse.sh`` run from the family directory;
all of them source the one ``hpc_env.sh``. The scripts are read as text, parsed
with ``bash -n``, and runCase.sh is run end to end against a stand-in
environment and a stand-in ``Allrun``. :mod:`plumetools.cai2012ext.hpc` is run
on real ``case.yaml`` files.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest
import yaml

from plumetools.cai2012.config import load_case_config
from plumetools.cai2012ext import feasibility, hpc, policy

REPO = Path(__file__).resolve().parents[2]
FAMILY = REPO / "cases" / "cai2012-kn-extension"

CASE_SCRIPTS = [FAMILY / "baseCase" / name for name in
                ("runMesh.sh", "runCase.sh", "runPost.sh", "runViz.sh")]
FAMILY_SCRIPTS = [FAMILY / name for name in ("runMeshAll.sh", "runAnalyse.sh")]
JOB_SCRIPTS = CASE_SCRIPTS + FAMILY_SCRIPTS

def bash_executable():
    """Prefer a native shell on Windows, rather than the WSL launch shim."""
    if sys.platform == "win32":
        git_bash = Path(os.environ.get("ProgramFiles", "C:/Program Files")) / "Git" / "bin" / "bash.exe"
        if git_bash.is_file():
            return str(git_bash)
    return shutil.which("bash")


BASH = bash_executable()
needs_bash = pytest.mark.skipif(BASH is None, reason="no bash")
needs_flock = pytest.mark.skipif(sys.platform == "win32", reason="no fcntl")


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def code_of(path: Path) -> str:
    """The script without its comment lines."""
    return "\n".join(line for line in read(path).splitlines()
                     if not line.lstrip().startswith("#"))


def write_case(directory: Path, knudsen: float) -> Path:
    """A case directory holding the template's case.yaml at ``knudsen``."""
    directory.mkdir(parents=True)
    data = yaml.safe_load(read(FAMILY / "baseCase" / "case.yaml"))
    data["exit"]["knudsen"] = knudsen
    (directory / "case.yaml").write_text(yaml.safe_dump(data, sort_keys=False),
                                         encoding="utf-8")
    return directory


def write_control(case: Path, end: float, delta: float) -> None:
    (case / "system").mkdir(exist_ok=True)
    (case / "system" / "controlDict").write_text(
        "FoamFile\n{\n    object      controlDict;\n}\n"
        "stopAt          endTime;\n"
        f"endTime         {end:g};\n"
        f"deltaT          {delta:g};\n", encoding="utf-8")


def python_env(tmp_path: Path) -> dict:
    """PATH with a ``python`` that is this interpreter, and plumetools on it."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    wrapper = bin_dir / "python"
    wrapper.write_text(f'#!/bin/sh\nexec "{Path(sys.executable).as_posix()}" "$@"\n',
                       encoding="utf-8")
    wrapper.chmod(0o755)
    env = {key: value for key, value in os.environ.items()
           if not key.startswith("PBS_")}
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    env["PYTHONPATH"] = os.pathsep.join(
        filter(None, [str(REPO), env.get("PYTHONPATH")]))
    return env


# --------------------------------------------------------------------------- #
# the job scripts, as text
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("script", JOB_SCRIPTS, ids=lambda p: p.name)
def test_every_job_starts_in_its_submission_directory(script):
    body = read(script)
    assert 'cd "$PBS_O_WORKDIR"' in body
    assert '[ -z "${PBS_O_WORKDIR:-}" ]' in body
    # hpc_env.sh: ./, ../ or ../../ from there, and nowhere else.
    assert "for dir in . .. ../..; do" in body
    assert '[ -f "$dir/hpc_env.sh" ]' in body
    assert 'source "$HPC_ENV"' in body
    assert "set -eo pipefail" in body and "set -u" not in code_of(script)
    assert re.search(r"^trap .* ERR$", body, re.MULTILINE)


@pytest.mark.parametrize("script", JOB_SCRIPTS, ids=lambda p: p.name)
def test_no_job_names_allviz_a_rank_count_or_allruncases(script):
    body = read(script)
    assert "Allviz" not in body
    assert not re.search(r"(^|\s)-np\b", body)
    assert "AllrunCases" not in body


@pytest.mark.parametrize("script", JOB_SCRIPTS, ids=lambda p: p.name)
def test_the_environment_lives_in_hpc_env_only(script):
    code = code_of(script)
    for setting in ("module ", "spack", "etc/bashrc", "activate", "PVPYTHON_PATH",
                    "VTK_DEFAULT_OPENGL_WINDOW"):
        assert setting not in code, setting


@pytest.mark.parametrize("script,queue,walltime,nodes,name", [
    ("baseCase/runMesh.sh", "shortq", "2:00:00", "nodes=1:ppn=1", "knx-mesh"),
    ("baseCase/runCase.sh", "longq", "120:00:00", "nodes=1:ppn=192", "knx-run"),
    ("baseCase/runPost.sh", "longq", "120:00:00", "nodes=1:ppn=192", "knx-post"),
    ("baseCase/runViz.sh", "longq", "120:00:00", "nodes=1:ppn=192", "knx-viz"),
    ("runMeshAll.sh", "shortq", "4:00:00", "nodes=1:ppn=1", "knx-meshall"),
    ("runAnalyse.sh", "shortq", "1:00:00", "nodes=1:ppn=1", "knx-analyse"),
])
def test_the_requested_resources(script, queue, walltime, nodes, name):
    directives = [line for line in read(FAMILY / script).splitlines()
                  if line.startswith("#PBS")]
    assert directives[0] == "#PBS -S /bin/bash"
    assert f"#PBS -q {queue}" in directives
    assert f"#PBS -l walltime={walltime}" in directives
    assert f"#PBS -l {nodes}" in directives
    assert f"#PBS -N {name}" in directives


def test_the_post_job_runs_the_case_chain_without_imagery_or_analysis():
    assert 'bash ../../AllpostCases --case "$CASE" --no-viz --no-analyse' in \
        code_of(FAMILY / "baseCase" / "runPost.sh")


def test_the_viz_job_draws_this_cases_images_only():
    assert 'bash ../../AllpostCases --case "$CASE" --viz-only' in \
        code_of(FAMILY / "baseCase" / "runViz.sh")


def test_the_family_jobs_generate_validate_mesh_and_analyse_in_order():
    code = code_of(FAMILY / "runMeshAll.sh")
    steps = ["python generate_cases.py\n",
             "python validate_cases.py --inputs-only\n",
             "bash ./AllmeshCases\n",
             "python validate_cases.py\n"]
    positions = [code.index(step) for step in steps]
    assert positions == sorted(positions)
    assert "bash ./AllpostCases --analyse-only" in code_of(FAMILY / "runAnalyse.sh")


def test_the_run_job_guards_labels_runs_and_records_in_that_order():
    code = code_of(FAMILY / "baseCase" / "runCase.sh")
    steps = ["hpc guard .", "hpc label .", "hpc finished .", "bash ./Allrun ",
             "hpc record ."]
    positions = [code.index(step) for step in steps]
    assert positions == sorted(positions)
    assert "bash ./Allrun >> output.log 2> error.log" in code


def test_the_environment_file_uses_paraview_6p2_and_shields_the_bashrc():
    env = read(FAMILY / "hpc_env.sh")
    assert ('export PVPYTHON_PATH="${PVPYTHON_PATH:-/home2/torresfa/src/paraview/'
            'ParaView-6.2.0-MPI-Linux-Python3.12-x86_64/bin/pvpython}"') in env
    assert "VTK_DEFAULT_OPENGL_WINDOW=vtkOSOpenGLRenderWindow" in env
    assert "ParaView-5" not in env
    assert "if source /home2/torresfa/scratch/OpenFOAM/OpenFOAM-v2512/etc/bashrc; then" in env
    assert env.index("set --\n") < env.index("etc/bashrc; then")


def test_the_template_runs_on_192_ranks():
    assert load_case_config(FAMILY / "baseCase").dsmc.n_subdomains == 192


# --------------------------------------------------------------------------- #
# the shell, parsing and running
# --------------------------------------------------------------------------- #

def _shell_scripts():
    found = []
    for directory in (FAMILY, FAMILY / "baseCase"):
        for path in sorted(directory.iterdir()):
            if not path.is_file():
                continue
            with open(path, "rb") as handle:
                first = handle.readline()
            if path.suffix == ".sh" or first.startswith(b"#!/bin/bash"):
                found.append(path)
    return found


SHELL_SCRIPTS = _shell_scripts()


def test_the_shell_script_list_is_complete():
    names = {path.name for path in SHELL_SCRIPTS}
    assert {"hpc_env.sh", "runMeshAll.sh", "runAnalyse.sh", "runCase.sh",
            "AllmeshCases", "AllpostCases", "Allrun"} <= names
    assert "AllrunCases" not in {p.name for p in FAMILY.iterdir()}


@needs_bash
@pytest.mark.parametrize("script", SHELL_SCRIPTS,
                         ids=lambda p: p.relative_to(FAMILY).as_posix())
def test_every_shell_script_parses(script):
    result = subprocess.run([BASH, "-n", script.as_posix()], capture_output=True,
                            text=True)
    assert result.returncode == 0, result.stderr


@pytest.fixture
def bare_family(tmp_path):
    """A family with its marker files and one case, but no hpc_env.sh."""
    family = tmp_path / "fam"
    case = family / "Cases" / "Kn1_np1x"
    case.mkdir(parents=True)
    (case / "case.yaml").write_text("{}\n", encoding="utf-8")
    for name in ("study.yaml", "generate_cases.py", "AllpostCases"):
        (family / name).write_text("", encoding="utf-8")
    return family


@needs_bash
@pytest.mark.parametrize("script", JOB_SCRIPTS, ids=lambda p: p.name)
def test_a_job_fails_loudly_without_its_directory_or_environment(
        script, bare_family, tmp_path):
    env = {key: value for key, value in os.environ.items()
           if not key.startswith("PBS_")}
    run = lambda workdir: subprocess.run(  # noqa: E731
        [BASH, script.as_posix()], capture_output=True, text=True, cwd=tmp_path,
        env=env if workdir is None else {**env, "PBS_O_WORKDIR": workdir.as_posix()})

    result = run(None)
    assert result.returncode == 1
    assert "PBS_O_WORKDIR is not set" in result.stderr

    empty = tmp_path / "elsewhere"
    empty.mkdir()
    result = run(empty)
    assert result.returncode == 1
    assert ("no case.yaml" in result.stderr
            or "not a generated case" in result.stderr
            or "not the family directory" in result.stderr), result.stderr

    workdir = (bare_family if script in FAMILY_SCRIPTS
               else bare_family / "Cases" / "Kn1_np1x")
    result = run(workdir)
    assert result.returncode == 1
    assert "no hpc_env.sh" in result.stderr


@needs_bash
def test_hpc_env_refuses_set_u(tmp_path):
    result = subprocess.run(
        [BASH, "-c", f'set -u; source "{(FAMILY / "hpc_env.sh").as_posix()}" || echo "refused $?"'],
        capture_output=True, text=True, cwd=tmp_path)
    assert "set -u" in result.stderr
    assert result.stdout == "refused 1\n"           # before any module is touched


@needs_bash
def test_hpc_env_hides_the_callers_arguments_from_the_openfoam_bashrc(tmp_path):
    """The site paths swapped for stubs; the argument handling is the file's own."""
    stubs = tmp_path / "stubs"
    (stubs / "lmod" / "lmod" / "init").mkdir(parents=True)
    (stubs / "lmod" / "lmod" / "init" / "bash").write_text("", encoding="utf-8")
    for name in ("setup-env.sh", "activate"):
        (stubs / name).write_text("", encoding="utf-8")
    seen = tmp_path / "bashrc-saw"
    (stubs / "bashrc").write_text(f'echo "$#:$*" > "{seen.as_posix()}"\n', encoding="utf-8")

    text = read(FAMILY / "hpc_env.sh")
    for site, stub in (
            ("/home2/appman/spack/share/spack/setup-env.sh", "setup-env.sh"),
            ("/home2/torresfa/scratch/OpenFOAM/OpenFOAM-v2512/etc/bashrc", "bashrc"),
            ("/home2/torresfa/scratch/OpenFOAM/rarefiedPlume/.venv/bin/activate",
             "activate")):
        assert site in text
        text = text.replace(site, (stubs / stub).as_posix())
    env_file = tmp_path / "hpc_env.sh"
    env_file.write_text(text, encoding="utf-8")

    script = (f'set -eo pipefail; spack() {{ echo "{stubs.as_posix()}"; }}; module() {{ :; }}; '
              f'source "{env_file.as_posix()}"; echo "args:$#:$*"; echo "pv:$PVPYTHON_PATH"')
    result = subprocess.run([BASH, "-c", script, "job", "case.yaml", "other"],
                            capture_output=True, text=True, cwd=tmp_path,
                            env={**os.environ, "PVPYTHON_PATH": "/my/pvpython"})
    assert result.returncode == 0, result.stderr
    assert read(seen).strip() == "0:"
    assert "args:2:case.yaml other" in result.stdout
    assert "pv:/my/pvpython" in result.stdout
    assert "<<< Completed: load OpenFOAM environment" in result.stdout


@needs_bash
@needs_flock
def test_run_case_labels_runs_records_and_then_does_not_record_again(tmp_path):
    family = tmp_path / "fam"
    case = write_case(family / "Cases" / "Kn0p001_np1x", 1.0e-3)
    (family / "study.yaml").write_text("", encoding="utf-8")
    (family / "hpc_env.sh").write_text("printf 'stand-in environment\\n'\n",
                                       encoding="utf-8")
    shutil.copy2(FAMILY / "baseCase" / "runCase.sh", case / "runCase.sh")
    write_control(case, end=0.5, delta=0.001)
    (case / "Allrun").write_text(
        "#!/bin/bash\necho \"stand-in Allrun $#\"\nmkdir -p 0.5\n",
        encoding="utf-8")
    env = {**python_env(tmp_path), "PBS_O_WORKDIR": str(case)}

    result = subprocess.run([BASH, (case / "runCase.sh").as_posix()], capture_output=True,
                            text=True, cwd=tmp_path, env=env, timeout=600)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "stand-in environment" in result.stdout
    log = read(case / "output.log")
    assert log.startswith("This case was run by its runCase.sh job.")
    assert log.index("UNSUITABLE FOR PRODUCTION") < log.index("stand-in Allrun 0")
    assert (case / hpc.MARKER).is_file()
    model = family / "results" / "cost-model.yaml"
    measurements = yaml.safe_load(read(model))["measurements"]
    assert [m["case"] for m in measurements] == ["Kn0p001_np1x"]
    assert measurements[0]["production_status"] == feasibility.STATUS_UNSUITABLE

    model.unlink()
    result = subprocess.run([BASH, (case / "runCase.sh").as_posix()], capture_output=True,
                            text=True, cwd=tmp_path, env=env, timeout=600)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "had already reached endTime" in result.stdout
    assert not model.exists()


# --------------------------------------------------------------------------- #
# AllpostCases --no-analyse
# --------------------------------------------------------------------------- #

@needs_bash
def test_allpostcases_no_analyse_skips_the_study_analysis(tmp_path):
    family = tmp_path / "fam"
    (family / "Cases_full_suite" / "Kn1_np1x").mkdir(parents=True)
    shutil.copy2(FAMILY / "AllpostCases", family / "AllpostCases")
    shutil.copy2(FAMILY / "study.yaml", family / "study.yaml")
    (family / "manifest-full-suite.yaml").write_text(yaml.safe_dump({
        "meta": {"knudsen_suite": policy.FULL_SUITE}, "cases": [
        {"name": "Kn1_np1x", "path": "Cases_full_suite/Kn1_np1x", "Kn": 1.0,
         "feasibility": {"status": "production-candidate"}}]}), encoding="utf-8")
    ran = family / "analysed"
    (family / "analyse.py").write_text(
        "from pathlib import Path\nPath('analysed').write_text('yes')\n",
        encoding="utf-8")
    env = python_env(tmp_path)
    run = lambda *flags: subprocess.run(  # noqa: E731
        [BASH, (family / "AllpostCases").as_posix(), *flags], capture_output=True,
        text=True, cwd=tmp_path, env=env, timeout=300)

    result = run("--case", "Kn1_np1x", "--no-viz", "--no-analyse")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "no solver output yet" in result.stdout
    assert not ran.exists()

    result = run("--case", "Kn1_np1x", "--no-viz")
    assert result.returncode == 0, result.stdout + result.stderr
    assert ran.is_file()

    result = run("--analyse-only", "--no-analyse")
    assert result.returncode == 1
    assert "opposite" in result.stderr


@needs_bash
def test_cleanup_dry_run_scopes_active_suite_and_preserves_legacy_results(tmp_path):
    family = tmp_path / "family"
    family.mkdir()
    shutil.copy2(FAMILY / "AllcleanCases", family / "AllcleanCases")
    shutil.copy2(FAMILY / "study.yaml", family / "study.yaml")
    active = family / "Cases_full_suite" / "Kn100_np1x" / "sentinel"
    legacy = family / "Cases" / "Kn1_np1x" / "results" / "sentinel"
    shared = family / "results" / "cost-model.yaml"
    for item in (active, legacy, shared):
        item.parent.mkdir(parents=True, exist_ok=True)
        item.write_text("preserved", encoding="utf-8")
    result = subprocess.run([BASH, (family / "AllcleanCases").as_posix(), "--purge", "--dry-run"],
                            capture_output=True, text=True, cwd=tmp_path,
                            env=python_env(tmp_path), timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Cases_full_suite" in result.stdout and "manifest-full-suite.yaml" in result.stdout
    assert "--dry-run: nothing removed" in result.stdout
    assert all(item.read_text(encoding="utf-8") == "preserved" for item in (active, legacy, shared))


# --------------------------------------------------------------------------- #
# plumetools.cai2012ext.hpc
# --------------------------------------------------------------------------- #

def test_guard_refuses_an_excluded_kn(tmp_path, capsys):
    case = write_case(tmp_path / "Kn0p01_np1x", 0.01)
    assert hpc.guard(case) is not None
    assert hpc.main(["guard", str(case)]) == 2
    assert "Cai's own Knudsen" in capsys.readouterr().err


def test_guard_accepts_an_extension_kn(tmp_path):
    case = write_case(tmp_path / "Kn1_np1x", 1.0)
    assert hpc.guard(case) is None
    assert hpc.main(["guard", str(case)]) == 0


def test_a_directory_without_a_case_is_refused(tmp_path, capsys):
    assert hpc.main(["guard", str(tmp_path)]) == 1
    assert "no case.yaml" in capsys.readouterr().err


def test_label_marks_a_dense_case_non_production(tmp_path, capsys):
    case = write_case(tmp_path / "Kn0p001_np1x", 1.0e-3)
    assert hpc.main(["label", str(case)]) == 0
    marker = (case / hpc.MARKER).read_text(encoding="utf-8")
    assert "UNSUITABLE FOR PRODUCTION" in marker
    unresolved = [line for line in marker.splitlines() if line.startswith("  - ")]
    assert {line[4:].split(":")[0] for line in unresolved} == {
        "mesh_budget", "cell_over_mean_free_path", "time_step_collision_time"}
    assert capsys.readouterr().out.startswith(marker)


def test_label_leaves_a_candidate_unmarked(tmp_path, capsys):
    case = write_case(tmp_path / "Kn10_np1x", 10.0)
    assert hpc.main(["label", str(case)]) == 0
    assert not (case / hpc.MARKER).exists()
    assert feasibility.STATUS_CANDIDATE in capsys.readouterr().out


def test_finished_reads_processor_times_with_half_a_step_of_slack(tmp_path):
    case = tmp_path / "case"
    case.mkdir()
    (case / "case.yaml").write_text("{}\n", encoding="utf-8")
    assert not hpc.reached_end(case)                  # no controlDict
    write_control(case, end=0.0628, delta=1.0e-4)
    (case / "0").mkdir()
    assert not hpc.reached_end(case)                  # time 0 is not a result
    (case / "processor0" / "0.06276").mkdir(parents=True)
    assert hpc.reached_end(case)                      # 0.4 dt short: finished
    shutil.rmtree(case / "processor0")
    (case / "processor0" / "0.06274").mkdir(parents=True)
    assert not hpc.reached_end(case)                  # 0.6 dt short: not yet
    assert hpc.main(["finished", str(case)]) == 1
    (case / "0.0628").mkdir()
    assert hpc.main(["finished", str(case)]) == 0     # reconstructed counts too


@needs_flock
def test_record_replaces_the_cases_previous_entry(tmp_path):
    family = tmp_path / "fam"
    case = write_case(family / "Cases" / "Kn1_np1x", 1.0)
    (family / "study.yaml").write_text("", encoding="utf-8")
    model = family / hpc.COST_MODEL
    other = {"case": "Kn10_np1x", "seconds": 50.0, "parcels": 2.0e6,
             "cells": 1.3e6, "steps": 45924}
    feasibility.record_measurement(model, other)

    assert hpc.main(["record", str(case), "100"]) == 0
    assert hpc.main(["record", str(case), "250"]) == 0
    measurements = yaml.safe_load(read(model))["measurements"]
    assert [m["case"] for m in measurements] == ["Kn10_np1x", "Kn1_np1x"]
    mine = measurements[1]
    assert mine["seconds"] == 250.0
    assert mine["path"] == "Cases/Kn1_np1x"
    assert mine["Kn"] == 1.0
    assert mine["cells"] == 1_310_720
    assert mine["production_status"] == feasibility.STATUS_CANDIDATE
    assert read(model).startswith("# GENERATED by each case's runCase.sh")


def test_record_outside_a_family_is_not_timed(tmp_path, capsys):
    case = write_case(tmp_path / "loose" / "Kn1_np1x", 1.0)
    assert hpc.main(["record", str(case), "100"]) == 0
    assert "not timed" in capsys.readouterr().out
    assert not list(tmp_path.rglob("cost-model.yaml"))


WRITER = """
import sys, time
from pathlib import Path
from plumetools.cai2012ext import hpc
path, tag, flags = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
(flags / f"ready-{tag}").touch()
while not (flags / "go").exists():
    time.sleep(0.001)
for i in range(int(sys.argv[4])):
    hpc.record_entry(path, {"case": f"{tag}{i:03d}", "seconds": 100.0 + i,
                            "parcels": 1.0e6 * (1 + i), "cells": 1.0e6 * (2 + i % 3),
                            "steps": 1000})
"""


@needs_flock
def test_concurrent_records_lose_no_entry(tmp_path):
    path = tmp_path / "results" / "cost-model.yaml"
    flags = tmp_path / "flags"
    flags.mkdir()
    tags, per_writer = ("a", "b", "c"), 20
    env = python_env(tmp_path)
    writers = [subprocess.Popen([sys.executable, "-c", WRITER, str(path), tag,
                                 str(flags), str(per_writer)], env=env)
               for tag in tags]
    deadline = time.monotonic() + 300
    while not all((flags / f"ready-{tag}").exists() for tag in tags):
        assert time.monotonic() < deadline, "the writers never started"
        assert all(w.poll() is None for w in writers), "a writer died"
        time.sleep(0.01)
    (flags / "go").touch()
    assert [w.wait(timeout=300) for w in writers] == [0] * len(tags)

    measurements = yaml.safe_load(read(path))["measurements"]
    assert sorted(m["case"] for m in measurements) == sorted(
        f"{tag}{i:03d}" for tag in tags for i in range(per_writer))
