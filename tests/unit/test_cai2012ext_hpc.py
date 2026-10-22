"""What a case's own PBS job needs from the extension.

``AllpostCases --no-analyse`` runs one case's post-processing chain without
the study-level analysis, which waits for every case. It is run here against a
stand-in ``analyse.py``. :mod:`plumetools.cai2012ext.hpc` -- the Kn guard, the
NON-PRODUCTION label, the finished test and the locked wall-clock record -- is
run on real ``case.yaml`` files.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest
import yaml

from plumetools.cai2012ext import feasibility, hpc

REPO = Path(__file__).resolve().parents[2]
FAMILY = REPO / "cases" / "cai2012-kn-extension"

needs_bash = pytest.mark.skipif(shutil.which("bash") is None, reason="no bash")
needs_flock = pytest.mark.skipif(sys.platform == "win32", reason="no fcntl")


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


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
    wrapper.write_text(f'#!/bin/sh\nexec "{sys.executable}" "$@"\n',
                       encoding="utf-8")
    wrapper.chmod(0o755)
    env = {key: value for key, value in os.environ.items()
           if not key.startswith("PBS_")}
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    env["PYTHONPATH"] = os.pathsep.join(
        filter(None, [str(REPO), env.get("PYTHONPATH")]))
    return env


# --------------------------------------------------------------------------- #
# AllpostCases --no-analyse
# --------------------------------------------------------------------------- #

@needs_bash
def test_allpostcases_no_analyse_skips_the_study_analysis(tmp_path):
    family = tmp_path / "fam"
    (family / "Cases" / "Kn1_np1x").mkdir(parents=True)
    shutil.copy2(FAMILY / "AllpostCases", family / "AllpostCases")
    (family / "manifest.yaml").write_text(yaml.safe_dump({"cases": [
        {"name": "Kn1_np1x", "path": "Cases/Kn1_np1x", "Kn": 1.0,
         "feasibility": {"status": "production-candidate"}}]}), encoding="utf-8")
    ran = family / "analysed"
    (family / "analyse.py").write_text(
        "from pathlib import Path\nPath('analysed').write_text('yes')\n",
        encoding="utf-8")
    env = python_env(tmp_path)
    run = lambda *flags: subprocess.run(  # noqa: E731
        ["bash", str(family / "AllpostCases"), *flags], capture_output=True,
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
