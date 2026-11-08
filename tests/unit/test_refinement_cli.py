from __future__ import annotations

import json
import csv
import os
from pathlib import Path
import subprocess
import sys
import pytest

from _refinement import write_cartesian_case, tree_hash

REPO = Path(__file__).resolve().parents[2]


def run(cwd, *args):
    env = dict(os.environ, PYTHONPATH=str(REPO), PYTHONDONTWRITEBYTECODE="1")
    return subprocess.run([sys.executable, "-B", "-m", "plumetools.refinement", *map(str, args)],
                          cwd=cwd, env=env, capture_output=True, text=True, timeout=30)


def test_standalone_assessment_and_explicit_report_output(tmp_path):
    case = write_cartesian_case(tmp_path / "unrelated case")
    before = tree_hash(case)
    r = run(tmp_path, "assess", case, "--settings", case / "refinement.yaml")
    assert r.returncode == 0, r.stderr
    assert "needs_mesh" in r.stdout
    assert tree_hash(case) == before
    out = tmp_path / "assessment"
    r = run(tmp_path, "assess", case, "--settings", case / "refinement.yaml", "--output", out)
    assert r.returncode == 0, r.stderr
    report = json.loads((out / "assessment.json").read_text())
    assert report["summary"]["mesh_failed_cells"] == 1
    assert (out / "cells.csv").is_file()


def test_help_and_dry_run_write_nothing(tmp_path):
    case = write_cartesian_case(tmp_path / "case")
    out = tmp_path / "not created"
    before = tree_hash(tmp_path)
    assert run(tmp_path, "--help").returncode == 0
    result = run(tmp_path, "plan", case, "--settings", case / "refinement.yaml", "--output", out, "--dry-run")
    assert result.returncode == 0, result.stderr
    assert "needs_mesh" in result.stdout
    assert tree_hash(tmp_path) == before
    assert not out.exists()


def test_plan_command_exports_selected_cell_and_source_data(tmp_path):
    case = write_cartesian_case(tmp_path / "case")
    out = tmp_path / "plan"
    r = run(tmp_path, "plan", case, "--settings", case / "refinement.yaml", "--output", out)
    assert r.returncode == 0, r.stderr
    assert (out / "plan.json").is_file()
    with (out / "cells.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert float(rows[0]["target_edge_m"]) == pytest.approx(.6)


def test_apply_of_passing_plan_creates_no_mesh_or_solver_work(tmp_path):
    case = write_cartesian_case(tmp_path / "case", lambdas=(2, 2))
    artifact = tmp_path / "plan"
    planned = run(tmp_path, "plan", case, "--settings", case / "refinement.yaml", "--output", artifact)
    assert planned.returncode == 0, planned.stderr
    out = tmp_path / "not created"
    result = run(tmp_path, "apply", case, "--plan", artifact / "plan.json", "--backend", "hex", "--output", out)
    assert result.returncode == 0, result.stderr
    assert not out.exists()
