from __future__ import annotations

import importlib
from pathlib import Path

import pytest

from _refinement import write_cartesian_case, tree_hash
from plumetools.refinement.assess import assess
from plumetools.refinement.config import load_settings
from plumetools.refinement.foam import read_snapshot
from plumetools.refinement.plan import plan_refinement


def setup_case(tmp_path, *, lambdas=(.75, 2)):
    source = write_cartesian_case(tmp_path / "source", lambdas=lambdas)
    snap = read_snapshot(source, load_settings(source / "refinement.yaml"))
    plan = plan_refinement(snap.geometry, snap.measured, assess(snap.geometry, snap.measured, snap.species))
    backend = importlib.import_module("plumetools.refinement.openfoam")
    return source, plan, backend


def test_passing_plan_needs_no_tools_or_new_directory(tmp_path):
    source, plan, module = setup_case(tmp_path, lambdas=(2, 2))
    before = tree_hash(source)
    output = tmp_path / "not created"
    result = module.HexRefiner().refine(source, plan, output)
    assert result.destination is None
    assert result.requires_flow_rerun is False
    assert not output.exists()
    assert tree_hash(source) == before


def test_existing_output_and_nested_output_are_refused(tmp_path):
    source, plan, module = setup_case(tmp_path)
    output = tmp_path / "existing"; output.mkdir(); (output / "keep").write_text("keep")
    with pytest.raises(ValueError, match="exists"):
        module.HexRefiner().refine(source, plan, output)
    with pytest.raises(ValueError, match="source"):
        module.HexRefiner().refine(source, plan, source / "child")
    assert (output / "keep").read_text() == "keep"


def test_stale_plan_is_refused_before_work(tmp_path):
    source, plan, module = setup_case(tmp_path)
    p = source / "1/TObserved"; p.write_text(p.read_text().replace("300", "301"))
    before = tree_hash(source)
    with pytest.raises(ValueError, match="changed|stale"):
        module.HexRefiner().refine(source, plan, tmp_path / "out")
    assert tree_hash(source) == before


def test_budget_does_not_silently_coarsen_the_request(tmp_path):
    source, plan, module = setup_case(tmp_path)
    with pytest.raises(ValueError, match="budget|cells"):
        module.HexRefiner(max_cells=8).refine(source, plan, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_failed_backend_keeps_source_and_destination_intact(tmp_path, monkeypatch):
    source, plan, module = setup_case(tmp_path)
    before = tree_hash(source)
    monkeypatch.setenv("WM_PROJECT_VERSION", "v2512")
    monkeypatch.setattr(module.shutil, "which", lambda name: name)
    def fail(command, cwd, log):
        Path(log).write_text("intentional tool failure")
        raise RuntimeError("intentional tool failure")
    monkeypatch.setattr(module, "_run", fail)
    with pytest.raises(RuntimeError, match="intentional tool failure"):
        module.HexRefiner().refine(source, plan, tmp_path / "out")
    assert not (tmp_path / "out").exists()
    assert tree_hash(source) == before


def test_non_3d_topology_gets_an_actionable_error(tmp_path):
    source = write_cartesian_case(tmp_path / "source")
    p = source / "constant/polyMesh/boundary"
    p.write_text(p.read_text().replace("type wall", "type empty"))
    snap = read_snapshot(source, load_settings(source / "refinement.yaml"))
    plan = plan_refinement(snap.geometry, snap.measured, assess(snap.geometry, snap.measured, snap.species))
    module = importlib.import_module("plumetools.refinement.openfoam")
    before = tree_hash(source)
    with pytest.raises(ValueError, match="3D|three-dimensional|topology"):
        module.HexRefiner().refine(source, plan, tmp_path / "out")
    assert tree_hash(source) == before


def test_quality_failure_with_zero_process_status_is_rejected(tmp_path):
    module = importlib.import_module("plumetools.refinement.openfoam")
    log = tmp_path / "checkMesh.log"
    log.write_text("Mesh check output\nFailed 1 mesh checks.\nEnd\n")
    with pytest.raises(ValueError, match="mesh check|quality"):
        module._require_mesh_ok(log)
    log.write_text("Mesh check output\nMesh OK.\nEnd\n")
    module._require_mesh_ok(log)
