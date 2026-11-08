"""Actual portable refinement on a case without any Cai or radial metadata."""
from __future__ import annotations

import importlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "unit"))
from _refinement import write_cartesian_case, tree_hash
from plumetools.refinement.assess import assess
from plumetools.refinement.config import load_settings
from plumetools.refinement.foam import read_snapshot
from plumetools.refinement.geometry import read_geometry
from plumetools.refinement.plan import plan_refinement
from plumetools.foamio.volfield import read_internal

pytestmark = pytest.mark.needs_openfoam


@pytest.fixture
def tools():
    missing = [n for n in ("refineHexMesh", "postProcess", "checkMesh") if not shutil.which(n)]
    if missing:
        pytest.skip("requires OpenFOAM v2512: " + ", ".join(missing))


def build(tmp_path, lambdas):
    source = write_cartesian_case(tmp_path / "unrelated case", lambdas=lambdas)
    snap = read_snapshot(source, load_settings(source / "refinement.yaml"))
    p = plan_refinement(snap.geometry, snap.measured, assess(snap.geometry, snap.measured, snap.species))
    return source, p


def test_refine_unrelated_cartesian_case(tmp_path, tools):
    source, plan = build(tmp_path, (.75, 2))
    before = tree_hash(source)
    output = tmp_path / "refined mesh"
    result = importlib.import_module("plumetools.refinement.openfoam").HexRefiner().refine(source, plan, output)
    assert result.n_cells == 9
    assert result.requires_flow_rerun is True
    assert result.target_status == "inherited_targets_met"
    g = read_geometry(output, time="0")
    assert g.volumes_m3.sum() == pytest.approx(2.0)
    assert sorted(g.volumes_m3) == pytest.approx([.125]*8 + [1])
    assert tree_hash(source) == before
    assert not (output / "1").exists()
    assert not (output / "0/lagrangian").exists()
    metadata = json.loads((output / "refinement-result.json").read_text())
    assert metadata["requires_flow_rerun"] is True


def test_multiple_passes_preserve_targets_and_extend_for_conformity(tmp_path, tools):
    source, plan = build(tmp_path, (.3, 2))
    output = tmp_path / "refined"
    result = importlib.import_module("plumetools.refinement.openfoam").HexRefiner().refine(source, plan, output)
    g = read_geometry(output, time="0")
    target = read_internal(output / "0/refinementTargetEdge").over(result.n_cells, "cells")
    assert np.all(g.max_edge_m <= target * (1 + 1e-12))
    assert g.volumes_m3.sum() == pytest.approx(2.0)
    metadata = json.loads((output / "refinement-result.json").read_text())
    assert metadata["stages"] == 3
    assert any(step["conformity_cells"] > 0 for step in metadata["history"])


def test_stage_limit_never_publishes_an_incomplete_mesh(tmp_path, tools):
    source, plan = build(tmp_path, (.3, 2))
    before = tree_hash(source)
    module = importlib.import_module("plumetools.refinement.openfoam")
    with pytest.raises(ValueError, match="stage|limit"):
        module.HexRefiner(max_stages=1).refine(source, plan, tmp_path / "not published")
    assert not (tmp_path / "not published").exists()
    assert tree_hash(source) == before


def test_native_hex_refinement_under_spaced_ancestor(tmp_path,tools):
    source,plan = build(tmp_path,(.75,2))
    output = tmp_path/'ancestor with spaces'/'mesh'
    result = importlib.import_module('plumetools.refinement.openfoam').HexRefiner().refine(source,plan,output)
    assert result.n_cells == 9
    assert read_geometry(output,time='0').volumes_m3.sum() == pytest.approx(2.)
