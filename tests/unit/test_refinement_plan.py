from __future__ import annotations

from dataclasses import replace
import importlib
import shutil

import numpy as np
import pytest

from _refinement import unit_inputs, write_cartesian_case
from plumetools.refinement.assess import assess
from plumetools.refinement.config import load_settings
from plumetools.refinement.foam import read_snapshot


def api():
    return importlib.import_module("plumetools.refinement.plan")


def case_plan(case):
    snap = read_snapshot(case, load_settings(case / "refinement.yaml"))
    a = assess(snap.geometry, snap.measured, snap.species)
    return api().plan_refinement(snap.geometry, snap.measured, a)


def test_plan_uses_measured_targets_and_does_not_refine_passing_cells():
    g, m, gas = unit_inputs(edges=[1, .5], lambdas=[.75, 2], parcels=[40, 40])
    p = api().plan_refinement(g, m, assess(g, m, gas))
    assert p.selected_cells.tolist() == [0]
    assert p.target_edge_m == pytest.approx([.6, 1.6])
    assert p.recommended_weight is None
    good = replace(g, max_edge_m=np.array([.5, .5]))
    passing = api().plan_refinement(good, m, assess(good, m, gas))
    assert passing.action == "criteria_satisfied"
    assert passing.selected_cells.size == 0


def test_serialized_plan_preserves_targets_and_source_geometry(tmp_path):
    case = write_cartesian_case(tmp_path / "source")
    p = case_plan(case)
    file = api().write_plan(p, tmp_path / "artifact")
    loaded = api().read_plan(file)
    assert loaded.mesh_sha256 == p.mesh_sha256
    assert loaded.selected_cells.tolist() == [0]
    assert loaded.source_geometry.centres_m.tolist() == p.source_geometry.centres_m.tolist()
    assert loaded.target_edge_m == pytest.approx([.6, 1.6])
    assert (file.parent / "targets.npz").is_file()
    assert (file.parent / "cells.csv").is_file()
    assert api().write_plan(p, file.parent) == file


def test_plan_directory_is_accepted_by_the_documented_workflow(tmp_path):
    source = write_cartesian_case(tmp_path/'source')
    original = case_plan(source)
    file = api().write_plan(original,tmp_path/'targets')
    loaded = api().read_plan(file.parent)
    assert api().plan_identity(loaded) == api().plan_identity(original)


@pytest.mark.parametrize("kind", ["field", "mesh"])
def test_changed_source_invalidates_plan(kind, tmp_path):
    case = write_cartesian_case(tmp_path / "source")
    p = case_plan(case)
    path = case / ("1/TObserved" if kind == "field" else "constant/polyMesh/points")
    path.write_text(path.read_text().replace("300", "301") if kind == "field"
                    else path.read_text().replace("(0 0 0)", "(-0.1 0 0)"))
    with pytest.raises(ValueError, match="changed|stale"):
        api().verify_sources(case, p)


def test_plan_is_relocatable(tmp_path):
    source = write_cartesian_case(tmp_path / "source")
    p = case_plan(source)
    copied = tmp_path / "copied"; shutil.copytree(source, copied)
    api().verify_sources(copied, p)


def test_low_population_gets_weight_estimate_but_no_guessed_mesh_plan():
    g, m, gas = unit_inputs(edges=[2], lambdas=[1], parcels=[10])
    p = api().plan_refinement(g, m, assess(g, m, gas), particle_weight=100)
    assert p.action == "needs_parcels"
    assert p.selected_cells.size == 0
    assert p.recommended_weight == pytest.approx(50)
    unknown = api().plan_refinement(g, m, assess(g, m, gas))
    assert unknown.recommended_weight is None


def test_unresolved_cells_prevent_executable_mesh_plan():
    g, m, gas = unit_inputs(edges=[2], lambdas=[1], parcels=[40])
    m = replace(m, temperature_K=np.array([np.nan]))
    p = api().plan_refinement(g, m, assess(g, m, gas), particle_weight=100)
    assert p.action == "unresolved"
    assert p.selected_cells.size == 0
    assert p.recommended_weight is None


def test_assessment_from_other_measurements_is_rejected():
    g, m, gas = unit_inputs(edges=[2], lambdas=[1], parcels=[40])
    a = assess(g, m, gas)
    with pytest.raises(ValueError, match="assessment"):
        api().plan_refinement(g, replace(m, source_sha256="different"), a)


def test_artifact_collision_and_corrupt_arrays_are_rejected(tmp_path):
    case = write_cartesian_case(tmp_path / "source")
    p = case_plan(case)
    out = tmp_path / "unrelated"; out.mkdir(); (out / "keep.txt").write_text("keep")
    with pytest.raises(ValueError, match="exist"):
        api().write_plan(p, out)
    file = api().write_plan(p, tmp_path / "artifact")
    (file.parent / "targets.npz").write_bytes(b"bad data")
    with pytest.raises(ValueError, match="hash|corrupt|checksum"):
        api().read_plan(file)
