from __future__ import annotations

import importlib
import math
import numpy as np
import pytest
import yaml

from _refinement import write_cartesian_case, tree_hash, write_field


def read(case, *, time=None):
    cfg = importlib.import_module("plumetools.refinement.config").load_settings(case / "refinement.yaml")
    return importlib.import_module("plumetools.refinement.foam").read_snapshot(case, cfg, time=time)


def test_custom_fields_outside_study(tmp_path):
    from plumetools.refinement.assess import assess
    case = write_cartesian_case(tmp_path / "unrelated case")
    before = tree_hash(case)
    snap = read(case)
    assert snap.geometry.max_edge_m == pytest.approx([1, 1])
    assert snap.measured.temperature_K == pytest.approx([300, 300])
    assert assess(snap.geometry, snap.measured, snap.species).mesh_status.tolist() == [0, 1]
    assert not (case / "case.yaml").exists()
    assert tree_hash(case) == before


@pytest.mark.parametrize("profile", ["standard", "hystrath"])
def test_native_profiles_use_the_correct_temperature_and_count(profile, tmp_path):
    case = write_cartesian_case(tmp_path / "case", profile=profile)
    snap = read(case)
    assert snap.measured.temperature_K == pytest.approx([300, 300])
    assert snap.measured.mean_parcels.tolist() == [40, 40]
    assert snap.measured.window["source"] == ("verified" if profile == "standard" else "declared")


def test_uniform_fields_are_broadcast_only_when_declared_uniform(tmp_path):
    case = write_cartesian_case(tmp_path / "case")
    write_field(case / "1/TObserved", [350], [0, 0, 0, 1, 0, 0, 0], uniform=True)
    assert read(case).measured.temperature_K.tolist() == [350, 350]
    write_field(case / "1/TObserved", [350], [0, 0, 0, 1, 0, 0, 0])
    with pytest.raises(ValueError, match="cells|mesh"):
        read(case)


def test_sampling_window_required_and_consistent(tmp_path):
    from plumetools.refinement.assess import assess
    case = write_cartesian_case(tmp_path / "case")
    path = case / "refinement.yaml"
    raw = yaml.safe_load(path.read_text()); raw["sampling"] = {}
    path.write_text(yaml.safe_dump(raw))
    snap = read(case)
    assert assess(snap.geometry, snap.measured, snap.species).summary["action"] == "unresolved"
    standard = write_cartesian_case(tmp_path / "standard", profile="standard")
    state = standard / "1/uniform/functionObjects/functionObjectProperties"
    state.write_text(state.read_text().replace("rhoM\n{ totalIter 10; totalTime 1; }", "rhoM\n{ totalIter 9; totalTime 0.9; }"))
    with pytest.raises(ValueError, match="different windows"):
        read(standard)


def test_wrong_field_dimensions_are_not_silently_interpreted(tmp_path):
    case = write_cartesian_case(tmp_path / "case")
    write_field(case / "1/TObserved", [300, 300], [0]*7)
    with pytest.raises(ValueError, match="dimensions"):
        read(case)


@pytest.mark.parametrize("kind", ["binary", "decomposed", "dynamic"])
def test_unsupported_input_is_read_only(kind, tmp_path):
    case = write_cartesian_case(tmp_path / "case")
    if kind == "binary":
        p = case / "1/TObserved"; p.write_text(p.read_text().replace("format ascii", "format binary"))
    elif kind == "decomposed":
        (case / "processor0").mkdir()
        (case / "constant/polyMesh").rename(case / "processor0/polyMesh")
    else:
        (case / "constant/dynamicMeshDict").write_text("dynamicFvMesh dynamicRefineFvMesh;\n")
    before = tree_hash(case)
    with pytest.raises(ValueError, match="ASCII|reconstruct|static|dynamic"):
        read(case)
    assert tree_hash(case) == before


def test_missing_geometry_explains_how_to_export_it(tmp_path):
    case = write_cartesian_case(tmp_path / "case")
    (case / "1/V").unlink()
    with pytest.raises(ValueError, match="writeCellVolumes"):
        read(case)


def test_cell_geometry_field_names_can_be_overridden(tmp_path):
    case = write_cartesian_case(tmp_path / "case")
    (case / "1/C").rename(case / "1/centresCustom")
    (case / "1/V").rename(case / "1/volumesCustom")
    path = case / "refinement.yaml"
    raw = yaml.safe_load(path.read_text())
    raw["fields"].update(cell_centres="centresCustom", cell_volumes="volumesCustom")
    path.write_text(yaml.safe_dump(raw))
    assert read(case).geometry.volumes_m3.tolist() == [1, 1]


def test_commented_weighting_models_do_not_disable_constant_weight(tmp_path):
    case = write_cartesian_case(tmp_path / "case", lambdas=(1, 1))
    weight = 1 / (math.sqrt(2) * math.pi * 40)
    (case / "constant/dsmcProperties").write_text(
        f"// coordinateSystem dsmcAxisymmetric;\n// typeIdList (Ar N2);\nnEquivalentParticles {weight:.17g};\n")
    assert read(case).particle_weight == pytest.approx(weight)


def test_nonpositive_configured_particle_weight_is_rejected(tmp_path):
    case = write_cartesian_case(tmp_path / "case", lambdas=(1, 1))
    (case / "constant/dsmcProperties").write_text("nEquivalentParticles -1;\n")
    with pytest.raises(ValueError, match="nEquivalentParticles"):
        read(case)


@pytest.mark.parametrize('role,old,new', [
    ('number_density','rhoNMean','nObserved'),
    ('mean_parcels','dsmcRhoNMean','countObserved'),
])
def test_renamed_fields_need_complete_window_evidence(tmp_path,role,old,new):
    from plumetools.refinement.assess import assess
    case = write_cartesian_case(tmp_path/'case',profile='standard',lambdas=(2,2))
    (case/'1'/old).rename(case/'1'/new)
    settings = case/'refinement.yaml'
    raw = yaml.safe_load(settings.read_text()); raw['fields'][role] = new
    settings.write_text(yaml.safe_dump(raw))
    before = tree_hash(case)
    snap = read(case)
    assert snap.measured.window.get('source') != 'verified'
    assert assess(snap.geometry,snap.measured,snap.species).summary['action'] == 'unresolved'
    assert tree_hash(case) == before


def test_explicit_window_for_renamed_field_must_match_available_native_data(tmp_path):
    case = write_cartesian_case(tmp_path/'case',profile='standard')
    (case/'1/rhoNMean').rename(case/'1/nObserved')
    settings = case/'refinement.yaml'
    raw = yaml.safe_load(settings.read_text())
    raw['fields']['number_density'] = 'nObserved'
    raw['sampling'] = {'start_s':0,'end_s':1,'source':'declared','basis':'mean_per_cell','n_steps':10}
    settings.write_text(yaml.safe_dump(raw))
    assert read(case).measured.window['source'] == 'declared'
    raw['sampling']['start_s'] = .5
    settings.write_text(yaml.safe_dump(raw))
    with pytest.raises(ValueError,match='window'):
        read(case)


def test_renamed_field_native_metadata_is_checked_when_present(tmp_path):
    case = write_cartesian_case(tmp_path/'case',profile='standard')
    (case/'1/rhoNMean').rename(case/'1/nObserved')
    settings = case/'refinement.yaml'
    raw = yaml.safe_load(settings.read_text()); raw['fields']['number_density'] = 'nObserved'
    settings.write_text(yaml.safe_dump(raw))
    state = case/'1/uniform/functionObjects/functionObjectProperties'
    state.write_text(state.read_text()+'\nnObserved { totalIter 5; totalTime 0.5; }\n')
    with pytest.raises(ValueError,match='different windows'):
        read(case)
