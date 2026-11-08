from dataclasses import replace
from pathlib import Path

import pytest
import yaml


STUDY = Path(__file__).resolve().parents[2] / 'cases/cai2012-radial-resolution/study.yaml'


def test_default_study_inherits_only_reference_physics():
    from plumetools.radialstudy.config import load_study
    cfg = load_study(STUDY)
    assert cfg.knudsen_numbers == (100., 10., 1., .1, .01, .001, .0001)
    assert cfg.reference.nozzle.diameter_m == .2
    assert cfg.reference.exit.T0_K == 300
    assert cfg.reference.exit.speed_ratio == 2
    assert cfg.reference.gas.species_name == 'Ar'
    assert cfg.mesh.predicted_edge_over_mfp == 1
    assert cfg.targets.edge_over_mfp == .8
    assert cfg.targets.min_mean_parcels == 20
    assert cfg.run.cores == 1


def test_adaptive_defaults_and_initial_cells():
    from plumetools.radialstudy.config import load_study
    from plumetools.radialstudy.physics import source_for_kn,predict_breakdown
    from plumetools.radialstudy.mesh import initial_mesh
    cfg = load_study(STUDY)
    assert cfg.domain.angle_deg==2.5 and cfg.domain.outer_over_D==10
    assert cfg.mesh.initial_cells==100 and cfg.controller.enabled
    assert cfg.controller.min_windows==4 and cfg.controller.stability_relative==.05
    assert cfg.controller.retry_parcels==25 and cfg.controller.max_retries==5
    m = source_for_kn(cfg,.0001)
    mesh = initial_mesh(m)
    assert mesh.n_cells==100 and mesh.transverse==1
    assert mesh.boundaries[-1]/m.source.diameter_m==pytest.approx(10)
    assert predict_breakdown(m)['status']=='NOT_REACHED_IN_DOMAIN'


@pytest.mark.parametrize('section,value',[
    ('mesh',{'initial_cells':0}),('mesh',{'initial_cells':True}),
    ('controller',{'enabled':'yes'}),('controller',{'min_windows':3}),
    ('controller',{'stability_relative':0}),('controller',{'max_retries':-1})])
def test_invalid_controller_settings(tmp_path,section,value):
    from plumetools.radialstudy.config import load_study
    data = yaml.safe_load(STUDY.read_text())
    data['reference_case'] = str(STUDY.parent.parent/'cai2012-kn-extension/baseCase/case.yaml')
    data[section] = value
    path = tmp_path/'study.yaml'
    path.write_text(yaml.safe_dump(data))
    with pytest.raises(ValueError):
        load_study(path)


@pytest.mark.parametrize('section,value', [('surprise', {}), ('mesh', {'max_cells': 0}),
    ('knudsen_numbers', [0]), ('time', {'blocks': 1}), ('run', {'cores': 2})])
def test_invalid_study_is_rejected(tmp_path, section, value):
    from plumetools.radialstudy.config import load_study
    data = yaml.safe_load(STUDY.read_text())
    data['reference_case'] = str(STUDY.parent.parent / 'cai2012-kn-extension/baseCase/case.yaml')
    data[section] = value
    p = tmp_path / 'study.yaml'
    p.write_text(yaml.safe_dump(data))
    with pytest.raises(ValueError):
        load_study(p)


def test_numerical_options_do_not_change_source_identity():
    from plumetools.radialstudy.config import load_study
    from plumetools.radialstudy.physics import source_for_kn
    c = load_study(STUDY)
    other = replace(c, targets=replace(c.targets, min_mean_parcels=40),
                    mesh=replace(c.mesh, predicted_edge_over_mfp=.9))
    assert source_for_kn(c, .001).physical_hash == source_for_kn(other, .001).physical_hash


def test_incompatible_physical_reference(tmp_path):
    from plumetools.radialstudy.config import load_study
    data = yaml.safe_load(STUDY.read_text())
    ref = tmp_path / 'reference.yaml'
    ref.write_text('model: cai2012_circular_plume\nexit:\n  characteristic_length: radius\n')
    data['reference_case'] = str(ref)
    p = tmp_path / 'study.yaml'
    p.write_text(yaml.safe_dump(data))
    with pytest.raises(ValueError, match='diameter'):
        load_study(p)
