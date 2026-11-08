import shutil

import pytest

from test_radialstudy_run import prepared
from test_radialstudy_post import synthetic


def measured(tmp_path,monkeypatch,windows=4,**kwargs):
    from plumetools.radialstudy.mesh import read_radial_mesh
    from plumetools.foamio.times import time_directories
    path,record = prepared(tmp_path)
    synthetic(path,record,read_radial_mesh(path),monkeypatch,**kwargs)
    for _,name in time_directories(path)[windows:]:
        shutil.rmtree(path/name)
    return path


def test_four_stable_windows_accept_without_claiming_independence(tmp_path,monkeypatch):
    from plumetools.radialstudy.post import checkpoint_assessment,assess_case
    p = measured(tmp_path,monkeypatch)
    result,plan = checkpoint_assessment(p)
    assert result['status']=='criteria_satisfied'
    assert plan.action=='criteria_satisfied'
    assert result['sampling']['basis']=='adaptive_sizing'
    assert not result['sampling']['independence_supported']
    assert not result['sampling']['complete']
    assert assess_case(p)['status']=='unresolved'


@pytest.mark.parametrize('kwargs',[{'windows':3},{'drift':True},{'parcels':0}])
def test_unreliable_partial_data_cannot_pass_or_refine(tmp_path,monkeypatch,kwargs):
    from plumetools.radialstudy.post import checkpoint_assessment
    p = measured(tmp_path,monkeypatch,**kwargs)
    result,plan = checkpoint_assessment(p)
    assert result['status']=='unresolved' and plan is None


@pytest.mark.parametrize('kwargs,action',[
    ({'parcels':10},'needs_parcels'),({'density_scale':1e7},'needs_mesh'),
    ({'velocity':5000},'needs_timestep')])
def test_incremental_checks_select_measured_correction(tmp_path,monkeypatch,kwargs,action):
    from plumetools.radialstudy.post import checkpoint_assessment
    p = measured(tmp_path,monkeypatch,**kwargs)
    result,plan = checkpoint_assessment(p)
    assert result['status']==action and plan is not None
    if action=='needs_timestep':
        assert result['recommended_delta_t_s']>0
    else:
        assert plan.action==action
