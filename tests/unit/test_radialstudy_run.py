from dataclasses import replace
import json
from pathlib import Path

import pytest

from test_radialstudy_physics import model


def test_transit_integrates_radial_speed_and_seeds_are_distinct(tmp_path):
    from scipy.integrate import quad
    from plumetools.radialstudy.case import sizing,write_case
    from plumetools.radialstudy.mesh import initial_mesh
    a,b = model(100),model(10)
    ma,mb = initial_mesh(a),initial_mesh(b)
    _,s = sizing(a,ma)
    reference = quad(lambda x: 1/float(a.state(x)['U']),a.r_in,a.r_out,epsabs=1e-13)[0]
    assert s['transit_s'] == pytest.approx(reference,rel=1e-8)
    ra = write_case(tmp_path/'seed-a',a,ma)
    rb = write_case(tmp_path/'seed-b',b,mb)
    assert ra['seed'] != rb['seed']
    assert write_case(tmp_path/'seed-a',a,ma)['seed'] == ra['seed']


def prepared(tmp_path):
    from plumetools.radialstudy.case import write_case
    from plumetools.radialstudy.mesh import initial_mesh
    m = model(100)
    m = replace(m,config=replace(m.config,domain=replace(m.config.domain,outer_over_D=.82,angle_deg=.5)))
    from plumetools.radialstudy.physics import source_for_kn
    m = source_for_kn(m.config,100)
    p = tmp_path/'case with spaces'
    r = write_case(p,m,initial_mesh(m))
    return p,r


def test_environment_detection_is_actionable(monkeypatch):
    from plumetools.radialstudy.run import require_environment
    monkeypatch.delenv('WM_PROJECT_VERSION',raising=False)
    with pytest.raises(ValueError,match='source'):
        require_environment()
    monkeypatch.setenv('WM_PROJECT_VERSION','v2512')
    monkeypatch.setattr('shutil.which',lambda x: None)
    with pytest.raises(ValueError,match='dsmcFoam'):
        require_environment()


def test_production_solver_has_no_five_minute_utility_timeout(tmp_path,monkeypatch):
    from plumetools.radialstudy.run import _native_run
    seen = []
    monkeypatch.setattr('plumetools.radialstudy.run._run',lambda *a,**k:seen.append(k))
    _native_run(['dsmcFoam+','-case',str(tmp_path)],tmp_path,tmp_path/'log')
    assert seen[0]['timeout'] is None


def test_fresh_and_resume_keep_cloud_and_use_argv(tmp_path,monkeypatch):
    from plumetools.radialstudy.run import run_case
    from plumetools.foamio.checkpoints import END_DIVIDER
    p,r = prepared(tmp_path)
    monkeypatch.setattr('plumetools.radialstudy.run.require_environment',lambda: {})
    monkeypatch.setattr('plumetools.radialstudy.run.mesh_case',lambda *a,**k: {'mesh_sha256':'mesh'})
    calls = []
    def runner(argv,cwd,log):
        calls.append(argv)
        Path(log).write_text('command log')
    run_case(p,runner=runner)
    assert [a[0] for a in calls] == ['dsmcInitialise+','dsmcFoam+']
    assert str(p) in calls[-1]
    with pytest.raises(ValueError,match='fresh'):
        run_case(p,runner=runner)
    t = r['sampling']['cadence_s']
    checkpoint = p/f'{t:.17g}'
    (checkpoint/'uniform').mkdir(parents=True)
    (checkpoint/'uniform/time').write_bytes(END_DIVIDER)
    cloud = checkpoint/'lagrangian/dsmc/positions'
    cloud.parent.mkdir(parents=True)
    cloud.write_bytes(END_DIVIDER)
    calls.clear()
    run_case(p,resume=True,runner=runner)
    assert [a[0] for a in calls] == ['dsmcFoam+']
    assert cloud.read_bytes() == END_DIVIDER
    bad = p/f'{2*t:.17g}'
    bad.mkdir()
    calls.clear()
    with pytest.raises(ValueError,match='incomplete|cannot be continued'):
        run_case(p,resume=True,runner=runner)
    assert not calls


def test_command_failure_has_log_and_preserves_failure_state(tmp_path,monkeypatch):
    from plumetools.radialstudy.run import run_case
    p,r = prepared(tmp_path)
    monkeypatch.setattr('plumetools.radialstudy.run.require_environment',lambda: {})
    monkeypatch.setattr('plumetools.radialstudy.run.mesh_case',lambda *a,**k: {'mesh_sha256':'mesh'})
    def fail(argv,cwd,log):
        Path(log).write_text('initializer failed visibly')
        raise RuntimeError('initializer failed visibly')
    with pytest.raises(RuntimeError,match='failed'):
        run_case(p,runner=fail)
    assert (p/'log.dsmcInitialise').read_text() == 'initializer failed visibly'
    assert json.loads((p/'run-state.json').read_text())['segments'][-1]['status'] == 'failed'
