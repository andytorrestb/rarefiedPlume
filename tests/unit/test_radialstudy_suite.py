"""Serial study commands reuse each case's initialization/restart lifecycle."""
from dataclasses import replace
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from _refinement import tree_hash
from test_radialstudy_physics import model

REPO = Path(__file__).resolve().parents[2]


def case(root,kn,revision=0,*,over_budget=False):
    from plumetools.radialstudy.case import write_case
    from plumetools.radialstudy.physics import source_for_kn
    from plumetools.radialstudy.mesh import initial_mesh
    m = model(kn)
    cfg = replace(m.config,domain=replace(m.config.domain,outer_over_D=.82,angle_deg=.5),
                  mesh=replace(m.config.mesh,geometric_cap_over_D=.04))
    if over_budget:
        cfg = replace(cfg,time=replace(cfg.time,max_steps=1))
    m = source_for_kn(cfg,kn)
    name = f'{kn:g}'.replace('.','p')
    path = root/f'Kn{name}'/f'r{revision}_ppc20'
    record = write_case(path,m,initial_mesh(m))
    return path,record


def started(path,record,*,complete=False):
    from plumetools.foamio.checkpoints import END_DIVIDER
    time = record['sampling']['end_s'] if complete else record['sampling']['cadence_s']
    checkpoint = path/f'{time:.17g}'
    (checkpoint/'uniform').mkdir(parents=True)
    (checkpoint/'uniform/time').write_bytes(END_DIVIDER)
    cloud = checkpoint/'lagrangian/dsmc/positions'
    cloud.parent.mkdir(parents=True)
    cloud.write_bytes(END_DIVIDER)
    (path/'run-state.json').write_text(json.dumps({'case_id':record['case_id'],
        'mesh_sha256':'mesh','sampling':record['sampling'],'segments':[]}))
    return checkpoint


def test_all_cases_run_sequentially_in_kn_and_revision_order(tmp_path,monkeypatch):
    from plumetools.radialstudy.run import run_study
    root = tmp_path/'Cases with spaces'
    lower,_ = case(root,10)
    later,_ = case(root,100,revision=2)
    first,_ = case(root,100)
    # An in-progress generation must not become part of the queue.
    shutil.copytree(first,first.with_name('r0_ppc20.stage-incomplete'))
    calls = []
    def run(path,**kwargs):
        calls.append((path,kwargs))
        return {'status':'completed'}
    monkeypatch.setattr('plumetools.radialstudy.run.run_case',run)
    result = run_study(root)
    assert [p for p,_ in calls]==[first,later,lower]
    assert all(k['resume'] is False for _,k in calls)
    assert result['completed']==3


def test_restart_skips_finished_resumes_interrupted_and_starts_untouched(tmp_path,monkeypatch):
    from plumetools.radialstudy.run import run_study
    root = tmp_path/'Cases'
    done,r = case(root,100,over_budget=True)
    started(done,r,complete=True)
    interrupted,r = case(root,10)
    started(interrupted,r)
    fresh,_ = case(root,1)
    monkeypatch.setattr('plumetools.radialstudy.run.mesh_fingerprint',lambda _: 'mesh')
    calls = []
    def run(path,**kwargs):
        calls.append((path,kwargs['resume']))
        return {'status':'completed'}
    monkeypatch.setattr('plumetools.radialstudy.run.run_case',run)
    before = tree_hash(done)
    result = run_study(root,resume=True)
    assert calls==[(interrupted,True),(fresh,False)]
    assert result['completed']==2 and result['skipped']==1
    assert tree_hash(done)==before


def test_fresh_mode_refuses_started_cases_before_running_any(tmp_path,monkeypatch):
    from plumetools.radialstudy.run import run_study
    root = tmp_path/'Cases'
    case(root,100)
    p,r = case(root,10)
    started(p,r)
    monkeypatch.setattr('plumetools.radialstudy.run.run_case',lambda *a,**k:pytest.fail('started a case before preflight finished'))
    with pytest.raises(ValueError,match='restart|resume'):
        run_study(root)


def test_bad_checkpoint_blocks_restart_without_touching_other_cases(tmp_path,monkeypatch):
    from plumetools.radialstudy.run import run_study
    root = tmp_path/'Cases'
    case(root,100)
    p,r = case(root,10)
    checkpoint = started(p,r)
    (checkpoint/'uniform/time').unlink()
    monkeypatch.setattr('plumetools.radialstudy.run.mesh_fingerprint',lambda _: 'mesh')
    monkeypatch.setattr('plumetools.radialstudy.run.run_case',lambda *a,**k:pytest.fail('ran a case despite an incomplete checkpoint'))
    before = tree_hash(root)
    with pytest.raises(ValueError,match='incomplete|cannot be continued'):
        run_study(root,resume=True)
    assert tree_hash(root)==before


def test_budget_gaps_are_reported_in_dry_run_and_checked_before_execution(tmp_path,monkeypatch):
    from plumetools.radialstudy.run import run_study
    root = tmp_path/'Cases'
    case(root,100)
    case(root,10,over_budget=True)
    calls = []
    monkeypatch.setattr('plumetools.radialstudy.run.run_case',lambda p,**k:calls.append(k) or {'status':'completed'})
    before = tree_hash(root)
    result = run_study(root,dry_run=True)
    assert result['dry_run'] and len(result['cases'])==2
    assert result['cases'][1]['budget_gaps']
    assert not calls and tree_hash(root)==before
    with pytest.raises(ValueError,match='budget'):
        run_study(root)
    assert not calls
    run_study(root,allow_budget_exceedance=True)
    assert len(calls)==2 and all(k['allow_budget_exceedance'] for k in calls)


def test_failure_stops_the_serial_queue(tmp_path,monkeypatch):
    from plumetools.radialstudy.run import run_study
    root = tmp_path/'Cases'
    case(root,100)
    case(root,10)
    calls = []
    def failed(path,**kwargs):
        calls.append(path)
        raise RuntimeError('native solver failed; see log.dsmcFoam')
    monkeypatch.setattr('plumetools.radialstudy.run.run_case',failed)
    with pytest.raises(RuntimeError,match='Kn100.*native solver failed'):
        run_study(root)
    assert len(calls)==1


def test_empty_case_tree_explains_generation(tmp_path):
    from plumetools.radialstudy.run import run_study
    with pytest.raises(ValueError,match='[Gg]enerate|generate_cases'):
        run_study(tmp_path/'missing')


@pytest.mark.parametrize('script',['Allrun','Allrestart'])
def test_study_shell_scripts_dry_run_from_another_directory(tmp_path,script):
    if not shutil.which('bash') or os.name!='posix':
        pytest.skip('requires a POSIX shell')
    root = tmp_path/'Cases with spaces'
    case(root,100)
    before = tree_hash(root)
    env = dict(os.environ,PYTHON=sys.executable,PYTHONPATH=str(REPO),PYTHONDONTWRITEBYTECODE='1')
    path = REPO/'cases/cai2012-radial-resolution'/script
    result = subprocess.run(['bash',str(path),'--cases',str(root),'--dry-run'],
                            cwd=tmp_path,env=env,capture_output=True,text=True,timeout=30)
    assert result.returncode==0,result.stderr
    report = json.loads(result.stdout)
    assert report['dry_run'] and len(report['cases'])==1
    assert report['mode']==('restart' if script=='Allrestart' else 'fresh')
    assert tree_hash(root)==before
    assert b'\r' not in path.read_bytes()
