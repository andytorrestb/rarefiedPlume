from dataclasses import replace
import json
from pathlib import Path
import re

import numpy as np
import pytest

from _refinement import tree_hash,write_field
from test_radialstudy_physics import model


def configured(tmp_path,monkeypatch,*,first_parcels=40,max_retries=5,max_steps=10000000):
    from plumetools.radialstudy.case import write_case,load_case
    from plumetools.radialstudy.mesh import initial_mesh,read_radial_mesh
    from plumetools.radialstudy.physics import source_for_kn
    from plumetools.foamio.checkpoints import END_DIVIDER
    import plumetools.radialstudy.segments as segments
    import plumetools.radialstudy.adaptive as adaptive
    m = model(100)
    cfg = replace(m.config,domain=replace(m.config.domain,outer_over_D=.82,angle_deg=.5),
                  mesh=replace(m.config.mesh,initial_cells=3),
                  controller=replace(m.config.controller,enabled=True,max_retries=max_retries),
                  time=replace(m.config.time,max_steps=max_steps))
    m = source_for_kn(cfg,100)
    p = tmp_path/'Cases'/'Kn100'/'r0_ppc20'
    r = write_case(p,m,initial_mesh(m))
    geometry = lambda path:read_radial_mesh(path).geometry
    def mesh(path,**kwargs):
        g = geometry(path)
        write_field(path/'0/C',g.centres_m,[0,1,0,0,0,0,0])
        write_field(path/'0/V',g.volumes_m3,[0,3,0,0,0,0,0])
        return {'mesh_sha256':g.mesh_sha256}
    monkeypatch.setattr(segments,'require_environment',lambda:{})
    monkeypatch.setattr(segments,'mesh_case',mesh)
    for target in ('plumetools.radialstudy.run.mesh_fingerprint',
                   'plumetools.refinement.plan.mesh_fingerprint',
                   'plumetools.radialstudy.adaptive.mesh_fingerprint'):
        monkeypatch.setattr(target,lambda path:geometry(path).mesh_sha256)
    monkeypatch.setattr('plumetools.radialstudy.post.read_geometry',lambda path,**kwargs:geometry(path))
    calls = []
    def runner(argv,path,log):
        calls.append((path.name,argv[0]))
        Path(log).write_text('native log\n')
        if argv[0]!='dsmcFoam+':
            return
        m,mesh,r = load_case(path)
        end = float(re.search(r'endTime\s+([^;]+)',(path/'system/controlDict').read_text())[1])
        folder = path/f'{end:.17g}'
        (folder/'uniform').mkdir(parents=True,exist_ok=True)
        (folder/'uniform/time').write_bytes(END_DIVIDER)
        (folder/'lagrangian/dsmc').mkdir(parents=True,exist_ok=True)
        (folder/'lagrangian/dsmc/positions').write_bytes(END_DIVIDER)
        values = {'rhoN_Ar':(mesh.density_estimate,[0,-3,0,0,0,0,0]),
                  'Ttra_Ar':(np.full(mesh.n_cells,300.),[0,0,0,1,0,0,0]),
                  'U_Ar':(np.tile([500.,0,0],(mesh.n_cells,1)),[0,1,-1,0,0,0,0]),
                  'dsmcNMean_Ar':(np.full(mesh.n_cells,first_parcels if path.name=='r0_ppc20' else 40),[0]*7)}
        for name,(value,dims) in values.items():
            write_field(folder/name,value,dims)
    monkeypatch.setattr(adaptive,'run_segment',lambda path,end:segments.run_segment(path,end,runner=runner))
    return p,r,calls


def test_accept_early_and_restart_is_read_only(tmp_path,monkeypatch):
    from plumetools.radialstudy.adaptive import run_adaptive,adaptive_preview
    p,r,calls = configured(tmp_path,monkeypatch)
    before = tree_hash(tmp_path)
    assert adaptive_preview(p)['action']=='fresh'
    assert tree_hash(tmp_path)==before
    result = run_adaptive(p)
    assert result['status']=='accepted'
    state = json.loads((p.parent/'controller-state.json').read_text())
    assert state['accepted_time_s']<r['sampling']['end_s']
    assert len([c for c in calls if c[1]=='dsmcInitialise+'])==1
    log = p.parent.parent/'adaptive-study.log'
    assert 'accepted' in log.read_text() and 'min_mean_parcels=40' in log.read_text()
    before = tree_hash(tmp_path)
    assert run_adaptive(p,resume=True)['status']=='already_finished'
    assert tree_hash(tmp_path)==before
    (p/'system/controlDict').write_text('edited')
    with pytest.raises(ValueError,match='changed'):
        adaptive_preview(p,resume=True)


def test_retry_prepares_child_cleans_parent_and_retains_summary(tmp_path,monkeypatch):
    from plumetools.radialstudy.adaptive import run_adaptive
    from plumetools.radialstudy.case import load_case
    p,r,calls = configured(tmp_path,monkeypatch,first_parcels=10)
    note = p/'my-notes.txt'
    note.write_text('keep this')
    result = run_adaptive(p)
    child = p.parent/'r1_ppc20'
    assert result['status']=='accepted' and child.exists()
    assert note.read_text()=='keep this' and not list(p.glob('0.0*'))
    assert not (p/'0/C').exists() and (p/'case.json').exists()
    assert load_case(child)[2]['numerics']['particle_weight']==pytest.approx(r['numerics']['particle_weight']*10/25)
    assert [name for name,cmd in calls if cmd=='dsmcInitialise+']==['r0_ppc20','r1_ppc20']
    log = (p.parent.parent/'adaptive-study.log').read_text()
    assert 'retry' in log and 'cleanup' in log and 'accepted' in log
    assert 'needs_parcels' in log and 'old_settings' in log


def test_retry_and_cumulative_step_limits_preserve_evidence(tmp_path,monkeypatch):
    from plumetools.radialstudy.adaptive import run_adaptive
    p,r,_ = configured(tmp_path,monkeypatch,first_parcels=10,max_retries=0)
    with pytest.raises(RuntimeError,match='retry limit'):
        run_adaptive(p,allow_budget_exceedance=True)
    assert (p/'0/C').exists() and not (p.parent/'r1_ppc20').exists()
    assert 'budget_limited' in (p.parent.parent/'adaptive-study.log').read_text()


def test_step_budget_is_checked_before_launch(tmp_path,monkeypatch):
    from plumetools.radialstudy.adaptive import run_adaptive
    p,r,calls = configured(tmp_path,monkeypatch,max_steps=1)
    with pytest.raises(RuntimeError,match='step budget'):
        run_adaptive(p)
    assert not calls


def test_unresolved_data_never_passes_and_is_not_cleaned(tmp_path,monkeypatch):
    from plumetools.radialstudy.adaptive import run_adaptive
    p,_,_ = configured(tmp_path,monkeypatch)
    monkeypatch.setattr('plumetools.radialstudy.adaptive.checkpoint_assessment',
                        lambda path:({'status':'unresolved','reason':'not stable'},None))
    with pytest.raises(RuntimeError,match='sampling limit'):
        run_adaptive(p)
    assert (p/'0/C').exists() and not (p.parent/'r1_ppc20').exists()


def test_cleanup_refuses_unowned_and_linked_outputs(tmp_path,monkeypatch):
    from plumetools.radialstudy.adaptive import clean_attempt
    p,_,_ = configured(tmp_path,monkeypatch)
    with pytest.raises(ValueError,match='owned|owner'):
        clean_attempt(p,p.parent,'unknown')


def test_restart_finishes_prepared_retry_without_duplicate_attempt(tmp_path,monkeypatch):
    import plumetools.radialstudy.adaptive as adaptive
    p,_,_ = configured(tmp_path,monkeypatch,first_parcels=10)
    clean = adaptive.clean_attempt
    monkeypatch.setattr(adaptive,'clean_attempt',lambda *a:(_ for _ in ()).throw(RuntimeError('interrupted cleanup')))
    with pytest.raises(RuntimeError,match='interrupted cleanup'):
        adaptive.run_adaptive(p)
    prefix = (p.parent.parent/'adaptive-study.log').read_text()
    assert (p.parent/'r1_ppc20/case.json').exists()
    monkeypatch.setattr(adaptive,'clean_attempt',clean)
    assert adaptive.run_adaptive(p,resume=True)['status']=='accepted'
    assert not (p.parent/'r2_ppc20').exists()
    assert (p.parent.parent/'adaptive-study.log').read_text().startswith(prefix)


def test_cleanup_detects_outside_link_before_deleting_any_outputs(tmp_path,monkeypatch):
    import plumetools.radialstudy.adaptive as adaptive
    p,_,_ = configured(tmp_path,monkeypatch,first_parcels=10)
    outside = tmp_path/'outside'
    outside.mkdir()
    (outside/'precious').write_text('preserve')
    clean = adaptive.clean_attempt
    def linked(path,root,owner):
        try:
            (path/'fieldMeasurements').symlink_to(outside,target_is_directory=True)
        except OSError:
            pytest.skip('creating symlinks requires operating-system permission')
        clean(path,root,owner)
    monkeypatch.setattr(adaptive,'clean_attempt',linked)
    with pytest.raises(ValueError,match='linked|outside'):
        adaptive.run_adaptive(p)
    assert (outside/'precious').read_text()=='preserve'
    assert (p/'0/C').exists()


def test_cumulative_budget_does_not_reset_at_retry(tmp_path,monkeypatch):
    from plumetools.radialstudy.adaptive import run_adaptive
    _,r,_ = configured(tmp_path/'reference',monkeypatch)
    budget = r['numerics']['steps_per_window']*12
    p,_,_ = configured(tmp_path/'bounded',monkeypatch,first_parcels=10,max_steps=budget)
    with pytest.raises(RuntimeError,match='step budget'):
        run_adaptive(p)
    state = json.loads((p.parent/'controller-state.json').read_text())
    assert state['current']=='r1_ppc20' and state['total_steps']==budget


def test_post_reports_sizing_basis_and_analysis_excludes_cleaned_attempt(tmp_path,monkeypatch):
    from plumetools.radialstudy.adaptive import run_adaptive
    from plumetools.radialstudy.post import post_case,analyse_study
    p,_,_ = configured(tmp_path,monkeypatch,first_parcels=10)
    result = run_adaptive(p)
    post = post_case(result['path'],plots=False)
    assert post['status']=='criteria_satisfied' and post['sampling']['basis']=='adaptive_sizing'
    analysis = analyse_study(p.parent.parent)
    assert len(analysis['cases'])==1 and analysis['cases'][0]['path']==result['path']


def test_interrupted_child_audit_publication_can_restart(tmp_path,monkeypatch):
    import plumetools.radialstudy.adaptive as adaptive
    import plumetools.radialstudy.generation as generation
    p,_,_ = configured(tmp_path,monkeypatch,first_parcels=10)
    write = generation.write_plan
    monkeypatch.setattr(generation,'write_plan',lambda *a:(_ for _ in ()).throw(RuntimeError('publication interrupted')))
    with pytest.raises(RuntimeError,match='publication interrupted'):
        adaptive.run_adaptive(p)
    monkeypatch.setattr(generation,'write_plan',write)
    assert adaptive.run_adaptive(p,resume=True)['status']=='accepted'
    assert (p.parent/'r1_ppc20/revision.json').exists()
    assert not (p.parent/'r2_ppc20').exists()


def test_acceptance_decision_is_not_lost_at_last_sampling_checkpoint(tmp_path,monkeypatch):
    import plumetools.radialstudy.adaptive as adaptive
    from plumetools.foamio.times import time_directories
    p,r,_ = configured(tmp_path,monkeypatch)
    assess = adaptive.checkpoint_assessment
    def final_only(path):
        return (assess(path) if time_directories(path)[-1][0]>=r['sampling']['end_s']*.999999 else
                ({'status':'unresolved','reason':'wait for final fixture checkpoint'},None))
    monkeypatch.setattr(adaptive,'checkpoint_assessment',final_only)
    save = adaptive._save
    interrupted = []
    def interrupt(path,data):
        save(path,data)
        if (Path(path).name=='controller-state.json' and not interrupted and
                data['last_checked_s']>=r['sampling']['end_s']*.999999 and data['outcome']=='running'
                and not data.get('pending_retry')):
            interrupted.append(True)
            raise RuntimeError('decision interrupted')
    monkeypatch.setattr(adaptive,'_save',interrupt)
    try:
        result = adaptive.run_adaptive(p)
    except RuntimeError:
        result = adaptive.run_adaptive(p,resume=True)
    assert result['status']=='accepted'


def test_cleanup_preserves_unknown_files_inside_output_directories(tmp_path,monkeypatch):
    import plumetools.radialstudy.adaptive as adaptive
    p,_,_ = configured(tmp_path,monkeypatch,first_parcels=10)
    clean = adaptive.clean_attempt
    note = p/'postProcessing/my-unrelated-analysis.txt'
    def with_note(path,root,owner):
        note.parent.mkdir(exist_ok=True)
        note.write_text('my analysis')
        clean(path,root,owner)
    monkeypatch.setattr(adaptive,'clean_attempt',with_note)
    assert adaptive.run_adaptive(p)['status']=='accepted'
    assert note.read_text()=='my analysis'


@pytest.mark.parametrize('legacy',[False,True])
def test_concurrent_controller_refuses_before_any_mutation(tmp_path,monkeypatch,legacy):
    from plumetools.radialstudy.adaptive import run_adaptive
    from plumetools.radialstudy.locking import exclusive_lock
    p,_,calls = configured(tmp_path,monkeypatch)
    lockpath = p.parent.parent/'.adaptive-study.lock' if legacy else p.parent/'.adaptive-controller.lock'
    with exclusive_lock(lockpath):
        before = tree_hash(p)
        with pytest.raises(RuntimeError,match='already running|locked'):
            run_adaptive(p,resume=True)
        assert not calls and tree_hash(p)==before
        assert not (p.parent/'controller-state.json').exists()
        assert not (p.parent.parent/'adaptive-study.log').exists()
