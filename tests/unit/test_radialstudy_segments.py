import hashlib
import json
import re

import pytest

from test_radialstudy_run import prepared


def setup_native(monkeypatch):
    monkeypatch.setattr('plumetools.radialstudy.segments.require_environment',lambda: {})
    monkeypatch.setattr('plumetools.radialstudy.segments.mesh_case',lambda *a,**k:{'mesh_sha256':'mesh'})


def test_segments_stop_at_checkpoints_restore_inputs_and_initialize_once(tmp_path,monkeypatch):
    from plumetools.radialstudy.segments import run_segment
    from plumetools.foamio.checkpoints import END_DIVIDER
    p,r = prepared(tmp_path)
    setup_native(monkeypatch)
    canonical = (p/'system/controlDict').read_bytes()
    calls = []
    def runner(argv,cwd,log):
        calls.append(argv[0])
        if argv[0]=='dsmcFoam+':
            end = float(re.search(r'endTime\s+([^;]+)',(p/'system/controlDict').read_text())[1])
            folder = p/f'{end:.17g}'
            (folder/'uniform').mkdir(parents=True)
            (folder/'uniform/time').write_bytes(END_DIVIDER)
            (folder/'lagrangian/dsmc').mkdir(parents=True)
            (folder/'lagrangian/dsmc/positions').write_bytes(END_DIVIDER)
    dt = r['sampling']['cadence_s']
    run_segment(p,dt,runner=runner)
    result = run_segment(p,2*dt,runner=runner)
    assert calls==['dsmcInitialise+','dsmcFoam+','dsmcFoam+']
    assert result['start_s']==pytest.approx(dt)
    assert (p/'system/controlDict').read_bytes()==canonical
    assert not (p/'controller-control.json').exists()
    assert json.loads((p/'run-state.json').read_text())['segments'][-1]['status']=='completed'


def test_failed_segment_restores_control_and_preserves_failure(tmp_path,monkeypatch):
    from plumetools.radialstudy.segments import run_segment
    p,r = prepared(tmp_path)
    setup_native(monkeypatch)
    canonical = (p/'system/controlDict').read_bytes()
    def failed(*args):
        raise RuntimeError('native failure')
    with pytest.raises(RuntimeError,match='native failure'):
        run_segment(p,r['sampling']['cadence_s'],runner=failed)
    assert (p/'system/controlDict').read_bytes()==canonical
    assert json.loads((p/'run-state.json').read_text())['segments'][-1]['status']=='failed'


def test_override_recovery_is_exact_and_refuses_unrelated_edits(tmp_path):
    from plumetools.radialstudy.segments import recover_control
    from plumetools.refinement.geometry import file_sha256
    p,r = prepared(tmp_path)
    control = p/'system/controlDict'
    original = control.read_bytes()
    active = re.sub(rb'endTime\s+[^;]+;',b'endTime 1e-8;',original)
    manifest = {'canonical':original.decode(),'canonical_sha256':file_sha256(control),
                'active_sha256':hashlib.sha256(active).hexdigest()}
    state = p/'controller-control.json'
    state.write_text(json.dumps(manifest))
    control.write_bytes(active+b'\n// unrelated edit')
    with pytest.raises(ValueError,match='changed|override'):
        recover_control(p)
    assert state.exists()
    control.write_bytes(active)
    recover_control(p)
    assert control.read_bytes()==original and not state.exists()


def test_successful_exit_without_target_checkpoint_is_not_completion(tmp_path,monkeypatch):
    from plumetools.radialstudy.segments import run_segment
    p,r = prepared(tmp_path)
    setup_native(monkeypatch)
    with pytest.raises(ValueError,match='checkpoint'):
        run_segment(p,r['sampling']['cadence_s'],runner=lambda *a:None)
