import threading
import time

import pytest

from test_radialstudy_suite import case


def test_parallel_cases_overlap_and_respect_worker_limit(tmp_path,monkeypatch):
    from plumetools.radialstudy.run import run_study
    for kn in (100,10,1):
        case(tmp_path,kn)
    gate = threading.Barrier(2)
    lock = threading.Lock()
    active = peak = count = 0
    def run(path,**kwargs):
        nonlocal active,peak,count
        with lock:
            active += 1
            peak = max(peak,active)
            count += 1
            ordinal = count
        if ordinal<=2:
            gate.wait(timeout=10)
        time.sleep(.02)
        with lock:
            active -= 1
        return {'status':'completed'}
    monkeypatch.setattr('plumetools.radialstudy.run.run_case',run)
    result = run_study(tmp_path,jobs=2)
    assert peak==2 and result['completed']==3
    assert result['jobs']==2


def test_parallel_failure_does_not_launch_remaining_cases(tmp_path,monkeypatch):
    from plumetools.radialstudy.run import run_study
    for kn in (100,10,1):
        case(tmp_path,kn)
    calls = []
    def run(path,**kwargs):
        calls.append(path.parent.name)
        if path.parent.name=='Kn100':
            raise RuntimeError('native failure')
        time.sleep(.05)
        return {'status':'completed'}
    monkeypatch.setattr('plumetools.radialstudy.run.run_case',run)
    with pytest.raises(RuntimeError,match='native failure'):
        run_study(tmp_path,jobs=2)
    assert 'Kn1' not in calls


@pytest.mark.parametrize('jobs',[0,-1,True,1.5])
def test_invalid_worker_counts_are_rejected_without_running(tmp_path,jobs):
    from plumetools.radialstudy.run import run_study
    with pytest.raises(ValueError,match='jobs'):
        run_study(tmp_path,jobs=jobs)
