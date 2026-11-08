import json

import pytest

from test_radialstudy_run import prepared


def test_six_blocks_two_disjoint_windows_and_hashes(tmp_path,monkeypatch):
    from plumetools.radialstudy.run import collect_case,read_windows
    p,r = prepared(tmp_path)
    monkeypatch.setattr('plumetools.radialstudy.run.mesh_fingerprint',lambda p:'mesh')
    (p/'run-state.json').write_text(json.dumps({'case_id':r['case_id'],'mesh_sha256':'mesh',
        'segments':[{'start_s':0,'status':'completed'}]}))
    cadence = r['sampling']['cadence_s']
    for i in range(1,17):
        time = i*cadence
        folder = p/f'{time:.17g}'
        folder.mkdir()
        for name in ('rhoN_Ar','Ttra_Ar','UMean_Ar','dsmcNMean_Ar'):
            (folder/name).write_text(f'data-{i}')
    windows = collect_case(p)
    assert len(windows) == 12
    assert {w['block'] for w in windows} == set(range(6))
    assert all(w['steps'] == r['numerics']['steps_per_window'] for w in windows)
    assert all(a['end_s'] == pytest.approx(b['start_s']) for a,b in zip(windows,windows[1:]))
    assert all(w['basis'] == 'disjoint_native_windows' for w in windows)
    assert read_windows(p) == windows
    # A previous Allpost may have inventoried only the first half of a run.
    (p/'sampling-windows.json').write_text(json.dumps(windows[:6]))
    assert len(read_windows(p)) == 12
    first = windows[0]
    (p/next(iter(first['source_files']))).write_text('changed')
    with pytest.raises(ValueError,match='changed'):
        read_windows(p)


def test_overlapping_windows_are_not_independent(tmp_path):
    from plumetools.radialstudy.run import validate_windows
    with pytest.raises(ValueError,match='overlap|disjoint'):
        validate_windows([{'start_s':0,'end_s':1,'basis':'disjoint_native_windows'},
                          {'start_s':0,'end_s':2,'basis':'disjoint_native_windows'}])
    with pytest.raises(ValueError,match='disjoint'):
        validate_windows([{'start_s':0,'end_s':1,'basis':'cumulative'}])
