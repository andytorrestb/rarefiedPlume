from dataclasses import replace
import json
import sys

import numpy as np
import pytest

from _refinement import tree_hash
from test_radialstudy_physics import model


def test_case_weight_native_dicts_and_idempotence(tmp_path):
    from plumetools.radialstudy.case import write_case
    from plumetools.radialstudy.mesh import initial_mesh
    m = model(100)
    mesh = initial_mesh(m)
    p = tmp_path/'r0_ppc20'
    r = write_case(p,m,mesh)
    assert r['numerics']['particle_weight'] == pytest.approx(np.min(mesh.density_estimate*mesh.geometry.volumes_m3)/20)
    assert r['numerics']['predicted_min_parcels'] == pytest.approx(20)
    assert r['inlet']['n'] == pytest.approx(float(m.state(m.r_in)['n']))
    assert r['inlet']['T'] == pytest.approx(299.845, abs=.001)
    assert 'dsmcFoam+' in (p/'system/controlDict').read_text()
    assert f"seedNumber {r['seed']};" in (p/'constant/dsmcProperties').read_text()
    assert r['seed_policy']['base_seed'] == 1729
    assert 'dsmcFreeStreamInflowPatch' in (p/'system/boundariesDict').read_text()
    assert 'dsmcSpecularWallPatch' in (p/'system/boundariesDict').read_text()
    assert 'dsmcSpherical' not in (p/'constant/dsmcProperties').read_text()
    assert 'resetAtOutputUntilTime' not in (p/'system/fieldPropertiesDict').read_text()
    assert 'resetAtOutput on;' in (p/'system/fieldPropertiesDict').read_text()
    (p/'0.01').mkdir()
    (p/'0.01/sentinel').write_text('do not overwrite results')
    before = tree_hash(p)
    assert write_case(p,m,mesh)['case_id'] == r['case_id']
    assert tree_hash(p) == before
    with pytest.raises(ValueError, match='conflict'):
        write_case(p,m,mesh,target_parcels=40)


def test_explicit_weight_timestep_and_physical_hashes(tmp_path):
    from plumetools.radialstudy.case import write_case, load_case
    from plumetools.radialstudy.mesh import initial_mesh
    m = model(.001)
    mesh = initial_mesh(m)
    a = write_case(tmp_path/'a',m,mesh)
    b = write_case(tmp_path/'b',m,mesh,particle_weight=a['numerics']['particle_weight']/2,
                   delta_t_s=a['numerics']['delta_t_s']/2,target_parcels=40)
    assert a['physical_identity'] == b['physical_identity']
    assert a['sampling'] == b['sampling']
    loaded, readmesh, record = load_case(tmp_path/'b')
    assert loaded.physical_hash == m.physical_hash
    assert record['numerics']['particle_weight'] == b['numerics']['particle_weight']
    assert readmesh.identity == mesh.identity


@pytest.mark.skipif(sys.platform != 'win32',reason='exercises actual Windows newline translation')
def test_windows_lf_scripts(tmp_path):
    from plumetools.radialstudy.case import write_case
    from plumetools.radialstudy.mesh import initial_mesh
    m = model(100)
    path = tmp_path/'r0'
    write_case(path,m,initial_mesh(m))
    for name in ('Allmesh','Allrun','Allpost'):
        assert b'\r' not in (path/name).read_bytes()


def test_loading_freezes_recorded_source_inlet_domain_and_sampling(tmp_path,monkeypatch):
    from dataclasses import replace
    from plumetools.radialstudy.case import write_case,load_case
    from plumetools.radialstudy.mesh import initial_mesh
    from plumetools.radialstudy.physics import source_for_kn
    m = model(100)
    mesh = initial_mesh(m)
    parent = tmp_path/'parent'
    old = write_case(parent,m,mesh)
    def other_platform(config,kn):
        calculated = source_for_kn(config,kn)
        # One representable step is the observed Windows/WSL anchor difference.
        return replace(calculated,r_in=np.nextafter(calculated.r_in,np.inf),
            physical_hash='different-derived-hash',
            source=replace(calculated.source,r_a=np.nextafter(calculated.source.r_a,np.inf),
                           T_t=np.nextafter(calculated.source.T_t,np.inf)))
    monkeypatch.setattr('plumetools.radialstudy.case.source_for_kn',other_platform)
    loaded,frozen,_ = load_case(parent)
    assert loaded.r_in == old['physical_identity']['r_in_m']
    assert loaded.physical_hash == old['physical_identity']['source_domain_hash']
    assert loaded.source.T_t == m.source.T_t
    assert loaded.source.r_a == m.source.r_a
    assert {k:float(v) for k,v in loaded.state(loaded.r_in).items()} == old['inlet']
    child = write_case(tmp_path/'child',loaded,frozen,target_parcels=40,parent=parent)
    assert child['sampling'] == old['sampling']
    assert child['physical_identity'] == old['physical_identity']
    assert child['inlet'] == old['inlet']


def test_changed_case_record_cannot_supply_frozen_physics(tmp_path):
    from plumetools.radialstudy.case import write_case,load_case
    from plumetools.radialstudy.mesh import initial_mesh
    m = model(100)
    case = tmp_path/'case'
    record = write_case(case,m,initial_mesh(m))
    record['inlet']['T'] = 123.
    (case/'case.json').write_text(json.dumps(record))
    with pytest.raises(ValueError,match='case.*identity|case.*record'):
        load_case(case)


@pytest.mark.parametrize('winerror',[5,32])
def test_transient_windows_publish_lock_keeps_case_atomic(tmp_path,monkeypatch,winerror):
    from pathlib import Path
    from plumetools.radialstudy.case import write_case,load_case
    from plumetools.radialstudy.mesh import initial_mesh
    original = Path.rename
    attempts = []
    def sharing_lock(stage,destination):
        attempts.append(stage)
        if len(attempts)<3:
            error = PermissionError('transient Windows directory lock')
            error.winerror = winerror
            raise error
        return original(stage,destination)
    monkeypatch.setattr(Path,'rename',sharing_lock)
    m = model(100)
    record = write_case(tmp_path/'case',m,initial_mesh(m))
    assert len(attempts)==3
    assert load_case(tmp_path/'case')[2]['case_id']==record['case_id']
    assert not list(tmp_path.glob('*.stage-*'))


def test_persistent_windows_publish_lock_is_bounded_and_leaves_no_partial_case(tmp_path,monkeypatch):
    from pathlib import Path
    import time
    from plumetools.radialstudy.case import write_case
    from plumetools.radialstudy.mesh import initial_mesh
    attempts = []
    def locked(stage,destination):
        attempts.append(stage)
        error = PermissionError('persistent Windows directory lock')
        error.winerror = 5
        raise error
    monkeypatch.setattr(Path,'rename',locked)
    monkeypatch.setattr(time,'sleep',lambda _:None)
    m = model(100)
    with pytest.raises(PermissionError,match='persistent'):
        write_case(tmp_path/'case',m,initial_mesh(m))
    assert 1<len(attempts)<=10
    assert not (tmp_path/'case').exists()
    assert not list(tmp_path.glob('*.stage-*'))
def test_native_measurement_interval_is_physical(tmp_path):
    import re
    from plumetools.radialstudy.case import write_case
    from plumetools.radialstudy.mesh import initial_mesh
    from test_radialstudy_physics import model
    m = model(100)
    path = tmp_path/'timing'
    record = write_case(path,m,initial_mesh(m))
    text = (path/'system/controlDict').read_text()
    assert re.search(r'writeControl\s+(\w+)',text)[1]=='runTime'
    interval = float(re.search(r'writeInterval\s+([^;]+)',text)[1])
    assert interval==pytest.approx(record['sampling']['cadence_s'])
    assert interval/record['numerics']['delta_t_s']==pytest.approx(record['numerics']['steps_per_window'])
