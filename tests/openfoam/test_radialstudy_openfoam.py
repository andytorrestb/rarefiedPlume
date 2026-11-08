"""Built native mesh geometry, later extended with a bounded solver mechanics run."""
from dataclasses import replace
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'unit'))
from test_radialstudy_physics import model

pytestmark = pytest.mark.needs_openfoam


def test_actual_radial_mesh_geometry(tmp_path):
    if not shutil.which('blockMesh'):
        pytest.skip('requires OpenFOAM v2512')
    from plumetools.radialstudy.mesh import initial_mesh, materialize
    from plumetools.refinement.geometry import read_geometry
    m = model(100)
    m = replace(m, r_out=.164)
    mesh = initial_mesh(m)
    # Deliberately exercise transverse subdivisions too.
    mesh = replace(mesh, transverse=3)
    destination = tmp_path/'ancestor with spaces'/'mesh'
    materialize(mesh, destination)
    actual = read_geometry(destination, time='0')
    np.testing.assert_allclose(actual.centres_m, mesh.geometry.centres_m, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(actual.volumes_m3, mesh.geometry.volumes_m3, rtol=1e-9)
    np.testing.assert_allclose(actual.max_edge_m, mesh.geometry.max_edge_m, rtol=1e-9)


def test_bounded_plus_mechanics_run(tmp_path):
    if not shutil.which('dsmcFoam+'):
        pytest.skip('requires the v2512 hyStrath plus port')
    from plumetools.radialstudy.case import write_case
    from plumetools.radialstudy.physics import source_for_kn
    from plumetools.radialstudy.mesh import initial_mesh
    from plumetools.radialstudy.run import mesh_case,run_case,collect_case
    m = model(100)
    cfg = replace(m.config,domain=replace(m.config.domain,outer_over_D=.82,angle_deg=.5),
                  mesh=replace(m.config.mesh,geometric_cap_over_D=.04))
    m = source_for_kn(cfg,100)
    mesh = initial_mesh(m)
    record = write_case(tmp_path/'mechanics only',m,mesh)
    assert mesh.n_cells <= 100
    assert record['numerics']['planned_steps'] <= 2000
    assert record['numerics']['predicted_total_parcels'] <= 5000
    case = tmp_path/'mechanics only'
    assert mesh_case(case)['n_cells'] == mesh.n_cells
    result = run_case(case)
    assert result['status'] == 'completed'
    windows = collect_case(case)
    assert len(windows) == 12
    for w in windows:
        assert all((case/w['time']/name).is_file() for name in ('rhoN_Ar','Ttra_Ar','U_Ar','dsmcNMean_Ar','mfp_Ar','mct_Ar'))
        assert (case/w['time']/'lagrangian/dsmc/positions').is_file()
    assert run_case(case,resume=True)['status'] == 'already_finished'
    from plumetools.radialstudy.post import case_plan
    from plumetools.refinement.plan import verify_sources
    plan,_ = case_plan(case)
    for relative in ('0/C','0/V','run-state.json','sampling-windows.json'):
        assert relative in plan.provenance['source_files']
        path = case/relative
        original = path.read_bytes()
        path.write_bytes(original+b'\n')
        with pytest.raises(ValueError,match='changed|stale'):
            verify_sources(case,plan)
        path.write_bytes(original)
    verify_sources(case,plan)


def test_serial_scripts_run_and_restart_native_cases(tmp_path):
    if not shutil.which('dsmcFoam+'):
        pytest.skip('requires the v2512 hyStrath plus port')
    from _refinement import tree_hash
    from test_radialstudy_suite import case
    from plumetools.foamio.times import time_directories
    root = tmp_path/'Cases with spaces'
    done,_ = case(root,100)
    interrupted,_ = case(root,10)
    repo = Path(__file__).resolve().parents[2]
    env = dict(os.environ,PYTHON=sys.executable,PYTHONPATH=str(repo),PYTHONDONTWRITEBYTECODE='1')

    def script(name):
        result = subprocess.run(['bash',str(repo/'cases/cai2012-radial-resolution'/name),
                                 '--cases',str(root)],cwd=tmp_path,env=env,
                                capture_output=True,text=True,timeout=120)
        assert result.returncode==0,result.stderr
        return json.loads(result.stdout)

    report = script('Allrun')
    assert report['completed']==2 and report['skipped']==0
    before = tree_hash(done)
    # Restore the second case's first saved checkpoint, as from an interrupted run.
    times = time_directories(interrupted)
    assert len(times)>1
    for _,name in times[1:]:
        shutil.rmtree(interrupted/name)
    initializer = (interrupted/'log.dsmcInitialise').read_bytes()
    initial_fields = tree_hash(interrupted/'0')
    fresh,_ = case(root,1)
    report = script('Allrestart')
    assert report['completed']==2 and report['skipped']==1
    assert [row['action'] for row in report['cases']]==['skip','resume','fresh']
    assert tree_hash(done)==before
    assert (interrupted/'log.dsmcInitialise').read_bytes()==initializer
    assert tree_hash(interrupted/'0')==initial_fields
    state = json.loads((interrupted/'run-state.json').read_text())
    assert state['segments'][-1]['start_s']==times[0][0]
    assert state['segments'][-1]['status']=='completed'
    assert time_directories(fresh)


def test_native_write_boundary_segments(tmp_path):
    if not shutil.which('dsmcFoam+'):
        pytest.skip('requires the v2512 hyStrath plus port')
    from test_radialstudy_suite import case
    from plumetools.radialstudy.segments import run_segment
    from plumetools.radialstudy.run import verify_inputs
    path,record = case(tmp_path/'Cases',100)
    canonical = (path/'system/controlDict').read_bytes()
    cadence = record['sampling']['cadence_s']
    first = run_segment(path,cadence)
    initializer = (path/'log.dsmcInitialise').read_bytes()
    second = run_segment(path,2*cadence)
    assert first['end_s']==pytest.approx(second['start_s'],rel=1e-12)
    assert (path/'log.dsmcInitialise').read_bytes()==initializer
    assert (path/'system/controlDict').read_bytes()==canonical
    verify_inputs(path,record)


def test_native_adaptive_clean_retry_and_persistent_summary(tmp_path):
    if not shutil.which('dsmcFoam+'):
        pytest.skip('requires the v2512 hyStrath plus port')
    from plumetools.radialstudy.case import write_case,sizing
    from plumetools.radialstudy.physics import source_for_kn
    from plumetools.radialstudy.mesh import initial_mesh
    from plumetools.radialstudy.run import run_study
    from _refinement import tree_hash
    m = model(100)
    cfg = replace(m.config,domain=replace(m.config.domain,outer_over_D=.82,angle_deg=.5),
                  mesh=replace(m.config.mesh,initial_cells=3),
                  controller=replace(m.config.controller,enabled=True,stability_relative=.5))
    m = source_for_kn(cfg,100)
    mesh = initial_mesh(m)
    num,_ = sizing(m,mesh)
    root = tmp_path/'Cases-adaptive'
    parent = root/'Kn100/r0_ppc20'
    write_case(parent,m,mesh,particle_weight=num['particle_weight']*4,delta_t_s=num['delta_t_s']/2)
    result = run_study(root)
    assert result['completed']==1
    state = json.loads((parent.parent/'controller-state.json').read_text())
    assert state['outcome']=='accepted' and len(state['attempts'])>=2
    assert not (parent/'0/lagrangian').exists() and not (parent/'constant/polyMesh').exists()
    log = (root/'adaptive-study.log').read_text()
    assert 'retry' in log and 'cleanup' in log and 'accepted' in log
    before = tree_hash(root)
    assert run_study(root,resume=True)['skipped']==1
    assert tree_hash(root)==before


def test_native_parallel_adaptive_cases_and_restart(tmp_path):
    if not shutil.which('dsmcFoam+'):
        pytest.skip('requires the v2512 hyStrath plus port')
    from plumetools.radialstudy.case import write_case,sizing
    from plumetools.radialstudy.physics import source_for_kn
    from plumetools.radialstudy.mesh import initial_mesh
    from plumetools.radialstudy.run import run_study
    from _refinement import tree_hash
    root = tmp_path/'parallel cases'
    for kn in (100,10):
        m = model(kn)
        cfg = replace(m.config,domain=replace(m.config.domain,outer_over_D=.82,angle_deg=.5),
                      mesh=replace(m.config.mesh,initial_cells=3),
                      controller=replace(m.config.controller,enabled=True,stability_relative=.5))
        m = source_for_kn(cfg,kn)
        mesh = initial_mesh(m)
        num,_ = sizing(m,mesh)
        write_case(root/f'Kn{kn}/r0_ppc20',m,mesh,
                   particle_weight=num['particle_weight']/2,delta_t_s=num['delta_t_s']/2)
    result = run_study(root,jobs=2)
    assert result['completed']==2 and result['jobs']==2
    lines = (root/'adaptive-study.log').read_text().splitlines()
    assert sum(' accepted ' in line for line in lines)==2
    assert all(re.match(r'^\d{4}-\d\d-\d\dT',line) for line in lines)
    before = tree_hash(root)
    assert run_study(root,resume=True,jobs=2)['skipped']==2
    assert tree_hash(root)==before
