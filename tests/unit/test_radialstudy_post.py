from dataclasses import replace
import json

import numpy as np
import pytest

from _refinement import write_field
from test_radialstudy_run import prepared


def synthetic(case,record,mesh,monkeypatch,*,density_scale=1,parcels=40,drift=False,velocity=500.):
    from plumetools.radialstudy.case import load_case
    m,_,_ = load_case(case)
    g = mesh.geometry
    write_field(case/'0/C',g.centres_m,[0,1,0,0,0,0,0])
    write_field(case/'0/V',g.volumes_m3,[0,3,0,0,0,0,0])
    monkeypatch.setattr('plumetools.radialstudy.run.mesh_fingerprint',lambda p:g.mesh_sha256)
    monkeypatch.setattr('plumetools.radialstudy.post.read_geometry',lambda *a,**k:g)
    (case/'run-state.json').write_text(json.dumps({'case_id':record['case_id'],'mesh_sha256':g.mesh_sha256,
        'segments':[{'start_s':0,'status':'completed'}]}))
    for i in range(12):
        time = record['sampling']['start_s']+(i+1)*record['sampling']['cadence_s']
        folder = case/f'{time:.17g}'
        folder.mkdir(exist_ok=True)
        n = mesh.density_estimate*density_scale*(1+i*.1 if drift else 1)
        write_field(folder/'rhoN_Ar',n,[0,-3,0,0,0,0,0])
        write_field(folder/'Ttra_Ar',np.full(mesh.n_cells,300.),[0,0,0,1,0,0,0])
        write_field(folder/'U_Ar',np.tile([velocity,0,0],(mesh.n_cells,1)),[0,1,-1,0,0,0,0])
        write_field(folder/'dsmcNMean_Ar',np.full(mesh.n_cells,parcels),[0,0,0,0,0,0,0])


def test_measured_mfp_controls_acceptance_no_extra_sweep(tmp_path,monkeypatch):
    from plumetools.radialstudy.post import assess_case,case_plan
    from plumetools.radialstudy.mesh import read_radial_mesh
    p,r = prepared(tmp_path)
    mesh = read_radial_mesh(p)
    synthetic(p,r,mesh,monkeypatch,parcels=20)
    result = assess_case(p)
    assert result['status'] == 'criteria_satisfied'
    assert result['assessment']['action'] == 'criteria_satisfied'
    assert result['comparison_required'] is False
    assert result['profile']['lambda_predicted_m'] != []
    plan,_ = case_plan(p)
    assert plan.action == 'criteria_satisfied'
    synthetic(p,r,mesh,monkeypatch,density_scale=1e7)
    # Recollection is explicit after replacing test measurements.
    (p/'sampling-windows.json').unlink()
    new = assess_case(p)
    assert new['assessment']['action'] == 'needs_mesh'
    assert new['recommended_delta_t_s'] < r['numerics']['delta_t_s']
    assert new['profile']['lambda_measured_m'][0] < result['profile']['lambda_measured_m'][0]


def test_timestep_only_failure_survives_passing_spatial_assessment(tmp_path,monkeypatch):
    from plumetools.radialstudy.post import assess_case,case_plan
    from plumetools.radialstudy.mesh import read_radial_mesh
    p,r = prepared(tmp_path)
    synthetic(p,r,read_radial_mesh(p),monkeypatch,velocity=5000.)
    result = assess_case(p)
    assert result['assessment']['action'] == 'criteria_satisfied'
    assert result['status'] == 'needs_timestep'
    plan,dt = case_plan(p)
    assert plan.action == 'criteria_satisfied'
    assert dt < r['numerics']['delta_t_s']


def test_missing_underoccupied_and_drifting_data(tmp_path,monkeypatch):
    from plumetools.radialstudy.post import assess_case,case_plan
    from plumetools.radialstudy.mesh import read_radial_mesh
    p,r = prepared(tmp_path)
    mesh = read_radial_mesh(p)
    assert assess_case(p)['status'] == 'missing_measurements'
    synthetic(p,r,mesh,monkeypatch,parcels=10)
    assert assess_case(p)['assessment']['action'] == 'needs_parcels'
    plan,_ = case_plan(p)
    assert plan.recommended_weight == pytest.approx(r['numerics']['particle_weight']/2)
    synthetic(p,r,mesh,monkeypatch,drift=True)
    (p/'sampling-windows.json').unlink()
    assert assess_case(p)['sampling']['stationary'] is False
    assert assess_case(p)['status'] == 'unresolved'


def test_projection_conserves_moments_and_recovers_thermal_energy():
    from plumetools.radialstudy.post import project_moments,primitive
    from plumetools.cai2012.gas import ARGON
    # Equal and opposite velocities yield additional thermal energy after averaging.
    rho = np.ones(2)
    p = np.array([[2.,0,0],[-2.,0,0]])
    energy = np.full(2, 1.5*ARGON.gas_constant_j_per_kg_k*300+2)
    raw = {'rho':rho,'momentum':p,'energy':energy}
    out,coverage = project_moments(np.array([1.,2.,3.]),raw,np.array([1.,3.]))
    assert out['rho'][0] == pytest.approx(1)
    # Frustum weighting is proportional to r^3 differences.
    assert out['momentum'][0,0] == pytest.approx((2*7-2*19)/26)
    assert coverage[0] == 1
    states = primitive(out,ARGON)
    assert states['T'][0] > 300
    assert out['energy'][0]*26 == pytest.approx(np.sum(energy*np.array([7,19])))


def test_optional_measured_breakdown_needs_gradient_coverage():
    from plumetools.radialstudy.post import measured_breakdown
    result = measured_breakdown(np.array([1.,2.]),np.array([1.,np.nan]),np.ones(2),
        np.ones(2),np.ones(2),np.array([True,False]))
    assert result['status'] == 'UNAVAILABLE'


def test_effective_target_survives_generation_assessment_and_planning(tmp_path,monkeypatch):
    from plumetools.radialstudy.case import load_case,write_case
    from plumetools.radialstudy.post import assess_case,case_plan
    original,r = prepared(tmp_path)
    m,mesh,_ = load_case(original)
    child = tmp_path/'ppc40'
    record = write_case(child,m,mesh,target_parcels=40,particle_weight=r['numerics']['particle_weight']/2)
    synthetic(child,record,mesh,monkeypatch,parcels=30)
    assert load_case(child)[0].config.targets.min_mean_parcels == 40
    assert assess_case(child)['assessment']['action'] == 'needs_parcels'
    plan,_ = case_plan(child)
    assert plan.target_parcels == 40
    assert plan.recommended_weight == pytest.approx(record['numerics']['particle_weight']*.75)
