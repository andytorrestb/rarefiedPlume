from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from _refinement import tree_hash

LOCAL = Path(__file__).resolve().parents[2]/'cases/cai2012-radial-resolution/study-local.yaml'


@pytest.mark.parametrize('kn',[100,10,1,.1,.01,.001,.0001])
def test_local_initial_mesh_keeps_100_cells_and_meets_predicted_edges(kn):
    from plumetools.radialstudy.config import load_study
    from plumetools.radialstudy.local import initial_local_model
    from plumetools.radialstudy.mesh import initial_mesh
    from plumetools.radialstudy.physics import source_for_kn
    cfg = load_study(LOCAL)
    m = initial_local_model(cfg,kn)
    mesh = initial_mesh(m)
    assert mesh.n_cells==100 and mesh.transverse==1
    assert m.r_out-m.r_in<=m.source.diameter_m*(1+1e-12)
    assert 0<m.config.domain.angle_deg<=2.5
    assert m.r_in==source_for_kn(cfg,kn).r_in
    ratio = mesh.geometry.max_edge_m/m.state(mesh.geometry.centres_m[:,0])['lambda']
    assert np.max(ratio)<=.8


def test_local_generation_is_read_only_in_preview_and_small(tmp_path):
    from plumetools.radialstudy.config import load_study
    from plumetools.radialstudy.generation import generate_initial
    from plumetools.radialstudy.case import load_case
    cfg = load_study(LOCAL)
    before = tree_hash(tmp_path)
    rows = generate_initial(cfg,tmp_path/'Cases-local',dry_run=True)
    assert tree_hash(tmp_path)==before
    assert len(rows)==7
    assert all(r['numerics']['n_cells']==100 for r in rows)
    assert all(r['numerics']['predicted_total_parcels']<=5000 for r in rows)
    rows = generate_initial(cfg,tmp_path/'Cases-local')
    for row in rows:
        m,mesh,record = load_case(row['path'])
        assert mesh.n_cells==100 and m.config.local.enabled
        assert record['numerics']['predicted_min_parcels']==pytest.approx(25)
        assert 'local_study_identity' in record
    assert rows[-1]['record']['reference_breakdown']['status']=='CROSSING_IN_DOMAIN'
    assert rows[-1]['record']['reference_breakdown']['r_over_D']>10


def test_measured_local_contraction_shortens_sampling(tmp_path):
    from plumetools.radialstudy.config import load_study
    from plumetools.radialstudy.local import initial_local_model,contract_local_model
    from plumetools.radialstudy.mesh import initial_mesh
    from plumetools.radialstudy.case import sizing
    m = initial_local_model(load_study(LOCAL),.01)
    mesh = initial_mesh(m)
    child,newmesh = contract_local_model(m,mesh,np.min(mesh.geometry.max_edge_m)*.1)
    assert newmesh.n_cells==100 and child.r_in==m.r_in
    assert child.r_out<m.r_out and newmesh.angle_rad<mesh.angle_rad
    assert np.max(newmesh.geometry.max_edge_m)<=np.min(mesh.geometry.max_edge_m)*.1
    assert child.state(child.r_in)==m.state(m.r_in)
    assert sizing(child,newmesh)[1]['end_s']<sizing(m,mesh)[1]['end_s']


def parent_plan(tmp_path,monkeypatch,action='needs_mesh'):
    from plumetools.radialstudy.config import load_study
    from plumetools.radialstudy.generation import generate_initial
    from plumetools.radialstudy.case import load_case
    from plumetools.refinement.geometry import file_sha256
    from plumetools.refinement.model import RefinementPlan
    row = generate_initial(load_study(LOCAL),tmp_path/'Cases-local',kn_values=[.01])[0]
    parent = Path(row['path'])
    m,mesh,record = load_case(parent)
    monkeypatch.setattr('plumetools.refinement.plan.mesh_fingerprint',lambda path:mesh.geometry.mesh_sha256)
    plan = RefinementPlan(action,mesh.geometry.mesh_sha256,'measured',mesh.geometry,
        mesh.geometry.max_edge_m*.2,mesh.density_estimate*2,
        np.arange(100) if action=='needs_mesh' else np.array([],dtype=int),20,
        record['numerics']['particle_weight']/2,{'source_files':{'case.json':file_sha256(parent/'case.json')}})
    return parent,m,mesh,record,plan


def test_local_revision_contracts_not_refines_and_keeps_common_identity(tmp_path,monkeypatch):
    from plumetools.radialstudy.generation import revise_case
    from plumetools.radialstudy.case import load_case,study_identity
    p,m,mesh,old,plan = parent_plan(tmp_path,monkeypatch)
    before = tree_hash(p)
    preview = revise_case(p,plan,dry_run=True)
    assert tree_hash(p)==before and preview['numerics']['n_cells']==100
    assert preview['sampling']['end_s']<old['sampling']['end_s']
    row = revise_case(p,plan)
    child,newmesh,record = load_case(row['path'])
    assert study_identity(record)==study_identity(old)
    assert record['physical_identity']!=old['physical_identity']
    assert record['inlet']==old['inlet'] and newmesh.n_cells==100
    assert child.r_out<m.r_out and record['sampling']['end_s']<old['sampling']['end_s']
    assert 'domain' in row['corrections']
    assert tree_hash(p)==before


def test_local_population_retry_keeps_geometry_and_sampling(tmp_path,monkeypatch):
    from plumetools.radialstudy.generation import revise_case
    p,m,mesh,old,plan = parent_plan(tmp_path,monkeypatch,'needs_parcels')
    new = revise_case(p,plan)['record']
    assert new['physical_identity']==old['physical_identity']
    assert new['sampling']==old['sampling'] and new['numerics']['n_cells']==100
    assert new['numerics']['particle_weight']==plan.recommended_weight


def test_local_parent_rejects_a_different_kn(tmp_path,monkeypatch):
    from plumetools.radialstudy.local import initial_local_model
    from plumetools.radialstudy.mesh import initial_mesh
    from plumetools.radialstudy.case import write_case
    p,m,mesh,_,_ = parent_plan(tmp_path,monkeypatch)
    other = initial_local_model(m.config,.001)
    with pytest.raises(ValueError,match='local revision changes'):
        write_case(p.parent/'r1_ppc20',other,initial_mesh(other),parent=p)


def test_dense_local_patch_uses_available_transverse_width():
    from plumetools.radialstudy.config import load_study
    from plumetools.radialstudy.local import initial_local_model
    from plumetools.radialstudy.mesh import initial_mesh
    m = initial_local_model(load_study(LOCAL),.0001)
    mesh = initial_mesh(m)
    width = 2*m.r_out*np.tan(mesh.angle_rad)
    assert width>=.9*.8*float(m.state(m.r_in)['lambda'])


def test_local_controller_retries_smaller_patch_and_restart_groups_once(tmp_path,monkeypatch):
    from test_radialstudy_adaptive import configured
    from plumetools.radialstudy.run import run_study
    from plumetools.radialstudy.case import load_case
    p,old,_ = configured(tmp_path,monkeypatch,local=True)
    result = run_study(p.parent.parent)
    assert result['completed']==1
    _,mesh,new = load_case(result['cases'][0]['path'])
    assert mesh.n_cells==100
    assert new['sampling']['end_s']<old['sampling']['end_s']
    before = tree_hash(tmp_path)
    result = run_study(p.parent.parent,resume=True)
    assert result['skipped']==1 and len(result['cases'])==1
    assert tree_hash(tmp_path)==before


def test_local_measured_population_budget_stops_before_continuing(tmp_path,monkeypatch):
    from test_radialstudy_adaptive import configured
    from plumetools.radialstudy.adaptive import run_adaptive
    import json
    p,record,_ = configured(tmp_path,monkeypatch,local=True,first_parcels=80,max_parcels=5000)
    assert record['numerics']['predicted_total_parcels']<5000
    with pytest.raises(RuntimeError,match='measured parcel budget'):
        run_adaptive(p)
    state = json.loads((p.parent/'controller-state.json').read_text())
    assert state['outcome']=='budget_limited' and 'certificate' not in state
    assert (p/'0/C').exists() and not (p.parent/'r1_ppc20').exists()
