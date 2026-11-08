from dataclasses import replace
from pathlib import Path
import shutil

import numpy as np
import pytest

from _refinement import tree_hash
from test_radialstudy_config import STUDY
from test_radialstudy_physics import model


def test_generate_seven_only_and_dry_run_writes_nothing(tmp_path):
    from plumetools.radialstudy.config import load_study
    from plumetools.radialstudy.generation import generate_initial
    cfg = load_study(STUDY)
    dest = tmp_path/'Cases'
    assert len(generate_initial(cfg,dest,dry_run=True)) == 7
    assert not dest.exists()
    rows = generate_initial(cfg,dest)
    assert len(rows) == 7
    assert all(Path(row['path']).name == 'r0_ppc20' for row in rows)
    assert len(list(dest.glob('Kn*/r*/case.json'))) == 7
    before = tree_hash(dest)
    generate_initial(cfg,dest)
    assert tree_hash(dest) == before
    moved = tmp_path/'moved'
    shutil.copytree(dest,moved)
    from plumetools.radialstudy.case import load_case
    m, mesh, r = load_case(next(moved.glob('Kn*/r0_ppc20')))
    assert r['physical_identity'] == mesh.identity
    assert mesh.identity == m.identity(angle_rad=mesh.angle_rad)


def test_independent_numerical_revisions_and_noop(tmp_path, monkeypatch):
    from plumetools.radialstudy.case import write_case
    from plumetools.radialstudy.mesh import initial_mesh
    from plumetools.radialstudy.generation import revise_case
    from plumetools.refinement.model import RefinementPlan
    # Source verification is separately exercised with actual meshes/field hashes.
    monkeypatch.setattr('plumetools.radialstudy.generation.verify_sources', lambda *a: None)
    m = model(100)
    mesh = initial_mesh(m)
    parent = tmp_path/'r0_ppc20'
    old = write_case(parent,m,mesh)
    g = mesh.geometry
    target = g.max_edge_m*2
    target[0] *= .2
    p = RefinementPlan('needs_mesh',g.mesh_sha256,'fields',g,target,
        mesh.density_estimate,np.array([0]),20,old['numerics']['particle_weight']/2,{})
    before = tree_hash(parent)
    row = revise_case(parent,p,correction='mesh')
    assert row['record']['numerics']['n_cells'] > mesh.n_cells
    assert row['record']['physical_identity'] == old['physical_identity']
    assert row['record']['numerics']['target_parcels'] == 20
    parcel = revise_case(parent,replace(p,action='needs_parcels',selected_cells=np.array([],dtype=int)),correction='parcels')
    assert parcel['record']['numerics']['mesh_id'] == old['numerics']['mesh_id']
    assert parcel['record']['numerics']['particle_weight'] == p.recommended_weight
    step = revise_case(parent,p,correction='timestep',delta_t_s=old['numerics']['delta_t_s']/2)
    assert step['record']['numerics']['particle_weight'] == old['numerics']['particle_weight']
    assert step['record']['numerics']['mesh_id'] == old['numerics']['mesh_id']
    assert tree_hash(parent) == before
    satisfied = replace(p,action='criteria_satisfied',selected_cells=np.array([],dtype=int))
    assert revise_case(parent,satisfied)['action'] == 'no_revision_needed'
    assert len(list(tmp_path.glob('r*/case.json'))) == 4
