from dataclasses import replace

import numpy as np
import pytest

from _refinement import unit_inputs
from test_radialstudy_physics import model


def test_mesh_uses_real_edges_and_exact_frustum_geometry():
    from plumetools.radialstudy.mesh import initial_mesh
    m = model(100)
    mesh = initial_mesh(m)
    assert mesh.n_cells < 1000
    assert np.max(mesh.geometry.max_edge_m) <= .01 * (1 + 1e-10)
    expected = 4 * np.tan(mesh.angle_rad)**2 * (m.r_out**3 - m.r_in**3) / 3
    assert mesh.geometry.volumes_m3.sum() == pytest.approx(expected, rel=1e-12)
    edges = mesh.vertices[:, [1,2,3,0,5,6,7,4,4,5,6,7]] - mesh.vertices[:, [0,1,2,3,4,5,6,7,0,1,2,3]]
    np.testing.assert_allclose(mesh.geometry.max_edge_m, np.linalg.norm(edges, axis=2).max(axis=1))
    assert np.all(mesh.transport_length_m > 0)
    assert np.all(mesh.geometry.centres_m[:, 0] > mesh.boundaries[:-1].repeat(mesh.transverse**2))


def test_low_kn_mesh_has_distinct_precise_vertices():
    from plumetools.radialstudy.mesh import initial_mesh, render_block_mesh
    m = model(.0001)
    mesh = initial_mesh(m)
    assert np.all(np.diff(mesh.boundaries) > 0)
    assert mesh.boundaries[0] == m.r_in
    assert mesh.boundaries[-1] == m.r_out
    assert np.all(mesh.geometry.max_edge_m <= m.state(mesh.boundaries[:-1])['lambda'] * (1+1e-10))
    assert 'inlet' in render_block_mesh(mesh)


def test_refinement_is_local_and_keeps_fixed_geometry(tmp_path):
    from plumetools.radialstudy.mesh import initial_mesh, refine_mesh, write_radial_mesh, read_radial_mesh
    from plumetools.refinement.model import RefinementPlan
    parent = initial_mesh(model(100))
    target = parent.geometry.max_edge_m * 1.2
    target[0] *= .31
    density = np.full(parent.n_cells, 1e15)
    p = RefinementPlan('needs_mesh', parent.geometry.mesh_sha256, 'source', parent.geometry,
                       target, density, np.array([0]), 20, None, {})
    child = refine_mesh(parent, p)
    assert child.identity == parent.identity
    assert len(set(parent.boundaries[2:]) - set(child.boundaries)) == 0
    assert parent.n_cells < child.n_cells < 2 * parent.n_cells
    inherited = np.searchsorted(parent.boundaries, child.geometry.centres_m[:, 0], side='right') - 1
    assert np.all(child.geometry.max_edge_m <= target[inherited] * (1+1e-12))
    assert child.density_kind == 'parent-measured'
    write_radial_mesh(child, tmp_path)
    reread = read_radial_mesh(tmp_path)
    np.testing.assert_array_equal(reread.boundaries, child.boundaries)
    np.testing.assert_array_equal(reread.geometry.volumes_m3, child.geometry.volumes_m3)
    assert reread.identity == child.identity


def test_passing_plan_does_not_create_mesh(tmp_path):
    from plumetools.radialstudy.mesh import initial_mesh, refine_mesh, RadialRefiner
    from plumetools.refinement.model import RefinementPlan
    a = initial_mesh(model(100))
    p = RefinementPlan('criteria_satisfied', a.geometry.mesh_sha256, 'src', a.geometry,
        a.geometry.max_edge_m*2, a.density_estimate, np.array([], dtype=int), 20, None, {})
    assert refine_mesh(a, p) is a
    assert RadialRefiner().refine(tmp_path, p, tmp_path/'unused').destination is None
    assert not (tmp_path/'unused').exists()


def test_budget_never_coarsens():
    from plumetools.radialstudy.mesh import initial_mesh
    m = model(100)
    with pytest.raises(ValueError, match='budget'):
        initial_mesh(replace(m, config=replace(m.config, mesh=replace(m.config.mesh, max_cells=10))))
