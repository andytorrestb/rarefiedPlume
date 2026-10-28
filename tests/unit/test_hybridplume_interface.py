"""hybridplume: mesh geometry, breakdown connectivity and interface extraction."""

from __future__ import annotations

import numpy as np
import pytest

from _hybridmesh import box_mesh
from plumetools.hybridplume import breakdown as bd
from plumetools.hybridplume import interface as itf
from plumetools.hybridplume.config import from_document
from plumetools.hybridplume.meshgeom import read_mesh
from plumetools.hybridplume.physics import BIRD_ARGON, SourceState

R = 0.1


@pytest.fixture(scope="module")
def mesh(tmp_path_factory):
    case = tmp_path_factory.mktemp("box")
    h = 0.025
    box_mesh(case, np.arange(0, 0.2 + 1e-9, h), np.arange(-0.15, 0.15 + 1e-9, h),
             np.arange(-0.15, 0.15 + 1e-9, h), R)
    return read_mesh(case)


def test_geometry_is_exact_for_boxes(mesh):
    assert np.allclose(mesh.cell_volumes, 0.025 ** 3, rtol=1e-12)
    assert mesh.cell_volumes.sum() == pytest.approx(0.2 * 0.3 * 0.3, rel=1e-12)
    # outward boundary area vectors of a closed box sum to zero
    total = sum(mesh.face_areas[mesh.patch_slice(p)].sum(axis=0) for p in mesh.patches)
    assert np.allclose(total, 0.0, atol=1e-14)
    src = mesh.face_areas[mesh.patch_slice("source")]
    assert np.allclose(src / np.linalg.norm(src, axis=1)[:, None], [-1, 0, 0])


def test_source_staircase_area_error_is_reported_not_hidden():
    from plumetools.hybridplume import mesh as hm
    doc = {"schema": "hybridplume/1", "base": {
        "source": {"diameter_m": 0.2, "speed_ratio": 2.0, "temperature_K": 300.0,
                   "knudsen": 0.01},
        "gas": {"mass_kg": 6.63e-26, "diameter_m": 4.17e-10, "t_ref_K": 273.0, "omega": 0.81},
        "domain": {"x_max_over_D": 1.0, "half_width_over_D": 1.0},
        "mesh": {"core_cell_over_D": 0.1}, "dsmc": {"backend": "dsmcFoam"}}}
    errs = []
    for h in (0.1, 0.05, 0.025):
        cfg = from_document(doc, overrides={"mesh": {"core_cell_over_D": h}})
        src = SourceState(BIRD_ARGON)
        pred = hm.predicted_source(hm.plan(cfg, src), src)
        errs.append(abs(pred.area_error))
        assert pred.exact_area_m2 == pytest.approx(np.pi * 0.01)
    assert errs[2] < errs[0]          # refinement reduces the staircase error


def _uniform_state(n_cells, rho=1e-4, t=300.0, u=(700.0, 0.0, 0.0)):
    return bd.CfdState(np.full(n_cells, rho), np.full(n_cells, t),
                       np.tile(u, (n_cells, 1)), np.full(n_cells, rho * BIRD_ARGON.gas_constant * t),
                       "1")


def test_uniform_field_is_continuum_everywhere_and_connected(mesh):
    st = _uniform_state(mesh.n_cells)
    ind = bd.indicators(st, bd.gradients(mesh, st), BIRD_ARGON)
    assert np.nanmax(ind["Kn_GLL"]) < 1e-12
    mask = bd.continuum_mask(ind, {"Kn_GLL": 0.05, "B_exp": 0.05})
    reg = bd.source_connected(mesh, mask, "source")
    assert reg["retained"].all()


def test_invalid_cells_are_never_continuum(mesh):
    st = _uniform_state(mesh.n_cells)
    rho = st.rho.copy()
    rho[5] = -1.0
    rho[7] = np.nan
    st = bd.CfdState(rho, st.T, st.U, st.p, "1")
    ind = bd.indicators(st, bd.gradients(mesh, st), BIRD_ARGON)
    mask = bd.continuum_mask(ind, {"Kn_GLL": 0.05, "B_exp": 0.05})
    assert not mask[5] and not mask[7]
    assert np.isnan(ind["Kn_GLL"][5])


def test_disconnected_pocket_is_handed_to_dsmc(mesh):
    c = mesh.cell_centres
    near = np.linalg.norm(c, axis=1) < 0.08
    pocket = np.linalg.norm(c - [0.17, 0.0, 0.0], axis=1) < 0.03
    reg = bd.source_connected(mesh, near | pocket, "source")
    assert reg["retained"][near].all() and not reg["retained"][pocket].any()
    assert reg["n_disconnected_continuum_components"] == 1


def test_interface_closes_and_points_cfd_to_dsmc(mesh):
    c = mesh.cell_centres
    retained = itf.connected_to_source(mesh, np.linalg.norm(c, axis=1) < 0.09)
    surf = itf.extract(mesh, retained)
    cl = itf.closure(mesh, surf)
    assert cl["relative_to_boundary_area"] < 1e-12
    # the normal points from the CFD cell to the DSMC cell
    d = c[surf.dsmc_cell] - c[surf.cfd_cell]
    assert np.all(np.einsum("ij,ij->i", d, surf.normal) > 0)
    assert retained[surf.cfd_cell].all() and not retained[surf.dsmc_cell].any()
    # a ball of radius 0.09 does not cover the r = 0.1 disk fully
    assert surf.source_uncovered.size > 0 and surf.source_covered.size > 0
    assert surf.outer_contact.size == 0


def test_region_touching_the_far_field_is_detected(mesh):
    retained = itf.connected_to_source(mesh, np.abs(mesh.cell_centres[:, 1]) < 0.06)
    surf = itf.extract(mesh, retained)
    assert surf.outer_contact.size > 0


def test_linear_face_interpolation_is_exact_for_linear_fields(mesh):
    c = mesh.cell_centres
    retained = itf.connected_to_source(mesh, np.linalg.norm(c, axis=1) < 0.09)
    surf = itf.extract(mesh, retained)
    phi = 3.0 + 2.0 * c[:, 0] - 5.0 * c[:, 1] + 0.5 * c[:, 2]
    fv = itf.face_values(mesh, surf, phi)
    fc = surf.centre
    assert np.allclose(fv, 3.0 + 2.0 * fc[:, 0] - 5.0 * fc[:, 1] + 0.5 * fc[:, 2], rtol=1e-12)


def _cfg(mode="automatic", **criteria):
    doc = {"schema": "hybridplume/1", "base": {
        "source": {"diameter_m": 0.2, "speed_ratio": 2.0, "temperature_K": 300.0,
                   "knudsen": 0.01},
        "gas": {"mass_kg": 6.63e-26, "diameter_m": 4.17e-10, "t_ref_K": 273.0, "omega": 0.81},
        "domain": {"x_max_over_D": 1.0, "half_width_over_D": 0.75},
        "mesh": {"core_cell_over_D": 0.125}, "dsmc": {"backend": "dsmcFoam"},
        "breakdown": {"thresholds": {"Kn_GLL": 0.05, "B_exp": 0.05}},
        "interface": {"mode": mode, "criteria": criteria}}}
    return from_document(doc)


def test_supersonic_outflow_surface_passes_and_reversed_flow_fails(mesh):
    src = SourceState(BIRD_ARGON, knudsen=0.01)
    c = mesh.cell_centres
    retained = itf.connected_to_source(mesh, c[:, 0] < 0.1)   # a slab: outflow face +x
    surf = itf.extract(mesh, retained)
    cfg = _cfg(max_outer_boundary_contact_fraction=1.0, min_overlap_layers=0)
    st = _uniform_state(mesh.n_cells, u=(3000.0, 0.0, 0.0))
    g = bd.gradients(mesh, st)
    ind = bd.indicators(st, g, BIRD_ARGON)
    doc = itf.assess(mesh, surf, st, g, ind, cfg, src, mode="automatic")
    assert doc["criteria"]["return_mass_fraction"]["pass"]
    assert doc["normal_mach"]["subsonic_area_fraction"] == 0.0
    back = _uniform_state(mesh.n_cells, u=(-3000.0, 0.0, 0.0))
    g2 = bd.gradients(mesh, back)
    doc2 = itf.assess(mesh, surf, back, g2, bd.indicators(back, g2, BIRD_ARGON), cfg, src,
                      mode="automatic")
    assert doc2["status"] == itf.STATUS_UNSUITABLE
    assert any("return_mass_fraction" in r for r in doc2["reasons"])


def test_empty_region_is_an_explicit_outcome(mesh):
    src = SourceState(BIRD_ARGON, knudsen=0.01)
    surf = itf.extract(mesh, np.zeros(mesh.n_cells, bool))
    st = _uniform_state(mesh.n_cells)
    g = bd.gradients(mesh, st)
    doc = itf.assess(mesh, surf, st, g, bd.indicators(st, g, BIRD_ARGON), _cfg(), src,
                     mode="automatic")
    assert doc["status"] == itf.STATUS_NONE
