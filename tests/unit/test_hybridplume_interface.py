"""hybridplume: mesh geometry and breakdown connectivity, on a Python-written hex box."""

from __future__ import annotations

import numpy as np
import pytest

from _hybridmesh import box_mesh
from plumetools.hybridplume import breakdown as bd
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
