"""hybridplume: mesh geometry, on a Python-written 3-D hex box."""

from __future__ import annotations

import numpy as np
import pytest

from _hybridmesh import box_mesh
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
