"""The collisionless plume-plus-plate-gas density (plumetools.plateugf.collisionless)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from plumetools.nozzleplume.config import load_case_config
from plumetools.plateugf import collisionless

ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "cases" / "plate-unigas" / "baseCase"


@pytest.fixture(scope="module")
def model(tmp_path_factory):
    case = tmp_path_factory.mktemp("case")
    cfg = yaml.safe_load((TEMPLATE / "case.yaml").read_text(encoding="utf-8"))
    cfg["plate"]["standoff_m"] = 5.0
    cfg["exit"]["knudsen"] = 100.0
    cfg["domain"]["x_max_m"] = 7.0
    (case / "case.yaml").write_text(yaml.safe_dump(cfg), encoding="utf-8")
    # a coarse quadrature keeps the test quick; the physics is the same
    return collisionless.Model(load_case_config(case), patch_m=1.0)


def test_on_the_wall_the_plate_gas_is_the_half_range_maxwellian_density(model):
    # Over the middle of a patch, just off the wall: Gamma sqrt(pi) / c_w.
    centre_patch = model.gamma[model.gamma.shape[0] // 2, model.gamma.shape[1] // 2]
    point = np.array([[model.standoff_m - 1.0e-4, 0.5, 0.5]])
    n = model.plate_gas(point)[0]
    assert n == pytest.approx(centre_patch * model.wall_density_factor, rel=0.01)


def test_the_plate_shadows_the_plume_and_emits_nothing_behind_itself(model):
    s = model.standoff_m
    behind = np.array([[s + 1.0, 0.0, 0.0]])          # straight behind the centre
    beside = np.array([[s + 1.0, 30.0, 0.0]])         # past the plate's edge
    assert model.plume(behind)[0] == 0.0
    assert model.plume(beside)[0] > 0.0
    assert model.plate_gas(np.vstack([behind, beside])).tolist() == [0.0, 0.0]
    # in front, the plate's gas adds to the plume
    front = np.array([[s - 0.5, 0.0, 0.0]])
    assert model.plate_gas(front)[0] > 0.5 * model.plume(front)[0]


def test_the_grid_reproduces_the_exact_sum_away_from_the_plate_edges(model):
    rng = np.random.default_rng(1)
    s = model.standoff_m
    points = np.column_stack([rng.uniform(0.2, s - 0.05, 400),
                              rng.uniform(0.0, 0.8 * model.half_y_m, 400),
                              rng.uniform(0.0, 0.8 * model.half_z_m, 400)])
    exact = model.plate_gas(points, grid=False)
    tabulated = model.plate_gas(points, grid=True)
    assert np.allclose(tabulated, exact, rtol=0.03)
