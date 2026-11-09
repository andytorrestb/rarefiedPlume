"""Small numerical and setup checks for the opt-in planar hybrid prototype."""
from pathlib import Path

import numpy as np
import pytest


def ellipse_grid():
    x, y = np.meshgrid(np.arange(20) / 10, np.arange(-10, 10) / 10, indexing="ij")
    lower = np.column_stack((x.ravel(), y.ravel()))
    upper = lower + 0.1
    corners = np.maximum(np.abs(lower), np.abs(upper))
    retained = (corners[:, 0] / 1.4)**2 + (corners[:, 1] / 0.45)**2 < 1
    return lower, upper, retained


def test_fit_respects_resolved_ellipse():
    from plumetools.hybridplume.ellipse import fit_ellipse, ellipse_points
    lower, upper, retained = ellipse_grid()
    fit = fit_ellipse(lower, upper, retained, source_radius=0.5, clearance=0.01)
    assert fit["status"] == "fitted"
    assert fit["area_m2"] > 0.45
    arc = ellipse_points(fit, 301)
    assert np.all((arc[:, 0] / 1.4)**2 + (arc[:, 1] / 0.45)**2 <= 1 + 1e-7)
    assert np.allclose(arc[[0, -1], 0], 0, atol=1e-9)
    assert np.max(np.abs(arc[[0, -1], 1])) <= 0.5


def test_excluded_interior_pocket_limits_the_whole_ellipse():
    from plumetools.hybridplume.ellipse import fit_ellipse
    lower, upper, retained = ellipse_grid()
    hole = (np.isclose(lower[:, 0], 0.6)
            & (lower[:, 1] >= -0.1 - 1e-8) & (upper[:, 1] <= 0.1 + 1e-8))
    retained[hole] = False
    fit = fit_ellipse(lower, upper, retained, source_radius=0.5, clearance=0.01)
    assert fit["status"] == "fitted"
    assert fit["xc"] + fit["a"] <= 0.59 + 1e-6


def test_no_source_connected_ellipse_is_reported():
    from plumetools.hybridplume.ellipse import fit_ellipse
    lower, upper, retained = ellipse_grid()
    retained[lower[:, 0] < 0.3] = False
    fit = fit_ellipse(lower, upper, retained, source_radius=0.5, clearance=0.01)
    assert fit["status"] == "no_feasible_ellipse"
    assert "a" not in fit


ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / "cases/cai2012-hybrid/hybrid.yaml"


def test_planar_profile_generates_one_layer_and_slit_flux(tmp_path):
    from plumetools.hybridplume import config, mesh, cfd
    cfg = config.load(STUDY, "planar-laptop-kn0p01")
    source = config.build_source(cfg)
    plan = mesh.plan(cfg, source)
    assert plan.divisions == (40, 40, 1)
    assert source.area_m2 == pytest.approx(0.002)
    assert source.with_knudsen(0.1).area_m2 == pytest.approx(0.002)
    assert mesh.predicted_source(plan, source).area_m2 == pytest.approx(0.002)
    dictionary = mesh.render_block_mesh_dict(plan)
    assert dictionary.count("type empty;") == 2
    assert all(text.count("type empty;") == 2 for text in cfd.render_fields(cfg, source).values())


@pytest.mark.parametrize("overrides", [
    {"mesh": {"geometry": "cartesian"}},
    {"mesh": {"depth_over_D": 0}},
    {"interface": {"fit": {"clearance_cells": -0.1}}},
])
def test_incompatible_planar_settings_are_refused(overrides):
    from plumetools.hybridplume import config
    with pytest.raises(config.ConfigError, match="geometry|depth|clearance|requires"):
        config.load(STUDY, "planar-laptop-kn0p01", overrides)


def test_curved_mesh_keeps_partial_slit_and_native_ellipsoid():
    from plumetools.hybridplume import config
    from plumetools.hybridplume.planar_mesh import render_dsmc_mesh
    cfg = config.load(STUDY, "planar-laptop-kn0p01")
    fit = {"xc": 0.025, "a": 0.11, "b": 0.07}
    dictionary = render_dsmc_mesh(cfg, config.build_source(cfg), fit)
    assert "type searchableSphere;" in dictionary
    assert dictionary.count("type empty;") == 2
    assert "source" in dictionary and "interface" in dictionary
