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


def test_containing_cell_mapping_refuses_invalid_and_outside_points(tmp_path):
    from _hybridmesh import box_mesh
    from plumetools.hybridplume.meshgeom import read_mesh
    from plumetools.hybridplume.planar import sample_containing_cells
    from plumetools.hybridplume.workflow import WorkflowError
    box_mesh(tmp_path, [0, .1, .2], [-.1, 0, .1], [-.005, .005], .1)
    mesh = read_mesh(tmp_path)
    allowed = np.ones(4, bool)
    points = np.array([[.03, -.03, 0], [.16, .02, 0]])
    assert sample_containing_cells(mesh, points, allowed).tolist() == [0, 3]
    allowed[3] = False
    with pytest.raises(WorkflowError, match="accepted|allowed"):
        sample_containing_cells(mesh, points, allowed)
    with pytest.raises(WorkflowError, match="outside"):
        sample_containing_cells(mesh, [[.3, 0, 0]], np.ones(4, bool))


def test_empty_depth_planes_add_no_reconstructed_flux(tmp_path):
    from types import SimpleNamespace
    from plumetools.hybridplume.budgets import reconstructed_patch_fluxes
    from plumetools.hybridplume.physics import BIRD_ARGON
    mesh = SimpleNamespace(patches={"front": SimpleNamespace(n_faces=1, type="empty")})
    assert reconstructed_patch_fluxes(tmp_path, "1", mesh, BIRD_ARGON) == {}


def test_planar_timestep_uses_actual_inplane_width():
    from plumetools.hybridplume import config, resources
    cfg = config.load(STUDY, "planar-laptop-kn0p01", {"dsmc": {"collisions": False}})
    source = config.build_source(cfg)
    est = resources.dsmc_estimate(
        cfg, source, number_density=[1e15], temperature=[300], speed=[10],
        volumes=[1e-6], sizes=[[.1, .1]], tracking_sizes=[.001],
        length_scale_m=.4, basis="synthetic planar cell")
    expected = .2*.001/(10 + 3*np.sqrt(source.gas.gas_constant*300))
    assert est["dt_s"] == pytest.approx(expected)
    assert est["densest_cell"]["cell_size_m"] == pytest.approx(.1)


def test_planar_workflow_reports_infeasible_fit_before_dsmc(tmp_path):
    from _hybridmesh import box_mesh
    from plumetools.hybridplume import workflow
    case = tmp_path/"planar"
    setup = workflow.setup(case, STUDY, "planar-laptop-kn0p01")
    assert setup["n_ranks"] == 1
    box_mesh(case/"cfd", [0, .1, .2], [-.1, 0, .1], [-.005, .005], .1)
    np.savez(case/"results/breakdown.npz", retained=np.zeros(4, bool), time="0.001")
    result = workflow.interface(case)
    assert result["status"] == "no_continuum_region"
    with pytest.raises(workflow.WorkflowError, match="no feasible ellipse"):
        workflow.dsmc_setup(case)


def test_mapping_tolerates_mesh_writer_coordinate_roundoff(tmp_path):
    from _hybridmesh import box_mesh
    from plumetools.hybridplume.meshgeom import read_mesh
    from plumetools.hybridplume.planar import sample_containing_cells
    box_mesh(tmp_path, [0, .1, .2], [-.1, 0, .1], [-.005, .005], .1)
    mesh = read_mesh(tmp_path)
    mesh.points[(mesh.points[:, 1] == 0) & (mesh.points[:, 2] < 0), 1] += 1e-16
    points = [[.03, -.03, 0], [.16, .02, 0]]
    assert sample_containing_cells(mesh, points, np.ones(4, bool)).tolist() == [0, 3]
