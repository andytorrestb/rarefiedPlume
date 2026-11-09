"""Focused full-3D fitted-interface checks; no production case or solver."""
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / "cases/cai2012-hybrid/hybrid.yaml"


def voxels():
    x, y, z = np.meshgrid(np.arange(10)/5, np.arange(-4, 4)/5,
                          np.arange(-4, 4)/5, indexing="ij")
    lower = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    return lower, lower+0.2


def test_fit_excludes_whole_axis_and_off_axis_voxels():
    from plumetools.hybridplume.ellipsoid import fit_ellipsoid, verify_ellipsoid
    lo, hi = voxels()
    allowed = np.ones(len(lo), bool)
    pocket = (np.isclose(lo[:, 0], .6) & (lo[:, 1] <= 0) & (hi[:, 1] >= 0)
              & (lo[:, 2] <= 0) & (hi[:, 2] >= 0))
    allowed[pocket] = False
    fit = fit_ellipsoid(lo, hi, allowed, source_radius=.5, clearance=.01)
    assert fit["status"] == "fitted"
    assert fit["xc"]+fit["a"] <= .59+1e-6
    assert fit["volume_m3"] > .01
    allowed[:] = True
    off_axis = (np.isclose(lo[:, 0], .4) & np.isclose(lo[:, 1], .2)
                & np.isclose(lo[:, 2], 0))
    allowed[off_axis] = False
    fit = fit_ellipsoid(lo, hi, allowed, source_radius=.5, clearance=.01)
    closest_x = np.clip(fit["xc"], .39, .61)
    assert ((closest_x-fit["xc"])/fit["a"])**2 + (.19/fit["b"])**2 >= 1-1e-7
    assert verify_ellipsoid(fit, lo, hi, allowed)["feasible"]


def test_sampler_is_full_3d_and_rejects_excluded_cells(tmp_path):
    from _hybridmesh import box_mesh
    from plumetools.hybridplume.meshgeom import read_mesh
    from plumetools.hybridplume.fitted import sample_containing_cells
    from plumetools.hybridplume.workflow import WorkflowError
    box_mesh(tmp_path, [0, .1, .2], [-.1, 0, .1], [-.1, 0, .1], .1)
    mesh = read_mesh(tmp_path)
    mask = np.ones(8, bool)
    points = [[.025, -.025, -.025], [.175, .025, .025]]
    assert sample_containing_cells(mesh, points, mask, axes=(0, 1, 2)).tolist() == [0, 7]
    mask[7] = False
    with pytest.raises(WorkflowError, match="accepted|allowed"):
        sample_containing_cells(mesh, points, mask, axes=(0, 1, 2))


def test_hpc_profile_is_3d_kn001_with_predicted_budget():
    from plumetools.hybridplume import config, mesh
    cfg = config.load(STUDY, "hpc-kn1em3-ellipsoid")
    source = config.build_source(cfg)
    plan = mesh.plan(cfg, source)
    assert source.knudsen == .001 and source.depth_m is None
    assert plan.divisions == (160, 240, 240)
    assert cfg.get("interface.mode") == "fitted-ellipsoid"
    assert cfg.get("resources.default_ranks") == 192
