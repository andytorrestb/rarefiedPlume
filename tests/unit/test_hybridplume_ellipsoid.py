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
    assert cfg.get("resources.max_memory_GB")*1024**3 < 1200*10**9


def test_facet_enclosure_refuses_excluded_voxel():
    from plumetools.hybridplume.spatial import certify_facets
    from plumetools.hybridplume.workflow import WorkflowError
    lo = np.array([[.4, -.1, -.1]])
    hi = np.array([[.6, .1, .1]])
    # Vertices lie on a large enclosing surface whose interior crosses the hole.
    with pytest.raises(WorkflowError, match="excluded|continuum"):
        certify_facets(np.array([[1, 0, 0], [0, .5, 0], [0, 0, .5]]),
                       dict(xc=0, a=1, b=.5), lo, hi, np.array([False]), .1)


def test_peak_budget_counts_parcels_before_initialization():
    from plumetools.hybridplume import config, resources
    cfg = config.load(STUDY, "hpc-kn1em3-ellipsoid")
    result = resources.peak_memory_estimate(cfg, 9216000,
                                          dict(cells=1000000, parcels=1e9), ranks=192)
    assert not result["fits_memory_budget"]
    assert result["peak_GiB"] > 1200


def test_3d_precfd_bound_uses_exit_state_without_disk_quadrature(monkeypatch):
    from plumetools.cai2012 import analytical
    from plumetools.hybridplume import config, mesh, resources
    cfg = config.load(STUDY, "hpc-kn1em3-ellipsoid")
    source = config.build_source(cfg)
    monkeypatch.setattr(analytical, "density_ratio", lambda *args:
                        pytest.fail("a fitted 3D collisional sizing bound must not use disk quadrature"))
    estimate = resources.assess_configuration(cfg, source, mesh.plan(cfg, source))["full_dsmc"]
    assert estimate["densest_cell"]["number_density"] == source.number_density


def test_mesh_publication_never_publishes_mesher_processors(tmp_path):
    from _hybridmesh import box_mesh
    from plumetools.hybridplume.spatial import publish_mesh
    box_mesh(tmp_path/"mesh-build", [0, .1], [-.1, .1], [-.1, .1], .1)
    (tmp_path/"mesh-build/processor0").mkdir()
    (tmp_path/"dsmc/constant").mkdir(parents=True)
    publish_mesh(tmp_path)
    assert (tmp_path/"dsmc/constant/polyMesh/owner").is_file()
    assert not (tmp_path/"dsmc/processor0").exists()


def test_ellipsoid_post_uses_own_mesh_without_cfd_addressing(tmp_path):
    import json
    import yaml
    from _hybridmesh import box_mesh
    from plumetools.hybridplume import config, workflow
    cfg = config.load(STUDY, "hpc-kn1em3-ellipsoid")
    source = config.build_source(cfg)
    (tmp_path/workflow.RESOLVED).write_text(yaml.safe_dump(
        {"profile": cfg.profile, "digest": cfg.digest, "config": cfg.doc}))
    res = tmp_path/"results"
    res.mkdir()
    box_mesh(tmp_path/"dsmc", [.1, .2, .3], [-.005, 0, .005], [-.005, 0, .005], .1)
    td = tmp_path/"dsmc/0.000002"
    td.mkdir()
    n, u, temp = source.number_density, source.velocity, 300.0
    rho = n*source.gas.mass_kg
    for name, kind, row in (
        ("rhoNMean", "scalar", str(n)), ("rhoMMean", "scalar", str(rho)),
        ("momentumMean", "vector", f"({rho*u} 0 0)"),
        ("linearKEMean", "scalar", str(.5*rho*(u*u+3*source.gas.gas_constant*temp))),
        ("dsmcRhoNMean", "scalar", "12")):
        (td/name).write_text(f"internalField nonuniform List<{kind}>\n8\n(\n"
                            + "\n".join([row]*8)+"\n);\n")
    (res/"dsmc_setup.json").write_text(json.dumps({
        "run": {"average_start_s": 0, "end_time_s": 2e-6, "dt_s": 1e-6,
                "weight": n*2.5e-6/12}, "inflow_patches": [], "ranks": 1,
        "smoke_only": True, "interface_status": "unsuitable"}))
    (tmp_path/"dsmc/log.dsmcFoam").touch()
    from plumetools.hybridplume.meshgeom import read_mesh
    centres = read_mesh(tmp_path/"dsmc").cell_centres
    cloud = td/"lagrangian/dsmc"
    cloud.mkdir(parents=True)
    speeds = np.vstack((np.eye(3), -np.eye(3)))*np.sqrt(3*source.gas.gas_constant*temp)
    positions, velocities = [], []
    for cell, centre in enumerate(centres):
        for vel in np.tile(speeds, (2, 1)) + [u, 0, 0]:
            positions.append(f"({' '.join(map(str, centre))}) {cell}")
            velocities.append(f"({' '.join(map(str, vel))})")
    for name, rows in (("positions", positions), ("U", velocities)):
        (cloud/name).write_text(f"{len(rows)}\n(\n"+"\n".join(rows)+"\n)\n")
    doc = workflow.dsmc_post(tmp_path)
    assert doc["geometry"] == "cartesian" and doc["cells"] == 8
    assert doc["valid_cells"] == 8 and doc["smoke_only"]
    assert "collisionless limit" in doc["centreline_comparison"]["note"]
    assert not (res/"dsmc_to_cfd_cells.npy").exists()
    with np.load(res/"spatial_fields.npz") as fields:
        assert fields["centres"].shape == (8, 3)
        assert np.allclose(fields["T"], 300)
    assert (td/"spatial_n_over_ne").is_file()
    assert (tmp_path/"dsmc/spatial.foam").is_file()


def test_retry_refuses_a_changed_dsmc_mesh_before_reusing_cloud(tmp_path, monkeypatch):
    import json
    import yaml
    from _hybridmesh import box_mesh
    from plumetools.hybridplume import config, fitted, spatial, workflow
    cfg = config.load(STUDY, "hpc-kn1em3-ellipsoid")
    (tmp_path/workflow.RESOLVED).write_text(yaml.safe_dump(
        {"profile": cfg.profile, "digest": cfg.digest, "config": cfg.doc}))
    res = tmp_path/"results"
    res.mkdir()
    box_mesh(tmp_path/"dsmc", [0, .1], [-.1, .1], [-.1, .1], .1)
    mesh_dir = tmp_path/"dsmc/constant/polyMesh"
    identity = {name: workflow.provenance.sha256(mesh_dir/name)
                for name in ("points", "faces", "owner", "neighbour", "boundary")}
    (res/"ellipsoid.json").write_text(json.dumps({"status": "fitted"}))
    (res/"interface.json").write_text(json.dumps({"fitted_identity": "same", "time": "1"}))
    (res/"dsmc_setup.json").write_text(json.dumps(
        {"fitted_identity": "same", "mesh_identity": identity}))
    monkeypatch.setattr(fitted, "_identity", lambda *args: "same")
    monkeypatch.setattr(spatial.resources, "assess_configuration",
                        lambda *args: pytest.fail("mesh identity must be checked before resource analysis"))
    with (mesh_dir/"points").open("a") as stream:
        stream.write("\n// changed after initialization\n")
    with pytest.raises(workflow.WorkflowError, match="mesh.*changed|mesh.*identity"):
        spatial.dsmc_setup(tmp_path)
