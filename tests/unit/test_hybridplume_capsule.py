"""Exact conservative capsule candidate sweep and flux screen; synthetic data only."""
import math
from types import SimpleNamespace

import numpy as np
import pytest

H = 0.05          # cell size: every limit used below lies on a cell edge
CLEAR = 0.01


def voxels():
    x, y, z = np.meshgrid(np.arange(40)*H, np.arange(-12, 12)*H,
                          np.arange(-12, 12)*H, indexing="ij")
    lower = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    return lower, lower+H


def radial(lo, hi):
    """Closest radial distance of each voxel footprint to the x-axis."""
    return np.hypot(*np.maximum(np.maximum(lo[:, 1:], -hi[:, 1:]), 0).T)


def sweep(lo, hi, retained, source_radius=.5):
    from plumetools.hybridplume.capsule import feasible_candidates
    return feasible_candidates(lo, hi, retained, source_radius=source_radius,
                               clearance=CLEAR, radius_samples=12, length_samples=24)


def test_unbounded_column_reaches_source_and_domain_limits():
    lo, hi = voxels()
    found = sweep(lo, hi, np.ones(len(lo), bool))
    assert found
    best = max(found, key=lambda c: c["volume_m3"])
    domain_pad = max(CLEAR, .25*H)
    assert abs(best["R"]-(.5-CLEAR)) < 1e-12
    assert abs(best["L"]+best["a"]-(2.0-domain_pad)) < 1e-9
    assert best["a"] > 0 and best["L"] > 0


def test_step_sets_length_limit():
    lo, hi = voxels()
    slab = (lo[:, 0] >= .6-1e-9) & (hi[:, 0] <= .8+1e-9)
    # Disk r < 0.3 removed: nothing may cross x = 0.6.
    retained = ~(slab & (radial(lo, hi) < .3-1e-9))
    found = sweep(lo, hi, retained)
    assert found
    assert all(c["L"] <= .6 and c["L"]+c["a"] <= .6 for c in found)
    # Annulus r >= 0.3 removed: wide capsules stop, narrow ones pass through.
    retained = ~(slab & (radial(lo, hi) >= .3-1e-9))
    found = sweep(lo, hi, retained)
    assert all(c["L"] <= .6 for c in found if c["R"] > .3-CLEAR)
    assert any(c["L"] > .8 for c in found if c["R"] < .3-CLEAR)


def test_tip_only_exclusion_limits_a():
    lo, hi = voxels()
    box = ((lo[:, 0] >= 1.2-1e-9) & (hi[:, 0] <= 1.4+1e-9)
           & (radial(lo, hi) < 1e-9))
    found = sweep(lo, hi, ~box)
    assert found
    assert all(c["L"]+c["a"] <= 1.2 for c in found)
    assert max(c["L"]+c["a"] for c in found) > 1.15


def test_fit_handles_zero_length_cylinder():
    from plumetools.hybridplume.shapes import Capsule
    lo, hi = voxels()
    # Exclusion right at the exit, outside r = 0.2: the exit disc itself is
    # blocked beyond r = 0.19, so only narrow capsules survive.
    at_exit = (hi[:, 0] <= .05+1e-9) & (radial(lo, hi) >= .2-1e-9)
    found = sweep(lo, hi, ~at_exit, source_radius=.3)
    assert found and all(c["R"] <= .2-CLEAR for c in found)
    tips = [c for c in found if c["L"] == 0]
    assert tips and all(c["a"] > 0 for c in tips)
    # Exclusion just downstream of the exit leaves room only for short shapes;
    # the L = 0 candidates are pure tips sitting on the source.
    near = ((lo[:, 0] >= .05-1e-9) & (hi[:, 0] <= .1+1e-9)
            & (radial(lo, hi) >= .2-1e-9))
    found = sweep(lo, hi, ~near)
    tips = [c for c in found if c["L"] == 0]
    assert tips and all(c["a"] > 0 for c in tips)
    wide = max(tips, key=lambda c: c["R"])
    rho = .2-CLEAR
    expected = (.05-CLEAR)/np.sqrt(1-(rho/wide["R"])**2)
    assert abs(wide["a"]-expected) < 1e-9
    for c in tips:
        assert Capsule(c["R"], 0.0, c["a"]).verify(lo, hi, ~near, CLEAR)["feasible"]
    assert len({(c["R"], c["L"]) for c in found}) == len(found)


def test_every_candidate_is_certified():
    from plumetools.hybridplume.shapes import Capsule
    lo, hi = voxels()
    rng = np.random.default_rng(3)
    retained = np.ones(len(lo), bool)
    downstream = np.flatnonzero(lo[:, 0] >= .5-1e-9)
    retained[rng.choice(downstream, 400, replace=False)] = False
    found = sweep(lo, hi, retained)
    assert found
    for c in found:
        assert set(c) == {"R", "L", "a", "volume_m3", "verification"}
        assert c["verification"]["feasible"]
        shape = Capsule(c["R"], c["L"], c["a"])
        assert shape.verify(lo, hi, retained, CLEAR)["feasible"]
        assert abs(c["volume_m3"]-shape.volume()) < 1e-15
    # Nothing feasible: empty list, not an exception.
    assert sweep(lo, hi, np.zeros(len(lo), bool)) == []


# --- pre-mesh half-range flux screen and selection -------------------------
D = 1.0


def half_spheroid_area(a, R):
    """Half the surface of a spheroid with polar semi-axis a and equatorial R."""
    if abs(a-R) < 1e-12:
        return 2*math.pi*R**2
    if a < R:
        e = math.sqrt(1-a*a/(R*R))
        return math.pi*R**2*(1+(1-e*e)/e*math.atanh(e))
    e = math.sqrt(1-R*R/(a*a))
    return math.pi*R**2*(1+a/(R*e)*math.asin(e))


def test_surface_area_matches_analytic():
    from plumetools.hybridplume.capsule import surface_samples
    from plumetools.hybridplume.shapes import Capsule
    for R, L, a in [(.45, 1.2, .6), (.3, .5, .9), (.4, 0.0, .4), (.2, 2.0, .05)]:
        shape = Capsule(R, L, a)
        points, normals, areas = surface_samples(shape, 40, 90, D)
        exact = 2*math.pi*R*L+half_spheroid_area(a, R)
        assert abs(areas.sum()/exact-1) < 5e-3
        assert points.shape == normals.shape == (len(areas), 3)
        assert np.allclose(np.linalg.norm(normals, axis=1), 1)
        # Samples lie on the surface with normals along the outward gradient.
        assert np.allclose(shape.q(points), 1)
        g = shape.outward_gradient(points)
        assert np.allclose(g/np.linalg.norm(g, axis=1)[:, None], normals)
        assert np.all(points[:, 0] > 0)  # no base disc; L = 0 is a pure tip


@pytest.fixture
def box(tmp_path):
    from _hybridmesh import box_mesh
    from plumetools.hybridplume.meshgeom import read_mesh
    from plumetools.hybridplume.fitted import RectilinearLocator
    edges = np.linspace(-.5, .5, 11)
    box_mesh(tmp_path, np.linspace(0, 1, 11), edges, edges, .5)
    mesh = read_mesh(tmp_path)
    return mesh, RectilinearLocator(mesh, axes=(0, 1, 2))


def cfg(screen=None, criteria=None):
    data = {"interface.fit.screen": screen or {}, "interface.criteria": criteria or {}}
    return SimpleNamespace(get=lambda key, default=None: data.get(key, default))


def state(mesh, velocity):
    n = mesh.n_cells
    return SimpleNamespace(rho=np.full(n, 1e-3), T=np.full(n, 300.0),
                           U=np.asarray(velocity(mesh.cell_centres), float))


def candidate(R, L, a):
    from plumetools.hybridplume.shapes import Capsule
    return {"R": R, "L": L, "a": a, "volume_m3": Capsule(R, L, a).volume(),
            "verification": {"feasible": True}}


def test_screen_return_is_small_for_radial_outflow(box):
    from plumetools.hybridplume.capsule import screen
    from plumetools.hybridplume.physics import Gas
    mesh, locator = box
    # Xenon at 300 K: 500 m/s is a normal speed ratio of about 2.6 (argon's
    # 1.4 cannot give a return fraction below 1e-3 even for exactly normal flow).
    gas = Gas(name="Xe", mass_kg=2.18e-25)
    st = state(mesh, lambda c: 500*c/np.linalg.norm(c, axis=1)[:, None])
    out = screen([candidate(.4, .1, .4)], locator, st, np.ones(mesh.n_cells, bool),
                 gas, D, cfg())
    assert out[0]["screen_feasible"] is True
    assert out[0]["return_mass_fraction"] < 1e-3
    assert out[0]["subsonic_area_fraction"] == 0.0
    assert out[0]["R"] == .4 and "verification" in out[0]


def test_screen_return_is_large_for_axial_flow_on_cylinder(box):
    from plumetools.hybridplume.capsule import screen
    from plumetools.hybridplume.physics import Gas
    mesh, locator = box
    st = state(mesh, lambda c: np.tile([500.0, 0, 0], (len(c), 1)))
    out = screen([candidate(.3, .6, .2)], locator, st, np.ones(mesh.n_cells, bool),
                 Gas(), D, cfg())
    assert out[0]["screen_feasible"] is True
    assert out[0]["return_mass_fraction"] > .3
    assert out[0]["subsonic_area_fraction"] > .5


def test_screen_marks_unsampleable_candidate_infeasible(box):
    from plumetools.hybridplume.capsule import screen
    from plumetools.hybridplume.physics import Gas
    mesh, locator = box
    st = state(mesh, lambda c: np.tile([500.0, 0, 0], (len(c), 1)))
    retained = mesh.cell_centres[:, 0] < .8
    # Tip touching excluded cells, a clean candidate, one leaving the domain.
    pool = [candidate(.3, .5, .4), candidate(.3, .3, .3), candidate(.3, .9, .3)]
    out = screen(pool, locator, st, retained, Gas(), D,
                 cfg({"axial_points_per_D": 20, "azimuthal_points": 36}))
    assert [c["screen_feasible"] for c in out] == [False, True, False]
    for c in (out[0], out[2]):
        assert c["return_mass_fraction"] is None and c["subsonic_area_fraction"] is None
    assert out[1]["return_mass_fraction"] is not None


def test_screen_raises_on_corrupt_inputs(box):
    from plumetools.hybridplume.capsule import screen
    from plumetools.hybridplume.physics import Gas
    from plumetools.hybridplume.workflow import WorkflowError
    mesh, locator = box
    st = state(mesh, lambda c: np.tile([500.0, 0, 0], (len(c), 1)))
    pool = [candidate(.3, .3, .3)]
    for bad in (np.ones(mesh.n_cells-1, bool), np.ones(mesh.n_cells)):
        with pytest.raises(WorkflowError, match="retained"):
            screen(pool, locator, st, bad, Gas(), D, cfg())
    retained = np.ones(mesh.n_cells, bool)
    wall = locator.locate([[.05, .25, .05]], retained)[0]  # under the cylinder wall
    for field, value in (("rho", 0.0), ("T", -1.0), ("U", np.nan)):
        corrupt = SimpleNamespace(rho=st.rho.copy(), T=st.T.copy(), U=st.U.copy())
        getattr(corrupt, field)[wall] = value
        with pytest.raises(WorkflowError, match="invalid CFD state"):
            screen(pool, locator, corrupt, retained, Gas(), D, cfg())


def screened(R, L, a, ret, sub, ok=True):
    c = candidate(R, L, a)
    c.update(return_mass_fraction=ret if ok else None,
             subsonic_area_fraction=sub if ok else None, screen_feasible=ok)
    return c


def test_select_prefers_largest_passing_volume():
    from plumetools.hybridplume.capsule import select
    pool = [screened(.4, 1.0, .5, .01, .1), screened(.45, 1.5, .5, .02, .2),
            screened(.45, 2.0, .5, .2, .1),          # larger, fails return
            screened(.45, 1.8, .5, .01, .5),         # larger, fails subsonic
            screened(.49, 2.5, .5, 0, 0, ok=False)]  # largest, unsampleable
    out = select(pool, {})
    assert out["status"] == "fitted"
    assert out["shape"] == {"shape": "capsule", "R": .45, "L": 1.5, "a": .5}
    assert out["screen"]["status"] == "passed"
    assert out["screen"]["candidates"] == pool
    # Explicit limits are honoured.
    out = select(pool, {"max_return_mass_fraction": .3, "max_subsonic_normal_area_fraction": .25})
    assert out["shape"]["L"] == 2.0


def test_select_reports_no_passing_candidate():
    from plumetools.hybridplume.capsule import select
    pool = [screened(.4, 1.0, .5, .3, .1), screened(.45, 1.5, .5, .2, .9),
            screened(.45, 2.0, .5, .4, .1), screened(.3, .1, .1, 0, 0, ok=False)]
    out = select(pool, {})
    assert out["status"] == "fitted"
    assert out["shape"]["L"] == 1.5  # lowest return among screen-feasible
    assert out["screen"]["status"] == "no_passing_candidate"
    for empty in ([], [screened(.3, .1, .1, 0, 0, ok=False)]):
        out = select(empty, {})
        assert out["status"] == "no_feasible_capsule" and out["reason"]
        assert "shape" not in out


def test_sample_containing_cells_unchanged(box):
    from plumetools.hybridplume.fitted import sample_containing_cells
    from plumetools.hybridplume.workflow import WorkflowError
    mesh, locator = box
    points = [[.05, -.45, -.45], [.95, .45, .45]]
    mask = np.ones(mesh.n_cells, bool)
    expected = sample_containing_cells(mesh, points, mask, axes=(0, 1, 2))
    assert expected.tolist() == [0, mesh.n_cells-1]
    assert np.array_equal(locator.locate(points, mask), expected)
    mask[-1] = False
    for call in (lambda: locator.locate(points, mask),
                 lambda: sample_containing_cells(mesh, points, mask, axes=(0, 1, 2))):
        with pytest.raises(WorkflowError, match="accepted/allowed"):
            call()
    with pytest.raises(WorkflowError, match="outside the CFD domain"):
        locator.locate([[1.5, 0, 0]], np.ones(mesh.n_cells, bool))
    with pytest.raises(WorkflowError, match="finite"):
        locator.locate([[np.nan, 0, 0]], np.ones(mesh.n_cells, bool))


def _stl_shape():
    from plumetools.hybridplume import shapes
    return shapes.Capsule(.5, .8, .6)


def test_stl_is_closed_and_outward():
    from collections import Counter
    from plumetools.hybridplume.capsule_mesh import capsule_stl
    shape, back = _stl_shape(), .05
    v, t = capsule_stl(shape, .002, back)
    edges = Counter()
    for tri in t:
        for i in range(3):
            edges[tuple(sorted((tri[i], tri[(i+1) % 3])))] += 1
    assert set(edges.values()) == {2}
    p = v[t]
    volume = np.einsum("ij,ij->i", p[:, 0], np.cross(p[:, 1], p[:, 2])).sum()/6
    expected = shape.volume()+math.pi*shape.R**2*back
    assert volume > 0 and volume == pytest.approx(expected, rel=.01)


@pytest.mark.parametrize("tol", [.002, .01])
@pytest.mark.parametrize("length", [.8, 0.0])
def test_stl_chord_error_within_tolerance(length, tol):
    from plumetools.hybridplume import shapes
    from plumetools.hybridplume.capsule_mesh import capsule_stl
    shape = shapes.Capsule(.5, length, .6)
    v, t = capsule_stl(shape, tol, .05)
    n = 12   # independent, finer grid than the generator's, includes edge midpoints
    w = np.array([(i/n, j/n, (n-i-j)/n) for i in range(n+1) for j in range(n+1-i)])
    pts = np.einsum("sk,tkd->tsd", w, v[t]).reshape(-1, 3)
    pts = pts[pts[:, 0] >= 0]
    assert np.max(np.abs(np.sqrt(shape.q(pts))-1)) <= tol


def test_snappy_dict_uses_trisurface_and_features(tmp_path):
    from plumetools.hybridplume import config
    from plumetools.hybridplume.capsule_mesh import background, write_mesh_inputs
    from pathlib import Path
    study = Path(__file__).resolve().parents[2]/"cases/cai2012-hybrid/hybrid.yaml"
    cfg = config.load(study, "hpc-kn1em3-ellipsoid")
    source = config.build_source(cfg)
    from plumetools.hybridplume import shapes
    shape = shapes.Capsule(.1, .3, .15)   # R matches the source radius
    assert shape.q(np.array([background(cfg, source)["seed"]]))[0] > 1
    paths = write_mesh_inputs(tmp_path, cfg, source, shape)
    snappy = (tmp_path/"mesh-build/system/snappyHexMeshDict").read_text()
    assert "interface.stl\n    {\n        type triSurfaceMesh;" in snappy
    assert "name interface;" in snappy and 'file "interface_rim.obj";' in snappy
    assert "level (3 3);" in snappy and "            level 3;" in snappy
    assert "explicitFeatureSnap     true;" in snappy and "implicitFeatureSnap     false;" in snappy
    assert "No feature edges" not in snappy and ".eMesh" not in snappy
    assert (tmp_path/"mesh-build/constant/triSurface/interface.stl").read_text().startswith("solid interface")
    assert not (tmp_path/"mesh-build/system/surfaceFeatureExtractDict").exists()
    assert len(paths) == 8 and all(p.is_file() for p in paths)
    # The exit-plane rim: a closed loop at x = 0, radius R, within the STL chord tolerance.
    tol = .1*float(cfg.get("interface.fit.geometry_relative_tolerance", .02))
    lines = (tmp_path/"mesh-build/constant/triSurface/interface_rim.obj").read_text().splitlines()
    v = np.array([[float(t) for t in ln.split()[1:]] for ln in lines if ln.startswith("v ")])
    segs = [tuple(int(t)-1 for t in ln.split()[1:]) for ln in lines if ln.startswith("l ")]
    n = len(v)
    assert n >= 8 and len(segs) == n
    assert sorted(segs) == sorted((i, (i+1) % n) for i in range(n))   # one closed loop
    assert np.all(v[:, 0] == 0) and np.allclose(np.hypot(v[:, 1], v[:, 2]), shape.R, rtol=1e-11)
    assert 1-math.cos(math.pi/n) <= tol                                  # sagitta / R
    assert 1-math.cos(math.pi/(n-1)) > tol                               # and not over-refined


def test_capsule_snappy_edits_refuse_a_changed_renderer(tmp_path, monkeypatch):
    from plumetools.hybridplume import capsule_mesh, config, shapes
    cfg = config.load(STUDY, "hpc-kn1em3-ellipsoid")
    source = config.build_source(cfg)
    real = capsule_mesh.render_snappy_dict
    for marker in ("    interface\n    {", "features ();", "explicitFeatureSnap     false;"):
        monkeypatch.setattr(capsule_mesh, "render_snappy_dict",
                            lambda *a, m=marker, **k: real(*a, **k).replace(m, "", 1))
        with pytest.raises(RuntimeError, match="not found"):
            capsule_mesh.write_mesh_inputs(tmp_path, cfg, source, shapes.Capsule(.1, .3, .15))


def test_ellipsoid_snappy_dict_keeps_feature_snapping_off(tmp_path):
    from plumetools.hybridplume import config
    from plumetools.hybridplume.ellipsoid_mesh import write_mesh_inputs
    cfg = config.load(STUDY, "hpc-kn1em3-ellipsoid")
    source = config.build_source(cfg)
    write_mesh_inputs(tmp_path, cfg, source, {"xc": 0.0, "a": .3, "b": .12, "source_radius_m": .1})
    snappy = (tmp_path/"mesh-build/system/snappyHexMeshDict").read_text()
    assert "explicitFeatureSnap     false;" in snappy and "    features ();" in snappy
    assert "interface_rim.obj" not in snappy


# --- Task 9: capsule wired into spatial.interface and Allrun ------------------

STUDY = __import__("pathlib").Path(__file__).resolve().parents[2]/"cases/cai2012-hybrid/hybrid.yaml"
CFD_TIME = "0.001"


def capsule_case(tmp_path, monkeypatch, *, retained=None):
    """Synthetic fitted-capsule case: box CFD mesh, current breakdown, radial outflow."""
    import yaml
    from _hybridmesh import box_mesh
    from plumetools.hybridplume import config, fitted, spatial, workflow
    from plumetools.hybridplume.meshgeom import read_mesh
    base = config.load(STUDY, "hpc-kn1em3-capsule")
    doc = config.deep_merge(base.doc, {"interface": {"fit": {"screen": {
        "radius_samples": 3, "length_samples": 4, "axial_points_per_D": 10, "azimuthal_points": 16}}}})
    (tmp_path/workflow.RESOLVED).write_text(yaml.safe_dump({"profile": base.profile, "config": doc}))
    cfg = workflow.load_case(tmp_path)
    edges = np.linspace(-.3, .3, 13)
    box_mesh(tmp_path/"cfd", np.linspace(0, .6, 13), edges, edges, .1)
    (tmp_path/"cfd"/CFD_TIME).mkdir()
    for name in ("p", "rho", "T", "U"):
        (tmp_path/"cfd"/CFD_TIME/name).write_text(f"field {name}\n")
    mesh = read_mesh(tmp_path/"cfd")
    if retained is None:
        retained = mesh.cell_centres[:, 0] < .4   # continuum region ends at x = 0.4
    np.savez(workflow.results_dir(tmp_path)/"breakdown.npz", time=CFD_TIME, retained=retained,
             Kn_GLL=np.full(mesh.n_cells, .01), input_identity=fitted.breakdown_identity(tmp_path, cfg, CFD_TIME))
    centres = mesh.cell_centres
    velocity = 3000*centres/np.linalg.norm(centres, axis=1)[:, None]
    monkeypatch.setattr(spatial.bd, "read_state", lambda *args: SimpleNamespace(
        rho=np.full(mesh.n_cells, 1e-3), T=np.full(mesh.n_cells, 300.0), U=velocity))
    return cfg, mesh


def source_diameter(cfg):
    from plumetools.hybridplume import config
    return config.build_source(cfg).diameter_m


def test_capsule_interface_writes_fit_screen_and_mesh_inputs(tmp_path, monkeypatch):
    import json
    from plumetools.hybridplume import fitted, shapes, spatial
    cfg, _ = capsule_case(tmp_path, monkeypatch)
    doc = spatial.interface(tmp_path)
    res = tmp_path/"results"
    fit = json.loads((res/"capsule.json").read_text())
    assert fit["status"] == "fitted" and fit["shape"] == "capsule"
    assert fit["screen"]["status"] in ("passed", "no_passing_candidate")
    assert fit["screen"]["candidates"] and fit["time"] == CFD_TIME
    assert fit["clearance_m"] == pytest.approx(.25*.05*source_diameter(cfg))
    # capsule.json round-trips through load_shape to the selected capsule.
    _, shape = spatial.load_shape(tmp_path, cfg)
    assert shape == shapes.Capsule(fit["R"], fit["L"], fit["a"])
    assert isinstance(shape, shapes.Capsule) and shape.L + shape.a <= .4
    on_disk = json.loads((res/"interface.json").read_text())
    assert on_disk == json.loads(json.dumps(doc))
    assert doc["mode"] == "fitted-capsule" and doc["status"] == "mesh_pending"
    assert doc["fit"]["R"] == shape.R and "screen" not in doc["fit"]
    selected = doc["screen"]["selected"]
    assert doc["screen"]["status"] == fit["screen"]["status"]
    assert (selected["R"], selected["L"], selected["a"]) == (shape.R, shape.L, shape.a)
    assert selected["screen_feasible"] and selected["return_mass_fraction"] < .05
    assert doc["fitted_identity"] == fitted._identity(tmp_path, cfg, shape, CFD_TIME)
    build = tmp_path/"mesh-build"
    assert (build/"constant/triSurface/interface_rim.obj").is_file()
    assert not (build/"system/surfaceFeatureExtractDict").exists()
    assert (build/"constant/triSurface/interface.stl").read_text().startswith("solid interface")
    assert "interface.stl" in (build/"system/snappyHexMeshDict").read_text()


def test_capsule_interface_reports_no_feasible_capsule(tmp_path, monkeypatch):
    import json
    from plumetools.hybridplume import interface as itf, spatial, workflow
    capsule_case(tmp_path, monkeypatch, retained=np.zeros(12**3, bool))
    doc = spatial.interface(tmp_path)
    assert doc["mode"] == "fitted-capsule" and doc["status"] == itf.STATUS_NONE
    assert doc["reasons"] and "capsule" in doc["reasons"][0]
    assert json.loads((tmp_path/"results/capsule.json").read_text())["status"] == "no_feasible_capsule"
    assert not (tmp_path/"mesh-build").exists()
    with pytest.raises(workflow.WorkflowError, match="no feasible capsule"):
        spatial.load_shape(tmp_path, workflow.load_case(tmp_path))


def test_capsule_interface_refuses_existing_ellipsoid_mesh(tmp_path, monkeypatch):
    import json
    from _hybridmesh import box_mesh
    from plumetools.hybridplume import spatial, workflow
    capsule_case(tmp_path, monkeypatch)
    # A DSMC mesh already built for an ellipsoid fit of this case.
    box_mesh(tmp_path/"dsmc", [0, .1], [-.1, .1], [-.1, .1], .1)
    (tmp_path/"results/interface.json").write_text(json.dumps(
        {"mode": "fitted-ellipsoid", "fitted_identity": "ellipsoid-identity", "time": CFD_TIME}))
    with pytest.raises(workflow.WorkflowError, match="existing DSMC mesh"):
        spatial.interface(tmp_path)
    assert not (tmp_path/"results/capsule.json").exists()
    assert not (tmp_path/"mesh-build").exists()


def test_dsmc_setup_rewrites_capsule_mesh_inputs(tmp_path, monkeypatch):
    from plumetools.hybridplume import shapes, spatial
    capsule_case(tmp_path, monkeypatch)
    spatial.interface(tmp_path)
    calls = []
    monkeypatch.setattr(spatial.capsule_mesh, "write_mesh_inputs",
                        lambda case, cfg, source, shape: calls.append(shape))
    monkeypatch.setattr(spatial, "write_mesh_inputs",
                        lambda *args: pytest.fail("ellipsoid mesh writer used for a capsule"))
    monkeypatch.setattr(spatial.resources, "assess_configuration",
                        lambda *args: {"full_dsmc": None})
    monkeypatch.setattr(spatial.dsmc, "plan_run", lambda *args: SimpleNamespace(describe=lambda: {}))
    monkeypatch.setattr(spatial.dsmc, "write_case", lambda *args: None)
    (tmp_path/"dsmc").mkdir()
    spatial.dsmc_setup(tmp_path, np_ranks=1)
    assert len(calls) == 1 and isinstance(calls[0], shapes.Capsule)


def test_plots_draw_capsule_section_under_capsule_name(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    from plumetools.hybridplume import spatial
    capsule_case(tmp_path, monkeypatch)
    spatial.interface(tmp_path)
    out = spatial.plots(tmp_path)
    assert [p.name for p in out] == ["breakdown_capsule.png"]
    assert not (tmp_path/"results/breakdown_ellipsoid.png").exists()


def test_allrun_snaps_without_surface_feature_extract():
    # The rim edge mesh is written by capsule_mesh; there is no extraction step.
    text = (STUDY.parent/"Allrun").read_text()
    start = text.index('cd "$CASE/mesh-build"')
    block = text[start:text.index("run topoSet topoSet", start)]
    assert "surfaceFeatureExtract" not in text
    assert block.index("run blockMesh blockMesh") < block.index("snappyHexMesh")


def test_allrun_spatial_fit_modes_match_config():
    import re
    from plumetools.hybridplume import config
    text = (STUDY.parent/"Allrun").read_text()
    body = text[text.index("spatial_fit() {"):]
    body = body[:body.index("\n}\n")]
    pattern = re.search(r"^\s*([\w|-]+)\)\s*return 0", body, re.M).group(1)
    assert set(pattern.split("|")) == set(config.SPATIAL_FIT_MODES)


def test_certify_facets_accepts_capsule_surface_and_refuses_outward_vertex():
    from plumetools.hybridplume import shapes, spatial, workflow
    shape = shapes.Capsule(.1, .3, .15)
    lo, hi = voxels()
    retained = np.ones(len(lo), bool)
    phi = np.linspace(0, 2*math.pi, 12, endpoint=False)
    t = np.linspace(0, math.pi/2, 5)
    x = np.concatenate((np.zeros(12), np.full(12, .15), .3+.15*np.cos(np.repeat(t, 12))))
    r = np.concatenate((np.full(24, .1), .1*np.sin(np.repeat(t, 12))))
    p = np.tile(phi, 7)
    vertices = np.column_stack((x, r*np.cos(p), r*np.sin(p)))
    checks = spatial.certify_facets(vertices, shape, lo, hi, retained, .02)
    assert checks["max_relative_radius_deviation"] < 1e-12
    pushed = vertices.copy()
    pushed[3, 1:] *= 1.03                       # 3 % outward on the exit-plane rim
    with pytest.raises(workflow.WorkflowError, match="tolerance"):
        spatial.certify_facets(pushed, shape, lo, hi, retained, .02)


@pytest.mark.parametrize("sign,match", [(1, "reversed normals"), (-1, "source coverage")])
def test_geometry_checks_capsule_normals_with_outward_gradient(monkeypatch, sign, match):
    # One cylinder-side face at y = R: the interface area vector points into
    # the DSMC region (-y), so the outward normal -face_area must be +y.
    from types import SimpleNamespace as NS
    from plumetools.hybridplume import config, shapes, spatial, workflow
    cfg = config.load(STUDY, "hpc-kn1em3-capsule")
    source = config.build_source(cfg)
    mesh = NS(cell_volumes=np.ones(1), patches={"interface": NS(type="patch")},
              patch_slice=lambda name: slice(0, 1), offsets=np.array([0, 3]),
              labels=np.arange(3), points=np.zeros((3, 3)),
              face_areas=np.array([[0, sign*1e-4, 0]]), face_centres=np.array([[.15, .1, 0]]))
    cmesh = NS(cell_centres=np.array([[.5, 0, 0]]), cell_size=lambda: np.ones((1, 3)))
    monkeypatch.setattr(spatial, "certify_facets", lambda *args: {})
    with pytest.raises(workflow.WorkflowError, match=match):
        spatial._geometry(mesh, cmesh, {"retained": np.ones(1, bool)},
                          shapes.Capsule(.1, .3, .15), cfg, source)
