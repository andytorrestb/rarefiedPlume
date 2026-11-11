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
