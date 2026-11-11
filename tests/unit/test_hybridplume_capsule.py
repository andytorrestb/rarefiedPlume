"""Exact conservative capsule candidate sweep; synthetic voxels only."""
import numpy as np

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
