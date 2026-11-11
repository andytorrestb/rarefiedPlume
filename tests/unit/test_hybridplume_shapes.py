"""Shared fitted-shape interface checks (spheroid and capsule)."""
import math
import numpy as np
import pytest

from plumetools.hybridplume import shapes
from plumetools.hybridplume.ellipsoid import verify_ellipsoid
from test_hybridplume_ellipsoid import voxels

SPHEROID = shapes.Spheroid(xc=0.1, a=0.5, b=0.3)
CAPSULE = shapes.Capsule(R=0.3, L=0.4, a=0.35)


def surface_points(shape, n=200):
    pts = []
    for x, y in shape.section(n):
        for phi in (0.0, 0.7, 2.1):
            pts.append((x, y*math.cos(phi), y*math.sin(phi)))
    return np.array(pts)


@pytest.mark.parametrize("shape", [SPHEROID, CAPSULE])
def test_q_is_one_on_surface(shape):
    pts = surface_points(shape)
    assert np.allclose(shape.q(pts), 1.0, atol=1e-9)
    assert shape.deviation(pts) < 1e-9


@pytest.mark.parametrize("shape,inside", [(SPHEROID, (0.1, 0, 0)), (CAPSULE, (0.2, 0, 0))])
def test_outward_gradient_points_away_from_interior(shape, inside):
    pts = surface_points(shape)
    g = shape.outward_gradient(pts)
    assert np.all(np.einsum("ij,ij->i", g, pts-np.array(inside)) > 0)


def test_capsule_volume_and_base():
    c = CAPSULE
    assert c.volume() == pytest.approx(math.pi*c.R**2*c.L + 2/3*math.pi*c.R**2*c.a)
    assert c.base_radius() == c.R
    assert c.params() == [c.R, c.L, c.a]


def test_spheroid_volume_and_base():
    from plumetools.hybridplume.ellipsoid import clipped_volume
    s = SPHEROID
    assert s.volume() == pytest.approx(clipped_volume([s.xc, s.a, s.b]))
    assert s.base_radius() == pytest.approx(s.b*math.sqrt(1-(s.xc/s.a)**2))
    assert s.params() == [s.xc, s.a, s.b]


def test_spheroid_verify_matches_verify_ellipsoid():
    lo, hi = voxels()
    allowed = np.ones(len(lo), bool)
    allowed[np.isclose(lo[:, 0], .6) & (lo[:, 1] <= 0) & (hi[:, 1] >= 0)] = False
    s = shapes.Spheroid(xc=0.1, a=0.4, b=0.3)
    expected = verify_ellipsoid({"xc": .1, "a": .4, "b": .3, "clearance_m": .01},
                                lo, hi, allowed)
    assert s.verify(lo, hi, allowed, 0.01) == expected


def test_capsule_grown_fails_where_unit_passes():
    lo, hi = voxels()
    allowed = np.ones(len(lo), bool)
    blocked = np.isclose(lo[:, 0], .6) & (lo[:, 1] <= 0) & (hi[:, 1] >= 0)
    allowed[blocked] = False
    c = shapes.Capsule(R=0.2, L=0.3, a=0.3-1e-3)
    assert c.verify(lo, hi, allowed, 0.0)["feasible"]
    assert not c.grown(1.01).verify(lo, hi, allowed, 0.0)["feasible"]
    assert c.grown(1.01).L == c.L
    assert c.grown(1.01).R == pytest.approx(c.R*1.01)


def test_capsule_cylinder_blocked_by_side_voxel():
    lo, hi = voxels()
    allowed = np.ones(len(lo), bool)
    side = np.isclose(lo[:, 0], .2) & np.isclose(lo[:, 1], .2) & np.isclose(lo[:, 2], 0)
    allowed[side] = False
    out = shapes.Capsule(R=0.3, L=0.4, a=0.1).verify(lo, hi, allowed, 0.0)
    assert not out["feasible"] and out["intersected_excluded_voxels"] >= 1
    assert set(out) == {"feasible", "intersected_excluded_voxels",
                        "excluded_voxels_checked", "min_scaled_distance_squared"}
    assert shapes.Capsule(R=0.15, L=0.4, a=0.1).verify(lo, hi, allowed, 0.0)["feasible"]


def test_capsule_zero_length_is_pure_tip():
    lo, hi = voxels()
    allowed = np.ones(len(lo), bool)
    c = shapes.Capsule(R=0.3, L=0.0, a=0.3)
    assert c.verify(lo, hi, allowed, 0.0)["feasible"]
    assert c.volume() == pytest.approx(2/3*math.pi*0.3**2*0.3)
    sec = c.section(20)
    assert np.isfinite(sec).all() and sec[:, 0].min() >= -1e-12
    assert np.allclose(c.q(surface_points(c)), 1.0)
    allowed[np.isclose(lo[:, 0], 0) & np.isclose(lo[:, 1], 0) & np.isclose(lo[:, 2], 0)] = False
    assert not c.verify(lo, hi, allowed, 0.0)["feasible"]


def test_json_round_trip_and_legacy():
    for s in (SPHEROID, CAPSULE):
        doc = s.to_json()
        assert doc["shape"] == ("spheroid" if s is SPHEROID else "capsule")
        assert shapes.from_json(doc) == s
    legacy = {"status": "fitted", "xc": 0.1, "a": 0.5, "b": 0.3, "volume_m3": 1.0,
              "source_radius_m": 0.2, "clearance_m": 0.01}
    assert shapes.from_json(legacy) == shapes.Spheroid(0.1, 0.5, 0.3)
    with pytest.raises(ValueError):
        shapes.from_json({"shape": "cube"})
