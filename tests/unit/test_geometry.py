"""Face centroids, normals, and the spherical-coordinate convention."""

from __future__ import annotations

import math

import numpy as np
import pytest

from plumetools.geometry import centroid, centroids, normal, normals, spherical

UNIT_TRIANGLE = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
UNIT_SQUARE = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])


# --------------------------------------------------------------------------- #
# centroid -- finding GP-01
# --------------------------------------------------------------------------- #

def test_triangle_centroid():
    np.testing.assert_allclose(centroid(UNIT_TRIANGLE), [1 / 3, 1 / 3, 0.0])


def test_square_centroid():
    np.testing.assert_allclose(centroid(UNIT_SQUARE), [0.5, 0.5, 0.0])


def test_regular_polygon_centroid_is_its_centre():
    n = 7
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    poly = np.stack([np.cos(ang), np.sin(ang), np.zeros(n)], axis=1)
    np.testing.assert_allclose(centroid(poly), [0, 0, 0], atol=1e-15)


def test_fix_is_bit_identical_on_triangles():
    """The GP-01 fix must not perturb the golden.

    Dividing by len(vertices) is the original expression with 3.0 replaced by
    len(), so on a 3-vertex face the two are the same computation. Asserted
    exactly, not approximately -- this is what lets the fix land inside a
    behaviour-preserving extraction.
    """
    rng = np.random.default_rng(0)
    for _ in range(200):
        tri = rng.normal(size=(3, 3))
        assert centroid(tri).tolist() == centroid(tri, legacy_triangle=True).tolist()


def test_legacy_divisor_inflates_quads_by_four_thirds():
    """Finding AR-01, reproduced exactly.

    The archived wake-cylinder inflow faces are quads, and the legacy /3.0 put
    every centroid at 4/3 of its true radius -- measured mean |r| 0.3946 m
    against a true 0.2960 m, ratio 1.3331.
    """
    true = centroid(UNIT_SQUARE)
    legacy = centroid(UNIT_SQUARE, legacy_triangle=True)
    np.testing.assert_allclose(legacy, np.asarray(true) * 4.0 / 3.0)


def test_centroid_rejects_degenerate_input():
    with pytest.raises(ValueError, match="at least 3"):
        centroid(np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
    with pytest.raises(ValueError, match=r"\(n, 3\)"):
        centroid(np.zeros((3, 2)))


def test_centroids_maps_over_faces():
    points = np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
    got = centroids([[0, 1, 2], [0, 2, 3]], points)
    assert got.shape == (2, 3)
    np.testing.assert_allclose(got[0], [2 / 3, 1 / 3, 0])


# --------------------------------------------------------------------------- #
# normal -- finding GP-02
# --------------------------------------------------------------------------- #

def test_normal_is_a_unit_vector():
    assert np.linalg.norm(normal(UNIT_TRIANGLE)) == pytest.approx(1.0)


def test_normal_flips_with_winding():
    reversed_tri = UNIT_TRIANGLE[::-1]
    np.testing.assert_allclose(normal(reversed_tri), -normal(UNIT_TRIANGLE), atol=1e-15)


def test_normal_orientation_matches_commit_b22f573():
    """n = (p0 - p2) x (p0 - p1), the orientation that commit settled on."""
    p = UNIT_TRIANGLE
    expected = np.cross(p[0] - p[2], p[0] - p[1])
    np.testing.assert_allclose(normal(p), expected / np.linalg.norm(expected))


def test_normal_rejects_collinear_face():
    line = np.array([[0.0, 0, 0], [1, 0, 0], [2, 0, 0]])
    with pytest.raises(ValueError, match="collinear"):
        normal(line)


def test_normals_maps_over_faces():
    points = np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
    assert normals([[0, 1, 2], [0, 2, 3]], points).shape == (2, 3)


# --------------------------------------------------------------------------- #
# spherical -- findings SM-04, SM-09
# --------------------------------------------------------------------------- #

def test_theta_is_measured_from_plus_z():
    """Not from the plume axis. See SM-04."""
    theta, _ = spherical(np.array([[0.0, 0, 1], [0, 0, -1], [1, 0, 0]]))
    np.testing.assert_allclose(theta, [0.0, math.pi, math.pi / 2], atol=1e-12)


def test_phi_is_azimuth_from_plus_x():
    _, phi = spherical(np.array([[1.0, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]))
    np.testing.assert_allclose(phi, [0, math.pi / 2, math.pi, -math.pi / 2], atol=1e-12)


def test_round_trip_through_cartesian():
    rng = np.random.default_rng(1)
    pts = rng.normal(size=(50, 3))
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    theta, phi = spherical(pts)
    back = np.stack([np.cos(phi) * np.sin(theta),
                     np.sin(phi) * np.sin(theta),
                     np.cos(theta)], axis=1)
    np.testing.assert_allclose(back, pts, atol=1e-12)


def test_fixed_radius_differs_from_per_point_radius():
    """The fixed radius is deliberate, not an approximation (SM-09).

    A centroid inside the nominal sphere gets a different polar angle under the
    two conventions; that difference is what tilts U off radial by 3.1e-4 in the
    reference case.
    """
    pt = np.array([[0.1, 0.1, 0.4]])
    fixed, _ = spherical(pt, radius=0.5)
    per_point, _ = spherical(pt, radius=None)
    assert fixed != pytest.approx(per_point)


def test_fixed_radius_smaller_than_the_points_raises():
    with pytest.raises(ValueError, match="exceeds 1"):
        spherical(np.array([[0.0, 0.0, 1.0]]), radius=0.5)


def test_unimplemented_polar_axis_raises():
    """The model as written only supports theta from +z; fail rather than pretend."""
    with pytest.raises(NotImplementedError, match="polar_axis"):
        spherical(np.array([[1.0, 0.0, 0.0]]), polar_axis="x")
