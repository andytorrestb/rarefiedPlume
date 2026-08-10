"""The Cai 2012 geometry: the exit disk, the vacuum box, and X/D normalisation."""

from __future__ import annotations

import math

import numpy as np
import pytest

from _cai2012 import load_case
from plumetools.cai2012.geometry import CaiGeometry, from_config

D = 0.2


@pytest.fixture
def geom(tmp_path):
    return from_config(load_case(tmp_path))


# --------------------------------------------------------------------------- #
# the exit disk
# --------------------------------------------------------------------------- #

def test_radius_is_half_the_diameter(geom):
    assert geom.radius_m == pytest.approx(0.1)


def test_exact_nozzle_area_is_pi_r_squared(geom):
    """The number the staircase patch is measured against."""
    assert geom.exact_nozzle_area_m2 == pytest.approx(math.pi * 0.1 ** 2)


def test_half_domain_halves_the_exit_area(tmp_path):
    geom = from_config(load_case(tmp_path, geometry={"symmetry_mode": "half_y"}))
    assert geom.exact_nozzle_area_m2 == pytest.approx(0.5 * math.pi * 0.1 ** 2)


def test_the_exit_plane_is_the_origin_of_x_over_d(geom):
    assert geom.x_min_m == 0.0


# --------------------------------------------------------------------------- #
# membership -- the rule topoSet applies
# --------------------------------------------------------------------------- #

def test_the_disk_centre_is_on_the_nozzle(geom):
    assert bool(geom.is_on_nozzle(0.0, 0.0))


def test_a_point_just_inside_the_rim_is_on_the_nozzle(geom):
    assert bool(geom.is_on_nozzle(0.0999, 0.0))


def test_a_point_just_outside_the_rim_is_not(geom):
    assert not bool(geom.is_on_nozzle(0.1001, 0.0))


def test_the_rim_itself_is_inclusive(geom):
    """r <= R0, matching topoSet's cylinderToFace."""
    assert bool(geom.is_on_nozzle(0.1, 0.0))


def test_membership_is_radial_not_square(geom):
    """(0.08, 0.08) is inside the bounding square and outside the disk."""
    assert not bool(geom.is_on_nozzle(0.08, 0.08))
    assert bool(geom.is_on_nozzle(0.07, 0.07))


def test_membership_broadcasts(geom):
    y = np.array([0.0, 0.05, 0.2])
    z = np.zeros(3)
    assert list(geom.is_on_nozzle(y, z)) == [True, True, False]


def test_half_domain_excludes_the_negative_y_half(tmp_path):
    geom = from_config(load_case(tmp_path, geometry={"symmetry_mode": "half_y"}))
    assert bool(geom.is_on_nozzle(0.05, 0.0))
    assert not bool(geom.is_on_nozzle(-0.05, 0.0))


def test_contains_covers_the_box(geom):
    assert bool(geom.contains(1.0, 0.0, 0.0))
    assert not bool(geom.contains(-0.1, 0.0, 0.0))
    assert not bool(geom.contains(1.0, 5.0, 0.0))


# --------------------------------------------------------------------------- #
# the domain
# --------------------------------------------------------------------------- #

def test_extents_are_diameters_times_D(geom):
    """geometry.*_over_D = 10 with D = 0.2 gives a 2 m half-extent."""
    assert geom.x_max_m == pytest.approx(2.0)
    assert (geom.y_min_m, geom.y_max_m) == pytest.approx((-2.0, 2.0))
    assert (geom.z_min_m, geom.z_max_m) == pytest.approx((-2.0, 2.0))


def test_half_domain_starts_at_y_zero(tmp_path):
    geom = from_config(load_case(tmp_path, geometry={"symmetry_mode": "half_y"}))
    assert geom.y_min_m == 0.0
    assert geom.y_max_m == pytest.approx(2.0)


def test_domain_volume(geom):
    assert geom.domain_volume_m3 == pytest.approx(2.0 * 4.0 * 4.0)


def test_bounds_are_in_the_documented_order(geom):
    assert geom.bounds == pytest.approx((0.0, 2.0, -2.0, 2.0, -2.0, 2.0))


# --------------------------------------------------------------------------- #
# normalisation
# --------------------------------------------------------------------------- #

def test_over_diameter_normalises_by_D_not_R0(geom):
    """Cai plots X/D and Z/D. Using the radius would double every coordinate."""
    assert float(geom.over_diameter(0.2)) == pytest.approx(1.0)
    assert float(geom.over_diameter(2.0)) == pytest.approx(10.0)


def test_normalisation_round_trips(geom):
    values = np.array([0.0, 0.35, 1.9])
    assert geom.from_diameter(geom.over_diameter(values)) == pytest.approx(values)


def test_changing_the_diameter_rescales_the_whole_domain(tmp_path):
    """Every extent is in diameters, so D is the only length that sets a scale."""
    big = from_config(load_case(tmp_path, nozzle={"diameter_m": 0.4}))
    assert big.x_max_m == pytest.approx(4.0)
    assert big.radius_m == pytest.approx(0.2)
    assert float(big.over_diameter(big.x_max_m)) == pytest.approx(10.0)


def test_describe_reports_the_dimensional_and_normalised_extents(geom):
    lines = "\n".join(geom.describe())
    assert "nozzle diameter D" in lines
    assert "pi R0^2" in lines
    assert "10 D" in lines


def test_geometry_is_frozen(geom):
    with pytest.raises(Exception):
        geom.diameter_m = 1.0


def test_a_geometry_can_be_built_directly():
    """Not everything has to go through a config -- the mesh tests need this."""
    geom = CaiGeometry(diameter_m=1.0, x_max_m=5.0, y_min_m=-3.0, y_max_m=3.0,
                       z_min_m=-3.0, z_max_m=3.0)
    assert geom.radius_m == 0.5
    assert geom.domain_volume_m3 == pytest.approx(5.0 * 6.0 * 6.0)
