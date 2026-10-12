"""The geometry: the exit disk, the box on its lattice, the plate, the symmetry.

The properties protected here are the ones everything downstream leans on
without re-checking: the box is a whole number of background cells, ``y = 0``
and ``z = 0`` are lattice planes, the plate is strictly inside, and the fraction
of the flow a symmetric domain models is the fraction every total is scaled by.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from _nozzleplume import R4D_CASE, load_case
from plumetools.nozzleplume.config import load_case_config
from plumetools.nozzleplume.geometry import GeometryError, from_config


# --------------------------------------------------------------------------- #
# the shipped case
# --------------------------------------------------------------------------- #

def test_the_r4d_plate_is_the_iss_panel_normal_to_the_plume():
    """34.14 m along y, 11.88 m along z, 0.25 m along the plume axis -- the
    footprint of the `panel` patch in cases/iss_solar_panels_3."""
    plate = from_config(load_case_config(R4D_CASE)).plate
    assert (plate.y_min_m, plate.y_max_m) == (-17.07, 17.07)
    assert (plate.z_min_m, plate.z_max_m) == (-5.94, 5.94)
    assert plate.x_downstream_m - plate.x_upstream_m == pytest.approx(0.25)
    assert plate.face_area_m2 == pytest.approx(34.14 * 11.88)


def test_the_r4d_domain_is_full_three_dimensional_by_default():
    geom = from_config(load_case_config(R4D_CASE))
    assert geom.symmetry_multiplicity == 1
    assert geom.y_min_m == -geom.y_max_m and geom.z_min_m == -geom.z_max_m


# --------------------------------------------------------------------------- #
# the box sits on the background lattice
# --------------------------------------------------------------------------- #

def test_extents_that_are_whole_cells_are_left_alone(tmp_path):
    geom = from_config(load_case(tmp_path))
    assert (geom.x_max_m, geom.y_max_m, geom.z_max_m) == (6.0, 3.0, 2.0)


def test_extents_are_rounded_up_to_whole_background_cells(tmp_path):
    """Up, never down: a domain that shrank could cut the plate's clearance."""
    geom = from_config(load_case(
        tmp_path, domain={"x_max_m": 6.1, "y_half_m": 3.26, "z_half_m": 2.4}))
    assert (geom.x_max_m, geom.y_max_m, geom.z_max_m) == (6.5, 3.5, 2.5)
    assert geom.requested_extents_m == (6.1, 3.26, 2.4)


def test_the_rounding_is_reported(tmp_path):
    geom = from_config(load_case(tmp_path, domain={"x_max_m": 6.1}))
    assert any("rounded up" in line for line in geom.describe())
    assert not any("rounded up" in line
                   for line in from_config(load_case(tmp_path / "b")).describe())


def test_the_plume_axis_lies_on_lattice_planes(tmp_path):
    """y = 0 and z = 0 are a whole number of cells from the low boundary, which
    is what makes the staircase nozzle four-fold symmetric."""
    geom = from_config(load_case(tmp_path, domain={"y_half_m": 3.3, "z_half_m": 2.2}))
    for low in (geom.y_min_m, geom.z_min_m):
        assert (0.0 - low) / 0.5 == pytest.approx(round((0.0 - low) / 0.5))


# --------------------------------------------------------------------------- #
# symmetry
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("planes,multiplicity", [([], 1), (["y"], 2), (["z"], 2),
                                                 (["y", "z"], 4)])
def test_each_symmetry_plane_halves_what_is_modelled(tmp_path, planes, multiplicity):
    geom = from_config(load_case(tmp_path, domain={"symmetry_planes": planes}))
    assert geom.symmetry_multiplicity == multiplicity
    assert geom.exact_nozzle_area_m2 == pytest.approx(
        math.pi * 0.2 ** 2 / multiplicity)
    assert geom.full_nozzle_area_m2 == pytest.approx(math.pi * 0.2 ** 2)


def test_a_symmetry_plane_puts_the_low_boundary_at_zero(tmp_path):
    geom = from_config(load_case(tmp_path, domain={"symmetry_planes": ["y"]}))
    assert geom.y_min_m == 0.0
    assert geom.z_min_m == -geom.z_max_m


def test_the_modelled_plate_is_the_part_on_the_modelled_side(tmp_path):
    """The Plate is always the whole plate; the mesh holds part of it."""
    geom = from_config(load_case(tmp_path, domain={"symmetry_planes": ["y", "z"]}))
    assert geom.plate.bounds == (4.0, 4.25, -2.0, 2.0, -1.0, 1.0)
    assert geom.modelled_plate_bounds == (4.0, 4.25, 0.0, 2.0, 0.0, 1.0)


def test_only_the_modelled_part_of_the_disk_is_the_nozzle(tmp_path):
    geom = from_config(load_case(tmp_path, domain={"symmetry_planes": ["y"]}))
    y = np.array([0.1, -0.1, 0.1, 0.3])
    z = np.array([0.1, 0.1, -0.1, 0.0])
    assert geom.is_on_nozzle(y, z).tolist() == [True, False, True, False]


# --------------------------------------------------------------------------- #
# the plate
# --------------------------------------------------------------------------- #

def test_the_plate_is_placed_by_its_upstream_face(tmp_path):
    plate = from_config(load_case(tmp_path)).plate
    assert plate.x_upstream_m == 4.0
    assert plate.x_downstream_m == 4.25


def test_an_offset_moves_the_plate_off_the_plume_axis(tmp_path):
    plate = from_config(load_case(
        tmp_path, plate={"centre_y_m": 0.5, "centre_z_m": -0.25})).plate
    assert (plate.y_min_m, plate.y_max_m) == (-1.5, 2.5)
    assert (plate.z_min_m, plate.z_max_m) == (-1.25, 0.75)


def test_without_a_plate_the_jet_runs_the_length_of_the_box(tmp_path):
    geom = from_config(load_case(tmp_path, plate=None))
    assert geom.plate is None
    assert geom.modelled_plate_bounds is None
    assert geom.flow_length_m == geom.x_max_m


def test_with_a_plate_the_jet_runs_to_its_upstream_face(tmp_path):
    assert from_config(load_case(tmp_path)).flow_length_m == 4.0


@pytest.mark.parametrize("changes,edge", [
    ({"plate": {"standoff_m": 5.5}}, "lee face"),
    ({"plate": {"size_y_m": 5.5}}, "y edge"),
    ({"plate": {"size_z_m": 3.5}}, "z edge"),
    ({"plate": {"centre_y_m": -0.75}}, "y edge"),
])
def test_a_plate_within_a_cell_of_an_open_boundary_is_rejected(tmp_path, changes,
                                                               edge):
    """Flush with the boundary it would be meshed as though it ran on past it,
    and the edge the gas escapes round would be gone."""
    with pytest.raises(GeometryError, match=edge):
        from_config(load_case(tmp_path, **changes))


def test_a_plate_may_cross_a_symmetry_plane(tmp_path):
    """That is not an open boundary: the plate is meant to straddle it."""
    geom = from_config(load_case(tmp_path, domain={"symmetry_planes": ["y", "z"]}))
    assert geom.modelled_plate_bounds[2] == 0.0


def test_a_box_narrower_than_the_exit_is_rejected(tmp_path):
    with pytest.raises(GeometryError, match="exit disk"):
        from_config(load_case(tmp_path, nozzle={"exit_diameter_m": 5.0},
                              plate=None))


def test_plate_membership_grows_with_the_margin(tmp_path):
    plate = from_config(load_case(tmp_path)).plate
    assert bool(plate.contains(4.1, 0.0, 0.0))
    assert not bool(plate.contains(3.9, 0.0, 0.0))
    assert bool(plate.contains(3.9, 0.0, 0.0, margin_m=0.2))
