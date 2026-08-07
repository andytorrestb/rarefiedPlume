"""AIAA 99-3455 geometry: unit conversions, derived coordinates, clearances.

The point of this module is that no coordinate is a literal. So the tests check
the *arithmetic that produces* each coordinate, not just its value -- moving the
cylinder or the gap has to move the plate.
"""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
import pytest

from plumetools.markelov1999.constants import INCH_TO_M
from plumetools.markelov1999.geometry import (
    ASSUMED_DIMENSIONS,
    PAPER_DIMENSIONS_IN,
    GeometryError,
    MarkelovGeometry,
    cylinder_surface_angles,
    report_paper_conversions,
)


@pytest.fixture
def geom():
    return MarkelovGeometry()


# --------------------------------------------------------------------------- #
# unit conversion of the paper's dimensions
# --------------------------------------------------------------------------- #

def test_inch_conversion_is_exact_by_definition():
    assert INCH_TO_M == 0.0254


@pytest.mark.parametrize("inches,metres", [
    (6.0, 0.1524),      # cylinder diameter, plate width, gap
    (3.0, 0.0762),      # cylinder radius
    (18.0, 0.4572),     # cylinder length
    (11.75, 0.29845),   # cylinder centre x
    (15.0, 0.3810),     # plate height
])
def test_paper_inches_convert_to_the_stated_si_values(inches, metres):
    assert inches * INCH_TO_M == pytest.approx(metres, rel=1e-12)


def test_orifice_diameter_converts_to_the_paper_radius():
    """0.8255 mm diameter -> 0.00041275 m radius. Requirement 3 and 14 both
    single this out, because the legacy config has it ten times too large."""
    assert PAPER_DIMENSIONS_IN["orifice_diameter_mm"] == 0.8255
    assert 0.8255e-3 / 2.0 == pytest.approx(0.00041275, rel=1e-15)
    assert MarkelovGeometry().orifice_radius_m == pytest.approx(0.00041275, rel=1e-15)


def test_baseline_matches_every_paper_dimension(geom):
    assert geom.cylinder_radius_m == pytest.approx(0.0762, rel=1e-12)
    assert geom.cylinder_length_m == pytest.approx(0.4572, rel=1e-12)
    assert geom.cylinder_centre_x_m == pytest.approx(0.29845, rel=1e-12)
    assert geom.gap_m == pytest.approx(0.1524, rel=1e-12)
    assert geom.plate_width_y_m == pytest.approx(0.1524, rel=1e-12)
    assert geom.plate_height_z_m == pytest.approx(0.3810, rel=1e-12)


def test_the_two_assumptions_are_declared_as_such():
    """Requirement 15: paper values and modern choices must be distinguishable
    mechanically, not only in prose."""
    assert set(ASSUMED_DIMENSIONS) == {"plate_thickness_in", "inflow_radius_in"}
    assert not set(ASSUMED_DIMENSIONS) & set(PAPER_DIMENSIONS_IN)

    metadata = MarkelovGeometry().to_metadata()
    assert "plate_thickness_in" in metadata["assumed_dimensions"]
    assert "plate_thickness_in" not in metadata["paper_dimensions"]


# --------------------------------------------------------------------------- #
# derived x coordinates
# --------------------------------------------------------------------------- #

def test_cylinder_surfaces_straddle_its_centre(geom):
    assert geom.cylinder_upstream_x_m == pytest.approx(0.29845 - 0.0762, rel=1e-12)
    assert geom.cylinder_downstream_x_m == pytest.approx(0.29845 + 0.0762, rel=1e-12)
    assert geom.cylinder_upstream_x_m == pytest.approx(0.22225, rel=1e-12)
    assert geom.cylinder_downstream_x_m == pytest.approx(0.37465, rel=1e-12)


def test_plate_x_is_computed_from_cylinder_radius_and_gap(geom):
    """0.29845 + 0.0762 + 0.1524 = 0.52705 m. Requirement 3: not a literal."""
    expected = geom.cylinder_centre_x_m + geom.cylinder_radius_m + geom.gap_m
    assert geom.plate_upstream_x_m == pytest.approx(expected, rel=1e-15)
    assert geom.plate_upstream_x_m == pytest.approx(0.52705, rel=1e-12)


def test_plate_downstream_face_adds_the_thickness(geom):
    assert geom.plate_downstream_x_m == pytest.approx(
        geom.plate_upstream_x_m + geom.plate_thickness_m, rel=1e-15)
    assert geom.plate_centre_x_m == pytest.approx(
        0.5 * (geom.plate_upstream_x_m + geom.plate_downstream_x_m), rel=1e-15)


def test_moving_the_gap_moves_the_plate_and_nothing_else(geom):
    """The 12 in gap case is out of scope, but the machinery that would produce
    it must already be there -- which is what "derived, not buried" buys."""
    twelve = replace(geom, gap_m=12.0 * INCH_TO_M)
    assert twelve.plate_upstream_x_m == pytest.approx(0.37465 + 0.3048, rel=1e-12)
    assert twelve.plate_upstream_x_m - geom.plate_upstream_x_m == pytest.approx(
        0.1524, rel=1e-12)
    assert twelve.cylinder_centre_x_m == geom.cylinder_centre_x_m
    assert twelve.cylinder_downstream_x_m == geom.cylinder_downstream_x_m


def test_moving_the_cylinder_moves_the_plate(geom):
    shifted = replace(geom, cylinder_centre_x_m=geom.cylinder_centre_x_m + 0.05)
    assert shifted.plate_upstream_x_m == pytest.approx(
        geom.plate_upstream_x_m + 0.05, rel=1e-12)


def test_changing_the_plate_thickness_leaves_the_gap_alone(geom):
    """The gap is measured to the *upstream* face, so thickness must not move it."""
    thick = replace(geom, plate_thickness_m=0.0254)
    assert thick.plate_upstream_x_m == pytest.approx(geom.plate_upstream_x_m, rel=1e-15)
    assert thick.plate_downstream_x_m > geom.plate_downstream_x_m


def test_gap_is_recoverable_from_the_derived_coordinates(geom):
    """The six-inch gap, read back out of the geometry rather than out of the
    field it was stored in."""
    measured = geom.plate_upstream_x_m - geom.cylinder_downstream_x_m
    assert measured == pytest.approx(0.1524, rel=1e-12)
    assert measured / INCH_TO_M == pytest.approx(6.0, rel=1e-12)


# --------------------------------------------------------------------------- #
# extents and the symmetry plane
# --------------------------------------------------------------------------- #

def test_cylinder_spans_its_full_length_in_z(geom):
    assert geom.cylinder_z_min_m == pytest.approx(-0.2286, rel=1e-12)
    assert geom.cylinder_z_max_m == pytest.approx(+0.2286, rel=1e-12)
    assert geom.cylinder_z_max_m - geom.cylinder_z_min_m == pytest.approx(0.4572, rel=1e-12)


def test_plate_spans_its_full_height_in_z(geom):
    assert geom.plate_z_min_m == pytest.approx(-0.1905, rel=1e-12)
    assert geom.plate_z_max_m == pytest.approx(+0.1905, rel=1e-12)
    assert geom.plate_z_max_m - geom.plate_z_min_m == pytest.approx(0.3810, rel=1e-12)


def test_the_only_symmetry_plane_is_y_equals_zero(geom):
    """Requirement 3. z = 0 is NOT a symmetry plane: the cylinder ends and the
    plate edges are inside the domain, so halving z would delete the
    three-dimensional wake this case exists to study."""
    assert geom.symmetry_plane_y_m == 0.0
    # The modelled half keeps the full z extent of both bodies.
    assert geom.cylinder_z_min_m < 0.0 < geom.cylinder_z_max_m
    assert geom.plate_z_min_m < 0.0 < geom.plate_z_max_m
    # ... and half the plate width in y.
    assert geom.plate_y_max_m == pytest.approx(0.5 * geom.plate_width_y_m, rel=1e-15)


def test_cylinder_axis_endpoints_lie_on_the_axis(geom):
    lo = np.array(geom.cylinder_axis_point)
    hi = np.array(geom.cylinder_axis_end_point)
    assert lo[0] == hi[0] == pytest.approx(geom.cylinder_centre_x_m)
    assert lo[1] == hi[1] == 0.0
    # The axis is along z, so the endpoint difference has no x or y component.
    np.testing.assert_allclose((hi - lo)[:2], 0.0, atol=1e-15)
    assert np.linalg.norm(hi - lo) == pytest.approx(geom.cylinder_length_m, rel=1e-12)


# --------------------------------------------------------------------------- #
# clearances
# --------------------------------------------------------------------------- #

def test_inflow_hemisphere_clears_the_cylinder_at_baseline(geom):
    """0.22225 - 0.1524 = 0.06985 m = 2.75 in."""
    assert geom.inflow_to_cylinder_clearance_m() == pytest.approx(0.06985, rel=1e-9)
    assert geom.inflow_to_cylinder_clearance_m() > 0.0


def test_inflow_hemisphere_clears_the_plate_at_baseline(geom):
    assert geom.inflow_to_plate_clearance_m() == pytest.approx(
        0.52705 - 0.1524, rel=1e-9)


def test_an_inflow_radius_reaching_the_cylinder_is_an_error(geom):
    """Requirement 5 and 10: this must fail, not warn."""
    too_big = replace(geom, inflow_radius_m=0.25)
    assert too_big.inflow_to_cylinder_clearance_m() < 0.0
    with pytest.raises(GeometryError, match="reaches the cylinder"):
        too_big.validate()


def test_an_inflow_radius_exactly_touching_the_cylinder_is_an_error(geom):
    touching = replace(geom, inflow_radius_m=geom.cylinder_upstream_x_m)
    assert touching.inflow_to_cylinder_clearance_m() == pytest.approx(0.0, abs=1e-15)
    with pytest.raises(GeometryError, match="reaches the cylinder"):
        touching.validate()


def test_clearance_accounts_for_an_axially_offset_cylinder(geom):
    """If the cylinder is pushed far along z, the nearest point moves from its
    lateral surface to an end-cap rim and the distance grows."""
    offset = replace(geom, cylinder_centre_z_m=0.5)
    assert offset.inflow_to_cylinder_clearance_m() > geom.inflow_to_cylinder_clearance_m()


def test_baseline_geometry_validates(geom):
    geom.validate()  # must not raise


def test_orifice_must_be_inside_the_inflow_surface(geom):
    bad = replace(geom, orifice_radius_m=0.2)
    with pytest.raises(GeometryError, match="smaller than"):
        bad.validate()


@pytest.mark.parametrize("field", [
    "inflow_radius_m", "cylinder_radius_m", "cylinder_length_m", "gap_m",
    "plate_width_y_m", "plate_height_z_m", "plate_thickness_m",
])
def test_every_length_must_be_positive_and_finite(geom, field):
    with pytest.raises(GeometryError, match=f"{field} must be finite and positive"):
        replace(geom, **{field: 0.0}).validate()
    with pytest.raises(GeometryError, match=f"{field} must be finite and positive"):
        replace(geom, **{field: float("nan")}).validate()


def test_a_cylinder_straddling_the_source_plane_is_an_error(geom):
    bad = replace(geom, cylinder_centre_x_m=0.05, inflow_radius_m=0.001)
    with pytest.raises(GeometryError, match="straddles the source plane"):
        bad.validate()


# --------------------------------------------------------------------------- #
# cylinder azimuth -- windward vs leeward
# --------------------------------------------------------------------------- #

def test_windward_is_at_pi_and_leeward_at_zero(geom):
    """The plume travels +x, so the source-facing generator is the one at -x
    relative to the cylinder centre. Swapping these silently swaps the reported
    pressures, which is the whole result."""
    r = geom.cylinder_radius_m
    xc = geom.cylinder_centre_x_m
    points = np.array([
        [xc - r, 0.0, 0.0],   # windward / stagnation
        [xc + r, 0.0, 0.0],   # leeward / base
        [xc, r, 0.0],         # side, +y
    ])
    psi = cylinder_surface_angles(points, geom)
    assert abs(psi[0]) == pytest.approx(math.pi, rel=1e-12)
    assert psi[1] == pytest.approx(0.0, abs=1e-12)
    assert psi[2] == pytest.approx(math.pi / 2, rel=1e-12)


def test_cylinder_azimuth_is_independent_of_z(geom):
    r, xc = geom.cylinder_radius_m, geom.cylinder_centre_x_m
    at_mid = cylinder_surface_angles(np.array([[xc - r, 0.0, 0.0]]), geom)
    at_end = cylinder_surface_angles(np.array([[xc - r, 0.0, 0.22]]), geom)
    assert at_mid[0] == pytest.approx(at_end[0])


def test_cylinder_azimuth_rejects_bad_shapes(geom):
    with pytest.raises(ValueError, match=r"\(n, 3\)"):
        cylinder_surface_angles(np.zeros(3), geom)


# --------------------------------------------------------------------------- #
# reporting
# --------------------------------------------------------------------------- #

def test_describe_shows_the_arithmetic_not_just_the_answer(geom):
    text = "\n".join(geom.describe())
    assert "0.527050" in text                     # the derived plate face
    assert "= 0.374650 + 0.152400" in text        # and where it came from
    assert "[ASSUMPTION]" in text
    assert "[PAPER]" in text


def test_report_paper_conversions_covers_every_dimension():
    lines = report_paper_conversions()
    text = "\n".join(lines)
    assert "0.8255 mm" in text
    for key in PAPER_DIMENSIONS_IN:
        if key.endswith("_in"):
            assert key[:-3] in text
    assert text.count("[ASSUMPTION") == len(ASSUMED_DIMENSIONS)


def test_metadata_is_json_serialisable(geom):
    import json
    json.dumps(geom.to_metadata())  # must not raise


def test_metadata_records_the_gap_and_the_clearance(geom):
    si = geom.to_metadata()["si"]
    assert si["gap_m"] == pytest.approx(0.1524, rel=1e-12)
    assert si["plate_upstream_x_m"] - si["cylinder_downstream_x_m"] == pytest.approx(
        0.1524, rel=1e-9)
    assert si["inflow_to_cylinder_clearance_m"] > 0.0
