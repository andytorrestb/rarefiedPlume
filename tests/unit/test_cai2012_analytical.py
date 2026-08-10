"""Cai's collisionless circular-exit solution.

The closed forms in :mod:`plumetools.cai2012.analytical` are re-derived rather
than transcribed from the paper, so they are checked against every limit that can
be established independently:

* the exit plane, where the answer is the outgoing half of a Maxwellian and is
  known exactly for any speed ratio;
* the stationary case ``S0 = 0``, where the half-Maxwellian moments are textbook
  (``n/n0 = 1/2``, ``<|v_x|> = sqrt(2kT/pi m)``, ``T/T0 = 1 - 2/(3 pi)``);
* the disk quadrature, which is the same solution evaluated a completely
  different way and must agree on the axis to machine precision;
* monotonicity, symmetry and the far field.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from plumetools.cai2012 import analytical

R0 = 0.1     # nozzle radius [m]
D = 0.2      # nozzle diameter [m]
S0 = 2.0     # [PAPER] Cai's speed ratio


# --------------------------------------------------------------------------- #
# t0, the cone the exit disk subtends
# --------------------------------------------------------------------------- #

def test_cos_theta_max_is_zero_in_the_exit_plane():
    """At x = 0 the disk fills a hemisphere: theta_max = 90 deg."""
    assert analytical.cos_theta_max(0.0, R0) == pytest.approx(0.0)


def test_cos_theta_max_approaches_one_far_downstream():
    assert analytical.cos_theta_max(1.0e6, R0) == pytest.approx(1.0, abs=1e-12)


def test_cos_theta_max_is_monotone():
    x = np.linspace(0.0, 10.0 * D, 50)
    t0 = analytical.cos_theta_max(x, R0)
    assert np.all(np.diff(t0) > 0.0)


def test_upstream_x_is_rejected():
    """The solution is one-sided: there is no plume upstream of the exit plane."""
    with pytest.raises(ValueError, match="x >= 0"):
        analytical.cos_theta_max(-0.01, R0)


def test_non_positive_radius_is_rejected():
    with pytest.raises(ValueError, match="radius"):
        analytical.cos_theta_max(1.0, 0.0)


# --------------------------------------------------------------------------- #
# A -- centreline density
# --------------------------------------------------------------------------- #

def test_exit_centreline_density_is_the_outgoing_half_of_the_distribution():
    """n/n0 -> (1 + erf S0)/2 at x = 0, NOT 1.

    Only the forward-moving half of the exit distribution is present in the
    plume. The limit is 1 as S0 -> infinity and exactly 1/2 with no drift.
    """
    from scipy.special import erf
    value = float(analytical.centerline_density_ratio(0.0, R0, S0))
    assert value == pytest.approx(0.5 * (1.0 + erf(S0)))


def test_exit_centreline_density_is_within_a_quarter_percent_of_one_at_S0_2():
    """At Cai's S0 = 2 the exit density is 0.99766 n0 -- 'approaches 1'."""
    value = float(analytical.centerline_density_ratio(0.0, R0, S0))
    assert value == pytest.approx(1.0, rel=3e-3)
    assert value < 1.0


def test_exit_centreline_density_is_one_half_with_no_drift():
    assert float(analytical.centerline_density_ratio(0.0, R0, 0.0)) == (
        pytest.approx(0.5))


def test_exit_centreline_density_tends_to_one_for_a_large_speed_ratio():
    assert float(analytical.centerline_density_ratio(0.0, R0, 6.0)) == (
        pytest.approx(1.0, abs=1e-12))


def test_centreline_density_decreases_downstream():
    x = np.linspace(0.0, 10.0 * D, 200)
    n = analytical.centerline_density_ratio(x, R0, S0)
    assert np.all(np.diff(n) < 0.0)


def test_centreline_density_tends_to_zero_far_downstream():
    """1/x^2, so X/D = 1e4 is already 1e-8 and X/D = 1e5 is 1e-10.

    The closed form loses digits to cancellation out here -- both terms of (C1)
    approach (1 + erf S0)/2 -- but at X/D = 1e5 that costs about ten of the
    sixteen available, which still leaves the answer meaningful.
    """
    assert float(analytical.centerline_density_ratio(1.0e5 * D, R0, S0)) == (
        pytest.approx(0.0, abs=1e-9))


def test_far_field_density_falls_as_one_over_x_squared():
    """A point source seen from far away: doubling X quarters the density."""
    near = float(analytical.centerline_density_ratio(200.0 * D, R0, S0))
    far = float(analytical.centerline_density_ratio(400.0 * D, R0, S0))
    assert near / far == pytest.approx(4.0, rel=1e-2)


def test_centreline_density_is_bounded_by_zero_and_one():
    x = np.linspace(0.0, 10.0 * D, 200)
    n = analytical.centerline_density_ratio(x, R0, S0)
    assert np.all(n > 0.0) and np.all(n <= 1.0)


# --------------------------------------------------------------------------- #
# B -- centreline axial velocity
# --------------------------------------------------------------------------- #

def test_exit_axial_velocity_with_no_drift_is_the_half_maxwellian_mean():
    """<|v_x|> = sqrt(2kT/pi m), i.e. U sqrt(beta0) = 1/sqrt(pi)."""
    value = float(analytical.centerline_axial_speed_ratio(0.0, R0, 0.0))
    assert value == pytest.approx(1.0 / math.sqrt(math.pi))


def test_exit_axial_velocity_is_close_to_the_exit_speed_ratio():
    """At x = 0 the thermal spread barely moves the mean off S0."""
    value = float(analytical.centerline_axial_speed_ratio(0.0, R0, S0))
    assert value == pytest.approx(S0, rel=5e-3)


def test_axial_velocity_increases_downstream():
    """The visible cone narrows and the slow molecules leave it first, so the
    mean axial velocity RISES. This is the signature of a collisionless plume --
    no collisions are available to bring it back to the bulk speed."""
    x = np.linspace(0.0, 10.0 * D, 200)
    u = analytical.centerline_axial_speed_ratio(x, R0, S0)
    assert np.all(np.diff(u) > 0.0)


def test_axial_velocity_approaches_a_finite_limit():
    near = float(analytical.centerline_axial_speed_ratio(50.0 * D, R0, S0))
    far = float(analytical.centerline_axial_speed_ratio(5000.0 * D, R0, S0))
    assert far == pytest.approx(near, rel=1e-3)
    assert 2.0 < far < 3.0


# --------------------------------------------------------------------------- #
# C -- centreline temperature
# --------------------------------------------------------------------------- #

def test_exit_temperature_with_no_drift_is_the_half_maxwellian_value():
    """T_x = T(1 - 2/pi), T_y = T_z = T, so T/T0 = 1 - 2/(3 pi) = 0.78779."""
    value = float(analytical.centerline_temperature_ratio(0.0, R0, 0.0))
    assert value == pytest.approx(1.0 - 2.0 / (3.0 * math.pi))


def test_exit_temperature_is_close_to_T0_at_S0_2():
    value = float(analytical.centerline_temperature_ratio(0.0, R0, S0))
    assert value == pytest.approx(1.0, rel=1e-2)


def test_temperature_falls_downstream():
    x = np.linspace(0.0, 10.0 * D, 200)
    t = analytical.centerline_temperature_ratio(x, R0, S0)
    assert np.all(np.diff(t) < 0.0)


def test_temperature_stays_positive():
    """A collisionless expansion cools but never reaches zero: the axial
    velocity spread survives however far downstream the point is."""
    x = np.logspace(-3, 4, 60) * D
    t = analytical.centerline_temperature_ratio(x, R0, S0)
    assert np.all(t > 0.0)


# --------------------------------------------------------------------------- #
# the combined profile
# --------------------------------------------------------------------------- #

def test_profile_agrees_with_the_individual_functions():
    x = np.linspace(0.0, 10.0 * D, 25)
    profile = analytical.centerline_profile(x, R0, S0)
    assert profile["n_over_n0"] == pytest.approx(
        analytical.centerline_density_ratio(x, R0, S0))
    assert profile["U_sqrt_beta0"] == pytest.approx(
        analytical.centerline_axial_speed_ratio(x, R0, S0))
    assert profile["T_over_T0"] == pytest.approx(
        analytical.centerline_temperature_ratio(x, R0, S0))


def test_profile_normalises_x_by_the_diameter():
    profile = analytical.centerline_profile(np.array([D, 2.0 * D]), R0, S0)
    assert profile["x_over_D"] == pytest.approx([1.0, 2.0])


# --------------------------------------------------------------------------- #
# D -- the full field
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("x_over_d", [0.1, 0.5, 1.0, 5.0, 10.0])
def test_disk_quadrature_reproduces_the_closed_form_on_the_axis(x_over_d):
    """Two independent evaluations of the same solution: Cai Eq. 5 integrated
    numerically over the disk, against the closed form (C1). They must agree to
    quadrature precision or one of them is wrong."""
    x = x_over_d * D
    quadrature = float(analytical.density_ratio(x, 0.0, 0.0, R0, S0))
    closed = float(analytical.centerline_density_ratio(x, R0, S0))
    assert quadrature == pytest.approx(closed, rel=1e-10)


def test_field_is_symmetric_about_the_plume_axis():
    x = 2.0 * D
    values = [float(analytical.density_ratio(x, y, z, R0, S0))
              for y, z in ((0.3, 0.0), (-0.3, 0.0), (0.0, 0.3), (0.0, -0.3))]
    assert values[0] == pytest.approx(values[1])
    assert values[0] == pytest.approx(values[2])
    assert values[0] == pytest.approx(values[3])


def test_field_depends_only_on_the_radius_off_axis():
    """Axisymmetry: two points at the same (x, r) but different azimuth agree."""
    x, r = 2.0 * D, 0.25
    a = float(analytical.density_ratio(x, r, 0.0, R0, S0))
    b = float(analytical.density_ratio(
        x, r * math.cos(0.7), r * math.sin(0.7), R0, S0))
    assert a == pytest.approx(b, rel=1e-10)


def test_field_decreases_away_from_the_axis():
    x = 2.0 * D
    radii = np.array([0.0, 0.2, 0.4, 0.8, 1.6])
    values = np.array([float(analytical.density_ratio(x, r, 0.0, R0, S0))
                       for r in radii])
    assert np.all(np.diff(values) < 0.0)


def test_field_broadcasts_over_a_plane():
    x = np.linspace(0.05, 2.0, 7)[:, None]
    z = np.linspace(-0.5, 0.5, 5)[None, :]
    field = analytical.density_ratio(x, 0.0, z, R0, S0)
    assert field.shape == (7, 5)
    assert np.all(field > 0.0)


def test_field_in_the_exit_plane_is_rejected():
    """The disk integrand is singular at the field point itself when x = 0."""
    with pytest.raises(ValueError, match="x > 0"):
        analytical.density_ratio(0.0, 0.0, 0.0, R0, S0)


def test_quadrature_is_converged_at_the_default_order():
    x, y, z = 0.5 * D, 0.15, 0.0
    coarse = float(analytical.density_ratio(x, y, z, R0, S0,
                                            n_radial=32, n_azimuthal=64))
    fine = float(analytical.density_ratio(x, y, z, R0, S0,
                                          n_radial=128, n_azimuthal=256))
    assert coarse == pytest.approx(fine, rel=1e-8)


@pytest.mark.parametrize("n_radial,n_azimuthal", [(1, 64), (16, 2)])
def test_degenerate_quadrature_is_rejected(n_radial, n_azimuthal):
    with pytest.raises(ValueError, match="quadrature"):
        analytical.density_ratio(0.5, 0.0, 0.0, R0, S0,
                                 n_radial=n_radial, n_azimuthal=n_azimuthal)


# --------------------------------------------------------------------------- #
# the full moments -- density, velocity and temperature off the axis
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("x_over_d", [0.05, 0.5, 1.0, 5.0, 10.0])
def test_moments_reproduce_all_three_closed_forms_on_the_axis(x_over_d):
    """The disk quadrature of the second, third and fourth velocity moments,
    against (C1), (C2) and (C4). Three independent checks of one derivation."""
    x = x_over_d * D
    field = analytical.moments(x, 0.0, 0.0, R0, S0)
    assert float(field["n_over_n0"]) == pytest.approx(
        float(analytical.centerline_density_ratio(x, R0, S0)), rel=1e-8)
    assert float(field["U_sqrt_beta0"][..., 0]) == pytest.approx(
        float(analytical.centerline_axial_speed_ratio(x, R0, S0)), rel=1e-8)
    assert float(field["T_over_T0"]) == pytest.approx(
        float(analytical.centerline_temperature_ratio(x, R0, S0)), rel=1e-8)


def test_moments_density_agrees_with_density_ratio_off_axis():
    a = float(analytical.moments(0.4, 0.15, 0.05, R0, S0)["n_over_n0"])
    b = float(analytical.density_ratio(0.4, 0.15, 0.05, R0, S0))
    assert a == pytest.approx(b, rel=1e-8)


def test_the_transverse_velocity_vanishes_on_the_axis():
    field = analytical.moments(0.4, 0.0, 0.0, R0, S0)
    assert float(field["U_sqrt_beta0"][..., 1]) == pytest.approx(0.0, abs=1e-12)
    assert float(field["U_sqrt_beta0"][..., 2]) == pytest.approx(0.0, abs=1e-12)


def test_the_velocity_points_away_from_the_exit_off_axis():
    """Collisionless streaming: the mean velocity leans outward, and its
    transverse component must be dropped from nothing -- T is defined about the
    LOCAL mean, so ignoring it would overstate the temperature."""
    field = analytical.moments(0.4, 0.3, 0.0, R0, S0)
    assert float(field["U_sqrt_beta0"][..., 1]) > 0.0
    field_negative = analytical.moments(0.4, -0.3, 0.0, R0, S0)
    assert float(field_negative["U_sqrt_beta0"][..., 1]) < 0.0


def test_moments_broadcast_and_carry_a_vector_axis():
    x = np.linspace(0.05, 1.0, 6)
    field = analytical.moments(x, 0.05, 0.0, R0, S0)
    assert field["n_over_n0"].shape == (6,)
    assert field["U_sqrt_beta0"].shape == (6, 3)
    assert field["T_over_T0"].shape == (6,)


def test_moments_in_the_exit_plane_are_rejected():
    with pytest.raises(ValueError, match="x > 0"):
        analytical.moments(0.0, 0.0, 0.0, R0, S0)


# --------------------------------------------------------------------------- #
# quadrature refinement near the exit plane
# --------------------------------------------------------------------------- #

def test_the_quadrature_order_rises_as_the_exit_plane_is_approached():
    far = analytical._refined_order(R0, 1.0, 48, 96)
    near = analytical._refined_order(R0, 0.005, 48, 96)
    assert near[0] > far[0] and near[1] > far[1]
    assert far == (48, 96)      # never lowered


def test_the_quadrature_order_is_capped():
    order = analytical._refined_order(R0, 1e-12, 48, 96)
    assert order == (analytical.MAX_RADIAL, analytical.MAX_AZIMUTHAL)


def test_a_degenerate_requested_order_is_still_rejected():
    """Refinement must not hide a nonsensical request behind a bump."""
    with pytest.raises(ValueError, match="quadrature"):
        analytical._refined_order(R0, 1.0, 1, 96)


def test_close_to_the_exit_plane_the_density_stays_under_the_physical_bound():
    """Half the exit distribution is the ceiling anywhere in the plume.

    Without the refinement a fixed 64 x 128 rule returns 0.99799 at
    (x, r) = (0.05, 0.21) R0 -- above the bound, which is how an under-resolved
    quadrature announces itself.
    """
    from scipy.special import erf
    bound = 0.5 * (1.0 + erf(S0))
    for x in (0.002, 0.005, 0.01, 0.02):
        for r in (0.0, 0.015, 0.0212, 0.05, 0.099):
            assert float(analytical.density_ratio(x, r, 0.0, R0, S0)) <= bound
            assert float(analytical.moments(x, r, 0.0, R0, S0)["n_over_n0"]) <= bound


# --------------------------------------------------------------------------- #
# numerical robustness
# --------------------------------------------------------------------------- #

def test_large_speed_ratio_does_not_overflow():
    """exp(S0^2) alone overflows around S0 = 27; the moments are written so the
    exponentials combine before they are evaluated."""
    value = float(analytical.centerline_density_ratio(D, R0, 30.0))
    assert math.isfinite(value) and 0.0 < value <= 1.0


def test_everything_is_finite_across_ten_diameters():
    x = np.linspace(0.0, 10.0 * D, 500)
    profile = analytical.centerline_profile(x, R0, S0)
    for key in ("n_over_n0", "U_sqrt_beta0", "T_over_T0"):
        assert np.all(np.isfinite(profile[key])), key


def test_describe_reports_the_stations_it_is_given():
    lines = analytical.describe(R0, S0, x_over_d=(0.0, 1.0))
    assert any("X/D" in line for line in lines)
    assert len(lines) == 4  # title, header, two stations
