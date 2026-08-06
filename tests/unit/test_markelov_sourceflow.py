"""The paper-faithful source-flow model, equation by equation.

Every reference value here was computed **independently** of the implementation:
by hand from the printed formula, or by an alternative numerical route (adaptive
quadrature against trapezoid, a closed-form integral against the code's). A test
that only re-evaluates the same expression proves nothing, so none does.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import integrate

from plumetools.markelov1999 import sourceflow as sf
from plumetools.markelov1999.constants import (
    BOLTZMANN_J_PER_K,
    PSI_TO_PA,
    molecular_mass_kg,
)

GAMMA = 1.4
T0 = 300.0
M_N2 = molecular_mass_kg(28.0134)


# --------------------------------------------------------------------------- #
# P1  limiting angle
# --------------------------------------------------------------------------- #

def test_limiting_angle_matches_hand_evaluation():
    """theta_L = (pi/2)(sqrt(2.4/0.4) - 1) = (pi/2)(sqrt(6) - 1).

    sqrt(6) = 2.449489742783178, so theta_L = 1.5707963267948966 * 1.449489742783178.
    """
    expected = 1.5707963267948966 * (2.449489742783178 - 1.0)
    assert sf.limiting_angle(GAMMA) == pytest.approx(expected, rel=1e-15)
    assert math.degrees(sf.limiting_angle(GAMMA)) == pytest.approx(130.4540, abs=1e-4)


def test_limiting_angle_exceeds_ninety_degrees_exactly_when_gamma_is_below_five_thirds():
    """Which is why a hemispherical source never reaches the clip for N2. Stated
    as a property so nobody later "simplifies" the clipping away as unreachable:
    at gamma = 5/3 the cone closes to exactly 90 degrees, and a full-sphere source
    boundary would cross it."""
    for gamma in (1.2, 1.3, 1.4):
        assert sf.limiting_angle(gamma) > math.pi / 2
    assert sf.limiting_angle(5.0 / 3.0) == pytest.approx(math.pi / 2, rel=1e-15)
    assert sf.limiting_angle(2.0) < math.pi / 2


def test_limiting_angle_rejects_gamma_at_or_below_one():
    for bad in (1.0, 0.9, -1.0):
        with pytest.raises(ValueError, match="gamma must exceed 1"):
            sf.limiting_angle(bad)


# --------------------------------------------------------------------------- #
# P2  limiting velocity
# --------------------------------------------------------------------------- #

def test_limiting_velocity_matches_hand_evaluation():
    """V = sqrt(2*1.4*k*300 / (0.4*m)).

    Numerator 2*1.4*1.380649e-23*300 = 1.1597451600e-20 J.
    Denominator 0.4 * 4.6517346...e-26 = 1.86069...e-26 kg.
    """
    numerator = 2.0 * 1.4 * 1.380649e-23 * 300.0
    denominator = 0.4 * M_N2
    assert sf.limiting_velocity(GAMMA, T0, M_N2) == pytest.approx(
        math.sqrt(numerator / denominator), rel=1e-15)
    assert sf.limiting_velocity(GAMMA, T0, M_N2) == pytest.approx(789.4849, abs=1e-4)


def test_limiting_velocity_scales_as_sqrt_temperature():
    v300 = sf.limiting_velocity(GAMMA, 300.0, M_N2)
    v1200 = sf.limiting_velocity(GAMMA, 1200.0, M_N2)
    assert v1200 / v300 == pytest.approx(2.0, rel=1e-12)


def test_limiting_velocity_rejects_nonphysical_inputs():
    with pytest.raises(ValueError, match="T0_K"):
        sf.limiting_velocity(GAMMA, 0.0, M_N2)
    with pytest.raises(ValueError, match="mass_kg"):
        sf.limiting_velocity(GAMMA, T0, -1.0)


def test_molecular_mass_of_nitrogen():
    """m = 0.0280134 / 6.02214076e23 kg. The single conversion the package uses."""
    assert M_N2 == pytest.approx(4.651735e-26, rel=1e-6)


# --------------------------------------------------------------------------- #
# P3  angular profile
# --------------------------------------------------------------------------- #

def test_angular_exponent_is_four_point_five_two_five():
    """(1.4 + 0.41)/(1.4 - 1) = 1.81/0.4 = 4.525.

    Checked through the function rather than by reading the constant back: at
    theta = theta_L/2 the cosine argument is pi/4, so f = cos(pi/4)**4.525.
    """
    theta_l = sf.limiting_angle(GAMMA)
    expected = math.cos(math.pi / 4.0) ** 4.525
    assert float(sf.angular(GAMMA, theta_l / 2.0)) == pytest.approx(expected, rel=1e-14)


def test_angular_is_one_on_axis_and_zero_at_the_limit():
    assert float(sf.angular(GAMMA, 0.0)) == pytest.approx(1.0, rel=1e-15)
    assert float(sf.angular(GAMMA, sf.limiting_angle(GAMMA))) == 0.0


def test_angular_is_clipped_not_folded_beyond_the_limiting_angle():
    """The legacy model takes abs() of the cosine, so its profile folds and rises
    again past theta_L; the phi factor takes no abs() at all and returns NaN.
    Both are SF-6. Here it is a clean, finite zero."""
    theta_l = sf.limiting_angle(GAMMA)
    beyond = np.array([theta_l, theta_l + 0.1, math.pi, 2.0 * math.pi])
    values = sf.angular(GAMMA, beyond)
    assert np.all(np.isfinite(values)), "the legacy phi factor gives NaN here"
    assert np.all(values == 0.0)


def test_angular_is_symmetric_about_the_axis():
    theta = np.linspace(0.0, 2.0, 25)
    np.testing.assert_allclose(sf.angular(GAMMA, theta), sf.angular(GAMMA, -theta))


def test_angular_decreases_monotonically_inside_the_cone():
    theta = np.linspace(0.0, sf.limiting_angle(GAMMA), 200)
    values = sf.angular(GAMMA, theta)
    assert np.all(np.diff(values) <= 0.0)


def test_angular_takes_one_angle_not_two():
    """SF-1. The signature is the assertion: there is no phi. A separable
    f(theta)*f(phi) cannot be expressed through this API by accident."""
    import inspect
    params = list(inspect.signature(sf.angular).parameters)
    assert params == ["gamma", "theta", "exponent_offset"]


# --------------------------------------------------------------------------- #
# P4  normalisation
# --------------------------------------------------------------------------- #

def test_normalisation_integral_against_independent_quadrature():
    """Adaptive Gauss-Kronrod (the implementation) against a fine trapezoid rule
    built here from the printed formula, not from sf.angular."""
    theta_l = 1.5707963267948966 * (math.sqrt(6.0) - 1.0)
    theta = np.linspace(0.0, theta_l, 2_000_001)
    integrand = np.sin(theta) * np.cos(math.pi * theta / (2.0 * theta_l)) ** 4.525
    reference = integrate.trapezoid(integrand, theta)

    assert sf.normalization_integral(GAMMA) == pytest.approx(reference, rel=1e-9)
    assert sf.normalization_integral(GAMMA) == pytest.approx(0.357656358, rel=1e-8)


def test_normalisation_integral_trapezoid_option_agrees_with_adaptive():
    adaptive = sf.normalization_integral(GAMMA)
    trapezoid = sf.normalization_integral(GAMMA, quadrature_points=200_000)
    assert trapezoid == pytest.approx(adaptive, rel=1e-9)


def test_normalisation_coefficient_matches_hand_evaluation():
    """A = 0.5*sqrt(0.4/2.4) / 0.357656358 = 0.204124145 / 0.357656358."""
    numerator = 0.5 * math.sqrt(0.4 / 2.4)
    assert numerator == pytest.approx(0.2041241452, rel=1e-9)
    assert sf.normalization_coefficient(GAMMA) == pytest.approx(
        numerator / 0.357656358, rel=1e-8)
    assert sf.normalization_coefficient(GAMMA) == pytest.approx(0.570727014, rel=1e-8)


def test_normalisation_uses_one_gamma_throughout():
    """SF-5. The legacy API has a separate `integrand_gamma` that silently
    defaults to 1.4 whatever gamma the caller passed, so A did not normalise the
    profile it scaled. Here a different gamma must move A."""
    import inspect
    assert "integrand_gamma" not in inspect.signature(
        sf.normalization_coefficient).parameters

    a_14 = sf.normalization_coefficient(1.4)
    a_167 = sf.normalization_coefficient(5.0 / 3.0)
    assert a_14 != pytest.approx(a_167, rel=1e-6)


def test_normalisation_makes_the_profile_integrate_to_the_intended_value():
    """A is defined so that A * integral(sin f) = 0.5*sqrt((g-1)/(g+1)). Check the
    round trip, which is the property the constant exists to provide."""
    A = sf.normalization_coefficient(GAMMA)
    integral = sf.normalization_integral(GAMMA)
    assert A * integral == pytest.approx(0.5 * math.sqrt(0.4 / 2.4), rel=1e-12)


# --------------------------------------------------------------------------- #
# P5  density
# --------------------------------------------------------------------------- #

def test_sonic_density_ratio_matches_hand_evaluation():
    """(2/2.4)**(1/0.4) = 0.8333333333**2.5."""
    assert sf.sonic_density_ratio(GAMMA) == pytest.approx(
        (2.0 / 2.4) ** 2.5, rel=1e-15)
    assert sf.sonic_density_ratio(GAMMA) == pytest.approx(0.633938145, rel=1e-8)


def _rho_by_hand(p0_pa, r, theta):
    """P5 assembled here from the printed formula, with no call into sf."""
    theta_l = 1.5707963267948966 * (math.sqrt(6.0) - 1.0)
    V = math.sqrt(2.0 * 1.4 * BOLTZMANN_J_PER_K * 300.0 / (0.4 * M_N2))
    A = 0.5 * math.sqrt(0.4 / 2.4) / 0.357656358
    f = math.cos(math.pi * theta / (2.0 * theta_l)) ** 4.525 if theta < theta_l else 0.0
    return (2.0 * A * p0_pa / V ** 2) * (2.0 / 2.4) ** 2.5 * (4.1275e-4 / r) ** 2 * f


def test_mass_density_matches_an_independently_assembled_p5():
    for theta in (0.0, 0.3, 1.0, 1.5):
        got = float(sf.mass_density(theta, 0.1524, gamma=GAMMA, T0_K=T0,
                                    mass_kg=M_N2, p0_pa=475 * PSI_TO_PA))
        assert got == pytest.approx(_rho_by_hand(475 * PSI_TO_PA, 0.1524, theta), rel=1e-8)


def test_density_is_exactly_linear_in_reservoir_pressure():
    """SF-4. The legacy calculateRhoN hard-coded 475 psi and ignored the case's
    own value, so the whole archived pressure sweep varied nothing (AR-02)."""
    kwargs = dict(gamma=GAMMA, T0_K=T0, mass_kg=M_N2)
    base = float(sf.number_density(0.2, 0.1524, p0_pa=5 * PSI_TO_PA, **kwargs))
    for psi in (25, 100, 475):
        scaled = float(sf.number_density(0.2, 0.1524, p0_pa=psi * PSI_TO_PA, **kwargs))
        assert scaled / base == pytest.approx(psi / 5.0, rel=1e-12)


def test_density_falls_as_one_over_r_squared():
    kwargs = dict(gamma=GAMMA, T0_K=T0, mass_kg=M_N2, p0_pa=100 * PSI_TO_PA)
    near = float(sf.number_density(0.0, 0.1524, **kwargs))
    far = float(sf.number_density(0.0, 0.3048, **kwargs))
    assert near / far == pytest.approx(4.0, rel=1e-12)


def test_orifice_radius_defaults_to_the_paper_value():
    """SF-3. 0.8255 mm diameter -> 4.1275e-4 m radius. The legacy case.yaml
    carries 0.0041275, ten times too large, so its density is 100x too high."""
    assert sf.ORIFICE_DIAMETER_M == pytest.approx(8.255e-4, rel=1e-15)
    assert sf.ORIFICE_RADIUS_M == pytest.approx(4.1275e-4, rel=1e-15)

    kwargs = dict(gamma=GAMMA, T0_K=T0, mass_kg=M_N2, p0_pa=100 * PSI_TO_PA)
    default = float(sf.number_density(0.0, 0.1524, **kwargs))
    legacy = float(sf.number_density(0.0, 0.1524, orifice_radius_m=0.0041275, **kwargs))
    assert legacy / default == pytest.approx(100.0, rel=1e-12)


def test_density_scales_as_orifice_radius_squared():
    kwargs = dict(gamma=GAMMA, T0_K=T0, mass_kg=M_N2, p0_pa=100 * PSI_TO_PA)
    a = float(sf.number_density(0.0, 0.1524, orifice_radius_m=4.1275e-4, **kwargs))
    b = float(sf.number_density(0.0, 0.1524, orifice_radius_m=8.2550e-4, **kwargs))
    assert b / a == pytest.approx(4.0, rel=1e-12)


def test_density_is_zero_beyond_the_limiting_angle_and_never_nan():
    theta = np.array([0.0, 1.0, 2.0, sf.limiting_angle(GAMMA), 2.5, 3.0])
    n = sf.number_density(theta, 0.1524, gamma=GAMMA, T0_K=T0, mass_kg=M_N2,
                          p0_pa=100 * PSI_TO_PA)
    assert np.all(np.isfinite(n))
    assert np.all(n[theta >= sf.limiting_angle(GAMMA)] == 0.0)
    assert np.all(n[theta < sf.limiting_angle(GAMMA)] > 0.0)


def test_number_density_is_mass_density_over_molecular_mass():
    """The single documented conversion, checked as an identity."""
    kwargs = dict(gamma=GAMMA, T0_K=T0, mass_kg=M_N2, p0_pa=25 * PSI_TO_PA)
    theta = np.linspace(0.0, 1.4, 11)
    np.testing.assert_allclose(
        sf.number_density(theta, 0.2, **kwargs),
        sf.mass_density(theta, 0.2, **kwargs) / M_N2, rtol=1e-15)


def test_density_rejects_nonphysical_inputs():
    kwargs = dict(gamma=GAMMA, T0_K=T0, mass_kg=M_N2)
    with pytest.raises(ValueError, match="p0_pa"):
        sf.mass_density(0.0, 0.1, p0_pa=0.0, **kwargs)
    with pytest.raises(ValueError, match="orifice_radius_m"):
        sf.mass_density(0.0, 0.1, p0_pa=1e5, orifice_radius_m=-1.0, **kwargs)
    with pytest.raises(ValueError, match="radius_m"):
        sf.mass_density(0.0, 0.0, p0_pa=1e5, **kwargs)


def test_density_broadcasts_theta_and_radius_together():
    theta = np.linspace(0.0, 1.0, 5)
    radius = np.linspace(0.1524, 0.5, 5)
    n = sf.number_density(theta, radius, gamma=GAMMA, T0_K=T0, mass_kg=M_N2,
                          p0_pa=25 * PSI_TO_PA)
    assert n.shape == (5,)


# --------------------------------------------------------------------------- #
# angle convention and P6
# --------------------------------------------------------------------------- #

def test_off_axis_angle_is_measured_from_plus_x():
    """SF-2. The legacy convention measures theta from +z, so a point on the
    plume axis would come out at 90 degrees instead of 0."""
    points = np.array([
        [1.0, 0.0, 0.0],    # straight down the plume axis
        [0.0, 1.0, 0.0],    # 90 deg, in the symmetry plane's normal direction
        [0.0, 0.0, 1.0],    # 90 deg, along the cylinder axis
        [-1.0, 0.0, 0.0],   # straight back
    ])
    theta = sf.off_axis_angle(points)
    np.testing.assert_allclose(
        theta, [0.0, math.pi / 2, math.pi / 2, math.pi], atol=1e-15)


def test_off_axis_angle_treats_y_and_z_alike():
    """An axisymmetric model must not distinguish the two transverse directions.
    The legacy separable form does, which is the substance of SF-1."""
    for angle in (0.1, 0.5, 1.2):
        in_y = np.array([[math.cos(angle), math.sin(angle), 0.0]])
        in_z = np.array([[math.cos(angle), 0.0, math.sin(angle)]])
        assert sf.off_axis_angle(in_y)[0] == pytest.approx(sf.off_axis_angle(in_z)[0])


def test_off_axis_angle_is_independent_of_radius():
    for scale in (0.01, 1.0, 100.0):
        p = scale * np.array([[1.0, 1.0, 0.0]])
        assert sf.off_axis_angle(p)[0] == pytest.approx(math.pi / 4, rel=1e-12)


def test_off_axis_angle_survives_round_off_at_the_poles():
    """(p . a)/|p| can exceed 1 by ~1e-16; unclipped, arccos would give NaN."""
    p = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]) * (1.0 + 1e-16)
    assert np.all(np.isfinite(sf.off_axis_angle(p)))


def test_off_axis_angle_rejects_the_origin_and_bad_shapes():
    with pytest.raises(ValueError, match="undefined at the orifice centre"):
        sf.off_axis_angle(np.zeros((1, 3)))
    with pytest.raises(ValueError, match=r"\(n, 3\)"):
        sf.off_axis_angle(np.zeros((3,)))
    with pytest.raises(ValueError, match="non-zero"):
        sf.off_axis_angle(np.ones((1, 3)), plume_axis=(0.0, 0.0, 0.0))


def test_radial_velocity_has_magnitude_v_everywhere_and_points_outward():
    points = np.array([[1.0, 0.0, 0.0], [1.0, 1.0, 1.0], [0.1, -0.2, 0.3]])
    V = sf.limiting_velocity(GAMMA, T0, M_N2)
    U = sf.radial_velocity(V, points)

    np.testing.assert_allclose(np.linalg.norm(U, axis=1), V, rtol=1e-12)
    outward = np.einsum("ij,ij->i", U, points)
    assert np.all(outward > 0.0)


def test_radial_velocity_is_parallel_to_position():
    points = np.array([[0.3, -0.4, 0.5], [1.0, 2.0, -3.0]])
    U = sf.radial_velocity(500.0, points)
    cross = np.cross(U, points)
    np.testing.assert_allclose(cross, 0.0, atol=1e-10)


def test_model_name_is_explicit():
    """Requirement: a new, explicitly named variant, not a switch on the legacy
    model."""
    assert sf.MODEL_NAME == "markelov1999_axisymmetric"
