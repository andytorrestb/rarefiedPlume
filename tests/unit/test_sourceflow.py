"""The analytical source-flow model, equations E1-E8.

Reference values are the frozen ones: gamma = 1.4, T0 = 300 K, N2
(m = 4.6672e-26 kg) -- the combination the reference case actually runs, not the
argon its solver simulates (finding SM-03).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from plumetools.constants import AMU_TO_KG
from plumetools.sourceflow import (
    angular_separable,
    angular_theta_only,
    limiting_angle,
    limiting_velocity,
    normalization_coefficient,
    velocity,
)

GAMMA = 1.4
T0 = 300.0
M_N2 = 28.0134 * AMU_TO_KG
M_AR = 6.63e-26


# --------------------------------------------------------------------------- #
# E1 -- limiting angle
# --------------------------------------------------------------------------- #

def test_limiting_angle_reference_value():
    assert limiting_angle(GAMMA) == pytest.approx(2.27686, abs=1e-5)
    assert math.degrees(limiting_angle(GAMMA)) == pytest.approx(130.45, abs=0.01)


def test_limiting_angle_for_argon_is_exactly_pi_over_two():
    """At gamma = 5/3 the square root is exactly 2, so theta_l = pi/2.

    Worth pinning: it is the one value with a closed form, and it is the gamma
    the solver actually uses (SM-03).
    """
    assert limiting_angle(5 / 3) == pytest.approx(math.pi / 2, rel=1e-12)


def test_limiting_angle_decreases_with_gamma():
    values = [limiting_angle(g) for g in (1.1, 1.2, 1.4, 1.67, 2.0)]
    assert all(a > b for a, b in zip(values, values[1:]))


def test_limiting_angle_diverges_as_gamma_approaches_one():
    assert limiting_angle(1.001) > limiting_angle(1.1) > limiting_angle(1.4)


def test_limiting_angle_is_not_the_pre_22a5fdb_form():
    """Guard against reintroducing the misplaced parenthesis.

    Commit 22a5fdb corrected `sqrt((g+1)/(g-1) - 1)` to `sqrt((g+1)/(g-1)) - 1`.
    The old form gives 3.512 rad at gamma = 1.4, so results produced before
    2022-06-08 used a different limiting angle entirely.
    """
    broken = 0.5 * math.pi * math.sqrt((GAMMA + 1) / (GAMMA - 1) - 1)
    assert limiting_angle(GAMMA) != pytest.approx(broken, rel=1e-6)
    assert broken == pytest.approx(3.5124, abs=1e-4)


def test_limiting_angle_rejects_gamma_below_one():
    with pytest.raises(ValueError, match="gamma"):
        limiting_angle(1.0)


# --------------------------------------------------------------------------- #
# E2 -- limiting velocity
# --------------------------------------------------------------------------- #

def test_limiting_velocity_reference_value():
    """788.164 m/s -- measured from the captured golden, |U| on all 2044 faces."""
    assert limiting_velocity(GAMMA, T0, M_N2) == pytest.approx(788.164, abs=1e-3)


def test_argon_differs_from_nitrogen():
    """Makes SM-03 visible: the model's gas is not the solver's gas."""
    n2 = limiting_velocity(GAMMA, T0, M_N2)
    ar = limiting_velocity(5 / 3, T0, M_AR)
    assert ar == pytest.approx(558.885, abs=1e-3)
    assert n2 / ar == pytest.approx(1.410, abs=1e-3)


def test_scales_as_sqrt_temperature():
    assert (limiting_velocity(GAMMA, 4 * T0, M_N2)
            == pytest.approx(2 * limiting_velocity(GAMMA, T0, M_N2)))


def test_scales_as_inverse_sqrt_mass():
    assert (limiting_velocity(GAMMA, T0, 4 * M_N2)
            == pytest.approx(limiting_velocity(GAMMA, T0, M_N2) / 2))


def test_limiting_velocity_rejects_nonpositive_inputs():
    with pytest.raises(ValueError):
        limiting_velocity(GAMMA, -1.0, M_N2)
    with pytest.raises(ValueError):
        limiting_velocity(GAMMA, T0, 0.0)


# --------------------------------------------------------------------------- #
# E3 / E4 -- angular dependence
# --------------------------------------------------------------------------- #

def test_exponent_is_4_525_at_gamma_1_4():
    assert (GAMMA + 0.41) / (GAMMA - 1) == pytest.approx(4.525, abs=1e-3)


def test_theta_only_is_one_on_axis_and_zero_at_the_limit():
    assert angular_theta_only(GAMMA, 0.0) == pytest.approx(1.0)
    assert angular_theta_only(GAMMA, limiting_angle(GAMMA)) == pytest.approx(0.0, abs=1e-12)


def test_theta_only_decreases_monotonically():
    theta = np.linspace(0, limiting_angle(GAMMA), 50)
    f = angular_theta_only(GAMMA, theta)
    assert np.all(np.diff(f) < 0)


def test_exponent_offset_is_configurable():
    """The 0.41 is uncited (SM-07); a correction must be a config change."""
    assert (angular_theta_only(GAMMA, 0.5, exponent_offset=0.41)
            != pytest.approx(angular_theta_only(GAMMA, 0.5, exponent_offset=0.0)))


def test_separable_peaks_where_theta_is_pi_over_2_and_phi_is_zero():
    """The |theta - pi/2| shift means the peak sits on the equator, not the pole,
    because theta is measured from +z while the plume axis is +x (SM-04)."""
    peak = angular_separable(GAMMA, math.pi / 2, 0.0)
    assert peak == pytest.approx(1.0)
    assert angular_separable(GAMMA, 0.0, 0.0) < peak


def test_separable_produces_nan_beyond_the_limiting_azimuth():
    """Finding SM-06, frozen deliberately.

    f_phi has no abs(), so for |phi| > theta_l the base is negative and the
    fractional exponent yields NaN. Latent on the committed meshes (0 of 2044
    faces) but live on the archived wake-cylinder meshes (99 of 99).
    """
    theta_l = limiting_angle(GAMMA)
    with np.errstate(invalid="ignore"):
        assert np.isnan(angular_separable(GAMMA, math.pi / 2, math.pi))
    assert not np.isnan(angular_separable(GAMMA, math.pi / 2, 0.9 * theta_l))


def test_separable_theta_factor_folds_past_the_limit():
    """f_theta takes abs() of the cosine, so the profile rises again instead of
    going to zero past the limiting angle (SM-06)."""
    theta_l = limiting_angle(GAMMA)
    beyond = math.pi / 2 + 1.05 * theta_l
    assert angular_separable(GAMMA, beyond, 0.0) > 0.0


# --------------------------------------------------------------------------- #
# E5 -- normalisation
# --------------------------------------------------------------------------- #

def test_normalization_is_positive():
    assert normalization_coefficient(GAMMA, limiting_angle(GAMMA)) > 0


def test_normalization_converges_with_more_quadrature_points():
    theta_l = limiting_angle(GAMMA)
    coarse = normalization_coefficient(GAMMA, theta_l, n_points=500)
    fine = normalization_coefficient(GAMMA, theta_l, n_points=8000)
    assert coarse == pytest.approx(fine, rel=1e-4)


def test_integrand_gamma_is_a_separate_knob():
    """Documents SM-05: the original hard-coded 1.4 inside the integrand, so the
    gamma argument passed to calculateNormCoeff was silently ignored."""
    theta_l = limiting_angle(GAMMA)
    default = normalization_coefficient(GAMMA, theta_l)
    other = normalization_coefficient(GAMMA, theta_l, integrand_gamma=1.3)
    assert np.isfinite(default) and np.isfinite(other)
    assert default != pytest.approx(other, rel=1e-9)


def test_integrand_gamma_above_the_outer_gamma_yields_nan():
    """A second face of SM-05, worth knowing before anyone 'fixes' the default.

    The integration runs to theta_l(gamma), but the integrand uses
    theta_l(integrand_gamma). A larger integrand_gamma gives a smaller limiting
    angle, so the cosine goes negative inside the interval and the fractional
    exponent produces NaN.
    """
    with np.errstate(invalid="ignore"):
        result = normalization_coefficient(GAMMA, limiting_angle(GAMMA),
                                           integrand_gamma=1.67)
    assert np.isnan(result)


# --------------------------------------------------------------------------- #
# E7 -- velocity
# --------------------------------------------------------------------------- #

def test_velocity_magnitude_is_the_limiting_speed():
    theta = np.linspace(0.1, math.pi - 0.1, 20)
    phi = np.linspace(-math.pi, math.pi, 20)
    u = velocity(788.164, theta, phi)
    np.testing.assert_allclose(np.linalg.norm(u, axis=1), 788.164, rtol=1e-12)


def test_velocity_is_radial_for_consistent_angles():
    """When theta and phi come from the same point, U is exactly v_l * r_hat."""
    rng = np.random.default_rng(2)
    pts = rng.normal(size=(30, 3))
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    theta = np.arccos(pts[:, 2])
    phi = np.arctan2(pts[:, 1], pts[:, 0])
    u = velocity(1.0, theta, phi)
    np.testing.assert_allclose(u, pts, atol=1e-12)


def test_velocity_shape():
    assert velocity(1.0, np.zeros(5), np.zeros(5)).shape == (5, 3)
