"""Argon VHS properties and the Kn / S0 conversions the Cai matrix rests on.

Every number in the Cai study comes out of this module, so these tests are the
line of defence for the whole case family: a wrong mean-free-path convention or
a dropped factor of two in the characteristic length would produce a study that
runs cleanly at the wrong densities.
"""

from __future__ import annotations

import math

import pytest

from plumetools.cai2012 import gas
from plumetools.cai2012.constants import AVOGADRO_PER_MOL, BOLTZMANN_J_PER_K
from plumetools.cai2012.gas import ARGON, VhsSpecies

D = 0.2          # nozzle diameter [m]
T0 = 300.0       # exit temperature [K]


# --------------------------------------------------------------------------- #
# the species
# --------------------------------------------------------------------------- #

def test_argon_matches_the_repository_model():
    """The values cases/3d-inflow already runs argon with."""
    assert ARGON.name == "Ar"
    assert ARGON.mass_kg == 6.63e-26
    assert ARGON.diameter_m == 4.17e-10
    assert ARGON.omega == 0.74


def test_gas_constant_is_boltzmann_over_mass():
    assert ARGON.gas_constant_j_per_kg_k == pytest.approx(
        BOLTZMANN_J_PER_K / 6.63e-26)
    # Argon: ~208 J/(kg K).
    assert ARGON.gas_constant_j_per_kg_k == pytest.approx(208.24, rel=1e-3)


def test_monatomic_gamma_is_five_thirds():
    assert ARGON.gamma == pytest.approx(5.0 / 3.0)


def test_molar_mass_is_reported_not_configured():
    """6.63e-26 kg is 39.93 g/mol -- argon, to three figures."""
    assert ARGON.molar_mass_g_per_mol == pytest.approx(39.93, rel=1e-3)
    assert ARGON.molar_mass_g_per_mol == pytest.approx(
        6.63e-26 * AVOGADRO_PER_MOL * 1000.0)


@pytest.mark.parametrize("field,value", [
    ("mass_kg", 0.0), ("mass_kg", -1.0),
    ("diameter_m", 0.0), ("t_ref_K", -273.0),
])
def test_non_positive_species_constants_are_rejected(field, value):
    with pytest.raises(ValueError, match=field):
        VhsSpecies(**{field: value})


@pytest.mark.parametrize("omega", [0.4, 1.2])
def test_omega_outside_the_vhs_range_is_rejected(omega):
    with pytest.raises(ValueError, match="omega"):
        VhsSpecies(omega=omega)


def test_most_probable_speed_and_beta_are_inverses():
    """beta0 = 1/(2 R T0), so sqrt(beta0) * sqrt(2 R T0) = 1 exactly."""
    speed = ARGON.most_probable_speed(T0)
    assert speed * math.sqrt(ARGON.beta(T0)) == pytest.approx(1.0)


def test_mean_thermal_speed_is_the_maxwellian_mean():
    """c_bar = sqrt(8 R T / pi), and c_bar / sqrt(2RT) = sqrt(4/pi)."""
    ratio = ARGON.mean_thermal_speed(T0) / ARGON.most_probable_speed(T0)
    assert ratio == pytest.approx(math.sqrt(4.0 / math.pi))


# --------------------------------------------------------------------------- #
# Kn -> lambda
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("kn,expected", [(100.0, 20.0), (0.1, 0.02), (0.01, 0.002)])
def test_mean_free_path_from_kn_uses_the_diameter(kn, expected):
    """Cai's characteristic length is the nozzle DIAMETER, not the radius.

    Using the radius would halve every lambda0 and so double every n0 -- a
    factor of two hiding inside a definition.
    """
    assert gas.mean_free_path_from_kn(kn, D) == pytest.approx(expected)


def test_kn_on_the_radius_is_half_the_diameter_answer():
    assert gas.mean_free_path_from_kn(0.1, D / 2) == pytest.approx(
        0.5 * gas.mean_free_path_from_kn(0.1, D))


@pytest.mark.parametrize("kn", [0.0, -1.0])
def test_non_positive_knudsen_is_rejected(kn):
    with pytest.raises(ValueError, match="knudsen"):
        gas.mean_free_path_from_kn(kn, D)


# --------------------------------------------------------------------------- #
# lambda <-> n
# --------------------------------------------------------------------------- #

def test_hard_sphere_mean_free_path_is_the_textbook_expression():
    n = 1.0e20
    assert gas.mean_free_path(n, ARGON, T0, convention="hard_sphere") == (
        pytest.approx(1.0 / (math.sqrt(2.0) * math.pi * ARGON.diameter_m ** 2 * n)))


def test_vhs_equals_hard_sphere_at_the_reference_temperature():
    """Equation (G1) at T = T_ref: the correction factor is exactly 1.

    This is why the repository's existing hard-sphere expression is not a
    different model, only the VHS one evaluated at its reference state.
    """
    n = 1.0e20
    at_ref = gas.mean_free_path(n, ARGON, ARGON.t_ref_K, convention="vhs")
    hard = gas.mean_free_path(n, ARGON, ARGON.t_ref_K, convention="hard_sphere")
    assert at_ref == pytest.approx(hard, rel=1e-15)


def test_vhs_correction_at_the_exit_temperature_is_the_documented_2_3_percent():
    """(300/273)^(0.74 - 0.5) = 1.0229. Small, but not round-off."""
    factor = gas.vhs_temperature_factor(ARGON, T0)
    assert factor == pytest.approx((300.0 / 273.0) ** 0.24)
    assert factor == pytest.approx(1.0229, abs=5e-4)


def test_vhs_mean_free_path_grows_with_temperature():
    """omega > 1/2 means a hotter gas has a smaller cross-section."""
    n = 1.0e20
    cold = gas.mean_free_path(n, ARGON, 200.0)
    hot = gas.mean_free_path(n, ARGON, 400.0)
    assert hot > cold


@pytest.mark.parametrize("convention", ["vhs", "hard_sphere"])
def test_mean_free_path_and_number_density_round_trip(convention):
    n = 3.21e19
    lam = gas.mean_free_path(n, ARGON, T0, convention=convention)
    back = gas.number_density_from_mean_free_path(
        lam, ARGON, T0, convention=convention)
    assert back == pytest.approx(n, rel=1e-12)


def test_zero_density_is_an_infinite_mean_free_path():
    """A vacuum cell in a diagnostic loop must not raise."""
    assert gas.mean_free_path(0.0, ARGON, T0) == math.inf


def test_negative_density_is_rejected():
    with pytest.raises(ValueError, match="finite"):
        gas.mean_free_path(-1.0, ARGON, T0)


def test_unknown_convention_is_rejected():
    with pytest.raises(ValueError, match="convention"):
        gas.mean_free_path(1e20, ARGON, T0, convention="bird1994")


# --------------------------------------------------------------------------- #
# Kn -> n0, the conversion the whole study rests on
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("kn", [100.0, 0.1, 0.01])
def test_number_density_from_kn_inverts_back_to_kn(kn):
    n0 = gas.number_density_from_kn(kn, D, ARGON, T0)
    assert gas.knudsen_from_number_density(n0, D, ARGON, T0) == pytest.approx(kn)


def test_number_density_scales_inversely_with_knudsen():
    """A factor of ten in Kn is a factor of ten in n0, exactly."""
    dense = gas.number_density_from_kn(0.01, D, ARGON, T0)
    thin = gas.number_density_from_kn(0.1, D, ARGON, T0)
    assert dense / thin == pytest.approx(10.0)


def test_the_three_cai_cases_have_three_different_densities():
    """A configuration labelled Kn = 0.1 can never run another case's density."""
    densities = [gas.number_density_from_kn(kn, D, ARGON, T0)
                 for kn in (100.0, 0.1, 0.01)]
    assert len(set(densities)) == 3
    assert densities[0] < densities[1] < densities[2]


def test_kn100_exit_density_is_the_documented_value():
    """n0 = (T0/Tref)^(omega-1/2) / (sqrt(2) pi d^2 lambda0), lambda0 = 20 m."""
    n0 = gas.number_density_from_kn(100.0, D, ARGON, T0)
    expected = ((300.0 / 273.0) ** 0.24
                / (math.sqrt(2.0) * math.pi * 4.17e-10 ** 2 * 20.0))
    assert n0 == pytest.approx(expected, rel=1e-12)
    assert n0 == pytest.approx(6.62e16, rel=1e-3)


# --------------------------------------------------------------------------- #
# S0 -> U0
# --------------------------------------------------------------------------- #

def test_speed_from_speed_ratio_is_cais_definition():
    """U0 = S0 sqrt(2 R T0). No gamma, no gas-dynamic relation."""
    u0 = gas.speed_from_speed_ratio(2.0, ARGON, T0)
    assert u0 == pytest.approx(2.0 * math.sqrt(2.0 * ARGON.gas_constant_j_per_kg_k * T0))
    assert u0 == pytest.approx(706.95, rel=1e-4)


def test_speed_ratio_round_trips():
    u0 = gas.speed_from_speed_ratio(2.0, ARGON, T0)
    assert gas.speed_ratio_from_speed(u0, ARGON, T0) == pytest.approx(2.0)


def test_speed_ratio_zero_is_a_stationary_exit():
    assert gas.speed_from_speed_ratio(0.0, ARGON, T0) == 0.0


def test_speed_scales_with_the_square_root_of_temperature():
    hot = gas.speed_from_speed_ratio(2.0, ARGON, 4.0 * T0)
    cold = gas.speed_from_speed_ratio(2.0, ARGON, T0)
    assert hot / cold == pytest.approx(2.0)


# --------------------------------------------------------------------------- #
# the injected flux
# --------------------------------------------------------------------------- #

def test_number_flux_at_zero_speed_ratio_is_the_effusion_flux():
    """Gamma = n c_bar / 4 when there is no drift."""
    n = 1.0e20
    flux = gas.maxwellian_number_flux(n, 0.0, ARGON, T0)
    assert flux == pytest.approx(0.25 * n * ARGON.mean_thermal_speed(T0))


def test_number_flux_approaches_n_u_at_a_high_speed_ratio():
    """At S0 = 2 the thermal correction is already under a tenth of a percent."""
    n = 1.0e20
    u0 = gas.speed_from_speed_ratio(2.0, ARGON, T0)
    flux = gas.maxwellian_number_flux(n, 2.0, ARGON, T0)
    assert flux == pytest.approx(n * u0, rel=2e-3)


def test_number_flux_is_monotone_in_the_speed_ratio():
    n = 1.0e20
    fluxes = [gas.maxwellian_number_flux(n, s, ARGON, T0)
              for s in (0.0, 0.5, 1.0, 2.0, 4.0)]
    assert fluxes == sorted(fluxes)


def test_number_flux_is_linear_in_density():
    a = gas.maxwellian_number_flux(1.0e20, 2.0, ARGON, T0)
    b = gas.maxwellian_number_flux(3.0e20, 2.0, ARGON, T0)
    assert b == pytest.approx(3.0 * a)


# --------------------------------------------------------------------------- #
# collision time
# --------------------------------------------------------------------------- #

def test_collision_time_is_mean_free_path_over_mean_speed():
    n = 1.0e20
    expected = gas.mean_free_path(n, ARGON, T0) / ARGON.mean_thermal_speed(T0)
    assert gas.collision_time(n, ARGON, T0) == pytest.approx(expected)


def test_reference_collision_time_at_kn_0p01_is_about_five_microseconds():
    """lambda_ref = 2 mm at ~399 m/s. This is the t0 Cai's dt/t0 = 1 refers to."""
    n_ref = gas.number_density_from_kn(0.01, D, ARGON, T0)
    assert gas.collision_time(n_ref, ARGON, T0) == pytest.approx(5.01e-6, rel=1e-2)


def test_describe_reports_every_derived_gas_quantity():
    lines = "\n".join(gas.describe(ARGON, T0))
    for token in ("Ar", "VHS", "R = k_B/m", "beta0", "VHS factor"):
        assert token in lines
