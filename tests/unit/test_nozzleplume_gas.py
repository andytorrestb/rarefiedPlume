"""The VHS species and the kinetic-theory relations of the nozzle-plume family.

The three wall fluxes are the ones that matter most: the first is the injection
rate, the second the jet's thrust, the third its power, and every number the
plate post-processing is judged against is built from them. Each is checked at
both of its limits -- effusion and hypersonic -- where the answer is known
without the formula.
"""

from __future__ import annotations

import math

import pytest

from plumetools.nozzleplume import gas
from plumetools.nozzleplume.constants import BOLTZMANN_J_PER_K

EXHAUST = gas.VhsSpecies(name="exhaust", mass_kg=3.40055e-26, diameter_m=3.5e-10)
NITROGEN = gas.VhsSpecies(name="N2", mass_kg=4.65e-26, diameter_m=4.17e-10,
                          omega=0.74, t_ref_K=273.0,
                          internal_degrees_of_freedom=2)


# --------------------------------------------------------------------------- #
# the species
# --------------------------------------------------------------------------- #

def test_the_collision_cross_section_is_pi_d_squared():
    """The R-4D data sheet prints 3.84845e-19 m^2 for d = 3.5e-10 m."""
    assert EXHAUST.collision_cross_section_m2 == pytest.approx(3.84845e-19, rel=1e-5)


def test_molar_mass_follows_from_the_molecular_mass():
    assert EXHAUST.molar_mass_g_per_mol == pytest.approx(20.4786, rel=1e-5)


def test_gamma_follows_from_the_mode_count():
    assert EXHAUST.gamma == pytest.approx(5.0 / 3.0)
    assert NITROGEN.gamma == pytest.approx(1.4)


def test_the_three_thermal_speeds_are_in_their_fixed_ratios():
    """c_m : c_bar : sigma = sqrt(2) : sqrt(8/pi) : 1, whatever the gas."""
    sigma = EXHAUST.thermal_sigma(498.04)
    assert EXHAUST.most_probable_speed(498.04) == pytest.approx(math.sqrt(2.0) * sigma)
    assert EXHAUST.mean_thermal_speed(498.04) == pytest.approx(
        math.sqrt(8.0 / math.pi) * sigma)


@pytest.mark.parametrize("name", ["R-4D", "4d", "", "exhaust gas", "a/b"])
def test_a_species_name_that_is_not_an_identifier_is_rejected(name):
    """It becomes a typeId and part of a file name."""
    with pytest.raises(ValueError, match="not an identifier"):
        gas.VhsSpecies(name=name, mass_kg=1e-26, diameter_m=1e-10)


@pytest.mark.parametrize("omega", [0.49, 1.01, -1.0])
def test_omega_outside_the_vhs_range_is_rejected(omega):
    with pytest.raises(ValueError, match="omega"):
        gas.VhsSpecies(name="x", mass_kg=1e-26, diameter_m=1e-10, omega=omega)


@pytest.mark.parametrize("field", ["mass_kg", "diameter_m", "t_ref_K"])
def test_a_non_positive_property_is_rejected(field):
    values = {"name": "x", "mass_kg": 1e-26, "diameter_m": 1e-10, "t_ref_K": 273.0}
    values[field] = 0.0
    with pytest.raises(ValueError, match=field):
        gas.VhsSpecies(**values)


# --------------------------------------------------------------------------- #
# the mean free path
# --------------------------------------------------------------------------- #

def test_a_hard_sphere_mean_free_path_is_the_textbook_expression():
    n = 1.306e22
    expected = 1.0 / (math.sqrt(2.0) * math.pi * 3.5e-10 ** 2 * n)
    assert gas.mean_free_path(n, EXHAUST, 498.04) == pytest.approx(expected)


def test_at_omega_one_half_the_two_conventions_are_the_same_number():
    """Which is why a data sheet with one diameter and no T_ref is a hard sphere."""
    vhs = gas.mean_free_path(1e22, EXHAUST, 498.04, convention="vhs")
    hard = gas.mean_free_path(1e22, EXHAUST, 498.04, convention="hard_sphere")
    assert vhs == hard


def test_away_from_omega_one_half_the_temperature_correction_applies():
    vhs = gas.mean_free_path(1e20, NITROGEN, 546.0, convention="vhs")
    hard = gas.mean_free_path(1e20, NITROGEN, 546.0, convention="hard_sphere")
    assert vhs / hard == pytest.approx(2.0 ** 0.24)


def test_density_and_mean_free_path_are_inverses():
    lam = gas.mean_free_path(3.3e21, NITROGEN, 410.0)
    assert gas.number_density_from_mean_free_path(lam, NITROGEN, 410.0) == (
        pytest.approx(3.3e21))


def test_a_vacuum_has_an_infinite_mean_free_path_and_collision_time():
    assert gas.mean_free_path(0.0, EXHAUST, 300.0) == float("inf")
    assert gas.collision_time(0.0, EXHAUST, 300.0) == float("inf")


def test_an_unknown_convention_is_rejected():
    with pytest.raises(ValueError, match="convention"):
        gas.mean_free_path(1e20, EXHAUST, 300.0, convention="bird")


# --------------------------------------------------------------------------- #
# fluxes through a plane
# --------------------------------------------------------------------------- #

def test_at_rest_the_number_flux_is_the_effusion_flux():
    n, t = 1e21, 400.0
    assert gas.maxwellian_number_flux(n, 0.0, EXHAUST, t) == pytest.approx(
        0.25 * n * EXHAUST.mean_thermal_speed(t))


def test_at_rest_the_momentum_flux_is_half_the_pressure():
    """One side of the plane supplies half of ``n k T``."""
    n, t = 1e21, 400.0
    assert gas.maxwellian_momentum_flux(n, 0.0, EXHAUST, t) == pytest.approx(
        0.5 * n * BOLTZMANN_J_PER_K * t)


def test_at_rest_each_effusing_molecule_carries_two_kT():
    n, t = 1e21, 400.0
    flux = gas.maxwellian_number_flux(n, 0.0, EXHAUST, t)
    assert gas.maxwellian_energy_flux(n, 0.0, EXHAUST, t) == pytest.approx(
        flux * 2.0 * BOLTZMANN_J_PER_K * t)


def test_internal_modes_add_their_equipartition_energy_to_the_flux():
    n, t = 1e21, 400.0
    monatomic = gas.VhsSpecies(name="x", mass_kg=NITROGEN.mass_kg,
                               diameter_m=NITROGEN.diameter_m)
    extra = (gas.maxwellian_energy_flux(n, 3.0, NITROGEN, t)
             - gas.maxwellian_energy_flux(n, 3.0, monatomic, t))
    flux = gas.maxwellian_number_flux(n, 3.0, NITROGEN, t)
    assert extra == pytest.approx(flux * BOLTZMANN_J_PER_K * t)   # (2/2) k T each


@pytest.mark.parametrize("speed_ratio", [5.0, 8.0])
def test_at_high_speed_ratio_the_fluxes_are_the_continuum_ones(speed_ratio):
    """n U, rho (U^2 + R T) and rho U (U^2/2 + 5/2 R T)."""
    n, t = 1.306e22, 498.04
    r = EXHAUST.gas_constant_j_per_kg_k
    u = speed_ratio * EXHAUST.most_probable_speed(t)
    rho = n * EXHAUST.mass_kg

    assert gas.maxwellian_number_flux(n, speed_ratio, EXHAUST, t) == (
        pytest.approx(n * u, rel=1e-9))
    assert gas.maxwellian_momentum_flux(n, speed_ratio, EXHAUST, t) == (
        pytest.approx(rho * (u * u + r * t), rel=1e-9))
    assert gas.maxwellian_energy_flux(n, speed_ratio, EXHAUST, t) == (
        pytest.approx(rho * u * (0.5 * u * u + 2.5 * r * t), rel=1e-9))


def test_a_negative_density_is_rejected_by_every_flux():
    for flux in (gas.maxwellian_number_flux, gas.maxwellian_momentum_flux,
                 gas.maxwellian_energy_flux):
        with pytest.raises(ValueError, match="number density"):
            flux(-1.0, 1.0, EXHAUST, 300.0)
