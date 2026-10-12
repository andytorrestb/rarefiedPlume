"""The exit state: one object, however ``case.yaml`` chose to state it.

Three kinds of check:

* the shipped R-4D case reproduces **the rest of its own data sheet** -- the
  exit area, the mass density and the collision cross-section were supplied
  alongside the inputs, and are derived here rather than typed in, so each is a
  free cross-check of the others;
* the alternative forms of the density and the speed are the same state;
* given Cai & Wang's inputs, this package and ``plumetools.cai2012`` produce the
  same exit state. That is the regression link between the general family and
  the validated one it was generalised from.
"""

from __future__ import annotations

import math

import pytest

from _cai2012 import load_case as load_cai_case
from _nozzleplume import R4D_CASE, load_case
from plumetools.cai2012 import inflow as cai_inflow
from plumetools.nozzleplume import exitstate
from plumetools.nozzleplume.config import load_case_config
from plumetools.nozzleplume.constants import BOLTZMANN_J_PER_K


@pytest.fixture(scope="module")
def r4d():
    return exitstate.from_config(load_case_config(R4D_CASE))


# --------------------------------------------------------------------------- #
# the R-4D data sheet
# --------------------------------------------------------------------------- #

def test_r4d_exit_area_matches_the_data_sheet(r4d):
    assert r4d.area_m2 == pytest.approx(0.111590, abs=5e-7)


def test_r4d_mass_density_matches_the_data_sheet(r4d):
    """4.442e-4 kg/m^3 is quoted to four figures; n m is 4.4411e-4."""
    assert r4d.mass_density_kg_per_m3 == pytest.approx(4.442e-4, rel=3e-4)


def test_r4d_collision_cross_section_matches_the_data_sheet(r4d):
    assert r4d.species.collision_cross_section_m2 == pytest.approx(
        3.84845e-19, rel=1e-5)


def test_r4d_exit_knudsen_number_is_of_order_ten_to_the_minus_four(r4d):
    """The reason the case is expensive, and the reason the weight is pinned."""
    assert r4d.knudsen == pytest.approx(3.732e-4, rel=1e-3)
    assert r4d.mean_free_path_m == pytest.approx(1.4069e-4, rel=1e-4)


def test_r4d_is_hypersonic_at_the_exit(r4d):
    assert r4d.speed_ratio == pytest.approx(5.012, rel=1e-3)
    assert r4d.mach_number > 5.0


def test_r4d_thrust_and_mass_flow_are_those_of_a_500_newton_thruster(r4d):
    """A sanity check no data sheet row gives directly: an R-4D is a ~490 N
    engine burning ~0.16 kg/s, and the supplied exit state implies both."""
    assert r4d.mass_flow_kg_per_s() == pytest.approx(0.15795, rel=1e-4)
    assert r4d.thrust_n() == pytest.approx(513.4, rel=1e-3)
    assert 300.0 < r4d.specific_impulse_s < 350.0


def test_the_report_marks_which_values_were_given(r4d):
    lines = {line.split("  ")[1].strip(): line for line in r4d.describe()[1:]}
    assert lines["number density n"].endswith("[GIVEN]")
    assert lines["bulk speed U"].endswith("[GIVEN]")
    assert lines["mass density rho"].endswith("[DERIVED]")
    assert lines["speed ratio S"].endswith("[DERIVED]")
    assert lines["Knudsen number"].endswith("[DERIVED]")


# --------------------------------------------------------------------------- #
# derived quantities
# --------------------------------------------------------------------------- #

def test_thrust_is_momentum_plus_pressure_thrust(tmp_path):
    state = exitstate.from_config(load_case(tmp_path))
    rho, u = state.mass_density_kg_per_m3, state.speed_m_per_s
    expected = (rho * u * u + state.static_pressure_pa) * state.area_m2
    assert state.thrust_n() == pytest.approx(expected, rel=1e-9)


def test_static_pressure_is_n_k_T(tmp_path):
    state = exitstate.from_config(load_case(tmp_path))
    assert state.static_pressure_pa == pytest.approx(
        1.0e22 * BOLTZMANN_J_PER_K * 500.0)


def test_injection_through_a_given_area_scales_with_it(tmp_path):
    """The meshed patch is not the disk; the rate must follow the patch."""
    state = exitstate.from_config(load_case(tmp_path))
    assert state.injection_rate_per_s(0.5 * state.area_m2) == pytest.approx(
        0.5 * state.injection_rate_per_s())


def test_specific_impulse_does_not_depend_on_the_area(tmp_path):
    state = exitstate.from_config(load_case(tmp_path))
    assert state.thrust_n(0.01) / (state.mass_flow_kg_per_s(0.01) * 9.80665) == (
        pytest.approx(state.specific_impulse_s))


def test_the_time_step_speed_is_three_sigma_above_the_bulk_speed(tmp_path):
    state = exitstate.from_config(load_case(tmp_path))
    sigma = math.sqrt(state.species.gas_constant_j_per_kg_k * 500.0)
    assert state.characteristic_speed_m_per_s == pytest.approx(3000.0 + 3.0 * sigma)


# --------------------------------------------------------------------------- #
# the alternative forms are the same state
# --------------------------------------------------------------------------- #

def test_a_mass_density_gives_the_same_state_as_its_number_density(tmp_path):
    direct = exitstate.from_config(load_case(tmp_path / "n"))
    by_mass = exitstate.from_config(load_case(tmp_path / "rho", exit={
        "number_density_per_m3": None,
        "mass_density_kg_per_m3": direct.mass_density_kg_per_m3}))
    assert by_mass.number_density_per_m3 == pytest.approx(1.0e22)
    assert by_mass.density_given_as == "mass_density_kg_per_m3"


def test_a_knudsen_number_gives_the_same_state_as_its_number_density(tmp_path):
    direct = exitstate.from_config(load_case(tmp_path / "n"))
    by_kn = exitstate.from_config(load_case(tmp_path / "kn", exit={
        "number_density_per_m3": None, "knudsen": direct.knudsen}))
    assert by_kn.number_density_per_m3 == pytest.approx(1.0e22)
    assert by_kn.knudsen == pytest.approx(direct.knudsen)


def test_a_speed_ratio_gives_the_same_state_as_its_speed(tmp_path):
    direct = exitstate.from_config(load_case(tmp_path / "u"))
    by_ratio = exitstate.from_config(load_case(tmp_path / "s", exit={
        "speed_m_per_s": None, "speed_ratio": direct.speed_ratio}))
    assert by_ratio.speed_m_per_s == pytest.approx(3000.0)
    assert by_ratio.speed_given_as == "speed_ratio"


# --------------------------------------------------------------------------- #
# agreement with the Cai family
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("knudsen", [100.0, 0.1, 0.01])
def test_cai_inputs_give_the_cai_exit_state(tmp_path, knudsen):
    """Kn, S0, T0 and argon, through both packages.

    ``plumetools.cai2012`` is the validated record; this package was
    generalised from it without importing it. If the two ever disagree on the
    state a Knudsen number implies, one of them has changed its physics.
    """
    cai = cai_inflow.from_config(
        load_cai_case(tmp_path / "cai", exit={"knudsen": knudsen}))

    general = exitstate.from_config(load_case(
        tmp_path / "general",
        nozzle={"exit_diameter_m": 0.2},
        exit={"temperature_K": 300.0, "speed_m_per_s": None, "speed_ratio": 2.0,
              "number_density_per_m3": None, "knudsen": knudsen},
        gas={"species_name": "Ar", "molecular_mass_kg": 6.63e-26,
             "molecular_diameter_m": 4.17e-10, "omega": 0.74, "t_ref_K": 273.0},
        plate=None,
        mesh={"refinement_regions": []}))

    assert general.number_density_per_m3 == pytest.approx(
        cai.number_density_per_m3, rel=1e-12)
    assert general.speed_m_per_s == pytest.approx(cai.velocity_m_per_s, rel=1e-12)
    assert general.mean_free_path_m == pytest.approx(cai.mean_free_path_m, rel=1e-12)
    assert general.number_flux_per_m2_s == pytest.approx(
        cai.number_flux_per_m2_s, rel=1e-12)
    assert general.collision_time_s == pytest.approx(cai.collision_time_s(), rel=1e-12)
