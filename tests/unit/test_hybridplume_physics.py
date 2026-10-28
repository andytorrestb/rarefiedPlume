"""hybridplume: the physical state and its transport."""

from __future__ import annotations

import math

import pytest

from plumetools.hybridplume.physics import BIRD_ARGON, K_BOLTZMANN, Gas, SourceState


@pytest.fixture
def source():
    return SourceState(BIRD_ARGON, diameter_m=0.2, speed_ratio=2.0, temperature_K=300.0,
                       knudsen=0.01)


# --------------------------------------------------------------------------- #
# state
# --------------------------------------------------------------------------- #

def test_state_is_self_consistent(source):
    g = source.gas
    assert source.pressure == pytest.approx(source.number_density * K_BOLTZMANN * 300.0, rel=1e-14)
    assert source.mass_density == pytest.approx(g.mass_kg * source.number_density, rel=1e-14)
    assert source.velocity == pytest.approx(2.0 * math.sqrt(2 * g.gas_constant * 300), rel=1e-14)
    # ae is NOT the sound speed: Me = S0 sqrt(2/gamma)
    assert source.mach == pytest.approx(source.velocity / source.sound_speed, rel=1e-12)
    assert source.mach == pytest.approx(2.0 * math.sqrt(2.0 / (5.0 / 3.0)), rel=1e-14)
    assert source.area_m2 == pytest.approx(math.pi * 0.01, rel=1e-15)
    # The density reproduces the requested Knudsen number under its own convention.
    lam = float(g.mean_free_path(source.number_density, 300.0, "vhs"))
    assert lam / source.diameter_m == pytest.approx(0.01, rel=1e-12)


def test_density_scales_inversely_with_knudsen(source):
    assert source.with_knudsen(1e-3).number_density == pytest.approx(
        10.0 * source.number_density, rel=1e-12)


def test_vhs_viscosity_reproduces_bird_table_a1():
    # Bird (1994) Table A1, argon: d = 4.17e-10 m, omega = 0.81, mu(273 K) = 2.117e-5.
    # Independent of the code: the table, not a recomputation. 0.2 % covers the
    # table's 4-figure rounding and m = 6.63e-26 vs Bird's 66.3e-27.
    assert BIRD_ARGON.mu_ref == pytest.approx(2.117e-5, rel=2e-3)


def test_mean_free_path_conventions_ratio_is_the_closed_form():
    w = BIRD_ARGON.omega
    ratios = BIRD_ARGON.mean_free_path_conventions(BIRD_ARGON.t_ref_K)
    assert ratios["hard_sphere"] == pytest.approx(1.0, rel=1e-14)   # equal at T_ref
    assert ratios["viscosity_simple"] == pytest.approx(
        15 * math.pi / (2 * (5 - 2 * w) * (7 - 2 * w)), rel=1e-12)


def test_gas_rejects_non_vhs_exponent():
    with pytest.raises(ValueError):
        Gas(omega=1.2)
