"""hybridplume: the physical state, the half-range fluxes, the continuum fluxes.

Tolerances: identities that are exact in real arithmetic are checked to 1e-12
relative (a few hundred ulps after the special functions); quadrature
comparisons to 1e-10, the integrator's own requested accuracy being 1e-13.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import quad

from plumetools.hybridplume import halfrange
from plumetools.hybridplume.fluxes import continuum_flux, viscous_stress
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


# --------------------------------------------------------------------------- #
# half-range moments
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("s", [-30.0, -8.0, -3.0, -2.0, -1.999, -0.5, 0.0, 0.7, 2.0, 6.0])
def test_half_range_moments_match_quadrature(s):
    f = lambda v: math.exp(-(v - s) ** 2) / math.sqrt(math.pi)
    j = halfrange.normal_moments(s, 1.0)
    lo = max(0.0, s - 40.0)
    for k, jk in zip((1, 2, 3), j):
        ref = quad(lambda v: v ** k * f(v), lo, max(lo, s) + 40.0, epsabs=0,
                   epsrel=1e-13, limit=400)[0]
        if ref == 0.0:          # s = -30: below double precision; must not go negative
            assert float(jk) >= 0.0
        else:
            assert float(jk) == pytest.approx(ref, rel=1e-10)


def test_reverse_flow_is_finite_positive_and_monotone():
    s = -np.linspace(2.0, 25.0, 200)
    j1, j2, j3 = halfrange.normal_moments(s, 1.0)
    for j in (j1, j2, j3):
        assert np.all(np.isfinite(j)) and np.all(j >= 0)
    nz = j1 > 0
    assert np.all(np.diff(j1[nz]) <= 0)


def test_forward_minus_reverse_is_the_inviscid_flux():
    rng = np.random.default_rng(1)
    gas = BIRD_ARGON
    k = 50
    rho = rng.uniform(1e-6, 1e-2, k)
    u = rng.normal(0, 600, (k, 3))
    t = rng.uniform(20, 900, k)
    n = rng.normal(size=(k, 3))
    n /= np.linalg.norm(n, axis=1)[:, None]
    net = halfrange.net_flux(rho, u, t, n, gas.gas_constant)
    euler = continuum_flux(rho, u, t, n, gas, inviscid=True)
    scale_m = rho * (np.linalg.norm(u, axis=1) + np.sqrt(gas.gas_constant * t))
    assert np.all(np.abs(net["mass"] - euler["mass"]) <= 1e-12 * scale_m)
    scale_p = scale_m * (np.linalg.norm(u, axis=1) + np.sqrt(gas.gas_constant * t))
    assert np.all(np.abs(net["momentum"] - euler["momentum"]) <= 1e-11 * scale_p[:, None])
    assert np.allclose(net["energy"], euler["energy"], rtol=1e-11,
                       atol=1e-12 * np.max(np.abs(euler["energy"])))


def test_solver_shaped_bird_expression_is_n_j1():
    n = np.array([1e20, 3e21])
    un = np.array([700.0, -50.0])
    a = np.array([350.0, 120.0])
    j1, _, _ = halfrange.normal_moments(un, a)
    assert np.allclose(halfrange.bird_eq_4_22_number_flux(n, un, a), n * j1, rtol=1e-13)


def test_gross_source_injection_differs_from_continuum_flux(source):
    """At S0 = 2 the half-range injection exceeds n Ue by J1/Ue - 1 = 2.4e-4."""
    gross = halfrange.gross_source_number_flux(source.number_density, 2.0,
                                               source.thermal_speed)
    ratio = gross / (source.number_density * source.velocity)
    s = 2.0
    expected = (s * 0.5 * (1 + math.erf(s)) + math.exp(-s * s) / (2 * math.sqrt(math.pi))) / s
    assert ratio == pytest.approx(expected, rel=1e-13)
    assert ratio - 1.0 == pytest.approx(2.45e-4, rel=0.02)


def test_normals_must_be_unit():
    with pytest.raises(ValueError):
        halfrange.incoming_flux([1.0], [[1, 0, 0]], [300.0], [[2.0, 0, 0]], 208.0)


# --------------------------------------------------------------------------- #
# continuum (NSF) fluxes, by manufactured states
# --------------------------------------------------------------------------- #

def test_couette_shear_and_conduction_signs():
    gas = BIRD_ARGON
    a, b = 1000.0, 50.0                      # dU_x/dy [1/s], dT/dy [K/m]
    grad_u = np.zeros((1, 3, 3))
    grad_u[0, 1, 0] = a                      # [i, j] = dU_j/dx_i
    grad_t = np.array([[0.0, b, 0.0]])
    rho, t = np.array([1e-3]), np.array([300.0])
    u = np.array([[100.0, 0.0, 0.0]])
    n = np.array([[0.0, 1.0, 0.0]])
    f = continuum_flux(rho, u, t, n, gas, grad_u=grad_u, grad_t=grad_t)
    mu = float(gas.viscosity(300.0))
    # x-momentum diffuses towards -y: its +y flux is -mu a
    assert f["viscous_momentum"][0] == pytest.approx([-mu * a, 0.0, 0.0])
    # heat flows down the gradient: +y conductive flux is -kappa b
    assert f["conductive_energy"][0] == pytest.approx(-float(gas.conductivity(300.0)) * b)
    # viscous work: -(tau.U).n = -mu a U_x
    assert f["viscous_energy"][0] == pytest.approx(-mu * a * 100.0)
    p = rho[0] * gas.gas_constant * 300.0
    assert f["momentum"][0, 1] == pytest.approx(p)


def test_viscous_stress_is_symmetric_and_traceless_for_dilatation_free_flow():
    rng = np.random.default_rng(2)
    g = rng.normal(size=(5, 3, 3))
    g -= np.trace(g, axis1=1, axis2=2)[:, None, None] * np.eye(3) / 3
    tau = viscous_stress(g, np.ones(5))
    assert np.allclose(tau, np.transpose(tau, (0, 2, 1)))
    assert np.allclose(np.trace(tau, axis1=1, axis2=2), 0.0, atol=1e-12)
