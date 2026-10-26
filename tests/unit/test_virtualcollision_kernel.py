"""Independent physical limits and invariants of the verified elastic law."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from plumetools.cai2012.constants import BOLTZMANN_J_PER_K
from plumetools.cai2012.gas import ARGON, VhsSpecies
from plumetools.virtualcollision.kernel import (
    VHSKernel, ZeroKernel, UnsupportedCollisionLaw, elastic_pair,
)


def test_v2512_independently_checked_argon_reference():
    # Independent calculation from the installed v2512 sigmaTcR expression,
    # with mR=3.315e-26 kg and solver k=1.38065e-23. Includes Gamma(1.76), not
    # an inferred cross-section chosen to match the mean free path.
    assert VHSKernel(ARGON).rate(500.) == pytest.approx(
        2.8978791142188915e-16, rel=2e-14)


def test_hard_sphere_units_and_equilibrium_collision_rate():
    # Omega=1/2 is a geometric hard sphere. Integrating over the known
    # Maxwell relative-speed PDF must give pi*d^2 * 4 sqrt(kT/(pi*m)).
    sp = VhsSpecies(omega=0.5)
    kernel = VHSKernel(sp)
    temperature = 300.
    a = sp.mass_kg / (4 * BOLTZMANN_J_PER_K * temperature)
    relative_pdf = lambda g: 4 / math.sqrt(math.pi) * a**1.5 * g*g * math.exp(-a*g*g)
    result = quad(lambda g: float(kernel.rate(g)) * relative_pdf(g),
                  0, 5000., epsabs=1e-28)[0]
    mean_relative_speed = 4 * math.sqrt(BOLTZMANN_J_PER_K * temperature
                                      / (math.pi * sp.mass_kg))
    assert result == pytest.approx(math.pi * sp.diameter_m**2
                                   * mean_relative_speed, rel=1e-10)
    # Multiplying [m^3/s] by density [1/m^3] gives a collision frequency.
    assert 1e20 * result == pytest.approx(30814.3, rel=2e-4)


def test_pair_energy_momentum_and_known_right_angle_scatter():
    vp, wp = elastic_pair([3., 0., 0.], [-1., 0., 0.], [0., 1., 0.])
    np.testing.assert_allclose(vp, [1., 2., 0.], atol=0)
    np.testing.assert_allclose(wp, [1., -2., 0.], atol=0)
    rng = np.random.default_rng(713)
    v = rng.normal(size=(60, 3)) * 500
    w = rng.normal(size=(60, 3)) * 500
    d = rng.normal(size=(60, 3))
    d /= np.linalg.norm(d, axis=1)[:, None]
    vp, wp = elastic_pair(v, w, d)
    np.testing.assert_allclose(vp + wp, v + w, atol=4e-13)
    np.testing.assert_allclose(np.sum(vp*vp + wp*wp, axis=1),
                               np.sum(v*v + w*w, axis=1), rtol=8e-16)


def test_exact_zero_and_zero_speed_conventions():
    np.testing.assert_array_equal(ZeroKernel().rate([0., 2., 1e6]), [0., 0., 0.])
    assert VHSKernel().rate(0.) == 0
    assert VHSKernel(VhsSpecies(omega=1.)).rate(0.) == 0
    with pytest.raises(ValueError, match="unit norm"):
        elastic_pair([1, 0, 0], [-1, 0, 0], [0, 2, 0])


@pytest.mark.parametrize("law", ["VSS", "LarsenBorgnakkeVariableHardSphere", "Unknown"])
def test_unsupported_collision_laws_are_explicit(law):
    with pytest.raises(UnsupportedCollisionLaw):
        VHSKernel(law=law)


@pytest.mark.parametrize("g", [-1., float("nan"), float("inf")])
def test_invalid_relative_speeds(g):
    with pytest.raises(ValueError, match="relative speed"):
        VHSKernel().rate(g)
