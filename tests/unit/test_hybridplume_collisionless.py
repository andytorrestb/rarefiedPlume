"""hybridplume: the independent collisionless reference.

It is checked against :mod:`plumetools.cai2012.analytical` -- a different
formulation (closed-form speed integrals, Gauss-Legendre azimuth) -- so the
agreement is evidence, not a self-test. The tolerance, 1e-9, is set by the
coarser of the two quadratures at the default orders; the measured agreement
is ~1e-13.
"""

from __future__ import annotations

import numpy as np
import pytest

from plumetools.cai2012 import analytical
from plumetools.hybridplume import collisionless
from plumetools.hybridplume.physics import BIRD_ARGON, SourceState


@pytest.fixture(scope="module")
def source():
    return SourceState(BIRD_ARGON, knudsen=0.1)


POINTS = np.array([[0.05, 0.0, 0.0], [0.2, 0.0, 0.0], [0.4, 0.05, 0.02],
                   [1.0, 0.3, -0.1], [2.0, 0.0, 0.5]])


def test_agrees_with_cai2012_closed_form(source):
    ref = collisionless.reference(POINTS, source)
    cf = analytical.moments(POINTS[:, 0], POINTS[:, 1], POINTS[:, 2], source.radius_m,
                            source.speed_ratio)
    assert np.allclose(ref["n_over_ne"], cf["n_over_n0"], rtol=1e-9)
    assert np.allclose(ref["U_over_ae"], cf["U_sqrt_beta0"], rtol=1e-9, atol=1e-12)
    assert np.allclose(ref["T_over_Te"], cf["T_over_T0"], rtol=1e-9)


def test_quadrature_converges(source):
    ref = collisionless.converged(POINTS, source)
    assert np.all(ref["quadrature_rel_change"] < 1e-9)


def test_near_the_exit_plane_the_error_estimate_says_so(source):
    """At x = 1e-5 m the kernel x/l^3 is 1e-5 m wide and the fixed quadrature
    cannot resolve it. The reference must REPORT that, not return a plausible
    number silently: the doubled-order change is large there and tiny at x = D/2."""
    near = collisionless.converged(np.array([[1e-5, 0.0, 0.0]]), source)
    far = collisionless.converged(np.array([[0.1, 0.0, 0.0]]), source)
    assert near["quadrature_rel_change"][0] > 1e-2
    assert far["quadrature_rel_change"][0] < 1e-10


def test_approaches_the_half_range_limit_on_the_axis(source):
    """n/ne -> (1 + erf S0)/2 as x -> 0: only forward-moving molecules exist."""
    from math import erf
    ref = collisionless.converged(np.array([[2e-3, 0.0, 0.0]]), source, order=(96, 192, 64))
    assert ref["quadrature_rel_change"][0] < 1e-6
    # at x = 0.01 D the closed form itself is 2.3e-4 below the limit
    assert ref["n_over_ne"][0] == pytest.approx(0.5 * (1 + erf(2.0)), rel=1e-3)


def test_rotational_symmetry_is_not_imposed_but_holds(source):
    r = 0.07
    pts = np.array([[0.3, r * np.cos(a), r * np.sin(a)] for a in np.radians([0, 23, 77, 141])])
    ref = collisionless.reference(pts, source)
    assert np.ptp(ref["n"]) / ref["n"].mean() < 1e-10
    assert np.ptp(ref["T"]) / ref["T"].mean() < 1e-10


def test_trace_and_directional_temperatures(source):
    ref = collisionless.reference(POINTS, source)
    assert np.allclose(ref["T_directional"].mean(axis=1), ref["T"], rtol=1e-12)
    # the far plume is anisotropic: axial spread differs from transverse
    far = ref["T_directional"][-1]
    assert far[0] / far[1] > 1.5


def test_refuses_the_exit_plane_and_the_lip(source):
    with pytest.raises(ValueError):
        collisionless.reference(np.array([[0.0, 0.0, 0.0]]), source)
    with pytest.raises(ValueError):
        collisionless.reference(np.array([[1e-4, 0.1, 0.0]]), source, lip_exclusion=1e-3)
