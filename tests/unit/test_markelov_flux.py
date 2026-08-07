"""Inflow flux verification: mass, momentum and energy.

The half-range moments are checked three ways -- closed form, deterministic
quadrature, and Monte-Carlo sampling -- and the drift-only form is checked to be
*different*, since that is the whole reason the closed forms are needed.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from plumetools.markelov1999 import flux
from plumetools.markelov1999.constants import BOLTZMANN_J_PER_K, molecular_mass_kg

M_N2 = molecular_mass_kg(28.0134)
T = 300.0


@pytest.fixture
def one_face():
    """A single 1 m^2 face with its inward normal along +x and U along +x."""
    return dict(
        number_density=np.array([1.0e19]),
        U=np.array([[789.4849, 0.0, 0.0]]),
        T_K=np.array([T]),
        areas=np.array([1.0]),
        inward_normals=np.array([[1.0, 0.0, 0.0]]),
        mass_kg=M_N2,
        internal_dof=2,
    )


# --------------------------------------------------------------------------- #
# the flux coefficient itself
# --------------------------------------------------------------------------- #

def test_number_flux_coefficient_at_zero_speed_ratio():
    """At s = 0 the half-range flux is n*c/(2*sqrt(pi)), the classical effusion
    result. Checked against the closed form written out here, not against the
    implementation."""
    assert float(flux.number_flux_coefficient(0.0)) == pytest.approx(
        1.0 / (2.0 * math.sqrt(math.pi)), rel=1e-15)


def test_number_flux_coefficient_tends_to_the_speed_ratio():
    """As s grows the thermal terms become negligible and the flux tends to the
    drift value n*U. That limit is exactly what makes rho*U*A look adequate."""
    for s in (10.0, 100.0, 1000.0):
        assert float(flux.number_flux_coefficient(s)) == pytest.approx(s, rel=2e-2 / s)


def test_number_flux_coefficient_against_direct_quadrature():
    """F0 = integral_0^inf u exp(-(u-s)^2)/sqrt(pi) du, integrated numerically."""
    from scipy import integrate
    for s in (-1.0, 0.0, 0.5, 1.87, 4.0):
        reference, _ = integrate.quad(
            lambda u: u * math.exp(-(u - s) ** 2) / math.sqrt(math.pi),
            0.0, 40.0 + abs(s), limit=200)
        assert float(flux.number_flux_coefficient(s)) == pytest.approx(
            reference, rel=1e-9)


def test_normal_momentum_coefficient_against_direct_quadrature():
    from scipy import integrate
    for s in (0.0, 1.0, 1.87, 3.0):
        reference, _ = integrate.quad(
            lambda u: u * u * math.exp(-(u - s) ** 2) / math.sqrt(math.pi),
            0.0, 40.0 + abs(s), limit=200)
        assert float(flux.normal_momentum_coefficient(s)) == pytest.approx(
            reference, rel=1e-9)


def test_most_probable_speed_matches_the_solver_definition():
    """DSMCCloud::maxwellianMostProbableSpeed is sqrt(2 k T / m)."""
    assert float(flux.most_probable_speed(T, M_N2)) == pytest.approx(
        math.sqrt(2.0 * BOLTZMANN_J_PER_K * T / M_N2), rel=1e-15)


# --------------------------------------------------------------------------- #
# closed form against independent quadrature
# --------------------------------------------------------------------------- #

def test_analytical_and_quadrature_agree_on_one_face(one_face):
    analytical = flux.analytical_face_fluxes(**one_face)
    quadrature = flux.quadrature_face_fluxes(**one_face)

    assert analytical.number_per_s == pytest.approx(quadrature.number_per_s, rel=1e-9)
    assert analytical.mass_kg_per_s == pytest.approx(quadrature.mass_kg_per_s, rel=1e-9)
    np.testing.assert_allclose(analytical.momentum_N, quadrature.momentum_N, rtol=1e-9)
    assert analytical.energy_W == pytest.approx(quadrature.energy_W, rel=1e-9)


@pytest.mark.parametrize("speed", [0.0, 100.0, 789.4849, 2000.0])
def test_the_two_routes_agree_across_speed_ratios(one_face, speed):
    """Including s = 0, where the drift flux vanishes but the true flux does not."""
    case = dict(one_face, U=np.array([[speed, 0.0, 0.0]]))
    analytical = flux.analytical_face_fluxes(**case)
    quadrature = flux.quadrature_face_fluxes(**case)
    assert analytical.number_per_s == pytest.approx(quadrature.number_per_s, rel=1e-9)
    assert analytical.energy_W == pytest.approx(quadrature.energy_W, rel=1e-9)


def test_the_two_routes_agree_with_a_tangential_velocity_component(one_face):
    """The closed form for the energy flux assumes a normal drift, so the
    tangential drift kinetic energy is added separately. The quadrature includes
    it natively, so a disagreement here would mean the addition is wrong."""
    case = dict(one_face, U=np.array([[600.0, 400.0, 0.0]]))
    analytical = flux.analytical_face_fluxes(**case)
    quadrature = flux.quadrature_face_fluxes(**case)
    assert analytical.energy_W == pytest.approx(quadrature.energy_W, rel=1e-9)
    np.testing.assert_allclose(analytical.momentum_N, quadrature.momentum_N, rtol=1e-9)


def test_tangential_momentum_rides_only_on_the_number_flux(one_face):
    """The tangential velocity distribution is an undisplaced Maxwellian, so its
    thermal contribution to tangential momentum is zero by symmetry."""
    case = dict(one_face, U=np.array([[600.0, 400.0, 0.0]]))
    totals = flux.analytical_face_fluxes(**case)
    expected_y = M_N2 * 400.0 * totals.number_per_s
    assert float(totals.momentum_N[1]) == pytest.approx(expected_y, rel=1e-12)


# --------------------------------------------------------------------------- #
# why rho*U*A is not enough
# --------------------------------------------------------------------------- #

def test_drift_only_flux_is_zero_at_zero_drift_but_the_true_flux_is_not(one_face):
    """The clearest statement of the point: at s = 0 a drift-only check reports
    no inflow at all, while a real surface effuses n*c/(2*sqrt(pi))."""
    case = dict(one_face, U=np.zeros((1, 3)))
    drift = flux.drift_only_number_flux(
        case["number_density"], case["U"], case["inward_normals"])
    assert float(drift[0]) == 0.0

    true = flux.analytical_face_fluxes(**case)
    assert true.number_per_s > 0.0


def test_drift_only_underestimates_at_a_realistic_speed_ratio(one_face):
    """At this plume's speed ratio the thermal correction is small but real. The
    test records the size rather than asserting it is negligible."""
    drift = float(flux.drift_only_number_flux(
        one_face["number_density"], one_face["U"],
        one_face["inward_normals"])[0])
    true = flux.analytical_face_fluxes(**one_face).number_per_s
    assert drift < true
    assert 0.99 < drift / true < 1.0


def test_energy_flux_is_positive_even_with_no_drift(one_face):
    case = dict(one_face, U=np.zeros((1, 3)))
    assert flux.analytical_face_fluxes(**case).energy_W > 0.0


# --------------------------------------------------------------------------- #
# statistical check
# --------------------------------------------------------------------------- #

def test_sampled_number_flux_matches_the_closed_form():
    """Monte Carlo over the drifting Maxwellian, using no closed form."""
    n, u = 1.0e19, 789.4849
    sampled, stderr = flux.sampled_number_flux(
        number_density=n, U_normal=u, T_K=T, mass_kg=M_N2,
        n_samples=400_000, seed=7)

    c = flux.most_probable_speed(T, M_N2)
    expected = n * c * float(flux.number_flux_coefficient(u / c))
    assert abs(sampled - expected) < 5.0 * stderr


def test_sampled_number_flux_matches_at_zero_drift():
    n = 1.0e19
    sampled, stderr = flux.sampled_number_flux(
        number_density=n, U_normal=0.0, T_K=T, mass_kg=M_N2,
        n_samples=400_000, seed=11)
    c = flux.most_probable_speed(T, M_N2)
    expected = n * c / (2.0 * math.sqrt(math.pi))
    assert abs(sampled - expected) < 5.0 * stderr


def test_the_solvers_own_sampler_produces_the_predicted_mean_speed():
    """Bird eqn 12.5 acceptance-rejection, as plumeFieldInflow implements it.

    The accepted normal velocities follow the flux-weighted distribution, whose
    mean is F1/F0. This is what the two deterministic routes cannot check about
    each other: that the velocities the solver actually injects carry the
    momentum the closed form predicts.
    """
    c = flux.most_probable_speed(T, M_N2)
    for u_normal in (0.0, 400.0, 789.4849):
        s = u_normal / c
        mean, stderr = flux.sampled_mean_normal_speed(
            U_normal=u_normal, T_K=T, mass_kg=M_N2, n_samples=200_000, seed=3)
        expected = (float(flux.normal_momentum_coefficient(s))
                    / float(flux.number_flux_coefficient(s)))
        assert abs(mean - expected) < 5.0 * stderr, f"at s = {s:.3f}"


# --------------------------------------------------------------------------- #
# orientation
# --------------------------------------------------------------------------- #

def test_inward_normals_are_oriented_outward_from_the_source():
    """OpenFOAM winds a boundary face so its normal leaves the FLUID. The fluid is
    outside the inflow cavity, so the stored normal points at the origin and the
    into-the-domain direction is radially outward."""
    centroids = np.array([[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [-0.1, 0.0, 0.0]])
    stored = -centroids / np.linalg.norm(centroids, axis=1)[:, None]

    inward = flux.inward_normals_from_patch(stored, centroids)
    radial = centroids / np.linalg.norm(centroids, axis=1)[:, None]
    np.testing.assert_allclose(inward, radial, atol=1e-12)


def test_orientation_is_fixed_whichever_way_the_normals_were_stored():
    centroids = np.array([[0.1, 0.0, 0.0], [0.0, 0.1, 0.0]])
    radial = centroids / np.linalg.norm(centroids, axis=1)[:, None]
    for stored in (radial, -radial):
        np.testing.assert_allclose(
            flux.inward_normals_from_patch(stored, centroids), radial, atol=1e-12)


def test_a_reversed_normal_would_change_the_flux(one_face):
    """Stated so the sign convention is not taken on trust: the half-range
    integral is not antisymmetric, so getting it backwards gives a small positive
    flux rather than a negative one -- which would look plausible."""
    forward = flux.analytical_face_fluxes(**one_face)
    backward = flux.analytical_face_fluxes(
        **dict(one_face, inward_normals=-one_face["inward_normals"]))
    assert backward.number_per_s > 0.0
    assert backward.number_per_s < 0.01 * forward.number_per_s


# --------------------------------------------------------------------------- #
# integration over a patch
# --------------------------------------------------------------------------- #

def test_flux_is_additive_over_faces(one_face):
    """Two half-area faces must give the same total as one whole-area face."""
    single = flux.analytical_face_fluxes(**one_face)
    split = flux.analytical_face_fluxes(
        number_density=np.repeat(one_face["number_density"], 2),
        U=np.repeat(one_face["U"], 2, axis=0),
        T_K=np.repeat(one_face["T_K"], 2),
        areas=np.array([0.5, 0.5]),
        inward_normals=np.repeat(one_face["inward_normals"], 2, axis=0),
        mass_kg=M_N2, internal_dof=2)
    assert split.number_per_s == pytest.approx(single.number_per_s, rel=1e-14)
    assert split.energy_W == pytest.approx(single.energy_W, rel=1e-14)


def test_mass_flux_is_the_number_flux_times_the_molecular_mass(one_face):
    totals = flux.analytical_face_fluxes(**one_face)
    assert totals.mass_kg_per_s == pytest.approx(
        totals.number_per_s * M_N2, rel=1e-14)


def test_internal_degrees_of_freedom_add_energy_but_not_mass(one_face):
    monatomic = flux.analytical_face_fluxes(**dict(one_face, internal_dof=0))
    diatomic = flux.analytical_face_fluxes(**dict(one_face, internal_dof=2))

    assert diatomic.mass_kg_per_s == pytest.approx(monatomic.mass_kg_per_s, rel=1e-14)
    assert diatomic.energy_W > monatomic.energy_W

    # The difference is exactly (zeta/2) k T per injected molecule.
    extra = diatomic.energy_W - monatomic.energy_W
    assert extra == pytest.approx(
        BOLTZMANN_J_PER_K * T * diatomic.number_per_s, rel=1e-12)


def test_totals_serialise_for_the_case_summary(one_face):
    import json
    json.dumps(flux.analytical_face_fluxes(**one_face).as_dict())
