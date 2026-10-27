"""Independent physical limits and moment identities for reduced collisions.

These tests deliberately use known distributions, the existing Cai evaluator,
and conservation identities rather than reproducing the event implementation.
"""

import math

import numpy as np
import pytest

from plumetools.cai2012 import analytical
from plumetools.cai2012.constants import BOLTZMANN_J_PER_K as KB
from plumetools.cai2012.gas import ARGON
from plumetools.virtualcollision import DiskSource, NumericalOptions, VirtualCollisionModel
from plumetools.virtualcollision.kernel import VHSKernel, ZeroKernel
from plumetools.virtualcollision.moments import aggregate, derive
from plumetools.virtualcollision.model import ray_exit


def source(density=1.0e17, speed_ratio=2.0):
    return DiskSource(radius_m=0.1, number_density_m3=density,
                      temperature_K=300.0, speed_ratio=speed_ratio, species=ARGON)


def small_options(**changes):
    values = dict(samples=32, hazard_steps=8, partner_radial=4,
                  partner_azimuthal=8, partner_speed=12, seed=124,
                  domain_length_D=3.0, domain_radius_D=3.0, batches=4,
                  chunk_size=32)
    values.update(changes)
    return NumericalOptions(**values)


def test_disabled_mode_is_existing_cai_exactly_not_a_sample_limit():
    src = source()
    opts = small_options()
    points = np.array([[0.2, 0., 0.], [0.4, 0.08, 0.], [0.8, 0.12, -0.16]])
    result = VirtualCollisionModel(src, VHSKernel(ARGON), opts).evaluate(
        points, collisions=False)
    expected = analytical.moments(points[:, 0], points[:, 1], points[:, 2],
                                  src.radius_m, src.speed_ratio)
    speed = ARGON.most_probable_speed(src.temperature_K)
    np.testing.assert_array_equal(result["n"], expected["n_over_n0"] * src.number_density_m3)
    np.testing.assert_array_equal(result["U"], expected["U_sqrt_beta0"] * speed)
    np.testing.assert_array_equal(result["T_trans"], expected["T_over_T0"] * src.temperature_K)


def test_zero_kernel_uses_exact_cai_even_when_correction_enabled():
    src = source()
    model = VirtualCollisionModel(src, ZeroKernel(), small_options())
    points = np.array([[0.2, 0., 0.], [0.5, 0.1, 0.]])
    off = model.evaluate(points, collisions=False)
    zero = model.evaluate(points, half_width_m=(0.02, 0.02), collisions=True)
    for key in ("n", "U", "T_trans"):
        np.testing.assert_array_equal(zero[key], off[key])


def test_moment_tensor_recovers_independent_anisotropic_population():
    # Eight equiprobable corners have exactly diagonal covariance, with a drift.
    drift = np.array([2., -3., 5.])
    sigma = np.array([1., 2., 4.])
    signs = np.array([[i, j, k] for i in (-1., 1.) for j in (-1., 1.)
                      for k in (-1., 1.)])
    velocities = drift + signs * sigma
    n = 7.0
    result = derive(np.array([n]), n * velocities.mean(axis=0)[None, :],
                    n * np.einsum("ni,nj->ij", velocities, velocities)[None, :, :] / 8,
                    ARGON.mass_kg)
    np.testing.assert_allclose(result["U"][0], drift, rtol=0, atol=1e-14)
    expected = ARGON.mass_kg / KB * np.diag(sigma ** 2)
    np.testing.assert_allclose(result["T_tensor"][0], expected, rtol=1e-13, atol=1e-14)
    assert result["T_trans"][0] == pytest.approx(np.trace(expected) / 3)


def test_raw_aggregation_keeps_between_cell_velocity_dispersion():
    # Both cells are cold. The combined population is not cold: Ux=+-c.
    c = 17.0
    density = np.ones(2)
    velocity = np.array([[c, 0., 0.], [-c, 0., 0.]])
    second = np.einsum("ni,nj->nij", velocity, velocity)
    result = aggregate(density, velocity, second, np.ones(2), ARGON.mass_kg)
    np.testing.assert_allclose(result["U"], np.zeros(3), atol=1e-15)
    assert result["T_tensor"][0, 0] == pytest.approx(ARGON.mass_kg * c*c / KB)
    assert result["T_trans"] == pytest.approx(ARGON.mass_kg * c*c / (3 * KB))


def test_invalid_raw_density_and_tensor_are_preserved_and_flagged():
    density = np.array([-1., 1.])
    momentum = np.zeros((2, 3))
    second = np.array([np.eye(3), np.diag([-1., 1., 1.])])
    result = derive(density, momentum, second, ARGON.mass_kg)
    np.testing.assert_array_equal(result["M0"], density)
    np.testing.assert_array_equal(result["M2"], second)
    assert not np.any(result["valid"])
    assert result["T_tensor"][1, 0, 0] < 0


def test_flux_source_is_effusive_not_a_truncated_density_maxwellian():
    src = source(speed_ratio=0.)
    positions, velocities = src.sample_flux(16384, seed=781)
    a = ARGON.most_probable_speed(src.temperature_K)
    assert np.all(velocities[:, 0] > 0.)
    assert np.all(positions[:, 0] == 0.)
    assert np.max(np.linalg.norm(positions[:, 1:], axis=1)) <= src.radius_m
    assert np.mean(positions[:, 1:] ** 2) == pytest.approx(src.radius_m**2 / 4., rel=.01)
    assert np.mean(velocities[:, 0]) == pytest.approx(a * math.sqrt(math.pi) / 2., rel=.01)
    assert np.mean(velocities[:, 1:]) == pytest.approx(0., abs=.005*a)
    assert np.mean(velocities[:, 1:] ** 2) == pytest.approx(a*a / 2., rel=.01)


def test_phase_space_quadrature_converges_to_ballistic_density_and_moments():
    src = source(speed_ratio=0.)
    point = np.array([0.3, 0.08, 0.])
    exact = analytical.moments(*point, src.radius_m, src.speed_ratio)
    target = np.array([float(exact["n_over_n0"]),
                       float(exact["U_sqrt_beta0"][0]), float(exact["T_over_T0"])])
    errors = []
    for radial, azimuthal, speed in ((3, 6, 10), (20, 64, 32)):
        opts = small_options(partner_radial=radial, partner_azimuthal=azimuthal,
                             partner_speed=speed)
        velocities, weights = src.local_quadrature(point, opts)
        result = derive(np.array([weights.sum()]),
                        np.einsum("n,ni->i", weights, velocities)[None, :],
                        np.einsum("n,ni,nj->ij", weights, velocities, velocities)[None, :, :],
                        ARGON.mass_kg)
        actual = np.array([result["n"][0] / src.number_density_m3,
                           result["U"][0, 0] / ARGON.most_probable_speed(src.temperature_K),
                           result["T_trans"][0] / src.temperature_K])
        errors.append(np.linalg.norm(actual - target))
    assert errors[1] < errors[0]
    assert errors[1] < 1e-5


def test_representative_events_have_probability_and_two_particle_bookkeeping():
    model = VirtualCollisionModel(source(density=6.62e18), VHSKernel(ARGON),
        small_options(samples=16, hazard_steps=4, partner_radial=3,
                      partner_azimuthal=6, partner_speed=10))
    diagnostics = model.preprocess()
    exposure = np.array(diagnostics["exposure"])
    probabilities = np.array(diagnostics["event_probability"])
    uncollided = np.array(diagnostics["uncollided_weight"])
    assert np.all(exposure >= 0.)
    assert np.all((probabilities >= 0.) & (probabilities <= 1.))
    np.testing.assert_allclose(probabilities + uncollided, 1., atol=2e-16)
    np.testing.assert_allclose(uncollided, np.exp(-exposure), atol=0., rtol=2e-15)
    weights = np.array(diagnostics["symmetric_pair_weights"])
    assert np.all(weights > 0.)
    assert np.all(weights <= .5 * probabilities[probabilities > 0.])
    for batch in model.events:
        assert len(batch) % 4 == 0
        for offset in range(0, len(batch), 4):
            event = batch[offset:offset+4]
            positions = np.stack([row[0] for row in event])
            velocities = np.stack([row[1] for row in event])
            rates = np.array([row[2] for row in event])
            np.testing.assert_array_equal(positions, np.tile(positions[0], (4, 1)))
            np.testing.assert_allclose(rates, rates[0] * np.array([1., 1., -1., -1.]))
            np.testing.assert_allclose(rates @ velocities, 0., atol=2e-13 * abs(rates[0]) * 1000.)
            assert rates @ np.sum(velocities**2, axis=1) == pytest.approx(
                0., abs=2e-12 * abs(rates[0]) * 1e6)
    balance = diagnostics["closed_control_volume_flux"]
    assert abs(balance["number_relative_residual"]) < 1e-14
    assert np.linalg.norm(balance["momentum_residual_m_s"]) < 1e-10
    assert abs(balance["twice_energy_residual_m2_s2"]) < 1e-7
    assert sum(balance["outward_number_by_boundary_over_emission"].values()) == pytest.approx(1.)


def test_closed_control_volume_counts_return_side_and_downstream_exits():
    # Independent rays with analytically known first crossings of L=2, B=1.
    point = np.array([1., 0., 0.])
    assert ray_exit(point, np.array([2., 0., 0.]), 2., 1.) == (.5, "downstream")
    assert ray_exit(point, np.array([-2., 0., 0.]), 2., 1.) == (.5, "return")
    assert ray_exit(point, np.array([0., 2., 0.]), 2., 1.) == (.5, "side")
    # A single event redirects two incoming particles. All four signed rays
    # leave, even though the forward-only number correction is nonzero.
    velocities = np.array([[0., 2., 0.], [0., -2., 0.],
                           [2., 0., 0.], [-2., 0., 0.]])
    signs = np.array([1., 1., -1., -1.])
    crossings = [ray_exit(point, velocity, 2., 1.)[1] for velocity in velocities]
    assert np.sum(signs) == 0.
    np.testing.assert_array_equal(signs @ velocities, np.zeros(3))
    assert signs @ np.sum(velocities**2, axis=1) == 0.
    assert np.sum(signs[np.array(crossings) == "downstream"]) == -1.


def test_high_kn_limit_reduces_event_weight_without_fitted_coefficients():
    options = small_options(samples=16, hazard_steps=4, partner_radial=2,
                            partner_azimuthal=4, partner_speed=10)
    dilute = VirtualCollisionModel(source(density=1e10), VHSKernel(ARGON), options).preprocess()
    denser = VirtualCollisionModel(source(density=1e12), VHSKernel(ARGON), options).preprocess()
    dilute_probability = np.array(dilute["event_probability"])
    denser_probability = np.array(denser["event_probability"])
    assert np.max(dilute_probability) < 1e-7
    np.testing.assert_allclose(denser_probability / dilute_probability, 100., rtol=1e-5)
    np.testing.assert_allclose(np.array(denser["symmetric_pair_weights"])
                               / np.array(dilute["symmetric_pair_weights"]), 100., rtol=1e-4)


def test_ring_residence_moments_have_inverse_speed_not_crossing_weight():
    # One axial trajectory crosses an independently known dx=.4-m segment.
    # Its flux is its rate; its steady inventory is rate * dx/vx.
    model = VirtualCollisionModel(source(), ZeroKernel(), small_options())
    position = np.array([.1, .05, 0.])
    rate = 3.0
    bounds_x, bounds_r = np.array([[.2, .6]]), np.array([[.04, .06]])
    for speed in (2., 4.):
        velocity = np.array([speed, 0., 0.])
        zeroth, first, second = model._ray_moments(position, velocity, rate,
                                                  bounds_x, bounds_r)
        expected_inventory = rate * .4 / speed
        assert zeroth[0] == pytest.approx(expected_inventory, rel=2e-15)
        np.testing.assert_allclose(first[0], expected_inventory * velocity, rtol=2e-15)
        np.testing.assert_allclose(second[0], expected_inventory * np.outer(velocity, velocity), rtol=2e-15)
        assert first[0, 0] / .4 == pytest.approx(rate)
