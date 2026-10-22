"""The extension's off-axis comparisons: projections, limits, estimators, errors.

Every analytical value here is :func:`plumetools.cai2012.analytical.moments`,
already tested against its closed-form centreline in
``test_cai2012_analytical.py``. What is tested here is what the extension does
with it: the 3-D-to-meridional projections and their axis limits, the mass-flux
normalisation and its axisymmetric area factor (the collisionless flux through
the hemisphere must integrate to exactly one), pressure from the LOCAL
temperature, and that a field identical to the analytical one scores zero error
through every estimator.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from plumetools.cai2012 import analytical
from plumetools.cai2012.post import SampledField
from plumetools.cai2012ext import comparison

R0 = 0.1
D = 0.2
S0 = 2.0


class _Geom:
    diameter_m = D
    radius_m = R0


class _Exit:
    number_density_per_m3 = 6.62e18
    T0_K = 300.0
    speed_ratio = S0
    beta0 = 1.0 / (2.0 * (1.380649e-23 / 6.63e-26) * 300.0)


def synthetic_field(centres, *, volumes=None, scale_n=1.0, scale_t=1.0):
    """A 'sampled' field that IS the collisionless solution at its centres."""
    exact = comparison.analytical_fields(centres, R0, S0)
    velocity = np.column_stack([exact["ux"], np.zeros_like(exact["ux"]),
                                np.zeros_like(exact["ux"])])
    # Rebuild the Cartesian vector from the projections.
    r = exact["r"]
    safe = np.where(r > 0, r, 1.0)
    velocity[:, 1] = np.where(r > 0, exact["ur"] * centres[:, 1] / safe, 0.0)
    velocity[:, 2] = np.where(r > 0, exact["ur"] * centres[:, 2] / safe, 0.0)
    return SampledField(
        time="0.06",
        centres=np.asarray(centres, dtype=float),
        volumes=np.ones(len(centres)) if volumes is None else volumes,
        number_density=scale_n * exact["n_over_n0"] * _Exit.number_density_per_m3,
        velocity=velocity / math.sqrt(_Exit.beta0),
        temperature=scale_t * exact["T_over_T0"] * _Exit.T0_K,
    )


# --------------------------------------------------------------------------- #
# projections
# --------------------------------------------------------------------------- #

def test_cylindrical_radial_points_away_from_the_axis_on_every_side():
    centres = np.array([[1.0, 0.0, 0.5], [1.0, 0.0, -0.5], [1.0, 0.5, 0.0],
                        [1.0, -0.5, 0.0]])
    outward = centres.copy()
    outward[:, 0] = 0.0
    p = comparison.projections(centres, outward)
    assert np.allclose(p["ur"], 0.5)
    assert np.allclose(p["uphi"], 0.0)
    assert np.allclose(p["ux"], 0.0)


def test_swirl_is_the_azimuthal_component():
    centres = np.array([[1.0, 0.0, 1.0]])      # on +z: azimuthal is -y
    p = comparison.projections(centres, np.array([[0.0, -2.0, 0.0]]))
    assert p["uphi"][0] == pytest.approx(2.0)
    assert p["ur"][0] == pytest.approx(0.0)


def test_spherical_and_cylindrical_radial_are_different_quantities():
    # At 45 degrees a purely axial velocity has u_r = 0 but u_R = cos 45.
    centres = np.array([[1.0, 0.0, 1.0]])
    p = comparison.projections(centres, np.array([[1.0, 0.0, 0.0]]))
    assert p["ur"][0] == pytest.approx(0.0)
    assert p["uR"][0] == pytest.approx(math.cos(math.pi / 4))
    assert p["utheta"][0] == pytest.approx(-math.sin(math.pi / 4))
    assert math.degrees(p["theta"][0]) == pytest.approx(45.0)


def test_on_the_axis_the_transverse_components_take_their_limit_of_zero():
    p = comparison.projections(np.array([[2.0, 0.0, 0.0]]),
                               np.array([[1.0, 0.3, 0.4]]))
    assert p["ur"][0] == 0.0 and p["uphi"][0] == 0.0
    assert p["uR"][0] == pytest.approx(1.0)     # theta = 0: u_R = u_x


def test_the_analytical_radial_velocity_is_even_and_vanishes_on_the_axis():
    points = np.array([[0.4, 0.0, 0.3], [0.4, 0.0, -0.3], [0.4, 0.3, 0.0],
                       [0.4, 0.0, 0.0]])
    f = comparison.analytical_fields(points, R0, S0)
    assert f["ur"][0] == pytest.approx(f["ur"][1], rel=1e-12)
    assert f["ur"][0] == pytest.approx(f["ur"][2], rel=1e-12)
    assert f["ur"][0] > 0.0
    assert f["ur"][3] == 0.0
    assert np.max(np.abs(f["uphi"])) < 1e-12


def test_banding_reproduces_a_direct_call_near_the_exit_plane():
    points = np.array([[0.002, 0.0, 0.05], [0.05, 0.0, 0.05], [1.0, 0.0, 0.3]])
    banded = comparison.analytical_fields(points, R0, S0)
    for i, point in enumerate(points):
        direct = analytical.moments(point[0], point[1], point[2], R0, S0)
        assert banded["n_over_n0"][i] == pytest.approx(
            float(direct["n_over_n0"]), rel=1e-10)
        assert banded["T_over_T0"][i] == pytest.approx(
            float(direct["T_over_T0"]), rel=1e-10)


def test_points_in_the_exit_plane_are_nan_rather_than_an_error():
    f = comparison.analytical_fields(np.array([[0.0, 0.0, 0.5]]), R0, S0)
    assert np.isnan(f["n_over_n0"][0])


# --------------------------------------------------------------------------- #
# pressure and mass flux
# --------------------------------------------------------------------------- #

def test_pressure_uses_the_local_temperature_not_the_exit_one():
    centres = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.1]])
    sampled = synthetic_field(centres, scale_t=0.5)
    f = comparison.dsmc_fields(sampled, _Exit)
    assert np.allclose(f["p_over_p0"], f["n_over_n0"] * f["T_over_T0"])
    assert not np.allclose(f["p_over_p0"], f["n_over_n0"])


def test_the_exit_flux_factor_has_its_effusion_and_beam_limits():
    assert comparison.exit_flux_factor(0.0) == pytest.approx(
        1.0 / (2.0 * math.sqrt(math.pi)))
    assert comparison.exit_flux_factor(8.0) == pytest.approx(8.0, rel=1e-6)


@pytest.mark.parametrize("radius", [2.0, 1.0])
def test_the_collisionless_flux_through_the_hemisphere_integrates_to_one(radius):
    """Every molecule leaves the disk forward and crosses the sphere once.

    So the solid-angle integral of F, with the axisymmetric factor
    2 pi sin(theta) dtheta, is the exit flow: exactly one. A wrong Gamma0, a
    missing R^2 or a flat area factor would each miss by far more than the
    quadrature error at this resolution.
    """
    edges = np.linspace(0.0, math.pi / 2, 361)
    mid = 0.5 * (edges[:-1] + edges[1:])
    points = np.column_stack([radius * np.cos(mid), np.zeros_like(mid),
                              radius * np.sin(mid)])
    f = comparison.analytical_fields(points, R0, S0)
    flux = comparison.mass_flux_per_steradian(f["n_over_n0"], f["uR"], f["R"],
                                              R0, S0)
    omega = 2.0 * math.pi * (np.cos(edges[:-1]) - np.cos(edges[1:]))
    assert np.sum(flux * omega) == pytest.approx(1.0, abs=5e-5)


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #

def test_relative_errors_skip_the_excluded_and_count_them():
    measured = np.array([1.1, 2.0, 5.0, np.nan])
    exact = np.array([1.0, 2.0, 1.0e-9, 1.0])
    valid = np.array([True, True, False, True])
    m = comparison.masked_metrics(measured, exact, name="q", valid=valid,
                                  scale=1.0)
    assert m["q_n_valid"] == 2
    assert m["q_n_excluded"] == 2
    assert m["q_max_rel_error"] == pytest.approx(0.1)
    assert m["q_mean_rel_error"] == pytest.approx(0.05)
    assert m["q_max_abs_error"] == pytest.approx(0.1)


def test_a_quantity_through_zero_gets_no_relative_error():
    m = comparison.masked_metrics(np.array([0.1, -0.1]), np.array([0.0, 0.0]),
                                  name="ur", valid=np.array([True, True]),
                                  scale=2.0, relative=False)
    assert "ur_max_rel_error" not in m
    assert m["ur_max_abs_error"] == pytest.approx(0.05)
    assert m["ur_n_valid"] == 2


def test_nothing_valid_is_nan_not_a_crash():
    m = comparison.masked_metrics(np.ones(3), np.ones(3), name="q",
                                  valid=np.zeros(3, dtype=bool), scale=1.0)
    assert math.isnan(m["q_max_rel_error"]) and m["q_n_valid"] == 0


def test_the_centreline_pressure_is_the_product_of_the_binned_moments():
    profile = {key: np.array([0.5, 0.1]) for key in (
        "n_over_n0", "n_over_n0_analytical", "n_over_n0_analytical_axis")}
    profile.update({key: np.array([0.8, 0.4]) for key in (
        "T_over_T0", "T_over_T0_analytical", "T_over_T0_analytical_axis")})
    profile.update({key: np.array([2.0, 2.2]) for key in (
        "U_sqrt_beta0", "U_sqrt_beta0_analytical")})
    out, metrics = comparison.centerline_extended(profile, S0)
    assert np.allclose(out["p_over_p0"], [0.4, 0.04])
    assert metrics["pressure_max_rel_error"] == pytest.approx(0.0)
    assert metrics["velocity_abs_error_scale"] == S0


# --------------------------------------------------------------------------- #
# the estimators, on a field that is the analytical solution
# --------------------------------------------------------------------------- #

def _shell_cells():
    radii = (1.85, 2.0, 2.1)
    thetas = np.radians(np.arange(0.5, 90.0, 1.0))
    phis = np.radians(np.arange(0.0, 360.0, 45.0))
    rr, tt, pp = np.meshgrid(radii, thetas, phis, indexing="ij")
    return np.column_stack([(rr * np.cos(tt)).ravel(),
                            (rr * np.sin(tt) * np.cos(pp)).ravel(),
                            (rr * np.sin(tt) * np.sin(pp)).ravel()])


def _shell_field(**kwargs):
    centres = _shell_cells()
    return synthetic_field(centres, volumes=np.full(len(centres), 0.03 ** 3),
                           **kwargs)


def test_an_exact_field_scores_zero_on_the_sphere():
    profile = comparison.angular_profile(_shell_field(), _Geom, _Exit)
    metrics = comparison.angular_metrics(profile, S0)
    for name in ("angular_density", "angular_mass_flux",
                 "angular_radial_velocity", "angular_temperature"):
        assert metrics[f"{name}_max_abs_error"] < 1e-12, name
    # The floors are relative to the sphere's own peak (1e-3), so the
    # comparison reaches ~75 deg -- not the ~40 deg an n/n0 >= 1e-3 floor
    # would allow at r/D = 10.
    assert metrics["angular_density_n_valid"] >= 28
    assert metrics["angular_mass_flux_n_valid"] >= 28
    assert profile["integrals"]["F"] == pytest.approx(
        profile["integrals"]["F_analytical"], rel=1e-12)
    assert profile["integrals"]["F_exact"] == pytest.approx(1.0, abs=2e-3)
    # Every bin is covered, so the same-cells flux closes too.
    assert profile["integrals"]["bins_without_cells"] == 0
    assert profile["integrals"]["F_analytical"] == pytest.approx(1.0, abs=0.02)


def test_cells_wider_than_a_bin_are_shared_and_leave_no_bin_empty():
    """At r/D = 10 the graded mesh's cells are as wide as a 2.5 deg bin."""
    thetas = np.radians(np.arange(2.0, 90.0, 6.0))      # one cell per 6 deg
    centres = np.column_stack([2.0 * np.cos(thetas), np.zeros_like(thetas),
                               2.0 * np.sin(thetas)])
    edges = np.radians(np.arange(0.0, 92.5, 2.5))
    volumes = np.full(len(centres), (2.0 * math.radians(6.0)) ** 3)
    fractions = comparison.footprint_fractions(
        np.arctan2(np.hypot(centres[:, 1], centres[:, 2]), centres[:, 0]),
        np.linalg.norm(centres, axis=1), volumes, edges)
    assert np.allclose(fractions.sum(axis=1), 1.0)
    assert np.all(fractions.sum(axis=0) > 0.0)
    sampled = synthetic_field(centres, volumes=volumes)
    profile = comparison.angular_profile(sampled, _Geom, _Exit)
    assert profile["integrals"]["bins_without_cells"] == 0


def test_the_shell_estimator_is_compared_with_itself_not_with_the_sphere():
    """A thick shell sits off the sphere; the matched value absorbs that."""
    sampled = _shell_field()
    profile = comparison.angular_profile(sampled, _Geom, _Exit)
    # Inside the shell the density differs from the exact sphere value ...
    assert not np.allclose(profile["n_over_n0"], profile["n_over_n0_exact"],
                           rtol=1e-3)
    # ... and matches the same-cells analytical value exactly.
    assert np.allclose(profile["n_over_n0"], profile["n_over_n0_analytical"],
                       rtol=1e-12)


def test_a_uniform_deficit_shows_as_that_relative_error():
    sampled = _shell_field(scale_n=0.9)
    metrics = comparison.angular_metrics(
        comparison.angular_profile(sampled, _Geom, _Exit), S0)
    assert metrics["angular_density_mean_rel_error"] == pytest.approx(0.1)
    assert metrics["angular_mass_flux_mean_rel_error"] == pytest.approx(0.1)


class _Post:
    plane = "xz"
    plane_half_over_D = 5.0
    centerline_x_over_D_max = 10.0


class _Cfg:
    post = _Post


def test_the_plane_scores_zero_for_an_exact_field_and_counts_the_lip():
    dx = 0.02
    xs = np.arange(dx / 2, 1.0, dx)
    zs = np.arange(-0.5 + dx / 2, 0.5, dx)
    xx, zz = np.meshgrid(xs, zs, indexing="ij")
    layer = np.column_stack([xx.ravel(), np.full(xx.size, dx / 2), zz.ravel()])
    centres = np.vstack([layer, layer * np.array([1.0, -1.0, 1.0])])
    sampled = synthetic_field(centres, volumes=np.full(len(centres), dx ** 3))
    plane = comparison.meridional_plane(sampled, _Cfg, _Geom, _Exit)
    metrics = comparison.plane_metrics(plane, S0)
    assert metrics["plane_cells"] == len(centres)
    assert metrics["plane_excluded_nozzle_lip"] > 0
    for name in ("plane_density", "plane_axial_velocity",
                 "plane_radial_velocity", "plane_temperature", "plane_pressure"):
        assert metrics[f"{name}_max_abs_error"] < 1e-12, name
        assert metrics[f"{name}_n_valid"] > 100, name
    assert "plane_radial_velocity_max_rel_error" not in metrics
