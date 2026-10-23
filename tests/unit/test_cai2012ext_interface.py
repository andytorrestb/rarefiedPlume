"""The interface physics: flux targets, the equivalent source, Bird, the sphere.

Independent checks wherever there is one: numerical quadrature for the one-way
moments, finite differences for the density gradient, the frozen 1999
implementation for the angular law and the terminal density, the scalar gas
module for the mean free path.
"""

from __future__ import annotations

import dataclasses
import math
from pathlib import Path

import numpy as np
import pytest
from scipy import integrate

from plumetools.cai2012 import gas
from plumetools.cai2012.config import load_case_config
from plumetools.cai2012ext import interface as itf
from plumetools.markelov1999 import flux as markelov_flux
from plumetools.markelov1999 import sourceflow as sf

REPO = Path(__file__).resolve().parents[2]
TEMPLATE = REPO / "cases" / "cai2012-kn-extension" / "baseCase"


def config_at(knudsen: float, **geometry):
    base = load_case_config(TEMPLATE)
    cfg = dataclasses.replace(base, exit=dataclasses.replace(base.exit, knudsen=knudsen))
    if geometry:
        cfg = dataclasses.replace(cfg, geometry=dataclasses.replace(cfg.geometry, **geometry))
    return cfg


@pytest.fixture(scope="module")
def source():
    return itf.equivalent_source(config_at(1.0e-3))


def settings(**values):
    defaults = dict(enabled=True, purpose="development-smoke", run_mode="exploratory",
                    bird_threshold=0.05, core_half_angle_deg=30.0,
                    source_validity_min_radius_over_D=1.0)
    defaults.update(values)
    return itf.InterfaceConfig(**defaults)


def select(source, cfg, **values):
    s = settings(**values)
    r_max, why = itf.domain_radius_limit(cfg, s)
    return itf.select_sphere(source, s, r_max_m=r_max, r_max_reason=why)


# --------------------------------------------------------------------------- #
# 1. angular law, normalisation, core fraction
# --------------------------------------------------------------------------- #

def test_monatomic_limiting_angle_is_ninety_degrees(source):
    assert source.gamma == pytest.approx(5.0 / 3.0)
    assert source.theta_L == pytest.approx(0.5 * math.pi, abs=1e-15)


def test_angular_law_is_the_frozen_1999_function(source):
    theta = np.linspace(0.0, 0.5 * math.pi, 37)
    np.testing.assert_array_equal(source.angular(theta),
                                  sf.angular(source.gamma, theta, 0.41))
    assert source.I == sf.normalization_integral(source.gamma, 0.41)


def test_mass_per_steradian_integrates_to_the_target(source):
    total, _ = integrate.quad(lambda t: float(source.j(t)) * 2.0 * math.pi * math.sin(t),
                              0.0, source.theta_L, epsrel=1e-12)
    assert total == pytest.approx(source.targets.mdot_kg_per_s, rel=1e-10)


def test_core_fraction_is_unrenormalised_and_bounded(source):
    k = source.k_theta
    previous = 0.0
    for degrees in (0.0, 5.0, 15.0, 30.0, 45.0, 60.0, 89.0, 90.0, 120.0):
        fraction = source.core_fraction(math.radians(degrees))
        assert 0.0 <= fraction <= 1.0
        assert fraction >= previous
        previous = fraction
        # gamma = 5/3: f = cos^k, so the fraction is 1 - cos^(k+1) exactly.
        # 1e-10: the frozen normalization_integral is quad at its default
        # tolerances, 1.2e-11 from the exact 1/(k+1).
        exact = 1.0 - math.cos(math.radians(min(degrees, 90.0))) ** (k + 1.0)
        assert fraction == pytest.approx(exact, abs=1e-10)
    assert source.core_fraction(source.theta_L) == pytest.approx(1.0, abs=1e-10)


# --------------------------------------------------------------------------- #
# 2. the anchor, and conservation along the streamtubes
# --------------------------------------------------------------------------- #

def test_anchor_reproduces_the_exit_state(source):
    st = source.state(np.array([source.r_a]), np.array([0.0]))
    assert st["U"][0] == pytest.approx(source.U_e, rel=1e-12)
    assert st["T"][0] == pytest.approx(source.T_a, rel=1e-12)
    assert st["rho"][0] == pytest.approx(source.rho_e, rel=1e-12)
    assert st["M"][0] == pytest.approx(source.M_a, rel=1e-12)


@pytest.mark.parametrize("r_over_ra", [1.0, 1.7, 4.0, 15.0, 80.0])
@pytest.mark.parametrize("theta_deg", [0.0, 20.0, 55.0, 85.0])
def test_mass_and_energy_are_conserved_along_every_ray(source, r_over_ra, theta_deg):
    r = np.array([r_over_ra * source.r_a])
    theta = np.array([math.radians(theta_deg)])
    st = source.state(r, theta)
    assert (st["rho"] * st["U"] * r * r)[0] == pytest.approx(float(source.j(theta)[0]),
                                                             rel=1e-11)
    h = source.cp * st["T"] + 0.5 * st["U"] ** 2
    assert h[0] == pytest.approx(source.targets.h_J_per_kg, rel=1e-12)


def test_gross_targets_hold_the_exit_flux_not_its_net_flow(source):
    t = source.targets
    # S0 = 2: the gross flux exceeds the net by the returning tail, 0.02 %.
    assert 1.0 < t.mdot_kg_per_s / t.mdot_net_kg_per_s < 1.001
    assert t.h_J_per_kg == pytest.approx(t.h_net_J_per_kg, rel=1e-3)
    assert source.T_t == pytest.approx(t.h_J_per_kg / source.cp)


# --------------------------------------------------------------------------- #
# 3. area-Mach, the supersonic branch, the terminal speed
# --------------------------------------------------------------------------- #

def test_area_mach_and_isentropic_density_agree(source):
    r = source.r_a * np.geomspace(1.0, 200.0, 25)
    st = source.state(r, np.zeros_like(r))
    np.testing.assert_allclose(itf.area_mach(st["M"], source.gamma),
                               itf.area_mach(source.M_a, source.gamma) * (r / source.r_a) ** 2,
                               rtol=1e-12)
    np.testing.assert_allclose(st["rho"] / source.rho_e,
                               (st["T"] / source.T_a) ** (1.0 / (source.gamma - 1.0)),
                               rtol=1e-10)


def test_the_solver_never_takes_the_subsonic_root(source):
    r = source.r_a * np.geomspace(1.0, 50.0, 40)
    mach = source.mach(r)
    assert np.all(mach >= source.M_a * (1 - 1e-12))
    assert np.all(np.diff(mach) > 0.0)
    # Below the sonic area ratio there is no supersonic root at all.
    assert np.isnan(itf.supersonic_mach(0.5, source.gamma))


def test_speed_approaches_the_markelov_limiting_speed(source):
    far = source.state(np.array([1.0e5 * source.r_a]), np.array([0.0]))
    assert source.V_inf == pytest.approx(math.sqrt(2.0 * source.targets.h_J_per_kg),
                                         rel=1e-12)
    assert source.V_inf == pytest.approx(
        sf.limiting_velocity(source.gamma, source.T_t, source.mass_kg), rel=1e-15)
    assert far["U"][0] == pytest.approx(source.V_inf, rel=1e-5)
    assert far["U"][0] < source.V_inf


def test_terminal_model_is_the_closure_at_infinite_mach(source):
    """The frozen 1999 density, through p_t r_o^2, is j / (V_inf r^2)."""
    r = source.r_a * np.array([1.0, 3.0, 30.0])
    theta = np.radians([0.0, 25.0, 60.0])
    terminal = itf.markelov_terminal(source, r, theta)
    np.testing.assert_allclose(terminal["rho"], source.j(theta) / (source.V_inf * r * r),
                               rtol=1e-12)
    rows = itf.compare_with_terminal_model(source, [2.0, 10.0])
    for row in rows:
        assert row["h_closure_over_target"] == pytest.approx(1.0, rel=1e-12)
        assert row["h_terminal_over_target"] == pytest.approx(2.0, rel=1e-3)
        assert row["rho_closure_over_terminal"] == pytest.approx(
            source.V_inf / row["U_closure"], rel=1e-12)


def test_stewart_spherical_regression():
    """Stewart & Lumpkin Eq. (2) at r* = 1 m, T* = 1000 K, gamma = 5/3."""
    gamma, r_star, t_star, d = 5.0 / 3.0, 1.0, 1000.0, 4.17e-10
    mach = np.array([1.0, 1.5, 3.0, 10.0])
    out = itf.stewart_spherical(mach, gamma)
    assert out["T_over_T_star"][0] == pytest.approx(1.0)
    assert out["r_over_r_star"][0] == pytest.approx(1.0)
    np.testing.assert_allclose(out["n_over_n_star"],
                               out["T_over_T_star"] ** (1.0 / (gamma - 1.0)))
    # The same relations as the closure's: F(M) is r^2 / r*^2.
    np.testing.assert_allclose(out["r_over_r_star"] ** 2, itf.area_mach(mach, gamma))
    temperature = t_star * out["T_over_T_star"]
    # Eq. (4) inverts the radial Bird form exactly, hard sphere only.
    radius = r_star * out["r_over_r_star"]
    target = 0.05
    n = itf.stewart_hard_sphere_density(mach, radius, target, gamma, d)
    species = gas.VhsSpecies(name="Ar", mass_kg=6.63e-26, diameter_m=d, omega=0.5,
                             t_ref_K=273.0)
    big_r = species.gas_constant_j_per_kg_k
    state = {"M": mach, "lambda": itf.mean_free_path(n, temperature, species, "hard_sphere"),
             "U": mach * np.sqrt(gamma * big_r * temperature),
             "c_bar": np.sqrt(8.0 * big_r * temperature / math.pi)}
    np.testing.assert_allclose(itf.bird_parameter(state, radius, "radial_approx"), target,
                               rtol=1e-12)


# --------------------------------------------------------------------------- #
# 4. Bird's parameter
# --------------------------------------------------------------------------- #

def test_vectorised_mean_free_path_matches_the_gas_module(source):
    for n, t in ((1e20, 300.0), (3e18, 40.0), (7e22, 900.0)):
        for convention in ("vhs", "hard_sphere"):
            assert float(itf.mean_free_path(n, t, source.species, convention)) == \
                pytest.approx(gas.mean_free_path(n, source.species, t,
                                                 convention=convention), rel=1e-14)


@pytest.mark.parametrize("r_over_ra", [1.05, 2.0, 7.0, 30.0])
@pytest.mark.parametrize("theta_deg", [0.0, 40.0])
def test_bird_general_is_the_finite_difference_density_gradient(source, r_over_ra,
                                                                 theta_deg):
    r = r_over_ra * source.r_a
    theta = math.radians(theta_deg)
    h = 1.0e-6 * r
    rho = lambda x: source.state(np.array([x]), np.array([theta]))["rho"][0]
    dlnrho = (math.log(rho(r + h)) - math.log(rho(r - h))) / (2.0 * h)
    st = source.state(np.array([r]), np.array([theta]))
    # B = -(lambda/(rho c_bar)) U . grad(rho), U radial: -(lambda U/c_bar) dln(rho)/dr.
    finite = -(st["lambda"][0] * st["U"][0] / st["c_bar"][0]) * dlnrho
    assert finite > 0.0  # signed: an expansion
    assert float(source.bird(np.array([r]), np.array([theta]))[0]) == pytest.approx(
        finite, rel=1e-6)


def test_general_over_radial_is_the_acceleration_factor(source):
    r = source.r_a * np.geomspace(1.0, 40.0, 12)
    theta = np.zeros_like(r)
    general = source.bird(r, theta, "general")
    radial = source.bird(r, theta, "radial_approx")
    mach2 = source.mach(r) ** 2
    np.testing.assert_allclose(general / radial, mach2 / (mach2 - 1.0), rtol=1e-12)
    with pytest.raises(ValueError):
        source.bird(r, theta, "knudsen")


def test_bird_scales_with_knudsen_and_nothing_else():
    a = itf.equivalent_source(config_at(1.0e-3))
    b = itf.equivalent_source(config_at(1.0e-4))
    r = a.r_a * np.array([1.0, 3.0, 9.0])
    np.testing.assert_allclose(a.bird(r, 0 * r) / b.bird(r, 0 * r), 10.0, rtol=1e-12)


# --------------------------------------------------------------------------- #
# 5. one-way fluxes
# --------------------------------------------------------------------------- #

def _quadrature(u_n, a):
    def weight(v):
        return math.exp(-((v - u_n) / a) ** 2) / (math.sqrt(math.pi) * a)
    upper = max(u_n, 0.0) + 12.0 * a
    j1, _ = integrate.quad(lambda v: v * weight(v), 0.0, upper, epsabs=0.0, epsrel=1e-13,
                           limit=400)
    j3, _ = integrate.quad(lambda v: v ** 3 * weight(v), 0.0, upper, epsabs=0.0,
                           epsrel=1e-13, limit=400)
    return j1, j3


@pytest.mark.parametrize("s", [-8.0, -4.0, -1.0, 0.0, 0.3, 2.0, 5.0, 12.0])
def test_one_way_moments_match_independent_quadrature(s):
    big_r, temperature = 208.0, 300.0
    a = math.sqrt(2.0 * big_r * temperature)
    j1, j3 = itf.one_way_moments(s * a, temperature, big_r)
    q1, q3 = _quadrature(s * a, a)
    assert float(j1) == pytest.approx(q1, rel=1e-9)
    assert float(j3) == pytest.approx(q3, rel=1e-9)


def test_one_way_moments_have_the_drift_and_effusion_limits():
    big_r, temperature = 208.0, 300.0
    a = math.sqrt(2.0 * big_r * temperature)
    j1, j3 = itf.one_way_moments(0.0, temperature, big_r)
    assert float(j1) == pytest.approx(a / (2.0 * math.sqrt(math.pi)), rel=1e-15)
    assert float(j3) == pytest.approx(a ** 3 / (2.0 * math.sqrt(math.pi)), rel=1e-15)
    u = 40.0 * a
    j1, j3 = itf.one_way_moments(u, temperature, big_r)
    assert float(j1) == pytest.approx(u, rel=1e-15)
    assert float(j3) == pytest.approx(u ** 3 + 1.5 * u * a * a, rel=1e-15)


def test_one_way_moments_are_stable_far_upstream():
    """No cancellation garbage: positive, finite, decaying, exact zero at underflow."""
    big_r, temperature = 208.0, 300.0
    a = math.sqrt(2.0 * big_r * temperature)
    s = -np.array([5.0, 10.0, 20.0, 26.0, 40.0, 1e4])
    j1, j3 = itf.one_way_moments(s * a, temperature, big_r)
    assert np.all(np.isfinite(j1)) and np.all(np.isfinite(j3))
    assert np.all(j1 >= 0.0) and np.all(j3 >= 0.0)
    assert np.all(np.diff(j1) <= 0.0)
    assert j1[-1] == 0.0
    # The leading asymptotic term a exp(-s^2) / (4 sqrt(pi) s^2).
    expected = a * math.exp(-400.0) / (4.0 * math.sqrt(math.pi) * 400.0)
    assert float(j1[2]) == pytest.approx(expected, rel=5e-3)


def test_face_fluxes_agree_with_the_tested_1999_flux_module(source):
    centres, normals, areas = itf.ideal_sphere_faces(3.0 * source.r_a, 0.5 * math.pi,
                                                      n_theta=24, n_phi=48)
    faces = itf.face_inflow(source, settings(tails="inject_exploratory"),
                            centres=centres, outward_normals=normals, areas=areas)
    balance = itf.flux_balance(source, faces, 1.0)
    reference = markelov_flux.analytical_face_fluxes(
        number_density=faces.number_density, U=faces.velocity, T_K=faces.temperature,
        areas=areas, inward_normals=faces.inward_normals, mass_kg=source.mass_kg,
        internal_dof=0)
    assert balance["gross_injected_mass_kg_per_s"] == pytest.approx(
        reference.mass_kg_per_s, rel=1e-12)
    assert balance["gross_injected_energy_W"] == pytest.approx(reference.energy_W,
                                                               rel=1e-12)


# --------------------------------------------------------------------------- #
# 6. Kn scales the density and the flows, nothing else
# --------------------------------------------------------------------------- #

def test_knudsen_scales_density_and_flow_only():
    a = itf.equivalent_source(config_at(1.0e-3))
    b = itf.equivalent_source(config_at(1.0e-4))
    assert b.n_e / a.n_e == pytest.approx(10.0, rel=1e-12)
    assert b.targets.mdot_kg_per_s / a.targets.mdot_kg_per_s == pytest.approx(10.0,
                                                                              rel=1e-12)
    assert b.targets.edot_W / a.targets.edot_W == pytest.approx(10.0, rel=1e-12)
    for name in ("diameter_m", "T_e", "U_e", "T_t", "T_a", "M_a", "r_a", "r_star",
                 "V_inf", "I", "A_M", "theta_L"):
        assert getattr(a, name) == pytest.approx(getattr(b, name), rel=1e-12), name
    theta = np.radians([0.0, 30.0, 70.0])
    r = np.full(3, 4.0 * a.r_a)
    np.testing.assert_allclose(a.state(r, theta)["n"] * 10.0, b.state(r, theta)["n"],
                               rtol=1e-12)


# --------------------------------------------------------------------------- #
# 7. choosing the sphere
# --------------------------------------------------------------------------- #

def test_a_valid_case_finds_the_bird_root(source):
    cfg = config_at(1.0e-3)
    pick = select(source, cfg)
    assert pick.status == itf.STATUS_ADMISSIBLE and pick.limited_by == "bird"
    assert pick.max_core_bird_at_radius == pytest.approx(0.05, rel=1e-9)
    assert pick.radius_m > max(source.r_a, 1.0 * source.diameter_m)
    # The core EDGE decides, not the axis: the axis is still well below.
    assert pick.axis_bird_at_radius < 0.05
    assert not pick.non_monotonic


def test_angular_edge_breakdown_shrinks_the_sphere(source):
    cfg = config_at(1.0e-3)
    narrow = select(source, cfg, core_half_angle_deg=10.0)
    wide = select(source, cfg, core_half_angle_deg=45.0,
                  source_validity_min_radius_over_D=0.75)
    assert wide.radius_m < narrow.radius_m
    # Wide enough that the edge already breaks down at the smallest radius.
    edge = select(source, cfg, core_half_angle_deg=60.0)
    assert edge.status == itf.STATUS_NO_INTERFACE
    assert edge.axis_bird_at_r_min < 0.05 < edge.max_core_bird_at_r_min


def test_no_root_at_kn_001_for_any_core_angle():
    cfg = config_at(0.01)
    source = itf.equivalent_source(cfg)
    for threshold in (0.02, 0.03, 0.05):
        pick = select(source, cfg, bird_threshold=threshold, core_half_angle_deg=None,
                      source_validity_min_radius_over_D=None)
        assert pick.status == itf.STATUS_NO_INTERFACE
        assert pick.radius_m is None
        assert pick.axis_bird_at_r_min == pytest.approx(0.0624, abs=5e-4)


def test_insufficient_domain_is_reported(source):
    cfg = config_at(1.0e-3, x_max_over_D=1.2, y_half_over_D=1.2, z_half_over_D=1.2)
    pick = select(source, cfg)
    assert pick.status == itf.STATUS_INSUFFICIENT_DOMAIN


def test_domain_limit_is_reported_not_hidden():
    cfg = config_at(1.0e-4)
    source = itf.equivalent_source(cfg)
    pick = select(source, cfg)
    assert pick.status == itf.STATUS_ADMISSIBLE
    assert pick.limited_by in ("domain box", "diagnostic sphere")
    assert pick.max_core_bird_at_radius < 0.05


def test_invalid_and_incomplete_states(source):
    cfg = config_at(1.0e-3)
    s = dataclasses.replace(settings(), core_half_angle_deg=None)
    r_max, _ = itf.domain_radius_limit(cfg, s)
    assert itf.select_sphere(source, s, r_max_m=r_max).status == itf.STATUS_INCOMPLETE
    s = dataclasses.replace(settings(), core_half_angle_deg=95.0)
    assert itf.select_sphere(source, s, r_max_m=r_max).status == itf.STATUS_INVALID_STATE
    s = dataclasses.replace(settings(), purpose="production",
                            source_validity_min_radius_over_D=None)
    assert itf.select_sphere(source, s, r_max_m=r_max).status == itf.STATUS_INCOMPLETE


# --------------------------------------------------------------------------- #
# 8. truncation without renormalisation
# --------------------------------------------------------------------------- #

def test_core_truncation_injects_exactly_the_core_fraction(source):
    radius = 3.0 * source.r_a
    centres, normals, areas = itf.ideal_sphere_faces(radius, 0.5 * math.pi,
                                                      n_theta=720, n_phi=64)
    s = settings(core_half_angle_deg=30.0)
    faces = itf.face_inflow(source, s, centres=centres, outward_normals=normals,
                            areas=areas)
    assert np.all(faces.number_density[faces.theta > math.radians(30.0) + 1e-12] == 0.0)
    assert np.all(faces.temperature > 0.0)
    fraction = source.core_fraction(math.radians(30.0))
    balance = itf.flux_balance(source, faces, fraction)
    assert balance["net_mass_residual"] == pytest.approx(0.0, abs=2e-3)
    assert abs(balance["gross_minus_net_mass_fraction"]) < 1e-9
    assert balance["mdot_core_target_kg_per_s"] == pytest.approx(
        fraction * source.targets.mdot_kg_per_s, rel=1e-15)


def test_faces_oriented_away_from_the_flow_are_refused(source):
    centres, normals, areas = itf.ideal_sphere_faces(3.0 * source.r_a, 0.5 * math.pi,
                                                      n_theta=8, n_phi=16)
    with pytest.raises(ValueError, match="pointing back"):
        itf.face_inflow(source, settings(), centres=centres, outward_normals=-normals,
                        areas=areas)


def test_the_angular_law_mismatch_is_measured_not_tuned(source):
    out = itf.angular_shape_comparison(source, angles_deg=(30.0,))
    row = out["rows"][0]
    assert row["markelov_fraction"] == pytest.approx(source.core_fraction(math.radians(30)))
    assert row["cai_collisionless_fraction"] > row["markelov_fraction"]
    assert out["prandtl_meyer_limit_from_anchor_deg"] < out["markelov_limiting_angle_deg"]


# --------------------------------------------------------------------------- #
# configuration
# --------------------------------------------------------------------------- #

def test_malformed_interface_settings_are_refused():
    good = {"enabled": True, "purpose": "production", "bird_threshold": 0.05}
    itf.parse_interface(good)
    for bad, match in (({**good, "bird_treshold": 0.05}, "unknown"),
                       ({"purpose": "production"}, "enabled"),
                       ({**good, "bird_threshold": -1}, "positive"),
                       ({**good, "bird_threshold": True}, "number"),
                       ({**good, "purpose": "demo"}, "purpose"),
                       ({**good, "core_half_angle_deg": 120.0}, "<= 90"),
                       ({**good, "tails": "inject_exploratory"}, "different model"),
                       ({**good, "source_offset_over_D": [0.1, 0.0, 0.0]}, "CP-1"),
                       ({**good, "purpose": "development-smoke"}, "run_mode")):
        with pytest.raises(itf.InterfaceConfigError, match=match):
            itf.parse_interface(bad)


def test_provisional_values_are_for_smoke_runs_only():
    s = itf.parse_interface({"enabled": True, "purpose": "interface-validation",
                             "core_half_angle_deg": 30.0, "core_angle_provisional": True,
                             "source_validity_min_radius_over_D": 1.0})
    assert any("provisional" in m for m in s.missing_decisions())
    smoke = dataclasses.replace(s, purpose="development-smoke", run_mode="exploratory")
    assert smoke.scientific is False


def test_mode_selection(tmp_path):
    assert itf.load_interface(tmp_path) is None
    (tmp_path / itf.SETTINGS_FILE).write_text("interface:\n  enabled: false\n")
    assert itf.load_interface(tmp_path) is None
    (tmp_path / itf.SETTINGS_FILE).write_text("interface:\n  enabled: true\n")
    assert itf.load_interface(tmp_path).enabled
    (tmp_path / itf.SETTINGS_FILE).write_text("interfaces:\n  enabled: true\n")
    with pytest.raises(itf.InterfaceConfigError):
        itf.load_interface(tmp_path)
