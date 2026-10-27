"""Small numerical checks for empirical coordinate maps, not DSMC evidence."""

import numpy as np
import pytest

from plumetools.cai2012 import analytical
from plumetools.cai2012ext.scaling import CaiPredictor, CoordinateMap, kn_feature


POINTS = np.asarray([[.5, 0.], [1., .3], [2., 1.], [8., 2.]])


@pytest.mark.parametrize("family", ["baseline", "axial", "axial_transverse", "radial_angular"])
def test_identity_coordinates_and_density_are_exact_direct_evaluations(family):
    mapping = CoordinateMap(family)
    assert np.array_equal(mapping.map_points(POINTS, .01), POINTS)
    predictor = CaiPredictor(mapping, n_radial=24, n_azimuthal=48)
    expected = analytical.density_ratio(POINTS[:, 0], POINTS[:, 1], 0., .5, 2.,
                                        n_radial=24, n_azimuthal=48)
    assert np.array_equal(predictor.density(POINTS, .01), expected)
    assert np.array_equal(mapping.jacobian(POINTS, .01), np.ones(len(POINTS)))


@pytest.mark.parametrize("family,coefficients", [
    ("axial", [-.3]), ("axial_transverse", [-.3, .2]),
    ("radial_angular", [.1, -.1, .05, .3, -.1, .03])])
def test_high_kn_limit_is_exact_identity_and_large_kn_converges(family, coefficients):
    mapping = CoordinateMap(family, coefficients)
    assert np.array_equal(mapping.map_points(POINTS, np.inf), POINTS)
    np.testing.assert_allclose(mapping.map_points(POINTS, 1e12), POINTS, rtol=1e-12, atol=1e-12)


def test_affine_direction_positive_scales_and_volume_jacobian():
    # a_x=.5: displayed axial contours compress; W samples twice as far away.
    mapping = CoordinateMap("axial_transverse", [-2*np.log(2), 2*np.log(1.5)])
    np.testing.assert_allclose(mapping.map_points(POINTS, 1), POINTS/[.5, 1.5])
    scales = mapping.scale_factors(1)
    assert scales == pytest.approx({"axial": .5, "transverse": 1.5})
    assert np.all(mapping.jacobian(POINTS, 1) == pytest.approx(1/(.5*1.5**2)))


def test_radial_map_monotonicity_and_deterministic_quadrature_convergence():
    mapping = CoordinateMap("radial_angular", [.15, -.1, .04, .4, .1, -.03], quadrature_order=12)
    xi = np.geomspace(1e-4, 1e4, 80)
    theta = np.full_like(xi, .7)
    eta, derivative = mapping.radial_coordinates(xi, theta, .01)
    assert np.all(np.diff(eta) > 0)
    assert np.all(derivative > 0)
    assert np.all(derivative >= np.exp(-mapping.log_bound))
    points = np.column_stack((xi*np.cos(theta), xi*np.sin(theta)))
    check = mapping.quadrature_convergence(points, .01)
    assert check["q_vs_4q_max_relative"] < 1e-6
    assert check["2q_vs_4q_max_relative"] < 1e-11
    assert np.array_equal(mapping.map_points(points, .01), mapping.map_points(points, .01))


def test_radial_volume_jacobian_matches_independent_cartesian_differences():
    mapping = CoordinateMap("radial_angular", [.1, -.04, .02, .3, .05, -.02])
    point = np.asarray([1.2, .5, .7])
    step = 1e-5
    numerical = np.column_stack([
        (mapping.map_cartesian(point+step*np.eye(3)[i], .1)
         - mapping.map_cartesian(point-step*np.eye(3)[i], .1))/(2*step)
        for i in range(3)])
    expected = mapping.jacobian([[point[0], np.hypot(*point[1:])]], .1)[0]
    assert np.linalg.det(numerical) == pytest.approx(expected, rel=1e-8)


def test_origin_and_source_lip_are_reported_separately():
    angular_origin = CoordinateMap("radial_angular", [0., .1, 0., .2, 0., 0.])
    assert np.isnan(angular_origin.jacobian([[0., 0.]], 1)[0])
    assert not angular_origin.source_diagnostics(1)["origin_differentiable"]
    smooth_origin = CoordinateMap("radial_angular", [.1, 0., 0., .2, .1, 0.])
    assert smooth_origin.jacobian([[0., 0.]], 1)[0] == pytest.approx(np.exp(.15))
    diagnostics = smooth_origin.source_diagnostics(1)
    assert diagnostics["exit_lip_mapped"][0] == 0
    assert not diagnostics["emitting_disk_boundary_preserved"]
    assert diagnostics["lip_radial_displacement_over_D"] != 0


def test_azimuth_preserved_for_solver_axis_and_physical_units():
    mapping = CoordinateMap("axial_transverse", [.2, -.3])
    origin = np.asarray([2., 3., 4.])
    point = origin + np.asarray([.04, .08, .03])
    mapped = mapping.map_cartesian(point, .5, axis=1, origin=origin, diameter=.2)
    relative = mapped-origin
    assert relative[0]/relative[2] == pytest.approx(4/3)
    assert relative[1] == pytest.approx(.08/np.exp(.2/1.5))
    assert np.hypot(relative[0], relative[2]) == pytest.approx(.05/np.exp(-.3/1.5))
    assert np.array_equal(CoordinateMap().map_cartesian(point, .5, axis=1,
                                                      origin=origin, diameter=.2), point)


def test_global_bounds_check_polynomial_interior_and_angular_extrema():
    # c1*h+c2*h^2 is zero at both endpoints but reaches one in the interior.
    with pytest.raises(ValueError, match="log_bound"):
        CoordinateMap("axial", [4., -4.], degree=2, log_bound=.9)
    mapping = CoordinateMap("axial", [4., -4.], degree=2, log_bound=1.)
    assert mapping.bound_diagnostics()["bound_active"]
    # cos2 - cos4 has its largest value at cos2=1/4, away from angle endpoints.
    with pytest.raises(ValueError, match="log_bound"):
        CoordinateMap("radial_angular", [0., 1., -1., 0., 0., 0.], log_bound=1.01)
    with pytest.raises(ValueError, match="log_bound"):
        CoordinateMap("axial_transverse", [0., 3.])


def test_temperature_transfer_uses_same_density_map_without_refitting():
    mapping = CoordinateMap("axial_transverse", [-.2, .1])
    predictor = CaiPredictor(mapping, n_radial=24, n_azimuthal=48)
    original = mapping.coefficients.copy()
    mapped = mapping.map_points(POINTS, 1.)
    expected = analytical.moments(mapped[:, 0], mapped[:, 1], 0., .5, 2.,
                                 n_radial=24, n_azimuthal=48)["T_over_T0"]
    assert np.array_equal(predictor.temperature(POINTS, 1), expected)
    assert np.array_equal(mapping.coefficients, original)
    assert predictor.mapping is mapping


def test_serialized_predictor_reproduces_new_point_predictions(tmp_path):
    predictor = CaiPredictor(CoordinateMap("radial_angular", [.1, 0., -.02, .3, -.05, .01]),
                             speed_ratio=2.3, n_radial=20, n_azimuthal=40, chunk_size=3)
    path = tmp_path / "predictor.json"
    predictor.save(path)
    restored = CaiPredictor.load(path)
    independent = np.asarray([[.8, .12], [3.5, .9], [7., 3.], [11., 2.]])
    assert restored.to_dict() == predictor.to_dict()
    assert np.array_equal(restored.density(independent, .3), predictor.density(independent, .3))
    assert np.array_equal(restored.temperature(independent, .3), predictor.temperature(independent, .3))


def test_invalid_physical_and_mapped_coordinates_are_not_silently_dropped():
    predictor = CaiPredictor(CoordinateMap())
    with pytest.raises(ValueError, match="source-plane"):
        predictor.density([[0., .5], [1., .2]], 1)
    with pytest.raises(ValueError, match="downstream"):
        predictor.density([[-.1, .2]], 1)
    with pytest.raises(ValueError, match="finite"):
        predictor.density([[np.nan, .2]], 1)
    with pytest.raises(ValueError, match="Kn"):
        kn_feature(-1)


def test_alternative_feature_is_explicit_and_serialized():
    mapping = CoordinateMap("axial", [.2], feature="inverse_kn_squared")
    assert mapping.map_points([[1., .2]], 1)[0, 0] == pytest.approx(np.exp(-.05))
    assert CoordinateMap.from_dict(mapping.to_dict()).feature == "inverse_kn_squared"


def test_direct_disk_quadrature_caps_are_flagged_without_hiding_near_exit_points():
    predictor = CaiPredictor(CoordinateMap(), n_radial=12, n_azimuthal=24)
    diagnostics = predictor.evaluation_diagnostics([[.0001, .5], [.5, .2], [0., .5]], 1)
    assert diagnostics["n_points"] == 3
    assert diagnostics["n_adaptive_disk_order_capped"] == 1
    assert diagnostics["n_undefined_source_plane"] == 1
    check = predictor.quadrature_convergence(POINTS, 1)
    assert check["q_vs_4q_max_relative"] < 1e-5
    assert check["2q_vs_4q_max_relative"] < 1e-10


def test_large_radial_batch_agrees_with_independent_point_evaluations():
    mapping = CoordinateMap("radial_angular", [.1, -.04, .02, .2, .03, -.01])
    radii = np.geomspace(.01, 100., 50000)
    points = np.column_stack((radii*np.cos(.6), radii*np.sin(.6)))
    full = mapping.map_points(points, .1)
    selected = [0, 19001, 41665, 41666, 49999]
    independent = np.vstack([mapping.map_points(points[index], .1) for index in selected])
    np.testing.assert_allclose(full[selected], independent, rtol=1e-14, atol=1e-14)
    assert np.all(np.diff(np.linalg.norm(full, axis=1)) > 0)
