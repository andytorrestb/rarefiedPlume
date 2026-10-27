"""Small synthetic numerical checks; none of these fields is DSMC evidence."""

import copy

import numpy as np
import pytest

from plumetools.cai2012ext.scaling import CaiPredictor, CoordinateMap
from plumetools.cai2012ext.scaling_fit import (
    ScalingCase, block_refits, contour_diagnostics, evaluate, extract_crossing, fit,
    prepare_case, select_cases,
)


class SyntheticEvaluator:
    """Known smooth field independent of the Cai quadrature implementation."""

    def __init__(self, mapping, speed_ratio=2., **kwargs):
        self.mapping = mapping

    def density(self, points, kn):
        x, p = self.mapping.map_points(points, kn).T
        return 1. / (1. + x*x + 2*p*p)

    def temperature(self, points, kn):
        x, p = self.mapping.map_points(points, kn).T
        return .6 / (1. + np.hypot(x, p))


ANGLES = np.deg2rad([4., 19., 39., 64., 84.])
EDGES = np.deg2rad([0., 10., 30., 50., 75., 90.])


def ray_points(n=18):
    radius, angle = np.meshgrid(np.geomspace(.12, 100., n), ANGLES)
    return np.column_stack([(radius*np.cos(angle)).ravel(),
                            (radius*np.sin(angle)).ravel()])


def synthetic_cases(mapping, *, points=None, with_temperature=False):
    points = ray_points() if points is None else points
    cases = []
    for i, kn in enumerate([.1, 1., 10.]):
        predictor = SyntheticEvaluator(mapping)
        cases.append(ScalingCase(
            f"synthetic-Kn{kn:g}", kn, points,
            predictor.density(points, kn), np.linspace(1., 7., len(points)),
            temperature=predictor.temperature(points, kn) if with_temperature else None,
            suitable=i != 1, metadata={"source": "synthetic numerical verification"}))
    return cases


def test_fixed_full_domain_mask_and_coverage_without_lip_cut():
    points = np.array([[.01, .5], [0., .4], [-.1, .2], [.3, .4],
                       [.4, .1], [.5, .1], [.6, .1], [.7, .1]])
    density = np.array([.2, .2, .2, np.nan, .2, .0, .0, .3])
    volumes = np.ones(8)
    volumes[4] = 0
    occupied = np.array([True, True, True, True, True, False, True, True])
    case = ScalingCase("synthetic-mask", 1., points, density, volumes, occupied=occupied,
                       metadata={"domain_valid": [1, 1, 1, 1, 1, 1, 1, 0]})
    prepared = prepare_case(case, evaluator_factory=SyntheticEvaluator)
    assert prepared.mask.tolist() == [True, False, False, False, False, False, True, False]
    assert prepared.coverage["near_lip_samples"] == 1
    assert prepared.coverage["lip_excluded"] == 0
    assert prepared.coverage["resolved_zero_samples"] == 1
    assert sum(prepared.coverage["excluded"].values()) + prepared.mask.sum() == 8
    again = prepare_case(case, floor=.25, evaluator_factory=SyntheticEvaluator)
    np.testing.assert_array_equal(again.mask, prepared.mask)
    np.testing.assert_array_equal(again.density_bands, prepared.density_bands)
    assert np.isnan(case.density[3])  # raw fields were not overwritten


def test_equal_occupied_angle_band_weights_and_volume_weights_within_stratum():
    base = ray_points(n=12)
    # Duplicate one physical sample with different real volume. Its stratum's
    # total influence must remain the same, irrespective of this cell count.
    points = np.concatenate([base, base[:1]])
    volumes = np.ones(len(points))
    volumes[-1] = 17
    case = ScalingCase("synthetic-weighting", 1., points,
                       SyntheticEvaluator(CoordinateMap()).density(points, 1.), volumes)
    p = prepare_case(case, angle_edges=EDGES, evaluator_factory=SyntheticEvaluator)
    totals = np.array([s["weight"] for s in p.coverage["strata"]])
    assert np.allclose(totals, 1 / len(totals))
    assert np.isclose(p.weights.sum(), 1.)
    assert p.weights[-1] == pytest.approx(17 * p.weights[0])
    np.testing.assert_array_equal(p.density_bands,
                                  np.searchsorted([.001, .01, .1], p.baseline, side="right"))


@pytest.mark.parametrize("r,d,q,status", [
    ([1, 2, 3], [.5, .05, .005], .1, "ok"),
    ([1, 2, 3], [.5, .4, .3], .1, "domain_truncated_outer"),
    ([1, 2, 3], [.05, .04, .03], .1, "domain_truncated_inner"),
    ([1, 2, 3, 4], [.5, .05, .3, .01], .1, "multiple"),
    ([1, 2, 3, 4], [.5, .1, .1, .01], .1, "ambiguous"),
    ([1], [.5], .1, "missing"),
    ([1, 2, 3], [.01, .1, .5], .1, "ambiguous_increasing"),
    ([1, 2, 3], [.5, .3, .1], .1, "domain_truncated_outer"),
])
def test_contour_branch_flags(r, d, q, status):
    result = extract_crossing(r, d, q)
    assert result["status"] == status
    assert (result["radius"] is not None) == (status == "ok")


def test_contour_crossings_do_not_bridge_missing_radial_shells():
    result = extract_crossing([1., 2., 4.], [.5, .02, .001], .1,
                               gaps_after=[True, False])
    assert result["status"] == "missing_gap"
    assert result["radius"] is None


@pytest.mark.parametrize("family,coefficients", [
    ("axial", [-.3]), ("axial_transverse", [-.3, .15]),
    ("radial_angular", [.1, -.07, .04, -.12, .08, .03]),
])
def test_joint_known_map_recovery_and_independent_points(family, coefficients):
    truth = CoordinateMap(family, coefficients)
    cases = synthetic_cases(truth)
    result = fit(cases, CoordinateMap(family), evaluator_factory=SyntheticEvaluator,
                 angle_edges=EDGES, radial_bins=18, max_nfev=60,
                 evaluate_temperature=False)
    assert result.success, result.message
    assert result.identifiable
    np.testing.assert_allclose(result.mapping.coefficients, coefficients, atol=2e-5)
    independent = np.array([[.16, .51], [2.1, 4.3], [9.1, 3.8], [35., 45.]])
    expected = SyntheticEvaluator(truth).density(independent, .4)
    actual = SyntheticEvaluator(result.mapping).density(independent, .4)
    np.testing.assert_allclose(actual, expected, rtol=2e-6)
    assert result.training_cases == [c.name for c in cases]
    assert not result.diagnostics["predictive_validation"]
    metrics = result.diagnostics["training_metrics"]
    assert metrics["L_F"] < 1e-12
    assert metrics["L_C"] < 1e-12


def test_contour_track_separate_objective_and_both_diagnostics_reported():
    truth = CoordinateMap("axial_transverse", [-.3, .2])
    cases = synthetic_cases(truth)
    result = fit(cases, CoordinateMap("axial_transverse"), track="C",
                 evaluator_factory=SyntheticEvaluator, angle_edges=EDGES,
                 radial_bins=18, max_nfev=50, evaluate_temperature=False)
    assert result.success, result.message
    assert result.track == "C"
    np.testing.assert_allclose(result.mapping.coefficients, truth.coefficients, atol=1e-4)
    assert result.diagnostics["training_metrics"]["L_F"] < 1e-10
    assert result.diagnostics["training_metrics"]["L_C"] < 1e-10


def test_optional_whole_kn_holdout_before_any_learned_preprocessing():
    truth = CoordinateMap("axial", [-.2])
    cases = synthetic_cases(truth)
    # A repeated Kn at different mesh/particle level must also be held out.
    duplicate = copy.copy(cases[1])
    duplicate.name = "synthetic-Kn1-second-level"
    cases.append(duplicate)

    class RejectHeldKn(SyntheticEvaluator):
        def density(self, points, kn):
            assert kn != 1., "held-out DSMC reached preprocessing or tuning"
            return super().density(points, kn)

    result = fit(cases, CoordinateMap("axial"), holdout_kn=[1.],
                 evaluator_factory=RejectHeldKn, angle_edges=EDGES, radial_bins=18,
                 max_nfev=20, evaluate_temperature=False)
    assert result.success
    assert result.holdout_cases == [cases[1].name, duplicate.name]
    assert {row["regime"] for row in result.diagnostics["holdout"]} == {"interpolation"}
    assert not any("Kn1" == name for name in result.training_cases)
    assert result.diagnostics["preprocessing_policy"].startswith("whole-Kn split before")


def test_quality_filtered_and_exploratory_selection_and_empty_subset():
    cases = synthetic_cases(CoordinateMap())
    train, held, excluded = select_cases(cases)
    assert len(train) == 3 and not held and not excluded
    train, held, excluded = select_cases(cases, suitable_only=True)
    assert len(train) == 2 and excluded == [cases[1]]
    for case in cases:
        case.suitable = False
    result = fit(cases, CoordinateMap("axial"), suitable_only=True,
                 evaluator_factory=SyntheticEvaluator)
    assert not result.success and not result.identifiable
    assert result.training_cases == []
    assert result.diagnostics["training_metrics"]["L_F"] is None


def test_joint_kn_only_fit_refuses_incompatible_physics():
    cases = synthetic_cases(CoordinateMap())
    cases[1].speed_ratio = 3.
    with pytest.raises(ValueError, match="speed ratio"):
        fit(cases, CoordinateMap("axial"), evaluator_factory=SyntheticEvaluator)
    cases[1].speed_ratio = 2.
    cases[0].metadata["stratum"] = "synthetic-gas-law-1"
    cases[1].metadata["stratum"] = "synthetic-gas-law-2"
    with pytest.raises(ValueError, match="strata"):
        fit(cases, CoordinateMap("axial"), evaluator_factory=SyntheticEvaluator)


def test_temperature_transfer_keeps_density_coefficients_and_uses_own_mask():
    mapping = CoordinateMap("axial", [-.25])
    cases = synthetic_cases(mapping, with_temperature=True)
    cases[0].temperature[0] = np.nan
    coefficients = mapping.coefficients.copy()
    report = evaluate(cases, mapping, evaluator_factory=SyntheticEvaluator,
                       angle_edges=EDGES, radial_bins=18)
    np.testing.assert_array_equal(mapping.coefficients, coefficients)
    for case in report["cases"]:
        transfer = case["temperature_transfer"]
        assert transfer["available"]
        assert not transfer["coefficients_refitted"]
        assert transfer["metrics"]["rmse"] == 0
        assert case["uncertainty"]["statistical"] == "unavailable"
        assert not case["uncertainty"]["occupancy_is_confidence_interval"]
    assert report["cases"][0]["temperature_transfer"]["excluded_temperature_samples"] == 1


def test_model_invalid_values_cannot_improve_score_by_disappearing():
    cases = synthetic_cases(CoordinateMap())

    class InvalidMappedEvaluator(SyntheticEvaluator):
        def density(self, points, kn):
            values = super().density(points, kn)
            if self.mapping.family != "baseline":
                values[0] = np.nan
            return values

    report = evaluate(cases, CoordinateMap("axial", [.1]),
                      evaluator_factory=InvalidMappedEvaluator, angle_edges=EDGES,
                      radial_bins=18, temperature=False)
    assert report["L_F"] == 1e6
    assert all(not c["model_evaluation_valid"] for c in report["cases"])
    assert all(c["coverage"]["valid_samples"] == len(cases[0].density) for c in report["cases"])


def test_empty_results_and_unidentifiable_contour_objective():
    assert not fit([], CoordinateMap("axial")).identifiable
    points = np.array([[10., 0.], [11., 0.]])
    case = ScalingCase("synthetic-missing-contours", 1, points, [.7, .6], [1, 1])
    result = fit([case], CoordinateMap("axial"), track="C",
                 evaluator_factory=SyntheticEvaluator, radial_bins=2,
                 evaluate_temperature=False)
    assert not result.identifiable
    assert result.diagnostics["omitted_no_computable_objective"] == [case.name]
    assert result.diagnostics["training_metrics"]["L_C"] is None


def test_baseline_has_no_coefficients_and_is_identifiable_on_computable_data():
    result = fit(synthetic_cases(CoordinateMap()), CoordinateMap(),
                 evaluator_factory=SyntheticEvaluator, angle_edges=EDGES,
                 radial_bins=18, evaluate_temperature=False)
    assert result.success and result.identifiable
    assert result.mapping.parameter_count == 0
    assert result.diagnostics["training_metrics"]["L_F"] == 0


def test_actual_cai_joint_fit_recovers_coordinate_map_on_independent_points():
    options = {"n_radial": 12, "n_azimuthal": 24, "chunk_size": 32}
    points = np.array([[1., 0.], [1., .4], [2., .3], [3., 1.], [8., 3.], [15., 2.]])
    truth = CoordinateMap("axial", [-.18])
    cases = [ScalingCase(f"synthetic-Cai-Kn{k}", k, points,
                          CaiPredictor(truth, **options).density(points, k), np.ones(len(points)))
             for k in [.1, 10.]]
    result = fit(cases, CoordinateMap("axial"), max_nfev=20,
                 evaluator_options=options, evaluate_temperature=False)
    assert result.success, result.message
    assert result.mapping.coefficients[0] == pytest.approx(-.18, abs=2e-6)
    test_points = np.array([[1.4, .3], [7., .8], [20., 5.]])
    np.testing.assert_allclose(CaiPredictor(result.mapping, **options).density(test_points, 1.),
                               CaiPredictor(truth, **options).density(test_points, 1.), rtol=1e-6)


def test_block_refits_require_temporal_independence_and_report_parameter_spread():
    truth = CoordinateMap("axial", [-.2])
    cases = synthetic_cases(truth)
    for case in cases:
        case.metadata["density_blocks"] = np.stack([case.density*.99, case.density*1.01])
    unavailable = block_refits(cases, CoordinateMap("axial"), n_resamples=2,
                               evaluator_factory=SyntheticEvaluator, radial_bins=18,
                               angle_edges=EDGES)
    assert not unavailable["available"]
    assert unavailable["resamples_run"] == 0
    for case in cases:
        case.metadata.update(independent_blocks=True,
                             temporal_correlation_method="synthetic declared test correlation",
                             block_duration_s=10., integrated_autocorrelation_time_s=1.)
    report = block_refits(cases, CoordinateMap("axial"), n_resamples=3,
                          max_nfev=20, seed=7, evaluator_factory=SyntheticEvaluator,
                          radial_bins=18, angle_edges=EDGES)
    assert report["available"]
    assert report["successful_refits"] == 3
    assert report["coefficient_sample_std"][0] > 0
    assert not report["occupancy_is_confidence_interval"]
    assert "not a confidence interval" in report["interpretation"]
    cases[0].metadata["block_duration_s"] = 1.
    assert not block_refits(cases, CoordinateMap("axial"), n_resamples=2)["available"]


def test_model_evaluation_exception_is_reported_with_fixed_coverage():
    class FailingEvaluator(SyntheticEvaluator):
        def density(self, points, kn):
            if self.mapping.family != "baseline":
                raise ValueError("synthetic unavailable mapped field")
            return super().density(points, kn)

    report = evaluate(synthetic_cases(CoordinateMap()), CoordinateMap("axial", [.1]),
                       evaluator_factory=FailingEvaluator, temperature=False)
    assert report["L_F"] == 1e6
    assert all(c["model_evaluation_failure"] for c in report["cases"])
    assert all(c["linear"]["normalized_L2"] is None for c in report["cases"])


def test_independent_velocity_and_hybrid_flux_are_explicitly_labeled():
    from plumetools.cai2012 import analytical

    options = {"n_radial": 8, "n_azimuthal": 16}
    theta = np.deg2rad([4., 19., 39., 64.])
    points = np.column_stack((10*np.cos(theta), 10*np.sin(theta)))
    baseline = analytical.moments(points[:, 0], points[:, 1], 0., .5, 2., **options)
    case = ScalingCase("synthetic-velocity", 1, points, baseline["n_over_n0"],
                       np.ones(len(points)), metadata={
                           "cylindrical_velocity_sqrt_beta0": baseline["U_sqrt_beta0"][:, :2]})
    report = evaluate([case], CoordinateMap("axial", [-.1]), evaluator_options=options,
                       temperature=False, angle_edges=EDGES)
    physics = report["cases"][0]["independent_physics"]
    assert physics["available"] and not physics["velocity_corrected"]
    assert physics["axial_velocity"]["rmse"] == 0
    assert physics["angular_flux"]["available"]
    assert not physics["angular_flux"]["hybrid_is_conservative_prediction"]
    assert "unchanged Cai velocity" in physics["angular_flux"]["hybrid_definition"]
