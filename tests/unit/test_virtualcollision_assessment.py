"""Independent assessment semantics, skip policy, and uncertainty checks."""
import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from plumetools.virtualcollision.assessment import (
    aggregate_reference_bins, angular_ring_profiles, dsmc_reference_metrics, inventory,
    run_reference_assessment, time_block_standard_error, write_json,
)


def test_dsmc_reference_denominator_weights_and_negative_improvement():
    # Two physical volumes: reference squared norm = 1 + 3*4 = 13.
    result = dsmc_reference_metrics([1, 2], [2, 2], [3, 2], weights=[1, 3],
                                    common_mask=[True, True])
    assert result["baseline_error"] == pytest.approx(1 / np.sqrt(13))
    assert result["corrected_error"] == pytest.approx(2 / np.sqrt(13))
    assert result["relative_improvement"] == pytest.approx(-1.)
    assert result["corrected_signed_bias"] == pytest.approx(.5)


def test_failed_corrected_point_cannot_choose_comparison_mask():
    result = dsmc_reference_metrics([1, 1], [2, 2], [1, np.nan], weights=[1, 1],
                                    common_mask=[True, True])
    assert result["baseline_error"] == 1
    assert result["corrected_error"] is None
    assert result["requested_points"] == 2
    assert result["status"] == "failed-prediction-on-common-mask"


def test_zero_crossing_uses_fixed_scale_and_zero_baseline_has_no_improvement():
    result = dsmc_reference_metrics([-1, 0, 1], [-1, 0, 1], [-1, 1, 1],
                                    weights=[1, 1, 1], common_mask=[True] * 3, scale=2)
    assert result["corrected_error"] == pytest.approx(1 / np.sqrt(12))
    assert result["relative_improvement"] is None


def test_no_reference_points_is_explicit():
    result = dsmc_reference_metrics([1], [2], [3], weights=[1], common_mask=[False])
    assert result["status"] == "empty-common-mask"


def _case(tmp_path, folder, *, kn=1., law="VariableHardSphere", interface=False):
    path = tmp_path / "cases" / "cai2012-kn-extension" / folder
    path.mkdir(parents=True)
    config = {"nozzle": {"diameter_m": .2},
              "exit": {"knudsen": kn, "T0_K": 300, "speed_ratio": 2,
                       "characteristic_length": "diameter"},
              "gas": {"species_name": "Ar", "mass_kg": 6.63e-26,
                      "diameter_m": 4.17e-10, "omega": .74, "t_ref_K": 273.,
                      "mean_free_path_convention": "vhs", "internal_degrees_of_freedom": 0},
              "dsmc": {"binary_collision_model": law, "collisions_enabled": True}}
    (path / "case.yaml").write_text(yaml.safe_dump(config))
    if interface:
        (path / "interface.yaml").write_text("enabled: true\n")
    return path


def test_missing_and_unsupported_cases_are_explicit(tmp_path):
    _case(tmp_path, "Cases/missing")
    _case(tmp_path, "Cases/unsupported", law="VSS")
    _case(tmp_path, "Interface/source", interface=True)
    records = inventory(tmp_path)["cases"]
    missing = next(r for r in records if r["path"] and r["path"].endswith("missing"))
    bad = next(r for r in records if r["path"] and r["path"].endswith("unsupported"))
    interface = next(r for r in records if r["source_type"] == "interface")
    assert missing["quality"] == "unavailable"
    assert "missing-dsmc-results" in missing["reasons"]
    assert bad["quality"] == "unsupported"
    assert "unsupported-collision-law" in bad["reasons"]
    assert not interface["primary_included"]
    assert "interface-source-changes-boundary-value-problem" in interface["reasons"]
    extra = next(r for r in records if r["nominal_kn"] == .01)
    assert extra["additional_cai_family_comparison"]
    assert not extra["original_six_point_inventory"]


def test_known_dense_direct_exit_remains_nonproduction(tmp_path):
    path = _case(tmp_path, "Cases/dense", kn=.001)
    (path / "results").mkdir()
    (path / "results" / "centerline.csv").write_text("x,n\n1,2\n")
    result = next(r for r in inventory(tmp_path)["cases"] if r["path"] and r["path"].endswith("dense"))
    assert result["quality"] == "exploratory"
    assert not result["primary_included"]
    assert "known-inadequate-direct-exit-resolution-and-timestep" in result["reasons"]


def test_block_estimator_acknowledges_temporal_correlation():
    # Duplicated adjacent observations cannot provide eight independent values.
    data = np.repeat(np.array([-1, 1, -1, 1.] * 8), 8)
    result = time_block_standard_error(data)
    iid_se = data.std(ddof=1) / np.sqrt(len(data))
    assert result["standard_error"] > 2 * iid_se
    assert time_block_standard_error([1, 2])["standard_error"] is None


def test_json_preserves_missing_as_null_without_fake_errorbars(tmp_path):
    path = tmp_path / "report.json"
    write_json(path, {"missing": np.nan, "raw": np.array([1., np.inf])})
    assert json.loads(path.read_text()) == {"missing": None, "raw": [1., None]}


def test_reference_aggregation_preserves_between_cell_velocity_dispersion():
    from plumetools.cai2012.post import SampledField
    from plumetools.cai2012.constants import BOLTZMANN_J_PER_K
    mass = 6.63e-26
    sample = SampledField("1", np.array([[.4, .1, 0], [.6, .1, 0]]),
                          np.array([1., 3.]), np.array([2., 1.]),
                          np.array([[100., 0, 0], [-100., 0, 0]]), np.array([10., 20.]))
    result = aggregate_reference_bins(sample, [0, 1], [0, .5], mass)
    assert result["n"][0] == pytest.approx(1.25)
    assert result["U"][0, 0] == pytest.approx(-20.)
    # Independently: within-cell thermal energy gives 16 K; bulk-velocity
    # variance = E[vx²]-E[vx]² = 10000-400, adding to the trace temperature.
    assert result["T_trans"][0] == pytest.approx(16. + mass/BOLTZMANN_J_PER_K * 9600/3)
    assert result["T_trans"][0] > 16.
    assert result["tensor_status"].startswith("unavailable")


def test_interface_kn_does_not_satisfy_missing_direct_inventory(tmp_path):
    _case(tmp_path, "Interface/validation", kn=.01, interface=True)
    rows = inventory(tmp_path)["cases"]
    assert any(r["nominal_kn"] == .01 and r["source_type"] == "direct-nozzle"
               and r["kind"] == "missing" for r in rows)


def test_angular_same_ring_estimator_and_missing_bins_remain_missing():
    from plumetools.virtualcollision.source import DiskSource
    from plumetools.cai2012.constants import BOLTZMANN_J_PER_K
    source = DiskSource(.1, 1e18, 300., 2.)
    angle = np.radians([1., 45.])
    points = np.column_stack([2*np.cos(angle), 2*np.sin(angle), np.zeros(2)])
    n = np.array([1e15, 2e14])
    u = np.column_stack([100*np.cos(angle), 100*np.sin(angle), np.zeros(2)])
    data = {"n": n, "M1": n[:, None]*u,
            "M2_trace": n*(100**2+3*BOLTZMANN_J_PER_K/source.species.mass_kg*20)}
    profiles = angular_ring_profiles(points, np.full(2, .001), data, data, data, source)
    finite = [row for row in profiles if np.isfinite(row["DSMC_n_over_n0"])]
    assert finite
    assert len(finite) < len(profiles)
    for row in finite:
        assert row["DSMC_n_over_n0"] == row["Cai_n_over_n0"] == row["Virtual-collision model_n_over_n0"]
        assert row["DSMC_uR_sqrt_beta0"] == pytest.approx(100/source.thermal_speed)
        assert row["DSMC_T_over_T0"] == pytest.approx(20/300)


def test_full_reference_adapter_end_to_end_with_known_synthetic_field(tmp_path, monkeypatch):
    """Exercise API, raw output preservation and honest empty angular bins.

    The synthetic fixture is deliberately not scientific reference evidence.
    """
    from plumetools.cai2012 import post
    case = _case(tmp_path, "Cases/fixture", kn=10.)
    document = inventory(tmp_path)
    row = next(r for r in document["cases"] if r.get("path", "") and r["path"].endswith("fixture"))
    row["available_field_times"] = ["1"]
    row["available_geometry_times"] = ["0"]
    row["quality"] = "exploratory"
    points = np.array([[.05, .05, 0], [.05, .15, 0], [.15, .05, 0], [.15, .15, 0]])
    sampled = post.SampledField("1", points, np.full(4, 1e-4),
                               np.full(4, 1e16), np.tile([600., 0, 0], (4, 1)), np.full(4, 100.))
    monkeypatch.setattr(post, "read_case", lambda *a, **kw: sampled)
    config = {"_root": str(tmp_path),
              "smoke": {"kn": [10.], "samples": 8, "batches": 2, "hazard_steps": 2,
                        "partner_radial": 2, "partner_azimuthal": 4, "partner_speed": 8},
              "assessment_bins": {"x_count": 2, "r_count": 2, "x_max_over_D": 1., "r_max_over_D": 1.}}
    output = tmp_path / "synthetic-output"
    output.mkdir()
    import importlib.util
    can_plot = importlib.util.find_spec('matplotlib') is not None
    report = run_reference_assessment(config, document, output, smoke=True, plots=can_plot)
    assert report["processed_cases"] == 1
    assert report["primary_cases"] == 0
    target = output / row["path"].replace("/", "__")
    arrays = np.load(target / "raw.npz")
    assert arrays["VC_valid"].shape == (4,)
    assert arrays["VC_M2"].shape == (4, 3, 3)
    assert arrays["common_density_mask"].shape == (4,)
    scores = json.loads((target / "metrics.json").read_text())
    assert all(not score["headline_eligible"] for score in scores)
    angular = [score for score in scores if score["observable"].startswith("angular")]
    assert angular and all(score["status"] == "empty-common-mask" for score in angular)
    if can_plot:
        for name in ('density_three_way.png', 'density_errors.png',
                     'centerline_three_way.png', 'angular_three_way.png'):
            assert (target/name).is_file()
