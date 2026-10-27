"""Bounded orchestration checks; synthetic data are never scientific evidence."""

import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from plumetools.cai2012ext import scaling_fit as sf
from plumetools.cai2012ext import scaling_workflow as workflow
from plumetools.cai2012ext.scaling import CaiPredictor, CoordinateMap


FAMILY = Path(__file__).resolve().parents[2] / "cases" / "cai2012-kn-extension"


def configuration(tmp_path):
    cfg = yaml.safe_load((FAMILY / "scaling.yaml").read_text(encoding="utf-8"))
    cfg["resources"].update(smoke_max_points_per_case=24, smoke_max_nfev=3)
    cfg["fit"].update(families=["baseline"], tracks=["F"], subsets=["exploratory", "suitable"])
    cfg["evaluation"].update(n_radial=12, n_azimuthal=24)
    path = tmp_path / "scaling.yaml"
    path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    return cfg, path


def case(name, kn=1., suitable=True, stratum="same-source"):
    points = np.asarray([[.5, .1], [1., .2], [2., 1.], [5., 2.], [12., 3.]])
    predictor = CaiPredictor(CoordinateMap(), n_radial=12, n_azimuthal=24)
    return sf.ScalingCase(name=name, kn=kn, points=points,
                          density=predictor.density(points, kn), volumes=np.ones(len(points)),
                          temperature=predictor.temperature(points, kn), suitable=suitable,
                          metadata={"stratum": stratum, "synthetic": True})


def test_default_cli_only_inventories_and_skips_missing_results(tmp_path, monkeypatch):
    _, path = configuration(tmp_path)
    inventory = {"cases": [{"name": "missing", "Kn": .01, "compatible": True,
                             "fit_ready": False, "result_status": "missing", "warnings": ["missing field"]}],
                 "missing_results": [.01]}
    monkeypatch.setattr(workflow.si, "inventory", lambda *_: inventory)
    monkeypatch.setattr(workflow.sf, "fit", lambda *_args, **_kwargs: pytest.fail("inventory must not fit"))
    out = tmp_path / "output"
    assert workflow.main(["--config", str(path), "--output", str(out)], study_root=tmp_path) == 0
    document = json.loads((out / "inventory.json").read_text(encoding="utf-8"))
    assert document["missing_results"] == [.01]
    assert not (out / "assessment.json").exists()


def test_production_fit_requires_explicit_full_opt_in(tmp_path, monkeypatch):
    _, path = configuration(tmp_path)
    monkeypatch.setattr(workflow.si, "inventory", lambda *_: {"cases": []})
    monkeypatch.setattr(workflow, "load_available", lambda *_: pytest.fail("must refuse before loading fields"))
    with pytest.raises(SystemExit) as error:
        workflow.main(["fit", "--config", str(path), "--output", str(tmp_path / "out")], study_root=tmp_path)
    assert error.value.code == 2


def test_local_core_cap_cannot_be_exceeded(tmp_path):
    cfg, path = configuration(tmp_path)
    cfg["resources"]["cores"] = 9
    path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    with pytest.raises(ValueError, match="eight"):
        workflow.load_config(tmp_path, path)


def test_exploratory_suitable_and_source_strata_remain_separate(tmp_path, monkeypatch):
    cfg, _ = configuration(tmp_path)
    cfg["fit"]["tracks"] = ["F", "C"]
    observed = []

    def fake_fit(cases, mapping, *, track, **kwargs):
        observed.append((tuple(c.name for c in cases), track))
        return sf.FitResult(mapping, track, False, "Synthetic selection test", [], [], False, {}, {})

    monkeypatch.setattr(workflow.sf, "fit", fake_fit)
    cases = [case("adequate", suitable=True), case("questionable", suitable=False),
             case("different-gas", suitable=True, stratum="different-gas")]
    report = workflow.run_assessment(cases, cfg, tmp_path / "report", smoke=True)
    assert len(report["fits"]) == 8
    for names, _ in observed:
        assert not ("different-gas" in names and len(names) > 1)
    assert (("adequate", "questionable"), "F") in observed
    assert (("adequate",), "F") in observed
    assert report["synthetic"] is True
    assert report["production_dsmc_executed"] is False


def test_no_available_fields_produces_unidentifiable_report_without_predictor(tmp_path):
    cfg, _ = configuration(tmp_path)
    out = tmp_path / "report"
    report = workflow.run_assessment([], cfg, out, inv={"missing_results": [100., .1]})
    assert len(report["fits"]) == 1
    assert report["fits"][0]["identifiable"] is False
    assert not list(out.rglob("predictor.json"))
    text = (out / "README.md").read_text(encoding="utf-8")
    assert "100.0" in text and "0.1" in text
    assert "cannot establish model effectiveness" in text


def test_bounded_synthetic_pipeline_serializes_predictor_and_cross_diagnostics(tmp_path):
    cfg, _ = configuration(tmp_path)
    cases = [case("SYNTHETIC_A", kn=1.), case("SYNTHETIC_B", kn=10.)]
    out = tmp_path / "synthetic"
    report = workflow.run_assessment(cases, cfg, out, smoke=True)
    assert report["synthetic"] is True
    assert len(report["fits"]) == 2
    record = report["fits"][0]
    assert record["identifiable"] is True
    assert record["evaluation"]["L_F"] < 1e-20
    assert "L_C" in record["evaluation"]
    assert all(row["temperature_transfer"]["coefficients_refitted"] is False
               for row in record["evaluation"]["cases"])
    stored = CaiPredictor.load(out / record["predictor"])
    assert stored.mapping.family == "baseline"
    assert (out / "coefficients.csv").is_file()
    assert "SYNTHETIC" in (out / "README.md").read_text(encoding="utf-8")


def test_predictor_evaluate_cli_needs_no_dsmc_and_preserves_coordinates(tmp_path):
    _, path = configuration(tmp_path)
    predictor = CaiPredictor(CoordinateMap("axial", [-.2]), n_radial=12, n_azimuthal=24)
    model = tmp_path / "predictor.json"
    predictor.save(model)
    points = tmp_path / "points.csv"
    points.write_text("X,P\n1.0,0.2\n", encoding="utf-8")
    output = tmp_path / "predictions"
    assert workflow.main(["evaluate", "--config", str(path), "--output", str(output),
                          "--predictor", str(model), "--points", str(points), "--kn", ".1"],
                         study_root=tmp_path) == 0
    table = np.genfromtxt(output / "predictions.csv", delimiter=",", names=True)
    assert float(table["X"]) == 1.
    assert float(table["P"]) == .2
    assert float(table["n_over_ne"]) == pytest.approx(predictor.density([[1., .2]], .1)[0])


def test_strict_json_encodes_unavailable_diagnostics_as_null(tmp_path):
    path = tmp_path / "diagnostics.json"
    workflow.write_json(path, {"unavailable": np.nan, "infinite": np.inf,
                               "count": np.int64(2), "flags": np.asarray([True, False])})
    assert json.loads(path.read_text(encoding="utf-8")) == {
        "unavailable": None, "infinite": None, "count": 2, "flags": [True, False]}


def test_whole_kn_holdout_evaluation_stays_separate_from_training(tmp_path):
    cfg, _ = configuration(tmp_path)
    cfg["fit"].update(subsets=["exploratory"], holdout_kn=[10.])
    train = case("training_Kn1", kn=1.)
    held_a = case("held_Kn10_A", kn=10.)
    held_b = case("held_Kn10_B", kn=10.)
    held_a.density *= 1.5
    held_b.density *= 1.5
    report = workflow.run_assessment([train, held_a, held_b], cfg, tmp_path / "held", smoke=True)
    record = report["fits"][0]
    assert record["training_cases"] == ["training_Kn1"]
    assert record["holdout_cases"] == ["held_Kn10_A", "held_Kn10_B"]
    assert record["evaluation_role"] == "training descriptive agreement"
    assert [row["name"] for row in record["evaluation"]["cases"]] == ["training_Kn1"]
    assert [row["name"] for row in record["holdout_evaluation"]["cases"]] == ["held_Kn10_A", "held_Kn10_B"]
    assert record["evaluation"]["L_F"] < 1e-20
    assert record["holdout_evaluation"]["L_F"] == pytest.approx(np.log(1.5)**2)
    assert all(row["regime"] == "extrapolation" for row in record["diagnostics"]["holdout"])


def test_total_input_memory_cap_aborts_without_silently_dropping_case(tmp_path, monkeypatch):
    cfg, _ = configuration(tmp_path)
    cfg["resources"]["max_total_input_bytes"] = 1
    raw = case("too-large")
    monkeypatch.setattr(workflow.si, "load_case", lambda *_args, **_kwargs: raw)
    inv = {"cases": [{"name": raw.name, "compatible": True, "fit_ready": True, "stratum": "same"}]}
    with pytest.raises(MemoryError, match="Total input memory cap"):
        workflow.load_available(inv, cfg)


def test_tiny_synthetic_report_renders_contour_overlays_and_paired_maps(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    from plumetools.cai2012ext import scaling_report

    cfg, _ = configuration(tmp_path)
    cfg["fit"]["subsets"] = ["exploratory"]
    cfg["resources"]["smoke_max_points_per_case"] = 48
    cases = workflow.synthetic_cases(cfg)[:1]
    output = tmp_path / "rendered"
    report = workflow.run_assessment(cases, cfg, output, smoke=True)
    monkeypatch.setattr(workflow.sf, "fit", lambda *_args, **_kwargs: pytest.fail("report must not refit"))
    scaling_report.render_report(output, report, cases=cases)
    figures = output / "figures" / cases[0].name / "exploratory" / "F"
    for name in ("cartesian-density.png", "cartesian-residual.png", "source-centered-density.png",
                 "source-centered-residual.png", "density-profiles.png", "contour-displacement.png",
                 "contour-overlays.png"):
        assert (figures / name).stat().st_size > 1000
    assert (output / "mapping-versus-kn.png").stat().st_size > 1000
    diagnostic = json.loads((figures / "baseline-source-consistency.json").read_text(encoding="utf-8"))
    assert diagnostic["source"]["emitting_disk_boundary_preserved"] is True
    assert diagnostic["evaluator"]["direct_evaluation"] is True
    assert "q_vs_4q_max_relative" in diagnostic["map_quadrature"]
    assert "q_vs_4q_max_relative" in diagnostic["disk_quadrature"]
    rows = np.genfromtxt(figures / "contour-displacement.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    assert len(rows) == 3*(len(cfg["fit"]["angle_edges_deg"])-1)
    assert np.any(rows["status"] == "ok")


def test_report_records_empty_case_coverage_and_continues(tmp_path, monkeypatch):
    from plumetools.cai2012ext import figures, scaling_report

    cfg, _ = configuration(tmp_path)
    cfg["fit"]["subsets"] = ["exploratory"]
    empty = case("SYNTHETIC_no_estimators")
    empty.occupied[:] = False
    cases = [case("SYNTHETIC_valid"), empty]
    output = tmp_path / "empty-report"
    report = workflow.run_assessment(cases, cfg, output, smoke=True)
    assert report["fits"][0]["identifiable"] is True
    monkeypatch.setattr(figures, "_pyplot", lambda: None)
    scaling_report.render_report(output, report, cases=cases)
    coverage = json.loads((output / "figures" / empty.name / "empty-coverage.json").read_text(encoding="utf-8"))
    assert coverage["valid_samples"] == 0
    assert coverage["excluded"]["no_meaningful_dsmc_estimator"] == len(empty.density)
    assert not list((output / "figures" / empty.name).rglob("*.png"))
    assert (output / "figures" / "SYNTHETIC_valid" / "exploratory" / "F" / "baseline-field.csv").is_file()
