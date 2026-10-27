"""Registration cannot import another study, fabricate volumes, or relax guards."""

import dataclasses
from pathlib import Path
import shutil

import numpy as np
import pytest
import yaml

from plumetools.cai2012 import inflow
from plumetools.cai2012.config import load_case_config
from plumetools.cai2012.geometry import from_config
from plumetools.cai2012ext import policy, provenance
from plumetools.cai2012ext import scaling_inventory as inventory


REPO = Path(__file__).resolve().parents[2]
FAMILY = REPO / "cases" / "cai2012-kn-extension"


@pytest.fixture
def family(tmp_path, monkeypatch):
    root = tmp_path / "cai2012-kn-extension"
    (root / "baseCase").mkdir(parents=True)
    shutil.copyfile(FAMILY / "baseCase" / "case.yaml", root / "baseCase" / "case.yaml")
    # Preserve regression coverage of the original four-case default policy,
    # independently of which explicit profile the shipped study now selects.
    legacy = yaml.safe_load((FAMILY / "study.yaml").read_text(encoding="utf-8"))
    legacy.setdefault("meta", {}).pop("knudsen_suite", None)
    legacy["kn_cases"] = [row for row in legacy["kn_cases"]
                          if float(row["Kn"]) in (10., 1., .001, .0001)]
    legacy.update(cases_dir="Cases", manifest_name="manifest.yaml")
    (root / "study.yaml").write_text(yaml.safe_dump(legacy), encoding="utf-8")

    def bounded_assessment(cfg):
        return provenance.exit_conditions(cfg, from_config(cfg), inflow.from_config(cfg)), {
            "production_ready": True, "status": "production-candidate", "unresolved": []}
    monkeypatch.setattr(inventory, "_assess_config", bounded_assessment)
    return root


def register(root, cases=None, **kwargs):
    value = {"scope": inventory.REGISTRATION_SCOPE, "enabled": True,
             "purpose": "analysis-and-specification-only", "allow_paper_decades": True,
             "prepared_decades": list(inventory.TARGET_DECADES), "cases": cases or []}
    value.update(kwargs)
    path = root / "registration.yaml"
    path.write_text(yaml.safe_dump(value), encoding="utf-8")
    return path


def archive(root, name="archive", kn=1.0, speed_ratio=2.0):
    directory = root / name
    directory.mkdir()
    value = yaml.safe_load((root / "baseCase" / "case.yaml").read_text(encoding="utf-8"))
    value["exit"].update(knudsen=kn, speed_ratio=speed_ratio)
    (directory / "case.yaml").write_text(yaml.safe_dump(value), encoding="utf-8")
    return {"name": name, "path": name, "Kn": kn}


def full_study(root):
    value = yaml.safe_load((root / "study.yaml").read_text(encoding="utf-8"))
    value["meta"]["knudsen_suite"] = "full_decades"
    value.update(cases_dir="Cases_full_suite", manifest_name="manifest-full-suite.yaml")
    value["kn_cases"].extend({"name": name, "Kn": kn, "enabled": True}
                             for name, kn in [("Kn100", 100.), ("Kn0p1", .1), ("Kn0p01", .01)])
    (root / "study.yaml").write_text(yaml.safe_dump(value), encoding="utf-8")
    return value


def full_case(root, name, kn):
    directory = root / "Cases_full_suite" / name
    directory.mkdir(parents=True)
    cfg = yaml.safe_load((root / "baseCase" / "case.yaml").read_text(encoding="utf-8"))
    cfg["exit"]["knudsen"] = kn
    cfg["meta"]["knudsen_suite"] = "full_decades"
    (directory / "case.yaml").write_text(yaml.safe_dump(cfg), encoding="utf-8")
    return {"name": name, "path": directory.relative_to(root).as_posix(), "Kn": kn}


def test_default_inventory_skips_missing_results_and_does_not_opt_in(family):
    archive(family)
    result = inventory.inventory(family)
    assert len(result["cases"]) == 4
    assert result["missing_results"] == list(inventory.TARGET_DECADES)
    assert not any(row["fit_ready"] or row["suitable"] for row in result["cases"])
    assert all(inventory.load_case(row) is None for row in result["cases"])


def test_explicit_registration_prepares_all_decades_but_legacy_guards_still_refuse(family):
    path = register(family)
    result = inventory.inventory(family, path)
    assert {s["Kn"] for s in result["prepared_specs"]} == set(inventory.TARGET_DECADES)
    assert len(result["cases"]) == 7
    assert policy.check_manifest({"cases": [{"name": "Kn100", "Kn": 100.0}]})
    cfg = load_case_config(family / "baseCase")
    cfg = dataclasses.replace(cfg, exit=dataclasses.replace(cfg.exit, knudsen=100.0))
    with pytest.raises(policy.ExtensionError):
        policy.check_template(cfg)


def test_registration_cannot_refer_to_another_study(family):
    path = register(family, [{"name": "other", "path": "../cai2012/Cases/Kn100", "Kn": 100.0}])
    with pytest.raises(inventory.InventoryError, match="inside this study"):
        inventory.inventory(family, path)


def test_paper_decade_needs_separate_opt_in(family):
    path = register(family, allow_paper_decades=False)
    with pytest.raises(inventory.InventoryError, match="allow_paper_decades"):
        inventory.inventory(family, path)


def test_registration_rejects_production_purpose_and_duplicate_paths(family):
    path = register(family, purpose="production")
    with pytest.raises(inventory.InventoryError, match="never authorizes"):
        inventory.read_registration(family, path)
    entry = archive(family)
    path = register(family, [entry, dict(entry, name="duplicate")])
    with pytest.raises(inventory.InventoryError, match="Duplicate"):
        inventory.read_registration(family, path)


def test_partial_contour_export_does_not_become_whole_field_or_equal_volumes(family):
    entry = archive(family)
    results = family / entry["path"] / "results"
    results.mkdir()
    (results / "density_plane.csv").write_text("x_over_D,z_over_D,n_over_n0\n1,1,.02\n", encoding="utf-8")
    document = inventory.inventory(family, register(family, [entry]))
    row = next(r for r in document["cases"] if r["name"] == entry["name"])
    assert row["result_status"] == "partial-exports-only"
    assert row["partial_exports"] == ["density_plane.csv"]
    assert not row["fit_ready"]
    assert inventory.load_case(row) is None


def test_npz_adapter_preserves_actual_volumes_temperature_and_no_sample_mask(family):
    entry = archive(family)
    values = family / entry["path"] / "field.npz"
    np.savez(values, points=[[.1, .2], [2., .1]], density=[1e-303, .02],
             volumes=[.0001, .003], temperature=[0., .4], parcel_count=[0., 12.],
             velocity=[[0., 0., 0.], [300., 0., 0.]])
    entry["field_npz"] = f"{entry['path']}/field.npz"
    document = inventory.inventory(family, register(family, [entry]))
    row = next(r for r in document["cases"] if r["name"] == entry["name"])
    assert row["fit_ready"]
    assert not row["suitable"]  # Absent quality evidence never passes.
    sampled = inventory.load_case(row)
    np.testing.assert_array_equal(sampled.volumes, [.0001, .003])
    np.testing.assert_array_equal(sampled.temperature, [0., .4])
    np.testing.assert_array_equal(sampled.occupied, [False, True])
    assert sampled.metadata["velocity_sqrt_beta0"].shape == (2, 3)
    with pytest.raises(inventory.InventoryError, match="exceed limit"):
        inventory.load_case(row, max_bytes=1)


def test_npz_without_volumes_is_not_fit_ready(family):
    entry = archive(family)
    values = family / entry["path"] / "field.npz"
    np.savez(values, points=[[1., 0.]], density=[.1])
    entry["field_npz"] = f"{entry['path']}/field.npz"
    document = inventory.inventory(family, register(family, [entry]))
    row = next(r for r in document["cases"] if r["name"] == entry["name"])
    assert not row["fit_ready"]
    assert "needs" in " ".join(row["warnings"])


def test_differing_source_parameters_form_separate_strata(family):
    first = archive(family, "first", 1., 2.)
    second = archive(family, "second", 10., 3.)
    document = inventory.inventory(family, register(family, [first, second]))
    selected = {r["name"]: r for r in document["cases"]}
    assert selected["first"]["stratum"] != selected["second"]["stratum"]
    assert len(document["compatibility_strata"]) == 2


def test_quality_requires_every_stated_criterion():
    assessment = {"production_ready": True}
    diagnostics = {"resolution": {"target_cell_over_mean_free_path": 1/3,
                                   "bounding_box_max": {"max": .1}}}
    occupancy = {"regions": [{"name": "plume", "fraction_below_floor": 0.}],
                 "weight_check": {"consistent": True}}
    stationarity = {"assessment": {"verdict": "stable-window-demonstrated"}}
    assert inventory.reference_suitability(assessment, diagnostics, occupancy, stationarity, True)["suitable"]
    stationarity["assessment"]["verdict"] = "inconclusive"
    assert not inventory.reference_suitability(assessment, diagnostics, occupancy, stationarity, True)["suitable"]
    assert not inventory.reference_suitability({}, {}, {}, {}, True)["suitable"]


def test_preparation_only_writes_small_configurations_and_preserves_existing_inputs(family):
    destination = family / "prepared"
    paths = inventory.prepare_missing_specs(family, destination, register(family))
    assert len(paths) == 14
    assert {p.name for p in destination.rglob("*") if p.is_file()} == {"case.yaml", "specification.yaml"}
    assert not list(destination.rglob("polyMesh"))
    assert not list(destination.rglob("Allrun"))
    first = next(p for p in paths if p.name == "case.yaml")
    value = yaml.safe_load(first.read_text(encoding="utf-8"))
    assert value["meta"]["production_execution_authorized"] is False
    first.write_text("changed: true\n", encoding="utf-8")
    with pytest.raises(inventory.InventoryError, match="Refusing to overwrite"):
        inventory.prepare_missing_specs(family, destination, family / "registration.yaml")


def test_external_configurations_are_provenance_only(family):
    other = family.parent / "cai2012"
    other.mkdir()
    (other / "study.yaml").write_text("kn_cases:\n- {name: Kn100, Kn: 100}\n", encoding="utf-8")
    result = inventory.inventory(family)
    assert result["external_equivalents"][0]["Kn"] == 100
    assert result["external_equivalents"][0]["included"] is False
    assert not any(r["Kn"] == 100 for r in result["cases"])


def test_inaccessible_configuration_is_reported_and_other_cases_remain_visible(family, monkeypatch):
    entry = archive(family)
    original = inventory.load_case_config

    def inaccessible(directory):
        if Path(directory).name == "archive":
            raise PermissionError("configuration denied")
        return original(directory)

    monkeypatch.setattr(inventory, "load_case_config", inaccessible)
    document = inventory.inventory(family, register(family, [entry]))
    row = next(r for r in document["cases"] if r["name"] == entry["name"])
    assert row["result_status"] == "missing-configuration"
    assert not row["fit_ready"]
    assert "denied" in " ".join(row["warnings"])
    assert len(document["cases"]) == 8


def test_interface_manifest_is_visible_but_not_a_compatible_source(family, monkeypatch):
    entry = archive(family, "interface")
    (family / "interface-manifest.yaml").write_text(yaml.safe_dump({"cases": [entry]}), encoding="utf-8")
    monkeypatch.setattr(inventory.interface, "load_interface", lambda path: object() if Path(path).name == "interface" else None)
    document = inventory.inventory(family)
    row = next(r for r in document["cases"] if r["name"] == "interface")
    assert row["source_type"] == "downstream-interface"
    assert not row["compatible"]
    assert not row["fit_ready"]


def test_full_suite_is_registered_from_study_without_analysis_opt_in(family):
    full_study(family)
    result = inventory.inventory(family)
    assert len(result["cases"]) == 7
    assert {row["Kn"] for row in result["cases"]} == set(inventory.TARGET_DECADES)
    assert all("Cases_full_suite" in row["path"] for row in result["cases"])
    assert result["knudsen_suite"] == "full_decades"
    assert Path(result["active_manifest"]).name == "manifest-full-suite.yaml"
    assert result["missing_results"] == list(inventory.TARGET_DECADES)
    assert all(inventory.load_case(row) is None for row in result["cases"])


def test_full_suite_manifest_auto_wires_paper_kn_field_and_ignores_old_manifest(family):
    full_study(family)
    entry = full_case(family, "Kn100_np1x", 100.)
    field = family / entry["path"] / "full-field.npz"
    np.savez(field, points=[[.2, .1], [2., .4]], density=[.4, .03],
             volumes=[1e-4, 5e-4], occupied=[True, True])
    entry["field_npz"] = field.relative_to(family).as_posix()
    active = {"meta": {"knudsen_suite": "full_decades"}, "cases": [entry]}
    (family / "manifest-full-suite.yaml").write_text(yaml.safe_dump(active), encoding="utf-8")
    # A stale or incompatible old manifest must not override the active source.
    (family / "manifest.yaml").write_text("cases: [{name: obsolete, Kn: 100}]\n", encoding="utf-8")
    result = inventory.inventory(family)
    selected = next(row for row in result["cases"] if row["Kn"] == 100.)
    assert selected["fit_ready"]
    assert selected["registration"] == "manifest-full-suite.yaml + study.yaml"
    assert 100. not in result["missing_results"]
    assert len(result["cases"]) == 7
    loaded = inventory.load_case(selected)
    assert loaded.kn == 100.
    np.testing.assert_array_equal(loaded.volumes, [1e-4, 5e-4])


def test_active_manifest_profile_must_match_full_suite_study(family):
    full_study(family)
    (family / "manifest-full-suite.yaml").write_text(yaml.safe_dump({
        "cases": [{"name": "old", "path": "Cases/old", "Kn": 1.}]}), encoding="utf-8")
    with pytest.raises(inventory.InventoryError, match="profile disagrees"):
        inventory.inventory(family)


def test_known_low_kn_surface_is_exploratory_even_if_other_evidence_passes(family):
    full_study(family)
    entry = full_case(family, "Kn0p001_np1x", .001)
    field = family / entry["path"] / "full-field.npz"
    np.savez(field, points=[[1., .1], [2., .2]], density=[.1, .02], volumes=[.01, .02])
    entry["field_npz"] = field.relative_to(family).as_posix()
    results = family / entry["path"] / "results"
    results.mkdir()
    evidence = {
        "diagnostics_summary.yaml": {"resolution": {"target_cell_over_mean_free_path": 1.,
                                                      "bounding_box_max": {"max": .1}}},
        "occupancy.yaml": {"regions": [{"name": "plume", "fraction_below_floor": 0.}],
                           "weight_check": {"consistent": True}},
        "stationarity.yaml": {"assessment": {"verdict": "stable-window-demonstrated"}}}
    for name, value in evidence.items():
        (results / name).write_text(yaml.safe_dump(value), encoding="utf-8")
    (family / "manifest-full-suite.yaml").write_text(yaml.safe_dump({
        "meta": {"knudsen_suite": "full_decades"}, "cases": [entry]}), encoding="utf-8")
    result = inventory.inventory(family)
    selected = next(row for row in result["cases"] if row["Kn"] == .001)
    assert selected["compatible"] and selected["fit_ready"]
    assert not selected["suitable"]
    assert "known_low_kn_inflow_surface" in " ".join(selected["known_source_limitations"])
    assert inventory.load_case(selected) is not None
