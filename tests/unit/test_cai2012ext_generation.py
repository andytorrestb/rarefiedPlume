"""Generate the selected Kn suite idempotently, never over existing results.

The generator script runs once per module, into a scratch copy of the shipped
family, exactly as it runs in ``cases/cai2012-kn-extension``. Everything after
that inspects or perturbs the generated tree -- conflicts, results, strays,
tampering -- through the same functions the script uses.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest
import yaml

from plumetools.cai2012.study import load_study
from plumetools.cai2012 import dictionaries, mesh
from plumetools.cai2012.config import load_case_config
from plumetools.cai2012ext import generation, policy, summary, vizspec
from plumetools.cai2012ext.verification import verify_extension
from plumetools.viz.spec import load_spec

REPO = Path(__file__).resolve().parents[2]
FAMILY = REPO / "cases" / "cai2012-kn-extension"


def make_family(root: Path) -> Path:
    """A copy of the shipped family's inputs -- template, study, viz spec."""
    root.mkdir(parents=True, exist_ok=True)
    shutil.copytree(FAMILY / "baseCase", root / "baseCase",
                    ignore=shutil.ignore_patterns("system", "constant",
                                                  "results", "case-summary.json"))
    for name in ("study.yaml", "viz.yaml"):
        shutil.copy2(FAMILY / name, root / name)
    return root


def run_generator(family: Path, *args) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(FAMILY / "generate_cases.py"),
         "--study", str(family / "study.yaml"), *args],
        capture_output=True, text=True, timeout=900)


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    family = make_family(tmp_path_factory.mktemp("ext") / "family")
    result = run_generator(family)
    assert result.returncode == 0, result.stdout + result.stderr
    return family


def manifest_of(family: Path) -> dict:
    study = load_study(family / "study.yaml")
    return yaml.safe_load((family / study.manifest_name).read_text(encoding="utf-8"))


def cases_root(family: Path) -> Path:
    return family / load_study(family / "study.yaml").cases_dir


def clone(generated: Path, tmp_path: Path) -> Path:
    target = tmp_path / "family"
    shutil.copytree(generated, target)
    return target


# --------------------------------------------------------------------------- #
# what is generated
# --------------------------------------------------------------------------- #

def test_exactly_seven_baseline_cases_with_explicit_full_suite(generated):
    manifest = manifest_of(generated)
    assert manifest["n_generated"] == 7
    assert manifest["meta"]["knudsen_suite"] == "full_decades"
    pairs = sorted((c["Kn"], c["numerical_particle_multiplier"])
                   for c in manifest["cases"])
    assert pairs == [(1e-4, 1.0), (1e-3, 1.0), (1e-2, 1.0), (.1, 1.0),
                     (1.0, 1.0), (10.0, 1.0), (100.0, 1.0)]
    on_disk = sorted(p.name for p in cases_root(generated).iterdir())
    assert on_disk == ["Kn0p0001_np1x", "Kn0p001_np1x", "Kn0p01_np1x", "Kn0p1_np1x",
                       "Kn100_np1x", "Kn10_np1x", "Kn1_np1x"]
    for name in on_disk:
        # The diagnostics' settings travel with every case, like the scripts.
        assert (cases_root(generated) / name / "diagnostics.yaml").is_file()
        assert (cases_root(generated) / name / "Alldiagnose.py").is_file()
        cfg = load_case_config(cases_root(generated) / name)
        assert cfg.meta["knudsen_suite"] == "full_decades"
        # ... and the PBS job scripts, submitted from the case itself.
        for job in ("runMesh.sh", "runCase.sh", "runPost.sh", "runViz.sh"):
            assert (cases_root(generated) / name / job).read_bytes() == (
                generated / "baseCase" / job).read_bytes()


def test_d_s0_t0_and_the_gas_are_held_fixed(generated):
    held = manifest_of(generated)["held_fixed"]
    assert held["exit.D_m"] == 0.2
    assert held["exit.S0"] == 2.0
    assert held["exit.T0_K"] == 300.0
    assert held["gas.omega"] == 0.74
    assert held["exit.mean_free_path_convention"] == "vhs"
    assert not any(isinstance(v, str) and v.startswith("NOT SHARED")
                   for v in held.values())
    assert manifest_of(generated)["held_fixed_problems"] == []


def test_within_a_kn_only_the_weight_moves(generated):
    per_kn = manifest_of(generated)["per_knudsen"]
    assert set(per_kn) == {"Kn100", "Kn10", "Kn1", "Kn0p1", "Kn0p01", "Kn0p001", "Kn0p0001"}
    for group in per_kn.values():
        for key in ("mesh_id", "deltaT_s", "average_start_s", "end_time_s",
                    "write_interval_steps"):
            assert not str(group[key]).startswith("NOT SHARED"), key


def test_every_row_records_the_promised_provenance(generated):
    for entry in manifest_of(generated)["cases"]:
        exit_ = entry["exit"]
        for key in ("Kn_target", "Kn_recomputed", "D_m", "S0", "T0_K",
                    "U0_m_per_s", "n0_per_m3", "rho0_kg_per_m3", "p0_Pa",
                    "lambda0_m"):
            assert key in exit_, key
        assert abs(exit_["Kn_relative_deviation"]) < 1e-12
        assert entry["weight"]["n_equivalent_particles"] == pytest.approx(
            entry["weight"]["baseline_n_equivalent_particles"]
            / entry["numerical_particle_multiplier"], rel=1e-15)
        assert entry["startup"]["achieved_reference_collision_times"] == \
            pytest.approx(1.0e4)
        assert entry["sampling"]["requested_sampling_domain_transits"] == 4.5
        assert entry["sources"]["case"]["sha256"]
        assert entry["feasibility"]["status"] in (
            "production-candidate", "unsuitable-for-production")
        assert entry["production_status"] == entry["feasibility"]["status"]


def test_the_startup_is_the_same_instant_in_all_seven(generated):
    starts = {c["average_start_s"] for c in manifest_of(generated)["cases"]}
    assert len(starts) == 1


def test_the_verdicts_split_on_kn(generated):
    status = {c["name"]: c["production_status"]
              for c in manifest_of(generated)["cases"]}
    for name, verdict in status.items():
        dense = name.startswith(("Kn0p01_", "Kn0p001_", "Kn0p0001_"))
        assert verdict == ("unsuitable-for-production" if dense
                           else "production-candidate"), name


def test_the_inputs_verify_independently(generated):
    report = verify_extension(generated, manifest_of(generated),
                              with_dictionaries=False)
    assert report.passed, report.failures


def test_all_seven_generated_dictionary_sets_match_the_manifest_without_meshing(generated, tmp_path):
    """Exercise Allmesh's Python writers without any OpenFOAM executable."""
    family = clone(generated, tmp_path)
    manifest = manifest_of(family)
    for entry in manifest["cases"]:
        case_dir = family / entry["path"]
        cfg, geom, state, plan, run, _, _, checks = generation.derive(case_dir)
        checks.raise_if_failed()
        paths = mesh.write_mesh_setup(case_dir, cfg, geom, plan)
        paths += dictionaries.write_case_dictionaries(case_dir, cfg, state, run)
        assert len(paths) == 10
        assert all(path.is_file() for path in paths)
        assert not (case_dir / "constant" / "polyMesh").exists()
        assert not generation.has_solver_output(case_dir)
        assert "nozzle" in (case_dir / "system" / "createPatchDict").read_text(encoding="utf-8")
    report = verify_extension(family, manifest, with_dictionaries=True)
    assert report.passed, report.failures


def test_very_low_kn_source_limitation_is_retained_without_changing_the_source(generated):
    manifest = manifest_of(generated)
    for entry in manifest["cases"]:
        case_dir = generated / entry["path"]
        cfg = load_case_config(case_dir)
        expected = policy.source_limitations(cfg)
        assert entry["known_source_limitations"] == expected
        assert cfg.meta["known_source_limitations"] == expected
        assert not (case_dir / "interface.yaml").exists()
        if entry["Kn"] <= .001:
            assert expected and "inflow surface" in expected[0]
            assert entry["feasibility"]["production_ready"] is False
        else:
            assert expected == []


def test_declared_per_kn_occupancy_targets_are_checked_per_case(generated):
    manifest = manifest_of(generated)
    targets = {entry["Kn"]: entry["weight"]["baseline_exit_particles_per_cell"]
               for entry in manifest["cases"]}
    assert targets[10.] == 2000
    assert targets[1.] == 2000
    assert targets[.001] == 40
    report = verify_extension(generated, manifest, with_dictionaries=False)
    assert report.passed, report.failures
    altered = yaml.safe_load(yaml.safe_dump(manifest))
    altered["cases"][0]["weight"]["baseline_exit_particles_per_cell"] *= 2
    report = verify_extension(generated, altered, with_dictionaries=False)
    assert any("configured target" in line for line in report.failures)


def test_an_unmarked_historical_study_still_generates_four_cases(tmp_path):
    family = make_family(tmp_path / "legacy-family")
    document = yaml.safe_load((family / "study.yaml").read_text(encoding="utf-8"))
    document["meta"].pop("knudsen_suite", None)
    document["kn_cases"] = [row for row in document["kn_cases"]
                            if policy.is_extension_knudsen(row["Kn"])]
    document.update(cases_dir="LegacyCases", manifest_name="manifest-legacy.yaml")
    (family / "study.yaml").write_text(yaml.safe_dump(document), encoding="utf-8")
    result = run_generator(family)
    assert result.returncode == 0, result.stdout + result.stderr
    manifest = manifest_of(family)
    assert manifest["n_generated"] == 4
    assert not any(policy.is_excluded(row["Kn"]) for row in manifest["cases"])
    report = verify_extension(family, manifest, with_dictionaries=False)
    assert report.passed, report.failures


def test_full_suite_generation_preserves_prior_case_results_and_manifest(tmp_path):
    family = make_family(tmp_path / "family")
    original_case = family / "Cases" / "Kn0p001_np1x" / "0.06"
    original_case.mkdir(parents=True)
    field = original_case / "rhoNMean"
    field.write_bytes(b"existing DSMC result; must not be copied or overwritten\n")
    prior_manifest = family / "manifest.yaml"
    prior_manifest.write_bytes(b"historical manifest\n")
    result = run_generator(family)
    assert result.returncode == 0, result.stdout + result.stderr
    assert field.read_bytes() == b"existing DSMC result; must not be copied or overwritten\n"
    assert prior_manifest.read_bytes() == b"historical manifest\n"
    assert manifest_of(family)["n_generated"] == 7
    assert not any((family / row["path"] / "0.06").exists()
                   for row in manifest_of(family)["cases"])


def test_losing_the_explicit_case_profile_fails_full_suite_validation(generated, tmp_path):
    family = clone(generated, tmp_path)
    entry = next(row for row in manifest_of(family)["cases"] if row["Kn"] == 100)
    path = family / entry["path"] / "case.yaml"
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    config["meta"].pop("knudsen_suite")
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    report = verify_extension(family, manifest_of(family), with_dictionaries=False)
    assert any("different Knudsen suites" in line for line in report.failures)


# --------------------------------------------------------------------------- #
# idempotency and conflicts
# --------------------------------------------------------------------------- #

def test_regenerating_keeps_every_case_untouched(generated, tmp_path):
    family = clone(generated, tmp_path)
    case_yaml = cases_root(family) / "Kn1_np1x" / "case.yaml"
    (cases_root(family) / "Kn1_np1x" / "constant").mkdir()
    marker = cases_root(family) / "Kn1_np1x" / "constant" / "mesh-marker"
    marker.write_text("a mesh would be here", encoding="utf-8")
    before = case_yaml.stat().st_mtime_ns
    result = run_generator(family)
    assert result.returncode == 0, result.stderr
    assert result.stdout.count("keep (identical)") == 7
    assert case_yaml.stat().st_mtime_ns == before
    assert marker.is_file()


def _states(family: Path, tmp_path: Path, names, **flags):
    study = load_study(family / "study.yaml")
    cases = [c for c in study.expand() if c.name in names]
    scratch = Path(tempfile.mkdtemp(dir=tmp_path))
    return study, generation.decide(
        generation.classify(family, study, scratch, cases), **flags)


def _edit_input(family: Path, name: str):
    path = cases_root(family) / name / "case.yaml"
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    data["mesh"]["max_cells"] = 123456789
    path.write_text(yaml.safe_dump(data), encoding="utf-8")


def _add_results(family: Path, name: str):
    (cases_root(family) / name / "0.06").mkdir()


def test_an_input_conflict_is_refused_not_overwritten(generated, tmp_path):
    family = clone(generated, tmp_path)
    _edit_input(family, "Kn10_np1x")
    _, states = _states(family, tmp_path, {"Kn10_np1x", "Kn1_np1x"})
    by_name = {s.case.name: s for s in states}
    assert by_name["Kn10_np1x"].status == generation.STATUS_CONFLICT
    assert by_name["Kn10_np1x"].action == "refuse"
    assert any("mesh.max_cells" in line for line in by_name["Kn10_np1x"].inputs)
    assert by_name["Kn1_np1x"].action == "keep"
    with pytest.raises(policy.ExtensionError):
        generation.apply(states)


def test_regenerate_replaces_a_conflict_only_without_results(generated, tmp_path):
    family = clone(generated, tmp_path)
    _edit_input(family, "Kn10_np1x")
    _edit_input(family, "Kn1_np1x")
    _add_results(family, "Kn1_np1x")
    _, states = _states(family, tmp_path, {"Kn10_np1x", "Kn1_np1x"},
                        regenerate=True)
    by_name = {s.case.name: s for s in states}
    assert by_name["Kn10_np1x"].action == "replace"
    assert by_name["Kn1_np1x"].has_results
    assert by_name["Kn1_np1x"].action == "refuse"


def test_results_under_processor_directories_count_as_results(generated, tmp_path):
    family = clone(generated, tmp_path)
    (cases_root(family) / "Kn1_np1x" / "processor0" / "0.01").mkdir(parents=True)
    assert generation.has_solver_output(cases_root(family) / "Kn1_np1x")
    assert not generation.has_solver_output(cases_root(family) / "Kn10_np1x")


def test_only_overwrite_replaces_a_case_with_results(generated, tmp_path):
    family = clone(generated, tmp_path)
    _edit_input(family, "Kn1_np1x")
    _add_results(family, "Kn1_np1x")
    _, states = _states(family, tmp_path, {"Kn1_np1x"}, overwrite=True)
    assert states[0].action == "replace"
    generation.apply(states)
    restored = yaml.safe_load(
        (cases_root(family) / "Kn1_np1x" / "case.yaml").read_text(encoding="utf-8"))
    expected = load_study(family / "study.yaml")
    target_case = next(c for c in expected.expand() if c.name == "Kn1_np1x")
    assert restored["mesh"]["max_cells"] == target_case.overrides["mesh"]["max_cells"]
    assert not (cases_root(family) / "Kn1_np1x" / "0.06").exists()


def test_a_changed_runner_script_is_stale_not_a_conflict(generated, tmp_path):
    family = clone(generated, tmp_path)
    (cases_root(family) / "Kn0p001_np1x" / "Allpost").write_text(
        "#!/bin/bash\n", encoding="utf-8")
    (cases_root(family) / "Kn0p001_np1x" / "diagnostics.yaml").unlink()
    _add_results(family, "Kn0p001_np1x")
    _, states = _states(family, tmp_path, {"Kn0p001_np1x"})
    assert states[0].status == generation.STATUS_STALE
    assert states[0].action == "keep-stale"
    assert any("diagnostics.yaml: missing" in line for line in states[0].auxiliary)
    _, states = _states(family, tmp_path, {"Kn0p001_np1x"}, refresh=True)
    assert states[0].action == "refresh"
    generation.apply(states)
    case = cases_root(family) / "Kn0p001_np1x"
    assert (case / "Allpost").read_bytes() == (
        family / "baseCase" / "Allpost").read_bytes()
    assert (case / "diagnostics.yaml").is_file()
    assert (case / "0.06").is_dir()


# --------------------------------------------------------------------------- #
# strays and tampering
# --------------------------------------------------------------------------- #

def test_a_stray_case_outside_the_full_suite_fails_validation(generated, tmp_path):
    family = clone(generated, tmp_path)
    stray = cases_root(family) / "Kn0p02_np1x"
    shutil.copytree(cases_root(family) / "Kn1_np1x", stray)
    data = yaml.safe_load((stray / "case.yaml").read_text(encoding="utf-8"))
    data["exit"]["knudsen"] = 0.02
    (stray / "case.yaml").write_text(yaml.safe_dump(data), encoding="utf-8")

    study = load_study(family / "study.yaml")
    strays = generation.stray_directories(family, study, study.expand())
    assert strays == [{"path": f"{study.cases_dir}/Kn0p02_np1x", "Kn": 0.02, "excluded": True}]
    report = verify_extension(family, manifest_of(family), with_dictionaries=False)
    assert any("excluded Kn = 0.02" in f for f in report.failures)
    # ... and it is never scheduled, because nothing reads the directory.
    runnable, _ = policy.schedulable(manifest_of(family), allow_unsuitable=True)
    assert f"{study.cases_dir}/Kn0p02_np1x" not in [e["path"] for e in runnable]


def test_a_hand_edited_case_fails_validation(generated, tmp_path):
    family = clone(generated, tmp_path)
    _edit_input(family, "Kn0p0001_np1x")
    report = verify_extension(family, manifest_of(family), with_dictionaries=False)
    assert any("Kn0p0001_np1x" in f and "changed since" in f
               for f in report.failures)


# --------------------------------------------------------------------------- #
# downstream of the manifest
# --------------------------------------------------------------------------- #

def test_the_viz_specs_carry_the_new_averaging_start_and_n0(generated, tmp_path):
    manifest = manifest_of(generated)
    written = vizspec.write_group_specs(manifest, generated / "viz.yaml",
                                        tmp_path / "specs")
    assert set(written) == {"Kn100", "Kn10", "Kn1", "Kn0p1", "Kn0p01", "Kn0p001", "Kn0p0001"}
    start = manifest["cases"][0]["average_start_s"]
    assert start > 0.05      # Cai's 10^4 t0, not the health family's 5.658e-3
    for name, (path, cases) in written.items():
        spec = load_spec(path)
        assert spec.sampling.time_min == pytest.approx(start)
        assert spec.sampling.time == "all"
        n0 = next(c["n0_per_m3"] for c in manifest["cases"]
                  if c["cai_case"] == name)
        derived = {d.name: d for d in spec.derived}
        assert set(derived) >= {"U", "Ttra", "n_over_n0"}
        assert f"{n0:.9e}" in derived["n_over_n0"].expression
        fields = {f.name: f for f in spec.fields}
        assert tuple(fields["n_over_n0"].range) == vizspec.N_OVER_N0_RANGE
        assert len(cases) == 1


def test_the_diagnostics_spec_draws_the_latest_time_with_pinned_criteria():
    spec = load_spec(FAMILY / "viz-diagnostics.yaml")
    assert spec.sampling.time == "latest"
    fields = {f.name: f for f in spec.fields}
    assert {"mfpResolved", "KnCellMax", "mfp", "KnD", "cellOverMfpMax",
            "parcelsPerStep", "diagValid", "diagStatus"} <= set(fields)
    assert tuple(fields["KnCellMax"].range) == (0.003, 3000.0)
    assert tuple(fields["mfpResolved"].range) == (-1.0, 1.0)


def test_the_study_table_lists_cases_that_have_not_run(generated, tmp_path):
    manifest = manifest_of(generated)
    rows = [summary.study_row(generated, entry) for entry in manifest["cases"]]
    assert len(rows) == 7
    assert {row["state"] for row in rows} == {"not run"}
    assert all(row["reference"] == summary.REFERENCE for row in rows)
    path = summary.write_table(tmp_path / "table.csv", rows)
    header = path.read_text(encoding="utf-8").splitlines()[0].split(",")
    assert tuple(header) == summary.COLUMNS
    assert "cai_reported" not in ",".join(header)
