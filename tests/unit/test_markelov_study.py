"""study.yaml parsing, the case matrix, and the generation manifest."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest
import yaml

from plumetools.markelov1999.constants import PSI_TO_PA
from plumetools.markelov1999.study import (
    CASE_TEMPLATE_ENTRIES,
    PressureCase,
    StudyError,
    apply_case_parameters,
    case_is_complete,
    clone_base_case,
    filter_by_pressure,
    generation_backend,
    load_study,
    write_manifest,
)

STUDY_YAML = """\
base_case: baseCase
cases_dir: Cases
gap_in: 6.0
gap_dir: gap06in
pressure_cases:
  - {name: p005psi, pressure_psi: 5,   enabled: true}
  - {name: p025psi, pressure_psi: 25,  enabled: true}
  - {name: p100psi, pressure_psi: 100, enabled: true}
  - {name: p475psi, pressure_psi: 475, enabled: true}
"""


def write_study(tmp_path: Path, text: str = STUDY_YAML) -> Path:
    path = tmp_path / "study.yaml"
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture
def study(tmp_path):
    return load_study(write_study(tmp_path))


# --------------------------------------------------------------------------- #
# the four-case matrix
# --------------------------------------------------------------------------- #

def test_the_four_paper_pressures_are_generated(study):
    cases = study.enabled_cases()
    assert len(cases) == 4
    assert [c.name for c in cases] == ["p005psi", "p025psi", "p100psi", "p475psi"]
    assert [c.pressure_psi for c in cases] == [5.0, 25.0, 100.0, 475.0]


def test_pressures_convert_to_pascals_exactly(study):
    for case in study.enabled_cases():
        assert case.pressure_pa == pytest.approx(case.pressure_psi * PSI_TO_PA, rel=1e-15)
    assert study.enabled_cases()[0].pressure_pa == pytest.approx(34473.786466, rel=1e-9)
    assert study.enabled_cases()[-1].pressure_pa == pytest.approx(3275009.714255, rel=1e-9)


def test_case_ordering_is_deterministic_regardless_of_file_order(tmp_path):
    """Reordering study.yaml must not change the generated tree, or two people
    editing the same matrix produce different-looking output from it."""
    shuffled = textwrap.dedent("""\
        pressure_cases:
          - {name: p475psi, pressure_psi: 475}
          - {name: p005psi, pressure_psi: 5}
          - {name: p100psi, pressure_psi: 100}
          - {name: p025psi, pressure_psi: 25}
        """)
    study = load_study(write_study(tmp_path, shuffled))
    assert [c.name for c in study.enabled_cases()] == [
        "p005psi", "p025psi", "p100psi", "p475psi"]


def test_case_paths_form_the_expected_hierarchy(study):
    paths = [str(study.case_path(c)).replace("\\", "/") for c in study.enabled_cases()]
    assert paths == [
        "Cases/gap06in/p005psi",
        "Cases/gap06in/p025psi",
        "Cases/gap06in/p100psi",
        "Cases/gap06in/p475psi",
    ]


# --------------------------------------------------------------------------- #
# enabling and disabling
# --------------------------------------------------------------------------- #

def test_a_disabled_case_is_excluded_but_still_recorded(tmp_path):
    """Requirement 1: individual pressure cases must be easy to disable. The
    entry stays in study.yaml, so the intended matrix is not lost."""
    text = STUDY_YAML.replace(
        "{name: p025psi, pressure_psi: 25,  enabled: true}",
        "{name: p025psi, pressure_psi: 25,  enabled: false}")
    study = load_study(write_study(tmp_path, text))

    assert [c.name for c in study.enabled_cases()] == [
        "p005psi", "p100psi", "p475psi"]
    assert len(study.all_cases()) == 4
    assert [c.name for c in study.all_cases() if not c.enabled] == ["p025psi"]


def test_enabled_defaults_to_true(tmp_path):
    text = "pressure_cases:\n  - {name: p005psi, pressure_psi: 5}\n"
    assert load_study(write_study(tmp_path, text)).enabled_cases()[0].enabled


def test_every_case_disabled_gives_an_empty_matrix(tmp_path):
    text = STUDY_YAML.replace("enabled: true", "enabled: false")
    assert load_study(write_study(tmp_path, text)).enabled_cases() == []


# --------------------------------------------------------------------------- #
# the CLI pressure filter
# --------------------------------------------------------------------------- #

def test_pressure_filter_selects_a_subset(study):
    kept = filter_by_pressure(study.enabled_cases(), [5, 25])
    assert [c.name for c in kept] == ["p005psi", "p025psi"]


def test_no_filter_keeps_everything(study):
    assert len(filter_by_pressure(study.enabled_cases(), None)) == 4
    assert len(filter_by_pressure(study.enabled_cases(), [])) == 4


def test_a_pressure_that_matches_nothing_is_an_error(study):
    """A typo on the command line must stop the run, not silently generate
    nothing."""
    with pytest.raises(StudyError, match=r"no case at \[7.0\] psi"):
        filter_by_pressure(study.enabled_cases(), [7])


# --------------------------------------------------------------------------- #
# validation
# --------------------------------------------------------------------------- #

def test_a_missing_study_file_is_an_error(tmp_path):
    with pytest.raises(StudyError, match="no study.yaml"):
        load_study(tmp_path / "nope.yaml")


def test_a_directory_is_accepted_and_resolves_to_study_yaml(tmp_path):
    write_study(tmp_path)
    assert len(load_study(tmp_path).enabled_cases()) == 4


def test_an_unknown_top_level_key_is_an_error(tmp_path):
    with pytest.raises(StudyError, match="unknown top-level key"):
        load_study(write_study(tmp_path, STUDY_YAML + "gap_mm: 152.4\n"))


def test_an_unknown_case_key_is_an_error(tmp_path):
    text = "pressure_cases:\n  - {name: p005psi, pressure_psi: 5, enabld: true}\n"
    with pytest.raises(StudyError, match="unknown key"):
        load_study(write_study(tmp_path, text))


def test_a_missing_case_key_is_an_error(tmp_path):
    text = "pressure_cases:\n  - {name: p005psi}\n"
    with pytest.raises(StudyError, match="missing 'pressure_psi'"):
        load_study(write_study(tmp_path, text))


def test_an_empty_matrix_is_an_error(tmp_path):
    with pytest.raises(StudyError, match="non-empty list"):
        load_study(write_study(tmp_path, "pressure_cases: []\n"))


def test_a_nonpositive_pressure_is_an_error(tmp_path):
    text = "pressure_cases:\n  - {name: p000psi, pressure_psi: 0}\n"
    with pytest.raises(StudyError, match="must be positive"):
        load_study(write_study(tmp_path, text))


def test_duplicate_case_names_are_an_error(tmp_path):
    text = ("pressure_cases:\n"
            "  - {name: p005psi, pressure_psi: 5}\n"
            "  - {name: p005psi, pressure_psi: 25}\n")
    with pytest.raises(StudyError, match="duplicate case name"):
        load_study(write_study(tmp_path, text))


def test_duplicate_pressures_are_an_error(tmp_path):
    text = ("pressure_cases:\n"
            "  - {name: a, pressure_psi: 5}\n"
            "  - {name: b, pressure_psi: 5}\n")
    with pytest.raises(StudyError, match="duplicate pressure"):
        load_study(write_study(tmp_path, text))


# --------------------------------------------------------------------------- #
# cloning
# --------------------------------------------------------------------------- #

@pytest.fixture
def base_case(tmp_path):
    base = tmp_path / "baseCase"
    (base / "system").mkdir(parents=True)
    (base / "constant" / "polyMesh").mkdir(parents=True)
    (base / "case.yaml").write_text(
        "model: markelov1999_axisymmetric\n"
        "stagnation: {p0_pa: 34473.79}\n"
        "dsmc: {n_equivalent_particles: null}\n"
        "meta: {reference: AIAA 99-3455}\n", encoding="utf-8")
    (base / "Allrun").write_text("#!/bin/bash\n", encoding="utf-8")
    (base / "runInflow.py").write_text("# inflow\n", encoding="utf-8")
    (base / "system" / "fvSchemes").write_text("ddtSchemes{}\n", encoding="utf-8")
    (base / "constant" / "polyMesh" / "points").write_text("0\n", encoding="utf-8")
    return base


def test_clone_copies_the_template_files(base_case, tmp_path):
    destination = tmp_path / "out" / "p005psi"
    clone_base_case(base_case, destination)

    assert (destination / "case.yaml").is_file()
    assert (destination / "Allrun").is_file()
    assert (destination / "runInflow.py").is_file()
    assert (destination / "system" / "fvSchemes").is_file()


def test_clone_does_not_inherit_the_template_mesh(base_case, tmp_path):
    """Each case meshes itself from its own case.yaml. Inheriting a mesh would
    silently give every case the template's geometry."""
    destination = tmp_path / "out" / "p005psi"
    clone_base_case(base_case, destination)
    assert not (destination / "constant" / "polyMesh").exists()


def test_clone_rejects_a_missing_or_incomplete_base(tmp_path):
    with pytest.raises(StudyError, match="does not exist"):
        clone_base_case(tmp_path / "nope", tmp_path / "out")

    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(StudyError, match="no case.yaml"):
        clone_base_case(empty, tmp_path / "out")


def test_the_template_entry_list_names_no_generated_dictionary():
    """Generated dictionaries must NOT be copied: they would be a second,
    silently divergent source of truth for the physics."""
    for generated in ("dsmcProperties", "controlDict", "blockMeshDict",
                      "snappyHexMeshDict", "dsmcInitialiseDict"):
        assert generated not in CASE_TEMPLATE_ENTRIES


# --------------------------------------------------------------------------- #
# applying case parameters
# --------------------------------------------------------------------------- #

def test_parameters_are_applied_structurally(base_case, tmp_path, study):
    """Structured YAML editing, not string substitution: the value is replaced
    whatever its original formatting."""
    destination = tmp_path / "out" / "p025psi"
    clone_base_case(base_case, destination)
    case = PressureCase(name="p025psi", pressure_psi=25.0)

    applied = apply_case_parameters(destination, case, study,
                                    n_equivalent_particles=1.5e11)

    data = yaml.safe_load((destination / "case.yaml").read_text(encoding="utf-8"))
    assert data["stagnation"]["p0_pa"] == pytest.approx(25.0 * PSI_TO_PA, rel=1e-15)
    assert data["dsmc"]["n_equivalent_particles"] == pytest.approx(1.5e11)
    assert applied["pressure_psi"] == 25.0


def test_case_metadata_records_the_pressure_and_the_gap(base_case, tmp_path, study):
    destination = tmp_path / "out" / "p100psi"
    clone_base_case(base_case, destination)
    apply_case_parameters(destination, PressureCase("p100psi", 100.0), study)

    meta = yaml.safe_load((destination / "case.yaml").read_text(encoding="utf-8"))["meta"]
    assert meta["case_name"] == "p100psi"
    assert meta["pressure_psi"] == 100.0
    assert meta["gap_in"] == 6.0
    assert meta["reference"] == "AIAA 99-3455"   # inherited, not overwritten


def test_the_generated_case_yaml_says_it_is_generated(base_case, tmp_path, study):
    destination = tmp_path / "out" / "p005psi"
    clone_base_case(base_case, destination)
    apply_case_parameters(destination, PressureCase("p005psi", 5.0), study)
    text = (destination / "case.yaml").read_text(encoding="utf-8")
    assert text.startswith("# GENERATED by")
    assert "do not edit" in text


def test_applying_to_a_file_without_the_expected_sections_is_an_error(tmp_path, study):
    """Rather than growing the sections silently, which would hide that the file
    is not this case family's template."""
    destination = tmp_path / "wrong"
    destination.mkdir()
    (destination / "case.yaml").write_text("model: source_flow\n", encoding="utf-8")
    with pytest.raises(StudyError, match="no 'stagnation' section"):
        apply_case_parameters(destination, PressureCase("x", 5.0), study)


def test_leaving_the_weight_unset_keeps_it_unset(base_case, tmp_path, study):
    destination = tmp_path / "out" / "p005psi"
    clone_base_case(base_case, destination)
    apply_case_parameters(destination, PressureCase("p005psi", 5.0), study)
    data = yaml.safe_load((destination / "case.yaml").read_text(encoding="utf-8"))
    assert data["dsmc"]["n_equivalent_particles"] is None


# --------------------------------------------------------------------------- #
# overwrite protection
# --------------------------------------------------------------------------- #

def test_a_case_with_results_is_detected(tmp_path):
    case = tmp_path / "case"
    (case / "0").mkdir(parents=True)
    assert not case_is_complete(case), "time 0 alone is not a result"

    (case / "0.004").mkdir()
    assert case_is_complete(case)


def test_a_missing_or_empty_case_is_not_complete(tmp_path):
    assert not case_is_complete(tmp_path / "nope")
    empty = tmp_path / "empty"
    empty.mkdir()
    assert not case_is_complete(empty)


def test_non_numeric_directories_are_not_mistaken_for_results(tmp_path):
    case = tmp_path / "case"
    (case / "constant").mkdir(parents=True)
    (case / "system").mkdir()
    (case / "postProcessing").mkdir()
    assert not case_is_complete(case)


# --------------------------------------------------------------------------- #
# the manifest
# --------------------------------------------------------------------------- #

def test_manifest_records_every_generated_case(tmp_path, study):
    entries = [
        {"name": c.name, "path": str(study.case_path(c)).replace("\\", "/"),
         "pressure_psi": c.pressure_psi, "pressure_pa": c.pressure_pa,
         "gap_in": 6.0, "n_equivalent_particles": 1e11}
        for c in study.enabled_cases()
    ]
    path = write_manifest(tmp_path / "manifest.yaml", study, entries, "builtin")

    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert document["n_generated"] == 4
    assert document["reference"] == "AIAA 99-3455"
    assert document["gap_in"] == 6.0
    assert document["backend"] == "builtin"
    assert [c["name"] for c in document["cases"]] == [
        "p005psi", "p025psi", "p100psi", "p475psi"]
    for case in document["cases"]:
        assert case["pressure_pa"] == pytest.approx(
            case["pressure_psi"] * PSI_TO_PA, rel=1e-12)
        assert case["n_equivalent_particles"] > 0


def test_manifest_records_disabled_cases_as_a_decision(tmp_path):
    """An absent case should be visible as a choice, not as a gap."""
    text = STUDY_YAML.replace(
        "{name: p475psi, pressure_psi: 475, enabled: true}",
        "{name: p475psi, pressure_psi: 475, enabled: false, note: too expensive}")
    study = load_study(write_study(tmp_path, text))
    path = write_manifest(tmp_path / "manifest.yaml", study, [], "builtin")

    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert document["disabled_cases"] == [
        {"name": "p475psi", "pressure_psi": 475.0, "note": "too expensive"}]


def test_manifest_says_it_is_generated(tmp_path, study):
    path = write_manifest(tmp_path / "manifest.yaml", study, [], "builtin")
    assert path.read_text(encoding="utf-8").startswith("# GENERATED by")


def test_the_backend_is_recorded(tmp_path, study):
    """CaseFoam is optional; the built-in clone produces the same tree, so this
    is provenance rather than a behaviour switch."""
    assert generation_backend() in ("casefoam", "builtin")


# --------------------------------------------------------------------------- #
# the shipped study
# --------------------------------------------------------------------------- #

def test_the_shipped_study_yaml_defines_the_four_paper_pressures():
    root = Path(__file__).resolve().parents[2] / "cases" / "markelov1999"
    study = load_study(root / "study.yaml")
    assert study.gap_in == 6.0
    assert study.gap_dir == "gap06in"
    assert [c.pressure_psi for c in study.enabled_cases()] == [5.0, 25.0, 100.0, 475.0]


def test_the_shipped_study_covers_only_the_six_inch_gap():
    """Requirement 16: the 12 in cases are out of scope."""
    root = Path(__file__).resolve().parents[2] / "cases" / "markelov1999"
    study = load_study(root / "study.yaml")
    assert study.gap_in == 6.0
