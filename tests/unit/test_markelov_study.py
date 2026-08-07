"""study.yaml parsing, the case matrix, and the generation manifest."""

from __future__ import annotations

import shutil
import textwrap
from pathlib import Path

import pytest
import yaml

from plumetools.markelov1999.constants import PSI_TO_PA
from plumetools.markelov1999.study import (
    PressureCase,
    StudyError,
    apply_case_parameters,
    case_is_complete,
    clone_cases,
    filter_by_pressure,
    load_study,
    prune_generated_dictionaries,
    require_casefoam,
    write_manifest,
)


def _casefoam_installed() -> bool:
    try:
        import casefoam  # noqa: F401
    except ImportError:
        return False
    return True

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
# cloning -- CaseFoam is a dependency, so these need it
# --------------------------------------------------------------------------- #

needs_casefoam = pytest.mark.skipif(
    not _casefoam_installed(),
    reason='casefoam is not installed (pip install -e ".[cases]")')


@pytest.fixture
def base_case(tmp_path):
    """A minimal template, shaped like cases/markelov1999/baseCase."""
    base = tmp_path / "baseCase"
    (base / "system").mkdir(parents=True)
    (base / "case.yaml").write_text(
        "model: markelov1999_axisymmetric\n"
        "stagnation: {p0_pa: 34473.79}\n"
        "dsmc: {n_equivalent_particles: null}\n"
        "meta: {reference: AIAA 99-3455}\n", encoding="utf-8")
    (base / "Allrun").write_text("#!/bin/bash\n", encoding="utf-8")
    (base / "runInflow.py").write_text("# inflow\n", encoding="utf-8")
    (base / "system" / "fvSchemes").write_text("ddtSchemes{}\n", encoding="utf-8")
    return base


@pytest.fixture
def study_root(tmp_path, base_case):
    """A study directory: the template plus the study-level files."""
    root = tmp_path / "study"
    root.mkdir()
    shutil.copytree(base_case, root / "baseCase")
    (root / "study.yaml").write_text(STUDY_YAML, encoding="utf-8")
    (root / "AllrunCases").write_text("#!/bin/bash\n", encoding="utf-8")
    (root / "README.md").write_text("# study\n", encoding="utf-8")
    (root / "generate_cases.py").write_text("# generator\n", encoding="utf-8")
    return root


def test_casefoam_is_required_not_optional():
    """There is no built-in substitute. Reimplementing a case generator the
    project already depends on would mean maintaining two, and the one that is
    not exercised is the one that drifts."""
    import plumetools.markelov1999.study as study_module
    assert not hasattr(study_module, "clone_base_case")
    assert not hasattr(study_module, "generation_backend")


def test_a_missing_casefoam_names_the_install_command(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def no_casefoam(name, *args, **kwargs):
        if name == "casefoam":
            raise ImportError("no casefoam")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_casefoam)
    with pytest.raises(StudyError, match=r'pip install -e "\.\[cases\]"'):
        require_casefoam()


@needs_casefoam
def test_casefoam_builds_the_expected_hierarchy(study_root, study):
    generated = clone_cases(study_root, study, study.enabled_cases())

    assert [p.name for p in generated] == [
        "p005psi", "p025psi", "p100psi", "p475psi"]
    for path in generated:
        assert path.is_dir()
        assert path.parent.name == "gap06in"
        assert path.parent.parent.name == "Cases"
        assert (path / "case.yaml").is_file()
        assert (path / "Allrun").is_file()
        assert (path / "system" / "fvSchemes").is_file()


@needs_casefoam
def test_the_study_directory_is_not_disturbed(study_root, study):
    """mkCases is pointed at the baseCase DIRECTORY, not the study root. Pointed
    at the root it copies study.yaml, generate_cases.py and the drivers into
    Cases/ -- which is what the first attempt at this integration got wrong."""
    clone_cases(study_root, study, study.enabled_cases())

    for name in ("study.yaml", "AllrunCases", "README.md", "generate_cases.py"):
        assert (study_root / name).is_file(), f"{name} disappeared"

    leaked = {p.name for p in (study_root / "Cases").rglob("*")}
    for name in ("study.yaml", "AllrunCases", "generate_cases.py"):
        assert name not in leaked, f"{name} leaked into Cases/"


@needs_casefoam
def test_casefoams_own_artefacts_are_left_alone(study_root, study):
    """Cases/baseCase, and the rmCases/Allrun/Allclean CaseFoam writes into
    Cases/, are its own layout -- the same one cases/caseFoamEx has. They are
    not cleaned up."""
    clone_cases(study_root, study, study.enabled_cases())
    assert (study_root / "Cases" / "baseCase").is_dir()


@needs_casefoam
def test_regenerating_replaces_the_previous_tree(study_root, study):
    """mkCases swallows FileExistsError from its own copytree, so a stale Cases/
    would otherwise be restructured in place."""
    clone_cases(study_root, study, study.enabled_cases())
    (study_root / "Cases" / "gap06in" / "p005psi" / "stale.txt").write_text(
        "left over\n", encoding="utf-8")

    clone_cases(study_root, study, study.enabled_cases())
    assert not (study_root / "Cases" / "gap06in" / "p005psi" / "stale.txt").exists()


@needs_casefoam
def test_generated_cases_do_not_inherit_a_dictionary_from_the_template(
        study_root, study):
    """CaseFoam copies the template faithfully, as it should. A dictionary left
    in baseCase by someone running ./Allmesh there would otherwise let a case run
    against the template's physics instead of its own."""
    (study_root / "baseCase" / "system" / "controlDict").write_text(
        "// stale\n", encoding="utf-8")
    (study_root / "baseCase" / "constant").mkdir(exist_ok=True)
    (study_root / "baseCase" / "constant" / "dsmcProperties").write_text(
        "// stale\n", encoding="utf-8")

    generated = clone_cases(study_root, study, study.enabled_cases())

    for path in generated:
        assert not (path / "system" / "controlDict").exists()
        assert not (path / "constant" / "dsmcProperties").exists()
        # A static dictionary the case does not generate is kept.
        assert (path / "system" / "fvSchemes").is_file()


@needs_casefoam
def test_a_missing_template_is_an_error(tmp_path, study):
    root = tmp_path / "empty"
    root.mkdir()
    with pytest.raises(StudyError, match="does not exist"):
        clone_cases(root, study, study.enabled_cases())

    (root / "baseCase").mkdir()
    with pytest.raises(StudyError, match="no case.yaml"):
        clone_cases(root, study, study.enabled_cases())


def test_pruning_removes_what_a_case_must_not_inherit(tmp_path):
    """The hygiene step, tested directly so it needs no CaseFoam."""
    case = tmp_path / "case"
    (case / "system").mkdir(parents=True)
    (case / "constant" / "polyMesh").mkdir(parents=True)
    (case / "0").mkdir()
    (case / "case.yaml").write_text("model: x\n", encoding="utf-8")
    (case / "system" / "controlDict").write_text("// generated\n", encoding="utf-8")
    (case / "system" / "fvSchemes").write_text("// static\n", encoding="utf-8")
    (case / "constant" / "dsmcProperties").write_text("// generated\n", encoding="utf-8")

    prune_generated_dictionaries(case)

    assert not (case / "system" / "controlDict").exists()
    assert not (case / "constant" / "dsmcProperties").exists()
    assert not (case / "constant" / "polyMesh").exists()
    assert not (case / "0").exists()
    assert (case / "case.yaml").is_file()
    assert (case / "system" / "fvSchemes").is_file()


# --------------------------------------------------------------------------- #
# applying case parameters
# --------------------------------------------------------------------------- #

def test_parameters_are_applied_structurally(tmp_path, study):
    """Structured YAML editing, not string substitution: the value is replaced
    whatever its original formatting."""
    destination = tmp_path / "p025psi"
    destination.mkdir()
    (destination / "case.yaml").write_text(
        "model: markelov1999_axisymmetric\n"
        "stagnation: {p0_pa: 34473.79}\n"
        "dsmc: {n_equivalent_particles: null}\n"
        "meta: {reference: AIAA 99-3455}\n", encoding="utf-8")

    applied = apply_case_parameters(
        destination, PressureCase("p025psi", 25.0), study,
        n_equivalent_particles=1.5e11)

    data = yaml.safe_load((destination / "case.yaml").read_text(encoding="utf-8"))
    assert data["stagnation"]["p0_pa"] == pytest.approx(25.0 * PSI_TO_PA, rel=1e-15)
    assert data["dsmc"]["n_equivalent_particles"] == pytest.approx(1.5e11)
    assert applied["pressure_psi"] == 25.0

    meta = data["meta"]
    assert meta["case_name"] == "p025psi"
    assert meta["gap_in"] == 6.0
    assert meta["reference"] == "AIAA 99-3455"   # inherited, not overwritten

    assert (destination / "case.yaml").read_text(
        encoding="utf-8").startswith("# GENERATED by")


def test_applying_to_a_file_without_the_expected_sections_is_an_error(tmp_path, study):
    """Rather than growing the sections silently, which would hide that the file
    is not this case family's template."""
    destination = tmp_path / "wrong"
    destination.mkdir()
    (destination / "case.yaml").write_text("model: source_flow\n", encoding="utf-8")
    with pytest.raises(StudyError, match="no 'stagnation' section"):
        apply_case_parameters(destination, PressureCase("x", 5.0), study)


def test_leaving_the_weight_unset_keeps_it_unset(tmp_path, study):
    destination = tmp_path / "p005psi"
    destination.mkdir()
    (destination / "case.yaml").write_text(
        "stagnation: {p0_pa: 1}\ndsmc: {n_equivalent_particles: null}\nmeta: {}\n",
        encoding="utf-8")
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
    path = write_manifest(tmp_path / "manifest.yaml", study, entries)

    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert document["n_generated"] == 4
    assert document["reference"] == "AIAA 99-3455"
    assert document["gap_in"] == 6.0
    assert document["generator"].startswith("casefoam ")
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
    path = write_manifest(tmp_path / "manifest.yaml", study, [])

    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert document["disabled_cases"] == [
        {"name": "p475psi", "pressure_psi": 475.0, "note": "too expensive"}]


def test_manifest_says_it_is_generated(tmp_path, study):
    path = write_manifest(tmp_path / "manifest.yaml", study, [])
    assert path.read_text(encoding="utf-8").startswith("# GENERATED by")


def test_the_generator_version_is_recorded(tmp_path, study):
    """The generated tree is CaseFoam's output, so the manifest says which
    version produced what is on disk."""
    path = write_manifest(tmp_path / "manifest.yaml", study, [])
    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert document["generator"].startswith("casefoam ")


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
