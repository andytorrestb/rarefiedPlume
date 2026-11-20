"""Generation of the cai2012-unigas family and its equivalence to the kn-extension.

The family inherits everything physical and geometric from
cases/cai2012-kn-extension: the same case.yaml outside ``meta``, the same
study matrix, byte-identical mesh dictionaries and the same particle weight
(spec §3). Generation is idempotent and never touches a case with results.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
import yaml

from plumetools.cai2012.config import load_case_config
from plumetools.cai2012.study import load_study
from plumetools.cai2012ugf.generation import generate, stage_pilot, validate
from plumetools.cai2012ugf.policy import study_differences, template_differences
from plumetools.cai2012ugf.settings import load_settings

REPO = Path(__file__).resolve().parents[2]
UGF = REPO / "cases" / "cai2012-unigas"
EXT = REPO / "cases" / "cai2012-kn-extension"
SEVEN = ("Kn100", "Kn10", "Kn1", "Kn0p1", "Kn0p01", "Kn0p001", "Kn0p0001")


def copy_family(tmp_path: Path) -> Path:
    target = tmp_path / "cai2012-unigas"
    shutil.copytree(UGF, target, ignore=shutil.ignore_patterns(
        "Cases", "Pilot", "results", "manifest-*.yaml", "__pycache__"))
    return target


def edit(path: Path, old: str, new: str) -> None:
    text = path.read_text(encoding="utf-8")
    assert old in text, f"{old!r} not in {path}"
    path.write_text(text.replace(old, new, 1), encoding="utf-8")


@pytest.fixture(autouse=True)
def no_git(monkeypatch):
    """The manifest's git provenance is a ~30 s `git status` on a Windows mount."""
    from plumetools.cai2012ugf import generation
    monkeypatch.setattr(generation, "git_provenance", lambda root: {"stub": "test"})


@pytest.fixture
def generated(tmp_path):
    family = copy_family(tmp_path)
    states = generate(family)
    assert {s.status for s in states} == {"new"}
    return family


def test_template_matches_kn_extension_outside_meta():
    assert template_differences(load_case_config(UGF / "baseCase"),
                                load_case_config(EXT / "baseCase")) == []


def test_template_difference_is_named(tmp_path):             # Review Focus 4
    family = copy_family(tmp_path)
    edit(family / "baseCase" / "case.yaml", "courant_target: 0.2", "courant_target: 0.3")
    assert template_differences(load_case_config(family / "baseCase"),
                                load_case_config(EXT / "baseCase")) == ["dsmc.courant_target"]


def test_study_matches_kn_extension():
    assert study_differences(load_study(UGF / "study.yaml"), load_study(EXT / "study.yaml")) == []


def test_study_difference_is_named(tmp_path):
    family = copy_family(tmp_path)
    edit(family / "study.yaml", "core_cell_size_m: 0.02", "core_cell_size_m: 0.03")
    assert any("Kn100" in d for d in study_differences(load_study(family / "study.yaml"),
                                                       load_study(EXT / "study.yaml")))


def test_seven_cases_with_the_suite_names(generated):
    assert sorted(p.name for p in (generated / "Cases").iterdir()) == sorted(
        f"{k}_np1x" for k in SEVEN)


def test_dry_run_reports_the_inflow_limitation_at_low_kn(tmp_path):
    states = {s.case.knudsen_case.name: s for s in generate(copy_family(tmp_path), dry_run=True)}
    flagged = sorted(k for k, s in states.items() if s.known_limitations)
    assert flagged == ["Kn0p0001", "Kn0p001"]


def test_each_case_has_unigas_yaml(generated):
    for case in (generated / "Cases").iterdir():
        assert load_settings(case).cell_weighted is False


def test_case_meta_names_this_family_and_solver(generated):
    data = yaml.safe_load((generated / "Cases/Kn1_np1x/case.yaml").read_text(encoding="utf-8"))
    assert data["meta"]["family"] == "cases/cai2012-unigas"
    assert data["meta"]["solver"] == "uniGasFoam hybrid USP-DSMC"
    assert data["meta"]["generated_by"] == "plumetools.cai2012ugf.generation"


def test_manifest_records_solver_and_defers_the_verdict(generated):
    manifest = yaml.safe_load((generated / "manifest-unigas.yaml").read_text(encoding="utf-8"))
    assert manifest["solver"] == "uniGasFoam hybrid USP-DSMC"
    assert {e["solver"] for e in manifest["cases"]} == {"uniGasFoam hybrid USP-DSMC"}
    assert {e["production_status"] for e in manifest["cases"]} == {"post-run-verdict-pending"}
    assert all("dsmcfoam_prerun_feasibility" in e and "unigas_commit" in e for e in manifest["cases"])


def test_mesh_dictionaries_byte_identical_and_weight_equal(generated):
    assert validate(generated) == []          # includes both comparisons, all seven Kn


def test_hand_edited_case_yaml_is_named(generated):          # Review I3
    case_yaml = generated / "Cases" / "Kn1_np1x" / "case.yaml"
    data = yaml.safe_load(case_yaml.read_text(encoding="utf-8"))
    data["dsmc"]["sampling_domain_transits"] = 2.0
    case_yaml.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    assert any("Kn1_np1x" in p and "dsmc.sampling_domain_transits" in p for p in validate(generated))


def test_second_generation_keeps_everything(generated):
    states = generate(generated)
    assert {s.status for s in states} == {"identical"}


def test_unigas_yaml_edit_is_an_input_change(generated):
    edit(generated / "baseCase/unigas.yaml", "max_courant_number: 0.5", "max_courant_number: 0.4")
    assert {s.status for s in generate(generated, dry_run=True)} == {"conflict"}


def test_case_with_results_is_never_replaced(generated):
    case = generated / "Cases/Kn1_np1x"
    (case / "0.06").mkdir()
    (case / "0.06/rhoN_Ar").write_text("x")
    edit(generated / "baseCase/unigas.yaml", "theta: 0.2", "theta: 0.3")
    states = {s.case.name: s for s in generate(generated, regenerate=True)}
    assert states["Kn1_np1x"].status == "conflict" and states["Kn1_np1x"].action == "refuse"
    assert (case / "0.06/rhoN_Ar").read_text() == "x"


def test_pilot_is_outside_the_manifest(generated):
    pilot = stage_pilot(generated, "Kn0p001")
    assert (pilot / "PILOT").is_file() and pilot.parent.name == "Pilot"
    assert "Pilot" not in (generated / "manifest-unigas.yaml").read_text()
