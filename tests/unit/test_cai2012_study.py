"""``study.yaml`` parsing, case generation, and the manifest.

The property under test throughout is **traceability**: a case named for one
Knudsen number must not be able to run at another, and the manifest must record
every value the solver will use.
"""

from __future__ import annotations

import shutil
import textwrap
from pathlib import Path

import pytest
import yaml

from _cai2012 import derived, write_case
from plumetools.cai2012 import geometry as geometry_module
from plumetools.cai2012 import inflow as inflow_module
from plumetools.cai2012 import mesh as mesh_module
from plumetools.cai2012.config import load_case_config
from plumetools.cai2012.study import (
    KnudsenCase,
    StudyError,
    apply_case_parameters,
    case_is_complete,
    clone_cases,
    filter_by_knudsen,
    load_study,
    manifest_entry,
    prune_generated_dictionaries,
    write_manifest,
)

REPO = Path(__file__).resolve().parents[2]
STUDY_DIR = REPO / "cases" / "cai2012"

STUDY_YAML = """\
base_case: baseCase
cases_dir: Cases
kn_cases:
  - {name: Kn100,  Kn: 100.0, enabled: true}
  - {name: Kn0p1,  Kn: 0.1,   enabled: true}
  - {name: Kn0p01, Kn: 0.01,  enabled: true}
"""


def write_study(tmp_path: Path, text: str = STUDY_YAML) -> Path:
    path = tmp_path / "study.yaml"
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture
def study(tmp_path):
    return load_study(write_study(tmp_path))


# --------------------------------------------------------------------------- #
# the shipped study
# --------------------------------------------------------------------------- #

def test_the_shipped_study_defines_cais_three_knudsen_numbers():
    study = load_study(STUDY_DIR)
    assert sorted(c.knudsen for c in study.all_cases()) == [0.01, 0.1, 100.0]


def test_the_shipped_study_enables_all_three():
    assert len(load_study(STUDY_DIR).enabled_cases()) == 3


def test_the_shipped_study_records_the_paper_metadata():
    meta = load_study(STUDY_DIR).meta
    assert meta["reference"] == "Cai and Wang 2012"
    assert meta["nozzle_diameter_m"] == 0.2
    assert meta["speed_ratio"] == 2.0


# --------------------------------------------------------------------------- #
# parsing
# --------------------------------------------------------------------------- #

def test_a_directory_or_a_file_is_accepted(tmp_path):
    write_study(tmp_path)
    assert len(load_study(tmp_path).kn_cases) == 3
    assert len(load_study(tmp_path / "study.yaml").kn_cases) == 3


def test_a_missing_study_is_an_error(tmp_path):
    with pytest.raises(StudyError, match="no study.yaml"):
        load_study(tmp_path)


def test_an_unknown_top_level_key_is_rejected(tmp_path):
    with pytest.raises(StudyError, match="unknown top-level key"):
        load_study(write_study(tmp_path, STUDY_YAML + "pressure_cases: []\n"))


def test_an_unknown_case_key_is_rejected(tmp_path):
    text = "kn_cases:\n  - {name: a, Kn: 1.0, presure: 2}\n"
    with pytest.raises(StudyError, match="unknown key"):
        load_study(write_study(tmp_path, text))


def test_a_missing_kn_is_rejected(tmp_path):
    with pytest.raises(StudyError, match="missing 'Kn'"):
        load_study(write_study(tmp_path, "kn_cases:\n  - {name: a}\n"))


def test_an_empty_matrix_is_rejected(tmp_path):
    with pytest.raises(StudyError, match="non-empty list"):
        load_study(write_study(tmp_path, "kn_cases: []\n"))


@pytest.mark.parametrize("kn", ["0.0", "-1.0"])
def test_a_non_positive_knudsen_is_rejected(tmp_path, kn):
    """n0 = 1/(sqrt(2) pi d^2 Kn D) would be negative or infinite."""
    with pytest.raises(StudyError, match="Knudsen number is positive"):
        load_study(write_study(tmp_path, f"kn_cases:\n  - {{name: a, Kn: {kn}}}\n"))


def test_a_duplicate_case_name_is_rejected(tmp_path):
    text = "kn_cases:\n  - {name: a, Kn: 1.0}\n  - {name: a, Kn: 2.0}\n"
    with pytest.raises(StudyError, match="duplicate case name"):
        load_study(write_study(tmp_path, text))


def test_a_duplicate_knudsen_number_is_rejected(tmp_path):
    """Two cases at the same Kn would differ in name only."""
    text = "kn_cases:\n  - {name: a, Kn: 1.0}\n  - {name: b, Kn: 1.0}\n"
    with pytest.raises(StudyError, match="duplicate Knudsen"):
        load_study(write_study(tmp_path, text))


def test_overriding_the_knudsen_number_is_rejected(tmp_path):
    """The Knudsen number is what a case IS. Overriding it would let a case
    named for one Kn run at another."""
    text = textwrap.dedent("""\
        kn_cases:
          - name: a
            Kn: 1.0
            overrides:
              exit: {knudsen: 5.0}
        """)
    with pytest.raises(StudyError, match="may not set 'exit.knudsen'"):
        load_study(write_study(tmp_path, text))


def test_overriding_the_model_is_rejected(tmp_path):
    text = textwrap.dedent("""\
        kn_cases:
          - name: a
            Kn: 1.0
            overrides:
              model: something_else
        """)
    with pytest.raises(StudyError, match="may not set 'model'"):
        load_study(write_study(tmp_path, text))


def test_other_overrides_are_allowed(tmp_path):
    text = textwrap.dedent("""\
        kn_cases:
          - name: a
            Kn: 1.0
            overrides:
              mesh: {max_cells: 100}
        """)
    study = load_study(write_study(tmp_path, text))
    assert study.kn_cases[0].overrides == {"mesh": {"max_cells": 100}}


# --------------------------------------------------------------------------- #
# ordering and filtering
# --------------------------------------------------------------------------- #

def test_cases_come_out_most_rarefied_first(study):
    """Decreasing Kn: the collisionless benchmark first, which is the order the
    results have to be trusted in."""
    assert [c.name for c in study.enabled_cases()] == ["Kn100", "Kn0p1", "Kn0p01"]


def test_ordering_does_not_depend_on_file_order(tmp_path):
    shuffled = textwrap.dedent("""\
        kn_cases:
          - {name: Kn0p01, Kn: 0.01}
          - {name: Kn100,  Kn: 100.0}
          - {name: Kn0p1,  Kn: 0.1}
        """)
    study = load_study(write_study(tmp_path, shuffled))
    assert [c.name for c in study.enabled_cases()] == ["Kn100", "Kn0p1", "Kn0p01"]


def test_a_disabled_case_stays_in_the_matrix(tmp_path):
    text = "kn_cases:\n  - {name: a, Kn: 1.0, enabled: false}\n" \
           "  - {name: b, Kn: 2.0}\n"
    study = load_study(write_study(tmp_path, text))
    assert len(study.all_cases()) == 2
    assert [c.name for c in study.enabled_cases()] == ["b"]


def test_filtering_keeps_the_requested_cases(study):
    kept = filter_by_knudsen(study.enabled_cases(), [100.0])
    assert [c.name for c in kept] == ["Kn100"]


def test_filtering_by_nothing_keeps_everything(study):
    assert len(filter_by_knudsen(study.enabled_cases(), None)) == 3


def test_filtering_by_an_absent_knudsen_is_an_error(study):
    """A typo on the command line stops the run instead of generating nothing."""
    with pytest.raises(StudyError, match="no case at Kn"):
        filter_by_knudsen(study.enabled_cases(), [1.0])


def test_case_paths_are_flat(study):
    case = study.enabled_cases()[0]
    assert study.case_path(case) == Path("Cases") / "Kn100"


# --------------------------------------------------------------------------- #
# generation
# --------------------------------------------------------------------------- #

@pytest.fixture
def template(tmp_path):
    """A study directory with a baseCase and a study.yaml."""
    write_study(tmp_path)
    write_case(tmp_path / "baseCase")
    (tmp_path / "baseCase" / "Allmesh").write_text("#!/bin/bash\n", encoding="utf-8")
    return tmp_path


def test_clone_creates_one_directory_per_case(template):
    study = load_study(template)
    paths = clone_cases(template, study, study.enabled_cases())
    assert [p.name for p in paths] == ["Kn100", "Kn0p1", "Kn0p01"]
    for path in paths:
        assert (path / "case.yaml").is_file()
        assert (path / "Allmesh").is_file()


def test_clone_without_a_base_case_is_an_error(tmp_path):
    write_study(tmp_path)
    study = load_study(tmp_path)
    with pytest.raises(StudyError, match="does not exist"):
        clone_cases(tmp_path, study, study.enabled_cases())


def test_clone_prunes_inherited_dictionaries(template):
    """A dictionary left in baseCase by someone running ./Allmesh there would
    otherwise be inherited by every case, and ./Allrun only checks that
    constant/dsmcProperties EXISTS."""
    (template / "baseCase" / "constant").mkdir(parents=True, exist_ok=True)
    (template / "baseCase" / "constant" / "dsmcProperties").write_text(
        "stale", encoding="utf-8")
    (template / "baseCase" / "constant" / "polyMesh").mkdir()
    (template / "baseCase" / "0").mkdir()

    study = load_study(template)
    paths = clone_cases(template, study, study.enabled_cases())
    for path in paths:
        assert not (path / "constant" / "dsmcProperties").exists()
        assert not (path / "constant" / "polyMesh").exists()
        assert not (path / "0").exists()


def test_prune_removes_logs(tmp_path):
    case = tmp_path / "case"
    case.mkdir()
    (case / "log.dsmcFoam").write_text("x", encoding="utf-8")
    removed = prune_generated_dictionaries(case)
    assert not (case / "log.dsmcFoam").exists()
    assert any(p.name == "log.dsmcFoam" for p in removed)


def test_apply_writes_the_knudsen_number_structurally(template):
    study = load_study(template)
    case = study.enabled_cases()[2]      # Kn0p01
    clone_cases(template, study, study.enabled_cases())
    destination = template / study.case_path(case)

    applied = apply_case_parameters(destination, case, study)
    assert applied["knudsen"] == 0.01

    data = yaml.safe_load((destination / "case.yaml").read_text(encoding="utf-8"))
    assert data["exit"]["knudsen"] == 0.01
    assert data["meta"]["case_name"] == "Kn0p01"


def test_apply_leaves_a_loadable_case(template):
    study = load_study(template)
    clone_cases(template, study, study.enabled_cases())
    for case in study.enabled_cases():
        destination = template / study.case_path(case)
        apply_case_parameters(destination, case, study)
        cfg = load_case_config(destination)
        assert cfg.exit.knudsen == case.knudsen


def test_a_generated_case_derives_its_own_density(template):
    """Nothing stores n0. The three cases end up with three densities because
    they have three Knudsen numbers, and for no other reason."""
    study = load_study(template)
    clone_cases(template, study, study.enabled_cases())
    densities = []
    for case in study.enabled_cases():
        destination = template / study.case_path(case)
        apply_case_parameters(destination, case, study)
        cfg = load_case_config(destination)
        densities.append(inflow_module.from_config(cfg).number_density_per_m3)
    assert len(set(densities)) == 3


def test_overrides_are_applied(template):
    study = load_study(template)
    case = KnudsenCase(name="Kn100", knudsen=100.0,
                       overrides={"dsmc": {"n_subdomains": 9}})
    clone_cases(template, study, [case])
    destination = template / study.case_path(case)
    applied = apply_case_parameters(destination, case, study)
    assert applied["overrides"] == {"dsmc": {"n_subdomains": 9}}
    assert load_case_config(destination).dsmc.n_subdomains == 9


def test_an_override_of_a_key_the_template_lacks_is_rejected(template):
    """Better a StudyError naming study.yaml than a config error naming the
    generated file."""
    study = load_study(template)
    case = KnudsenCase(name="Kn100", knudsen=100.0,
                       overrides={"dsmc": {"nsubdomains": 9}})
    clone_cases(template, study, [case])
    with pytest.raises(StudyError, match="not in the base case.yaml"):
        apply_case_parameters(template / study.case_path(case), case, study)


def test_apply_refuses_a_foreign_template(tmp_path):
    """A case.yaml without an 'exit' section is not a Cai case, and silently
    growing one would hide that."""
    (tmp_path / "case.yaml").write_text("model: x\nmeta: {}\n", encoding="utf-8")
    study = load_study(write_study(tmp_path))
    with pytest.raises(StudyError, match="does not look like the cai2012 template"):
        apply_case_parameters(tmp_path, KnudsenCase("a", 1.0), study)


def test_case_is_complete_detects_solver_output(tmp_path):
    case = tmp_path / "case"
    (case / "0").mkdir(parents=True)
    assert not case_is_complete(case)
    (case / "0.05").mkdir()
    assert case_is_complete(case)


# --------------------------------------------------------------------------- #
# the manifest
# --------------------------------------------------------------------------- #

def test_manifest_entry_carries_every_required_value(tmp_path):
    """§7: Kn, D, lambda0, T0, U0, n0, particle weight, minimum cell size, deltaT."""
    cfg, geom, state, plan, run = derived(tmp_path)
    study = load_study(write_study(tmp_path))
    entry = manifest_entry(KnudsenCase("Kn100", 100.0), study, cfg, geom, state,
                           plan, run)
    for key in ("Kn", "D_m", "lambda0_m", "T0_K", "U0_m_per_s", "n0_per_m3",
                "n_equivalent_particles", "min_cell_size_m", "deltaT_s"):
        assert key in entry, key
    assert entry["Kn"] == 100.0
    assert entry["n0_per_m3"] == pytest.approx(state.number_density_per_m3)


def test_manifest_entry_records_a_coarsened_mesh(tmp_path):
    cfg, geom, state, plan, run = derived(
        tmp_path, exit={"knudsen": 0.001},
        mesh={"core_cell_size_m": None, "max_cells": 200000})
    study = load_study(write_study(tmp_path))
    entry = manifest_entry(KnudsenCase("x", 0.001), study, cfg, geom, state,
                           plan, run)
    assert entry["coarsened_for_budget"] is True
    assert entry["cell_over_mean_free_path"] > 1.0


def test_manifest_records_disabled_cases_as_a_decision(tmp_path):
    text = "kn_cases:\n  - {name: a, Kn: 1.0, enabled: false}\n" \
           "  - {name: b, Kn: 2.0}\n"
    study = load_study(write_study(tmp_path, text))
    path = write_manifest(tmp_path / "manifest.yaml", study, [{"name": "b"}])
    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert [c["name"] for c in document["disabled_cases"]] == ["a"]
    assert document["n_generated"] == 1


def test_manifest_is_lf_and_says_it_is_generated(tmp_path):
    study = load_study(write_study(tmp_path))
    path = write_manifest(tmp_path / "manifest.yaml", study, [])
    assert b"\r\n" not in path.read_bytes()
    assert "GENERATED" in path.read_text(encoding="utf-8")
