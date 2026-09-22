"""The numerical-particle axis: ``nEquivalentParticles`` and nothing else.

The property under test throughout is that the sweep is a **controlled**
experiment. A convergence study whose members also differ in the mesh, the time
step or the density does not measure statistical resolution, and it fails to
measure it after the HPC hours have been spent -- so the tests here are mostly
about what must NOT change, and the shipped study is checked against the matrix
it had before the axis existed.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest
import yaml

from _cai2012 import derived, write_case
from plumetools.cai2012 import geometry as geometry_module
from plumetools.cai2012 import inflow as inflow_module
from plumetools.cai2012 import mesh as mesh_module
from plumetools.cai2012.config import load_case_config
from plumetools.cai2012.dictionaries import (
    AVERAGED_FIELDS,
    render_control_dict,
    render_dsmc_properties,
)
from plumetools.cai2012.study import (
    ParticleLevel,
    StudyError,
    apply_case_parameters,
    case_is_complete,
    clone_cases,
    filter_by_particles,
    load_study,
    manifest_entry,
)
from plumetools.cai2012.verify import (
    field_average_block,
    foam_float,
    verify_tree,
)

REPO = Path(__file__).resolve().parents[2]
STUDY_DIR = REPO / "cases" / "cai2012"

#: The levels the study asks for.
LEVELS = (1.0, 2.0, 5.0, 10.0)


def derive_generated(case_dir: Path):
    """``(cfg, geom, exit_state, plan, run)`` from a case.yaml already on disk.

    ``_cai2012.derived`` *writes* a fresh minimal case.yaml first, which is what
    makes it convenient everywhere else and useless here: it would overwrite
    exactly the file :func:`apply_case_parameters` just produced. This loads
    what generation wrote, which is also what ``./Allmesh`` does.
    """
    cfg = load_case_config(case_dir)
    geom = geometry_module.from_config(cfg)
    exit_state = inflow_module.from_config(cfg)
    plan = mesh_module.plan(cfg, geom, exit_state)
    run = inflow_module.derive_run_settings(
        cfg, exit_state, geom,
        min_cell_size_m=plan.min_cell_size_m,
        exit_cell_volume_m3=plan.exit_cell_volume_m3)
    return cfg, geom, exit_state, plan, run


# --------------------------------------------------------------------------- #
# the shipped study: the physical matrix is UNCHANGED
# --------------------------------------------------------------------------- #

def test_the_shipped_study_still_defines_cais_three_knudsen_numbers():
    """§14.2. The physical matrix is the baseline for this work, and adding an
    axis to it must not have edited it."""
    study = load_study(STUDY_DIR)
    assert sorted(c.knudsen for c in study.all_cases()) == [0.01, 0.1, 100.0]
    assert [c.name for c in study.enabled_cases()] == ["Kn100", "Kn0p1", "Kn0p01"]


def test_the_shipped_physical_cases_keep_their_overrides():
    """The transient basis, the sampling duration and the subdomain count are
    per-case decisions with reasons recorded in study.yaml. The particle axis
    does not get to change them."""
    study = load_study(STUDY_DIR)
    for case in study.enabled_cases():
        assert case.overrides["dsmc"]["transient_basis"] == "transits"
        assert case.overrides["dsmc"]["transient_domain_transits"] == 2.0
        assert case.overrides["dsmc"]["sampling_domain_transits"] == 1.5
        assert case.overrides["dsmc"]["n_subdomains"] == 12
    kn0p01 = [c for c in study.enabled_cases() if c.knudsen == 0.01][0]
    assert kn0p01.overrides["mesh"]["max_cells"] == 4000000


def test_the_shipped_study_sweeps_the_four_particle_levels():
    study = load_study(STUDY_DIR)
    assert [p.multiplier for p in study.enabled_particle_levels()] == list(LEVELS)
    assert [p.name for p in study.enabled_particle_levels()] == [
        "np1x", "np2x", "np5x", "np10x"]


def test_the_shipped_study_asks_for_3x_run_time_and_4x_output():
    study = load_study(STUDY_DIR)
    assert study.run_time_multiplier == 3.0
    assert study.output_frequency_multiplier == 4.0


def test_every_physical_case_gets_four_variants():
    """§14.3. 3 physical x 4 levels = 12."""
    study = load_study(STUDY_DIR)
    expanded = study.expanded_cases()
    assert len(expanded) == 12
    for physical in study.enabled_cases():
        variants = [c for c in expanded if c.knudsen_case is physical]
        assert [c.multiplier for c in variants] == list(LEVELS)


def test_generated_names_and_paths_are_unique():
    """§14.13. Two cases sharing a directory would silently overwrite."""
    study = load_study(STUDY_DIR)
    expanded = study.expanded_cases()
    names = [c.name for c in expanded]
    paths = [str(study.case_path(c)) for c in expanded]
    assert len(set(names)) == len(names)
    assert len(set(paths)) == len(paths)


def test_a_generated_name_carries_both_axes():
    """§14.14. The physical case and the level are both readable from the name."""
    study = load_study(STUDY_DIR)
    assert "Kn0p01_np10x" in [c.name for c in study.expanded_cases()]


def test_variants_of_a_case_are_adjacent():
    """A case's four variants are compared with each other, so they are
    generated and run next to each other."""
    study = load_study(STUDY_DIR)
    names = [c.name for c in study.expanded_cases()]
    assert names[:4] == ["Kn100_np1x", "Kn100_np2x", "Kn100_np5x", "Kn100_np10x"]


# --------------------------------------------------------------------------- #
# parsing the axis
# --------------------------------------------------------------------------- #

STUDY_YAML = textwrap.dedent("""\
    base_case: baseCase
    cases_dir: Cases
    run_time_multiplier: 3.0
    output_frequency_multiplier: 4.0
    kn_cases:
      - {name: Kn100,  Kn: 100.0}
      - {name: Kn0p1,  Kn: 0.1}
    particle_levels:
      - {name: np1x,  multiplier: 1.0}
      - {name: np2x,  multiplier: 2.0}
      - {name: np5x,  multiplier: 5.0}
      - {name: np10x, multiplier: 10.0}
    """)


def write_study(tmp_path: Path, text: str = STUDY_YAML) -> Path:
    path = tmp_path / "study.yaml"
    path.write_text(text, encoding="utf-8")
    return path


def test_a_study_without_a_particle_axis_is_the_study_as_it_was(tmp_path):
    """Backwards compatibility, and the definition of the baseline: no axis
    means one 1x level, which divides the weight by one."""
    study = load_study(write_study(
        tmp_path, "kn_cases:\n  - {name: a, Kn: 1.0}\n"))
    assert [p.multiplier for p in study.enabled_particle_levels()] == [1.0]
    assert study.run_time_multiplier == 1.0
    assert study.output_frequency_multiplier == 1.0
    assert len(study.expanded_cases()) == 1


def test_a_matrix_without_a_1x_level_is_rejected(tmp_path):
    """Without the control there is nothing to read the other levels against,
    and nothing that reproduces the published results."""
    text = ("kn_cases:\n  - {name: a, Kn: 1.0}\n"
            "particle_levels:\n  - {name: np2x, multiplier: 2.0}\n")
    with pytest.raises(StudyError, match="no 1x level"):
        load_study(write_study(tmp_path, text))


@pytest.mark.parametrize("multiplier", ["0.0", "-2.0"])
def test_a_non_positive_multiplier_is_rejected(tmp_path, multiplier):
    """It DIVIDES the weight, so zero or negative is an infinite or negative
    number of molecules per parcel."""
    text = ("kn_cases:\n  - {name: a, Kn: 1.0}\n"
            "particle_levels:\n  - {name: np1x, multiplier: 1.0}\n"
            f"  - {{name: bad, multiplier: {multiplier}}}\n")
    with pytest.raises(StudyError, match="DIVIDES nEquivalentParticles"):
        load_study(write_study(tmp_path, text))


def test_a_duplicate_level_name_is_rejected(tmp_path):
    text = ("kn_cases:\n  - {name: a, Kn: 1.0}\n"
            "particle_levels:\n  - {name: np1x, multiplier: 1.0}\n"
            "  - {name: np1x, multiplier: 2.0}\n")
    with pytest.raises(StudyError, match="duplicate particle level name"):
        load_study(write_study(tmp_path, text))


def test_a_duplicate_multiplier_is_rejected(tmp_path):
    text = ("kn_cases:\n  - {name: a, Kn: 1.0}\n"
            "particle_levels:\n  - {name: np1x, multiplier: 1.0}\n"
            "  - {name: same, multiplier: 1.0}\n")
    with pytest.raises(StudyError, match="duplicate particle multiplier"):
        load_study(write_study(tmp_path, text))


def test_an_unknown_level_key_is_rejected(tmp_path):
    text = ("kn_cases:\n  - {name: a, Kn: 1.0}\n"
            "particle_levels:\n  - {name: np1x, multiplier: 1.0, weight: 3}\n")
    with pytest.raises(StudyError, match="unknown key"):
        load_study(write_study(tmp_path, text))


@pytest.mark.parametrize("key", ["numerical_particle_multiplier",
                                 "run_time_multiplier",
                                 "output_frequency_multiplier"])
def test_a_case_may_not_override_a_multiplier(tmp_path, key):
    """The multiplier is what the generated directory's suffix MEANS, and the
    other two are study-wide or the members are not comparable."""
    text = ("kn_cases:\n  - name: a\n    Kn: 1.0\n"
            f"    overrides:\n      dsmc: {{{key}: 7.0}}\n")
    with pytest.raises(StudyError, match="may not set"):
        load_study(write_study(tmp_path, text))


def test_levels_come_out_cheapest_first(tmp_path):
    """The 1x control first: if the sweep is cut short, what survives is the
    part comparable with the published results."""
    text = ("kn_cases:\n  - {name: a, Kn: 1.0}\n"
            "particle_levels:\n  - {name: np10x, multiplier: 10.0}\n"
            "  - {name: np1x, multiplier: 1.0}\n"
            "  - {name: np2x, multiplier: 2.0}\n")
    study = load_study(write_study(tmp_path, text))
    assert [p.name for p in study.enabled_particle_levels()] == [
        "np1x", "np2x", "np10x"]


def test_a_disabled_level_stays_in_the_matrix(tmp_path):
    text = ("kn_cases:\n  - {name: a, Kn: 1.0}\n"
            "particle_levels:\n  - {name: np1x, multiplier: 1.0}\n"
            "  - {name: np2x, multiplier: 2.0, enabled: false}\n")
    study = load_study(write_study(tmp_path, text))
    assert len(study.all_particle_levels()) == 2
    assert [p.name for p in study.enabled_particle_levels()] == ["np1x"]


def test_filtering_levels_accepts_a_name_or_a_multiplier(tmp_path):
    study = load_study(write_study(tmp_path))
    levels = study.enabled_particle_levels()
    assert [p.name for p in filter_by_particles(levels, ["np2x"])] == ["np2x"]
    assert [p.name for p in filter_by_particles(levels, [5])] == ["np5x"]
    assert [p.name for p in filter_by_particles(levels, None)] == [
        "np1x", "np2x", "np5x", "np10x"]


def test_filtering_by_an_absent_level_is_an_error(tmp_path):
    study = load_study(write_study(tmp_path))
    with pytest.raises(StudyError, match="no particle level"):
        filter_by_particles(study.enabled_particle_levels(), [3])


# --------------------------------------------------------------------------- #
# the weight itself
# --------------------------------------------------------------------------- #

def test_the_1x_case_runs_at_the_baseline_weight(tmp_path):
    """§14.4. Not 'close to' the value the study has always used -- the same
    object, because 1x divides by one."""
    _, _, _, _, baseline = derived(tmp_path / "one")
    _, _, _, _, run = derived(
        tmp_path / "two", dsmc={"numerical_particle_multiplier": 1.0})
    assert run.n_equivalent_particles == baseline.n_equivalent_particles


@pytest.mark.parametrize("multiplier", LEVELS)
def test_the_weight_is_the_baseline_divided_by_the_multiplier(tmp_path,
                                                              multiplier):
    """§14.5-7."""
    _, _, _, _, baseline = derived(tmp_path / "baseline")
    _, _, _, _, run = derived(
        tmp_path / f"np{multiplier:g}",
        dsmc={"numerical_particle_multiplier": multiplier})
    assert run.n_equivalent_particles == pytest.approx(
        baseline.n_equivalent_particles / multiplier, rel=1e-12)
    assert run.baseline_n_equivalent_particles == pytest.approx(
        baseline.n_equivalent_particles, rel=1e-12)


@pytest.mark.parametrize("multiplier", LEVELS)
def test_the_parcel_population_goes_up_by_the_multiplier(tmp_path, multiplier):
    """The point of the axis. One parcel stands for nEquivalentParticles
    molecules, so a smaller weight is more parcels."""
    _, _, _, _, baseline = derived(tmp_path / "baseline")
    _, _, _, _, run = derived(
        tmp_path / f"np{multiplier:g}",
        dsmc={"numerical_particle_multiplier": multiplier})
    assert run.exit_particles_per_cell == pytest.approx(
        baseline.exit_particles_per_cell * multiplier, rel=1e-12)


def test_the_multiplier_applies_to_a_pinned_weight_too(tmp_path):
    """A case that pins its weight is still entitled to a 2x variant, and the
    level is defined against whatever the 1x case runs at."""
    _, _, _, _, run = derived(
        tmp_path / "pinned",
        dsmc={"n_equivalent_particles": 1.0e+12,
              "numerical_particle_multiplier": 5.0})
    assert run.baseline_n_equivalent_particles == 1.0e+12
    assert run.n_equivalent_particles == pytest.approx(2.0e+11)


@pytest.mark.parametrize("key", ["numerical_particle_multiplier",
                                 "run_time_multiplier",
                                 "output_frequency_multiplier"])
def test_a_zero_multiplier_is_refused_by_the_loader(tmp_path, key):
    """Caught at load time, naming the key, rather than as a division by zero
    somewhere inside the derivation."""
    with pytest.raises(ValueError, match="finite and positive"):
        derived(tmp_path / "zero", dsmc={key: 0.0})


# --------------------------------------------------------------------------- #
# what must NOT move
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("multiplier", LEVELS)
def test_the_time_step_is_unchanged(tmp_path, multiplier):
    """§14.8. This is the whole controlled-experiment claim: vary the parcel
    count, hold the discretisation fixed."""
    _, _, _, _, baseline = derived(tmp_path / "baseline")
    _, _, _, _, run = derived(
        tmp_path / f"np{multiplier:g}",
        dsmc={"numerical_particle_multiplier": multiplier})
    assert run.delta_t_s == baseline.delta_t_s
    assert run.courant == baseline.courant


@pytest.mark.parametrize("multiplier", LEVELS)
def test_the_physics_and_the_mesh_are_unchanged(tmp_path, multiplier):
    """§14.9, §15. The density, the mean free path, the exit velocity and the
    mesh plan are all derived from Kn, and the particle axis is not in that
    chain."""
    _, _, base_exit, base_plan, _ = derived(tmp_path / "baseline")
    _, _, exit_state, plan, _ = derived(
        tmp_path / f"np{multiplier:g}",
        dsmc={"numerical_particle_multiplier": multiplier})
    assert exit_state.as_dict() == base_exit.as_dict()
    assert plan.n_cells == base_plan.n_cells
    assert plan.min_cell_size_m == base_plan.min_cell_size_m
    assert plan.core_cell_size_m == base_plan.core_cell_size_m


@pytest.mark.parametrize("multiplier", LEVELS)
def test_the_averaging_start_is_unchanged(tmp_path, multiplier):
    """§10. The four variants have to be averaged over the same window or the
    comparison between them is not about parcels."""
    _, _, _, _, baseline = derived(tmp_path / "baseline")
    _, _, _, _, run = derived(
        tmp_path / f"np{multiplier:g}",
        dsmc={"numerical_particle_multiplier": multiplier})
    assert run.average_start_s == baseline.average_start_s
    assert run.end_time_s == baseline.end_time_s
    assert run.write_interval_steps == baseline.write_interval_steps


# --------------------------------------------------------------------------- #
# run length and write frequency
# --------------------------------------------------------------------------- #

def test_the_end_time_is_extended_by_the_multiplier(tmp_path):
    """§14.10. Snapped up to a whole number of write intervals, so at or just
    above 3x and never below it."""
    _, _, _, _, baseline = derived(tmp_path / "baseline")
    _, _, _, _, run = derived(tmp_path / "long",
                              dsmc={"run_time_multiplier": 3.0})
    ratio = run.end_time_s / baseline.end_time_s
    assert 3.0 <= ratio <= 3.0 + run.write_interval_s / baseline.end_time_s
    assert run.baseline_end_time_s == pytest.approx(
        baseline.baseline_end_time_s)


def test_the_transient_does_not_move_with_the_run_length(tmp_path):
    """§7, §10. The extra time is SAMPLING time. A longer startup would be more
    expense for no more statistics."""
    _, _, _, _, baseline = derived(tmp_path / "baseline")
    _, _, _, _, run = derived(tmp_path / "long",
                              dsmc={"run_time_multiplier": 3.0})
    assert run.average_start_s == baseline.average_start_s
    # Every second the multiplier added landed inside the averaging window, so
    # the sampled window grew by MORE than the run did.
    assert (run.sampling_time_s / baseline.sampling_time_s
            > run.end_time_s / baseline.end_time_s)


def test_the_run_length_does_not_change_the_time_step(tmp_path):
    """§7. 'Do not change deltaT' -- including as a side effect of running
    longer."""
    _, _, _, _, baseline = derived(tmp_path / "baseline")
    _, _, _, _, run = derived(tmp_path / "long",
                              dsmc={"run_time_multiplier": 3.0})
    assert run.delta_t_s == baseline.delta_t_s


def test_output_is_four_times_more_frequent(tmp_path):
    """§14.11. writeControl is timeStep, so the interval is an integer number
    of steps and 4x is met to within that rounding."""
    _, _, _, _, baseline = derived(tmp_path / "baseline")
    _, _, _, _, run = derived(tmp_path / "often",
                              dsmc={"output_frequency_multiplier": 4.0})
    assert run.baseline_write_interval_steps == baseline.write_interval_steps
    assert run.write_interval_steps == max(
        1, round(baseline.write_interval_steps / 4.0))
    assert run.achieved_output_frequency_multiplier == pytest.approx(4.0, rel=0.02)


def test_the_write_interval_is_an_integer_number_of_steps(tmp_path):
    _, _, _, _, run = derived(tmp_path / "often",
                              dsmc={"output_frequency_multiplier": 4.0})
    assert isinstance(run.write_interval_steps, int)
    assert run.write_interval_steps >= 1


def test_a_huge_output_multiplier_still_writes_something(tmp_path):
    """Rounding down to a zero-step interval would mean writing nothing."""
    _, _, _, _, run = derived(tmp_path / "silly",
                              dsmc={"output_frequency_multiplier": 1.0e+6})
    assert run.write_interval_steps == 1


def test_the_output_frequency_does_not_change_the_time_step(tmp_path):
    """§8. 'Do not accidentally modify deltaT'."""
    _, _, _, _, baseline = derived(tmp_path / "baseline")
    _, _, _, _, run = derived(tmp_path / "often",
                              dsmc={"output_frequency_multiplier": 4.0})
    assert run.delta_t_s == baseline.delta_t_s


def test_the_write_interval_is_derived_from_the_baseline_run_not_the_long_one(
        tmp_path):
    """Otherwise 'four times more often than before' would be four times a
    number that itself moved with the run length."""
    _, _, _, _, baseline = derived(tmp_path / "baseline")
    _, _, _, _, run = derived(tmp_path / "both",
                              dsmc={"run_time_multiplier": 3.0,
                                    "output_frequency_multiplier": 4.0})
    assert run.baseline_write_interval_steps == baseline.write_interval_steps


def test_running_longer_gives_more_sampled_writes(tmp_path):
    """The point of §7 and §8 together: more frames inside the averaging
    window, which is what a statistical-convergence comparison reads."""
    _, _, _, _, baseline = derived(tmp_path / "baseline")
    _, _, _, _, run = derived(tmp_path / "both",
                              dsmc={"run_time_multiplier": 3.0,
                                    "output_frequency_multiplier": 4.0})
    assert run.n_sampled_writes > 10 * baseline.n_sampled_writes


# --------------------------------------------------------------------------- #
# the generated dictionaries
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("multiplier", LEVELS)
def test_dsmc_properties_carries_the_divided_weight(tmp_path, multiplier):
    """§14.12. What the solver reads, not what Python computed."""
    cfg, _, exit_state, _, run = derived(
        tmp_path / f"np{multiplier:g}",
        dsmc={"numerical_particle_multiplier": multiplier})
    text = render_dsmc_properties(cfg, exit_state, run)
    path = tmp_path / "dsmcProperties"
    path.write_text(text, encoding="utf-8")
    assert foam_float(path, "nEquivalentParticles") == pytest.approx(
        run.baseline_n_equivalent_particles / multiplier, rel=1e-6)


def test_the_control_dict_keeps_the_studys_write_control(tmp_path):
    """§8. 'Preserve the existing writeControl mode' -- it is timeStep for a
    documented reason, and a resumed run under runTime writes nothing."""
    cfg, _, _, _, run = derived(tmp_path / "case",
                                dsmc={"output_frequency_multiplier": 4.0})
    path = tmp_path / "controlDict"
    path.write_text(render_control_dict(cfg, run), encoding="utf-8")
    lines = [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines()]
    assert "writeControl    timeStep;" in lines
    assert f"writeInterval   {run.write_interval_steps};" in lines


def test_the_control_dict_averages_from_the_transient_not_from_zero(tmp_path):
    """§10. Averaging from t = 0 would fold the startup into every mean."""
    cfg, _, _, _, run = derived(tmp_path / "case",
                                dsmc={"run_time_multiplier": 3.0})
    path = tmp_path / "controlDict"
    path.write_text(render_control_dict(cfg, run), encoding="utf-8")
    block = field_average_block(path)
    assert block["timeStart"] == pytest.approx(run.average_start_s, rel=1e-5)
    assert block["timeStart"] > 0.0


def test_the_control_dict_keeps_the_variance_of_both_number_densities(tmp_path):
    """§9. CV = sqrt(dsmcRhoNPrime2Mean) / dsmcRhoNMean is the metric this
    study exists to make evaluable, and it cannot be recovered later from a
    mean alone."""
    cfg, _, _, _, run = derived(tmp_path / "case")
    path = tmp_path / "controlDict"
    path.write_text(render_control_dict(cfg, run), encoding="utf-8")
    fields = field_average_block(path)["fields"]

    assert set(fields) == set(AVERAGED_FIELDS)
    for name in AVERAGED_FIELDS:
        assert fields[name]["mean"] is True
    for name in ("rhoN", "dsmcRhoN"):
        assert fields[name]["prime2Mean"] is True


def test_prime2mean_stays_off_for_the_vector_fields(tmp_path):
    """A symmTensor per cell per write, that nothing reads, on a study that
    writes 60 times per case."""
    cfg, _, _, _, run = derived(tmp_path / "case")
    path = tmp_path / "controlDict"
    path.write_text(render_control_dict(cfg, run), encoding="utf-8")
    fields = field_average_block(path)["fields"]
    for name in ("momentum", "fD", "q"):
        assert fields[name]["prime2Mean"] is False


# --------------------------------------------------------------------------- #
# generation, metadata, idempotency
# --------------------------------------------------------------------------- #

@pytest.fixture
def template(tmp_path):
    """A study directory with a baseCase and a two-axis study.yaml."""
    write_study(tmp_path)
    write_case(tmp_path / "baseCase")
    (tmp_path / "baseCase" / "Allmesh").write_text("#!/bin/bash\n",
                                                   encoding="utf-8")
    return tmp_path


def test_generation_makes_one_directory_per_pair(template):
    study = load_study(template)
    paths = clone_cases(template, study, study.expanded_cases())
    assert [p.name for p in paths] == [
        "Kn100_np1x", "Kn100_np2x", "Kn100_np5x", "Kn100_np10x",
        "Kn0p1_np1x", "Kn0p1_np2x", "Kn0p1_np5x", "Kn0p1_np10x"]


def test_a_generated_case_states_its_own_level(template):
    """§14.14. Both axes recoverable from the file, not only from the name."""
    study = load_study(template)
    clone_cases(template, study, study.expanded_cases())
    for case in study.expanded_cases():
        destination = template / study.case_path(case)
        apply_case_parameters(destination, case, study)
        data = yaml.safe_load(
            (destination / "case.yaml").read_text(encoding="utf-8"))
        assert data["meta"]["cai_case"] == case.knudsen_case.name
        assert data["meta"]["particle_level"] == case.particle_level.name
        assert data["meta"]["numerical_particle_multiplier"] == case.multiplier
        assert data["dsmc"]["numerical_particle_multiplier"] == case.multiplier


def test_a_generated_case_carries_the_study_wide_schedule(template):
    study = load_study(template)
    clone_cases(template, study, study.expanded_cases())
    for case in study.expanded_cases():
        destination = template / study.case_path(case)
        apply_case_parameters(destination, case, study)
        cfg = load_case_config(destination)
        assert cfg.dsmc.run_time_multiplier == 3.0
        assert cfg.dsmc.output_frequency_multiplier == 4.0


def test_generated_cases_differ_only_in_the_multiplier(template):
    """§15, and the experiment. Everything else in case.yaml is identical
    across a physical case's four variants."""
    study = load_study(template)
    clone_cases(template, study, study.expanded_cases())
    documents = {}
    for case in study.expanded_cases():
        if case.knudsen_case.name != "Kn100":
            continue
        destination = template / study.case_path(case)
        apply_case_parameters(destination, case, study)
        data = yaml.safe_load(
            (destination / "case.yaml").read_text(encoding="utf-8"))
        data.pop("meta")
        documents[case.name] = data

    reference = documents["Kn100_np1x"]
    for name, data in documents.items():
        differing = [k for k in data["dsmc"]
                     if data["dsmc"][k] != reference["dsmc"][k]]
        assert differing == ([] if name == "Kn100_np1x"
                             else ["numerical_particle_multiplier"])
        for section in ("exit", "gas", "geometry", "mesh", "nozzle",
                        "resolution", "checks", "output", "post"):
            assert data[section] == reference[section], (name, section)


def test_the_manifest_row_names_both_axes(template):
    """§12."""
    study = load_study(template)
    case = study.expanded_cases()[1]          # Kn100_np2x
    clone_cases(template, study, [case])
    destination = template / study.case_path(case)
    apply_case_parameters(destination, case, study)
    cfg, geom, exit_state, plan, run = derive_generated(destination)

    entry = manifest_entry(case, study, cfg, geom, exit_state, plan, run)
    assert entry["name"] == "Kn100_np2x"
    assert entry["cai_case"] == "Kn100"
    assert entry["particle_level"] == "np2x"
    assert entry["numerical_particle_multiplier"] == 2.0
    assert entry["template"] == "baseCase"
    assert entry["path"] == "Cases/Kn100_np2x"
    assert entry["write_control"] == "timeStep"
    assert entry["run_time_multiplier"] == 3.0
    assert entry["output_frequency_multiplier"] == 4.0
    assert entry["start_time_s"] == 0.0
    assert entry["mesh_id"]
    assert entry["averaging"]["fields"]["dsmcRhoN"]["prime2Mean"] is True
    assert entry["n_equivalent_particles"] == pytest.approx(
        entry["baseline_n_equivalent_particles"] / 2.0)


def test_the_four_variants_share_a_mesh_id(template):
    """The manifest alone can answer 'are these on the same mesh?'."""
    study = load_study(template)
    ids = set()
    for case in study.expanded_cases():
        if case.knudsen_case.name != "Kn100":
            continue
        clone_cases(template, study, [case])
        destination = template / study.case_path(case)
        apply_case_parameters(destination, case, study)
        cfg, geom, exit_state, plan, run = derive_generated(destination)
        ids.add(manifest_entry(case, study, cfg, geom, exit_state, plan,
                               run)["mesh_id"])
    assert len(ids) == 1


def test_regenerating_is_idempotent(template):
    """§14.15. Generation is a pure function of study.yaml and the template, so
    running it twice gives byte-identical cases."""
    study = load_study(template)
    first = {}
    for _ in range(2):
        clone_cases(template, study, study.expanded_cases())
        for case in study.expanded_cases():
            destination = template / study.case_path(case)
            apply_case_parameters(destination, case, study)
            text = (destination / "case.yaml").read_text(encoding="utf-8")
            if case.name in first:
                assert text == first[case.name]
            else:
                first[case.name] = text


def test_a_variant_with_results_is_detected_as_complete(template):
    """The restart guard is per generated case, so one finished variant does
    not make its siblings look finished."""
    study = load_study(template)
    clone_cases(template, study, study.expanded_cases())
    finished = template / study.case_path(study.expanded_cases()[0])
    (finished / "0.001").mkdir()
    assert case_is_complete(finished)
    assert not case_is_complete(
        template / study.case_path(study.expanded_cases()[1]))


# --------------------------------------------------------------------------- #
# the validator
# --------------------------------------------------------------------------- #

def _generate_tree(root: Path) -> dict:
    """Generate the two-axis study under ``root`` and write its dictionaries.

    Returns the manifest document, in memory -- ``verify_tree`` takes it that
    way, so the test does not have to round-trip through YAML.
    """
    from plumetools.cai2012.dictionaries import write_case_dictionaries

    study = load_study(root)
    clone_cases(root, study, study.expanded_cases())
    entries = []
    for case in study.expanded_cases():
        destination = root / study.case_path(case)
        apply_case_parameters(destination, case, study)
        cfg, geom, exit_state, plan, run = derive_generated(destination)
        write_case_dictionaries(destination, cfg, exit_state, run)
        # blockMeshDict is written by the mesh module in a real run; a stand-in
        # is enough here, and it has to be identical across a group.
        (destination / "system" / "blockMeshDict").write_text(
            f"cells {plan.n_cells}\n", encoding="utf-8")
        entries.append(manifest_entry(case, study, cfg, geom, exit_state, plan,
                                      run))
    return {"cases": entries}


def test_the_validator_passes_a_correctly_generated_tree(template):
    """§13. And it is not vacuous: the negative cases below fail it."""
    manifest = _generate_tree(template)
    report = verify_tree(template, manifest)
    assert report.passed, report.lines()
    assert report.passes


def test_the_validator_catches_a_mis_scaled_weight(template):
    """The failure this whole check exists for: a variant that is not the
    multiple of the baseline it claims to be."""
    manifest = _generate_tree(template)
    properties = template / "Cases" / "Kn100_np5x" / "constant" / "dsmcProperties"
    text = properties.read_text(encoding="utf-8")
    original = foam_float(properties, "nEquivalentParticles")
    properties.write_text(
        text.replace(f"{original:.6e}", f"{original * 1.5:.6e}"),
        encoding="utf-8")

    report = verify_tree(template, manifest)
    assert not report.passed
    assert any("nEquivalentParticles" in f and "Kn100_np5x" in f
               for f in report.failures)


def test_the_validator_catches_a_different_time_step(template):
    """A sweep in which one variant also moved deltaT measures nothing."""
    manifest = _generate_tree(template)
    control = template / "Cases" / "Kn100_np2x" / "system" / "controlDict"
    text = control.read_text(encoding="utf-8")
    delta_t = foam_float(control, "deltaT")
    control.write_text(
        text.replace(f"deltaT          {delta_t:g};",
                     f"deltaT          {delta_t * 0.5:g};"),
        encoding="utf-8")

    report = verify_tree(template, manifest)
    assert not report.passed
    assert any("deltaT" in f for f in report.failures)


def test_the_validator_catches_a_different_mesh(template):
    manifest = _generate_tree(template)
    block_mesh = template / "Cases" / "Kn100_np10x" / "system" / "blockMeshDict"
    block_mesh.write_text("cells 17\n", encoding="utf-8")

    report = verify_tree(template, manifest)
    assert not report.passed
    assert any("blockMeshDict" in f for f in report.failures)


def test_the_validator_catches_prime2mean_being_switched_off(template):
    """Losing the variance is silent at run time and unrecoverable afterwards."""
    manifest = _generate_tree(template)
    control = template / "Cases" / "Kn100_np1x" / "system" / "controlDict"
    text = control.read_text(encoding="utf-8")
    control.write_text(text.replace("prime2Mean  on;", "prime2Mean  off;"),
                       encoding="utf-8")

    report = verify_tree(template, manifest)
    assert not report.passed
    assert any("prime2Mean" in f for f in report.failures)


def test_the_validator_reports_a_case_with_no_dictionaries(template):
    """Rather than passing a tree ./Allmesh has not been run on."""
    manifest = _generate_tree(template)
    (template / "Cases" / "Kn100_np1x" / "constant" / "dsmcProperties").unlink()

    report = verify_tree(template, manifest)
    assert not report.passed
    assert any("AllmeshCases" in f for f in report.failures)


def test_the_validator_catches_a_duplicate_case_name(template):
    manifest = _generate_tree(template)
    manifest["cases"].append(dict(manifest["cases"][0]))

    report = verify_tree(template, manifest)
    assert not report.passed
    assert any("duplicate" in f for f in report.failures)


def test_a_level_knows_whether_it_is_the_baseline():
    assert ParticleLevel(name="np1x", multiplier=1.0).is_baseline
    assert not ParticleLevel(name="np2x", multiplier=2.0).is_baseline
