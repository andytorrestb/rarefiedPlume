"""The statistical-health sweep: the matrix, the cost model, and the audit.

Nothing here runs a solver, meshes anything, or needs CaseFoam. What is covered
is the part that decides *what the nine cases are* and *what their output
means*, which is where a wrong answer is invisible: a matrix whose cases differ
in something other than the two swept axes still generates, still runs, still
produces nine pictures, and the pictures still look like plumes.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

from _cai2012 import derived, write_case
from plumetools.cai2012 import audit, health
from plumetools.cai2012.post import PostError, SampledField
from plumetools.cai2012.study import StudyError

REPO = Path(__file__).resolve().parents[2]

#: The shipped study, which the tests below check the shape of as well as the
#: loader. A matrix that stopped being 3 x 3 would otherwise only show up in a
#: 12-hour run.
STUDY = REPO / "cases" / "cai2012-health" / "study.yaml"

MINIMAL_STUDY = {
    "weights": [
        {"name": "ppc005", "target_particles_per_cell": 5.0},
        {"name": "ppc020", "target_particles_per_cell": 20.0},
    ],
    "samplings": [
        {"name": "s0p5", "sampling_domain_transits": 0.5},
        {"name": "s1p5", "sampling_domain_transits": 1.5},
    ],
}


def write_study(tmp_path, document=None, **changes):
    path = tmp_path / "study.yaml"
    document = dict(document or MINIMAL_STUDY, **changes)
    path.write_text(yaml.safe_dump(document), encoding="utf-8")
    return path


# --------------------------------------------------------------------------- #
# the matrix
# --------------------------------------------------------------------------- #

def test_the_matrix_is_the_cross_product_of_the_two_axes(tmp_path):
    study = health.load_study(write_study(tmp_path))
    assert [case.name for case in study.enabled_cases()] == [
        "ppc005/s0p5", "ppc005/s1p5", "ppc020/s0p5", "ppc020/s1p5"]


def test_cases_come_out_cheapest_first_whatever_order_the_file_is_in(tmp_path):
    """Reordering study.yaml must not change the tree, and the top weight row is
    most of the bill -- so it is generated and run last, not first."""
    reversed_study = {
        "weights": list(reversed(MINIMAL_STUDY["weights"])),
        "samplings": list(reversed(MINIMAL_STUDY["samplings"])),
    }
    study = health.load_study(write_study(tmp_path, reversed_study))
    assert [case.name for case in study.enabled_cases()] == [
        "ppc005/s0p5", "ppc005/s1p5", "ppc020/s0p5", "ppc020/s1p5"]


def test_the_case_path_is_a_tree_not_a_flat_name(tmp_path):
    """Cases/<weight>/<sampling>. A flat name would work equally well as a
    directory and would lose the hierarchy CaseFoam is here to build."""
    study = health.load_study(write_study(tmp_path))
    case = study.enabled_cases()[-1]
    assert study.case_path(case) == Path("Cases") / "ppc020" / "s1p5"


def test_a_disabled_level_removes_a_whole_row(tmp_path):
    document = {
        "weights": [dict(MINIMAL_STUDY["weights"][0], enabled=False),
                    MINIMAL_STUDY["weights"][1]],
        "samplings": MINIMAL_STUDY["samplings"],
    }
    study = health.load_study(write_study(tmp_path, document))
    assert [c.name for c in study.enabled_cases()] == ["ppc020/s0p5",
                                                       "ppc020/s1p5"]
    # ...and the disabled cells are still on record as intended-but-skipped
    assert len(study.all_cases()) == 4


def test_disabling_a_whole_axis_is_rejected(tmp_path):
    document = {
        "weights": [dict(w, enabled=False) for w in MINIMAL_STUDY["weights"]],
        "samplings": MINIMAL_STUDY["samplings"],
    }
    with pytest.raises(StudyError, match="not a sweep"):
        health.load_study(write_study(tmp_path, document))


def test_there_is_no_override_mechanism_and_asking_for_one_says_why(tmp_path):
    """Everything except the two axes must be identical across the matrix, or
    the nine images are not comparable. A per-case override block is the way
    that would quietly stop being true, so it is refused by name."""
    with pytest.raises(StudyError, match="IDENTICAL across the matrix"):
        health.load_study(write_study(tmp_path, overrides={"dsmc": {}}))


def test_two_levels_at_the_same_value_are_rejected(tmp_path):
    """They would differ in directory name only, and the matrix would carry a
    column that measures nothing."""
    document = dict(MINIMAL_STUDY)
    document["samplings"] = [
        {"name": "s1p5", "sampling_domain_transits": 1.5},
        {"name": "also1p5", "sampling_domain_transits": 1.5},
    ]
    with pytest.raises(StudyError, match="differ in name only"):
        health.load_study(write_study(tmp_path, document))


def test_a_non_positive_sampling_duration_is_rejected(tmp_path):
    """endTime would land at or before average_start, so fieldAverage writes no
    *Mean field and there is nothing to measure -- several hours later."""
    document = dict(MINIMAL_STUDY)
    document["samplings"] = [{"name": "s0", "sampling_domain_transits": 0.0}]
    with pytest.raises(StudyError, match="nothing to measure"):
        health.load_study(write_study(tmp_path, document))


def test_the_baseline_cell_is_the_one_that_matches_cases_cai2012(tmp_path):
    """20 parcels/cell and 1.5 transits is Kn100 re-run. The pinned colour
    ranges and the control on the published numbers both come from it, so
    losing it silently would leave the sweep with no anchor."""
    study = health.load_study(write_study(tmp_path))
    assert study.baseline is not None
    assert study.baseline.name == "ppc020/s1p5"


# --------------------------------------------------------------------------- #
# writing a case
# --------------------------------------------------------------------------- #

def _case(weight=20.0, transits=1.5, weight_name="ppc020", sampling_name="s1p5"):
    return health.HealthCase(
        weight=health.WeightLevel(name=weight_name, particles_per_cell=weight),
        sampling=health.SamplingLevel(name=sampling_name,
                                      domain_transits=transits))


def test_exactly_two_values_are_written_into_a_case(tmp_path):
    """The whole comparability claim is that nothing else varies. If generation
    wrote the mesh or the time step per case, nine runs on nine meshes would
    still produce nine plausible pictures."""
    write_case(tmp_path)
    before = yaml.safe_load((tmp_path / "case.yaml").read_text(encoding="utf-8"))

    study = health.load_study(write_study(tmp_path))
    health.apply_case_parameters(tmp_path, _case(weight=40.0, transits=4.5),
                                 study)
    after = yaml.safe_load((tmp_path / "case.yaml").read_text(encoding="utf-8"))

    assert after["resolution"]["target_particles_per_cell"] == 40.0
    assert after["dsmc"]["sampling_domain_transits"] == 4.5
    for section in ("nozzle", "exit", "gas", "geometry", "mesh", "checks",
                    "output", "post"):
        assert after[section] == before[section], (
            f"{section} changed; only the two swept keys may")
    unchanged = {k: v for k, v in before["dsmc"].items()
                 if k != "sampling_domain_transits"}
    assert {k: v for k, v in after["dsmc"].items()
            if k != "sampling_domain_transits"} == unchanged


def test_a_template_at_another_knudsen_number_is_refused(tmp_path):
    """The sweep exists to explain the Kn = 100 imagery. At another Kn the
    density, the parcel count and the picture all change, and nothing in the
    matrix could be compared with cases/cai2012 at all."""
    write_case(tmp_path, exit={"knudsen": 0.1})
    study = health.load_study(write_study(tmp_path))
    with pytest.raises(StudyError, match="ONE Knudsen number"):
        health.apply_case_parameters(tmp_path, _case(), study)


def test_a_case_yaml_from_another_family_is_refused(tmp_path):
    (tmp_path / "case.yaml").write_text(
        yaml.safe_dump({"model": "markelov1999_plate", "stagnation": {}}),
        encoding="utf-8")
    study = health.load_study(write_study(tmp_path))
    with pytest.raises(StudyError, match="cai2012 template"):
        health.apply_case_parameters(tmp_path, _case(), study)


def test_the_generated_case_still_loads_with_the_cai2012_loader(tmp_path):
    """The whole reuse argument is that a health case IS a cai2012 case: the
    same Allmesh, Allrun, checks and post-processing, with two values changed."""
    from plumetools.cai2012.config import load_case_config

    write_case(tmp_path)
    study = health.load_study(write_study(tmp_path))
    health.apply_case_parameters(tmp_path, _case(weight=5.0, transits=0.5),
                                 study)
    cfg = load_case_config(tmp_path)
    assert cfg.resolution.target_particles_per_cell == 5.0
    assert cfg.dsmc.sampling_domain_transits == 0.5
    assert cfg.exit.knudsen == health.FAMILY_KNUDSEN


# --------------------------------------------------------------------------- #
# cost
# --------------------------------------------------------------------------- #

def test_the_prior_reproduces_the_run_it_was_calibrated_on():
    """cases/cai2012's Kn100: 1.2e6 parcels, 1.31e6 cells, 7212 steps, ~1 h."""
    hours = health.PRIOR_COST_MODEL.seconds(
        parcels=1.2e6, cells=1.31e6, steps=7212) / 3600.0
    assert hours == pytest.approx(1.0, rel=1e-9)


def test_the_two_terms_are_fitted_separately_from_two_runs():
    """A one-term model calibrated on the cheap row understates the expensive
    one, which is the direction that costs twelve hours to find out about."""
    truth = health.CostModel(per_parcel_s=2.0e-7, per_cell_s=5.0e-8)
    measurements = [
        {"case": "a", "parcels": 3.0e5, "cells": 1.3e6, "steps": 5000},
        {"case": "b", "parcels": 2.4e6, "cells": 1.3e6, "steps": 9000},
        {"case": "c", "parcels": 1.2e6, "cells": 3.7e6, "steps": 7000},
    ]
    for entry in measurements:
        entry["seconds"] = truth.seconds(**{k: entry[k] for k in
                                            ("parcels", "cells", "steps")})

    fitted = health.fit_cost_model(measurements)
    assert fitted.per_parcel_s == pytest.approx(truth.per_parcel_s, rel=1e-9)
    assert fitted.per_cell_s == pytest.approx(truth.per_cell_s, rel=1e-9)
    assert fitted.n_measurements == 3


def test_one_measurement_scales_the_prior_and_says_so():
    """Two coefficients cannot be fitted to one number. Presenting the result
    as a two-term fit would make an assumed split look measured."""
    model = health.fit_cost_model([
        {"case": "a", "parcels": 1.2e6, "cells": 1.31e6, "steps": 7212,
         "seconds": 7200.0}])
    assert model.n_measurements == 1
    assert "not separable" in model.source
    assert model.seconds(parcels=1.2e6, cells=1.31e6,
                         steps=7212) == pytest.approx(7200.0)


def test_collinear_measurements_fall_back_rather_than_fitting_a_negative_cost():
    """Two runs at the same weight and mesh do not separate the terms, and the
    least-squares answer can come out saying a bigger mesh runs faster."""
    measurements = [
        {"case": "a", "parcels": 1.0e6, "cells": 1.3e6, "steps": 1000,
         "seconds": 400.0},
        {"case": "b", "parcels": 1.0e6, "cells": 1.3e6, "steps": 2000,
         "seconds": 800.0},
    ]
    model = health.fit_cost_model(measurements)
    assert model.per_parcel_s > 0.0 and model.per_cell_s > 0.0


def test_no_measurements_gives_the_prior_not_a_crash():
    assert health.fit_cost_model([]) is health.PRIOR_COST_MODEL
    assert health.load_cost_model(Path("nowhere.yaml")) is health.PRIOR_COST_MODEL


def test_re_running_a_case_replaces_its_measurement(tmp_path):
    """A resumed leg's wall clock is not the cost of the whole case; keeping
    both and averaging would quietly halve every estimate."""
    path = tmp_path / "cost-model.yaml"
    entry = {"case": "ppc005/s0p5", "parcels": 3e5, "cells": 1.3e6,
             "steps": 5000, "seconds": 600.0}
    health.record_measurement(path, entry)
    health.record_measurement(path, dict(entry, seconds=900.0))

    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert len(document["measurements"]) == 1
    assert document["measurements"][0]["seconds"] == 900.0


def test_the_cost_report_marks_an_unmeasured_estimate_as_a_prior():
    costs = [health.estimate_case(_case(), parcels=1.2e6, cells=1310720,
                                  steps=7212, writes=28,
                                  bytes_per_write=2.7e8,
                                  model=health.PRIOR_COST_MODEL)]
    text = "\n".join(health.cost_report(costs, health.PRIOR_COST_MODEL))
    assert "PRIOR" in text
    assert "TOTAL" in text


def test_sampled_writes_counts_only_the_frames_that_carry_a_mean_field(tmp_path):
    """Frames written during the transient have no *Mean field at all. Counting
    them as sampled frames would overstate every convergence curve's length and
    quietly include instantaneous images in a series of averages."""
    _, _, _, _, run = derived(
        tmp_path / "case",
        dsmc={"transient_basis": "transits", "transient_domain_transits": 2.0,
              "sampling_domain_transits": 1.5, "write_interval_s": None})
    total = run.n_writes
    sampled = health.sampled_writes(run)
    assert 0 < sampled <= total
    first_sampled = run.write_interval_s * (total - sampled + 1)
    assert first_sampled >= run.average_start_s - 1e-15


# --------------------------------------------------------------------------- #
# the manifest
# --------------------------------------------------------------------------- #

def test_the_manifest_records_what_is_held_fixed(tmp_path):
    entries = [
        {"name": "ppc005/s0p5", "Kn": 100.0, "n_cells": 1310720,
         "deltaT_s": 1.373e-6},
        {"name": "ppc020/s0p5", "Kn": 100.0, "n_cells": 1310720,
         "deltaT_s": 1.373e-6},
    ]
    study = health.load_study(write_study(tmp_path))
    path = health.write_manifest(tmp_path / "manifest.yaml", study, entries)
    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert document["held_fixed"]["n_cells"] == 1310720
    assert document["held_fixed"]["deltaT_s"] == 1.373e-6


def test_a_matrix_that_stopped_being_comparable_says_so_in_the_manifest(tmp_path):
    """Two cases on two meshes still generate, still run, and still draw nine
    plausible pictures. The manifest is where that becomes visible."""
    entries = [
        {"name": "ppc005/s0p5", "Kn": 100.0, "n_cells": 1310720},
        {"name": "ppc020/s0p5", "Kn": 100.0, "n_cells": 3726000},
    ]
    study = health.load_study(write_study(tmp_path))
    path = health.write_manifest(tmp_path / "manifest.yaml", study, entries)
    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert "NOT SHARED" in str(document["held_fixed"]["n_cells"])
    assert "no longer comparable" in str(document["held_fixed"]["n_cells"])


# --------------------------------------------------------------------------- #
# the shipped study
# --------------------------------------------------------------------------- #

@pytest.mark.skipif(not STUDY.is_file(), reason="the family is not generated yet")
def test_the_shipped_study_is_a_three_by_three_matrix():
    study = health.load_study(STUDY)
    assert len(study.enabled_weights()) == 3
    assert len(study.enabled_samplings()) == 3
    assert len(study.enabled_cases()) == 9


@pytest.mark.skipif(not STUDY.is_file(), reason="the family is not generated yet")
def test_the_shipped_axes_bracket_the_baseline_and_the_floor():
    """5 is checks.min_particles_per_cell, 20 is what cases/cai2012 runs at.
    A sweep that did not span them could not say anything about either."""
    study = health.load_study(STUDY)
    weights = [w.particles_per_cell for w in study.enabled_weights()]
    transits = [s.domain_transits for s in study.enabled_samplings()]
    assert min(weights) <= 5.0 and 20.0 in weights and max(weights) > 20.0
    assert min(transits) < 1.5 and 1.5 in transits and max(transits) > 1.5
    assert study.baseline is not None


@pytest.mark.skipif(not STUDY.is_file(), reason="the family is not generated yet")
def test_the_shipped_template_is_at_the_family_knudsen_number():
    from plumetools.cai2012.config import load_case_config

    cfg = load_case_config(STUDY.parent / "baseCase")
    assert cfg.exit.knudsen == health.FAMILY_KNUDSEN
    assert cfg.dsmc.write_interval_s is not None, (
        "the health family pins its write interval: the frame times have to "
        "line up across cases for the overlap check and the contact sheet")


@pytest.mark.skipif(not STUDY.is_file(), reason="the family is not generated yet")
def test_the_write_interval_gives_the_shortest_case_a_usable_curve():
    """Two frames are not a convergence curve, and 100 frames is a study whose
    imagery step outlasts its solve. The interval has to serve both ends of a
    9x sampling axis."""
    from plumetools.cai2012 import inflow, mesh
    from plumetools.cai2012.config import load_case_config
    from plumetools.cai2012.geometry import from_config

    study = health.load_study(STUDY)
    cfg = load_case_config(STUDY.parent / "baseCase")
    geom = from_config(cfg)
    exit_state = inflow.from_config(cfg)
    plan = mesh.plan(cfg, geom, exit_state)

    counts = {}
    for level in study.enabled_samplings():
        run = inflow.derive_run_settings(
            _sampling(cfg, level.domain_transits), exit_state, geom,
            min_cell_size_m=plan.min_cell_size_m,
            exit_cell_volume_m3=plan.exit_cell_volume_m3)
        counts[level.name] = health.sampled_writes(run)

    assert min(counts.values()) >= 3, (
        f"the shortest case gets {min(counts.values())} sampled frame(s): "
        f"{counts}. Fewer than three is not a curve.")
    assert max(counts.values()) <= 60, (
        f"the longest case gets {max(counts.values())} sampled frames: "
        f"{counts}. Every frame is ~255 MB and ~20 s of rendering per field.")


def _sampling(cfg, transits):
    import dataclasses

    return dataclasses.replace(
        cfg, dsmc=dataclasses.replace(cfg.dsmc,
                                      sampling_domain_transits=transits))


@pytest.mark.skipif(not STUDY.is_file(), reason="the family is not generated yet")
def test_the_viz_window_matches_the_template_transient():
    """viz.yaml's time_min is a literal copy of dsmc.average_start_s.

    It has to be: fieldAverage writes no *Mean field before then, and with
    prefer_mean on those frames fall back to the instantaneous field and the
    series silently changes quantity part way through. If the template's
    transient is ever changed and this is not, the study starts drawing shot
    noise labelled as a running average -- and both look like a plume.
    """
    from plumetools.cai2012 import inflow, mesh
    from plumetools.cai2012.config import load_case_config
    from plumetools.cai2012.geometry import from_config
    from plumetools.viz import load_spec

    cfg = load_case_config(STUDY.parent / "baseCase")
    geom = from_config(cfg)
    exit_state = inflow.from_config(cfg)
    plan = mesh.plan(cfg, geom, exit_state)
    run = inflow.derive_run_settings(
        cfg, exit_state, geom, min_cell_size_m=plan.min_cell_size_m,
        exit_cell_volume_m3=plan.exit_cell_volume_m3)

    spec = load_spec(STUDY.parent / "viz.yaml")
    assert spec.sampling.time_min is not None, (
        "prefer_mean is on over a series here; without a window the transient "
        "frames fall back to the instantaneous field")
    assert spec.sampling.time_min == pytest.approx(run.average_start_s,
                                                   rel=1e-3), (
        f"viz.yaml time_min is {spec.sampling.time_min:g} but the template "
        f"starts averaging at {run.average_start_s:g}")


@pytest.mark.skipif(not STUDY.is_file(), reason="the family is not generated yet")
def test_every_drawn_field_has_a_pinned_range():
    """Ranges are auto-pinned across the FRAMES of one series but not across
    cases. Without an explicit range the nine cases get nine colour scales, the
    contact sheet compares nothing, and it looks entirely fine -- the cells
    differ visibly, and the differences are the colour maps."""
    from plumetools.viz import load_spec

    spec = load_spec(STUDY.parent / "viz.yaml")
    unpinned = [field.name for field in spec.fields if field.range is None]
    assert not unpinned, (
        f"{unpinned} would autoscale per case, so no two cells of the contact "
        f"sheet would be on the same scale")


@pytest.mark.skipif(not STUDY.is_file(), reason="the family is not generated yet")
def test_the_health_spec_averages_where_the_other_studies_do_not():
    """The one setting that separates this study's imagery from cases/cai2012's.
    Flipped by accident, the series becomes shot noise and the convergence
    curve measures nothing -- while still producing plausible plume pictures."""
    from plumetools.viz import load_spec

    health_spec = load_spec(STUDY.parent / "viz.yaml")
    cai_spec = load_spec(REPO / "cases" / "cai2012" / "viz.yaml")
    assert health_spec.sampling.prefer_mean is True
    assert cai_spec.sampling.prefer_mean is False
    assert "dsmcRhoN" in [f.name for f in health_spec.fields]


@pytest.mark.skipif(not STUDY.is_file(), reason="the family is not generated yet")
def test_the_shared_runner_scripts_have_not_drifted_from_cases_cai2012():
    """A health case IS a cai2012 case; the scripts are copies rather than a
    fork. A fix applied to one and not the other is the failure mode, and it
    would show up as two families that quietly run differently."""
    shared = ("Allclean", "Allmesh", "Allpost", "Allrun", "postProcess.py",
              "runInflow.py")
    origin = REPO / "cases" / "cai2012" / "baseCase"
    copy = STUDY.parent / "baseCase"
    drifted = [name for name in shared
               if (origin / name).read_bytes() != (copy / name).read_bytes()]
    assert not drifted, (
        f"{drifted} differ between cases/cai2012/baseCase and "
        f"cases/cai2012-health/baseCase. Re-copy, or if the change is "
        f"deliberate, say so here.")


# --------------------------------------------------------------------------- #
# occupancy
# --------------------------------------------------------------------------- #

def _sampled(n_cells=8, volume=1.0e-6, density=1.0e16):
    centres = np.zeros((n_cells, 3))
    centres[:, 0] = np.linspace(0.005, 1.9, n_cells)
    return SampledField(
        time="0.01",
        centres=centres,
        volumes=np.full(n_cells, volume),
        number_density=np.full(n_cells, density),
        velocity=np.zeros((n_cells, 3)),
        temperature=np.full(n_cells, 300.0),
    )


def test_occupancy_is_dsmcRhoN_itself_with_no_volume_factor(tmp_path):
    """dsmcRhoN IS the parcel count in the cell -- DSMCCloud::calculateFields
    adds 1 per parcel and never divides by volume.

    Multiplying by V, which the repository's own catalog and specs used to say
    to do, is out by 1/V -- a factor of 1e6 on the Cai mesh. On a logarithmic
    colour scale that produces an entirely plausible picture, and it would
    scale every occupancy figure in this study by the same constant.
    """
    cfg, geom, exit_state, plan, _ = derived(tmp_path / "case")
    sampled = _sampled(volume=1.0e-6)
    parcel_count = np.full(sampled.n_cells, 20.0)

    result = audit.occupancy_audit(
        sampled, cfg, geom, exit_state,
        n_equivalent_particles=1.0e6,
        parcel_count=parcel_count,
        averaged_steps=1000.0,
        core_cell_size_m=plan.core_cell_size_m)

    domain = [r for r in result["regions"] if r["name"] == "domain"][0]
    assert domain["median"] == pytest.approx(20.0)
    assert domain["fraction_below_floor"] == 0.0


def test_the_median_is_the_headline_not_the_mean(tmp_path):
    """Occupancy spans orders of magnitude between the exit and the plume edge.
    A mean over that reports the exit cell as though it were typical."""
    result = audit.region_occupancy(
        "r", "d", [1.0, 1.0, 1.0, 1.0, 1.0, 10000.0],
        floor=5.0, target=20.0, averaged_steps=1.0)
    assert result.median == pytest.approx(1.0)
    assert result.mean > 1000.0


def test_an_empty_region_is_zeros_not_an_exception():
    """A mesh can legitimately have no cell in a band, and a study that stopped
    on it would be reporting a mesh property as a failure."""
    result = audit.region_occupancy("r", "d", [], floor=5.0, target=20.0,
                                    averaged_steps=1.0)
    assert result.n_cells == 0 and result.median == 0.0


def test_a_parcel_field_from_another_mesh_is_refused(tmp_path):
    cfg, geom, exit_state, plan, _ = derived(tmp_path / "case")
    sampled = _sampled(n_cells=8)
    with pytest.raises(PostError, match="not the same mesh"):
        audit.occupancy_audit(sampled, cfg, geom, exit_state,
                              n_equivalent_particles=1.0e6,
                              parcel_count=np.ones(5),
                              averaged_steps=1000.0)


def test_the_particle_weight_identity_is_checked_not_assumed():
    """Every occupancy number rests on rhoN*V/dsmcRhoN being one global weight.

    Two things can break it silently: dsmcRhoN not being a count (which is how
    the repository had it), or a solver with radial weighting. Either way the
    ratio stops being constant, and every figure in the audit would be wrong
    while still looking entirely plausible.
    """
    volumes = [1.0e-6, 1.0e-6]
    # weight = rhoN * V / count, so 20 parcels in a 1e-6 m^3 cell at a 1e6
    # weight is rhoN = 2e13.
    good = audit.weight_consistency([2.0e13, 4.0e13], [20.0, 40.0], volumes,
                                    1.0e6)
    assert good["consistent"] and good["measured_weight"] == pytest.approx(1e6)

    bad = audit.weight_consistency([2.0e13, 4.0e13], [20.0, 80.0], volumes,
                                   1.0e6)
    assert not bad["consistent"]


def test_reading_dsmcRhoN_as_a_density_fails_the_weight_check():
    """The specific error this check exists to catch, out by exactly 1/V."""
    result = audit.weight_consistency([2.0e13], [20.0 / 1.0e-6], [1.0e-6],
                                      1.0e6)
    assert not result["consistent"]
    # out by exactly V
    assert result["measured_weight"] == pytest.approx(1.0e6 * 1.0e-6)


def test_cells_with_no_parcels_do_not_break_the_weight_check():
    """0/0 says nothing about the weight; the far field is full of such cells."""
    result = audit.weight_consistency([1.0e13, 0.0], [10.0, 0.0], [1.0e-6] * 2,
                                      1.0e6)
    assert result["consistent"] and result["n_compared"] == 1


def test_the_statistical_budget_is_the_product_of_the_two_axes():
    """The sweep's central claim is that the axes trade off. Both directions
    have to give the same budget, or the comparison tests nothing."""
    assert audit.statistical_budget(20.0, 1.5, 2060.0) == pytest.approx(
        audit.statistical_budget(5.0, 6.0, 2060.0))


# --------------------------------------------------------------------------- #
# the overlap
# --------------------------------------------------------------------------- #

def _write_field(path: Path, values):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = "\n".join(f"{v:.10g}" for v in values)
    path.write_text(
        "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
        "    class volScalarField;\n    object rhoNMean;\n}\n\n"
        "dimensions      [0 0 0 0 0 0 0];\n\n"
        f"internalField   nonuniform List<scalar>\n{len(values)}\n(\n{rows}\n)\n;\n\n"
        "boundaryField\n{\n}\n", encoding="utf-8", newline="\n")
    return path


def test_the_same_run_truncated_agrees_byte_for_byte(tmp_path):
    """Two cases differing only in endTime are the same trajectory. If this
    ever fails, the sweep's differences are not all statistical -- which
    outranks every other finding in the study."""
    values = [1.0e16, 2.0e16, 3.0e16]
    for case in ("long", "short"):
        _write_field(tmp_path / case / "0.006" / "rhoNMean", values)
        _write_field(tmp_path / case / "0.006" / "dsmcRhoNMean", values)

    result = audit.compare_overlap(tmp_path / "long", tmp_path / "short",
                                   fields=("rhoNMean", "dsmcRhoNMean"))
    assert result["all_identical"]
    assert result["n_compared"] == 2
    assert "IDENTICAL" in "\n".join(audit.overlap_report(result))


def test_a_disagreement_is_quantified_rather_than_just_flagged(tmp_path):
    """"They differ" is not actionable. One part in 1e-12 is a formatting
    artefact; one part in 1e-2 is a different run."""
    _write_field(tmp_path / "long" / "0.006" / "rhoNMean", [1.0e16, 2.0e16])
    _write_field(tmp_path / "short" / "0.006" / "rhoNMean", [1.0e16, 2.02e16])

    result = audit.compare_overlap(tmp_path / "long", tmp_path / "short",
                                   fields=("rhoNMean",))
    assert not result["all_identical"]
    # Symmetric: the denominator is max(|a|, |b|), so neither run is treated as
    # the reference. 0.02/2.02, not 0.02/2.00.
    assert result["worst"]["max_rel_diff"] == pytest.approx(0.02 / 2.02, rel=1e-6)
    assert result["worst"]["max_abs_diff"] == pytest.approx(2.0e14, rel=1e-6)
    assert "not deterministic" in "\n".join(audit.overlap_report(result))


def test_cases_that_share_no_time_are_reported_as_nothing_compared(tmp_path):
    """An empty comparison passing as "all identical" would be the worst
    possible outcome: the check would be green having tested nothing."""
    _write_field(tmp_path / "long" / "0.006" / "rhoNMean", [1.0])
    _write_field(tmp_path / "short" / "0.008" / "rhoNMean", [1.0])

    result = audit.compare_overlap(tmp_path / "long", tmp_path / "short",
                                   fields=("rhoNMean",))
    assert result["n_compared"] == 0
    assert "NOTHING COMPARED" in "\n".join(audit.overlap_report(result))


def test_common_times_are_matched_by_name_not_by_float_tolerance(tmp_path):
    """OpenFOAM named the directories; two runs at the same step and interval
    produce the same names. A float comparison would invent a tolerance where
    there is an exact answer."""
    for time in ("0.0021211", "0.00424357"):
        _write_field(tmp_path / "long" / time / "rhoNMean", [1.0])
    _write_field(tmp_path / "short" / "0.0021211" / "rhoNMean", [1.0])

    assert audit.common_times(tmp_path / "long",
                              tmp_path / "short") == ["0.0021211"]


def test_a_time_the_caller_named_but_neither_case_wrote_is_an_error(tmp_path):
    _write_field(tmp_path / "long" / "0.006" / "rhoNMean", [1.0])
    _write_field(tmp_path / "short" / "0.006" / "rhoNMean", [1.0])
    with pytest.raises(PostError, match="overlap check was told"):
        audit.compare_overlap(tmp_path / "long", tmp_path / "short",
                              fields=("rhoNMean",), times=["0.009"])


# --------------------------------------------------------------------------- #
# the sweep table
# --------------------------------------------------------------------------- #

def test_a_sweep_row_survives_a_case_with_no_audit_yet(tmp_path):
    """The table is built while the matrix is part way through. A missing
    occupancy file leaves blanks, not an exception."""
    entry = {"name": "ppc005/s0p5", "target_particles_per_cell": 5.0,
             "sampling_domain_transits": 0.5,
             "exit_particles_per_cell_estimated": 5.0}
    row = audit.sweep_row(entry, {"centerline": {"density_mean_rel_error": 0.02}},
                          None, steps_per_transit=2060.0)
    assert row["density_mean_rel_error_percent"] == pytest.approx(2.0)
    assert row["plume_median_occupancy"] is None
    assert row["statistical_budget"] == pytest.approx(5.0 * 0.5 * 2060.0)


def test_the_sweep_table_writes_every_column_even_when_empty(tmp_path):
    path = audit.write_sweep_table(tmp_path / "sweep-table.csv", [
        {"case": "ppc005/s0p5", "particles_per_cell": 5.0}])
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0].split(",") == list(audit.SWEEP_COLUMNS)
    assert len(lines[1].split(",")) == len(audit.SWEEP_COLUMNS)
