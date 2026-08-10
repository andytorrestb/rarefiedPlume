"""Tier 2: the Cai 2012 case against the real OpenFOAM toolchain.

Marked ``needs_openfoam`` and deselected by default::

    pytest -m "not needs_openfoam"    # the fast suite
    pytest -m needs_openfoam          # these

What the Tier-0 tests cannot answer
-----------------------------------
Everything here is a question about the *solver*, not about this repository's
Python:

* does ``topoSet``'s ``cylinderToFace`` select the same faces the Python
  prediction says it will, and does the resulting patch area approach
  :math:`\\pi R_0^2`?
* does ``plumeFieldInflow`` inject on the nozzle and **only** the nozzle?
* do the vacuum boundaries actually delete molecules?
* is the solver standard ``dsmcFoam``, with argon exactly as ``case.yaml``
  describes it?
* do the three Knudsen cases produce three different densities in the files the
  solver reads?
* does the measured injection rate match the theoretical Maxwellian flux?

The meshed case is built once per session and shared. The DSMC run is
deliberately short -- a few hundred steps, not a converged solution -- because
these are behavioural questions, not validation ones. The validation numbers come
from ``./AllpostCases`` on a real run.
"""

from __future__ import annotations

import math
import os
import re
import shutil
import subprocess

import numpy as np
import pytest

from conftest import (
    PLUME_LIB,
    REPO,
    SOLVER_COMMANDS,
    missing,
    plume_library_path,
    run_or_fail,
)
from plumetools.cai2012 import gas, inflow, mesh, post
from plumetools.cai2012.config import load_case_config
from plumetools.cai2012.geometry import from_config
from plumetools.mesh.boundary import read_boundary

pytestmark = pytest.mark.needs_openfoam

CAI_BASE_CASE = REPO / "cases" / "cai2012" / "baseCase"

#: A cheap mesh for the toolchain tests: ~30k cells rather than 1.3 million.
#: The physics is untouched -- only the resolution and the run length change --
#: and the questions these tests ask do not depend on either.
TEST_MESH = {
    "mesh": {
        "core_cell_size_m": 0.02,
        "outer_x_cells": 6,
        "outer_lateral_cells": 6,
        "max_cells": 4000000,
    },
    "geometry": {"x_max_over_D": 4.0, "y_half_over_D": 3.0, "z_half_over_D": 3.0},
    "post": {"centerline_x_over_D_max": 4.0, "plane_half_over_D": 2.0},
    "checks": {"min_nozzle_faces": 30, "max_nozzle_area_error": 0.08},
}


def _write_case(destination, knudsen: float = 100.0):
    """Copy baseCase, set the Knudsen number, and shrink the mesh."""
    import yaml

    shutil.copytree(CAI_BASE_CASE, destination)
    for script in destination.glob("All*"):
        script.chmod(0o755)
    for script in destination.glob("*.py"):
        script.chmod(0o755)

    path = destination / "case.yaml"
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    data["exit"]["knudsen"] = knudsen
    data["dsmc"]["transient_basis"] = "transits"
    data["dsmc"]["transient_domain_transits"] = 0.5
    data["dsmc"]["sampling_domain_transits"] = 0.5
    for section, values in TEST_MESH.items():
        data[section].update(values)
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return destination


def _env():
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(REPO), env.get("PYTHONPATH", "")]).rstrip(os.pathsep)
    return env


def _script(name, case, *args, timeout=7200):
    result = subprocess.run(
        ["bash", f"./{name}", *args], cwd=str(case), capture_output=True,
        text=True, timeout=timeout, env=_env(), check=False)
    if result.returncode != 0:
        tail = "\n".join((result.stdout + result.stderr).splitlines()[-40:])
        pytest.fail(f"./{name} failed:\n{tail}")
    return result


@pytest.fixture(scope="module")
def cai_meshed(openfoam, tmp_path_factory):
    """A meshed Kn = 100 case, built by its own ./Allmesh."""
    absent = missing(("blockMesh", "topoSet", "createPatch", "checkMesh"))
    if absent:
        pytest.skip(f"mesh utilities not on PATH: {absent}")
    case = _write_case(tmp_path_factory.mktemp("cai") / "Kn100")
    _script("Allmesh", case)
    return case


@pytest.fixture(scope="module")
def cai_solved(cai_meshed, tmp_path_factory):
    """A short DSMC run of the meshed case, post-processed.

    ``./Allpost`` is part of the fixture rather than of one test: it writes the
    cell centres and volumes that every field-reading test below needs, and
    running it here means the whole pipeline a user would run is exercised
    exactly once.
    """
    absent = missing(SOLVER_COMMANDS)
    if absent:
        pytest.skip(f"dsmcFoam utilities not on PATH: {absent}")
    if plume_library_path() is None:
        pytest.skip(
            f"{PLUME_LIB} is not built; run applications/dsmcBoundaryModels/Allwmake")

    case = tmp_path_factory.mktemp("cai-run") / "Kn100"
    shutil.copytree(cai_meshed, case)
    for script in list(case.glob("All*")) + list(case.glob("*.py")):
        script.chmod(0o755)

    _shorten(case, steps=300)
    _script("Allrun", case, "--serial")
    _script("Allpost", case)
    return case


def _shorten(case, *, steps: int, writes: int = 2):
    """Cut the run to ``steps`` timesteps, with ``writes`` evenly spaced.

    deltaT is untouched -- the Courant check depends on it -- so this only
    changes how many steps are taken. ``writeInterval`` is a STEP COUNT, because
    the generated controlDict uses ``writeControl timeStep``; ``steps`` is chosen
    to be a whole number of intervals so the last write lands exactly on endTime,
    which is the property the resume tests depend on.
    """
    assert steps % writes == 0, "steps must be a whole number of write intervals"
    end_time = float(_entry(case, "deltaT")) * steps
    for entry, value in (("endTime", f"{end_time:.12g}"),
                         ("writeInterval", str(steps // writes)),
                         ("functions/fieldAverage1/timeStart", "0")):
        run_or_fail(["foamDictionary", "system/controlDict",
                     "-entry", entry, "-set", value], case)


def _entry(case, key: str) -> str:
    result = subprocess.run(
        ["foamDictionary", "system/controlDict", "-entry", key, "-value"],
        cwd=str(case), capture_output=True, text=True, check=True)
    return result.stdout.strip()


# --------------------------------------------------------------------------- #
# the mesh and the nozzle patch
# --------------------------------------------------------------------------- #

def test_the_nozzle_patch_exists(cai_meshed):
    """It is created by topoSet + createPatch, not by blockMesh: a disk is not
    a block face."""
    patches = read_boundary(cai_meshed)
    assert "nozzle" in patches
    assert patches["nozzle"].n_faces > 0


def test_the_nozzle_patch_is_type_patch(cai_meshed):
    """plumeFieldInflow refuses to inject across a wallPolyPatch."""
    assert read_boundary(cai_meshed)["nozzle"].type == "patch"


def test_topo_set_selects_exactly_the_predicted_faces(cai_meshed):
    """The Python rule (face centre inside r <= R0) and topoSet's
    cylinderToFace are supposed to be the same rule. If they are not, the mesh
    that was checked is not the mesh that was built."""
    cfg = load_case_config(cai_meshed)
    geom = from_config(cfg)
    plan = mesh.plan(cfg, geom, inflow.from_config(cfg))
    predicted = mesh.predicted_nozzle_faces(plan, geom)
    assert read_boundary(cai_meshed)["nozzle"].n_faces == predicted.n_faces


def test_the_nozzle_area_approaches_pi_r_squared(cai_meshed):
    """A circle has no exact Cartesian representation, so this is a staircase.
    The area error is what the injected molecule flow scales with."""
    cfg = load_case_config(cai_meshed)
    geom = from_config(cfg)

    result = subprocess.run(
        ["checkMesh", "-writeAllFields", "-noTopology"], cwd=str(cai_meshed),
        capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout[-2000:]

    plan = mesh.plan(cfg, geom, inflow.from_config(cfg))
    area = mesh.predicted_nozzle_faces(plan, geom).area_m2
    exact = geom.exact_nozzle_area_m2
    assert abs(area - exact) / exact < 0.08


def test_the_nozzle_area_error_falls_as_the_mesh_is_refined(openfoam,
                                                            tmp_path_factory):
    """Convergence, not just closeness: the staircase must be a discretisation
    of the disk rather than a fixed shape that happens to be near it."""
    errors = []
    for cell in (0.03, 0.01):
        case = _write_case(
            tmp_path_factory.mktemp(f"cai-{cell}") / "Kn100")
        import yaml
        path = case / "case.yaml"
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        data["mesh"]["core_cell_size_m"] = cell
        data["checks"]["min_nozzle_faces"] = 4
        data["checks"]["max_nozzle_area_error"] = 0.5
        path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")

        _script("Allmesh", case)
        cfg = load_case_config(case)
        geom = from_config(cfg)
        plan = mesh.plan(cfg, geom, inflow.from_config(cfg))
        prediction = mesh.predicted_nozzle_faces(plan, geom)
        assert read_boundary(case)["nozzle"].n_faces == prediction.n_faces
        errors.append(abs(prediction.area_error))
    assert errors[1] < errors[0]


def test_no_boundary_is_a_wall(cai_meshed):
    """A wall REFLECTS. The plume must expand into vacuum, not into a box."""
    patches = read_boundary(cai_meshed)
    assert [name for name, info in patches.items() if info.type == "wall"] == []


def test_check_mesh_is_happy(cai_meshed):
    log = (cai_meshed / "log.checkMesh").read_text(encoding="utf-8")
    assert "Mesh OK" in log


# --------------------------------------------------------------------------- #
# the generated dictionaries
# --------------------------------------------------------------------------- #

def test_the_solver_is_standard_dsmcFoam(cai_meshed):
    """Not the MNF fork's dsmcFoam+."""
    text = (cai_meshed / "system" / "controlDict").read_text(encoding="utf-8")
    assert "application     dsmcFoam;" in text
    assert "dsmcFoam+" not in text


def test_the_argon_properties_match_case_yaml(cai_meshed):
    cfg = load_case_config(cai_meshed)
    text = (cai_meshed / "constant" / "dsmcProperties").read_text(encoding="utf-8")
    assert f"mass                            {cfg.gas.mass_kg:g};" in text
    assert f"diameter                        {cfg.gas.diameter_m:g};" in text
    assert f"omega                           {cfg.gas.omega:g};" in text
    assert "typeIdList                      (Ar);" in text


def test_only_the_nozzle_is_named_as_an_inflow(cai_meshed):
    text = (cai_meshed / "constant" / "dsmcProperties").read_text(encoding="utf-8")
    assert "patches ( nozzle );" in text
    assert "vacuum" not in text.split("plumeFieldInflowCoeffs", 1)[1]


def test_collisions_are_enabled_even_at_kn_100(cai_meshed):
    """The claim is that DSMC APPROACHES the collisionless solution, which
    requires the collision machinery to be running."""
    text = (cai_meshed / "constant" / "dsmcProperties").read_text(encoding="utf-8")
    assert "BinaryCollisionModel            VariableHardSphere;" in text


def test_the_three_knudsen_cases_produce_three_different_densities(openfoam,
                                                                   tmp_path_factory):
    """§12: a configuration labelled Kn = 0.1 must never be able to run the
    density from another case. Checked in the file the SOLVER reads."""
    densities = {}
    for knudsen in (100.0, 0.1, 0.01):
        case = _write_case(
            tmp_path_factory.mktemp(f"cai-kn{knudsen}") / "case", knudsen)
        _script("Allmesh", case)
        _script("Allrun", case, "--no-solve")
        text = (case / "0" / "boundaryNumberDensity_Ar").read_text(
            encoding="utf-8")
        block = text.split("    nozzle\n    {", 1)[1].split("}", 1)[0]
        densities[knudsen] = float(
            block.split("uniform", 1)[1].split(";", 1)[0])
    assert len(set(densities.values())) == 3
    assert densities[100.0] < densities[0.1] < densities[0.01]
    assert densities[0.1] / densities[100.0] == pytest.approx(1000.0, rel=1e-6)


# --------------------------------------------------------------------------- #
# the solver run
# --------------------------------------------------------------------------- #

def test_the_custom_inflow_model_loaded_and_chose_the_nozzle(cai_solved):
    """plumeFieldInflow announces the patches it will inject across."""
    log = (cai_solved / "log.dsmcFoam").read_text(encoding="utf-8")
    assert "plumeFieldInflow: injecting across" in log
    assert "(nozzle)" in log


def test_particles_were_inserted(cai_solved):
    log = (cai_solved / "log.dsmcFoam").read_text(encoding="utf-8")
    inserted = [int(line.split("=")[1]) for line in log.splitlines()
                if "Particles inserted" in line]
    assert inserted and max(inserted) > 0


def test_the_measured_inlet_flux_matches_the_maxwellian_flux(cai_solved):
    """Bird eq. 4.22 against what the solver actually injected.

    plumeFieldInflow accumulates a fractional particle count per face per step
    and injects the integer part, so the per-step count fluctuates by O(1) and
    the AVERAGE over the run is what can be compared. 5% is generous for a
    few-hundred-step average and still tight enough to catch a wrong density, a
    wrong area or a missing thermal correction -- each of which is tens of
    percent or more.
    """
    cfg = load_case_config(cai_solved)
    geom = from_config(cfg)
    state = inflow.from_config(cfg)
    plan = mesh.plan(cfg, geom, state)
    run = inflow.derive_run_settings(
        cfg, state, geom, min_cell_size_m=plan.min_cell_size_m,
        exit_cell_volume_m3=plan.exit_cell_volume_m3)

    log = (cai_solved / "log.dsmcFoam").read_text(encoding="utf-8")
    inserted = np.array([int(line.split("=")[1]) for line in log.splitlines()
                         if "Particles inserted" in line], dtype=float)
    assert inserted.size > 50, "too few steps to average"

    # The meshed patch area, which is what the solver injects through.
    nozzle = mesh.predicted_nozzle_faces(plan, geom)
    theoretical = (state.injection_rate_per_s(nozzle.area_m2)
                   * run.delta_t_s / run.n_equivalent_particles)

    assert inserted.mean() == pytest.approx(theoretical, rel=0.05)


def test_the_vacuum_boundaries_delete_outgoing_molecules(cai_solved):
    """If they reflected, the parcel count would grow without bound instead of
    levelling off, and molecules would appear upstream of the exit plane."""
    cfg = load_case_config(cai_solved)
    sampled = post.read_case(cai_solved, cfg.gas.mass_kg)
    assert np.all(sampled.centres[:, 0] >= 0.0)

    log = (cai_solved / "log.dsmcFoam").read_text(encoding="utf-8")
    counts = [int(line.split("=")[1]) for line in log.splitlines()
              if "Number of dsmc particles" in line]
    if len(counts) > 20:
        # Over a short run the count rises as the plume fills the domain, but a
        # reflecting boundary would keep it rising linearly with no sign of the
        # growth slowing once molecules start reaching the far field.
        assert counts[-1] > counts[0]


def test_the_solver_wrote_the_averaged_fields_the_validation_reads(cai_solved):
    time = post.latest_time(cai_solved)
    for name in post.REQUIRED_FIELDS:
        assert (cai_solved / time / name).is_file(), name


def test_allpost_produces_the_four_validation_outputs(cai_solved):
    """§9: centreline density, velocity, temperature, and a density contour."""
    results = cai_solved / "results"
    for name in ("centerline.csv", "density_plane.csv", "metrics.yaml"):
        assert (results / name).is_file(), name
    header = (results / "centerline.csv").read_text(
        encoding="utf-8").splitlines()[0].split(",")
    for column in ("x_over_D", "n_over_n0", "n_over_n0_analytical",
                   "U_sqrt_beta0", "T_over_T0"):
        assert column in header, column


def test_the_metrics_file_scores_all_three_centreline_quantities(cai_solved):
    import yaml

    document = yaml.safe_load(
        (cai_solved / "results" / "metrics.yaml").read_text(encoding="utf-8"))
    for name in ("density", "velocity", "temperature"):
        assert f"{name}_max_rel_error" in document["centerline"], name
    assert document["case"]["Kn"] == 100.0


def test_the_sampled_density_is_within_reach_of_the_analytical_solution(cai_solved):
    """A sanity check, NOT a validation: this run is a few hundred steps on a
    coarse mesh, so the statistics are poor and the transient is not over. It
    only asserts that the near field is the right order of magnitude -- a
    wrong density, a wrong area or an injecting vacuum boundary would all be
    factors, not percents.
    """
    cfg = load_case_config(cai_solved)
    geom = from_config(cfg)
    state = inflow.from_config(cfg)
    sampled = post.read_case(cai_solved, cfg.gas.mass_kg)
    profile = post.centerline(sampled, cfg, geom, state)

    near = profile["x_over_D"] < 1.0
    assert np.any(near)
    ratio = (profile["n_over_n0"][near]
             / profile["n_over_n0_analytical"][near])
    assert 0.5 < float(np.median(ratio)) < 2.0


# --------------------------------------------------------------------------- #
# continuing a run rather than overwriting it
# --------------------------------------------------------------------------- #

def _time_dirs(case):
    """Numeric time directories, ascending, excluding 0 -- as ./Allrun sees them."""
    times = []
    for entry in case.iterdir():
        if not entry.is_dir():
            continue
        try:
            value = float(entry.name)
        except ValueError:
            continue
        if value > 0.0:
            times.append((value, entry.name))
    return [name for _, name in sorted(times)]


def _total_iter(case, time: str) -> int:
    """fieldAverage's accumulated iteration count at a written time."""
    path = case / "processor0" / time / "uniform" / "functionObjects" \
        / "functionObjectProperties"
    if not path.is_file():
        path = case / time / "uniform" / "functionObjects" \
            / "functionObjectProperties"
    text = path.read_text(encoding="utf-8")
    return int(re.search(r"totalIter\s+(\d+);", text).group(1))


@pytest.fixture(scope="module")
def cai_resumable(cai_solved, tmp_path_factory):
    """A private copy of the solved case, for the resume tests to mutate."""
    case = tmp_path_factory.mktemp("cai-resume") / "Kn100"
    shutil.copytree(cai_solved, case, symlinks=True)
    for script in list(case.glob("All*")) + list(case.glob("*.py")):
        script.chmod(0o755)
    return case


def test_the_generated_control_dict_resumes_from_the_latest_time(cai_meshed):
    """The one-line change the whole workflow rests on."""
    text = (cai_meshed / "system" / "controlDict").read_text(encoding="utf-8")
    assert "startFrom       latestTime;" in text
    assert "purgeWrite      0;" in text
    # A STEP count, not seconds: the runTime schedule restarts from the resume
    # point, so a resumed leg shorter than one interval would write nothing.
    assert "writeControl    timeStep;" in text


def test_rerunning_a_finished_case_reports_instead_of_restarting(cai_resumable):
    """The case is at endTime, so there is nothing to do -- and in particular
    nothing to overwrite."""
    before = _time_dirs(cai_resumable)
    assert before, "the fixture should have written at least one time"

    result = _script("Allrun", cai_resumable)
    assert "already reached endTime" in result.stdout
    assert _time_dirs(cai_resumable) == before


def test_allmesh_refuses_to_remesh_a_case_with_results(cai_resumable):
    """blockMesh would replace the geometry the results were computed on while
    leaving them in place, addressing cells that no longer exist."""
    result = subprocess.run(
        ["bash", "./Allmesh"], cwd=str(cai_resumable), capture_output=True,
        text=True, timeout=600, env=_env(), check=False)
    assert result.returncode != 0
    assert "already holds results" in result.stderr
    assert _time_dirs(cai_resumable), "the results are still there"

    # ...but --dict-only is allowed, because that is how endTime is raised.
    # It regenerates controlDict from case.yaml, so the short-run settings this
    # module works with have to be put back afterwards.
    _script("Allmesh", cai_resumable, "--dict-only")
    _shorten(cai_resumable, steps=300)


def test_raising_the_end_time_continues_from_the_latest_time(cai_resumable):
    """The whole point: more run, not a second run.

    Checks the three things that make it a continuation rather than a restart --
    the solver starts past 0, the earlier writes survive, and fieldAverage picks
    up its accumulators instead of starting the average again.
    """
    before = _time_dirs(cai_resumable)
    resume = before[-1]
    iterations_before = _total_iter(cai_resumable, resume)

    # Set both, so this test does not depend on what the controlDict happened to
    # hold: 450 steps is three write intervals, the last exactly at endTime.
    _shorten(cai_resumable, steps=450)

    _script("Allrun", cai_resumable)

    log = (cai_resumable / "log.dsmcFoam").read_text(encoding="utf-8")
    first = float(re.search(r"^Time = (\S+)", log, re.M).group(1))
    assert first > float(resume), "the solver restarted instead of resuming"

    after = _time_dirs(cai_resumable)
    assert set(before) <= set(after), "an earlier write was destroyed"
    assert len(after) > len(before), "the resumed leg wrote nothing"

    # dsmcInitialise and decomposePar -force would each have wiped the state.
    stdout = (cai_resumable / "log.Allrun").read_text(encoding="utf-8") \
        if (cai_resumable / "log.Allrun").is_file() else ""
    del stdout

    assert _total_iter(cai_resumable, after[-1]) > iterations_before, (
        "fieldAverage restarted its average instead of continuing it")


def test_the_resumed_leg_writes_at_the_end_time(cai_resumable):
    """`writeControl timeStep` counts the GLOBAL step index, so the schedule is
    the same whether the run went straight through or was resumed part-way.
    With runTime control the schedule restarts at the resume point and a short
    leg writes nothing at all."""
    latest = _time_dirs(cai_resumable)[-1]
    end_time = float(_entry(cai_resumable, "endTime"))
    delta_t = float(_entry(cai_resumable, "deltaT"))
    assert float(latest) >= end_time - 0.5 * delta_t


# --------------------------------------------------------------------------- #
# the flux formula itself, against the solver's own constants
# --------------------------------------------------------------------------- #

def test_the_theoretical_flux_agrees_with_n_times_u_at_high_speed_ratio(cai_meshed):
    """A cross-check on the flux expression that needs no run: at S0 = 2 the
    thermal correction is a fraction of a percent, so a formula that were wrong
    by the usual factors (2, 4, sqrt(pi)) could not pass this."""
    cfg = load_case_config(cai_meshed)
    state = inflow.from_config(cfg)
    assert state.number_flux_per_m2_s == pytest.approx(
        state.number_density_per_m3 * state.velocity_m_per_s, rel=0.01)
    assert gas.maxwellian_number_flux(
        state.number_density_per_m3, 0.0, cfg.species(), state.T0_K) == (
        pytest.approx(0.25 * state.number_density_per_m3
                      * cfg.species().mean_thermal_speed(state.T0_K)))
