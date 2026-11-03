"""Restart checkpoints: plumetools.foamio.checkpoints.

The trees here mimic a real decomposed dsmcFoam case (cases/r-4d/baseCaseHPC,
8 ranks): every rank's time directory holds its fields, ``uniform/time`` and a
``lagrangian/dsmc`` cloud, and every file ends with OpenFOAM's end divider.

The defect these pin: ``Allrun`` took the resume time from ``processor0``
alone and ``hpc finished`` from the newest time on ANY rank, so a time some
rank had written only partly -- a job killed mid-write -- was resumed from or
counted as the end of the run.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plumetools.foamio import checkpoints

FOOTER = "\n// " + "*" * 73 + " //\n"
FIELDS = ("rhoNMean", "boundaryT", "boundaryU")
CLOUD = ("positions", "coordinates", "U", "typeId")


def write_time(base: Path, name: str, *, fields=FIELDS, cloud=CLOUD,
               time_state=True, truncate=()):
    directory = base / name
    directory.mkdir(parents=True)
    for field in fields:
        text = f"FoamFile {{ object {field}; }}\ninternalField uniform 0;\n"
        (directory / field).write_text(text + ("" if field in truncate else FOOTER))
    if cloud:
        (directory / "lagrangian" / "dsmc").mkdir(parents=True)
        for name_ in cloud:
            (directory / "lagrangian" / "dsmc" / name_).write_text("0()\n" + FOOTER)
    if time_state:
        (directory / "uniform").mkdir()
        (directory / "uniform" / "time").write_text(f"value {name};\n" + FOOTER)


def case_with(tmp_path, n_ranks, times, *, end_time=None, delta_t=1e-3):
    case = tmp_path / "case"
    (case / "system").mkdir(parents=True)
    (case / "system" / "controlDict").write_text(
        f"endTime {end_time if end_time is not None else 1.0};\n"
        f"deltaT {delta_t};\nstartFrom latestTime;\n")
    for rank in range(n_ranks):
        base = case / f"processor{rank}"
        (base / "constant").mkdir(parents=True)
        for name in times:
            write_time(base, name)
    return case


def test_a_complete_parallel_run_resumes_in_place(tmp_path):
    case = case_with(tmp_path, 4, ["0.1", "0.2"])
    state = checkpoints.scan(case)
    assert state.decomposed == (0.2, "0.2") and state.reconstructed is None
    plan = checkpoints.plan_restart(case, 4)
    assert (plan.resume, plan.decompose, plan.refusal) == ("0.2", "none", None)


def test_a_rank_without_the_newest_time_blocks_the_resume(tmp_path):
    """The old Allrun read processor0's 0.3 and resumed from it."""
    case = case_with(tmp_path, 4, ["0.1", "0.2"])
    for rank in range(3):                       # processor3 never wrote 0.3
        write_time(case / f"processor{rank}", "0.3")
    state = checkpoints.scan(case)
    assert state.decomposed == (0.2, "0.2")
    plan = checkpoints.plan_restart(case, 4)
    assert plan.refusal and "processor3 has no 0.3/" in plan.refusal
    assert "Nothing was changed" in plan.refusal


def test_a_truncated_field_is_incomplete(tmp_path):
    case = case_with(tmp_path, 3, ["0.1"])
    for rank in range(3):
        write_time(case / f"processor{rank}", "0.2",
                   truncate=("rhoNMean",) if rank == 1 else ())
    plan = checkpoints.plan_restart(case, 3)
    assert "processor1: cut off before the end: rhoNMean" in plan.refusal


def test_a_vanished_cloud_is_incomplete_but_an_empty_subdomain_is_not(tmp_path):
    case = case_with(tmp_path, 2, [])
    write_time(case / "processor0", "0.1")
    write_time(case / "processor1", "0.1", cloud=())     # never held a parcel
    write_time(case / "processor0", "0.2")
    write_time(case / "processor1", "0.2", cloud=())
    assert checkpoints.scan(case).decomposed == (0.2, "0.2")

    write_time(case / "processor0", "0.3", cloud=())      # lost its parcels
    write_time(case / "processor1", "0.3", cloud=())
    plan = checkpoints.plan_restart(case, 2)
    assert "wrote a 'dsmc' cloud at 0.2 but none at 0.3" in plan.refusal


def test_missing_time_state_and_missing_cloud_files_are_incomplete(tmp_path):
    case = case_with(tmp_path, 2, ["0.1"])
    write_time(case / "processor0", "0.2", time_state=False)
    write_time(case / "processor1", "0.2", cloud=("positions",))
    refusal = checkpoints.plan_restart(case, 2).refusal
    assert "processor0: no uniform/time" in refusal
    assert "processor1: cloud 'dsmc' lacks U, coordinates, typeId" in refusal


def test_finished_counts_only_complete_checkpoints(tmp_path):
    """The old 'hpc finished' took the newest time on ANY rank."""
    case = case_with(tmp_path, 2, ["0.5"], end_time=1.0, delta_t=0.01)
    write_time(case / "processor0", "1")                  # processor1 was killed
    assert not checkpoints.reached_end(case)
    write_time(case / "processor1", "1")
    assert checkpoints.reached_end(case)


def test_end_time_has_half_a_step_of_slack(tmp_path):
    case = case_with(tmp_path, 1, ["0.00990533"], end_time=0.00990534, delta_t=1e-6)
    assert checkpoints.reached_end(case)
    assert not checkpoints.at_end(0.0098, 0.00990534, 1e-6)


def test_a_new_decomposition_reconstructs_the_newer_state_first(tmp_path):
    """decomposePar -force deletes processor*/; their newest state must survive."""
    case = case_with(tmp_path, 4, ["0.1", "0.2"])
    write_time(case, "0.1")                                # reconstructed earlier
    plan = checkpoints.plan_restart(case, 8)
    assert (plan.resume, plan.decompose, plan.reconstruct_first) == ("0.2", "latest", True)


def test_a_state_already_reconstructed_is_not_reconstructed_again(tmp_path):
    """The usual end of a parallel run: reconstructPar has caught up."""
    case = case_with(tmp_path, 4, ["0.1", "0.2"])
    write_time(case, "0.2")
    plan = checkpoints.plan_restart(case, 8)
    assert (plan.resume, plan.decompose, plan.reconstruct_first) == ("0.2", "latest", False)
    assert checkpoints.plan_restart(case, 4).decompose == "none"


def test_a_newer_reconstructed_state_is_decomposed_without_reconstructing(tmp_path):
    case = case_with(tmp_path, 4, ["0.1"])
    write_time(case, "0.2")
    plan = checkpoints.plan_restart(case, 4)
    assert (plan.resume, plan.decompose, plan.reconstruct_first) == ("0.2", "latest", False)


def test_a_fresh_case_decomposes_everything(tmp_path):
    case = case_with(tmp_path, 0, [])
    plan = checkpoints.plan_restart(case, 4)
    assert (plan.resume, plan.decompose, plan.refusal) == (None, "all", None)


def test_a_serial_run_does_not_ignore_newer_processor_state(tmp_path):
    case = case_with(tmp_path, 2, ["0.1", "0.2"])
    write_time(case, "0.1")
    plan = checkpoints.plan_restart(case, 1)
    assert plan.resume == "0.1"
    assert "newer than the case directory's 0.1" in plan.refusal


def test_a_missing_rank_directory_blocks(tmp_path):
    case = case_with(tmp_path, 3, ["0.1"])
    import shutil
    shutil.rmtree(case / "processor1")
    plan = checkpoints.plan_restart(case, 3)
    assert "processor directories [1] are missing" in plan.refusal


def test_the_shell_plan_prints_only_safe_values(tmp_path):
    case = case_with(tmp_path, 2, ["0.1"])
    text = checkpoints.plan_restart(case, 2).shell()
    assert "RESUME=0.1\n" in text and "DECOMPOSE=none\n" in text
    assert "RECONSTRUCT_FIRST=0\n" in text and "PROC_COUNT=2\n" in text


@pytest.mark.parametrize("args,status", [(["finished"], 1), (["report"], 0)])
def test_the_command_line(tmp_path, capsys, args, status):
    case = case_with(tmp_path, 2, ["0.1"])
    assert checkpoints.main([args[0], str(case)]) == status


def test_plan_exits_3_with_the_reason(tmp_path, capsys):
    case = case_with(tmp_path, 2, ["0.1"])
    write_time(case / "processor0", "0.2")
    assert checkpoints.main(["plan", str(case), "--np", "2"]) == 3
    assert "processor1 has no 0.2/" in capsys.readouterr().err
