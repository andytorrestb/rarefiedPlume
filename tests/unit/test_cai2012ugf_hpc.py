"""Run legs, cloud-safe restarts, build checks and records (spec §6, §11).

uniGasFoam keeps only the newest particle cloud (``cleanLagrangian``) and,
resuming from a time without ``lagrangian/uniGas``, silently re-fills the
subdomain from the initialisation dictionary. A rank that holds no particle
legitimately writes no cloud (OpenFOAM ``Cloud::writeObject`` writes only when
``size() > 0``), so a missing cloud is a refusal only when the rank's own
density at that time says it held gas: then its cloud was deleted.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
import yaml

from plumetools.cai2012ext import direct_case
from plumetools.cai2012ugf import dictionaries
from plumetools.cai2012ugf.hpc import (
    check_build,
    leg_for,
    record,
    resume_plan,
    write_provenance,
)
from plumetools.cai2012ugf.settings import settings_from_dict

REPO = Path(__file__).resolve().parents[2]
TEMPLATE = REPO / "cases" / "cai2012-unigas" / "baseCase"
DIVIDER = "// " + "*" * 73 + " //\n"


def make_case(tmp_path: Path, *, with_dictionaries: bool = False, commit: str = "def") -> Path:
    family = tmp_path / "family"
    case = family / "Cases" / "Kn1_np1x"
    case.mkdir(parents=True)
    text = (TEMPLATE / "case.yaml").read_text(encoding="utf-8")
    data = yaml.safe_load(text)
    data["meta"]["unigas_commit"] = commit
    data["meta"]["solver"] = "uniGasFoam hybrid USP-DSMC"
    (case / "case.yaml").write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    shutil.copy(TEMPLATE / "unigas.yaml", case / "unigas.yaml")
    (case / "system").mkdir()
    (case / "system" / "controlDict").write_text(
        "FoamFile { object controlDict; }\nendTime 0.0628;\ndeltaT 1e-7;\n" + DIVIDER)
    if with_dictionaries:
        d = direct_case.derive(case)
        dictionaries.write_case_dictionaries(case, d.cfg, d.exit_state, d.run,
                                             settings_from_dict({}, "test"))
    return case


def write_time(base: Path, time: str, state: str, delta_t: str = "1e-09") -> None:
    """``state``: ``cloud`` (gas, cloud), ``none`` (no gas, no cloud),
    ``deleted`` (gas, cloud gone), ``partial`` (gas, cut off mid-write)."""
    t = base / time
    (t / "uniform").mkdir(parents=True, exist_ok=True)
    (t / "uniform" / "time").write_text(f"value {time};\ndeltaT {delta_t};\n" + DIVIDER)
    density = "0" if state == "none" else "1e18"
    (t / "rhoN_Ar").write_text(f"internalField uniform {density};\n"
                               + ("" if state == "partial" else DIVIDER))
    if state in ("cloud", "partial"):
        c = t / "lagrangian" / "uniGas"
        c.mkdir(parents=True, exist_ok=True)
        (c / "positions").write_text("0()\n" + DIVIDER)


def make_ranks(tmp_path: Path, *, n: int, times: dict, delta_t: str = "1e-09") -> Path:
    """``times``: {time: {rank: state}}; see :func:`write_time`."""
    case = make_case(tmp_path)
    for time, ranks in times.items():
        for rank in range(n):
            write_time(case / f"processor{rank}", time, ranks.get(rank, "none"), delta_t)
    return case


def test_legs():
    assert leg_for(None, 0.05, 0.063) == "startup"
    assert leg_for(0.049, 0.05, 0.063) == "startup"
    assert leg_for(0.05 * (1 - 1e-12), 0.05, 0.063) == "sampling"
    assert leg_for(0.063, 0.05, 0.063) == "done"


def test_legs_read_rounded_time_directory_names():
    # timePrecision 8 names the leg-end write 5.6580853e-05 for an endTime of
    # 5.65808533701079e-05: that write IS the end of the startup leg.
    assert leg_for(5.6580853e-05, 5.65808533701079e-05, 1.4e-4) == "sampling"
    assert leg_for(1.4145213e-04, 5.7e-05, 1.41452133425e-04) == "done"
    assert leg_for(5.65e-05, 5.65808533701079e-05, 1.4e-4) == "startup"


def test_legs_allow_half_a_step_short():
    # OpenFOAM writes within half a step of a write time and stops within half a
    # step of endTime, so a leg's last write can fall up to dt/2 short of its end.
    # 1e-7 short: beyond the 1e-6 relative tolerance (5e-8 here), within half a 2e-7 step.
    assert leg_for(0.0501435 - 1e-7, 0.0501435, 0.063, slack_s=1.02e-7) == "sampling"
    assert leg_for(0.0501435 - 1e-7, 0.0501435, 0.063) == "startup"
    assert leg_for(0.063 - 1e-7, 0.0501435, 0.063, slack_s=1.02e-7) == "done"


def test_resume_reads_the_step_from_uniform_time(tmp_path):
    # 1e-7 short of average_start_s = 0.0501435 (beyond the 1e-6 relative
    # tolerance), with a step of 2e-7: within half a step, so the leg is over.
    case = make_ranks(tmp_path, n=2, times={"0.0501434": {0: "cloud", 1: "cloud"}},
                      delta_t="2e-07")
    assert resume_plan(case, 2).leg == "sampling"


def test_resume_refuses_time_without_cloud_on_every_rank(tmp_path):    # Review Focus 1
    # Rank 0 finished writing 0.06 and deleted its 0.055 cloud; rank 1 was cut
    # off mid-write at 0.06, and 0.06 was then moved aside -- the frozen
    # planner's own advice. Rank 0 held gas at 0.055 but has no cloud there.
    case = make_ranks(tmp_path, n=2, times={"0.055": {0: "deleted", 1: "cloud"}})
    plan = resume_plan(case, 2)
    assert plan.refusal and "lagrangian/uniGas" in plan.refusal and "manual restart" in plan.refusal
    assert "processor0" in plan.refusal


def test_partial_newer_time_is_refused_and_names_the_cloud_loss(tmp_path):
    case = make_ranks(tmp_path, n=2, times={"0.055": {0: "deleted", 1: "cloud"},
                                            "0.06": {0: "cloud", 1: "partial"}})
    plan = resume_plan(case, 2)
    assert plan.refusal and "0.06" in plan.refusal and "lagrangian/uniGas" in plan.refusal


def test_resume_uses_newest_time_with_complete_clouds(tmp_path):
    case = make_ranks(tmp_path, n=2, times={"0.055": {0: "deleted", 1: "deleted"},
                                            "0.06": {0: "cloud", 1: "cloud"}})
    plan = resume_plan(case, 2)
    assert plan.refusal is None and plan.base.resume == "0.06" and plan.leg == "sampling"


def test_a_rank_that_never_held_particles_is_not_refused(tmp_path):
    case = make_ranks(tmp_path, n=2, times={"0.06": {0: "cloud", 1: "none"}})
    plan = resume_plan(case, 2)
    assert plan.refusal is None and plan.base.resume == "0.06"


def test_fresh_case_starts_in_startup(tmp_path):
    plan = resume_plan(make_case(tmp_path), 2)
    assert plan.leg == "startup" and plan.base.resume is None and plan.refusal is None


def test_shell_lines_carry_leg_and_refusal(tmp_path):
    lines = resume_plan(make_case(tmp_path), 2).shell()
    assert "LEG=startup\n" in lines and "REFUSAL=0\n" in lines


def test_check_build_mismatch_is_reported(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "uniGasFoam.commit").write_text("abc\n")
    monkeypatch.setenv("FOAM_USER_APPBIN", str(bin_dir))
    message = check_build(make_case(tmp_path, commit="def"))
    assert message and "abc" in message and "def" in message


def test_check_build_match_passes(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "uniGasFoam.commit").write_text("def\n")
    monkeypatch.setenv("FOAM_USER_APPBIN", str(bin_dir))
    assert check_build(make_case(tmp_path, commit="def")) is None


def test_record_appends_with_solver(tmp_path):
    case = make_case(tmp_path)
    record(case, 100.0)
    path = record(case, 123.0)
    rows = yaml.safe_load(path.read_text(encoding="utf-8"))["measurements"]
    assert len(rows) == 2 and rows[-1]["solver"] == "uniGasFoam hybrid USP-DSMC"
    assert rows[-1]["seconds"] == 123.0 and rows[-1]["unigas_commit"] == "def"
    assert path == tmp_path / "family" / "results" / "cost-model.yaml"


def test_provenance_records_solver_commit_and_patch_keys(tmp_path, monkeypatch):
    monkeypatch.setenv("WM_PROJECT_VERSION", "v2512")
    p = yaml.safe_load(write_provenance(make_case(tmp_path, with_dictionaries=True)).read_text())
    assert p["solver"] == "uniGasFoam hybrid USP-DSMC" and p["openfoam_version"] == "v2512"
    assert p["uniGasProperties"]["dsmcCellsTimeStepLimited"] == "true"
    assert p["uniGasProperties"]["cellWeightedSimulation"] == "false"
    assert p["hybridDecompositionDict"]["breakdownMax"] == "0.05"
    assert p["unigas_commit"] == "def"
