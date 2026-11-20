"""The pilot report (spec §8) and the dsmcFoam cross-check (spec §9)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from plumetools.cai2012ugf.crosscheck import report as crosscheck_report
from plumetools.cai2012ugf.pilot import parse_log, projected_hours


def write_log(tmp_path: Path, *, steps) -> Path:
    """``steps``: (time, dt, particles, imbalance %, clock s) per terminal output."""
    lines = ["Starting time loop", ""]
    for time, dt, particles, imbalance, clock in steps:
        lines += [f"Time = {time:g}", "",
                  f"    Time step                       = {dt:.10g}",
                  f"    Number of particles             = {particles}",
                  f"    Maximum imbalance               = {imbalance}",
                  f"ExecutionTime = {clock} s  ClockTime = {clock:g} s", ""]
    path = tmp_path / "log.uniGasFoam"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def two_families(tmp_path: Path, *, ugf_cases, dsmc_cases):
    roots = {}
    for name, cases_dir, cases, value in (("ugf", "Cases", ugf_cases, 1.10),
                                          ("dsmc", "Cases_full_suite", dsmc_cases, 1.00)):
        root = tmp_path / name
        root.mkdir()
        (root / "study.yaml").write_text(yaml.safe_dump({"cases_dir": cases_dir}))
        for case in cases:
            results = root / cases_dir / case / "results"
            results.mkdir(parents=True)
            (results / "comparison.yaml").write_text(yaml.safe_dump(
                {"centreline": {"n_over_n0": {"rel_rms": value, "valid": 120}},
                 "status": "non-production"}))
        roots[name] = root
    return roots["ugf"], roots["dsmc"]


def test_parse_log_reads_the_solver_lines(tmp_path):
    log = write_log(tmp_path, steps=[(1e-7, 1e-7, 1000, 3.0, 10.0), (2e-7, 1e-7, 1100, 4.0, 20.0)])
    out = parse_log(log)
    assert list(out["particles"]) == [1000, 1100] and out["clock"][-1] == 20.0
    assert list(out["time"]) == [1e-7, 2e-7] and list(out["imbalance"]) == [3.0, 4.0]
    assert list(out["dt"]) == [1e-7, 1e-7]


def test_parse_log_spans_appended_legs(tmp_path):
    log = write_log(tmp_path, steps=[(1e-7, 1e-7, 1000, 3.0, 10.0)])
    with open(log, "a", encoding="utf-8") as handle:
        handle.write("\n// ===== Allrun: leg sampling =====\nTime = 3e-07\n\n"
                     "    Time step                       = 5e-08\n"
                     "    Number of particles             = 1200\n"
                     "ExecutionTime = 5 s  ClockTime = 5 s\n")
    out = parse_log(log)
    assert list(out["time"]) == [1e-7, 3e-7] and list(out["particles"]) == [1000, 1200]


def test_projection_uses_recent_rate():
    assert projected_hours(t_now=0.01, end=0.02, recent_clock_per_sim_s=3.6e6,
                           elapsed_s=0) == pytest.approx(10.0)
    assert projected_hours(t_now=0.01, end=0.02, recent_clock_per_sim_s=3.6e6,
                           elapsed_s=7200) == pytest.approx(12.0)


def test_crosscheck_pairs_metrics_and_reports_missing(tmp_path):
    ugf, dsmc = two_families(tmp_path, ugf_cases=["Kn1_np1x"], dsmc_cases=["Kn1_np1x", "Kn0p1_np1x"])
    r = crosscheck_report(ugf, dsmc)
    row = r["cases"]["Kn1_np1x"]["metrics"]["centreline.n_over_n0.rel_rms"]
    assert row["dsmcFoam"] == 1.00 and row["uniGasFoam"] == 1.10
    assert row["difference"] == pytest.approx(0.10)
    assert r["cases"]["Kn0p1_np1x"]["missing"] == ["uniGasFoam"]
    assert "Co 0.5" in r["numerics"]


def test_crosscheck_writes_its_files(tmp_path):
    ugf, dsmc = two_families(tmp_path, ugf_cases=["Kn1_np1x"], dsmc_cases=["Kn1_np1x"])
    crosscheck_report(ugf, dsmc)
    assert (ugf / "results" / "crosscheck.yaml").is_file()
    text = (ugf / "results" / "crosscheck.md").read_text(encoding="utf-8")
    assert "not validated" in text.lower() and "Kn1_np1x" in text
