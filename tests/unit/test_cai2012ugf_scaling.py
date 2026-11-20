"""The scaling inventory of the cai2012-unigas family (spec §10).

The frozen kn-extension inventory builds every row; this family relabels each
one: the solver, a stratum of its own (never pooled with dsmcFoam), and the
suitability decided by the post-run verdict instead of the dsmcFoam pre-run
feasibility that USP cells fail by design.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
import yaml

from plumetools.cai2012ext import provenance
from plumetools.cai2012ugf import generation
from plumetools.cai2012ugf.scaling import inventory

REPO = Path(__file__).resolve().parents[2]
UGF = REPO / "cases" / "cai2012-unigas"
FIELD = ("FoamFile { version 2.0; format ascii; class volScalarField; object x; }\n"
         "dimensions [0 0 0 0 0 0 0];\ninternalField uniform 1;\nboundaryField { }\n")


@pytest.fixture
def family_with_two_results(tmp_path, monkeypatch):
    monkeypatch.setattr(generation, "git_provenance", lambda root: {"commit": "stub"})
    monkeypatch.setattr(provenance, "git_provenance", lambda root: {"commit": "stub"})
    family = tmp_path / "cai2012-unigas"
    shutil.copytree(UGF, family, ignore=shutil.ignore_patterns(
        "Cases", "Pilot", "results", "manifest-*.yaml", "__pycache__"))
    generation.generate(family)
    for name, ready in (("Kn1_np1x", True), ("Kn0p001_np1x", True)):
        case = family / "Cases" / name
        t = case / "0.06"
        t.mkdir()
        for field in ("C", "V", "rhoNMean"):
            (t / field).write_text(FIELD, encoding="utf-8")
        results = case / "results"
        results.mkdir()
        (results / "unigas_verdict.yaml").write_text(yaml.safe_dump(
            {"status": "production-candidate" if ready else "non-production",
             "production_ready": ready, "unresolved": []}), encoding="utf-8")
        for doc in ("diagnostics_summary.yaml", "occupancy.yaml", "stationarity.yaml"):
            (results / doc).write_text("{}\n", encoding="utf-8")
    return family


def test_rows_carry_the_solver_and_a_separate_stratum(family_with_two_results):
    rows = inventory(family_with_two_results)["cases"]
    assert {r["solver"] for r in rows} == {"uniGasFoam hybrid USP-DSMC"}
    assert all(r["stratum"].endswith("+uniGasFoam") for r in rows if r["stratum"])


def test_suitability_comes_from_the_post_run_verdict(family_with_two_results):
    rows = {r["name"]: r for r in inventory(family_with_two_results)["cases"]}
    assert rows["Kn1_np1x"]["fit_ready"] is True
    assert rows["Kn1_np1x"]["suitable"] is True          # verdict production_ready, fit-ready fields
    assert rows["Kn0p001_np1x"]["suitable"] is False     # inflow-surface limitation
    assert not any("Pre-run" in w or "feasibility" in w.lower() for w in rows["Kn1_np1x"]["warnings"])
    assert rows["Kn1_np1x"]["suitability"]["criteria"] == "cai2012ugf post-run verdict (spec §7)"


def test_a_case_without_a_verdict_is_not_suitable(family_with_two_results):
    rows = {r["name"]: r for r in inventory(family_with_two_results)["cases"]}
    assert rows["Kn10_np1x"]["suitable"] is False
    assert any("post-run verdict not computed" in r for r in rows["Kn10_np1x"]["suitability"]["reasons"])


def test_strata_are_regrouped_under_the_new_labels(family_with_two_results):
    inv = inventory(family_with_two_results)
    assert inv["solver"] == "uniGasFoam hybrid USP-DSMC"
    assert all(key.endswith("+uniGasFoam") for key in inv["compatibility_strata"])
