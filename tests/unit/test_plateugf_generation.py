"""The plate family's generator (spec §1, §3, §5; plan Task 5)."""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest
import yaml

from plumetools.nozzleplume import workflow as npw
from plumetools.plateugf import generation

ROOT = Path(__file__).resolve().parents[2]
FAMILY = ROOT / "cases" / "plate-unigas"


@pytest.fixture
def family(tmp_path):
    shutil.copy(FAMILY / "study.yaml", tmp_path / "study.yaml")
    shutil.copytree(FAMILY / "baseCase", tmp_path / "baseCase")
    return tmp_path


def test_21_cases_named_and_parameterised(family):
    cases = generation.plan_cases(generation.load_study(family))
    assert len(cases) == 21
    names = [c.name for c in cases]
    assert names[0] == "Kn100_zD25" and names[-1] == "Kn0p0001_zD100"
    assert len(set(names)) == 21
    by = {c.name: c for c in cases}
    assert by["Kn0p01_zD50"].knudsen == 0.01
    assert by["Kn0p01_zD50"].standoff_m == pytest.approx(10.0)
    assert by["Kn1_zD100"].x_max_m == pytest.approx(20.0 + 0.25 + 1.75)


def test_dry_run_writes_nothing(family):
    entries = generation.generate(family, dry_run=True, estimate=False)
    assert len(entries) == 21
    assert not (family / "Cases").exists()


def test_generated_case_carries_its_parameters(family):
    generation.generate(family, only=["Kn0p1_zD50"], estimate=False)
    case = family / "Cases" / "Kn0p1_zD50"
    cfg = yaml.safe_load((case / "case.yaml").read_text(encoding="utf-8"))
    assert cfg["exit"]["knudsen"] == 0.1
    assert cfg["plate"]["standoff_m"] == pytest.approx(10.0)
    assert cfg["domain"]["x_max_m"] == pytest.approx(12.0)
    assert cfg["meta"]["case_name"] == "Kn0p1_zD50"
    assert any(r["level"] == 6 for r in cfg["mesh"]["refinement_regions"])
    for name in ("unigas.yaml", "Allmesh", "Allrun", "runCase.sh", "runMesh.sh"):
        assert (case / name).is_file()


def test_refuses_case_with_output(family):
    generation.generate(family, only=["Kn100_zD25"], estimate=False)
    (family / "Cases" / "Kn100_zD25" / "0.0125").mkdir()
    with pytest.raises(generation.GenerationError, match="solver output"):
        generation.generate(family, only=["Kn100_zD25"], estimate=False)


def test_densest_farthest_case_loads_and_passes_checks(family):
    generation.generate(family, only=["Kn0p0001_zD100"], estimate=False)
    case = npw.load(family / "Cases" / "Kn0p0001_zD100")
    assert not case.checks().failures
    assert case.plan.exit_level == 6


def test_smoke_case_is_small(family):
    path = generation.stage_smoke(family)
    assert path == family / "Smoke" / "Kn100_zD25"
    case = npw.load(path)
    assert sum(case.plan.cells_by_level.values()) < 200_000
    assert not case.checks().failures
    assert case.cfg.dsmc.n_subdomains == 12
