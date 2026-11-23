"""The cai2012-unigas chain end to end, on a coarse Kn = 0.01 case (spec §8 step 1).

Generate -> Allmesh -> Allrun on two ranks (decompose, startup leg, sampling
leg, reconstruct, convert) -> Allpostugf (Allpost, Alldiagnose, Allaudit,
Allconverge, Allverdict, Allcompare). This proves the plumbing only: the mesh is
coarsened to ~38k cells, the run lasts 0.3 domain transits and the averaging
start is moved inside it. Nothing it produces is a result.

Needs OpenFOAM v2512 sourced and the patched uniGasFoam built (see
docs/tutorials/cai2012-unigas.md). Slow: minutes.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
import yaml

pytestmark = pytest.mark.needs_openfoam

REPO = Path(__file__).resolve().parents[2]
FAMILY = REPO / "cases" / "cai2012-unigas"
STOP_TRANSITS = 0.3     # long enough for localKnudsen to decompose and DSMC cells to form

if shutil.which("uniGasFoam") is None or shutil.which("blockMesh") is None:
    pytest.skip("uniGasFoam and OpenFOAM are not on PATH", allow_module_level=True)


def run(command, cwd: Path, timeout: int = 3600):
    result = subprocess.run(command, cwd=str(cwd), capture_output=True, text=True,
                            timeout=timeout, check=False)
    if result.returncode != 0:
        pytest.fail(f"{' '.join(command)} failed in {cwd}:\n{result.stdout[-4000:]}\n"
                    f"{result.stderr[-4000:]}")
    return result


@pytest.fixture(scope="module")
def smoke_case(tmp_path_factory):
    from plumetools.cai2012.study import load_study
    from plumetools.cai2012ext import direct_case
    from plumetools.cai2012ugf import generation

    root = tmp_path_factory.mktemp("ugf-smoke")
    study = load_study(FAMILY / "study.yaml")
    case = next(c for c in study.expand() if c.knudsen_case.name == "Kn0p01")
    (root / "Cases").mkdir()
    case_dir = generation.stage(FAMILY / study.base_case, root / "Cases", case, study)

    data = yaml.safe_load((case_dir / "case.yaml").read_text(encoding="utf-8"))
    data["mesh"].update(core_cell_size_m=0.025, max_core_cell_over_D=0.2,
                        outer_x_cells=8, outer_lateral_cells=8, outer_expansion=40.0)
    data["dsmc"]["n_subdomains"] = 2
    (case_dir / "case.yaml").write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    stop = direct_case.derive(case_dir).run.domain_transit_s * STOP_TRANSITS
    data["dsmc"]["average_start_s"] = 0.4 * stop
    # ~6 sampled writes, so Allconverge.py has more than its 4-frame minimum
    # (the frozen stationarity figure cannot draw zero regions).
    data["dsmc"]["write_interval_s"] = 0.1 * stop
    (case_dir / "case.yaml").write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")

    settings = yaml.safe_load((case_dir / "unigas.yaml").read_text(encoding="utf-8"))
    settings["decomposition"]["decomposition_interval"] = 10
    (case_dir / "unigas.yaml").write_text(yaml.safe_dump(settings, sort_keys=False), encoding="utf-8")

    run(["bash", "./Allmesh"], case_dir)
    run(["bash", "./Allrun", "--stop-at-transits", str(STOP_TRANSITS)], case_dir)
    run(["bash", "./Allpostugf"], case_dir)
    return case_dir


def latest(case_dir: Path) -> Path:
    from plumetools.cai2012ugf.convert import time_directories
    return time_directories(case_dir)[-1]


def test_converted_fields_exist(smoke_case):
    t = latest(smoke_case)
    for name in ("rhoNMean", "rhoMMean", "momentumMean", "linearKEMean", "dsmcRhoNMean", "rhoN"):
        assert (t / name).is_file(), name


def test_both_cell_types_appear(smoke_case):
    from plumetools.foamio.volfield import read_internal_field
    ids = read_internal_field(latest(smoke_case) / "uniGasCollisionModelId")
    assert set(np.unique(np.rint(ids))) == {0.0, 1.0}


def test_weight_check_holds_with_cell_weighting_off(smoke_case):
    occupancy = yaml.safe_load((smoke_case / "results" / "occupancy.yaml").read_text(encoding="utf-8"))
    assert occupancy["weight_check"]["consistent"] is True


def test_verdict_and_comparison_are_written(smoke_case):
    v = yaml.safe_load((smoke_case / "results" / "unigas_verdict.yaml").read_text(encoding="utf-8"))
    assert v["status"] in ("production-candidate", "non-production")
    assert (smoke_case / "results" / "comparison.yaml").is_file()
    provenance = yaml.safe_load(
        (smoke_case / "results" / "provenance-unigas.yaml").read_text(encoding="utf-8"))
    assert provenance["postrun_verdict"]["status"] == v["status"]


def test_both_legs_ran(smoke_case):
    log = (smoke_case / "log.uniGasFoam").read_text(encoding="utf-8", errors="replace")
    assert "leg startup" in log and "leg sampling" in log


def test_pilot_report_reads_the_run(smoke_case):
    from plumetools.cai2012ugf.pilot import report
    doc = report(smoke_case)
    assert doc["particles"]["last"] > 0 and doc["fields"]["dsmc_fraction"] > 0
