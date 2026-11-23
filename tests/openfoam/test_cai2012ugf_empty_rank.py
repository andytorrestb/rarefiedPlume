"""uniGasFoam survives a processor that holds no particle at an adaptation.

Every Cai case starts from near vacuum, so at the first adaptation (step 20)
every processor away from the nozzle is empty -- on the cluster's 192 ranks,
most of them. uniGasFoam filled an empty processor's empty cells with its own
minimum over occupied cells, which over no cells is VGREAT; the collision rate
of a VGREAT density at a VGREAT temperature overflows, and the run died with a
floating-point exception at step 20 (upstream, and fork commit 3cd0e92). Fork
commit bd8ba5c takes the minimum over all processors, as a serial run does.

Two ranks, split by hand along the plume axis: rank 1 holds only the cells
beyond 3/4 of the domain length. The first time step comes from the
near-vacuum start (Courant 0.5 on the 25 mm core cell at the 300 K most probable
speed, ~3.5e-5 s), so in 20 steps a parcel at U0 + 3 most probable speeds
(~1.8 km/s) covers ~1.25 m of the 1.5 m to rank 1. The test also reads, from
the solver's own imbalance line, that rank 1 really was empty at step 20.

Needs OpenFOAM v2512 sourced, the patched uniGasFoam built and the project's
Python first on PATH (Allmesh calls ``python``). About three minutes.
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
ADAPTATION_INTERVAL = 20
TERMINAL_OUTPUT_EVERY = 10      # nTerminalOutputs in the rendered controlDict
END_TIME_S = 2.5e-3             # > 60 steps at any time step the run can take

if shutil.which("uniGasFoam") is None or not os.environ.get("FOAM_APPBIN"):
    pytest.skip("uniGasFoam and OpenFOAM are not on PATH", allow_module_level=True)


def foam(utility: str) -> str:
    # OpenFOAM's own binaries: a user-bin decomposePar can shadow the release one.
    return str(Path(os.environ["FOAM_APPBIN"]) / utility)


def run(command, cwd: Path, timeout: int = 1800):
    result = subprocess.run(command, cwd=str(cwd), capture_output=True, text=True,
                            timeout=timeout, check=False)
    if result.returncode != 0:
        pytest.fail(f"{' '.join(command)} failed in {cwd}:\n{result.stdout[-3000:]}\n"
                    f"{result.stderr[-3000:]}")
    return result


def write_manual_split(case_dir: Path, x_cut: float) -> int:
    """Rank 1 gets every cell whose centre lies beyond ``x_cut``; returns its size."""
    from plumetools.foamio.volfield import read_internal_field

    run([foam("postProcess"), "-func", "writeCellCentres", "-time", "0"], case_dir)
    x = read_internal_field(case_dir / "0" / "Cx")
    for written in case_dir.glob("0/C*"):
        written.unlink()
    ranks = (x > x_cut).astype(int)
    header = "FoamFile { version 2.0; format ascii; class %s; object %s; }\n"
    (case_dir / "constant" / "cellDecomposition").write_text(
        header % ("labelList", "cellDecomposition")
        + f"{ranks.size}\n(\n" + "\n".join(map(str, ranks)) + "\n)\n", encoding="utf-8")
    (case_dir / "system" / "decomposeParDict").write_text(
        header % ("dictionary", "decomposeParDict")
        + "numberOfSubdomains 2;\nmethod manual;\n"
        + 'coeffs { dataFile "cellDecomposition"; }\n', encoding="utf-8")
    return int(ranks.sum())


@pytest.fixture(scope="module")
def empty_rank_run(tmp_path_factory):
    from plumetools.cai2012.study import load_study
    from plumetools.cai2012ugf import generation

    root = tmp_path_factory.mktemp("ugf-empty-rank")
    study = load_study(FAMILY / "study.yaml")
    case = next(c for c in study.expand() if c.knudsen_case.name == "Kn0p01")
    (root / "Cases").mkdir()
    case_dir = generation.stage(FAMILY / study.base_case, root / "Cases", case, study)

    data = yaml.safe_load((case_dir / "case.yaml").read_text(encoding="utf-8"))
    data["mesh"].update(core_cell_size_m=0.025, max_core_cell_over_D=0.2,
                        outer_x_cells=8, outer_lateral_cells=8, outer_expansion=40.0)
    data["dsmc"]["n_subdomains"] = 2
    (case_dir / "case.yaml").write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    settings = yaml.safe_load((case_dir / "unigas.yaml").read_text(encoding="utf-8"))
    settings["adaptive"]["adaptation_interval"] = ADAPTATION_INTERVAL
    (case_dir / "unigas.yaml").write_text(yaml.safe_dump(settings, sort_keys=False), encoding="utf-8")

    run(["bash", "./Allmesh"], case_dir)
    x_max = data["geometry"]["x_max_over_D"] * data["nozzle"]["diameter_m"]     # 2.0 m
    assert write_manual_split(case_dir, 0.75 * x_max) > 0, "rank 1 would hold no cell"
    run([foam("decomposePar"), "-force"], case_dir)
    for processor in case_dir.glob("processor*"):
        (processor / "0").mkdir(exist_ok=True)
    run([foam("foamDictionary"), "-entry", "endTime", "-set", str(END_TIME_S),
         "system/controlDict"], case_dir)

    solver = subprocess.run(["mpirun", "-np", "2", "uniGasFoam", "-parallel"],
                            cwd=str(case_dir), capture_output=True, text=True,
                            timeout=1800, check=False)
    log = case_dir / "log.uniGasFoam"
    log.write_text(solver.stdout + solver.stderr, encoding="utf-8")
    return solver.returncode, log


def test_an_empty_processor_survives_the_first_adaptations(empty_rank_run):
    from plumetools.cai2012ugf.pilot import parse_log

    returncode, log = empty_rank_run
    text = log.read_text(encoding="utf-8", errors="replace")
    assert returncode == 0 and "sigFpe" not in text, \
        f"uniGasFoam died (exit {returncode}):\n{text[-3000:]}"

    out = parse_log(log)
    assert len(out["time"]) * TERMINAL_OUTPUT_EVERY >= 3 * ADAPTATION_INTERVAL, \
        "the run ended before its third adaptation"
    # The scenario itself: at step 20, the first adaptation, rank 1 held no
    # particle. With two ranks the imbalance is 100 % only when one is empty.
    step_20 = ADAPTATION_INTERVAL // TERMINAL_OUTPUT_EVERY - 1
    assert out["imbalance"][step_20] == pytest.approx(100.0, abs=1e-9), \
        "rank 1 held particles at the first adaptation: move the cut downstream"
