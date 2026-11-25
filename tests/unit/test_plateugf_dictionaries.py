"""uniGasFoam dictionaries for a nozzleplume plate case (spec §2, §4, §6)."""
from __future__ import annotations

import re
import shutil
from pathlib import Path

import pytest
import yaml

from plumetools.cai2012ugf import dictionaries as ugf
from plumetools.cai2012ugf.settings import load_settings
from plumetools.nozzleplume import workflow
from plumetools.plateugf import dictionaries, sizing

ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "cases" / "plate-unigas" / "baseCase"


@pytest.fixture(scope="module")
def case_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("Kn10_zD25")
    shutil.copy(TEMPLATE / "case.yaml", d / "case.yaml")
    shutil.copy(TEMPLATE / "unigas.yaml", d / "unigas.yaml")
    cfg = yaml.safe_load((d / "case.yaml").read_text(encoding="utf-8"))
    cfg["exit"]["knudsen"] = 10.0
    table = sizing.density_table(7.04, 21.0, 0.1, 2.0)
    cfg["mesh"]["refinement_regions"] = sizing.refinement_regions(
        table, knudsen=10.0, diameter_m=0.2, x_limit_m=7.04, rho_limit_m=21.0)
    (d / "case.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
    return d


@pytest.fixture(scope="module")
def case(case_dir):
    return workflow.load(case_dir)


@pytest.fixture(scope="module")
def written(case, case_dir):
    dictionaries.write_case_dictionaries(case, load_settings(case_dir))
    return case_dir


def _boundary_blocks(text: str) -> dict:
    blocks = {}
    for m in re.finditer(r"patch\s+(\w+);\s*\}\s*boundaryModel\s+(\w+);\s*\w+Properties\s*\{([^}]*)\}", text):
        blocks.setdefault(m.group(1), []).append((m.group(2), m.group(3)))
    return blocks


def test_exit_state_matches_first_family(case):
    from plumetools.cai2012 import inflow
    from plumetools.cai2012.config import load_case_config
    first = inflow.from_config(load_case_config(ROOT / "cases" / "cai2012-unigas" / "baseCase"))
    assert first.knudsen == pytest.approx(10.0)
    _, exit_view, _ = dictionaries.ugf_views(case)
    assert exit_view.number_density_per_m3 == pytest.approx(first.number_density_per_m3, rel=1e-12)
    assert exit_view.velocity_m_per_s == pytest.approx(first.velocity_m_per_s, rel=1e-12)
    assert exit_view.T0_K == pytest.approx(300.0)


def test_plate_is_diffuse_wall_at_300K(case, case_dir):
    blocks = _boundary_blocks(dictionaries.render_boundaries_dict(case, load_settings(case_dir)))
    (model, body), = blocks["plate"]
    assert model == "uniGasDiffuseWallPatch"
    assert re.search(r"velocity\s+\(0 0 0\);", body)
    assert re.search(r"temperature\s+300;", body)


def test_symmetry_has_no_boundary_entry_and_open_patches_delete(case, case_dir):
    text = dictionaries.render_boundaries_dict(case, load_settings(case_dir))
    blocks = _boundary_blocks(text)
    assert "symmetry" not in text.split("uniGasPatchBoundaries", 1)[1]
    for patch in ("upstreamVacuum", "vacuum", "nozzle"):
        assert [m for m, _ in blocks[patch]][0] == "uniGasDeletionPatch"
    general = text.split("uniGasGeneralBoundaries", 1)[1]
    assert "uniGasFreeStreamInflowPatch" in general
    assert re.search(r"patch\s+nozzle;", general)
    assert re.search(r"velocity\s+\(706\.95\d* 0 0\);", general)


def test_other_dictionaries_equal_cai2012ugf_for_same_views(case, written):
    cfg, _, run = dictionaries.ugf_views(case)
    settings = load_settings(written)
    assert (written / "constant" / "uniGasProperties").read_text(encoding="utf-8") == \
        ugf.render_uni_gas_properties(cfg, run, settings)
    assert (written / "system" / "fieldPropertiesDict").read_text(encoding="utf-8") == \
        ugf.render_field_properties_dict(cfg, run, settings)
    for name in ("snappyHexMeshDict", "blockMeshDict", "decomposeParDict",
                 "controlDict.startup", "controlDict.sampling", "boundariesDict",
                 "hybridDecompositionDict", "uniGasInitialisationDict", "fvSchemes"):
        assert (written / "system" / name).is_file(), name
    assert not (written / "constant" / "dsmcProperties").exists()


def test_sampling_writes_four_per_transit(case, written):
    text = (written / "system" / "controlDict.sampling").read_text(encoding="utf-8")
    interval = float(re.search(r"^writeInterval\s+([^;]+);", text, re.M).group(1))
    assert interval == pytest.approx(case.run.transit_s / 4.0, rel=0.1)
    startup = (written / "system" / "controlDict.startup").read_text(encoding="utf-8")
    assert float(re.search(r"^endTime\s+([^;]+);", startup, re.M).group(1)) == \
        pytest.approx(case.run.average_start_s, rel=1e-9)


def test_inflow_and_deletion_match_cai2012ugf_for_the_same_exit(case):
    from plumetools.cai2012 import inflow
    from plumetools.cai2012.config import load_case_config
    base = ROOT / "cases" / "cai2012-unigas" / "baseCase"
    first_cfg = load_case_config(base)
    first = ugf.render_boundaries_dict(first_cfg, inflow.from_config(first_cfg), load_settings(base))
    ours = dictionaries.render_boundaries_dict(case, load_settings(base))
    general = lambda text: text.split("uniGasGeneralBoundaries", 1)[1]  # noqa: E731
    assert general(ours) == general(first)
    ours_blocks = _boundary_blocks(ours)
    for patch, entries in _boundary_blocks(first).items():
        assert ours_blocks[patch] == entries
    assert set(ours_blocks) - set(_boundary_blocks(first)) == {"plate"}
