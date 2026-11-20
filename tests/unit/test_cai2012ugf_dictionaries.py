"""The uniGasFoam dictionaries of a cai2012-unigas case.

Every value is either inherited from the kn-extension case (the particle
weight, the run window, the exit state, the gas) or fixed by the approved
design (spec §4, the owner's decisions 4 and 5).
"""

from __future__ import annotations

import dataclasses
import re
from pathlib import Path

import pytest

from plumetools.cai2012.config import load_case_config
from plumetools.cai2012ext import direct_case
from plumetools.cai2012ugf.dictionaries import (
    render_boundaries_dict,
    render_control_dict,
    render_field_properties_dict,
    render_hybrid_decomposition_dict,
    render_initialisation_dict,
    render_uni_gas_properties,
    write_case_dictionaries,
)
from plumetools.cai2012ugf.settings import settings_from_dict

REPO = Path(__file__).resolve().parents[2]
KN_EXTENSION = REPO / "cases" / "cai2012-kn-extension"
SETTINGS = settings_from_dict({}, "test")


@pytest.fixture(scope="module")
def case():
    cfg = load_case_config(KN_EXTENSION / "baseCase")
    cfg = dataclasses.replace(cfg, exit=dataclasses.replace(cfg.exit, knudsen=0.001))
    return direct_case.derive_config(cfg)


# --------------------------------------------------------------------------- #
# readers for the rendered text
# --------------------------------------------------------------------------- #

def entry(text: str, key: str) -> str:
    match = re.search(rf"^\s*{re.escape(key)}\s+([^;\n]+);", text, re.MULTILINE)
    assert match, f"no '{key}' entry"
    return match.group(1).strip()


def props(case) -> str:
    return render_uni_gas_properties(case.cfg, case.run, SETTINGS)


def patch_models(text: str, patch: str) -> set:
    pairs = re.findall(r"patch\s+(\w+);\s*\}\s*boundaryModel\s+(\w+);", text)
    return {model for name, model in pairs if name == patch}


def number_density(text: str, species: str) -> float:
    match = re.search(rf"numberDensities\s*\{{\s*{species}\s+([^;]+);", text)
    assert match
    return float(match.group(1))


def inflow_velocity(text: str) -> tuple:
    block = text[text.index("uniGasFreeStreamInflowPatchProperties"):]
    match = re.search(r"velocity\s*\(([^)]*)\)", block)
    return tuple(float(v) for v in match.group(1).split())


def field_block(text: str, name: str) -> str:
    blocks = re.split(r"\bfieldModel\b", text)[1:]
    chosen = [b for b in blocks if re.search(rf"^\s*field\s+{name};", b, re.MULTILINE)]
    assert len(chosen) == 1, f"one field block for {name}"
    return chosen[0]


def reset_until(text: str, name: str) -> float:
    return float(entry(field_block(text, name), "resetAtOutputUntilTime"))


def measure_mfp(text: str, name: str) -> str:
    return entry(field_block(text, name), "measureMeanFreePath")


# --------------------------------------------------------------------------- #
# the tests
# --------------------------------------------------------------------------- #

def test_weight_is_the_kn_extension_weight(case):
    assert float(entry(props(case), "nEquivalentParticles")) == pytest.approx(
        case.run.n_equivalent_particles, rel=1e-12)


def test_cell_weighting_off_and_adaptive_keys(case):
    t = props(case)
    for key, value in [("cellWeightedSimulation", "false"), ("cellWeightAdaptation", "false"),
                       ("axisymmetricSimulation", "false"), ("subCellAdaptation", "true"),
                       ("timeStepAdaptation", "true"), ("adaptationInterval", "20"),
                       ("maxSubCellSizeMFPRatio", "0.5"), ("maxTimeStepMCTRatio", "0.2"),
                       ("maxCourantNumber", "0.5"), ("dsmcCellsTimeStepLimited", "true"),
                       ("maxSubCellLevels", "10")]:
        assert entry(t, key) == value, key


def test_argon_is_the_repository_vhs_model(case):
    t = props(case)
    assert float(entry(t, "omega")) == 0.74 and "0.81" not in t
    assert float(entry(t, "mass")) == 6.63e-26 and float(entry(t, "diameter")) == 4.17e-10
    assert entry(t, "rotationalDegreesOfFreedom") == "0" and entry(t, "vibrationalModes") == "0"


def test_collision_block(case):
    t = props(case)
    assert entry(t, "collisionModel") == "hybrid"
    assert entry(t, "bgkCollisionModel") == "unifiedStochasticParticleSBGK"
    assert entry(t, "dsmcCollisionPartnerModel") == "noTimeCounter"
    assert entry(t, "dsmcCollisionModel") == "variableHardSphere"
    assert entry(t, "Tref") == "273" and entry(t, "macroInterpolation") == "true"
    assert float(entry(t, "theta")) == 1e-3


def test_nozzle_injects_and_deletes_vacuum_only_deletes(case):
    t = render_boundaries_dict(case.cfg, case.exit_state, SETTINGS)
    assert patch_models(t, "nozzle") == {"uniGasDeletionPatch", "uniGasFreeStreamInflowPatch"}
    assert patch_models(t, "upstreamVacuum") == {"uniGasDeletionPatch"}
    assert patch_models(t, "vacuum") == {"uniGasDeletionPatch"}
    assert number_density(t, "Ar") == pytest.approx(case.exit_state.number_density_per_m3, rel=1e-12)
    assert inflow_velocity(t) == pytest.approx((case.exit_state.velocity_m_per_s, 0.0, 0.0), rel=1e-12)
    assert float(entry(t, "translationalTemperature")) == case.exit_state.T0_K


def test_averaging_window(case):
    t = render_field_properties_dict(case.cfg, case.run, SETTINGS)
    assert reset_until(t, "Ar") == pytest.approx(case.run.average_start_s, rel=1e-12)
    assert reset_until(t, "ArWindow") > case.run.end_time_s
    assert measure_mfp(t, "Ar") == "true"
    assert entry(field_block(t, "Ar"), "averagingAcrossManyRuns") == "true"


def test_two_legs(case):
    s = render_control_dict(case.cfg, case.run, SETTINGS, "startup")
    p = render_control_dict(case.cfg, case.run, SETTINGS, "sampling")
    assert entry(s, "application") == entry(p, "application") == "uniGasFoam"
    assert float(entry(s, "endTime")) == pytest.approx(case.run.average_start_s, rel=1e-12)
    assert float(entry(s, "writeInterval")) == pytest.approx(case.run.domain_transit_s, rel=1e-12)
    assert float(entry(p, "endTime")) == pytest.approx(case.run.end_time_s, rel=1e-12)
    assert float(entry(p, "writeInterval")) == pytest.approx(case.run.write_interval_s, rel=1e-12)
    assert entry(p, "writeControl") == "adjustableRunTime" and entry(p, "purgeWrite") == "0"
    assert entry(p, "startFrom") == "latestTime" and entry(p, "writeFormat") == "ascii"


def test_unknown_leg_is_refused(case):
    with pytest.raises(ValueError, match="leg"):
        render_control_dict(case.cfg, case.run, SETTINGS, "transient")


def test_decomposition(case):
    t = render_hybrid_decomposition_dict(case.run, SETTINGS)
    assert entry(t, "decompositionModel") == "localKnudsen"
    assert (entry(t, "breakdownMax"), entry(t, "smoothingPasses"), entry(t, "decompositionInterval")) == ("0.05", "5", "100")
    assert float(entry(t, "theta")) == 0.2 and entry(t, "resetAtDecomposition") == "on"
    assert float(entry(t, "resetAtDecompositionUntilTime")) == pytest.approx(case.run.average_start_s, rel=1e-12)


def test_initialisation_is_the_suite_near_vacuum(case):
    t = render_initialisation_dict(case.cfg)
    assert entry(t, "type") == "uniGasMeshFill"
    assert number_density(t, "Ar") == case.cfg.dsmc.initial_number_density_per_m3
    assert float(entry(t, "translationalTemperature")) == case.cfg.dsmc.initial_temperature_K


def test_write_case_dictionaries_writes_the_set(case, tmp_path):
    names = {p.relative_to(tmp_path).as_posix()
             for p in write_case_dictionaries(tmp_path, case.cfg, case.exit_state, case.run, SETTINGS)}
    assert {"constant/uniGasProperties", "system/controlDict", "system/controlDict.startup",
            "system/controlDict.sampling", "system/chemReactDict", "system/decomposeParDict",
            "system/boundariesDict", "system/fieldPropertiesDict", "system/uniGasInitialisationDict",
            "system/hybridDecompositionDict", "system/fvSchemes", "system/fvSolution"} <= names
    assert (tmp_path / "system/controlDict").read_text() == (tmp_path / "system/controlDict.startup").read_text()
    assert (tmp_path / "0").is_dir()
