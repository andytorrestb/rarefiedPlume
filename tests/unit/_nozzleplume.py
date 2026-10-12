"""Helpers for the nozzle-plume unit tests: a minimal, valid case on disk.

Every test that needs a config goes through :func:`load_case`, which writes a
``case.yaml`` and loads it with the real loader. Constructing
:class:`~plumetools.nozzleplume.config.CaseConfig` directly would skip
``validate()``, and the validation is a substantial part of what these tests are
checking.

The template is deliberately the *smallest* file that loads, not a copy of
``cases/r-4d/baseCase/case.yaml``: a fixture that mirrored the shipped case
would pass even if the shipped case and the schema drifted apart together.
``tests/unit/test_nozzleplume_config.py`` loads the real one separately.
"""

from __future__ import annotations

import copy
from pathlib import Path

import yaml

from plumetools.nozzleplume import exitstate, geometry, mesh, run
from plumetools.nozzleplume.config import load_case_config

REPO = Path(__file__).resolve().parents[2]

#: The shipped R-4D case.
R4D_CASE = REPO / "cases" / "r-4d" / "baseCase"

#: The smallest case.yaml this family's loader accepts, as a data structure.
#:
#: Round numbers, so an expected value can be worked out by hand: a 0.4 m exit,
#: a 4 x 2 m plate 4 m downstream, a 6 x 6 x 4 m box of 0.5 m cells.
MINIMAL_CASE = {
    "model": "nozzle_exit_plume",
    "meta": {"case_name": "test"},
    "nozzle": {"exit_diameter_m": 0.4},
    "exit": {
        "temperature_K": 500.0,
        "speed_m_per_s": 3000.0,
        "number_density_per_m3": 1.0e+22,
    },
    "gas": {
        "species_name": "exhaust",
        "molecular_mass_kg": 3.4e-26,
        "molecular_diameter_m": 3.5e-10,
    },
    "plate": {
        "standoff_m": 4.0,
        "thickness_m": 0.25,
        "size_y_m": 4.0,
        "size_z_m": 2.0,
    },
    "domain": {"x_max_m": 6.0, "y_half_m": 3.0, "z_half_m": 2.0},
    "mesh": {
        "background_cell_size_m": 0.5,
        "plate_refinement_level": 2,
        "refinement_regions": [
            {"name": "exitCore", "level": 3, "x_max_m": 1.0, "radius_m": 0.5},
        ],
    },
    "dsmc": {"n_equivalent_particles": 1.0e+17},
}


def _merge(target: dict, updates: dict) -> dict:
    """Merge ``updates`` into a deep copy of ``target``, one section deep.

    A section given as ``None`` is removed, which is how a test asks for a case
    with no plate.
    """
    result = copy.deepcopy(target)
    for section, values in (updates or {}).items():
        if values is None:
            result.pop(section, None)
        elif isinstance(values, dict) and isinstance(result.get(section), dict):
            result[section] = {**result[section], **values}
        else:
            result[section] = values
    return result


def write_case(directory: Path, **sections) -> Path:
    """Write a ``case.yaml`` under ``directory``, with per-section overrides.

    ``write_case(tmp_path, exit={"temperature_K": 300.0})`` gives the minimal
    case with that one key changed; ``plate=None`` gives it with no plate.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "case.yaml"
    path.write_text(
        yaml.safe_dump(_merge(MINIMAL_CASE, sections), sort_keys=False),
        encoding="utf-8")
    return path


def load_case(directory: Path, **sections):
    """Write and load a case config in one step."""
    write_case(directory, **sections)
    return load_case_config(directory)


def derive(cfg):
    """``(geom, exit_state, plan, nozzle, run)`` for a loaded config."""
    geom = geometry.from_config(cfg)
    exit_state = exitstate.from_config(cfg)
    plan = mesh.plan(cfg, geom, exit_state)
    nozzle = mesh.predicted_nozzle_faces(plan, geom)
    settings = run.derive(cfg, exit_state, geom, plan, nozzle)
    return geom, exit_state, plan, nozzle, settings


def derived(directory: Path, **sections):
    """``(cfg, geom, exit_state, plan, nozzle, run)`` -- the whole chain."""
    cfg = load_case(directory, **sections)
    return (cfg, *derive(cfg))
