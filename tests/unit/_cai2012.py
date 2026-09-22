"""Helpers for the Cai 2012 unit tests: a minimal, valid case on disk.

Every test that needs a config goes through :func:`load_case`, which writes a
``case.yaml`` and loads it with the real loader. Constructing
:class:`~plumetools.cai2012.config.CaiCaseConfig` directly would skip
``validate()``, and the validation is a substantial part of what these tests are
checking.

The template is deliberately the *smallest* file that loads, not a copy of
``cases/cai2012/baseCase/case.yaml``: a test fixture that mirrors the shipped
case would pass even if the shipped case and the schema drifted apart together.
``tests/unit/test_cai2012_config.py`` loads the real one separately.
"""

from __future__ import annotations

import copy
from pathlib import Path

import yaml

from plumetools.cai2012 import geometry as geometry_module
from plumetools.cai2012 import inflow as inflow_module
from plumetools.cai2012 import mesh as mesh_module
from plumetools.cai2012.config import load_case_config

#: The smallest case.yaml this family's loader accepts, as a data structure.
#:
#: Mesh counts are small so a plan is cheap to build; the physics is Cai's.
MINIMAL_CASE = {
    "model": "cai2012_circular_plume",
    "meta": {"case_name": "test"},
    "nozzle": {"diameter_m": 0.2},
    "exit": {
        "speed_ratio": 2.0,
        "T0_K": 300.0,
        "knudsen": 100.0,
        "characteristic_length": "diameter",
    },
    "gas": {
        "species_name": "Ar",
        "mass_kg": 6.63e-26,
        "diameter_m": 4.17e-10,
        "omega": 0.74,
        "t_ref_K": 273.0,
        "internal_degrees_of_freedom": 0,
        "mean_free_path_convention": "vhs",
    },
    "geometry": {
        "symmetry_mode": "none",
        "x_max_over_D": 10.0,
        "y_half_over_D": 10.0,
        "z_half_over_D": 10.0,
    },
    "mesh": {
        "type": "graded_cartesian",
        "core_cell_size_m": 0.05,      # coarse on purpose: a fast test mesh
        "target_cell_over_mfp": 1.0,
        "max_core_cell_over_D": 0.05,
        "core_x_over_D": 2.0,
        "core_half_over_D": 1.5,
        "outer_x_cells": 6,
        "outer_lateral_cells": 5,
        "outer_expansion": 40.0,
        "max_cells": 4000000,
        "min_core_cell_size_m": None,
        "patch_names": {
            "nozzle": "nozzle",
            "upstream_vacuum": "upstreamVacuum",
            "outer": "vacuum",
            "symmetry": "symmetry",
        },
    },
    "dsmc": {
        "binary_collision_model": "VariableHardSphere",
        "collisions_enabled": True,
        "n_equivalent_particles": None,
        "numerical_particle_multiplier": 1.0,
        "run_time_multiplier": 1.0,
        "output_frequency_multiplier": 1.0,
        "delta_t_s": None,
        "courant_target": 0.2,
        "end_time_s": None,
        "average_start_s": None,
        "write_interval_s": None,
        "transient_basis": "transits",
        "transient_collision_times": 10000.0,
        "transient_domain_transits": 3.0,
        "sampling_domain_transits": 3.0,
        "initial_number_density_per_m3": 1.0e+10,
        "initial_temperature_K": 300.0,
        "n_subdomains": 4,
    },
    "resolution": {
        "target_particles_per_cell": 20.0,
        "reference_knudsen": 0.01,
        "report_only": False,
    },
    "checks": {
        "enabled": True,
        "max_cell_over_mfp": 1.0,
        "max_courant": 1.0,
        "warn_courant": 0.5,
        "min_particles_per_cell": 5.0,
        "min_nozzle_faces": 4,          # the test mesh is deliberately coarse
        "max_nozzle_area_error": 0.35,
    },
    "output": {"dialect": "standard", "inflow_model": "plumeFieldInflow"},
    "post": {
        "centerline_x_over_D_max": 10.0,
        "centerline_points": 40,
        "centerline_radius_over_D": 0.25,
        "plane": "xz",
        "plane_half_over_D": 5.0,
        "contour_levels": [0.1, 0.01, 0.001],
        "plane_points_x": 20,
        "plane_points_lateral": 20,
    },
}


def _merge(target: dict, updates: dict) -> dict:
    """Merge ``updates`` into a deep copy of ``target``, one section deep."""
    result = copy.deepcopy(target)
    for section, values in (updates or {}).items():
        if isinstance(values, dict) and isinstance(result.get(section), dict):
            result[section] = {**result[section], **values}
        else:
            result[section] = values
    return result


def write_case(directory: Path, **sections) -> Path:
    """Write a ``case.yaml`` under ``directory``, with per-section overrides.

    ``write_case(tmp_path, exit={"knudsen": 0.01})`` gives the minimal case with
    that one key changed.
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


def derived(directory: Path, **sections):
    """``(cfg, geom, exit_state, plan, run)`` -- the whole derivation chain."""
    cfg = load_case(directory, **sections)
    geom = geometry_module.from_config(cfg)
    exit_state = inflow_module.from_config(cfg)
    plan = mesh_module.plan(cfg, geom, exit_state)
    run = inflow_module.derive_run_settings(
        cfg, exit_state, geom,
        min_cell_size_m=plan.min_cell_size_m,
        exit_cell_volume_m3=plan.exit_cell_volume_m3)
    return cfg, geom, exit_state, plan, run
