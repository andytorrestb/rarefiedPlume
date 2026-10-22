"""The nozzle-plume family's resolution diagnostics, run on a cai2012 case.

:mod:`plumetools.rarefaction` is tested in ``test_rarefaction.py``. What is
tested here is the wiring: that a Kn-extension case -- a cai2012 ``case.yaml``,
not a nozzle-plume one -- feeds it the right gas, D and thresholds, that its
centreline follows this family's tube, and that the exit-plane check lines the
measured resolution up with what generation planned. The state is a uniform
exit Maxwellian, so every number has an exact value: at ``T0`` the VHS mean free
path of ``n0`` is ``Kn D`` by construction.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest
import yaml

from _cai2012 import write_case
from _nozzleplume import write_field
# The rarefaction tests' synthetic polyMesh, so both suites read one layout.
from test_rarefaction import write_box_mesh

from plumetools.cai2012 import inflow, mesh
from plumetools.cai2012.config import load_case_config
from plumetools.cai2012.geometry import from_config
from plumetools.cai2012.post import read_internal_field
from plumetools.cai2012ext import diagnostics, policy
from plumetools.rarefaction.settings import DiagnosticSettings, SettingsError

REPO = Path(__file__).resolve().parents[2]
FAMILY = REPO / "cases" / "cai2012-kn-extension"
TIME = "0.06"
CELL = 0.05


def cai_case(tmp_path: Path, *, knudsen: float = 10.0, parcels: float = 40.0):
    """A cai2012 case over an 8^3 box of 5 cm cells, holding the exit state."""
    write_case(tmp_path, exit={"knudsen": knudsen})
    shutil.copy(FAMILY / "baseCase" / "diagnostics.yaml",
                tmp_path / "diagnostics.yaml")
    centres = write_box_mesh(tmp_path, (8, 8, 8), (0.0, -0.2, -0.2),
                             (CELL, CELL, CELL))
    cfg = load_case_config(tmp_path)
    exit_state = inflow.from_config(cfg)
    n = np.full(len(centres), exit_state.number_density_per_m3)
    rho = n * cfg.gas.mass_kg
    u = exit_state.velocity_m_per_s
    r = exit_state.species.gas_constant_j_per_kg_k
    volume = CELL ** 3
    weight = exit_state.number_density_per_m3 * volume / parcels

    directory = tmp_path / TIME
    write_field(directory / "C", internal=[tuple(c) for c in centres],
                field_class="volVectorField")
    write_field(directory / "V", internal=[volume] * len(centres))
    write_field(directory / "rhoNMean", internal=n.tolist())
    write_field(directory / "rhoMMean", internal=rho.tolist())
    write_field(directory / "momentumMean",
                internal=[(m * u, 0.0, 0.0) for m in rho],
                field_class="volVectorField")
    write_field(directory / "linearKEMean",
                internal=(rho * (0.5 * u * u + 1.5 * r * exit_state.T0_K)).tolist())
    write_field(directory / "dsmcRhoNMean",
                internal=(n * volume / weight).tolist())
    state = directory / "uniform" / "functionObjects"
    state.mkdir(parents=True)
    (state / "functionObjectProperties").write_text(
        "fieldAverage1\n{\n" + "".join(
            f"    {name}\n    {{\n        totalIter 500;\n"
            f"        totalTime 5e-4;\n    }}\n"
            for name in ("rhoN", "rhoM", "momentum", "linearKE", "dsmcRhoN"))
        + "}\n", encoding="utf-8")
    return centres, cfg


def test_the_fields_carry_the_cases_own_gas_and_diameter(tmp_path):
    cai_case(tmp_path, knudsen=10.0)
    fields, written, document = diagnostics.run(tmp_path)
    names = {Path(p).name for p in fields}
    assert {"mfp", "KnD", "KnCellMax", "mfpResolved", "cellOverMfpMax",
            "diagValid", "parcelsPerStep"} <= names
    # n0 at T0 is Kn D by construction: KnD IS the case's Kn everywhere.
    assert read_internal_field(tmp_path / TIME / "mfp") == pytest.approx(
        10.0 * 0.2, rel=1e-9)
    assert read_internal_field(tmp_path / TIME / "KnD") == pytest.approx(
        10.0, rel=1e-9)
    assert (tmp_path / "results" / "diagnostics_summary.yaml").is_file()
    assert document["knudsen"] == 10.0
    assert document["mesh_type"] == "graded_cartesian"


def test_the_exit_plane_check_lines_measured_up_with_planned(tmp_path):
    _, cfg = cai_case(tmp_path, knudsen=1.0)
    _, _, document = diagnostics.run(tmp_path)
    exit_plane = document["exit_plane"]
    geom = from_config(cfg)
    plan = mesh.plan(cfg, geom, inflow.from_config(cfg))
    assert exit_plane["target_knudsen"] == 1.0
    assert exit_plane["measured_KnD_median"] == pytest.approx(1.0, rel=1e-9)
    assert exit_plane["planned_cell_over_mfp"] == pytest.approx(plan.cell_over_mfp)
    # The box's cells are CELL wide against lambda = Kn D.
    assert exit_plane["measured_cell_over_mfp_max_median"] == pytest.approx(
        CELL / (1.0 * 0.2), rel=1e-9)
    assert exit_plane["valid_cells"] == exit_plane["cells"] > 0


def test_the_thresholds_come_from_diagnostics_yaml(tmp_path):
    cai_case(tmp_path)
    _, _, document = diagnostics.run(tmp_path)
    settings = document["settings"]
    assert settings["mesh_budget_cells"] == 4_000_000
    assert settings["practical_min_cell_size_m"] == pytest.approx(0.010)
    assert any("graded Cartesian" in note for note in document["notes"])


def test_forty_parcels_pass_the_validity_rule_and_ten_do_not(tmp_path):
    cai_case(tmp_path / "forty", parcels=40.0)
    cai_case(tmp_path / "ten", parcels=10.0)
    _, _, forty = diagnostics.run(tmp_path / "forty")
    _, _, ten = diagnostics.run(tmp_path / "ten")
    assert forty["validity"]["valid_cell_fraction"] == 1.0
    assert ten["validity"]["valid_cell_fraction"] == 0.0
    assert ten["validity"]["invalid_by_reason"]["too_few_sampled_parcels"][
        "cells"] == 512


def test_the_centreline_is_this_familys_tube(tmp_path):
    cai_case(tmp_path)
    diagnostics.run(tmp_path)
    table = np.genfromtxt(tmp_path / "results" / "centreline_diagnostics.csv",
                          delimiter=",", names=True)
    # post.centerline_radius_over_D = 0.25 -> 5 cm: the four cells round the
    # axis at each 5 cm station of post.centerline_points = 40 over 2 m.
    assert table.size == 8
    assert np.all(table["n_cells"] == 4)
    assert table["mfp"] == pytest.approx(2.0, rel=1e-9)
    pytest.importorskip("matplotlib")
    assert (tmp_path / "results" / "centreline-diagnostics.png").is_file()


def test_a_case_without_the_settings_file_says_how_to_get_it(tmp_path):
    cai_case(tmp_path)
    (tmp_path / "diagnostics.yaml").unlink()
    with pytest.raises(SettingsError, match="--refresh"):
        diagnostics.run(tmp_path)


def test_a_settings_file_may_not_restate_the_gas(tmp_path):
    cai_case(tmp_path)
    (tmp_path / "diagnostics.yaml").write_text(
        yaml.safe_dump({"gas": {"molecular_mass_kg": 1.0}, "diagnostics": {}}),
        encoding="utf-8")
    with pytest.raises(SettingsError, match="case.yaml"):
        diagnostics.run(tmp_path)


def test_an_excluded_kn_is_never_diagnosed(tmp_path):
    write_case(tmp_path, exit={"knudsen": 0.1})
    with pytest.raises(policy.ExtensionError, match="Cai's own"):
        diagnostics.run(tmp_path)


def test_the_shipped_settings_match_the_template():
    settings = diagnostics.load_settings(FAMILY / "baseCase")
    template = load_case_config(FAMILY / "baseCase")
    assert isinstance(settings, DiagnosticSettings)
    assert settings.mesh_budget_cells == template.mesh.max_cells
    assert settings.practical_min_cell_size_m == pytest.approx(
        template.mesh.max_core_cell_over_D * template.nozzle.diameter_m)
    assert settings.cells_per_mean_free_path == 3.0
