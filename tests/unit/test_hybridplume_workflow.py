"""hybridplume: configuration, generated dictionaries, and workflow refusals."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from _hybridmesh import box_mesh
from plumetools.hybridplume import cfd, config, dsmc, workflow
from plumetools.hybridplume import mesh as hm
from plumetools.hybridplume.physics import BIRD_ARGON, SourceState
from plumetools.mesh.boundary import read_boundary

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / "cases" / "cai2012-hybrid" / "hybrid.yaml"


def test_every_profile_validates_and_is_full_3d():
    raw = yaml.safe_load(STUDY.read_text())
    for name in ["base"] + list(raw["profiles"]):
        cfg = config.load(STUDY, name)
        src = config.build_source(cfg)
        plan = hm.plan(cfg, src)
        text = hm.render_block_mesh_dict(plan)
        for forbidden in ("wedge", "empty", "symmetry"):
            assert f"type {forbidden}" not in text
        assert plan.divisions[1] == plan.divisions[2]   # no half domain


@pytest.mark.parametrize("key", ["mesh.axisymmetric", "mesh.wedge", "dsmc.radial_weighting",
                                 "domain.symmetry"])
def test_symmetry_reductions_are_refused(key):
    section, item = key.split(".")
    with pytest.raises(config.ConfigError):
        config.load(STUDY, "smoke-kn0p01", {section: {item: True}})


def test_background_floor_must_be_a_floor():
    with pytest.raises(config.ConfigError):
        config.load(STUDY, "smoke-kn0p01", {"cfd": {"background_pressure_ratio": 0.5}})


def test_departures_from_paper_are_labelled():
    cfg = config.load(STUDY, "smoke-kn1em3")
    dep = config.departures_from_paper(cfg)
    assert dep["source.knudsen"]["label"].startswith("extension")


def test_cfd_case_uses_vhs_transport_and_caps_the_viscous_step(tmp_path):
    cfg = config.load(STUDY, "smoke-kn0p01")
    src = config.build_source(cfg)
    plan = hm.plan(cfg, src)
    out = cfd.write_case(tmp_path, cfg, src, plan)
    assert out["transport_table_max_rel_error"] < 1e-4
    thermo = (tmp_path / "constant" / "thermophysicalProperties").read_text()
    assert "transport       tabulated;" in thermo and "sutherland" not in thermo
    ctrl = (tmp_path / "system" / "controlDict").read_text()
    assert 'libs ("libhybridPlumeThermo.so");' in ctrl
    tm = out["timing"]
    assert tm["max_delta_t_s"] <= tm["viscous_dt_limit_s"] * (1 + 1e-12)
    p = (tmp_path / "0" / "p").read_text()
    assert "waveTransmissive" in p and f"{src.pressure:.10g}" in p


def test_dsmc_dictionaries(tmp_path):
    cfg = config.load(STUDY, "smoke-collisionless")
    src = config.build_source(cfg)
    run = dsmc.RunPlan(1e12, 1e-6, 10, 20, 2e-4, 1e-4, 1e-3, {})
    text = dsmc.render_dsmc_properties(cfg, src, run, ["source"])
    assert "BinaryCollisionModel  none;" in text          # v2512's name, not NoBinaryCollision
    assert "patches ( source );" in text and "axisymmetric" not in text.lower()
    on = config.load(STUDY, "smoke-kn0p01")
    text = dsmc.render_dsmc_properties(on, config.build_source(on), run, ["interface", "source"])
    assert "VariableHardSphere" in text and "patches ( interface source );" in text
    ctrl = dsmc.render_control_dict(cfg, run)
    assert "writeControl    timeStep;" in ctrl and "timeStart 0.0001;" in ctrl


def test_dsmc_fields_carry_values_only_on_inflow_patches(tmp_path):
    box_mesh(tmp_path, np.linspace(0, 0.2, 5), np.linspace(-0.15, 0.15, 7),
             np.linspace(-0.15, 0.15, 7), 0.1)
    patches = read_boundary(tmp_path)
    src = SourceState(BIRD_ARGON)
    k = patches["source"].n_faces
    from plumetools.hybridplume import mapping
    dsmc.write_fields(tmp_path, patches, "Ar", {"source": mapping.source_inflow(k, src)}, src)
    t = (tmp_path / "0" / "boundaryT").read_text()
    assert t.count("nonuniform") == 1 and "class       volScalarField;" in t
    n = (tmp_path / "0" / "boundaryNumberDensity_Ar").read_text()
    assert f"{src.number_density:.12g}" in n
    for name, *_ in dsmc.MEASUREMENT_FIELDS:
        assert (tmp_path / "0" / name).is_file()


def test_zero_face_interface_patch_is_appended(tmp_path):
    box_mesh(tmp_path, np.linspace(0, 0.2, 3), np.linspace(-0.15, 0.15, 3),
             np.linspace(-0.15, 0.15, 3), 0.1)
    workflow._add_zero_face_patch(tmp_path, "interface")
    p = read_boundary(tmp_path)
    assert p["interface"].type == "patch" and p["interface"].n_faces == 0
    workflow._add_zero_face_patch(tmp_path, "interface")      # idempotent
    assert len(read_boundary(tmp_path)) == 4


def test_keep_set_format(tmp_path):
    path = dsmc.write_keep_set(tmp_path, np.array([5, 1, 3]))
    text = path.read_text()
    assert "class cellSet;" in text and "\n3\n(\n1\n3\n5\n)" in text


def test_dsmc_setup_refuses_an_unsuitable_interface(tmp_path):
    case = tmp_path / "case"
    (case / "results").mkdir(parents=True)
    cfg = config.load(STUDY, "smoke-kn0p01")
    (case / workflow.RESOLVED).write_text(yaml.safe_dump(
        {"profile": cfg.profile, "digest": cfg.digest, "config": cfg.doc}))
    (case / "results" / "interface.json").write_text(json.dumps(
        {"status": "unsuitable", "reasons": ["return_mass_fraction"], "allow_unsuitable": False}))
    with pytest.raises(workflow.WorkflowError, match="unsuitable"):
        workflow.dsmc_setup(case)


def test_interface_step_reports_full_dsmc_for_mode_none(tmp_path):
    case = tmp_path / "case"
    case.mkdir()
    cfg = config.load(STUDY, "smoke-collisionless")
    (case / workflow.RESOLVED).write_text(yaml.safe_dump(
        {"profile": cfg.profile, "digest": cfg.digest, "config": cfg.doc}))
    doc = workflow.interface(case)
    assert doc["status"] == "full_dsmc"


def test_resource_estimate_counts_parcels_not_only_cells():
    from plumetools.hybridplume import resources
    cfg = config.load(STUDY, "smoke-kn0p01")
    src = config.build_source(cfg)
    est = resources.assess_configuration(cfg, src, hm.plan(cfg, src))["full_dsmc"]
    assert est["parcels"] > 0 and est["memory_GB"] > 0
    assert est["dt_s"] <= est["dt_courant_s"] * (1 + 1e-12)
    assert "cell_over_mfp" in est and "occupancy" in est
