"""Kinematic-viscosity cap for the coarse CFD profile; no solver."""
from pathlib import Path

import pytest

from plumetools.hybridplume import cfd, config, mesh

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / "cases" / "cai2012-hybrid" / "hybrid.yaml"


def load(name, overlay=None):
    return config.load(STUDY, name, overlay)


def test_nu_max_equals_background_kinematic_viscosity():
    cfg = load("hpc-kn1em3-capsule")
    src = config.build_source(cfg)
    tm = cfd.timing(cfg, src, mesh.plan(cfg, src))
    bg = cfd.background_state(cfg, src)
    assert tm["nu_max_m2_s"] == pytest.approx(
        src.gas.viscosity(bg["temperature_K"]) / bg["density_kg_m3"], rel=1e-12)


def test_thermo_dict_selects_capped_transport_only_when_enabled():
    cfg = load("hpc-kn1em3-capsule")
    src = config.build_source(cfg)
    nu = cfd.timing(cfg, src, mesh.plan(cfg, src))["nu_max_m2_s"]
    capped = cfd.render_thermophysical_properties(src, nu_max=nu)
    assert "transport       cappedTabulated;" in capped
    assert "nuMax" in capped
    ell = load("hpc-kn1em3-ellipsoid")
    esrc = config.build_source(ell)
    assert cfd.timing(ell, esrc, mesh.plan(ell, esrc))["nu_max_m2_s"] is None
    plain = cfd.render_thermophysical_properties(esrc)
    assert "nuMax" not in plain and "cappedTabulated" not in plain
    assert cfd.render_thermophysical_properties(esrc, nu_max=None) == plain


@pytest.mark.parametrize("factor", [0, -1, float("nan"), float("inf")])
def test_viscosity_cap_factor_must_be_positive(factor):
    with pytest.raises(config.ConfigError):
        load("hpc-kn1em3-capsule", {"cfd": {"viscosity_cap": {"factor": factor}}})


# --------------------------------------------------------------------------- #
# accountability: the cap must never act inside the retained continuum region
# --------------------------------------------------------------------------- #

def _capsule_case(tmp_path, times=("0.001",), overlay=None):
    """A 2x2x2-cell capsule-profile case with uniform states at ``times``."""
    import yaml
    from _hybridmesh import box_mesh
    from plumetools.hybridplume import workflow
    case = tmp_path / "capsule"
    cfg = load("hpc-kn1em3-capsule", overlay)
    case.mkdir()
    (case / workflow.RESOLVED).write_text(yaml.safe_dump(
        {"profile": cfg.profile, "digest": cfg.digest, "config": cfg.doc}))
    box_mesh(case / "cfd", [0, .1, .2], [-.1, 0, .1], [-.1, 0, .1], .1)
    for t in times:
        td = case / "cfd" / t
        td.mkdir()
        for name, kind, row in (("p", "scalar", "2"), ("T", "scalar", "300"),
                                ("rho", "scalar", "0.0001"), ("U", "vector", "(150 0 0)")):
            (td / name).write_text(f"FoamFile {{ version 2.0; format ascii; object {name}; }}\n"
                                   f"internalField nonuniform List<{kind}>\n8\n(\n"
                                   + "\n".join([row] * 8) + "\n);\n")
    return case


def test_capped_cells_flags_only_low_density():
    import numpy as np
    from plumetools.hybridplume import breakdown as bd
    cfg = load("hpc-kn1em3-capsule")
    src = config.build_source(cfg)
    bg = cfd.background_state(cfg, src)
    nu = cfd.timing(cfg, src, mesh.plan(cfg, src))["nu_max_m2_s"]
    rho = np.array([bg["density_kg_m3"], bg["density_kg_m3"] / 100])
    state = bd.CfdState(rho=rho, T=np.full(2, 300.0), U=np.zeros((2, 3)),
                        p=rho * src.gas.gas_constant * 300.0, time="1")
    assert bd.capped_cells(state, src.gas, nu).tolist() == [False, True]


def test_breakdown_reports_uncapped_retained_region(tmp_path):
    from plumetools.hybridplume import workflow
    doc = workflow.breakdown(_capsule_case(tmp_path), time="0.001")
    cap = doc["viscosity_cap"]
    assert cap["nu_max_m2_s"] > 0
    assert cap["capped_cells"] == 0 and cap["capped_in_retained"] == 0
    assert doc["summary"]["retained_cells"] > 0
    assert (tmp_path / "capsule/cfd/0.001/viscosityCapped").is_file()


def test_breakdown_errors_when_capped_cell_is_retained(tmp_path, monkeypatch):
    from plumetools.hybridplume import breakdown as bd
    from plumetools.hybridplume import workflow
    case = _capsule_case(tmp_path)
    real = bd.read_state

    def thin(case_dir, time, gas):
        # The same uniform (hence retained) flow, at a density the cap acts on.
        s = real(case_dir, time, gas)
        return bd.CfdState(rho=s.rho * 1e-6, T=s.T, U=s.U, p=s.p * 1e-6, time=s.time)
    monkeypatch.setattr(bd, "read_state", thin)
    with pytest.raises(workflow.WorkflowError, match="viscosity cap active in retained continuum"):
        workflow.breakdown(case, time="0.001")


def test_breakdown_omits_cap_section_when_disabled(tmp_path):
    from plumetools.hybridplume import workflow
    case = _capsule_case(tmp_path, overlay={"cfd": {"viscosity_cap": {"enabled": False}}})
    doc = workflow.breakdown(case, time="0.001")
    assert "viscosity_cap" not in doc
    assert not (case / "cfd/0.001/viscosityCapped").exists()


# --------------------------------------------------------------------------- #
# retained-region stability between the last two writes
# --------------------------------------------------------------------------- #

def test_retained_change_is_volume_weighted():
    import numpy as np
    from plumetools.hybridplume import breakdown as bd
    v = np.array([1.0, 3.0, 6.0])
    prev = np.array([True, True, False])   # 4
    cur = np.array([True, False, True])    # 7
    assert bd.retained_change(prev, cur, v) == pytest.approx(3.0 / 7.0)
    assert bd.retained_change(cur, cur, v) == 0.0


def test_capsule_profile_sets_retained_stability_limit():
    assert load("hpc-kn1em3-capsule").get("breakdown.retained_stability_limit") == 0.02
    assert load("hpc-kn1em3-ellipsoid").get("breakdown.retained_stability_limit") is None


def test_retained_stability_passes_between_identical_writes(tmp_path):
    from plumetools.hybridplume import workflow
    doc = workflow.breakdown(_capsule_case(tmp_path, ("0.001", "0.002")), time="0.002")
    assert doc["retained_stability"] == {"previous_write": "0.001",
                                         "volume_fraction_change": 0.0, "limit": 0.02}


def test_retained_stability_errors_above_limit(tmp_path, monkeypatch):
    from plumetools.hybridplume import breakdown as bd
    from plumetools.hybridplume import workflow
    case = _capsule_case(tmp_path, ("0.001", "0.002"))
    monkeypatch.setattr(bd, "retained_change", lambda *a: 0.03)
    with pytest.raises(workflow.WorkflowError, match="retained volume fraction changed"):
        workflow.breakdown(case, time="0.002")


def test_failed_gate_withdraws_an_earlier_breakdown(tmp_path, monkeypatch):
    from plumetools.hybridplume import breakdown as bd
    from plumetools.hybridplume import fitted, workflow
    case = _capsule_case(tmp_path, ("0.001", "0.002"))
    workflow.breakdown(case, time="0.001")
    fitted._load_breakdown(case)   # published and current
    monkeypatch.setattr(bd, "retained_change", lambda *a: 0.03)
    with pytest.raises(workflow.WorkflowError, match="retained volume fraction changed"):
        workflow.breakdown(case, time="0.002")
    results = case / "results"
    assert not (results / "breakdown.npz").exists()
    assert not (results / "breakdown.json").exists()
    assert (results / "breakdown.failed.json").is_file()
    with pytest.raises(workflow.WorkflowError, match="run the breakdown step first"):
        fitted._load_breakdown(case)
    with pytest.raises(workflow.WorkflowError):
        workflow.interface(case)


def test_retained_stability_skips_without_previous_write(tmp_path):
    from plumetools.hybridplume import workflow
    doc = workflow.breakdown(_capsule_case(tmp_path), time="0.001")
    assert doc["retained_stability"] == {"previous_write": None,
                                         "volume_fraction_change": None, "limit": 0.02}


def test_retained_stability_section_absent_without_limit(tmp_path):
    from plumetools.hybridplume import workflow
    case = _capsule_case(tmp_path, ("0.001", "0.002"),
                         overlay={"breakdown": {"retained_stability_limit": None}})
    assert "retained_stability" not in workflow.breakdown(case, time="0.002")
