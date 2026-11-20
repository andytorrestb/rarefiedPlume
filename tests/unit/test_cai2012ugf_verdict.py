"""The post-run suitability verdict (spec §7) and the post context it feeds."""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import numpy as np
import pytest
import yaml

from plumetools.cai2012ugf.context import post_context
from plumetools.cai2012ugf.settings import VerdictSettings
from plumetools.cai2012ugf.verdict import assess_fields, effective_spacing_over_mfp, resolve_limits

REPO = Path(__file__).resolve().parents[2]
TEMPLATE = REPO / "cases" / "cai2012-kn-extension" / "baseCase" / "case.yaml"
N0 = 1.0e20

COMMON = dict(
    n0=N0,
    settings=VerdictSettings(),
    limits={"max_effective_spacing_over_mfp": 1.0, "min_particles_per_dsmc_cell": 5.0},
    stationarity_verdict="stable-window-demonstrated",
    known_limitations=[],
    window_complete=True,
)


def fields_with(*, broken=None, empty_dsmc_cell=False, mfp_in_empty=0.0, all_dsmc=False):
    """Four cells: three DSMC and one USP, every gate passing unless broken."""
    ids = np.array([1.0, 1.0, 0.0, 1.0])
    if all_dsmc:
        ids[:] = 1.0
    f = {
        "collision_model_id": ids,
        "sub_cell_levels": np.ones((4, 3)),
        "particles": np.array([100.0, 50.0, 40.0, 30.0]),
        "mfp": np.array([1.0, 1.0, 0.01, 2.0]),
        "extents": np.full((4, 3), 0.1),
        "time_step_mct_ratio": np.array([0.1, 0.1, 3.0, 0.05]),
        "kn_gll": np.array([0.5, 0.2, 0.01, 1.0]),
        "number_density": N0 * np.array([1.0, 0.5, 0.9, 0.01]),
        "centres": np.array([[0.01, 0, 0], [0.2, 0, 0], [0.05, 0, 0], [1.5, 0, 0]]),
    }
    if broken == "spacing":
        f["extents"][0] = 5.0
    elif broken == "mct":
        f["time_step_mct_ratio"][0] = 0.5
    elif broken == "occupancy":
        f["particles"][1] = 2.0
    if empty_dsmc_cell:
        f = {k: np.concatenate([v, v[-1:]]) for k, v in f.items()}
        f["collision_model_id"][-1] = 1.0
        f["particles"][-1] = 0.0
        f["mfp"][-1] = mfp_in_empty
        f["number_density"][-1] = 0.0
    return f


def common(**changes):
    return {**COMMON, **changes}


# --------------------------------------------------------------------------- #
# the effective collision spacing
# --------------------------------------------------------------------------- #

# noTimeCounter.C bins a cell's N particles into all S sub-cells; a particle
# alone in its sub-cell (probability P = (1 - 1/S)^(N-1)) draws its partner from
# the WHOLE cell. The expected partner spacing is P * cell + (1 - P) * sub-cell.

def test_dense_sub_cells_resolve():
    r = effective_spacing_over_mfp(np.array([[1.0, 1, 1]]), np.array([[5, 5, 5]]),
                                   np.array([1000.0]), np.array([1.0]))
    p = (1 - 1 / 125) ** 999
    assert r[0] == pytest.approx(p * 1.0 + (1 - p) * 0.2, rel=1e-12)
    assert r[0] == pytest.approx(0.2, rel=1e-2)


def test_sparse_sub_cells_fall_back_toward_cell_size():
    r = effective_spacing_over_mfp(np.array([[1.0, 1, 1]]), np.array([[10, 10, 10]]),
                                   np.array([40.0]), np.array([0.1]))
    p = (1 - 1 / 1000) ** 39                                       # 0.962
    assert r[0] == pytest.approx((p * 1.0 + (1 - p) * 0.1) / 0.1, rel=1e-12)   # 9.66


def test_kn0p01_lip_cell_fails_the_gate():
    # cell/lambda 2.7 -> 6 sub-cells per direction (lambda/2), 40 particles:
    # 83 % of partners come from the whole cell, so the spacing is ~2.3 lambda.
    r = effective_spacing_over_mfp(np.array([[2.7, 2.7, 2.7]]), np.array([[6, 6, 6]]),
                                   np.array([40.0]), np.array([1.0]))
    assert 2.2 < r[0] < 2.4


def test_one_particle_means_cell_size():
    r = effective_spacing_over_mfp(np.array([[2.0, 1, 1]]), np.array([[10, 10, 10]]),
                                   np.array([1.0]), np.array([1.0]))
    assert r[0] == pytest.approx(2.0)


def test_undefined_mean_free_path_is_nan():
    # uniGasFoam writes GREAT (1e15) where the mean free path is undefined.
    r = effective_spacing_over_mfp(np.ones((4, 3)), np.ones((4, 3)), np.full(4, 10.0),
                                   np.array([0.0, np.inf, np.nan, 1.0e15]))
    assert np.isnan(r).all()


def test_time_step_sentinel_is_undefined_not_a_failure():
    # uniGasFoam writes GREAT (1e15) for dt/tau where the collision rate is zero.
    f = fields_with()
    f["time_step_mct_ratio"][1] = 1.0e15
    v = assess_fields(f, **COMMON)
    assert v["production_ready"] is True
    assert v["dsmc"]["undefined_time_step_ratio"] == 1
    assert v["dsmc"]["time_step_mct_ratio"]["max"] == pytest.approx(0.1)


# --------------------------------------------------------------------------- #
# the verdict
# --------------------------------------------------------------------------- #

def test_verdict_excludes_unoccupied_cells():            # Review Focus 2
    v = assess_fields(fields_with(empty_dsmc_cell=True, mfp_in_empty=0.0), **COMMON)
    assert v["dsmc"]["excluded_unoccupied"] == 1
    assert np.isfinite(v["dsmc"]["effective_spacing_over_mfp"]["max"])
    assert v["production_ready"] is True


def test_verdict_with_no_usp_cells():                    # Review Focus 5
    v = assess_fields(fields_with(all_dsmc=True), **COMMON)
    assert v["usp"]["fraction"] == 0.0 and v["usp"]["x_range_m"] is None
    assert v["status"] in ("production-candidate", "non-production")


def test_each_gate_names_its_failure():
    for broken, phrase in [("spacing", "collision spacing"), ("mct", "time step"),
                           ("occupancy", "particles")]:
        v = assess_fields(fields_with(broken=broken), **COMMON)
        assert v["production_ready"] is False and any(phrase in u for u in v["unresolved"]), broken
    v = assess_fields(fields_with(), **common(stationarity_verdict="not-demonstrated"))
    assert v["production_ready"] is False and any("stationarity" in u for u in v["unresolved"])


def test_occupancy_counts_only_cells_above_the_floor():
    f = fields_with()
    f["particles"][3] = 1.0             # n/n0 = 1e-2 is above the 1e-3 floor: fails
    assert assess_fields(f, **COMMON)["production_ready"] is False
    f["number_density"][3] = N0 * 1e-4  # now below the floor: not counted
    assert assess_fields(f, **COMMON)["production_ready"] is True


def test_usp_cells_are_not_gated():
    v = assess_fields(fields_with(), **COMMON)      # the USP cell has mct 3.0, mfp 0.01
    assert v["usp"]["n_cells"] == 1 and v["production_ready"] is True


def test_incomplete_sampling_window_blocks_production():
    # A walltime-killed, reconstructed case must not be judged on part of its window.
    v = assess_fields(fields_with(), **common(window_complete=False))
    assert v["production_ready"] is False
    assert any("sampling window incomplete" in u for u in v["unresolved"])


class Checks:
    max_cell_over_mfp = 1.0
    min_particles_per_cell = 5.0


def test_limits_default_to_case_checks():
    assert resolve_limits(VerdictSettings(), Checks()) == {
        "max_effective_spacing_over_mfp": 1.0, "min_particles_per_dsmc_cell": 5.0}


def test_limits_may_tighten_but_never_loosen():
    tight = VerdictSettings(max_effective_spacing_over_mfp=0.5, min_particles_per_dsmc_cell=10.0)
    assert resolve_limits(tight, Checks())["max_effective_spacing_over_mfp"] == 0.5
    for loose in (VerdictSettings(max_effective_spacing_over_mfp=2.0),
                  VerdictSettings(min_particles_per_dsmc_cell=1.0)):
        with pytest.raises(ValueError, match="looser"):
            resolve_limits(loose, Checks())


def test_known_limitation_blocks_production():
    v = assess_fields(fields_with(), **common(known_limitations=["direct-exit inflow surface"]))
    assert v["status"] == "non-production"
    assert "direct-exit inflow surface" in v["unresolved"]


def test_all_gates_pass():
    v = assess_fields(fields_with(), **COMMON)
    assert v["status"] == "production-candidate" and v["production_ready"] is True
    assert v["unresolved"] == []
    assert set(v["dsmc"]["effective_spacing_over_mfp"]) >= {"max", "p50", "p95", "p99", "fraction_over_limit"}


# --------------------------------------------------------------------------- #
# the post context
# --------------------------------------------------------------------------- #

@pytest.fixture
def kn_case_dir(tmp_path):
    text = TEMPLATE.read_text(encoding="utf-8")
    text, count = re.subn(r"^(\s*knudsen:\s*)10\.0", r"\g<1>1.0", text, flags=re.MULTILINE)
    assert count == 1
    (tmp_path / "case.yaml").write_text(text, encoding="utf-8")
    return tmp_path


def test_context_pending_without_verdict(kn_case_dir):   # Review Focus 5
    c = post_context(kn_case_dir)
    assert c.verdict["status"] == "post-run-verdict-pending" and c.verdict["production_ready"] is False
    assert c.mode == "direct" and c.refusal is None
    assert c.cfg.exit.knudsen == 1.0


def test_context_reads_written_verdict(kn_case_dir):
    (kn_case_dir / "results").mkdir()
    (kn_case_dir / "results/unigas_verdict.yaml").write_text(
        yaml.safe_dump({"status": "non-production", "production_ready": False, "unresolved": ["x"]}))
    assert post_context(kn_case_dir).verdict["unresolved"] == ["x"]
