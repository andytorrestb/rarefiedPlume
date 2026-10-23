"""The overlap comparison: region, estimator, metric, interval, verdict.

Synthetic runs only. The block means are checked against blocks the test put
into the running means itself, so the recovery is exact, not plausible.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from plumetools.cai2012ext import overlap

D = 0.2
MASS = 6.63e-26

REGION = {"r_min_over_D": 3.0, "r_max_over_D": 5.0, "theta_max_deg": 20.0,
          "r_bin_over_D": 1.0, "theta_bin_deg": 10.0, "interface_buffer_over_D": 0.5}
RAW = {"baseline_case": "../full", "region": REGION,
       "targets": {"number_density": 0.05, "speed": 0.05, "temperature": 0.10},
       "uncertainty": {"blocks": 4, "confidence": 0.95}, "low_signal_relative_se": 0.2}


@pytest.fixture
def cfg():
    return overlap.parse_overlap(RAW)


def test_every_key_is_required_and_checked():
    overlap.parse_overlap(RAW)
    for bad, match in (({**RAW, "extra": 1}, "unknown"),
                       ({k: v for k, v in RAW.items() if k != "targets"}, "missing"),
                       ({**RAW, "region": {**REGION, "r_bin_over_D": 0.7}}, "divide"),
                       ({**RAW, "region": {**REGION, "theta_max_deg": 120.0}}, "theta_max"),
                       ({**RAW, "uncertainty": {"blocks": 2, "confidence": 0.95}}, "blocks"),
                       ({**RAW, "targets": {"number_density": 0.05}}, "exactly")):
        with pytest.raises(overlap.OverlapConfigError, match=match):
            overlap.parse_overlap(bad)


def test_the_region_must_be_downstream_of_every_interface(cfg):
    assert overlap.check_region(cfg, interface_radii_m=[2.0 * D, None], diameter_m=D,
                                domain_limits_m=(10 * D, 10 * D)) == []
    problems = overlap.check_region(cfg, interface_radii_m=[2.0 * D, 2.8 * D],
                                    diameter_m=D, domain_limits_m=(10 * D, 10 * D))
    assert problems and "downstream of every interface" in problems[0]
    assert overlap.check_region(cfg, interface_radii_m=[], diameter_m=D,
                                domain_limits_m=(4 * D, 10 * D))


def test_bin_weights_are_the_exact_sector_volumes(cfg):
    w = overlap.bin_volumes(cfg, D)
    assert w.size == 2 * 2
    exact = (2 * math.pi / 3) * ((5 * D) ** 3 - (3 * D) ** 3) * (1 - math.cos(math.radians(20)))
    assert w.sum() == pytest.approx(exact, rel=1e-14)


def test_the_metric_and_its_limits():
    w = np.array([1.0, 2.0, 3.0])
    ref = np.array([1.0, 2.0, 4.0])
    assert overlap.normalised_l2(ref, ref, w)["E"] == 0.0
    assert overlap.normalised_l2(1.1 * ref, ref, w)["E"] == pytest.approx(0.1, rel=1e-12)
    # An empty bin is NaN, excluded -- not a zero that would score 100 %.
    gap = overlap.normalised_l2(np.array([1.0, np.nan, 4.0]), ref, w)
    assert gap["E"] == 0.0 and gap["bins"] == 2


def test_the_interval_decides_and_noise_is_never_agreement():
    w = np.ones(50)
    ref = np.ones(50)
    tiny = np.full(50, 1e-4)
    close = overlap.normalised_l2(1.01 * ref, ref, w, tiny, tiny)
    assert overlap.verdict(close, 0.05) == overlap.PASS
    far = overlap.normalised_l2(1.20 * ref, ref, w, tiny, tiny)
    assert overlap.verdict(far, 0.05) == overlap.FAIL
    # Identical means but noisy runs: inconclusive, not pass.
    noisy = overlap.normalised_l2(ref, ref, w, np.full(50, 0.05), np.full(50, 0.05))
    assert noisy["noise_floor"] > 0.025
    assert overlap.verdict(noisy, 0.05) == overlap.INCONCLUSIVE
    # No uncertainty estimate at all: inconclusive.
    assert overlap.verdict(overlap.normalised_l2(ref, ref, w), 0.05) == overlap.INCONCLUSIVE


def _bins(n, speed, temperature, cells):
    return {"cells": np.asarray(cells), "molecules": np.asarray(n) * 1.0,
            "n": np.asarray(n, float), "speed": np.asarray(speed, float),
            "u_R": np.asarray(speed, float), "u_theta": np.zeros(len(n)),
            "T": np.asarray(temperature, float)}


def test_compare_excludes_and_counts_empty_bins(cfg):
    w_size = overlap.bin_volumes(cfg, D).size
    se = {q: np.full(w_size, 1e-6) for q in overlap.QUANTITIES}
    full = {"final": _bins([1, 1, 1, 1], [800] * 4, [20] * 4, [5, 5, 5, 5]), "se": se}
    trunc = {"final": _bins([1, 0, 1, 1], [800, 0, 800, 800], [20, 0, 20, 20],
                            [5, 0, 5, 5]), "se": se}
    result = overlap.compare(trunc, full, cfg, diameter_m=D)
    assert result["bins"]["excluded_empty"] == 1
    assert result["metrics"]["n"]["E"] == 0.0
    assert result["verdict"] == overlap.PASS


# --------------------------------------------------------------------------- #
# reading a run: block means recovered exactly from running means
# --------------------------------------------------------------------------- #

def _field(path: Path, values) -> None:
    values = np.asarray(values, dtype=np.float64)
    kind = "vector" if values.ndim == 2 else "scalar"
    rows = "\n".join("(" + " ".join(f"{v:.17g}" for v in row) + ")" if kind == "vector"
                     else f"{row:.17g}" for row in values)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("FoamFile\n{\n    class vol%sField;\n}\n"
                    "internalField   nonuniform List<%s> \n%d\n(\n%s\n)\n;\n"
                    "boundaryField\n{\n}\n" % (kind.capitalize(), kind, len(values), rows))


def synthetic_run(case: Path, cfg, block_densities, *, t_start=1.0, dt=1.0) -> np.ndarray:
    """Cells at bin centres; running means built from the given block values."""
    radial, polar = cfg.edges(D)
    centres = []
    for i in range(len(radial) - 1):
        for j in range(len(polar) - 1):
            r = 0.5 * (radial[i] + radial[i + 1])
            t = 0.5 * (polar[j] + polar[j + 1])
            centres.append((r * math.cos(t), r * math.sin(t), 0.0))
    centres = np.array(centres)
    (case / "system").mkdir(parents=True)
    (case / "system" / "controlDict").write_text(
        f"functions {{ fieldAverage1 {{ timeStart {t_start:g}; }} }}\n")
    n_cells = len(centres)
    blocks = np.asarray(block_densities, dtype=np.float64)
    speed = 800.0
    for k in range(1, len(blocks) + 1):
        t = t_start + k * dt
        mean_n = blocks[:k].mean(axis=0) * np.ones(n_cells)
        rho = MASS * mean_n
        momentum = rho[:, None] * speed * centres / np.linalg.norm(centres, axis=1)[:, None]
        temperature = 20.0
        kinetic = 0.5 * rho * (speed ** 2 + 3.0 * 1.380649e-23 / MASS * temperature)
        directory = case / f"{t:g}"
        for name, values in (("rhoNMean", mean_n), ("rhoMMean", rho),
                             ("momentumMean", momentum), ("linearKEMean", kinetic)):
            _field(directory / name, values)
        _field(directory / "C", centres)
        _field(directory / "V", np.full(n_cells, 1e-3))
    return centres


def test_block_means_are_recovered_from_running_means(cfg, tmp_path):
    rng = np.random.default_rng(3)
    blocks = 1e20 * (1.0 + 0.05 * rng.standard_normal((4, 1)))
    synthetic_run(tmp_path / "run", cfg, blocks)
    result = overlap.run_bins(tmp_path / "run", cfg, diameter_m=D, mass_kg=MASS)
    assert result["blocks"] == 4 and result["frames"] == 4
    expected_se = blocks[:, 0].std(ddof=1) / 2.0
    np.testing.assert_allclose(result["se"]["n"], expected_se, rtol=1e-9)
    np.testing.assert_allclose(result["final"]["n"], blocks.mean(), rtol=1e-12)
    np.testing.assert_allclose(result["final"]["speed"], 800.0, rtol=1e-12)
    np.testing.assert_allclose(result["final"]["T"], 20.0, rtol=1e-9)
    np.testing.assert_allclose(result["final"]["u_theta"], 0.0, atol=1e-9)


def test_two_runs_of_one_flow_agree_and_a_biased_one_does_not(cfg, tmp_path):
    rng = np.random.default_rng(11)
    base = 1e20 * (1.0 + 0.001 * rng.standard_normal((4, 1)))
    synthetic_run(tmp_path / "a", cfg, base)
    synthetic_run(tmp_path / "b", cfg, base[::-1])
    synthetic_run(tmp_path / "c", cfg, 1.2 * base)
    a, b, c = (overlap.run_bins(tmp_path / x, cfg, diameter_m=D, mass_kg=MASS)
               for x in "abc")
    assert overlap.compare(a, b, cfg, diameter_m=D)["metrics"]["n"]["verdict"] == overlap.PASS
    assert overlap.compare(c, b, cfg, diameter_m=D)["metrics"]["n"]["verdict"] == overlap.FAIL
