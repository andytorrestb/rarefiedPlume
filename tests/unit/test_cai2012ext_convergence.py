"""Stationarity of sampled quantities: a flat total is not convergence.

The case the module exists for is the one measured in cases/cai2012: a total
that is flat while the far field is still filling. These tests build exactly
that -- a steady near field, a far-field region still rising -- and require the
verdict to fail on the far field and to ignore the flat total.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from plumetools.cai2012ext import convergence


def frames_from(series: dict, times, *, mean_series=None) -> list:
    """Frames whose instantaneous values are ``series`` and whose running mean
    is their cumulative mean, unless given explicitly."""
    frames = []
    for i, t in enumerate(times):
        instant = {name: float(values[i]) for name, values in series.items()}
        if mean_series is None:
            mean = {name: float(np.mean(values[:i + 1]))
                    for name, values in series.items()}
        else:
            mean = {name: float(values[i]) for name, values in mean_series.items()}
        frames.append({"time": t, "instant": instant, "mean": mean})
    return frames


def all_regions(value_for) -> dict:
    names = [name for name, _, _ in convergence.REGIONS] + [
        convergence.CONTEXT_REGION]
    return {name: value_for(name) for name in names}


TIMES = np.linspace(0.05, 0.063, 37)


def test_the_trend_of_a_line_is_recovered_exactly():
    values = 2.0 + 10.0 * (TIMES - TIMES[0])
    trend = convergence.linear_trend(TIMES, values)
    span = TIMES[-1] - TIMES[0]
    assert trend["relative_change"] == pytest.approx(10.0 * span / values.mean())
    assert trend["standard_error"] == pytest.approx(0.0, abs=1e-12)


def test_noise_about_a_constant_is_stationary():
    rng = np.random.default_rng(1)
    series = all_regions(lambda _: 1.0 + 0.01 * rng.standard_normal(TIMES.size))
    result = convergence.assess(frames_from(series, TIMES))
    assert result["verdict"] == "stable-window-demonstrated"
    assert all(r["stable"] for r in result["regions"].values())


def test_a_far_field_still_filling_fails_while_the_total_looks_flat():
    rising = np.linspace(0.9, 1.0, TIMES.size)        # +11% across the window
    flat = np.ones(TIMES.size)
    series = all_regions(lambda name: rising if name.startswith(
        "centreline X/D 7-10") else flat)
    result = convergence.assess(frames_from(series, TIMES))
    assert result["verdict"] == "not-demonstrated"
    far = result["regions"]["centreline X/D 7-10"]
    assert not far["instant_stationary"]
    assert far["instant_relative_change"] == pytest.approx(0.1 / 0.95, rel=1e-6)
    assert result["regions"][convergence.CONTEXT_REGION]["stable"]


def test_the_domain_total_alone_never_decides_the_verdict():
    rising = np.linspace(0.5, 1.0, TIMES.size)
    flat = np.ones(TIMES.size)
    series = all_regions(lambda name: rising if name == convergence.CONTEXT_REGION
                         else flat)
    result = convergence.assess(frames_from(series, TIMES))
    assert result["verdict"] == "stable-window-demonstrated"
    assert not result["regions"][convergence.CONTEXT_REGION]["in_verdict"]


def test_a_running_mean_that_is_still_moving_fails():
    flat = np.ones(TIMES.size)
    drifting = np.linspace(0.8, 1.0, TIMES.size)
    result = convergence.assess(frames_from(
        all_regions(lambda _: flat), TIMES,
        mean_series=all_regions(lambda _: drifting)))
    assert result["verdict"] == "not-demonstrated"
    region = result["regions"]["centreline X/D 0-1"]
    assert not region["running_mean_settled"]
    assert region["running_mean_drift"] > 0.05


def test_a_trend_too_noisy_to_resolve_is_inconclusive_not_stable():
    """A 20% error band would pass a 10% drift; that is not evidence."""
    rng = np.random.default_rng(3)
    flat = np.ones(TIMES.size)
    noisy = 1.0 + 0.3 * rng.standard_normal(TIMES.size)
    series = all_regions(lambda name: noisy if name.startswith(
        "centreline X/D 7-10") else flat)
    means = all_regions(lambda _: flat)      # a settled running mean
    result = convergence.assess(frames_from(series, TIMES, mean_series=means))
    far = result["regions"]["centreline X/D 7-10"]
    assert not far["instant_trend_resolvable"]
    assert far["status"] == "inconclusive"
    assert result["verdict"] == "inconclusive"


def test_too_few_frames_is_never_a_pass():
    series = all_regions(lambda _: np.ones(3))
    result = convergence.assess(frames_from(series, TIMES[:3]))
    assert result["verdict"] == "insufficient-frames"


class _Post:
    centerline_radius_over_D = 0.10


class _Cfg:
    post = _Post


class _Geom:
    diameter_m = 0.2


def test_the_regions_select_the_tube_bands_and_the_shell_sectors():
    centres = np.array([
        [0.1, 0.0, 0.0],                 # tube, X/D 0.5
        [0.9, 0.0, 0.01],                # tube, X/D 4.5
        [1.9, 0.0, 0.0],                 # tube X/D 9.5 AND shell theta 0
        [1.9 * math.cos(math.radians(30)), 0.0,
         1.9 * math.sin(math.radians(30))],   # shell, 30 deg
        [0.5, 0.5, 0.5],                 # nowhere but the domain
    ])
    masks = convergence.region_masks(centres, _Cfg, _Geom)
    assert list(np.flatnonzero(masks["centreline X/D 0-1"])) == [0]
    assert list(np.flatnonzero(masks["centreline X/D 4-6"])) == [1]
    assert list(np.flatnonzero(masks["centreline X/D 7-10"])) == [2]
    assert list(np.flatnonzero(masks["shell r/D 9-10, theta 0-15 deg"])) == [2]
    assert list(np.flatnonzero(masks["shell r/D 9-10, theta 15-45 deg"])) == [3]
    assert masks[convergence.CONTEXT_REGION].all()


def test_regional_density_is_molecules_over_volume():
    masks = {"a": np.array([True, True, False])}
    out = convergence.regional_densities(np.array([2.0, 4.0, 100.0]),
                                         np.array([1.0, 3.0, 1.0]), masks, 2.0)
    assert out["a"] == pytest.approx((2.0 * 1.0 + 4.0 * 3.0) / 4.0 / 2.0)
