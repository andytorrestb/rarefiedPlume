r"""Overlap validation: a truncated interface run against the full-domain benchmark.

The question it answers is narrow: **downstream of the interface, does the
truncated run reproduce the full-domain DSMC solution at the same Kn?** If it
does, injecting the closure state at that sphere is an acceptable replacement
for simulating the dense core; if it does not, the coupling assumptions failed
at that threshold. Cai's analytical solution is NOT the reference here -- the
full-domain DSMC run is (``cases/cai2012``'s ``Kn0p01_np1x``, read-only).

The comparison region, decided before any result
------------------------------------------------
``overlap.region`` in the validation case's ``interface.yaml``: an exit-centred
spherical sector ``r_min <= R/D <= r_max``, ``theta <= theta_max``, all azimuths,
cut into bins of ``r_bin_over_D`` by ``theta_bin_deg``. It must lie outside
every interface of the validation set by ``interface_buffer_over_D`` and inside
both domains; :func:`check_region` refuses one that does not. It is persisted
in every result. It is never adjusted after the error is seen, and no density
floor is applied.

The estimator: a conservative, mesh-independent mapping
-------------------------------------------------------
Each run's cells are assigned to bins by their centre. Per bin, ``n`` is the
volume-weighted mean (molecules over volume) and ``U``, ``T`` are
molecule-weighted means of the cell values -- ``post.centerline``'s rule -- so
the bins hold the same physical quantities whatever the meshes underneath.
Bins with no cell, or no molecule, in either run are excluded and counted.
Removed upstream regions never enter: the region starts outside the interface,
and an empty bin is never scored as zero.

The metric, for q in ``n``, ``|U|``, ``T``, with the bins' physical volumes
``w_i`` (the exact sector volumes, the same for both runs)::

    E_q = sqrt( sum_i w_i (q_interface_i - q_full_i)^2 / sum_i w_i q_full_i^2 )

Targets ``E_n <= 0.05``, ``E_speed <= 0.05``, ``E_T <= 0.10`` are this
project's, not published tolerances; they are read from ``overlap.targets``.
``u_R`` and ``u_theta`` (exit-centred spherical components) and the profiles
along the axis bins and across the outermost shell are reported beside them.

Uncertainty, from block means
-----------------------------
``fieldAverage`` writes running means from the averaging start ``t_s``. The mean
over a block ``(t_a, t_b]`` is recovered exactly from two of them::

    M_ab = ((t_b - t_s) M_b - (t_a - t_s) M_a) / (t_b - t_a)

``uncertainty.blocks`` blocks give each bin a standard error
``SE = std(blocks)/sqrt(K)``. With ``sigma_i^2 = SE_if^2 + SE_full^2``, the
expected squared error carries a noise term ``N = sum w sigma^2`` and has
variance ``V = sum w^2 (4 d^2 sigma^2 + 2 sigma^4)``, so::

    E_low  = sqrt(max(0, S - N - z sqrt(V)) / Q)
    E_high = sqrt(max(0, S - N + z sqrt(V)) / Q)       S = sum w d^2, Q = sum w q_full^2

``pass`` needs ``E_high`` below the target, ``fail`` needs ``E_low`` above it;
anything else is ``inconclusive``. So is any metric whose noise floor
``sqrt(N/Q)`` alone reaches half the target, or whose runs have fewer than
three blocks: low precision is never evidence of agreement. The noise in ``Q``
is neglected (second order).

Temperature here is the translational ``3RT = <v^2> - |u|^2`` of each cell, as
everywhere else in the extension.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

from plumetools.cai2012 import post
from plumetools.cai2012ext import interface

K_BOLTZMANN = 1.380649e-23

#: The verdicts.
PASS, FAIL, INCONCLUSIVE = "pass", "fail", "inconclusive"

#: Below this many blocks a standard error is not an estimate.
MIN_BLOCKS = 3

#: Fields read per frame; the cumulative means fieldAverage writes.
FIELDS = ("rhoNMean", "rhoMMean", "momentumMean", "linearKEMean")


class OverlapConfigError(ValueError):
    """``overlap:`` in ``interface.yaml`` is malformed or incompatible."""


@dataclass(frozen=True)
class OverlapConfig:
    """The predeclared comparison."""

    baseline_case: str
    r_min_over_D: float
    r_max_over_D: float
    theta_max_deg: float
    r_bin_over_D: float
    theta_bin_deg: float
    interface_buffer_over_D: float
    target_number_density: float
    target_speed: float
    target_temperature: float
    blocks: int
    confidence: float
    low_signal_relative_se: float

    def as_dict(self) -> dict:
        return {
            "baseline_case": self.baseline_case,
            "region": {"r_min_over_D": self.r_min_over_D,
                       "r_max_over_D": self.r_max_over_D,
                       "theta_max_deg": self.theta_max_deg,
                       "r_bin_over_D": self.r_bin_over_D,
                       "theta_bin_deg": self.theta_bin_deg,
                       "interface_buffer_over_D": self.interface_buffer_over_D,
                       "centre": "the physical exit centre (exit-centred "
                                 "coordinates, as every Cai comparison)"},
            "targets": {"number_density": self.target_number_density,
                        "speed": self.target_speed,
                        "temperature": self.target_temperature,
                        "note": "project targets, not published tolerances"},
            "uncertainty": {"blocks": self.blocks, "confidence": self.confidence,
                            "method": "block means from running-mean "
                                      "differences; linearised E interval"},
            "low_signal_relative_se": self.low_signal_relative_se,
        }

    def edges(self, diameter_m: float):
        """``(radial_edges_m, polar_edges_rad)``."""
        n_r = int(round((self.r_max_over_D - self.r_min_over_D) / self.r_bin_over_D))
        n_t = int(round(self.theta_max_deg / self.theta_bin_deg))
        radial = np.linspace(self.r_min_over_D, self.r_max_over_D, n_r + 1) * diameter_m
        polar = np.radians(np.linspace(0.0, self.theta_max_deg, n_t + 1))
        return radial, polar


_KEYS = {
    "baseline_case": None,
    "region": ("r_min_over_D", "r_max_over_D", "theta_max_deg", "r_bin_over_D",
               "theta_bin_deg", "interface_buffer_over_D"),
    "targets": ("number_density", "speed", "temperature"),
    "uncertainty": ("blocks", "confidence"),
    "low_signal_relative_se": None,
}


def parse_overlap(raw, where: str = "interface.yaml") -> OverlapConfig:
    """Validate ``overlap:``. Every key is required; none is inferred."""
    if not isinstance(raw, dict):
        raise OverlapConfigError(f"{where}: 'overlap' must be a mapping")
    unknown = sorted(set(raw) - set(_KEYS))
    missing = sorted(set(_KEYS) - set(raw))
    if unknown or missing:
        raise OverlapConfigError(f"{where}: overlap keys unknown {unknown}, "
                                 f"missing {missing}")
    values = {}
    for section, keys in _KEYS.items():
        if keys is None:
            continue
        block = raw[section]
        if not isinstance(block, dict) or set(block) != set(keys):
            raise OverlapConfigError(
                f"{where}: overlap.{section} must have exactly {list(keys)}")
        for key in keys:
            value = block[key]
            if isinstance(value, bool) or not isinstance(value, (int, float)) \
                    or not math.isfinite(float(value)) or float(value) < 0:
                raise OverlapConfigError(
                    f"{where}: overlap.{section}.{key} must be a non-negative "
                    f"number, got {value!r}")
            values[f"{section}.{key}"] = value
    if not isinstance(raw["baseline_case"], str) or not raw["baseline_case"]:
        raise OverlapConfigError(f"{where}: overlap.baseline_case must be a path")
    low = raw["low_signal_relative_se"]
    if isinstance(low, bool) or not isinstance(low, (int, float)) or not low > 0:
        raise OverlapConfigError(
            f"{where}: overlap.low_signal_relative_se must be positive")
    cfg = OverlapConfig(
        baseline_case=raw["baseline_case"],
        r_min_over_D=float(values["region.r_min_over_D"]),
        r_max_over_D=float(values["region.r_max_over_D"]),
        theta_max_deg=float(values["region.theta_max_deg"]),
        r_bin_over_D=float(values["region.r_bin_over_D"]),
        theta_bin_deg=float(values["region.theta_bin_deg"]),
        interface_buffer_over_D=float(values["region.interface_buffer_over_D"]),
        target_number_density=float(values["targets.number_density"]),
        target_speed=float(values["targets.speed"]),
        target_temperature=float(values["targets.temperature"]),
        blocks=int(values["uncertainty.blocks"]),
        confidence=float(values["uncertainty.confidence"]),
        low_signal_relative_se=float(low))
    if not (0.0 < cfg.r_min_over_D < cfg.r_max_over_D):
        raise OverlapConfigError(f"{where}: need 0 < r_min_over_D < r_max_over_D")
    if not (0.0 < cfg.theta_max_deg <= 90.0):
        raise OverlapConfigError(f"{where}: theta_max_deg must be in (0, 90]")
    if not (cfg.r_bin_over_D > 0 and cfg.theta_bin_deg > 0):
        raise OverlapConfigError(f"{where}: bin sizes must be positive")
    for span, step, name in (
            (cfg.r_max_over_D - cfg.r_min_over_D, cfg.r_bin_over_D, "r_bin_over_D"),
            (cfg.theta_max_deg, cfg.theta_bin_deg, "theta_bin_deg")):
        count = span / step
        if abs(count - round(count)) > 1e-9:
            raise OverlapConfigError(
                f"{where}: overlap.region.{name} does not divide its span")
    if not (0.0 < cfg.confidence < 1.0):
        raise OverlapConfigError(f"{where}: confidence must be in (0, 1)")
    if cfg.blocks < MIN_BLOCKS:
        raise OverlapConfigError(f"{where}: at least {MIN_BLOCKS} blocks")
    return cfg


def check_region(cfg: OverlapConfig, *, interface_radii_m, diameter_m: float,
                 domain_limits_m) -> list[str]:
    """Why the region is not common to every run, or ``[]``.

    Args:
        interface_radii_m: every interface radius in the validation set (or
            ``None`` for one with no admissible sphere).
        domain_limits_m: ``(x_max, lateral_half)`` of the smaller domain.
    """
    problems = []
    radii = [r for r in interface_radii_m if r is not None]
    if radii:
        needed = max(radii) + cfg.interface_buffer_over_D * diameter_m
        if cfg.r_min_over_D * diameter_m < needed:
            problems.append(
                f"r_min_over_D = {cfg.r_min_over_D:g} is inside the largest "
                f"interface ({max(radii) / diameter_m:.3f} D) plus the "
                f"{cfg.interface_buffer_over_D:g} D buffer; the region must lie "
                f"downstream of every interface")
    x_max, lateral = domain_limits_m
    theta = math.radians(cfg.theta_max_deg)
    reach_x = cfg.r_max_over_D * diameter_m
    reach_lateral = cfg.r_max_over_D * diameter_m * math.sin(theta)
    if reach_x > x_max or reach_lateral > lateral:
        problems.append(f"r_max_over_D = {cfg.r_max_over_D:g} leaves the domain")
    return problems


# --------------------------------------------------------------------------- #
# reading a run, one frame at a time
# --------------------------------------------------------------------------- #

def average_start(case_dir: Path) -> float:
    """``fieldAverage1``'s ``timeStart`` from the case's own controlDict."""
    text = (Path(case_dir) / "system" / "controlDict").read_text(encoding="utf-8")
    match = re.search(r"timeStart\s+([-+0-9.eE]+)\s*;", text)
    if not match:
        raise post.PostError(f"{case_dir}: no fieldAverage timeStart in controlDict")
    return float(match.group(1))


def frames(case_dir: Path, t_start: float) -> list[str]:
    """Sampled time directories, oldest first."""
    out = []
    for entry in Path(case_dir).iterdir():
        try:
            value = float(entry.name)
        except ValueError:
            continue
        if entry.is_dir() and value > t_start * (1 + 1e-9) and (entry / "rhoNMean").is_file():
            out.append((value, entry.name))
    return [name for _, name in sorted(out)]


def geometry(case_dir: Path, names: list[str]):
    """``(centres, volumes)`` from the newest frame that has ``C``."""
    for name in reversed(names):
        if (Path(case_dir) / name / "C").is_file():
            centres = post.read_internal_field(Path(case_dir) / name / "C")
            v_path = Path(case_dir) / name / "V"
            if not v_path.is_file():
                raise post.PostError(f"{case_dir}/{name}: no V; ./Allpost writes it")
            return centres, post.read_internal_field(v_path)
    raise post.PostError(f"{case_dir}: no frame has cell centres; run ./Allpost")


def bin_index(centres, cfg: OverlapConfig, diameter_m: float):
    """``(selected, flat_bin)``: the cells in the region and their bin."""
    radial, polar = cfg.edges(diameter_m)
    c = np.asarray(centres, dtype=np.float64)
    big_r = np.linalg.norm(c, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        theta = np.arccos(np.clip(c[:, 0] / big_r, -1.0, 1.0))
    selected = ((big_r >= radial[0]) & (big_r < radial[-1]) & (theta <= polar[-1])
                & (c[:, 0] > 0.0))
    i_r = np.clip(np.searchsorted(radial, big_r[selected], side="right") - 1, 0,
                  len(radial) - 2)
    i_t = np.clip(np.searchsorted(polar, theta[selected], side="right") - 1, 0,
                  len(polar) - 2)
    return selected, i_r * (len(polar) - 1) + i_t


def bin_volumes(cfg: OverlapConfig, diameter_m: float) -> np.ndarray:
    """The exact sector volumes ``w_i`` -- the same for both runs."""
    radial, polar = cfg.edges(diameter_m)
    shell = (radial[1:] ** 3 - radial[:-1] ** 3) / 3.0
    cone = 2.0 * math.pi * (np.cos(polar[:-1]) - np.cos(polar[1:]))
    return (shell[:, None] * cone[None, :]).reshape(-1)


def _read(case_dir: Path, name: str, selected) -> dict:
    values = {}
    for field in FIELDS:
        array = post.read_internal_field(Path(case_dir) / name / field)
        if array.shape[0] == 1:
            array = np.broadcast_to(array, (selected.size,) + array.shape[1:])
        if array.shape[0] != selected.size:
            raise post.PostError(f"{case_dir}/{name}/{field}: {array.shape[0]} "
                                 f"values for {selected.size} cells")
        values[field] = np.array(array[selected], dtype=np.float64)
    return values


def binned_moments(cell: dict, volumes, index, n_bins: int, mass_kg: float) -> dict:
    """Per-bin ``n``, ``|U|``, ``u_R``, ``u_theta``, ``T`` from cell means.

    ``n`` volume-weighted; the velocity components and ``T`` weighted by the
    cell's molecules ``V n`` -- the rule of :func:`plumetools.cai2012.post.centerline`.
    """
    rho_m = cell["rhoMMean"]
    occupied = rho_m > 0.0
    velocity = np.zeros((rho_m.size, 3))
    velocity[occupied] = cell["momentumMean"][occupied] / rho_m[occupied, None]
    mean_sq = np.zeros(rho_m.size)
    mean_sq[occupied] = 2.0 * cell["linearKEMean"][occupied] / rho_m[occupied]
    gas = K_BOLTZMANN / mass_kg
    temperature = np.maximum(0.0, (mean_sq - (velocity ** 2).sum(axis=1)) / (3.0 * gas))
    centres = cell["centres"]
    r_hat = centres / np.linalg.norm(centres, axis=1)[:, None]
    lateral = np.hypot(centres[:, 1], centres[:, 2])
    theta_hat = np.zeros_like(r_hat)
    off = lateral > 0.0
    # e_theta = (-sin th, cos th cos ph, cos th sin ph), from the exit centre.
    cos_t = r_hat[:, 0]
    sin_t = np.sqrt(np.clip(1.0 - cos_t ** 2, 0.0, None))
    theta_hat[off, 0] = -sin_t[off]
    theta_hat[off, 1] = cos_t[off] * centres[off, 1] / lateral[off]
    theta_hat[off, 2] = cos_t[off] * centres[off, 2] / lateral[off]
    n = cell["rhoNMean"]
    molecules = volumes * n

    def total(values):
        return np.bincount(index, weights=values, minlength=n_bins)

    vol = total(volumes)
    mol = total(molecules)
    cells = np.bincount(index, minlength=n_bins)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = {
            "cells": cells,
            "molecules": mol,
            "n": np.where(vol > 0, mol / vol, np.nan),
            "speed": total(molecules * np.linalg.norm(velocity, axis=1)) / mol,
            "u_R": total(molecules * np.einsum("ij,ij->i", velocity, r_hat)) / mol,
            "u_theta": total(molecules * np.einsum("ij,ij->i", velocity, theta_hat)) / mol,
            "T": total(molecules * temperature) / mol,
        }
    for key in ("speed", "u_R", "u_theta", "T"):
        out[key] = np.where(mol > 0, out[key], np.nan)
    return out


QUANTITIES = ("n", "speed", "u_R", "u_theta", "T")


def run_bins(case_dir: Path, cfg: OverlapConfig, *, diameter_m: float,
             mass_kg: float) -> dict:
    """The final means and the block means of one run, per bin.

    Reads the geometry once, then ONLY the ``blocks + 1`` boundary frames (the
    first holding nothing: at ``t_s`` the running mean has zero weight), each
    reduced to the region's cells before the next is read.
    """
    case_dir = Path(case_dir)
    t_s = average_start(case_dir)
    names = frames(case_dir, t_s)
    if not names:
        raise post.PostError(f"{case_dir}: no sampled frame after {t_s:g} s")
    centres, volumes = geometry(case_dir, names)
    if volumes.size == 1:
        volumes = np.full(centres.shape[0], float(volumes[0]))
    selected, index = bin_index(centres, cfg, diameter_m)
    n_bins = bin_volumes(cfg, diameter_m).size
    vol = volumes[selected]

    def moments(name):
        cell = _read(case_dir, name, selected)
        cell["centres"] = centres[selected]
        return cell

    final = binned_moments(moments(names[-1]), vol, index, n_bins, mass_kg)
    # Block i ends at frame ceil(i N / K) - 1: K contiguous blocks covering
    # every frame, the last ending at the final one.
    k = min(cfg.blocks, len(names))
    boundaries = [names[int(math.ceil(i * len(names) / k)) - 1] for i in range(1, k + 1)]
    blocks = []
    previous, t_prev = None, t_s
    for name in boundaries:
        t = float(name)
        cumulative = moments(name)
        if previous is None:
            block = {f: cumulative[f] for f in FIELDS}
        else:
            block = {f: ((t - t_s) * cumulative[f] - (t_prev - t_s) * previous[f])
                     / (t - t_prev) for f in FIELDS}
        block["centres"] = cumulative["centres"]
        blocks.append(binned_moments(block, vol, index, n_bins, mass_kg))
        previous, t_prev = {f: cumulative[f] for f in FIELDS}, t
    se = {}
    for q in QUANTITIES:
        if len(blocks) >= MIN_BLOCKS:
            stack = np.vstack([b[q] for b in blocks])
            with warnings.catch_warnings():
                # A bin empty in every block has no standard error: NaN, by
                # design, and excluded downstream -- not a warning.
                warnings.simplefilter("ignore", RuntimeWarning)
                se[q] = np.nanstd(stack, axis=0, ddof=1) / math.sqrt(len(blocks))
        else:
            se[q] = np.full(n_bins, np.nan)
    return {"case": str(case_dir), "t_start": t_s, "frames": len(names),
            "final_time": names[-1], "blocks": len(blocks),
            "block_ends": boundaries, "final": final, "se": se,
            "cells_in_region": int(selected.sum())}


# --------------------------------------------------------------------------- #
# the metric
# --------------------------------------------------------------------------- #

def normalised_l2(measured, reference, weights, se_measured=None, se_reference=None,
                  *, confidence: float = 0.95) -> dict:
    """``E`` and its interval, over the bins both runs populate."""
    m = np.asarray(measured, dtype=np.float64)
    ref = np.asarray(reference, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    keep = np.isfinite(m) & np.isfinite(ref) & (w > 0)
    d = m[keep] - ref[keep]
    q = float(np.sum(w[keep] * ref[keep] ** 2))
    s = float(np.sum(w[keep] * d * d))
    out = {"E": math.sqrt(s / q) if q > 0 else float("nan"), "bins": int(keep.sum())}
    if se_measured is None or se_reference is None:
        out.update({"E_low": None, "E_high": None, "noise_floor": None})
        return out
    sigma2 = (np.nan_to_num(np.asarray(se_measured, dtype=np.float64)[keep], nan=np.inf) ** 2
              + np.nan_to_num(np.asarray(se_reference, dtype=np.float64)[keep], nan=np.inf) ** 2)
    if not np.all(np.isfinite(sigma2)):
        out.update({"E_low": None, "E_high": None, "noise_floor": None})
        return out
    from scipy import stats

    z = float(stats.norm.ppf(0.5 + 0.5 * confidence))
    noise = float(np.sum(w[keep] * sigma2))
    variance = float(np.sum(w[keep] ** 2 * (4.0 * d * d * sigma2 + 2.0 * sigma2 ** 2)))
    centre = s - noise
    out.update({
        "E_bias_corrected": math.sqrt(max(centre, 0.0) / q),
        "E_low": math.sqrt(max(centre - z * math.sqrt(variance), 0.0) / q),
        "E_high": math.sqrt(max(centre + z * math.sqrt(variance), 0.0) / q),
        "noise_floor": math.sqrt(noise / q),
    })
    return out


def verdict(metric: dict, target: float) -> str:
    """``pass`` / ``fail`` / ``inconclusive`` against one target."""
    if metric.get("E_high") is None:
        return INCONCLUSIVE
    if metric["noise_floor"] >= 0.5 * target:
        return INCONCLUSIVE
    if metric["E_high"] <= target:
        return PASS
    if metric["E_low"] > target:
        return FAIL
    return INCONCLUSIVE


def compare(truncated: dict, full: dict, cfg: OverlapConfig, *,
            diameter_m: float) -> dict:
    """Every overlap metric between two :func:`run_bins` results."""
    w = bin_volumes(cfg, diameter_m)
    a, b = truncated["final"], full["final"]
    populated = (a["cells"] > 0) & (b["cells"] > 0) & (a["molecules"] > 0) \
        & (b["molecules"] > 0)
    excluded = ~populated
    targets = {"n": cfg.target_number_density, "speed": cfg.target_speed,
               "T": cfg.target_temperature}
    metrics = {}
    for q in ("n", "speed", "T", "u_R"):
        mask = populated
        metric = normalised_l2(np.where(mask, a[q], np.nan), np.where(mask, b[q], np.nan),
                               w, truncated["se"][q], full["se"][q],
                               confidence=cfg.confidence)
        if q in targets:
            metric["target"] = targets[q]
            metric["verdict"] = verdict(metric, targets[q])
        metrics[q] = metric
    # u_theta is ~0 in both: scored against the local speed, never relatively.
    scale = np.where(populated, b["speed"], np.nan)
    with np.errstate(invalid="ignore"):
        swirl = np.abs(a["u_theta"] - b["u_theta"]) / scale
    metrics["u_theta_over_speed_max"] = float(np.nanmax(swirl)) if np.any(
        np.isfinite(swirl)) else None
    with np.errstate(invalid="ignore", divide="ignore"):
        low = populated & (full["se"]["n"] / b["n"] > cfg.low_signal_relative_se)
    primary = [metrics[q]["verdict"] for q in ("n", "speed", "T")]
    overall = FAIL if FAIL in primary else (PASS if all(v == PASS for v in primary)
                                            else INCONCLUSIVE)
    return {
        "verdict": overall,
        "metrics": metrics,
        "bins": {"total": int(w.size), "compared": int(populated.sum()),
                 "excluded_empty": int(excluded.sum()),
                 "excluded_volume_fraction": float(w[excluded].sum() / w.sum()),
                 "low_signal": int(low.sum()),
                 "low_signal_volume_fraction": float(w[low].sum() / w.sum()),
                 "low_signal_rule": f"full-run SE(n)/n > {cfg.low_signal_relative_se:g}; "
                                    f"kept in the metric, reported here"},
    }


def profiles(truncated: dict, full: dict, cfg: OverlapConfig, *,
             diameter_m: float) -> list[dict]:
    """Per-bin rows: both runs, their errors and standard errors."""
    radial, polar = cfg.edges(diameter_m)
    n_t = len(polar) - 1
    rows = []
    for flat in range(bin_volumes(cfg, diameter_m).size):
        i_r, i_t = divmod(flat, n_t)
        row = {"r_over_D": 0.5 * (radial[i_r] + radial[i_r + 1]) / diameter_m,
               "theta_deg": math.degrees(0.5 * (polar[i_t] + polar[i_t + 1])),
               "cells_interface": int(truncated["final"]["cells"][flat]),
               "cells_full": int(full["final"]["cells"][flat])}
        for q in QUANTITIES:
            x, y = truncated["final"][q][flat], full["final"][q][flat]
            row[f"{q}_interface"] = float(x)
            row[f"{q}_full"] = float(y)
            row[f"{q}_se_interface"] = float(truncated["se"][q][flat])
            row[f"{q}_se_full"] = float(full["se"][q][flat])
            row[f"{q}_relative_error"] = float((x - y) / y) if (
                q not in ("u_theta",) and np.isfinite(x) and np.isfinite(y)
                and y != 0) else float("nan")
        rows.append(row)
    return rows


def write_rows(path: Path, rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


# --------------------------------------------------------------------------- #
# one validation case
# --------------------------------------------------------------------------- #

def run(case_dir: Path, *, baseline: Path | None = None) -> dict:
    """Compare a validation case with its baseline; write ``results/overlap.*``."""
    from plumetools.cai2012ext import interface_case

    case_dir = Path(case_dir)
    derived = interface_case.derive(case_dir)
    if derived.settings.purpose != "interface-validation":
        raise OverlapConfigError(f"{case_dir}: purpose is {derived.settings.purpose}, "
                                 f"not interface-validation")
    cfg = parse_overlap(derived.raw.get("overlap"), str(case_dir / interface.SETTINGS_FILE))
    d = derived.source.diameter_m
    geom = derived.geom
    problems = check_region(
        cfg, interface_radii_m=[derived.selection.radius_m], diameter_m=d,
        domain_limits_m=(geom.x_max_m, min(geom.y_max_m, geom.z_max_m)))
    core = derived.settings.core_half_angle_deg
    if core is not None and cfg.theta_max_deg > core:
        problems.append(
            f"theta_max_deg = {cfg.theta_max_deg:g} exceeds the injected core "
            f"({core:g} deg): beyond it the truncated run lacks the excluded tail "
            f"by construction, and the metric would measure that, not the coupling")
    if problems:
        raise OverlapConfigError("; ".join(problems))
    baseline = Path(baseline) if baseline else (case_dir / cfg.baseline_case).resolve()
    if not baseline.is_dir():
        raise post.PostError(
            f"the full-domain baseline {baseline} is not on this machine. It may "
            f"be on the cluster: run this there, or copy its C, V and the "
            f"{len(FIELDS)} *Mean fields of the boundary frames.")
    mass = derived.source.mass_kg
    truncated = run_bins(case_dir, cfg, diameter_m=d, mass_kg=mass)
    full = run_bins(baseline, cfg, diameter_m=d, mass_kg=mass)
    result = compare(truncated, full, cfg, diameter_m=d)
    rows = profiles(truncated, full, cfg, diameter_m=d)
    written = write_rows(case_dir / "results" / "overlap_bins.csv", rows)
    document = {
        "verdict": result["verdict"],
        "interface": {"bird_threshold": derived.settings.bird_threshold,
                      "bird_definition": derived.settings.bird_definition,
                      "core_half_angle_deg": derived.settings.core_half_angle_deg,
                      "R_interface_over_D": derived.selection.radius_m / d,
                      "core_mass_fraction": derived.selection.core_fraction},
        "comparison": cfg.as_dict(),
        "baseline": {k: full[k] for k in ("case", "t_start", "frames", "final_time",
                                          "blocks", "cells_in_region")},
        "truncated": {k: truncated[k] for k in ("case", "t_start", "frames",
                                                "final_time", "blocks",
                                                "cells_in_region")},
        **result,
        "definitions": {
            "E_q": "sqrt(sum w (q_interface - q_full)^2 / sum w q_full^2) over "
                   "populated bins, w the exact sector volume",
            "estimator": "n volume-weighted; |U|, u_R, u_theta, T molecule-"
                         "weighted cell means per bin",
            "excluded": "bins empty in either run; never scored as zero",
            "reference": "the full-domain DSMC run, not Cai's analytical "
                         "solution (which stays the primary comparison of each "
                         "case, Allcompare.py)",
        },
        "bins_csv": str(written),
        "note": "a pass supports the coupling at this threshold over this region "
                "only; it is not a converged-result claim",
    }
    path = case_dir / "results" / "overlap.yaml"
    path.write_text("# GENERATED by plumetools.cai2012ext.overlap -- do not edit.\n"
                    + yaml.safe_dump(_plain(document), sort_keys=False),
                    encoding="utf-8", newline="\n")
    document["path"] = str(path)
    return document


def _plain(value):
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m plumetools.cai2012ext.overlap",
        description="Truncated interface run vs the full-domain benchmark.")
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--baseline", type=Path, default=None,
                        help="override overlap.baseline_case (read-only)")
    args = parser.parse_args(argv)
    try:
        document = run(args.case_dir, baseline=args.baseline)
    except (OverlapConfigError, post.PostError, ValueError) as exc:
        print(f"overlap: {exc}", file=sys.stderr)
        return 1
    print(f"overlap verdict: {document['verdict'].upper()}")
    for q in ("n", "speed", "T"):
        m = document["metrics"][q]
        print(f"  E_{q:<6} {m['E']:.4f}  interval "
              f"[{m.get('E_low')}, {m.get('E_high')}]  target {m['target']:g}  "
              f"-> {m['verdict']}")
    print(f"  bins compared {document['bins']['compared']} of "
          f"{document['bins']['total']} (excluded empty "
          f"{document['bins']['excluded_empty']})")
    print(f"    wrote {document['path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
