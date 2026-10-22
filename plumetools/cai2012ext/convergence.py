r"""Whether the sampling window was steady where it matters: in the far field.

A flat total parcel count is **not** evidence that the far field has converged.
The count is dominated by the dense near field, and it plateaus while the far
centreline is still being filled by the slow tail of the exit distribution --
measured in ``cases/cai2012`` at Kn = 100: the count was flat at 1.32 domain
transits and the centreline density at X/D = 9.7 was still 5.6% low after two
(``docs/cai2012-case.md`` section 5.1). So this module watches **regional
sampled quantities** across the written frames of a run and asks two questions
of each:

**Is the instantaneous field stationary over the averaging window?** Each frame
also carries the single-step ``rhoN``. Its regional mean is fitted against time
by least squares; the relative change the fit implies across the window is
compared with :data:`DEFAULT_TOLERANCE` and with twice its own standard error.
A field still filling in shows up as a trend the noise cannot explain.

**Has the running mean stopped moving?** ``rhoNMean`` at the last frame against
``rhoNMean`` at the frame nearest the middle of the window. A running mean of a
stationary signal moves less and less; one averaging a transient keeps moving.

A trend test on noisy data can pass by having a wide error band. So a region is
``stable`` only if the test could have resolved a drift of
:data:`RESOLUTION_FACTOR` tolerances; otherwise it is ``inconclusive``, which
is reported as such and never counted as a pass.

The regions are the centreline tube in bands of ``X/D`` (the estimator
``post.centerline`` uses) and the r/D = 9-10 shell in polar sectors (where the
angular profiles are taken). The whole-domain total is reported as context and
**excluded from the verdict**, for the reason above.

Both tests need frames: fewer than :data:`MIN_FRAMES` sampled frames gives
``insufficient-frames``, never a pass.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

#: Relative change across the window above which a region is not stationary.
#: An analysis parameter, not a physical constant: the default is the size of
#: the far-field deficit the short transient produced in cases/cai2012, so a
#: repeat of that would fail by a clear margin.
DEFAULT_TOLERANCE = 0.02

#: Fewest sampled frames either test is attempted on.
MIN_FRAMES = 4

#: The trend test must resolve a change of this many tolerances.
#:
#: A noisy series has a wide confidence band, and "no trend outside the band"
#: then says nothing: a 20% band passes a 10% drift. When twice the standard
#: error of the fitted change exceeds ``RESOLUTION_FACTOR * tolerance`` the
#: region is ``inconclusive`` -- the data cannot tell -- never ``stable``.
RESOLUTION_FACTOR = 2.0

#: ``(name, kind, (lo, hi) in D or degrees)``. ``tube`` bands are X/D along the
#: centreline tube; ``shell`` sectors are polar angle on r/D in [9, 10].
REGIONS = (
    ("centreline X/D 0-1", "tube", (0.0, 1.0)),
    ("centreline X/D 2-3", "tube", (2.0, 3.0)),
    ("centreline X/D 4-6", "tube", (4.0, 6.0)),
    ("centreline X/D 7-10", "tube", (7.0, 10.0)),
    ("shell r/D 9-10, theta 0-15 deg", "shell", (0.0, 15.0)),
    ("shell r/D 9-10, theta 15-45 deg", "shell", (15.0, 45.0)),
    ("shell r/D 9-10, theta 45-90 deg", "shell", (45.0, 90.0)),
)

#: Reported, never part of the verdict.
CONTEXT_REGION = "domain total (context only)"


def sampled_frames(case_dir: Path, average_start_s: float) -> list[str]:
    """Time directories at or after the averaging start that hold ``rhoNMean``.

    The rule ``cases/cai2012-health/baseCase/Allaudit.py`` applies, with the
    same one-part-in-1e9 tolerance: the write schedule is counted in steps and
    the first sampled frame can land a rounding below the start.
    """
    times = []
    for entry in Path(case_dir).iterdir():
        if not entry.is_dir():
            continue
        try:
            value = float(entry.name)
        except ValueError:
            continue
        if (value > 0.0 and value >= float(average_start_s) * (1.0 - 1.0e-9)
                and (entry / "rhoNMean").is_file()):
            times.append((value, entry.name))
    return [name for _, name in sorted(times)]


def geometry_time(case_dir: Path, times: list[str]) -> str | None:
    """The newest of ``times`` holding cell centres ``C``, else any that does."""
    for name in reversed(times):
        if (Path(case_dir) / name / "C").is_file():
            return name
    for entry in sorted(Path(case_dir).iterdir()):
        if entry.is_dir() and (entry / "C").is_file():
            return entry.name
    return None


def region_masks(centres, cfg, geom) -> dict:
    """Boolean cell selections for every region, plus the context region."""
    d = float(geom.diameter_m)
    centres = np.asarray(centres, dtype=np.float64)
    x = centres[:, 0]
    r = np.hypot(centres[:, 1], centres[:, 2])
    big_r = np.hypot(x, r)
    theta = np.degrees(np.arctan2(r, x))
    tube = r <= float(cfg.post.centerline_radius_over_D) * d
    shell = (big_r >= 9.0 * d) & (big_r <= 10.0 * d) & (x > 0.0)

    masks = {}
    for name, kind, (lo, hi) in REGIONS:
        if kind == "tube":
            masks[name] = tube & (x > lo * d) & (x <= hi * d)
        else:
            masks[name] = shell & (theta >= lo) & (theta < hi if hi < 90.0
                                                   else theta <= hi)
    masks[CONTEXT_REGION] = np.ones(x.size, dtype=bool)
    return masks


def regional_densities(number_density, volumes, masks, n0: float) -> dict:
    """Volume-weighted ``n/n0`` in each region: molecules over volume."""
    number_density = np.asarray(number_density, dtype=np.float64)
    volumes = np.asarray(volumes, dtype=np.float64)
    out = {}
    for name, mask in masks.items():
        volume = float(volumes[mask].sum())
        out[name] = (float((number_density[mask] * volumes[mask]).sum()) / volume
                     / float(n0)) if volume > 0.0 else float("nan")
    return out


def linear_trend(times, values) -> dict:
    """Least-squares line through ``values(times)``, as a relative change.

    Returns:
        ``{"relative_change", "standard_error", "mean"}``: the change the fitted
        slope implies across the whole window, over the mean, and the standard
        error of that number. NaN with fewer than three points.
    """
    t = np.asarray(times, dtype=np.float64)
    v = np.asarray(values, dtype=np.float64)
    keep = np.isfinite(t) & np.isfinite(v)
    t, v = t[keep], v[keep]
    if t.size < 3 or not np.ptp(t) > 0.0:
        return {"relative_change": float("nan"), "standard_error": float("nan"),
                "mean": float(np.mean(v)) if v.size else float("nan")}
    mean = float(v.mean())
    tc = t - t.mean()
    slope = float((tc * (v - mean)).sum() / (tc * tc).sum())
    residual = v - (mean + slope * tc)
    dof = t.size - 2
    slope_se = math.sqrt(float((residual ** 2).sum()) / dof / float((tc * tc).sum()))
    span = float(np.ptp(t))
    scale = abs(mean) if mean != 0.0 else float("nan")
    return {"relative_change": slope * span / scale,
            "standard_error": slope_se * span / scale,
            "mean": mean}


def assess(frames: list, *, tolerance: float = DEFAULT_TOLERANCE,
           min_frames: int = MIN_FRAMES) -> dict:
    """The verdict over one run's sampled frames.

    Args:
        frames: one dict per sampled frame, ascending in time, each
            ``{"time": t, "mean": {region: n/n0}, "instant": {region: n/n0}}``.
            ``instant`` may be absent for a frame that wrote no ``rhoN``.
        tolerance: see :data:`DEFAULT_TOLERANCE`.
        min_frames: see :data:`MIN_FRAMES`.

    Returns:
        ``{"verdict", "tolerance", "n_frames", "regions": {name: {...}}}``.
        Each region is ``stable``, ``drifting`` (its running mean is still
        moving, or a resolvable trend exceeds the tolerance) or
        ``inconclusive`` (too noisy to resolve one -- see
        :data:`RESOLUTION_FACTOR`). ``verdict`` is ``stable-window-demonstrated``
        only if every judged region is stable; ``not-demonstrated`` if any
        drifts; ``inconclusive`` otherwise; ``insufficient-frames`` below
        ``min_frames``.
    """
    frames = sorted(frames, key=lambda f: float(f["time"]))
    regions = {}
    names = [name for name, _, _ in REGIONS] + [CONTEXT_REGION]
    if len(frames) < min_frames:
        return {"verdict": "insufficient-frames", "tolerance": float(tolerance),
                "n_frames": len(frames), "min_frames": int(min_frames),
                "regions": regions}

    times = [float(f["time"]) for f in frames]
    middle = times[0] + 0.5 * (times[-1] - times[0])
    half_index = int(np.argmin([abs(t - middle) for t in times]))

    for name in names:
        means = [f["mean"].get(name, float("nan")) for f in frames]
        instants = [(f.get("instant") or {}).get(name, float("nan"))
                    for f in frames]
        last, halfway = means[-1], means[half_index]
        drift = (abs(last - halfway) / abs(last)
                 if np.isfinite(last) and last != 0.0 and np.isfinite(halfway)
                 else float("nan"))
        trend = linear_trend(times, instants)
        band = 2.0 * trend["standard_error"]
        resolvable = np.isfinite(band) and band <= RESOLUTION_FACTOR * float(tolerance)
        stationary = (np.isfinite(trend["relative_change"])
                      and abs(trend["relative_change"]) <= max(float(tolerance),
                                                               band))
        settled = np.isfinite(drift) and drift <= float(tolerance)
        if not settled or (resolvable and not stationary):
            status = "drifting"
        elif not resolvable:
            status = "inconclusive"
        else:
            status = "stable"
        regions[name] = {
            "running_mean_final": float(last),
            "running_mean_halfway": float(halfway),
            "running_mean_drift": float(drift),
            "running_mean_settled": bool(settled),
            "instant_mean": trend["mean"],
            "instant_relative_change": trend["relative_change"],
            "instant_trend_standard_error": trend["standard_error"],
            "instant_trend_resolvable": bool(resolvable),
            "instant_stationary": bool(stationary),
            "status": status,
            "stable": status == "stable",
            "in_verdict": name != CONTEXT_REGION,
        }

    judged = [r for r in regions.values() if r["in_verdict"]]
    if any(r["status"] == "drifting" for r in judged):
        verdict = "not-demonstrated"
    elif judged and all(r["status"] == "stable" for r in judged):
        verdict = "stable-window-demonstrated"
    else:
        verdict = "inconclusive"
    return {"verdict": verdict, "tolerance": float(tolerance),
            "resolution_factor": RESOLUTION_FACTOR,
            "n_frames": len(frames), "halfway_time": times[half_index],
            "regions": regions,
            "note": ("The domain total is context only: a flat parcel count is "
                     "not evidence that the far field has converged. "
                     "'inconclusive' means the instantaneous data are too noisy "
                     "to resolve a drift of the tolerance -- more parcels or "
                     "more frames, not a pass.")}
