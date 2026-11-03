r"""Post-processing: DSMC moments, profiles, azimuthal statistics, comparisons.

Distinct quantities, never collapsed into one "validation error":

* the prescribed source state (:class:`.physics.SourceState`);
* the CFD precursor solution (cells);
* the mapped interface state (``interface_fields.npz``);
* the expected Maxwellian injection (half-range integrals);
* the measured DSMC injection and return (log counts, crossings);
* the sampled downstream DSMC solution (``fieldAverage`` moments; tensor
  moments from parcel snapshots);
* the independent collisionless reference (:mod:`.collisionless`).

Profiles are taken from the actual 3-D cells. A "centreline" value is the
volume-weighted mean of the cells whose centres lie within one core cell of the
axis -- a tube of finite radius -- and the reference is averaged over **the same
cell centres** so the comparison carries no tube-averaging bias. Azimuthal
averaging is a reporting operation only: :func:`azimuthal_statistics` returns the
unaveraged spread as well.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from plumetools.foamio.volfield import read_internal_field
from plumetools.rarefaction.core import translational_temperature


def dsmc_fieldaverage(time_dir: Path, gas) -> dict:
    """Cell moments from fieldAverage means; temperature from the moments."""
    d = Path(time_dir)
    need = ("rhoNMean", "rhoMMean", "momentumMean", "linearKEMean", "dsmcRhoNMean")
    missing = [n for n in need if not (d / n).is_file()]
    if missing:
        raise FileNotFoundError(f"{d}: no {missing}; the time is before fieldAverage's "
                                "timeStart or the run did not write it")
    n = read_internal_field(d / "rhoNMean")
    rho = read_internal_field(d / "rhoMMean")
    mom = read_internal_field(d / "momentumMean").reshape(-1, 3)
    ke = read_internal_field(d / "linearKEMean")
    parcels = read_internal_field(d / "dsmcRhoNMean")
    with np.errstate(divide="ignore", invalid="ignore"):
        u = mom / rho[:, None]
    t = translational_temperature(rho, mom, ke, gas.gas_constant)
    out = {"n": n, "rho": rho, "U": u, "T": t, "parcels_per_step": parcels}
    if (d / "rhoNPrime2Mean").is_file():
        var = read_internal_field(d / "rhoNPrime2Mean")
        with np.errstate(divide="ignore", invalid="ignore"):
            out["n_cv_per_step"] = np.sqrt(np.clip(var, 0, None)) / n
    return out


def axis_cells(centres: np.ndarray, half_width: float) -> np.ndarray:
    return (np.abs(centres[:, 1]) <= half_width) & (np.abs(centres[:, 2]) <= half_width)


def line_profile(centres, volumes, select, fields: dict, axis: int = 0, decimals: int = 9) -> dict:
    """Volume-weighted means of ``fields`` over ``select`` cells, grouped by coordinate."""
    x = np.round(centres[select, axis], decimals)
    keys, inv = np.unique(x, return_inverse=True)
    w = volumes[select]
    out = {"x": keys}
    for name, values in fields.items():
        v = np.asarray(values)[select]
        if v.ndim == 1:
            ok = np.isfinite(v)
            num = np.bincount(inv, np.where(ok, v * w, 0.0), minlength=keys.size)
            den = np.bincount(inv, np.where(ok, w, 0.0), minlength=keys.size)
            with np.errstate(divide="ignore", invalid="ignore"):
                out[name] = num / den
        else:
            comps = []
            for k in range(v.shape[1]):
                ok = np.isfinite(v[:, k])
                num = np.bincount(inv, np.where(ok, v[:, k] * w, 0.0), minlength=keys.size)
                den = np.bincount(inv, np.where(ok, w, 0.0), minlength=keys.size)
                with np.errstate(divide="ignore", invalid="ignore"):
                    comps.append(num / den)
            out[name] = np.stack(comps, axis=1)
    return out


def azimuthal_statistics(centres, values, x_bin: float, r_bin: float, min_cells: int = 8,
                         mask=None) -> dict:
    """Mean and spread over azimuth in (x, r) bins; the spread is NOT averaged away.

    Returns per bin: x, r, mean, relative standard deviation, count, and the
    largest relative standard deviation over bins with ``min_cells`` or more.
    """
    c = centres
    v = np.asarray(values, float)
    ok = np.isfinite(v) & (True if mask is None else mask)
    ix = np.floor(c[:, 0] / x_bin).astype(np.int64)
    ir = np.floor(np.hypot(c[:, 1], c[:, 2]) / r_bin).astype(np.int64)
    key = ix * 100000 + ir
    keys, inv = np.unique(key[ok], return_inverse=True)
    vv = v[ok]
    cnt = np.bincount(inv)
    mean = np.bincount(inv, vv) / cnt
    var = np.bincount(inv, vv * vv) / cnt - mean ** 2
    rel = np.sqrt(np.clip(var, 0, None)) / np.abs(mean)
    good = cnt >= min_cells
    return {"x": (keys // 100000 + 0.5) * x_bin, "r": (keys % 100000 + 0.5) * r_bin,
            "mean": mean, "rel_std": rel, "count": cnt,
            "max_rel_std": float(np.nanmax(rel[good])) if good.any() else None,
            "median_rel_std": float(np.nanmedian(rel[good])) if good.any() else None}


def compare(a: np.ndarray, b: np.ndarray, ref_scale: float | None = None) -> dict:
    """Differences of ``a`` from ``b``: relative where ``b`` is a positive scale."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if not ok.any():
        return {"points": 0}
    if ref_scale is None:
        rel = (a[ok] - b[ok]) / np.abs(b[ok])
    else:
        rel = (a[ok] - b[ok]) / ref_scale
    return {"points": int(ok.sum()), "max_abs_rel": float(np.max(np.abs(rel))),
            "mean_rel": float(np.mean(rel)), "rms_rel": float(np.sqrt(np.mean(rel ** 2)))}


def write_csv(path: Path, columns: dict) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    names, cols = [], []
    for k, v in columns.items():
        v = np.asarray(v)
        if v.ndim == 2:
            for j in range(v.shape[1]):
                names.append(f"{k}_{'xyz'[j]}")
                cols.append(v[:, j])
        else:
            names.append(k)
            cols.append(v)
    arr = np.column_stack(cols)
    np.savetxt(path, arr, delimiter=",", header=",".join(names), comments="")
    return path
