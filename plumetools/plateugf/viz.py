"""Flow slices on the two symmetry planes of a plate case (phase B, B4).

The quarter domain's symmetry planes y = 0 and z = 0 are the plume's two
mid-planes. The cells touching each (centre within half a cell of it) are drawn
with ``tripcolor`` on their centres -- no resampling -- for the averaged number
density over the exit value, the speed and the translational temperature, on
FIXED ranges so every case of the sweep shares one colour scale. The plate is
drawn as an outline. Needs ``C`` and ``V`` in the time directory
(``postProcess -func writeCellCentres`` / ``writeCellVolumes``; ``Allpost``
writes them).

The density comparison: the collisionless density of the same geometry
(:mod:`plumetools.plateugf.collisionless`: Cai's plume, shadowed by the plate,
plus the gas the plate re-emits) is drawn on the same scale (``ncai``), the
simulated over the collisionless density beside it (``nratio``), and
``results/cai_density.yaml`` holds that ratio for the gas in four regions. At
large Kn_e it should sit near 1; the departures are the collisional effects.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from plumetools.cai2012ugf.settings import load_settings
from plumetools.foamio.volfield import read_internal_field

#: (title, unit, range, log, colour map) per drawn field.
FIELDS = {
    "n": ("number density / exit value", "n/n0", (1.0e-7, 1.0), True, "viridis"),
    "ncai": ("collisionless density (Cai + plate gas) / exit value", "n/n0", (1.0e-7, 1.0), True, "viridis"),
    "nratio": ("simulated / collisionless density", "n_sim / n_collisionless", (0.1, 10.0), True, "RdBu_r"),
    "U": ("speed", "|U| [m/s]", (0.0, 1500.0), False, "viridis"),
    "T": ("translational temperature", "T [K]", (0.0, 600.0), False, "viridis"),
}
_KEYS = {"n": "n_over_n0", "ncai": "ncai_over_n0", "nratio": "nratio", "U": "U", "T": "T"}
#: Below this collisionless n/n0 the ratio is not drawn (shown as 1).
RATIO_FLOOR = 1.0e-7


def regions(x: np.ndarray, standoff_m: float) -> dict:
    """Masks along the axis, as fractions of the stand-off ``s``."""
    s = float(standoff_m)
    return {
        "near_nozzle": (x > 0.0) & (x < 0.25 * s),
        "mid_plume": (x >= 0.25 * s) & (x < 0.75 * s),
        "before_plate": (x >= 0.75 * s) & (x < s),
        "beside_behind_plate": x >= s,
        "all": np.ones(x.shape, bool),
    }


def density_summary(centres, volumes, n_sim, n_cai, standoff_m: float) -> dict:
    """Simulated over collisionless gas held in each region: sum(n V) / sum(n V)."""
    out = {}
    for name, mask in regions(np.asarray(centres)[:, 0], standoff_m).items():
        sim = float((n_sim[mask] * volumes[mask]).sum())
        cai = float((n_cai[mask] * volumes[mask]).sum())
        out[name] = {"cells": int(mask.sum()), "molecules_sim": sim, "molecules_cai": cai,
                     "ratio": sim / cai if cai > 0 else None}
    return out


def touching(centres: np.ndarray, volumes: np.ndarray, axis: int) -> np.ndarray:
    """Cells whose centre lies within half a cell of the plane ``axis = 0``."""
    return np.abs(centres[:, axis]) <= 0.51 * np.cbrt(volumes)


def read(case_dir: Path, time: str) -> dict:
    """Cell centres, volumes and the drawn fields at ``time``."""
    from plumetools.nozzleplume import exitstate
    from plumetools.nozzleplume.config import load_case_config
    case_dir = Path(case_dir)
    t = case_dir / time
    field = load_settings(case_dir).fields.mean_field
    n0 = exitstate.from_config(load_case_config(case_dir)).number_density_per_m3
    from plumetools.plateugf.collisionless import Model
    rho = np.asarray(read_internal_field(t / f"rhoN_{field}")).ravel()
    speed = np.linalg.norm(np.asarray(read_internal_field(t / f"UMean_{field}")).reshape(-1, 3), axis=1)
    centres = np.asarray(read_internal_field(t / "C")).reshape(-1, 3)
    cai = Model(load_case_config(case_dir)).density(centres)
    shown = cai > RATIO_FLOOR * n0
    return {
        "C": centres,
        "V": np.asarray(read_internal_field(t / "V")).ravel(),
        "n_sim": rho, "n_cai": cai,
        "n_over_n0": rho / n0,
        "ncai_over_n0": cai / n0,
        "nratio": np.where(shown, rho / np.where(shown, cai, 1.0), 1.0),
        "U": np.where(rho > 0, speed, 0.0),
        "T": np.where(rho > 0, np.asarray(read_internal_field(t / f"translationalT_{field}")).ravel(), 0.0),
    }


def draw(data: dict, out_dir: Path, *, label: str, plate: tuple) -> list[Path]:
    """One image per ``FIELDS`` entry held in ``data``, on y = 0 and on z = 0.
    ``plate = (x0, x1, y_half, z_half)``."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize
    from matplotlib.patches import Rectangle

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    centres, volumes = data["C"], data["V"]
    x0, x1, y_half, z_half = plate
    written = []
    for plane, axis, other, other_name, extent in (("y0", 1, 2, "z", z_half), ("z0", 2, 1, "y", y_half)):
        cells = touching(centres, volumes, axis)
        if int(cells.sum()) < 3:
            continue
        x, w = centres[cells, 0], centres[cells, other]
        for key, (title, unit, (lo, hi), log, cmap) in FIELDS.items():
            if _KEYS[key] not in data:
                continue
            values = np.asarray(data[_KEYS[key]])[cells]
            norm = LogNorm(lo, hi) if log else Normalize(lo, hi)
            values = np.clip(np.where(values > 0, values, lo) if log else values, lo, hi)
            span_x, span_w = max(np.ptp(x), 1e-9), max(np.ptp(w), 1e-9)
            figure, axes = plt.subplots(figsize=(10.0, max(3.0, min(9.0, 10.0 * span_w / span_x)) + 1.0))
            image = axes.tripcolor(x, w, values, norm=norm, cmap=cmap, shading="gouraud")
            figure.colorbar(image, ax=axes, label=unit)
            axes.add_patch(Rectangle((x0, 0.0), x1 - x0, extent, fill=True,
                                     facecolor="0.85", edgecolor="k", lw=1.0))
            axes.set_xlabel("x [m] (plume axis)")
            axes.set_ylabel(f"{other_name} [m]")
            axes.set_aspect("equal", adjustable="box")
            axes.set_title(f"{title} on {plane[0]} = 0 -- {label}", fontsize=10)
            figure.tight_layout()
            path = out_dir / f"slice-{key}-{plane}.png"
            figure.savefig(path, dpi=130)
            plt.close(figure)
            written.append(path)
    return written


def slices(case_dir: Path, time: str | None = None) -> list[Path]:
    """Draw the slices of ``case_dir`` at ``time`` (default: newest with C and V) into results/."""
    from plumetools.nozzleplume.config import load_case_config
    from plumetools.plateugf.post import PostError, _times
    case_dir = Path(case_dir)
    if time is None:
        held = [name for _, name in _times(case_dir) if (case_dir / name / "C").is_file()]
        if not held:
            raise PostError(f"no time in {case_dir} has cell centres C: run "
                            f"'postProcess -func writeCellCentres' (Allpost does)")
        time = held[-1]
    cfg = load_case_config(case_dir)
    p = cfg.plate
    plate = (p.standoff_m, p.standoff_m + p.thickness_m, 0.5 * p.size_y_m, 0.5 * p.size_z_m)
    label = f"{cfg.meta.get('case_name') or case_dir.resolve().name}, t = {time} s"
    data = read(case_dir, time)
    written = draw(data, case_dir / "results", label=label, plate=plate)
    import yaml
    summary = {"time": time, "standoff_m": float(p.standoff_m),
               "regions_note": "x as fractions of the stand-off s: near_nozzle 0-0.25 s, "
                               "mid_plume 0.25-0.75 s, before_plate 0.75 s-s, then x >= s",
               "regions": density_summary(data["C"], data["V"], data["n_sim"], data["n_cai"],
                                          p.standoff_m)}
    path = case_dir / "results" / "cai_density.yaml"
    path.write_text("# Simulated vs collisionless (Cai's plume + the plate's re-emitted gas)\n"
                    "# molecules per region: plumetools.plateugf.viz. Reported, not judged.\n"
                    + yaml.safe_dump(summary, sort_keys=False), encoding="utf-8", newline="\n")
    return [*written, path]
