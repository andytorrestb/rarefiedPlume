r"""Plate loads from a uniGasFoam plate case, and the collisionless Cai comparison.

uniGasFoam's diffuse wall writes the time-averaged loads on the plate patch
(``fD_<field>``: the momentum the gas hands the wall per area and time;
``surfaceHeatTransfer_<field>``: the heat flux into it), with ``<field>`` the
running mean of ``unigas.yaml`` (``Ar``). They are read straight into
nozzleplume's :class:`~plumetools.nozzleplume.post.PlateLoads`, whose pressure
(``fD . n``), shear (the tangential part) and totals (scaled from the quarter
domain to the whole plate) are r-4d's estimators, unchanged.

The comparison is on the plate's UPSTREAM face only -- the face Cai's
collisionless solution (:mod:`plumetools.plateugf.caiplate`) describes; edge
and lee faces are left out. Differences are reported, never judged.

Maps are drawn normalised by the exit static pressure ``p0 = n0 k T0`` on fixed
logarithmic ranges, so every case of the sweep shares one colour scale.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import yaml

from plumetools.cai2012ugf.settings import load_settings
from plumetools.foamio.errors import FoamDataError
from plumetools.foamio.volfield import read_patch_field
from plumetools.mesh.patches import read_patches
from plumetools.nozzleplume import exitstate, geometry
from plumetools.nozzleplume.config import load_case_config
from plumetools.nozzleplume.constants import BOLTZMANN_J_PER_K
from plumetools.nozzleplume.post import (PlateLoads, plot_results, write_csv,
                                         write_plate_results, write_yaml)
from plumetools.plateugf import caiplate

#: Fixed colour range for loads normalised by p0 (log10 decades).
MAP_RANGE = (1.0e-9, 1.0e-1)
RATIO_RANGE = (0.25, 4.0)


class PostError(RuntimeError):
    """The case cannot be post-processed as asked."""


def _setup(case_dir: Path):
    cfg = load_case_config(Path(case_dir))
    return cfg, geometry.from_config(cfg), exitstate.from_config(cfg)


def _times(case_dir: Path) -> list[tuple[float, str]]:
    out = []
    for entry in Path(case_dir).iterdir():
        try:
            value = float(entry.name)
        except ValueError:
            continue
        if entry.is_dir() and value > 0.0:
            out.append((value, entry.name))
    return sorted(out)


def read_plate(case_dir: Path, time: str | None = None) -> PlateLoads:
    """The averaged plate loads at ``time`` (default: the newest write holding them)."""
    case_dir = Path(case_dir)
    cfg, geom, _ = _setup(case_dir)
    field = load_settings(case_dir).fields.mean_field
    names = (f"fD_{field}", f"surfaceHeatTransfer_{field}")
    if time is None:
        held = [name for _, name in _times(case_dir) if (case_dir / name / names[0]).is_file()]
        if not held:
            raise PostError(
                f"no reconstructed time in {case_dir} holds {names[0]}: run "
                f"reconstructPar (Allpost does), or the run has not written yet")
        time = held[-1]
    plate = cfg.mesh.patch_names["plate"]
    patch = read_patches(case_dir, [plate])[plate]
    values = {}
    for name in names:
        path = case_dir / time / name
        try:
            values[name] = np.asarray(read_patch_field(path, plate), dtype=float)
        except (OSError, FoamDataError) as exc:
            raise PostError(f"{path}: {exc}") from None
    return PlateLoads(
        time=time, centres=patch.centres, areas=patch.areas, normals=patch.normals,
        force_density=values[names[0]].reshape(-1, 3), heat_flux=values[names[1]].ravel(),
        multiplicity=geom.symmetry_multiplicity, symmetry_planes=tuple(geom.symmetry_planes))


def upstream(loads: PlateLoads) -> np.ndarray:
    """Mask of the upstream (+x facing, plume-side) faces."""
    return loads.side == 1


def profiles(loads: PlateLoads, values: dict) -> dict:
    """Values along the plate's two centrelines on the upstream face: the row of
    faces nearest z = 0 (``along_y``) and the column nearest y = 0 (``along_z``)."""
    up = upstream(loads)
    y, z = loads.centres[up, 1], loads.centres[up, 2]
    half = 0.51 * np.sqrt(loads.areas[up])
    out = {}
    for key, across, along in (("along_y", z, y), ("along_z", y, z)):
        # the faces that touch the plane across = 0: centre within half a face
        # of it. Robust to snappy's micrometre noise and to mixed face sizes.
        row = np.flatnonzero(np.abs(across) <= half)
        row = row[np.argsort(along[row])]
        out[key] = {"y_m": y[row], "z_m": z[row],
                    **{name: np.asarray(v)[up][row] for name, v in values.items()}}
    return out


def compare(case_dir: Path, time: str | None = None, loads: PlateLoads | None = None) -> dict:
    """Simulated vs collisionless loads on the upstream face."""
    case_dir = Path(case_dir)
    cfg, geom, es = _setup(case_dir)
    loads = loads or read_plate(case_dir, time)
    up = upstream(loads)
    centres, areas = loads.centres[up], loads.areas[up]
    cai = caiplate.plate_loads(
        centres, radius_m=0.5 * es.diameter_m, number_density=es.number_density_per_m3,
        speed=es.speed_m_per_s, temperature_K=es.temperature_K, mass_kg=es.species.mass_kg,
        wall_temperature_K=float(cfg.plate.temperature_K))
    p_sim, tau_sim, q_sim = loads.pressure[up], loads.shear[up], loads.heat_flux[up]
    mult = loads.multiplicity
    centre = int(np.argmin(np.hypot(centres[:, 1], centres[:, 2]))) if centres.size else 0

    def ratio(a, b):
        return float(a / b) if b else None

    force_sim, force_cai = float((p_sim * areas).sum()) * mult, float((cai.pressure * areas).sum()) * mult
    heat_sim, heat_cai = float((q_sim * areas).sum()) * mult, float((cai.heat_flux * areas).sum()) * mult
    summary = {
        "time": loads.time,
        "faces": int(up.sum()),
        "p0_pa": float(es.number_density_per_m3 * BOLTZMANN_J_PER_K * es.temperature_K),
        "whole_plate": {
            "normal_force_sim_n": force_sim, "normal_force_cai_n": force_cai,
            "normal_force_ratio": ratio(force_sim, force_cai),
            "heat_load_sim_w": heat_sim, "heat_load_cai_w": heat_cai,
            "heat_load_ratio": ratio(heat_sim, heat_cai),
        },
        "centre_face": {
            "y_m": float(centres[centre, 1]), "z_m": float(centres[centre, 2]),
            "p_sim_pa": float(p_sim[centre]), "p_cai_pa": float(cai.pressure[centre]),
            "p_ratio": ratio(p_sim[centre], cai.pressure[centre]),
        },
        "max_shear": {"sim_pa": float(tau_sim.max()), "cai_pa": float(cai.shear_magnitude.max()),
                      "ratio": ratio(tau_sim.max(), cai.shear_magnitude.max())},
    }
    table = {
        "y_m": centres[:, 1], "z_m": centres[:, 2], "area_m2": areas,
        "p_sim_pa": p_sim, "p_cai_pa": cai.pressure,
        "tau_sim_pa": tau_sim, "tau_cai_pa": cai.shear_magnitude,
        "q_sim_w_per_m2": q_sim, "q_cai_w_per_m2": cai.heat_flux,
    }
    return {"faces": int(up.sum()), "table": table, "summary": summary, "cai": cai, "loads": loads}


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #

def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    return plt, LogNorm


def _map(path: Path, y, z, values, title: str, label: str, norm, cmap: str) -> Path:
    plt, _ = _plt()
    span_y, span_z = max(np.ptp(y), 1e-9), max(np.ptp(z), 1e-9)
    figure, axes = plt.subplots(figsize=(8.0, max(2.5, min(8.0, 8.0 * span_z / span_y)) + 1.0))
    image = axes.tripcolor(y, z, np.clip(values, norm.vmin, norm.vmax), norm=norm,
                           cmap=cmap, shading="flat") if len(y) >= 3 else None
    if image is not None:
        figure.colorbar(image, ax=axes, label=label)
    axes.set_xlabel("y [m]")
    axes.set_ylabel("z [m]")
    axes.set_aspect("equal", adjustable="box")
    axes.set_title(title, fontsize=10)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return path


def _figures(results: Path, result: dict, label: str) -> list[Path]:
    try:
        plt, LogNorm = _plt()
    except ImportError:
        return []
    t, s = result["table"], result["summary"]
    p0 = s["p0_pa"]
    y, z = t["y_m"], t["z_m"]
    norm = LogNorm(*MAP_RANGE)
    written = [_map(results / "plate-shear.png", y, z, t["tau_sim_pa"] / p0,
                    f"Plate upstream face: shear / p0\n{label}, t = {s['time']} s",
                    "|tau| / p0", norm, "viridis")]
    # sim / Cai ratio maps
    figure, panels = plt.subplots(2, 1, figsize=(8.0, 7.0))
    rnorm = LogNorm(*RATIO_RANGE)
    for axes, key, name in zip(panels, ("p", "tau"), ("pressure", "shear")):
        sim, cai = t[f"{key}_sim_pa"], t[f"{key}_cai_pa"]
        ratio = np.where(cai > 0, sim / np.where(cai > 0, cai, 1.0), np.nan)
        if len(y) >= 3:
            image = axes.tripcolor(y, z, np.clip(np.nan_to_num(ratio, nan=1.0), *RATIO_RANGE),
                                   norm=rnorm, cmap="RdBu_r", shading="flat")
            figure.colorbar(image, ax=axes, label=f"{name}: simulation / Cai")
        axes.set_aspect("equal", adjustable="box")
        axes.set_xlabel("y [m]")
        axes.set_ylabel("z [m]")
    panels[0].set_title(f"Simulation over collisionless Cai, upstream face\n{label}", fontsize=10)
    figure.tight_layout()
    written.append(results / "cai-ratio-maps.png")
    figure.savefig(written[-1], dpi=150)
    plt.close(figure)
    # centreline profiles, simulation vs Cai
    loads = result["loads"]
    prof = profiles(loads, {"p_sim": loads.pressure, "tau_sim": loads.shear})
    up = upstream(loads)
    cai_p = np.zeros(len(loads.pressure))
    cai_t = np.zeros(len(loads.pressure))
    cai_p[up], cai_t[up] = t["p_cai_pa"], t["tau_cai_pa"]
    cprof = profiles(loads, {"p_cai": cai_p, "tau_cai": cai_t})
    figure, panels = plt.subplots(2, 2, figsize=(10.0, 7.0), sharey="row")
    for col, key, axis_name in ((0, "along_y", "y_m"), (1, "along_z", "z_m")):
        for row, q, name in ((0, "p", "pressure"), (1, "tau", "shear")):
            ax = panels[row][col]
            ax.semilogy(prof[key][axis_name], prof[key][f"{q}_sim"] / p0, "o", ms=3, label="simulation")
            ax.semilogy(cprof[key][axis_name], cprof[key][f"{q}_cai"] / p0, "-", label="collisionless Cai")
            ax.set_xlabel(axis_name.replace("_m", " [m]"))
            ax.set_ylabel(f"{name} / p0")
            ax.grid(True, which="both", alpha=0.3)
    panels[0][0].legend(fontsize=8)
    panels[0][0].set_title("along y (row nearest z = 0)", fontsize=9)
    panels[0][1].set_title("along z (column nearest y = 0)", fontsize=9)
    figure.suptitle(f"Plate upstream face, simulation vs collisionless Cai -- {label}", fontsize=10)
    figure.tight_layout()
    written.append(results / "cai-profiles.png")
    figure.savefig(written[-1], dpi=150)
    plt.close(figure)
    return written


def _window(case_dir: Path, time: str, allow_startup: bool) -> dict:
    """Which leg ``time`` is in, and how long its mean has been accumulating.

    uniGasFoam writes the plate fields at EVERY write. Before
    ``average_start_s`` each is a one-interval average of the startup, not the
    sampled mean, so it is refused unless ``allow_startup``.
    """
    from plumetools.plateugf.workflow import run_window
    try:
        average_start, _ = run_window(case_dir)
    except FoamDataError:
        return {"leg": "unknown", "averaging_window_s": None}
    t = float(time)
    if t < average_start * (1.0 - 1.0e-6):
        if not allow_startup:
            raise PostError(
                f"t = {time} s is before averaging starts ({average_start:.8g} s): its "
                f"plate fields average one write interval of the startup, not the "
                f"sampling window. Wait for the sampling leg's writes, or pass "
                f"--allow-startup to draw it anyway (it is labelled 'startup').")
        return {"leg": "startup", "averaging_window_s": None}
    return {"leg": "sampling", "averaging_window_s": t - average_start}


def run_post(case_dir: Path, time: str | None = None, allow_startup: bool = False) -> list[Path]:
    """Every plate result and figure of one case into ``results/``."""
    case_dir = Path(case_dir)
    cfg, _, es = _setup(case_dir)
    label = str(cfg.meta.get("case_name") or case_dir.resolve().name)
    loads = read_plate(case_dir, time)
    window = _window(case_dir, loads.time, allow_startup)
    if window["leg"] == "startup":
        label += " [STARTUP write: not the sampled mean]"
    result = compare(case_dir, loads=loads)
    result["summary"].update(window)
    results = case_dir / "results"
    results.mkdir(parents=True, exist_ok=True)
    written = list(write_plate_results(results, loads, es))
    prof = profiles(loads, {"pressure_pa": loads.pressure, "shear_pa": loads.shear,
                            "heat_flux_w_per_m2": loads.heat_flux})
    # along_z = 0: the row nearest z = 0 (a profile in y); 1: the column nearest y = 0.
    columns = {"along_z": [], "y_m": [], "z_m": [], "pressure_pa": [], "shear_pa": [],
               "heat_flux_w_per_m2": []}
    for flag, key in ((0, "along_y"), (1, "along_z")):
        data = prof[key]
        columns["along_z"] += [flag] * len(data["y_m"])
        for name in ("y_m", "z_m", "pressure_pa", "shear_pa", "heat_flux_w_per_m2"):
            columns[name] += list(data[name])
    written.append(write_csv(results / "plate_profiles.csv", columns))
    written.append(write_csv(results / "cai_comparison.csv", result["table"]))
    written.append(write_yaml(results / "cai_comparison.yaml", result["summary"],
                              "Plate upstream face: simulation vs collisionless Cai "
                              "(plumetools.plateugf.caiplate). Differences are reported, "
                              "not judged; nothing here is validated."))
    written += plot_results(results, loads, None, label)
    written += _figures(results, result, label)
    return written


def family_table(root: Path) -> list[str]:
    """One line per case under ``root`` that has ``results/cai_comparison.yaml``."""
    lines = [f"{'case':<18}{'t [s]':>12}{'p_c sim':>12}{'p_c Cai':>12}{'ratio':>8}"
             f"{'F sim [N]':>12}{'F Cai [N]':>12}{'ratio':>8}"]
    for path in sorted(Path(root).glob("*/results/cai_comparison.yaml")):
        s = yaml.safe_load(path.read_text(encoding="utf-8"))
        c, w = s["centre_face"], s["whole_plate"]
        lines.append(f"{path.parents[1].name:<18}{float(s['time']):>12.5g}"
                     f"{c['p_sim_pa']:>12.4g}{c['p_cai_pa']:>12.4g}{(c['p_ratio'] or 0):>8.3f}"
                     f"{w['normal_force_sim_n']:>12.4g}{w['normal_force_cai_n']:>12.4g}"
                     f"{(w['normal_force_ratio'] or 0):>8.3f}")
    return lines


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m plumetools.plateugf.post")
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--time")
    args = parser.parse_args(argv)
    for path in run_post(args.case_dir, args.time):
        print(f"    wrote {path}")
    return 0
