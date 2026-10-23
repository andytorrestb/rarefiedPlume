r"""Figures for the extension's comparisons, on fixed normalised colour scales.

Every case of the extension -- every Kn and every particle level -- is drawn on
the SAME colour range for the same quantity (:data:`PLANE_FIELDS`). Normalised
quantities are comparable across Kn by construction, and a figure auto-scaled
to its own data could not be held against another: images on different scales
differ in their colour maps, not in their physics.

Figures of a case that is unsuitable for production carry that verdict in the
title, in capitals. matplotlib is optional (``pip install -e ".[plots]"``); every
function returns ``[]`` without it and the numbers are already in the CSVs.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from plumetools.cai2012ext import comparison

#: ``key: (label, log, vmin, vmax, difference kind, difference half-range)``.
#: The difference panel is DSMC/analytical - 1 for "relative" and
#: DSMC - analytical for "absolute".
PLANE_FIELDS = {
    "n_over_n0": (r"$n/n_0$", True, 1.0e-4, 1.0, "relative", 0.5),
    "ux": (r"$u_x\sqrt{\beta_0}$  (axial)", False, 0.0, 2.6, "absolute", 0.5),
    "ur": (r"$u_r\sqrt{\beta_0}$  (cylindrical radial)", False, 0.0, 2.6,
           "absolute", 0.5),
    "T_over_T0": (r"$T/T_0$  (local translational)", False, 0.0, 1.1,
                  "absolute", 0.5),
    "p_over_p0": (r"$p/p_0 = (n/n_0)(T/T_0)$", True, 1.0e-5, 1.0, "relative", 0.5),
}

#: File name per plane quantity.
PLANE_FILES = {
    "n_over_n0": "plane-density.png",
    "ux": "plane-velocity-axial.png",
    "ur": "plane-velocity-radial.png",
    "T_over_T0": "plane-temperature.png",
    "p_over_p0": "plane-pressure.png",
}


def status_banner(production_ready: bool) -> str:
    """A title suffix stating the case's verdict."""
    return ("" if production_ready else
            "\nNON-PRODUCTION -- UNSUITABLE FOR PRODUCTION (see comparison.yaml)")


def _pyplot():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    return plt


def plot_centerline_pressure(results_dir: Path, profile: dict, label: str,
                             banner: str = "") -> list[Path]:
    """``p/p0`` along the axis: DSMC against the analytical solution."""
    plt = _pyplot()
    if plt is None:
        return []
    figure, axes = plt.subplots(figsize=(6.0, 4.2))
    axes.plot(profile["x_over_D"], profile["p_over_p0_analytical_axis"], "-",
              color="black", linewidth=1.6, label="Cai collisionless, on axis")
    axes.plot(profile["x_over_D"], profile["p_over_p0"], "o", markersize=3.0,
              color="tab:red", label=f"dsmcFoam, {label}")
    axes.set_yscale("log")
    axes.set_xlabel("$X/D$")
    axes.set_ylabel(r"$p/p_0 = (n/n_0)(T/T_0)$")
    axes.grid(True, which="both", alpha=0.3)
    axes.legend()
    axes.set_title(f"Centreline pressure -- {label}{banner}", fontsize=9)
    figure.tight_layout()
    path = Path(results_dir) / "centerline-pressure.png"
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return [path]


def plot_angular(results_dir: Path, profile: dict, label: str,
                 banner: str = "") -> list[Path]:
    """The r/D = 10 density and mass-flux profiles."""
    plt = _pyplot()
    if plt is None:
        return []
    results_dir = Path(results_dir)
    theta = profile["theta_deg"]
    radius = profile["shell"]["radius_over_D"]
    written = []

    figure, axes = plt.subplots(figsize=(6.4, 4.4))
    axes.plot(theta, profile["n_over_n0_exact"], "-", color="black",
              linewidth=1.4, label=f"Cai collisionless, on r/D = {radius:g}")
    axes.plot(theta, profile["n_over_n0_analytical"], "--", color="grey",
              linewidth=1.0, label="Cai collisionless, same shell cells")
    axes.plot(theta, profile["n_over_n0"], "o", markersize=3.5,
              color="tab:red", label=f"dsmcFoam, {label}")
    axes.set_yscale("log")
    axes.set_xlabel(r"polar angle $\theta$ from the axis [deg]")
    axes.set_ylabel(r"$n/n_0$")
    axes.grid(True, which="both", alpha=0.3)
    axes.legend(fontsize=8)
    axes.set_title(f"Angular density at r/D = {radius:g} -- {label}{banner}",
                   fontsize=9)
    figure.tight_layout()
    path = results_dir / "angular-density.png"
    figure.savefig(path, dpi=150)
    plt.close(figure)
    written.append(path)

    figure, axes = plt.subplots(figsize=(6.4, 4.4))
    axes.plot(theta, profile["F_exact"], "-", color="black", linewidth=1.4,
              label="Cai collisionless, on the sphere")
    axes.plot(theta, profile["F_analytical"], "--", color="grey", linewidth=1.0,
              label="Cai collisionless, same shell cells")
    axes.plot(theta, profile["F"], "o", markersize=3.5, color="tab:red",
              label=f"dsmcFoam, {label}")
    axes.set_xlabel(r"polar angle $\theta$ from the axis [deg]")
    axes.set_ylabel(r"$F = n u_R R^2 / (\Gamma_0 \pi R_0^2)$  [1/sr]")
    axes.grid(True, alpha=0.3)
    integrals = profile["integrals"]
    partial = ""
    if integrals["bins_without_cells"]:
        partial = (f" ({integrals['F_exact_occupied_bins']:.4f} over the "
                   f"covered bins; {integrals['bins_without_cells']} empty)")
    axes.set_title(
        f"Mass flux through r/D = {radius:g} -- {label}{banner}\n"
        rf"$\int F\,2\pi\sin\theta\,d\theta$: DSMC {integrals['F']:.4f}, "
        f"collisionless {integrals['F_exact']:.4f}{partial}", fontsize=9)
    axes.legend(fontsize=8)
    figure.tight_layout()
    path = results_dir / "angular-massflux.png"
    figure.savefig(path, dpi=150)
    plt.close(figure)
    written.append(path)
    return written


def plot_plane(results_dir: Path, plane: dict, grid: dict, label: str,
               banner: str = "") -> list[Path]:
    """DSMC | analytical | difference, in the centre plane, per quantity.

    One cell layer is drawn (the one on the positive side of the plane): the
    two straddling layers are mirror images and drawing both superimposes
    duplicate points. Both enter the metrics.
    """
    plt = _pyplot()
    if plt is None:
        return []
    from matplotlib.colors import LogNorm, Normalize
    from matplotlib.tri import Triangulation

    measured, exact = plane["measured"], plane["exact"]
    normal_index = 1 if plane["lateral_axis"] == "Z" else 2
    layer = measured["centres"][:, normal_index] >= 0.0
    xs = plane["x_over_D"][layer]
    ls = plane["lateral_over_D"][layer]
    if xs.size < 3:
        return []
    triangles = Triangulation(xs, ls)

    exact_n = exact["n_over_n0"][layer]
    usable = (np.isfinite(exact_n)
              & (exact_n >= comparison.DENSITY_FLOOR_OVER_N0)
              & ~plane["lip"][layer] & measured["occupied"][layer])

    written = []
    for key, (title, log, vmin, vmax, kind, span) in PLANE_FIELDS.items():
        figure, axes = plt.subplots(1, 3, figsize=(15.0, 4.6), sharey=True)
        norm = LogNorm(vmin=vmin, vmax=vmax) if log else Normalize(vmin, vmax)

        values = np.asarray(measured[key][layer], dtype=float)
        shown = np.where(values > 0.0, values, np.nan) if log else values
        mesh = axes[0].tripcolor(triangles, np.nan_to_num(
            shown, nan=vmin if not log else vmin * 1.0e-3), norm=norm,
            shading="gouraud", cmap="viridis")
        axes[0].set_title(f"dsmcFoam -- {label}", fontsize=9)

        reference = np.asarray(grid[key], dtype=float)
        reference = np.where(reference > 0.0, reference, np.nan) if log else reference
        axes[1].pcolormesh(grid["x_over_D"], grid["lateral_over_D"], reference,
                           norm=norm, shading="auto", cmap="viridis")
        axes[1].set_title("Cai collisionless (analytical)", fontsize=9)
        figure.colorbar(mesh, ax=axes[:2], label=title, shrink=0.85)

        exact_values = np.asarray(exact[key][layer], dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            if kind == "relative":
                difference = values / exact_values - 1.0
                dlabel = "DSMC / analytical - 1"
            else:
                difference = values - exact_values
                dlabel = "DSMC - analytical"
        difference = np.where(usable, difference, 0.0)
        diff_mesh = axes[2].tripcolor(triangles, difference, cmap="RdBu_r",
                                      vmin=-span, vmax=span, shading="gouraud")
        axes[2].set_title(f"{dlabel} (0 outside the compared plume)", fontsize=9)
        figure.colorbar(diff_mesh, ax=axes[2], shrink=0.85)

        for axis in axes:
            axis.set_xlabel("$X/D$")
            axis.set_aspect("equal", adjustable="box")
        axes[0].set_ylabel(f"${plane['lateral_axis']}/D$")
        figure.suptitle(f"{title} in the centre plane{banner}", fontsize=10)
        path = Path(results_dir) / PLANE_FILES[key]
        figure.savefig(path, dpi=130)
        plt.close(figure)
        written.append(path)
    return written


def plot_stationarity(results_dir: Path, frames: list, assessment: dict,
                      transit_s: float, start_s: float, label: str,
                      banner: str = "") -> list[Path]:
    """Regional n/n0 per frame: running mean and instantaneous, normalised."""
    plt = _pyplot()
    if plt is None or not frames:
        return []
    times = np.array([float(f["time"]) for f in frames])
    elapsed = (times - float(start_s)) / float(transit_s)
    names = list(assessment.get("regions") or {})
    figure, axes = plt.subplots(len(names), 1, figsize=(7.0, 1.7 * len(names)),
                                sharex=True, squeeze=False)
    for axis, name in zip(axes[:, 0], names):
        mean = np.array([f["mean"].get(name, np.nan) for f in frames])
        instant = np.array([(f.get("instant") or {}).get(name, np.nan)
                            for f in frames])
        reference = mean[-1] if np.isfinite(mean[-1]) and mean[-1] else 1.0
        axis.plot(elapsed, instant / reference, ".", color="tab:grey",
                  markersize=3, label="instantaneous")
        axis.plot(elapsed, mean / reference, "-", color="tab:red",
                  label="running mean")
        region = assessment["regions"][name]
        verdict = region["status"] if region["in_verdict"] else (
            f"{region['status']}, context only")
        axis.set_ylabel("n / final", fontsize=7)
        axis.set_title(f"{name}: {verdict}", fontsize=8)
        axis.grid(True, alpha=0.3)
    axes[-1, 0].set_xlabel("domain transits since the averaging start")
    axes[0, 0].legend(fontsize=7)
    figure.suptitle(f"Stationarity -- {label}: {assessment['verdict']}{banner}",
                    fontsize=10)
    figure.tight_layout()
    path = Path(results_dir) / "stationarity.png"
    figure.savefig(path, dpi=130)
    plt.close(figure)
    return [path]
