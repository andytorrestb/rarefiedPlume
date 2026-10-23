r"""The study table and the Knudsen comparison, from per-case results.

Reads what each case's scripts wrote -- ``metrics.yaml`` (``postProcess.py``),
``diagnostics_summary.yaml`` (``Alldiagnose.py``), ``comparison.yaml``
(``Allcompare.py``), ``occupancy.yaml`` (``Allaudit.py``) and
``stationarity.yaml`` (``Allconverge.py``) -- and the manifest, and writes one
row per case. A case that has not run, or not been post-processed, still
gets its row, with the missing columns empty and ``state`` saying why: a table
that silently lost a case would look complete.

Every row carries the production verdict and the reference it was compared
with. There is no "Cai reported" column: Cai published DSMC results at
Kn = 100, 0.1 and 0.01 only, none of which is in this study.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
import yaml

#: Columns of ``results/study-table.csv``, in order.
COLUMNS = (
    "case", "Kn", "particle_level", "multiplier", "exit_ppc_target",
    "production_status", "state", "cell_over_mfp", "dt_over_tau0", "n_cells",
    "estimated_parcels", "time", "averaged_transits",
    "centerline_density_max_rel_pct", "centerline_density_mean_rel_pct",
    "centerline_velocity_max_rel_pct", "centerline_temperature_max_rel_pct",
    "centerline_pressure_mean_rel_pct",
    "angular_density_mean_rel_pct", "angular_mass_flux_mean_rel_pct",
    "mass_flux_integral_dsmc", "mass_flux_integral_collisionless",
    "plane_density_mean_rel_pct", "plane_axial_velocity_mean_rel_pct",
    "plane_radial_velocity_max_abs", "plane_temperature_mean_rel_pct",
    "plane_pressure_mean_rel_pct",
    "exit_ppc_measured", "plume_ppc_median", "stationarity",
    "diag_valid_volume_fraction", "diag_exit_KnD", "diag_exit_cell_over_mfp",
    "diag_cell_over_mfp_median", "diag_cell_over_mfp_max",
    "diag_volume_fraction_coarser_than_target", "diag_cells_to_resolve",
    "reference",
)

REFERENCE = "Cai collisionless analytical; no published DSMC at this Kn"


def _read(path: Path):
    path = Path(path)
    if not path.is_file():
        return None
    try:
        return yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError:
        return None


def _pct(value):
    if value is None:
        return None
    value = float(value)
    return None if math.isnan(value) else 100.0 * value


def _region(occupancy: dict | None, name: str):
    for region in (occupancy or {}).get("regions") or []:
        if region.get("name") == name:
            return region.get("median")
    return None


def study_row(root: Path, entry: dict) -> dict:
    """One row of the study table."""
    results = Path(root) / entry["path"] / "results"
    metrics = _read(results / "metrics.yaml")
    comparison = _read(results / "comparison.yaml")
    occupancy = _read(results / "occupancy.yaml")
    stationarity = _read(results / "stationarity.yaml")
    diagnostics = _read(results / "diagnostics_summary.yaml") or {}
    feasibility = entry.get("feasibility") or {}
    resolution = (diagnostics.get("resolution") or {}).get("bounding_box_max") or {}
    exit_plane = diagnostics.get("exit_plane") or {}
    budget = (diagnostics.get("planning") or {}).get("mesh_budget") or {}

    if comparison is not None:
        state = "compared"
    elif metrics is not None:
        state = "post-processed (no Allcompare yet)"
    elif (Path(root) / entry["path"] / "log.dsmcFoam").is_file():
        state = "run, not post-processed"
    else:
        state = "not run"

    centerline = (comparison or {}).get("centerline") or {}
    angular = (comparison or {}).get("angular") or {}
    plane = (comparison or {}).get("plane") or {}
    case = (comparison or {}).get("case") or (metrics or {}).get("case") or {}
    return {
        "case": entry["name"],
        "Kn": entry["Kn"],
        "particle_level": entry.get("particle_level"),
        "multiplier": entry.get("numerical_particle_multiplier"),
        "exit_ppc_target": (entry.get("weight") or {}).get(
            "exit_particles_per_cell"),
        "production_status": entry.get("production_status"),
        "state": state,
        "cell_over_mfp": (feasibility.get("resolution") or {}).get(
            "cell_over_mean_free_path"),
        "dt_over_tau0": (feasibility.get("time_step") or {}).get(
            "dt_over_exit_collision_time"),
        "n_cells": entry.get("n_cells"),
        "estimated_parcels": entry.get("estimated_particles"),
        "time": case.get("time"),
        "averaged_transits": case.get("averaged_domain_transits"),
        "centerline_density_max_rel_pct": _pct(centerline.get("density_max_rel_error")),
        "centerline_density_mean_rel_pct": _pct(centerline.get("density_mean_rel_error")),
        "centerline_velocity_max_rel_pct": _pct(centerline.get("velocity_max_rel_error")),
        "centerline_temperature_max_rel_pct": _pct(
            centerline.get("temperature_max_rel_error")),
        "centerline_pressure_mean_rel_pct": _pct(centerline.get("pressure_mean_rel_error")),
        "angular_density_mean_rel_pct": _pct(angular.get("angular_density_mean_rel_error")),
        "angular_mass_flux_mean_rel_pct": _pct(
            angular.get("angular_mass_flux_mean_rel_error")),
        "mass_flux_integral_dsmc": angular.get("angular_mass_flux_integral"),
        "mass_flux_integral_collisionless": angular.get(
            "angular_mass_flux_integral_exact"),
        "plane_density_mean_rel_pct": _pct(plane.get("plane_density_mean_rel_error")),
        "plane_axial_velocity_mean_rel_pct": _pct(
            plane.get("plane_axial_velocity_mean_rel_error")),
        "plane_radial_velocity_max_abs": plane.get("plane_radial_velocity_max_abs_error"),
        "plane_temperature_mean_rel_pct": _pct(plane.get("plane_temperature_mean_rel_error")),
        "plane_pressure_mean_rel_pct": _pct(plane.get("plane_pressure_mean_rel_error")),
        "exit_ppc_measured": _region(occupancy, "exit"),
        "plume_ppc_median": _region(occupancy, "plume"),
        "stationarity": ((stationarity or {}).get("assessment") or {}).get("verdict"),
        "diag_valid_volume_fraction": (diagnostics.get("validity") or {}).get(
            "valid_volume_fraction"),
        "diag_exit_KnD": exit_plane.get("measured_KnD_median"),
        "diag_exit_cell_over_mfp": exit_plane.get(
            "measured_cell_over_mfp_max_median"),
        "diag_cell_over_mfp_median": resolution.get("median"),
        "diag_cell_over_mfp_max": resolution.get("max"),
        "diag_volume_fraction_coarser_than_target": (
            resolution.get("exceeding_target") or {}).get("volume_fraction"),
        "diag_cells_to_resolve": budget.get("refined_cells_uncapped"),
        "reference": REFERENCE,
    }


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return "" if math.isnan(value) else f"{value:.6g}"
    return str(value)


def write_table(path: Path, rows: list) -> Path:
    """``study-table.csv``, one row per case, :data:`COLUMNS` in order."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(COLUMNS)
        for row in rows:
            writer.writerow([_cell(row.get(column)) for column in COLUMNS])
    return path


def write_markdown(path: Path, rows: list) -> Path:
    """A readable subset of the table, for the README of a run."""
    columns = ("case", "production_status", "state", "cell_over_mfp",
               "dt_over_tau0", "diag_exit_cell_over_mfp",
               "diag_volume_fraction_coarser_than_target",
               "centerline_density_mean_rel_pct",
               "angular_mass_flux_mean_rel_pct", "plane_temperature_mean_rel_pct",
               "stationarity")
    lines = ["| " + " | ".join(columns) + " |",
             "|" + "---|" * len(columns)]
    for row in rows:
        lines.append("| " + " | ".join(_cell(row.get(c)) for c in columns) + " |")
    lines += ["", f"Reference for every row: {REFERENCE}.", ""]
    Path(path).write_text("\n".join(lines), encoding="utf-8", newline="\n")
    return Path(path)


def _csv_columns(path: Path) -> dict | None:
    path = Path(path)
    if not path.is_file():
        return None
    with open(path, encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if not header:
            return None
        data = [[float(v) if v else np.nan for v in row] for row in reader]
    array = np.asarray(data, dtype=float).reshape(-1, len(header))
    return {name: array[:, i] for i, name in enumerate(header)}


def plot_knudsen(root: Path, entries: list, output_dir: Path) -> list[Path]:
    """The four Knudsen numbers on one figure, against the collisionless curve.

    Normalised, the collisionless solution is the same for every Kn -- it
    depends on ``S0`` alone -- so one reference curve serves all four, and the
    spread between the cases is the collisional departure from it. The last
    panel is the mesh resolution along the axis from ``Alldiagnose.py``
    (``Delta_max / lambda`` against the ``1/k`` target), so each departure can
    be read beside how well that case resolved its mean free path. A case that
    is unsuitable for production is labelled so in the legend.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []

    def results(entry, name):
        return _csv_columns(Path(root) / entry["path"] / "results" / name)

    curves = [(entry, results(entry, "centerline_extended.csv"),
               results(entry, "angular_r10D.csv"),
               results(entry, "centreline_diagnostics.csv"))
              for entry in sorted(entries, key=lambda e: -float(e["Kn"]))]
    if not any(any(c is not None for c in curve[1:]) for curve in curves):
        return []

    figure, axes = plt.subplots(2, 3, figsize=(15.0, 8.0))
    panels = (("n_over_n0", "$n/n_0$", True),
              ("U_sqrt_beta0", r"$U_x\sqrt{\beta_0}$", False),
              ("T_over_T0", "$T/T_0$", False))
    reference_drawn = {"centerline": False, "angular": False, "target": False}
    for entry, centerline, angular, diagnostics in curves:
        label = f"Kn = {float(entry['Kn']):g}" + (
            "" if entry.get("production_status") == "production-candidate"
            else " (NON-PRODUCTION)")
        if centerline is not None:
            for axis, (key, ylabel, log) in zip(axes[0], panels):
                if not reference_drawn["centerline"]:
                    axis.plot(centerline["x_over_D"],
                              centerline[f"{key}_analytical_axis"], "-",
                              color="black", linewidth=1.4,
                              label="Cai collisionless")
                axis.plot(centerline["x_over_D"], centerline[key], ".",
                          markersize=3, label=label)
                axis.set_xlabel("$X/D$")
                axis.set_ylabel(ylabel)
                if log:
                    axis.set_yscale("log")
            reference_drawn["centerline"] = True
        if angular is not None:
            if not reference_drawn["angular"]:
                axes[1, 0].plot(angular["theta_deg"], angular["n_over_n0_exact"],
                                "-", color="black", linewidth=1.4,
                                label="Cai collisionless")
                axes[1, 1].plot(angular["theta_deg"], angular["F_exact"], "-",
                                color="black", linewidth=1.4,
                                label="Cai collisionless")
                reference_drawn["angular"] = True
            axes[1, 0].plot(angular["theta_deg"], angular["n_over_n0"], ".",
                            markersize=4, label=label)
            axes[1, 1].plot(angular["theta_deg"], angular["F"], ".",
                            markersize=4, label=label)
        if diagnostics is not None and "cellOverMfpMax" in diagnostics:
            shown = diagnostics["cellOverMfpMax"] > 0.0
            axes[1, 2].plot(diagnostics["x_over_D"][shown],
                            diagnostics["cellOverMfpMax"][shown], "-",
                            linewidth=1.0, label=label)
            if not reference_drawn["target"]:
                axes[1, 2].axhline(1.0 / 3.0, color="0.4", linestyle="--",
                                   linewidth=0.9, label="target 1/3")
                reference_drawn["target"] = True
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_xlabel(r"$\theta$ at r/D = 10 [deg]")
    axes[1, 0].set_ylabel("$n/n_0$")
    axes[1, 1].set_xlabel(r"$\theta$ at r/D = 10 [deg]")
    axes[1, 1].set_ylabel(r"$F$ [1/sr]")
    axes[1, 2].set_yscale("log")
    axes[1, 2].set_xlabel("$X/D$")
    axes[1, 2].set_ylabel(r"centreline $\Delta_{max}/\lambda$ (Alldiagnose.py)")
    for axis in axes.ravel():
        axis.grid(True, which="both", alpha=0.3)
        if axis.get_legend_handles_labels()[0]:
            axis.legend(fontsize=7)
    figure.suptitle("The Knudsen extension at the baseline weight -- reference: "
                    "Cai's collisionless ANALYTICAL solution; no published DSMC "
                    "at these Kn", fontsize=11)
    figure.tight_layout()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "knudsen-comparison.png"
    figure.savefig(path, dpi=130)
    plt.close(figure)
    return [path]
