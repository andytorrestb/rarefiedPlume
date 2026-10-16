r"""Rarefaction and mesh-resolution diagnostics for a nozzle-plume case.

:mod:`plumetools.rarefaction` does the work; this supplies what only a
``case.yaml`` of this family knows, and puts the results where the rest of the
family's post-processing puts its own:

=========================  ====================================================
from ``case.yaml``         used for
=========================  ====================================================
``gas``                    the species and ``mean_free_path_convention``
``nozzle.exit_diameter_m`` ``D``, the length ``KnD = lambda / D`` is built on
``mesh.background_cell_  the octree levels a practical minimum allows
size_m``
``post.centreline_*``      the same tube and stations as ``centreline.csv``
``diagnostics``            thresholds and the rest; optional
=========================  ====================================================

Outputs
-------
=====================================  ========================================
``<time>/mfp``, ``<time>/KnD``, ...    volume fields beside ``rhoNMean``; open
                                       ``open.foam`` in ParaView
``results/diagnostics_summary.yaml``   the summary
``results/centreline_diagnostics.csv`` the fields along the plume axis
``results/centreline-diagnostics.png`` the same, drawn (needs matplotlib)
=====================================  ========================================

A centreline value is the volume-weighted mean of the written field over the
cells of that station -- the cells ``centreline.csv`` averages -- that have a
value, valid or low-confidence; ``valid_cells`` says how many are trusted, and
the figure draws a station with none of them hollow. A station with no value
at all is ``nan``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml

from plumetools.nozzleplume import post
from plumetools.rarefaction import summary
from plumetools.rarefaction.run import diagnose, file_sha256

#: Diagnostic fields sampled along the centreline, in column order.
CENTRELINE_FIELDS = ("mfp", "KnD", "cellSizeV", "cellSizeMax", "cellSizeTarget",
                     "cellOverMfpV", "cellOverMfpMax", "KnGradN", "KnGradT",
                     "parcelsPerStep", "sampledParcels")


def reference_lengths(geom) -> dict:
    """``{"D": exit diameter}`` -- the case's own Knudsen length."""
    return {"D": float(geom.diameter_m)}


def _provenance(case_dir: Path) -> dict:
    """What identifies this case's inputs: the case.yaml hash, and the
    iteration it is, if :mod:`plumetools.iteration` prepared it."""
    document = {"case_yaml_sha256": file_sha256(Path(case_dir) / "case.yaml")}
    iteration = Path(case_dir) / "iteration.yaml"
    if iteration.is_file():
        record = yaml.safe_load(iteration.read_text(encoding="utf-8")) or {}
        document["iteration"] = {key: record.get(key)
                                 for key in ("name", "parent", "created_utc")}
    return document


def centreline_profile(case, diagnosis) -> dict:
    """The diagnostic fields along the plume axis. See the module docstring."""
    data = diagnosis.data
    centres = data.inputs.get("cell_centres")
    volumes = data.inputs.get("cell_volume")
    if centres is None or volumes is None:
        raise post.PostError("the centreline needs cell centres and volumes")
    inside, index, edges = post.centreline_tube(centres, volumes, case.cfg,
                                                case.geom)
    n_bins = edges.size - 1
    volume = volumes[inside]
    counts = np.bincount(index, minlength=n_bins)
    held = counts > 0

    valid = diagnosis.context["valid"][inside] if not diagnosis.context.missing(
        "valid") else np.zeros(inside.sum(), dtype=bool)
    valid_counts = np.bincount(index, weights=valid.astype(float), minlength=n_bins)
    station = 0.5 * (edges[:-1] + edges[1:])
    profile = {
        "x_m": station[held],
        "x_over_D": station[held] / case.geom.diameter_m,
        "n_cells": counts[held].astype(float),
        "valid_cells": valid_counts[held],
    }
    for name in CENTRELINE_FIELDS:
        field = diagnosis.outputs.get(name)
        if field is None:
            continue
        good = np.isfinite(field.values[inside])
        weight = np.where(good, volume, 0.0)
        total = np.bincount(index, weights=weight, minlength=n_bins)
        summed = np.bincount(index, weights=np.where(good, field.values[inside], 0.0)
                             * weight, minlength=n_bins)
        with np.errstate(divide="ignore", invalid="ignore"):
            profile[name] = np.where(total > 0.0, summed / total, np.nan)[held]
    profile["time"] = data.time
    return profile


def plot_centreline(path: Path, profile: dict, settings, label: str) -> Path | None:
    """Three panels: lengths, resolution, Knudsen numbers. ``None`` without
    matplotlib."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    x = profile["x_m"]
    trusted = profile["valid_cells"] > 0
    figure, panels = plt.subplots(3, 1, figsize=(6.5, 8.5), sharex=True)
    for axes, entries, ylabel in (
            (panels[0], (("mfp", "mean free path $\\lambda$"),
                         ("cellSizeMax", "cell $\\Delta_{max}$"),
                         ("cellSizeTarget", "target $\\lambda/k$")), "length [m]"),
            (panels[1], (("cellOverMfpMax", "$\\Delta_{max}/\\lambda$"),
                         ("cellOverMfpV", "$V^{1/3}/\\lambda$")),
             "cell / mean free path"),
            (panels[2], (("KnD", "$Kn_D = \\lambda/D$"),
                         ("KnGradN", "$\\lambda|\\nabla n|/n$"),
                         ("KnGradT", "$\\lambda|\\nabla T|/T$")), "Knudsen number")):
        for key, legend in entries:
            values = profile.get(key)
            if values is None:
                continue
            shown = np.isfinite(values) & (values > 0.0)
            if shown.any():
                # Hollow at every station; filled where it holds valid cells.
                line, = axes.plot(x[shown], values[shown], "o-", markersize=2.5,
                                  linewidth=0.8, mfc="none", label=legend)
                axes.plot(x[shown & trusted], values[shown & trusted], "o",
                          markersize=2.5, color=line.get_color())
        axes.set_yscale("log")
        axes.set_ylabel(ylabel)
        axes.grid(True, which="both", alpha=0.3)
    panels[0].axhline(settings.practical_min_cell_size_m, color="0.4",
                      linestyle="--", linewidth=0.9,
                      label=f"practical minimum {settings.practical_min_cell_size_m:g} m")
    panels[1].axhline(settings.target_ratio, color="0.4", linestyle="--",
                      linewidth=0.9, label=f"target $1/k$ = {settings.target_ratio:.3g}")
    for axes in panels:
        handles, _ = axes.get_legend_handles_labels()
        if handles:
            axes.legend(fontsize=7, loc="best")
    panels[2].set_xlabel("x [m] from the exit plane   (hollow: no valid cell at "
                         "that station, diagValid = 0)")
    panels[0].set_title(f"Centreline diagnostics -- {label}, t = {profile['time']} s",
                        fontsize=10)
    figure.tight_layout()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return path


def run(case, *, time: str | None = None,
        results: Path | None = None) -> tuple[list[Path], list[Path]]:
    """Compute, write and report the diagnostics of a loaded case.

    Args:
        case: a :class:`plumetools.nozzleplume.workflow.Case`.
        time: the time directory; the latest by default.
        results: where the summary and centreline go; ``<case>/results``.

    Returns:
        ``(fields, results)``: the field files written into the time
        directory, and the summary, centreline and figure.

    Raises:
        plumetools.nozzleplume.post.PostError: see
            :func:`plumetools.rarefaction.run.diagnose`.
    """
    cfg = case.cfg
    settings = cfg.diagnostics
    diagnosis = diagnose(
        case.case_dir, cfg.species(), settings,
        convention=cfg.gas.mean_free_path_convention, time=time,
        reference_lengths=reference_lengths(case.geom),
        background_cell_size_m=float(cfg.mesh.background_cell_size_m),
        header={"case": case.label, **_provenance(case.case_dir),
                "collision_model": cfg.dsmc.collision_model})

    document = diagnosis.document
    if cfg.dsmc.collision_model == "NoBinaryCollision":
        document.setdefault("notes", []).append(
            "This run was COLLISIONLESS (dsmc.collision_model: NoBinaryCollision). "
            "lambda is the gas's mean free path; the simulated molecules never "
            "collided, so no mesh resolves anything about collisions here.")
    print("\n".join(summary.report(document)))
    print()

    out = Path(results) if results else case.case_dir / "results"
    written = [summary.write_yaml(out / "diagnostics_summary.yaml", document)]
    try:
        profile = centreline_profile(case, diagnosis)
    except post.PostError as exc:
        print(f"    no diagnostic centreline: {exc}")
    else:
        columns = ["x_m", "x_over_D", "n_cells", "valid_cells",
                   *[name for name in CENTRELINE_FIELDS if name in profile]]
        written.append(post.write_csv(out / "centreline_diagnostics.csv",
                                      {key: profile[key] for key in columns}))
        figure = plot_centreline(out / "centreline-diagnostics.png", profile,
                                 settings, case.label)
        if figure is not None:
            written.append(figure)
    return list(diagnosis.written), written
