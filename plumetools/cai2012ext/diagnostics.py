r"""The nozzle-plume family's rarefaction and mesh-resolution diagnostics, here.

R-4D's ``./Allpost`` runs ``python -m plumetools.nozzleplume diagnostics``:
:func:`plumetools.rarefaction.run.diagnose` fed from a nozzle-plume
``case.yaml``. A Cai case has another schema, so this module feeds the SAME
function from a cai2012 ``case.yaml`` and writes the same outputs to the same
places. Nothing in :mod:`plumetools.rarefaction` or
:mod:`plumetools.nozzleplume` is changed or copied; the centreline figure is
:func:`plumetools.nozzleplume.diagnostics.plot_centreline` itself.

=============================  ==============================================
from the case                  used for
=============================  ==============================================
``case.yaml`` ``gas``          the species and ``mean_free_path_convention``
``case.yaml`` ``nozzle``       ``D``, the length ``KnD = lambda / D`` is on --
                               the length the case's own Kn is defined on
``case.yaml`` ``post``         the tube and stations of ``centerline.csv``
``diagnostics.yaml``           the nozzle-plume ``diagnostics:`` section --
                               thresholds, practical minimum, budget. A file of
                               its own because the cai2012 schema rejects
                               sections it does not know
=============================  ==============================================

Outputs, as in the nozzle-plume family
--------------------------------------
=====================================  ========================================
``<time>/mfp``, ``KnD``, ``KnCellMax``,  volume fields beside ``rhoNMean``: open
``mfpResolved``, ``diagValid``, ...    ``open.foam`` in ParaView. -1 where
                                       there is no value
``results/diagnostics_summary.yaml``   validity, Delta/lambda statistics, the
                                       worst cells, the lambda/k mesh estimate
``results/centreline_diagnostics.csv`` every diagnostic along the axis
``results/centreline-diagnostics.png`` the same, drawn
=====================================  ========================================

The summary gains one block of this family's, ``exit_plane``: the measured
``KnD`` and ``Delta_max/lambda`` in the exit cells beside the target Kn and the
pre-run ``cell/lambda0`` -- the check that what generation planned is what the
run produced.

There is no background cell: the mesh is graded Cartesian (``blockMesh``), not
an octree, so no refinement levels are listed. The planning block still gives
every cell's ``lambda/k`` target and lower bounds on a mesh that met it, with
and without ``practical_min_cell_size_m``. See ``docs/rarefaction-diagnostics.md``
for what every field and statistic means, and what none of them proves.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml

from plumetools.cai2012 import inflow, mesh
from plumetools.cai2012.config import load_case_config
from plumetools.cai2012.geometry import from_config
from plumetools.cai2012ext import FAMILY_DIRECTORY, policy
from plumetools.nozzleplume.diagnostics import CENTRELINE_FIELDS, plot_centreline
from plumetools.rarefaction import summary
from plumetools.rarefaction.run import diagnose, file_sha256
from plumetools.rarefaction.settings import DiagnosticSettings, SettingsError

#: The per-case settings file, copied from ``baseCase/`` by generation.
SETTINGS_FILE = "diagnostics.yaml"

#: Added to every summary's notes.
GRADED_MESH_NOTE = (
    "The mesh is graded Cartesian (blockMesh), not an octree: no refinement "
    "levels are listed, and the mesh estimates assume each cell could take its "
    "own lambda/k. A graded mesh cannot -- refining the core refines a whole "
    "lattice of planes -- so they are lower bounds on a lower bound.")


def load_settings(case_dir: Path) -> DiagnosticSettings:
    """The ``diagnostics:`` section of the case's ``diagnostics.yaml``.

    Raises:
        SettingsError: for a missing file, an unknown top-level key, or a bad
            setting -- naming it.
    """
    path = Path(case_dir) / SETTINGS_FILE
    if not path.is_file():
        raise SettingsError(
            f"no {SETTINGS_FILE} in {case_dir}. It is copied from baseCase/ by "
            f"generate_cases.py; ./generate_cases.py --refresh adds it to a case "
            f"generated before it existed.")
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    unknown = sorted(set(raw) - {"diagnostics"})
    if unknown:
        raise SettingsError(f"{path}: unknown key(s) {unknown}; the gas and D "
                            f"come from case.yaml, so only 'diagnostics' belongs "
                            f"here")
    return DiagnosticSettings.from_mapping(raw.get("diagnostics"),
                                           f"{path}:diagnostics")


def centreline_profile(diagnosis, cfg, geom) -> dict:
    """The diagnostic fields along the axis, on ``centerline.csv``'s stations.

    The tube and bins of :func:`plumetools.cai2012.post.centerline` -- cells
    within ``post.centerline_radius_over_D`` of the axis, ``post.centerline_
    points`` stations over ``[0, post.centerline_x_over_D_max]`` -- so a
    station here is the same cells as a station there. A value is the
    volume-weighted mean over the station's cells that have one;
    ``valid_cells`` says how many of them are trusted (``diagValid``).
    """
    data = diagnosis.data
    centres = data.inputs.get("cell_centres")
    volumes = data.inputs.get("cell_volume")
    if centres is None or volumes is None:
        raise ValueError("the centreline needs cell centres and volumes "
                         "(./Allpost writes them)")
    d = float(geom.diameter_m)
    x = centres[:, 0]
    radial = np.hypot(centres[:, 1], centres[:, 2])
    x_max = float(cfg.post.centerline_x_over_D_max) * d
    inside = ((radial <= float(cfg.post.centerline_radius_over_D) * d)
              & (x >= 0.0) & (x <= x_max))
    n_bins = int(cfg.post.centerline_points)
    edges = np.linspace(0.0, x_max, n_bins + 1)
    index = np.clip(np.digitize(x[inside], edges) - 1, 0, n_bins - 1)
    counts = np.bincount(index, minlength=n_bins)
    held = counts > 0

    valid = (diagnosis.context["valid"][inside]
             if not diagnosis.context.missing("valid")
             else np.zeros(int(inside.sum()), dtype=bool))
    station = 0.5 * (edges[:-1] + edges[1:])
    profile = {
        "x_m": station[held],
        "x_over_D": station[held] / d,
        "n_cells": counts[held].astype(float),
        "valid_cells": np.bincount(index, weights=valid.astype(float),
                                   minlength=n_bins)[held],
    }
    volume = volumes[inside]
    for name in CENTRELINE_FIELDS:
        field = diagnosis.outputs.get(name)
        if field is None:
            continue
        values = field.values[inside]
        good = np.isfinite(values) & (values >= 0.0)    # -1 is "no value"
        weight = np.where(good, volume, 0.0)
        total = np.bincount(index, weights=weight, minlength=n_bins)
        summed = np.bincount(index, weights=np.where(good, values, 0.0) * weight,
                             minlength=n_bins)
        with np.errstate(divide="ignore", invalid="ignore"):
            profile[name] = np.where(total > 0.0, summed / total, np.nan)[held]
    profile["time"] = data.time
    return profile


def exit_plane(diagnosis, cfg, geom, plan) -> dict:
    """Measured resolution in the exit cells, beside what generation planned.

    The exit cells are the first layer of the uniform core inside the nozzle
    disk -- the cells whose volume set the particle weight. Their ``KnD`` should
    sit near the case's target Kn (the gas there is the exit Maxwellian, with
    ``T`` a little below ``T0`` because only the outgoing half is present), and
    their ``cellOverMfpMax`` near the pre-run ``cell/lambda0``.
    """
    centres = diagnosis.data.inputs.get("cell_centres")
    if centres is None:
        return {"evaluated": False, "reason": "no cell centres"}
    radial = np.hypot(centres[:, 1], centres[:, 2])
    cells = (centres[:, 0] <= float(plan.x_segments[0].first_cell_m)) & (
        radial <= float(geom.radius_m))
    valid = (diagnosis.context["valid"] if not diagnosis.context.missing("valid")
             else np.zeros(centres.shape[0], dtype=bool))
    out = {"cells": int(cells.sum()), "valid_cells": int((cells & valid).sum()),
           "target_knudsen": float(cfg.exit.knudsen),
           "planned_cell_over_mfp": float(plan.cell_over_mfp)}
    for name, key in (("KnD", "measured_KnD_median"),
                      ("cellOverMfpMax", "measured_cell_over_mfp_max_median"),
                      ("parcelsPerStep", "measured_parcels_per_step_median")):
        field = diagnosis.outputs.get(name)
        if field is None:
            continue
        pick = cells & valid & np.isfinite(field.values)
        out[key] = float(np.median(field.values[pick])) if pick.any() else None
    return out


def run(case_dir: Path, *, time: str | None = None,
        results: Path | None = None, write_fields: bool = True) -> tuple:
    """Compute, write and report the diagnostics of one case.

    Args:
        case_dir: a reconstructed, post-processed case (``./Allpost`` has
            written ``C`` and ``V``).
        time: the time directory; the latest by default.
        results: where the summary and centreline go; ``<case>/results``.
        write_fields: write the volume fields into the time directory.

    Returns:
        ``(fields, written, document)``: the field files written into the time
        directory, the summary/centreline/figure paths, and the summary.

    Raises:
        policy.ExtensionError: for a case at an excluded Kn.
        SettingsError, plumetools.rarefaction.foam.PostError, ValueError: as
            :func:`plumetools.rarefaction.run.diagnose`.
    """
    from plumetools.cai2012ext import interface, interface_case

    case_dir = Path(case_dir)
    cfg = load_case_config(case_dir)
    sphere = interface.load_interface(case_dir)
    if sphere is not None:
        problem = interface_case.knudsen_problem(cfg, sphere, str(case_dir))
        if problem:
            raise policy.ExtensionError(problem)
    else:
        problem = policy.direct_exit_problem(cfg, str(case_dir))
        if problem:
            raise policy.ExtensionError(problem)
    geom = from_config(cfg)
    exit_state = inflow.from_config(cfg)
    # In interface mode the exit cells are not in the mesh: there is no exit
    # plane to check and no direct-exit plan to check it against.
    plan = mesh.plan(cfg, geom, exit_state) if sphere is None else None
    settings = load_settings(case_dir)
    label = cfg.meta.get("case_name") or case_dir.name

    diagnosis = diagnose(
        case_dir, cfg.species(), settings,
        convention=cfg.gas.mean_free_path_convention, time=time,
        reference_lengths={"D": float(geom.diameter_m)},
        background_cell_size_m=None, write_fields=write_fields,
        header={"case": label, "family": FAMILY_DIRECTORY,
                "case_yaml_sha256": file_sha256(case_dir / "case.yaml"),
                "diagnostics_yaml_sha256": file_sha256(case_dir / SETTINGS_FILE),
                "knudsen": float(cfg.exit.knudsen),
                "collision_model": cfg.dsmc.binary_collision_model,
                "mesh_type": cfg.mesh.type})
    document = diagnosis.document
    if plan is not None:
        document["exit_plane"] = exit_plane(diagnosis, cfg, geom, plan)
    else:
        document["exit_plane"] = {
            "evaluated": False,
            "reason": "interface mode: the exit cells were removed with the "
                      "dense core; the inlet is the downstream sphere"}
        document["interface"] = interface_case.post_context(case_dir).interface
    document.setdefault("notes", []).append(GRADED_MESH_NOTE)

    print("\n".join(summary.report(document)))
    measured = document["exit_plane"]
    if measured.get("measured_KnD_median") is not None:
        print(f"  exit cells: KnD {measured['measured_KnD_median']:.4g} against "
              f"Kn = {measured['target_knudsen']:g}; Delta_max/lambda "
              f"{measured.get('measured_cell_over_mfp_max_median', float('nan')):.4g}"
              f" against the planned {measured['planned_cell_over_mfp']:.4g}")
    print()

    out = Path(results) if results else case_dir / "results"
    written = [summary.write_yaml(out / "diagnostics_summary.yaml", document)]
    try:
        profile = centreline_profile(diagnosis, cfg, geom)
    except ValueError as exc:
        print(f"    no diagnostic centreline: {exc}")
    else:
        from plumetools.cai2012.post import write_csv

        columns = ["x_m", "x_over_D", "n_cells", "valid_cells",
                   *[name for name in CENTRELINE_FIELDS if name in profile]]
        written.append(write_csv(out / "centreline_diagnostics.csv",
                                 {key: profile[key] for key in columns}))
        figure = plot_centreline(out / "centreline-diagnostics.png", profile,
                                 settings, label)
        if figure is not None:
            written.append(figure)
    return list(diagnosis.written), written, document
