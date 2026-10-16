r"""Condense the diagnostic fields into what a remeshing decision needs.

Every statistic here is over **valid** cells only, weighted by volume where it
says so; an invalid cell is counted in the validity block and nowhere else.

Refinement levels
-----------------
An octree mesh offers only the sizes ``B / 2**L`` of its background cell ``B``.
A practical minimum ``h_min`` therefore allows the levels with
``B / 2**L >= h_min``; for ``B = 0.5 m`` and ``h_min = 10 mm`` that is up to
level 5, 15.625 mm, and level 6, 7.8125 mm, is excluded. The minimum is a
planning constraint. It decides which levels are listed as available; it never
clips the target ``lambda / k`` that is reported, and a cell whose target is
below it is reported as such rather than hidden.

The mesh estimate
-----------------
``refined_cells_*`` is what the mesh would hold if every valid cell were
refined -- never coarsened -- to its target, cell by cell:
``sum V / h**3``. ``within_practical_minimum`` stops each cell at the finest
permitted level; ``uncapped`` does not. Both leave invalid cells as they are and
both omit the buffer layers ``snappyHexMesh`` adds between levels, so they are
lower bounds, good for comparing a plan with the budget and not for more.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from plumetools.rarefaction.core import STATUS_FLAGS

#: Read into every summary: what the numbers in it are, and are not.
NOTES = [
    "cellOverMfp* (Delta/lambda) and KnCell* (lambda/Delta) measure MESH "
    "RESOLUTION: how many mean free paths one cell spans. KnD and the other "
    "reference-length numbers, and KnGradN / KnGradT, characterise RAREFACTION, "
    "on the nozzle scale and on the local density and temperature scales. They "
    "answer different questions.",
    "No one of these numbers, and no threshold on any of them, establishes that "
    "the DSMC solution is valid. The time step against the collision time, the "
    "parcels per cell the collisions are selected from, the sampling and the "
    "length of the transient all matter as well.",
    "lambda is the equilibrium mean free path (equation G1) at the local "
    "translational temperature, from time-averaged n and T. In a strongly "
    "non-equilibrium region -- the expanding plume, the shock layer, where "
    "re-emitted wall gas meets the jet -- it is a nominal length.",
    "Statistics exclude invalid cells. In the written fields those hold "
    "invalid_value; diagStatus says why.",
    "Gradients are least-squares fits of the averaged fields over valid face "
    "neighbours. They amplify sampling noise, most where cells are smaller than "
    "lambda, and are first-order at mesh-level transitions.",
]


def _plain(value):
    """numpy scalars and arrays to what yaml.safe_dump writes."""
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.ndarray):
        return [_plain(v) for v in value.tolist()]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        value = float(value)
        return value if math.isfinite(value) else None
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def statistics(values, valid, volumes) -> dict:
    """Min, median, 95th percentile, max and volume-weighted mean over valid cells."""
    valid = np.asarray(valid, dtype=bool)
    if not valid.any():
        return {"cells": 0}
    v = np.asarray(values)[valid]
    w = np.asarray(volumes)[valid]
    return {
        "cells": int(valid.sum()),
        "min": float(v.min()),
        "median": float(np.median(v)),
        "p95": float(np.percentile(v, 95.0)),
        "max": float(v.max()),
        "volume_weighted_mean": float((v * w).sum() / w.sum()),
    }


def fraction_above(values, valid, volumes, threshold: float) -> dict:
    """Cells and volume above ``threshold``, as fractions of the VALID ones."""
    valid = np.asarray(valid, dtype=bool)
    above = valid & (np.nan_to_num(values, nan=-np.inf) > threshold)
    total_volume = float(np.asarray(volumes)[valid].sum())
    return {
        "threshold": float(threshold),
        "cells": int(above.sum()),
        "cell_fraction": float(above.sum() / max(valid.sum(), 1)),
        "volume_fraction": (float(np.asarray(volumes)[above].sum() / total_volume)
                            if total_volume > 0 else 0.0),
    }


def bounds(centres, mask) -> list | None:
    """``[x0, x1, y0, y1, z0, z1]`` of the cell centres in ``mask``."""
    if not np.any(mask):
        return None
    points = np.asarray(centres)[mask]
    lo, hi = points.min(axis=0), points.max(axis=0)
    return [float(lo[0]), float(hi[0]), float(lo[1]), float(hi[1]),
            float(lo[2]), float(hi[2])]


def worst_cells(ratio, valid, centres, sizes, n: int,
                separation: float) -> list[int]:
    """Indices of the ``n`` largest valid ratios, spread out.

    A cell is skipped if it is within ``separation`` of its own size of one
    already chosen, so the list names ``n`` places rather than ``n`` neighbours.
    """
    valid = np.asarray(valid, dtype=bool) & np.isfinite(ratio)
    candidates = np.flatnonzero(valid)
    if candidates.size == 0 or n <= 0:
        return []
    # Enough candidates that the suppression cannot run out on a large mesh.
    take = min(candidates.size, max(1000, 200 * n))
    top = candidates[np.argpartition(-ratio[candidates], take - 1)[:take]]
    top = top[np.argsort(-ratio[top], kind="stable")]
    chosen: list[int] = []
    for index in top:
        reach = separation * sizes[index]
        if all(np.linalg.norm(centres[index] - centres[c]) > reach for c in chosen):
            chosen.append(int(index))
            if len(chosen) == n:
                break
    return chosen


def refinement_levels(background_cell_size_m: float,
                      practical_min_cell_size_m: float) -> dict:
    """The octree levels a practical minimum permits, and the first it does not."""
    b = float(background_cell_size_m)
    h_min = float(practical_min_cell_size_m)
    finest = int(math.floor(math.log2(b / h_min) + 1e-9)) if b >= h_min else -1
    rows = [{"level": level, "cell_size_m": b / 2 ** level,
             "permitted": level <= finest} for level in range(max(finest, 0) + 2)]
    return {
        "background_cell_size_m": b,
        "practical_min_cell_size_m": h_min,
        "finest_permitted_level": finest,
        "finest_permitted_cell_size_m": b / 2 ** finest if finest >= 0 else None,
        "first_excluded_level": finest + 1,
        "first_excluded_cell_size_m": b / 2 ** (finest + 1),
        "levels": rows,
    }


def required_level(target, background_cell_size_m: float) -> np.ndarray:
    """The coarsest level whose cell is no larger than ``target``; >= 0."""
    with np.errstate(divide="ignore", invalid="ignore"):
        level = np.ceil(np.log2(float(background_cell_size_m) / target) - 1e-9)
    return np.maximum(level, 0.0)


def mesh_estimate(volumes, current_size, new_size, valid) -> float:
    """``sum V / h**3``, each valid cell refined -- never coarsened -- to
    ``new_size``. See the module docstring."""
    h = np.where(valid, np.minimum(current_size, np.nan_to_num(new_size, nan=np.inf)),
                 current_size)
    return float((np.asarray(volumes) / h ** 3).sum())


def build(outputs: dict, ctx, data, settings, *, skipped: dict,
          background_cell_size_m: float | None = None,
          header: dict | None = None) -> dict:
    """The summary document. See the module docstring for its parts.

    Args:
        outputs: ``{name: Field}`` from :func:`plumetools.rarefaction.core.evaluate`.
        ctx: the :class:`~plumetools.rarefaction.core.Context` they came from.
        data: the :class:`~plumetools.rarefaction.foam.CaseData` read.
        settings: the :class:`~plumetools.rarefaction.settings.DiagnosticSettings`.
        skipped: ``{name: reason}`` for the diagnostics that could not run.
        background_cell_size_m: the level-0 cell, if the mesh is an octree
            refinement of a uniform one; enables the level planning.
        header: case identification, placed first.
    """
    document: dict = dict(header or {})
    document["time"] = data.time
    document["sampling"] = {
        "window": data.window.as_dict() if data.window else None,
        "fields": {role: settings.input_field(role) for role in (
            "number_density", "mass_density", "momentum_density",
            "kinetic_energy_density", "parcels_per_step")},
        "temperature": "translational, from the averaged moments: 3 R T = "
                       "2 linearKE/rhoM - |momentum/rhoM|^2",
    }
    species = ctx.species
    document["settings"] = {
        "gas": {"molecular_diameter_m": float(species.diameter_m),
                "omega": float(species.omega), "t_ref_K": float(species.t_ref_K),
                "mean_free_path_convention": ctx.convention},
        **settings.as_dict(),
    }
    document["checks"] = data.checks
    document["skipped"] = dict(skipped)

    volumes = data.inputs.get("cell_volume")
    n_cells = data.n_cells
    if volumes is None or ctx.missing("valid"):
        document["validity"] = {"cells": n_cells, "evaluated": False,
                                "reason": {**ctx.missing("valid"),
                                           **ctx.missing("cell_volume")}}
        document["notes"] = NOTES
        return _plain(document)

    valid = ctx["valid"]
    status = ctx["status"]
    centres = data.inputs.get("cell_centres")
    total_volume = float(volumes.sum())
    validity = {
        "cells": n_cells,
        "valid_cells": int(valid.sum()),
        "valid_cell_fraction": float(valid.mean()),
        "valid_volume_fraction": float(volumes[valid].sum() / total_volume),
        "invalid_cells": int((~valid).sum()),
        "invalid_cell_fraction": float((~valid).mean()),
        "invalid_volume_fraction": float(volumes[~valid].sum() / total_volume),
        # Where each kind of invalid cell is, so an excluded region -- the
        # exit, when it holds too few parcels -- is named rather than lost.
        "invalid_by_reason": {
            reason: {"cells": int(((status & flag) > 0).sum()),
                     "volume_fraction": float(
                         volumes[(status & flag) > 0].sum() / total_volume),
                     "bounds_m": (bounds(centres, (status & flag) > 0)
                                  if centres is not None else None)}
            for flag, reason in STATUS_FLAGS.items()},
    }
    for name, label in (("grad_number_density", "number_density"),
                        ("grad_temperature", "temperature")):
        if not ctx.missing(name):
            ok = ctx[name].valid
            validity.setdefault("gradient_fitted", {})[label] = {
                "cells": int(ok.sum()),
                "fraction_of_valid_cells": float(ok.sum() / max(valid.sum(), 1))}
    document["validity"] = validity

    per_step = np.asarray(data.inputs.get("parcels_per_step", np.zeros(n_cells)))
    particles = {"threshold": float(settings.min_sampled_particles),
                 "basis": settings.particle_count_basis,
                 "rule": "valid only with MORE sampled parcels than the threshold"}
    occupied = per_step > 0.0
    for basis, count in (("per_step", per_step),
                         ("accumulated", per_step * float(data.window.n_steps)
                          if data.window else None)):
        if count is None:
            continue
        particles[basis] = {
            "median_over_occupied_cells": float(np.median(count[occupied]))
            if occupied.any() else 0.0,
            "cells_at_or_below_threshold": int(
                (count <= settings.min_sampled_particles).sum()),
            "volume_fraction_at_or_below_threshold": float(
                volumes[count <= settings.min_sampled_particles].sum() / total_volume),
        }
    document["particles"] = particles

    # --- resolution ----------------------------------------------------------
    ratio_k = settings.target_ratio
    resolution: dict = {"target_cell_over_mean_free_path": ratio_k,
                        "cells_per_mean_free_path": float(
                            settings.cells_per_mean_free_path)}
    for measure, name in (("volume_equivalent", "cellOverMfpV"),
                          ("bounding_box_max", "cellOverMfpMax")):
        if name in outputs:
            field = outputs[name]
            resolution[measure] = {
                "field": name,
                **statistics(field.values, field.valid, volumes),
                "exceeding_target": fraction_above(field.values, field.valid,
                                                   volumes, ratio_k),
            }
    primary = outputs.get("cellOverMfpMax") or outputs.get("cellOverMfpV")
    size_field = outputs.get("cellSizeMax") or outputs.get("cellSizeV")
    if primary is not None and centres is not None and size_field is not None:
        sizes = size_field.values
        lam = outputs["mfp"].values if "mfp" in outputs else None
        target = outputs["cellSizeTarget"].values if "cellSizeTarget" in outputs else None
        chosen = worst_cells(primary.values, primary.valid, centres, sizes,
                             int(settings.n_worst_cells),
                             float(settings.worst_cells_separation))
        rows = []
        for index in chosen:
            row = {"at_m": [float(v) for v in centres[index]],
                   "cell_over_mean_free_path": float(primary.values[index]),
                   "cell_size_m": float(sizes[index])}
            if lam is not None:
                row["mean_free_path_m"] = float(lam[index])
            if target is not None:
                row["target_cell_size_m"] = float(target[index])
                if background_cell_size_m:
                    row["required_level"] = int(required_level(
                        target[index], background_cell_size_m))
            rows.append(row)
        resolution["worst_cells"] = rows
        resolution["by_cell_size"] = _by_cell_size(
            primary, sizes, centres, volumes, ratio_k)
    document["resolution"] = resolution

    # --- rarefaction ----------------------------------------------------------
    document["rarefaction"] = {
        name: {"description": ctx.registry[name].output.description,
               **statistics(field.values, field.valid, volumes)}
        for name, field in outputs.items()
        if name.startswith("Kn") and not name.startswith("KnCell")}

    # --- planning -------------------------------------------------------------
    if "cellSizeTarget" in outputs:
        document["planning"] = _planning(
            outputs, centres, volumes, n_cells, settings, background_cell_size_m)
    document["notes"] = NOTES
    return _plain(document)


def _by_cell_size(primary, sizes, centres, volumes, ratio_k) -> list[dict]:
    """One row per cell size -- per octree level on an octree mesh."""
    valid = primary.valid
    with np.errstate(divide="ignore", invalid="ignore"):
        key = np.round(np.log2(sizes))
    rows = []
    for value in np.unique(key[np.isfinite(key)]):
        group = key == value
        good = group & valid
        over = good & (np.nan_to_num(primary.values, nan=-np.inf) > ratio_k)
        row = {"cell_size_m": float(np.median(sizes[group])),
               "cells": int(group.sum()), "valid_cells": int(good.sum()),
               "exceeding_cells": int(over.sum()),
               "exceeding_volume_fraction": float(
                   volumes[over].sum() / volumes[good].sum()) if good.any() else 0.0}
        if good.any():
            worst = np.flatnonzero(good)[np.argmax(primary.values[good])]
            row["max_cell_over_mean_free_path"] = float(primary.values[worst])
            row["worst_at_m"] = [float(v) for v in centres[worst]]
            row["exceeding_bounds_m"] = bounds(centres, over)
        rows.append(row)
    return rows


def _planning(outputs, centres, volumes, n_cells, settings,
              background_cell_size_m) -> dict:
    target = outputs["cellSizeTarget"]
    valid = target.valid
    h_min = float(settings.practical_min_cell_size_m)
    current = (outputs["cellSizeV"].values if "cellSizeV" in outputs
               else np.cbrt(volumes))
    planning: dict = {}

    # What the next mesh could actually give each cell: on an octree, the
    # level sizes B / 2**L -- uncapped, or stopped at the finest permitted.
    if background_cell_size_m:
        b = float(background_cell_size_m)
        levels = refinement_levels(b, h_min)
        planning["refinement_levels"] = levels
        need = required_level(np.where(valid, target.values, np.inf), b)
        uncapped_size = b / 2.0 ** need
        capped_size = b / 2.0 ** np.minimum(
            need, max(levels["finest_permitted_level"], 0))
        histogram = {}
        for level in np.unique(need[valid]).astype(int):
            pick = valid & (need == level)
            histogram[int(level)] = {
                "cells": int(pick.sum()),
                "volume_fraction_of_valid": float(
                    volumes[pick].sum() / volumes[valid].sum()),
                "permitted": bool(level <= levels["finest_permitted_level"])}
        planning["required_level_histogram"] = histogram
    else:
        uncapped_size = target.values
        capped_size = np.maximum(target.values, h_min)

    below = valid & (np.nan_to_num(target.values, nan=np.inf) < h_min)
    smallest = (int(np.flatnonzero(valid)[np.argmin(target.values[valid])])
                if valid.any() else None)
    planning["target_below_practical_minimum"] = {
        "practical_min_cell_size_m": h_min,
        "cells": int(below.sum()),
        "cell_fraction_of_valid": float(below.sum() / max(valid.sum(), 1)),
        "volume_fraction_of_valid": float(
            volumes[below].sum() / volumes[valid].sum()) if valid.any() else 0.0,
        "bounds_m": bounds(centres, below) if centres is not None else None,
        "smallest_target_m": (float(target.values[smallest])
                              if smallest is not None else None),
        "smallest_target_at_m": ([float(v) for v in centres[smallest]]
                                 if smallest is not None and centres is not None
                                 else None),
    }
    budget = int(settings.mesh_budget_cells)
    capped = mesh_estimate(volumes, current, capped_size, valid)
    uncapped = mesh_estimate(volumes, current, uncapped_size, valid)
    planning["mesh_budget"] = {
        "budget_cells": budget,
        "current_cells": int(n_cells),
        "refined_cells_within_practical_minimum": capped,
        "refined_cells_uncapped": uncapped,
        "within_budget": bool(capped <= budget),
        "note": "lower bounds: every valid cell refined to the level that "
                "meets lambda/k (none coarsened), invalid cells unchanged, no "
                "snappyHexMesh buffer layers",
    }
    return planning


def report(document: dict) -> list[str]:
    """The console version: the handful of numbers worth reading first."""
    lines = [f"Rarefaction diagnostics at t = {document['time']}"]
    window = document["sampling"]["window"]
    if window:
        lines.append(f"  sampling window          {window['start_s']:.6g} .. "
                     f"{window['end_s']:.6g} s, {window['averaged_steps']:,} steps")
    else:
        lines.append("  sampling window          not recorded")
    weight = document.get("checks", {}).get("particle_weight")
    if weight and weight.get("checked"):
        verdict = "one number in every occupied cell" if weight["consistent"] \
            else "NOT one number: a field mapping is wrong"
        configured = weight.get("configured_particle_weight")
        lines.append(
            f"  rhoN V / dsmcRhoN        {weight['implied_particle_weight']:.6e}"
            f"   ({verdict}"
            + (f"; nEquivalentParticles {configured:.6e})" if configured
               else "; no constant/dsmcProperties to compare with)"))
    validity = document["validity"]
    if not validity.get("valid_cells") and validity.get("evaluated") is False:
        lines.append(f"  not evaluated: {validity['reason']}")
        return lines
    lines.append(
        f"  valid cells              {validity['valid_cells']:,} of "
        f"{validity['cells']:,} ({validity['valid_cell_fraction']:.1%} of cells, "
        f"{validity['valid_volume_fraction']:.1%} of volume)")
    particles = document["particles"]
    lines.append(f"  validity rule            more than {particles['threshold']:g} "
                 f"sampled parcels ({particles['basis']})")
    if validity["valid_cells"] == 0:
        median = particles.get(particles["basis"], {}).get(
            "median_over_occupied_cells", 0.0)
        lines += [
            f"  ** NO CELL IS VALID: the median occupied cell has {median:.3g}"
            f" sampled parcels ({particles['basis']}), against MORE than"
            f" {particles['threshold']:g}.",
            "     The fields still hold values wherever they can be computed, with"
            " diagValid = 0;",
            "     every statistic below is empty. For valid cells, lower the particle"
            " weight",
            "     (nEquivalentParticles), or relax diagnostics.min_sampled_particles"
            " or particle_count_basis.",
        ]
    for measure, label in (("bounding_box_max", "Delta_max / lambda"),
                           ("volume_equivalent", "V^(1/3) / lambda")):
        block = document["resolution"].get(measure)
        if not block or not block.get("cells"):
            continue
        over = block["exceeding_target"]
        lines.append(
            f"  {label:<24} median {block['median']:.3g}, max {block['max']:.3g};"
            f" above {over['threshold']:.3g} in {over['cell_fraction']:.1%} of"
            f" valid cells, {over['volume_fraction']:.1%} of their volume")
    worst = document["resolution"].get("worst_cells") or []
    if worst:
        x, y, z = worst[0]["at_m"]
        lines.append(f"  worst resolution         {worst[0]['cell_over_mean_free_path']:.3g}"
                     f" mean free paths per cell at ({x:.4g}, {y:.4g}, {z:.4g}) m")
    planning = document.get("planning", {})
    levels = planning.get("refinement_levels")
    if levels:
        lines.append(
            f"  finest permitted level   {levels['finest_permitted_level']} = "
            f"{levels['finest_permitted_cell_size_m'] * 1e3:.4g} mm  (level "
            f"{levels['first_excluded_level']} = "
            f"{levels['first_excluded_cell_size_m'] * 1e3:.4g} mm is below the "
            f"{levels['practical_min_cell_size_m'] * 1e3:.4g} mm minimum)")
    below = planning.get("target_below_practical_minimum")
    if below:
        lines.append(f"  target below minimum     {below['cells']:,} valid cells"
                     f" ({below['volume_fraction_of_valid']:.1%} of valid volume)")
    budget = planning.get("mesh_budget")
    if budget:
        lines.append(
            f"  mesh to reach lambda/k   >= {budget['refined_cells_within_practical_minimum']:.3g}"
            f" cells within the minimum, >= {budget['refined_cells_uncapped']:.3g}"
            f" uncapped (budget {budget['budget_cells']:.3g}, now "
            f"{budget['current_cells']:,})")
    for name, reason in document.get("skipped", {}).items():
        lines.append(f"  skipped {name}: {reason}")
    return lines


def write_yaml(path: Path, document: dict) -> Path:
    """Write the summary, with a provenance comment."""
    import yaml

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# GENERATED by plumetools.rarefaction -- do not edit.\n"
        "# Local mean free path, Knudsen numbers and mesh resolution, from one\n"
        "# sampled time. Statistics are over VALID cells only. See notes.\n\n"
        + yaml.safe_dump(_plain(document), sort_keys=False,
                         default_flow_style=None, width=100),
        encoding="utf-8", newline="\n")
    return path
