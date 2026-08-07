#!/usr/bin/env python3
"""Post-process one case: mesh stats, particle stats, wall pressures, audit.

    python postProcess.py [case_dir]

Writes ``case-summary.json`` -- the machine-readable record ``AllpostCases``
collects into the study table -- and prints a human-readable report.

Everything here is generic across the pressure cases. Nothing is keyed to a
particular reservoir pressure or mesh.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from plumetools.config import load_case_config
from plumetools.markelov1999 import postprocess as pp
from plumetools.markelov1999 import resolution as resolution_module
from plumetools.markelov1999.foamfields import read_patch_field
from plumetools.markelov1999.geometry import from_config


def read_internal_field(path: Path) -> np.ndarray | None:
    """Read the internalField of an OpenFOAM volScalarField.

    Returns ``None`` if the file is absent or uniform -- a uniform ``dsmcRhoN``
    means the solver never wrote per-cell values, which is a fact worth reporting
    rather than a number worth averaging.
    """
    import re

    path = Path(path)
    if not path.is_file():
        return None
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)?\s*\((.*?)\)\s*;",
        text, re.S)
    if not match:
        return None
    return np.asarray(
        [float(v) for v in re.findall(r"[-+0-9.eE]+", match.group(2))],
        dtype=np.float64)


def cell_centres(case_dir: Path) -> np.ndarray | None:
    """Cell centres [m], from a ``writeCellCentres`` output if one exists.

    Returns ``None`` when absent. Computing them here would mean reconstructing
    the cell-face connectivity, which is what OpenFOAM's own utility is for::

        postProcess -func writeCellCentres

    The audit degrades to domain-wide statistics without them rather than
    guessing.
    """
    for name in ("C", "cellCentres"):
        for time_dir in sorted(Path(case_dir).glob("[0-9]*")):
            candidate = time_dir / name
            if candidate.is_file():
                import re
                text = candidate.read_text(encoding="utf-8", errors="replace")
                match = re.search(
                    r"internalField\s+nonuniform\s+List<vector>\s*(\d+)?\s*\((.*?)\n\)\s*;",
                    text, re.S)
                if match:
                    values = [float(v) for v in re.findall(
                        r"[-+0-9.eE]+", match.group(2))]
                    return np.asarray(values, dtype=np.float64).reshape(-1, 3)
    return None


def main(argv) -> int:
    case_dir = Path(argv[1]) if len(argv) > 1 else Path(__file__).resolve().parent
    cfg = load_case_config(case_dir)
    geom = from_config(cfg)

    summary: dict = {}
    existing = case_dir / "case-summary.json"
    if existing.is_file():
        # runInflow.py wrote the inflow, flux and check sections. Keep them:
        # post-processing adds to the record rather than replacing it.
        summary = json.loads(existing.read_text(encoding="utf-8"))

    summary.setdefault("case", dict(cfg.meta))

    print(f"### {case_dir.name}")
    print()

    # 1. mesh statistics
    mesh_stats = pp.mesh_statistics(case_dir)
    summary["mesh"] = mesh_stats.as_dict()
    print("\n".join(mesh_stats.report()))
    print()

    # 2. particle statistics
    particle_stats = pp.parse_solver_log(case_dir / "log.dsmcFoam")
    summary["particles"] = particle_stats.as_dict()
    print("\n".join(particle_stats.report()))
    print()

    # 3-5. wall pressures and the ratio
    try:
        windows, extras = pp.wall_pressure(case_dir, cfg, geom)
        summary["pressure"] = {
            "windows": [w.as_dict() for w in windows],
            **extras,
        }
        print("\n".join(pp.pressure_report(windows, extras)))
    except (FileNotFoundError, ValueError) as exc:
        summary["pressure"] = {"error": str(exc)}
        print(f"Cylinder surface pressure: NOT AVAILABLE\n  {exc}")
    print()

    # post-run resolution audit
    time_dir = pp.latest_time_dir(case_dir)
    counts = read_internal_field(time_dir / "dsmcRhoN") if time_dir else None
    if counts is None:
        summary["resolution_audit"] = {
            "error": "no per-cell dsmcRhoN field found; the audit needs one to "
                     "say what the occupancy actually was"}
        print("Particle resolution AUDIT: NOT AVAILABLE")
        print("  dsmcRhoN holds the literal parcel count per cell, and none was")
        print("  found. Without it the only occupancy figure is the analytical")
        print("  estimate, which is not evidence.")
    else:
        centres = cell_centres(case_dir)
        if centres is None or len(centres) != len(counts):
            print("Particle resolution AUDIT (domain-wide only)")
            print(f"  cells {len(counts)}  mean {counts.mean():.2f}  "
                  f"median {np.median(counts):.2f}  max {counts.max():.0f} "
                  f"particles/cell")
            print("  Per-region occupancy needs cell centres. Generate them with:")
            print("      postProcess -func writeCellCentres")
            summary["resolution_audit"] = {
                "n_cells": int(len(counts)),
                "mean_particles_per_cell": float(counts.mean()),
                "median_particles_per_cell": float(np.median(counts)),
                "max_particles_per_cell": float(counts.max()),
                "total_particles": float(counts.sum()),
                "region_occupancy": {},
                "note": "no cell centres; per-region occupancy not measured",
            }
        else:
            result = resolution_module.audit(
                counts, centres, cfg, geom,
                time=time_dir.name, field="dsmcRhoN")
            summary["resolution_audit"] = {
                "time": result.time,
                "n_cells": result.n_cells,
                "mean_particles_per_cell": result.mean_particles_per_cell,
                "median_particles_per_cell": result.median_particles_per_cell,
                "max_particles_per_cell": result.max_particles_per_cell,
                "occupied_cells": result.occupied_cells,
                "total_particles": result.total_particles,
                "region_occupancy": result.region_occupancy,
                "target": result.target,
                "meets_target": result.meets_target(),
            }
            print("\n".join(result.report()))
    print()

    path = case_dir / "case-summary.json"
    path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
