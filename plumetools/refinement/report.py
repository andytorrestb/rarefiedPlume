"""Portable JSON/CSV reporting; no plotting or case-specific dependency."""
from __future__ import annotations

from dataclasses import asdict
import csv
import json
from pathlib import Path

import numpy as np

from .plan import json_safe


def write_cells_csv(path, geometry, *, target, density, selected=None, measured=None, assessment=None):
    """Export source cell identities/geometry and targets for another mesher."""
    chosen = set() if selected is None else set(np.asarray(selected).tolist())
    names = ["cell_id", "x_m", "y_m", "z_m", "volume_m3", "max_edge_m", "target_edge_m", "number_density_estimate", "selected"]
    if measured is not None:
        names += ["temperature_K", "mean_parcels", "lambda_measured_m", "mesh_status", "occupancy_status"]
    with Path(path).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(names)
        for i, (xyz, volume, edge, limit, n) in enumerate(zip(geometry.centres_m, geometry.volumes_m3, geometry.max_edge_m, target, density)):
            row = [i, *map(float, xyz), float(volume), float(edge), float(limit), float(n), int(i in chosen)]
            if measured is not None:
                row += [float(measured.temperature_K[i]), float(measured.mean_parcels[i]),
                        float(assessment.lambda_m[i]), int(assessment.mesh_status[i]), int(assessment.occupancy_status[i])]
            writer.writerow(row)


def write_assessment(geometry, measured, assessment, directory: Path) -> Path:
    """Write an explicitly selected output directory, leaving case fields alone."""
    directory = Path(directory)
    document = json_safe({"kind": "plumetools.refinement.assessment", "schema_version": 1,
                          "mesh_sha256": geometry.mesh_sha256, "source_sha256": measured.source_sha256,
                          "targets": asdict(assessment.targets), "summary": assessment.summary,
                          "window": measured.window, "provenance": measured.provenance})
    text = json.dumps(document, indent=2, allow_nan=False) + "\n"
    if directory.exists() and any(directory.iterdir()):
        existing = directory / "assessment.json"
        if existing.is_file() and existing.read_text(encoding="utf-8") == text and (directory / "cells.csv").is_file():
            return existing
        raise ValueError(f"output already exists: {directory}; choose a fresh directory")
    directory.mkdir(parents=True, exist_ok=True)
    write_cells_csv(directory / "cells.csv", geometry, target=assessment.target_edge_m,
                    density=measured.number_density, selected=np.flatnonzero(assessment.mesh_status == 0),
                    measured=measured, assessment=assessment)
    path = directory / "assessment.json"
    path.write_text(text, encoding="utf-8")
    return path
