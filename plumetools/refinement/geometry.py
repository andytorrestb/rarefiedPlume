"""Actual cell edges and source identities using shared OpenFOAM readers."""
from __future__ import annotations

import hashlib
from pathlib import Path
import re
from typing import Optional

import numpy as np

from plumetools.foamio.volfield import read_internal
from plumetools.mesh.polymesh import read_faces_flat, read_label_list, read_points_array
from .model import Geometry

MESH_FILES = ("points", "faces", "owner", "neighbour", "boundary")
HISTORY_FILES = ("cellLevel", "pointLevel", "refinementHistory", "level0Edge")
DIMS = {"cell_centres": (0, 1, 0, 0, 0, 0, 0), "cell_volumes": (0, 3, 0, 0, 0, 0, 0)}


def file_sha256(path: Path) -> str:
    """Stream a content hash without loading a second copy of a field."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_ascii(path: Path):
    """Reject binary before the shared ASCII parser reads the payload."""
    if not path.is_file():
        raise ValueError(f"missing {path}")
    with path.open("rb") as stream:
        header = stream.read(8192).decode("ascii", errors="replace")
    if re.search(r"\bformat\s+binary\s*;", header):
        raise ValueError(f"{path}: ASCII input required; convert a copy with foamFormatConvert")
    return header


def field_values(path: Path, n_cells: int, dimensions, *, vector=False):
    """Read a dimension-checked field with declared uniform broadcasting only."""
    header = require_ascii(path)
    match = re.search(r"\bdimensions\s*\[([^\]]+)\]\s*;", header)
    if not match:
        raise ValueError(f"{path}: missing dimensions")
    try:
        dims = tuple(float(v) for v in match.group(1).split())
    except ValueError:
        raise ValueError(f"{path}: malformed dimensions") from None
    if dims != tuple(dimensions):
        raise ValueError(f"{path}: dimensions {dims} do not match {tuple(dimensions)}")
    data = read_internal(path).over(n_cells, "cells in the mesh")
    shape = (n_cells, 3) if vector else (n_cells,)
    if data.shape != shape:
        raise ValueError(f"{path}: expected field shape {shape}, got {data.shape}")
    return data


def mesh_files(case: Path) -> dict:
    """Relative mesh inputs, including optional hex refinement history."""
    poly = Path(case) / "constant/polyMesh"
    if not poly.is_dir():
        if any(Path(case).glob("processor*")):
            raise ValueError("reconstructed mesh required; run reconstructParMesh/reconstructPar on a copy")
        raise ValueError(f"{poly}: no static mesh")
    found = {}
    for name in MESH_FILES + HISTORY_FILES:
        path = poly / name
        if name in MESH_FILES:
            require_ascii(path)
        if path.is_file():
            found[path.relative_to(case).as_posix()] = file_sha256(path)
    return found


def mesh_fingerprint(case: Path) -> str:
    """Content identity independent of absolute location and unused cell sets."""
    digest = hashlib.sha256()
    for name, value in sorted(mesh_files(Path(case)).items()):
        digest.update((name + "\0" + value + "\n").encode())
    return digest.hexdigest()


def read_geometry(case: Path, *, time: str, field_names: Optional[dict] = None) -> Geometry:
    """Read physical mesh data; configurable C/V names support outside cases."""
    case = Path(case)
    mesh_files(case)
    if (case / time / "polyMesh").exists():
        raise ValueError("time-specific/dynamic mesh not supported; provide a static reconstructed snapshot")
    dynamic = case / "constant/dynamicMeshDict"
    if dynamic.is_file():
        match = re.search(r"\bdynamicFvMesh\s+(\w+)\s*;", dynamic.read_text(errors="replace"))
        if match and match.group(1) != "staticFvMesh":
            raise ValueError("dynamic mesh input requires a static exported snapshot")
    points = read_points_array(case)
    offsets, labels = read_faces_flat(case)
    owner = read_label_list(case / "constant/polyMesh/owner")
    neighbour = read_label_list(case / "constant/polyMesh/neighbour")
    n_faces = offsets.size - 1
    if owner.size != n_faces or neighbour.size > n_faces or np.any(owner < 0) or np.any(neighbour < 0):
        raise ValueError("mesh owner/neighbour do not match face topology")
    if np.any(np.diff(offsets) < 3) or np.any(labels < 0) or np.any(labels >= len(points)):
        raise ValueError("mesh has invalid face vertices")
    n_cells = int(max(owner.max(initial=-1), neighbour.max(initial=-1))) + 1
    max_edge = np.zeros(n_cells)
    for first in range(0, n_faces, 65536):
        last = min(first + 65536, n_faces)
        local = offsets[first:last+1] - offsets[first]
        vertices = labels[offsets[first]:offsets[last]]
        following = np.arange(vertices.size) + 1
        following[local[1:]-1] = local[:-1]
        lengths = np.linalg.norm(points[vertices] - points[vertices[following]], axis=1)
        face_max = np.maximum.reduceat(lengths, local[:-1])
        np.maximum.at(max_edge, owner[first:last], face_max)
        internal_last = min(last, neighbour.size)
        if first < internal_last:
            np.maximum.at(max_edge, neighbour[first:internal_last], face_max[:internal_last-first])
    names = {"cell_centres": "C", "cell_volumes": "V", **(field_names or {})}
    arrays = {}
    for role, command in (("cell_centres", "writeCellCentres"), ("cell_volumes", "writeCellVolumes")):
        path = case / time / names[role]
        if not path.is_file():
            raise ValueError(f"missing {path}; postProcess -func {command} -time {time} exports it")
        arrays[role] = field_values(path, n_cells, DIMS[role], vector=(role == "cell_centres"))
    return Geometry(arrays["cell_centres"], arrays["cell_volumes"], max_edge, mesh_fingerprint(case))
