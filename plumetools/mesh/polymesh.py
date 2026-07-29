"""Parse ``constant/polyMesh/points`` and ``constant/polyMesh/faces``."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from plumetools.mesh._parse import find_list_block, read_text, strip_comments


def read_points(case_dir: Path) -> np.ndarray:
    """Read mesh point coordinates.

    Args:
        case_dir: an OpenFOAM case directory.

    Returns:
        ``(n_points, 3)`` float64 array of coordinates [m], indexed by point label.

    Raises:
        ValueError: if the list structure or a coordinate triple is malformed.
    """
    path = Path(case_dir) / "constant" / "polyMesh" / "points"
    text = strip_comments(read_text(path))
    count, first = find_list_block(text, path)
    lines = text.split("\n")

    points = np.empty((count, 3), dtype=np.float64)
    for i in range(count):
        raw = lines[first + i].strip()
        if not (raw.startswith("(") and raw.endswith(")")):
            raise ValueError(f"{path}: line {first + i + 1}: expected '(x y z)', got {raw!r}")
        parts = raw[1:-1].split()
        if len(parts) != 3:
            raise ValueError(f"{path}: line {first + i + 1}: expected 3 coordinates, got {len(parts)}")
        points[i] = [float(p) for p in parts]
    return points


def read_faces(case_dir: Path) -> list[list[int]]:
    """Read face-to-point connectivity.

    Args:
        case_dir: an OpenFOAM case directory.

    Returns:
        List of ``n_faces`` point-label lists, in winding order, indexed by global
        face label.

    Faces of any vertex count are accepted -- ``3(...)`` for the tetrahedral
    Pointwise meshes, ``4(...)`` for anything ``blockMesh`` produces, and larger
    polygons. The legacy reader loaded the whole file into a string-keyed dict and
    then parsed it a second time into an identical dict (finding GP-03).
    """
    path = Path(case_dir) / "constant" / "polyMesh" / "faces"
    text = strip_comments(read_text(path))
    count, first = find_list_block(text, path)
    lines = text.split("\n")

    faces: list[list[int]] = []
    for i in range(count):
        raw = lines[first + i].strip()
        open_paren = raw.find("(")
        close_paren = raw.rfind(")")
        if open_paren < 0 or close_paren < open_paren:
            raise ValueError(f"{path}: line {first + i + 1}: expected 'N(a b c ...)', got {raw!r}")
        labels = [int(v) for v in raw[open_paren + 1:close_paren].split()]
        declared = raw[:open_paren].strip()
        if declared and int(declared) != len(labels):
            raise ValueError(
                f"{path}: line {first + i + 1}: declared {declared} vertices, found {len(labels)}"
            )
        faces.append(labels)
    return faces


def read_patch_faces(case_dir: Path, start_face: int, n_faces: int) -> list[list[int]]:
    """Read only the faces belonging to one patch.

    Patch faces are contiguous in ``constant/polyMesh/faces``, beginning at the
    patch's ``startFace``.
    """
    return read_faces(case_dir)[start_face:start_face + n_faces]
