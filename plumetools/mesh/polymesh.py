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


def read_face_range(case_dir: Path, start_face: int, n_faces: int) -> list[list[int]]:
    """Read a contiguous range of faces, parsing **only** that range.

    Args:
        case_dir: an OpenFOAM case directory.
        start_face: global label of the first face wanted.
        n_faces: how many faces to read.

    Returns:
        ``n_faces`` point-label lists, in winding order.

    Raises:
        ValueError: if the range runs past the end of the file, or a line in it
            is malformed.

    The same answer as :func:`read_patch_faces`, without converting the rest of
    the file. Boundary faces are a few percent of a mesh, and on a mesh of a few
    million faces converting all of them to reach one patch takes longer than
    everything else a caller does with the result.

    Comments are stripped from the **header only** -- everything up to the
    list's opening parenthesis. :func:`strip_comments` walks its input a
    character at a time, which on a 50 MB ``faces`` file is most of the cost of
    reading it (8.5 s of 10.7 s, measured on a 530 000-cell mesh), and the
    payload of a file OpenFOAM wrote carries no comments to strip. A line in the
    requested range that is not a face is still an error, so a payload that did
    carry one would stop the read rather than be mis-parsed.
    """
    path = Path(case_dir) / "constant" / "polyMesh" / "faces"
    text = read_text(path)
    opening = text.find("\n(\n")
    if opening < 0:
        raise ValueError(f"{path}: no '<count> ( ... )' list block found")
    count, _ = find_list_block(strip_comments(text[:opening + 3]), path)
    if start_face < 0 or n_faces < 0 or start_face + n_faces > count:
        raise ValueError(
            f"{path}: faces {start_face} .. {start_face + n_faces} requested, but "
            f"the mesh has {count}")
    lines = text[opening + 3:].split("\n")

    faces: list[list[int]] = []
    for i in range(start_face, start_face + n_faces):
        raw = lines[i].strip()
        open_paren = raw.find("(")
        close_paren = raw.rfind(")")
        if open_paren < 0 or close_paren < open_paren:
            raise ValueError(
                f"{path}: face {i}: expected 'N(a b c ...)', got {raw!r}")
        faces.append([int(v) for v in raw[open_paren + 1:close_paren].split()])
    return faces
