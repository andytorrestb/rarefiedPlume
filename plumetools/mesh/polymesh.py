"""Parse ``constant/polyMesh/points`` and ``constant/polyMesh/faces``.

Two families of reader live here. :func:`read_points` and :func:`read_faces`
walk the file a line at a time and return Python lists -- fine for the boundary
patches and small meshes they were written for. :func:`read_points_array`,
:func:`read_faces_flat` and :func:`read_label_list` return flat numpy arrays
and do no per-line Python work, which is what a whole-mesh pass over a
multi-million-cell octree needs (see :mod:`plumetools.rarefaction.foam`).
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from plumetools.mesh._parse import find_list_block, read_text, strip_comments

#: OpenFOAM writes a list of ten entries or fewer on one line: ``3(0 1 2)``.
_INLINE_LIST = re.compile(r"^\s*(\d+)\s*\((.*)\)\s*$", re.S)

#: Turns every parenthesis into a space, so a payload parses as one stream.
_NO_PARENS = str.maketrans("()", "  ")


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


# --------------------------------------------------------------------------- #
# whole-mesh readers: flat numpy arrays, no per-line Python
# --------------------------------------------------------------------------- #

def _list_payload(path: Path) -> tuple[int, str]:
    """The declared count and the raw payload of an ASCII OpenFOAM list file.

    Handles both forms OpenFOAM writes: the multi-line one (``N``, then ``(``,
    one entry per line, ``)``) and, for ten entries or fewer, the inline
    ``N(...)``. Comments are stripped from the header only -- see
    :func:`read_face_range` for why.

    Raises:
        ValueError: for a binary file, or one with no list in it.
    """
    path = Path(path)
    text = read_text(path)
    opening = text.find("\n(\n")
    if opening >= 0:
        count, _ = find_list_block(strip_comments(text[:opening + 3]), path)
        closing = text.rfind(")")
        if closing <= opening:
            raise ValueError(f"{path}: unterminated list")
        return count, text[opening + 3:closing]

    body = strip_comments(text)
    header_end = body.find("}")            # the FoamFile dictionary
    # find_list_block's binary check applies to the inline form too.
    if re.search(r"^\s*format\s+binary\s*;", body[:header_end + 1], re.M):
        find_list_block(body, path)
    match = _INLINE_LIST.match(body[header_end + 1:].strip())
    if not match:
        raise ValueError(f"{path}: no '<count> ( ... )' list block found")
    return int(match.group(1)), match.group(2)


def read_label_list(path: Path) -> np.ndarray:
    """A ``labelList`` file -- ``owner`` or ``neighbour`` -- as an int64 array.

    Raises:
        ValueError: if the number of labels is not the declared count.
    """
    count, payload = _list_payload(path)
    labels = np.fromstring(payload, dtype=np.int64, sep=" ") if payload.strip() \
        else np.empty(0, dtype=np.int64)
    if labels.size != count:
        raise ValueError(f"{path}: declares {count} labels, holds {labels.size}")
    return labels


def read_points_array(case_dir: Path) -> np.ndarray:
    """``constant/polyMesh/points`` as an ``(n_points, 3)`` float64 array.

    The same values as :func:`read_points`, without a Python loop per point.
    """
    path = Path(case_dir) / "constant" / "polyMesh" / "points"
    count, payload = _list_payload(path)
    values = np.fromstring(payload.translate(_NO_PARENS), dtype=np.float64, sep=" ")
    if values.size != 3 * count:
        raise ValueError(
            f"{path}: declares {count} points, holds {values.size / 3:g}")
    return values.reshape(count, 3)


def read_faces_flat(case_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    """``constant/polyMesh/faces`` as CSR arrays.

    Returns:
        ``(offsets, labels)``: face ``i`` is ``labels[offsets[i]:offsets[i+1]]``,
        in winding order. ``offsets`` has ``n_faces + 1`` entries.

    The vertex count of each face is read from the digits just before its
    ``(`` -- one or two of them for any face a mesher makes -- so finding where
    each face starts needs no sequential walk. Every count is then checked
    against the label stream it was read from.
    """
    path = Path(case_dir) / "constant" / "polyMesh" / "faces"
    count, payload = _list_payload(path)
    raw = np.frombuffer(payload.encode("ascii"), dtype=np.uint8)
    opens = np.flatnonzero(raw == ord("("))
    if opens.size != count:
        raise ValueError(f"{path}: declares {count} faces, holds {opens.size}")

    # The count before each '(' has at most three digits; a non-digit (newline,
    # space, the start of the payload) ends it.
    padded = np.concatenate([np.full(3, ord(" "), dtype=np.uint8), raw])
    sizes = np.zeros(count, dtype=np.int64)
    place = np.ones(count, dtype=np.int64)
    going = np.ones(count, dtype=bool)
    for back in (1, 2, 3):
        digit = padded[opens + 3 - back].astype(np.int64) - ord("0")
        going &= (digit >= 0) & (digit <= 9)
        sizes += np.where(going, digit * place, 0)
        place *= 10
    if np.any(sizes < 3):
        raise ValueError(f"{path}: a face has fewer than 3 vertices, or no count")

    tokens = np.fromstring(payload.translate(_NO_PARENS), dtype=np.int64, sep=" ")
    starts = np.concatenate([[0], np.cumsum(sizes + 1)[:-1]])
    if tokens.size != int((sizes + 1).sum()) or not np.array_equal(tokens[starts], sizes):
        raise ValueError(f"{path}: face vertex counts do not match the labels")
    is_label = np.ones(tokens.size, dtype=bool)
    is_label[starts] = False
    offsets = np.concatenate([[0], np.cumsum(sizes)])
    return offsets, tokens[is_label]
