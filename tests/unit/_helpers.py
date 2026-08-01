"""Builders for synthetic OpenFOAM meshes used by the unit tests."""

from __future__ import annotations

import textwrap
from pathlib import Path

HEADER = """\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "constant/polyMesh";
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""

#: A unit square in the z = 0 plane, extended by one column of points.
SQUARE_POINTS = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (2, 0, 0), (2, 1, 0)]
TRI_FACES = [[0, 1, 2], [0, 2, 3]]
QUAD_FACES = [[0, 1, 2, 3], [1, 4, 5, 2]]
MIXED_FACES = [[0, 1, 2], [0, 2, 3], [1, 4, 5, 2]]


def write(path: Path, text: str, newline: str = "\n") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline=newline) as f:
        f.write(text)


def patch(name, n_faces, start_face, ptype="patch", extra=""):
    """One boundary entry. The body is written verbatim so a test can control
    key order and spelling."""
    return (name, f"type {ptype};\n{extra}nFaces {n_faces};\nstartFace {start_face};\n")


def make_case(root: Path, *, points, faces, patches, newline: str = "\n",
              sets: dict | None = None) -> Path:
    """Build a minimal polyMesh.

    Args:
        root: case directory to create.
        points: iterable of (x, y, z).
        faces: iterable of point-label lists.
        patches: iterable of ``(name, body_text)``.
        newline: line terminator, for exercising CRLF handling.
        sets: optional ``{set_name: [face labels]}``.
    """
    points, faces, patches = list(points), list(faces), list(patches)
    mesh = Path(root) / "constant" / "polyMesh"

    pts = "\n".join(f"({x} {y} {z})" for x, y, z in points)
    write(mesh / "points",
          HEADER.format(cls="vectorField", obj="points")
          + f"\n{len(points)}\n(\n{pts}\n)\n", newline)

    fcs = "\n".join(f"{len(f)}({' '.join(str(i) for i in f)})" for f in faces)
    write(mesh / "faces",
          HEADER.format(cls="faceList", obj="faces")
          + f"\n{len(faces)}\n(\n{fcs}\n)\n", newline)

    body = "\n".join(
        f"    {name}\n    {{\n"
        + textwrap.indent(text.strip("\n"), "        ")
        + "\n    }"
        for name, text in patches
    )
    write(mesh / "boundary",
          HEADER.format(cls="polyBoundaryMesh", obj="boundary")
          + f"\n{len(patches)}\n(\n{body}\n)\n", newline)

    for name, labels in (sets or {}).items():
        joined = "\n".join(str(i) for i in labels)
        write(mesh / "sets" / name,
              HEADER.format(cls="faceSet", obj=name)
              + f"\n{len(labels)}\n(\n{joined}\n)\n", newline)
    return Path(root)
