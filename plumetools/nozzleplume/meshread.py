"""Read boundary-patch geometry out of a built ``constant/polyMesh``.

Shared by :mod:`plumetools.nozzleplume.verify`, which measures the mesh against
``case.yaml``, and :mod:`plumetools.nozzleplume.post`, which integrates the
plate loads over it. Both need the same three things per patch -- where each
face is, how big it is, and which way it faces -- and nothing about the interior
of the mesh.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from plumetools.mesh.boundary import read_boundary
from plumetools.mesh.polymesh import read_face_range, read_points


@dataclass(frozen=True)
class PatchGeometry:
    """The faces of one boundary patch.

    Attributes:
        name: the patch name.
        centres: ``(n, 3)`` face centres [m] -- the vertex average, which is the
            centroid for the rectangles an octree mesh is made of.
        area_vectors: ``(n, 3)`` face area vectors [m^2]. OpenFOAM winds a
            boundary face so its normal points **out of the fluid**: into the
            plate on a wall, and towards ``-x`` on the nozzle.
        vertices: ``(m, 3)`` the distinct mesh points the patch uses [m].
    """

    name: str
    centres: np.ndarray
    area_vectors: np.ndarray
    vertices: np.ndarray

    @property
    def n_faces(self) -> int:
        return int(self.centres.shape[0])

    @property
    def areas(self) -> np.ndarray:
        """``(n,)`` face areas [m^2]."""
        return np.linalg.norm(self.area_vectors, axis=1)

    @property
    def normals(self) -> np.ndarray:
        """``(n, 3)`` unit normals, pointing out of the fluid."""
        return self.area_vectors / self.areas[:, None]

    @property
    def area_m2(self) -> float:
        """Total patch area [m^2]."""
        return float(self.areas.sum())


def face_geometry(faces: list, points: np.ndarray) -> tuple:
    """Centres and area vectors of many polygonal faces.

    Args:
        faces: point-label lists, one per face, in winding order.
        points: ``(n_points, 3)`` mesh coordinates [m].

    Returns:
        ``(centres, area_vectors)``, each ``(n_faces, 3)``.

    The area vector is the triangle-fan sum from the first vertex, exact for any
    planar polygon. Faces are grouped by vertex count so each group is one array
    operation: a plate patch is tens of thousands of quads, and a loop over them
    is most of the time this would otherwise take.
    """
    centres = np.empty((len(faces), 3), dtype=np.float64)
    area_vectors = np.empty((len(faces), 3), dtype=np.float64)

    by_size: dict = {}
    for index, face in enumerate(faces):
        by_size.setdefault(len(face), []).append(index)

    for size, indices in by_size.items():
        if size < 3:
            raise ValueError(f"a face needs at least 3 vertices, got {size}")
        vertices = points[np.asarray([faces[i] for i in indices], dtype=np.int64)]
        fan = np.cross(vertices[:, 1:-1] - vertices[:, :1],
                       vertices[:, 2:] - vertices[:, :1])
        centres[indices] = vertices.mean(axis=1)
        area_vectors[indices] = 0.5 * fan.sum(axis=1)
    return centres, area_vectors


def read_patches(case_dir: Path, names) -> dict:
    """Read the geometry of the named patches.

    Args:
        case_dir: an OpenFOAM case directory with ``constant/polyMesh``.
        names: the patch names wanted.

    Returns:
        ``{name: PatchGeometry}``, in the order asked for.

    Raises:
        KeyError: for a patch the mesh does not have, listing the ones it does.

    The points are read once and shared; the faces are read patch by patch
    through :func:`plumetools.mesh.polymesh.read_face_range`, which does not
    convert the interior faces it is not asked for.
    """
    case_dir = Path(case_dir)
    boundary = read_boundary(case_dir)
    absent = sorted(name for name in names if name not in boundary)
    if absent:
        raise KeyError(
            f"patch(es) {absent} not in the mesh; it has {sorted(boundary)}")

    points = read_points(case_dir)
    out = {}
    for name in names:
        info = boundary[name]
        faces = read_face_range(case_dir, info.start_face, info.n_faces)
        centres, area_vectors = face_geometry(faces, points)
        used = np.unique(np.asarray(
            [label for face in faces for label in face], dtype=np.int64))
        out[name] = PatchGeometry(
            name=name, centres=centres, area_vectors=area_vectors,
            vertices=points[used] if used.size else np.empty((0, 3)))
    return out


_N_CELLS = re.compile(r"nCells\s*:\s*(\d+)")


def read_cell_count(case_dir: Path) -> int | None:
    """The mesh's cell count, from the ``note`` in ``polyMesh/owner``.

    OpenFOAM writes ``note "nPoints:... nCells:... nFaces:..."`` into the header
    of ``owner``. Reading it there costs nothing; counting distinct owners would
    mean parsing a file with one line per face. ``None`` if the note is absent,
    which a mesh from another tool may be.
    """
    path = Path(case_dir) / "constant" / "polyMesh" / "owner"
    if not path.is_file():
        return None
    with open(path, encoding="utf-8", errors="replace") as handle:
        head = handle.read(4096)
    match = _N_CELLS.search(head)
    return int(match.group(1)) if match else None
