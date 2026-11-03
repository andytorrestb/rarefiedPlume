"""Whole-mesh geometry from ``constant/polyMesh``: faces, cells, patches, adjacency.

The breakdown diagnostics, the interface extraction and the budgets all need
the same things -- every face's centre and area vector, every cell's centre and
volume, owner/neighbour, and which patch each boundary face is on -- for the
actual 3-D cells. Reading goes through :mod:`plumetools.mesh` (the shared
streamed readers); this module only does geometry.

Face and cell geometry follow OpenFOAM's ``primitiveMesh`` algorithms:

* a face is split into triangles about its vertex average; the area vector is
  their sum and the centre their area-weighted centroid (exact for planar
  polygons, identical to the vertex average for rectangles);
* a cell is split into pyramids from the average of its face centres; the
  volume is ``sum(S_f . (C_f - c_est))/3`` and the centre is the pyramid-volume
  weighted centroid ``(3/4) C_f + (1/4) c_est``.

Area vectors point from owner to neighbour, and out of the domain on boundary
faces, exactly as OpenFOAM stores them.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from plumetools.mesh.boundary import read_boundary
from plumetools.mesh.polymesh import read_faces_flat, read_label_list, read_points_array


@dataclass(frozen=True)
class MeshGeometry:
    points: np.ndarray
    offsets: np.ndarray            # CSR into labels, n_faces + 1
    labels: np.ndarray
    owner: np.ndarray              # (n_faces,)
    neighbour: np.ndarray          # (n_internal,)
    face_centres: np.ndarray
    face_areas: np.ndarray         # area vectors
    cell_centres: np.ndarray
    cell_volumes: np.ndarray
    patches: dict                  # name -> PatchInfo

    @property
    def n_cells(self) -> int:
        return int(self.cell_volumes.size)

    @property
    def n_faces(self) -> int:
        return int(self.owner.size)

    @property
    def n_internal(self) -> int:
        return int(self.neighbour.size)

    def patch_slice(self, name: str) -> slice:
        info = self.patches[name]
        return slice(info.start_face, info.start_face + info.n_faces)

    def face_patch(self) -> np.ndarray:
        """Patch index of every face; -1 for internal faces."""
        out = np.full(self.n_faces, -1, dtype=np.int64)
        for i, name in enumerate(self.patches):
            out[self.patch_slice(name)] = i
        return out

    def face_vertices(self, face: int) -> np.ndarray:
        return self.points[self.labels[self.offsets[face]:self.offsets[face + 1]]]

    def cell_size(self) -> np.ndarray:
        """``(n_cells, 3)`` axis-aligned bounding-box extents [m]."""
        from plumetools.rarefaction.foam import cell_extents
        return cell_extents(self.points, self.offsets, self.labels, self.owner,
                            self.neighbour, self.n_cells)

    def cell_graph(self):
        """:class:`plumetools.rarefaction.core.CellGraph` for least-squares gradients."""
        from plumetools.rarefaction.foam import cell_graph
        return cell_graph(self.cell_centres, self.points, self.offsets, self.labels,
                          self.owner, self.neighbour, self.patches)


def face_geometry_flat(points, offsets, labels) -> tuple:
    """Centres and area vectors of every face of a CSR face list."""
    n_faces = offsets.size - 1
    sizes = np.diff(offsets)
    centres = np.empty((n_faces, 3))
    areas = np.empty((n_faces, 3))
    for size in np.unique(sizes):
        idx = np.nonzero(sizes == size)[0]
        verts = points[labels[offsets[idx][:, None] + np.arange(size)[None, :]]]
        c_est = verts.mean(axis=1)
        nxt = np.roll(verts, -1, axis=1)
        tri_n = 0.5 * np.cross(verts - c_est[:, None, :], nxt - c_est[:, None, :])
        tri_c = (verts + nxt + c_est[:, None, :]) / 3.0
        s = tri_n.sum(axis=1)
        smag = np.linalg.norm(s, axis=1)
        unit = s / np.where(smag > 0, smag, 1.0)[:, None]
        w = np.einsum("fti,fi->ft", tri_n, unit)
        wsum = w.sum(axis=1)
        centres[idx] = np.einsum("ft,fti->fi", w, tri_c) / np.where(
            wsum > 0, wsum, 1.0)[:, None]
        areas[idx] = s
    return centres, areas


def cell_geometry(face_centres, face_areas, owner, neighbour, n_cells) -> tuple:
    n_int = neighbour.size
    count = np.bincount(owner, minlength=n_cells) + np.bincount(neighbour, minlength=n_cells)
    c_est = np.zeros((n_cells, 3))
    for a in range(3):
        c_est[:, a] = (np.bincount(owner, face_centres[:, a], minlength=n_cells)
                       + np.bincount(neighbour, face_centres[:n_int, a], minlength=n_cells))
    c_est /= count[:, None]

    pyr_o = np.einsum("fi,fi->f", face_areas, face_centres - c_est[owner])
    pyr_n = -np.einsum("fi,fi->f", face_areas[:n_int],
                       face_centres[:n_int] - c_est[neighbour])
    vol3 = np.bincount(owner, pyr_o, minlength=n_cells) + np.bincount(
        neighbour, pyr_n, minlength=n_cells)
    centres = np.zeros((n_cells, 3))
    for a in range(3):
        pc_o = 0.75 * face_centres[:, a] + 0.25 * c_est[owner, a]
        pc_n = 0.75 * face_centres[:n_int, a] + 0.25 * c_est[neighbour, a]
        centres[:, a] = (np.bincount(owner, pyr_o * pc_o, minlength=n_cells)
                         + np.bincount(neighbour, pyr_n * pc_n, minlength=n_cells))
    centres /= vol3[:, None]
    return centres, vol3 / 3.0


def read_mesh(case_dir: Path) -> MeshGeometry:
    case_dir = Path(case_dir)
    poly = case_dir / "constant" / "polyMesh"
    points = read_points_array(case_dir)
    offsets, labels = read_faces_flat(case_dir)
    owner = read_label_list(poly / "owner").astype(np.int64)
    neighbour = read_label_list(poly / "neighbour").astype(np.int64)
    n_cells = int(max(owner.max(initial=-1), neighbour.max(initial=-1))) + 1
    fc, fa = face_geometry_flat(points, offsets.astype(np.int64), labels.astype(np.int64))
    cc, cv = cell_geometry(fc, fa, owner, neighbour, n_cells)
    return MeshGeometry(points, offsets.astype(np.int64), labels.astype(np.int64), owner,
                        neighbour, fc, fa, cc, cv, read_boundary(case_dir))
