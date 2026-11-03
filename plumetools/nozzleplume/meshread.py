"""Compatibility names: boundary-patch geometry now lives in :mod:`plumetools.mesh`.

:class:`PatchGeometry`, :func:`face_geometry` and :func:`read_patches` are in
:mod:`plumetools.mesh.patches`, and :func:`read_cell_count` in
:mod:`plumetools.mesh.polymesh`, so that code outside this family -- the
rarefaction diagnostics -- reads mesh geometry without going through it. The
names below are the same objects, for callers written against this module.
"""

from plumetools.mesh.patches import PatchGeometry, face_geometry, read_patches
from plumetools.mesh.polymesh import read_cell_count

__all__ = ["PatchGeometry", "face_geometry", "read_cell_count", "read_patches"]
