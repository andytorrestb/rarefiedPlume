# `plumetools.mesh`

`constant/polyMesh` parsing: the `boundary` file, points and faces (whole, by
range, by label — streamed), owner/neighbour, the cell count, face sets, and
the geometry of boundary patches (centres, area vectors, normals). Nothing else
in the package parses `polyMesh` itself.

Depends on: `numpy`. Used by every family that reads a built mesh.

::: plumetools.mesh
    options:
      show_root_heading: false
      members: false

## Modules

| Module | Purpose | Principal interfaces |
|---|---|---|
| [`plumetools.mesh.boundary`](boundary.md) | `constant/polyMesh/boundary` | `read_boundary`, `require_patch`, `PatchInfo` |
| [`plumetools.mesh.polymesh`](polymesh.md) | `points`, `faces`, labels, the cell count, streamed | `read_points`, `read_faces`, `read_face_range`, `read_face_ranges`, `read_faces_at`, `read_patch_faces`, `read_label_list`, `read_cell_count`, `read_points_array`, `read_faces_flat` |
| [`plumetools.mesh.patches`](patches.md) | boundary-patch geometry | `read_patches`, `face_geometry`, `PatchGeometry` |
| [`plumetools.mesh.sets`](sets.md) | `constant/polyMesh/sets/<name>` face sets | `read_face_set` |
| [`plumetools.mesh._parse`](_parse.md) | shared ASCII list-file helpers (internal) | `read_text`, `strip_comments`, `find_list_block` |

`read_boundary(case_dir)` returns `{name: PatchInfo}` with each patch's type,
face count and start face; `read_patches(case_dir, names)` returns a
`PatchGeometry` per named patch, with face centres, area vectors and normals as
arrays. Faces are read in one pass, keeping only the requested labels, which is
what made the 2026-10-27 refactor's memory reductions (44 MB `faces` file: 592
to 250 MiB for the nozzle patch, recorded in the
[refactor record](../../../refactor-2026-10-27.md#evidence)).

Two geometry routines are deliberately **not** replaced by
`patches.face_geometry`: `cai2012ext.interface_case.openfoam_normals` and the
legacy `plumetools.geometry` centroids — a different summation order could
change written fields in the last bit, which the byte-for-byte pins would catch.
