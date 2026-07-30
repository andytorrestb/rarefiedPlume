"""Verify a generated mesh against the geometry declared in ``case.yaml``.

Run after ``blockMesh``::

    python -m plumetools.verify_mesh cases/3d-inflow

Checks that need no OpenFOAM, so they also run in CI against a mesh produced
elsewhere. ``checkMesh`` covers mesh *quality*; this covers whether the mesh is
the geometry the model is about to be evaluated on -- which ``checkMesh`` cannot
know.

The most important check is the last one: ``geometry.sphere_radius_m`` is used as
a fixed radius by the source-flow model, so if the meshed sphere differs the
model is silently evaluated on the wrong surface (finding SM-09).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from plumetools import geometry
from plumetools.config import CaseConfig, load_case_config
from plumetools.inflow import read_patch_geometry
from plumetools.mesh.boundary import read_boundary


class MeshVerificationError(AssertionError):
    """A generated mesh does not match the geometry the config declares."""


def verify(case_dir: Path, cfg: CaseConfig | None = None, tol: float = 1e-9) -> list[str]:
    """Check a case's mesh against its config.

    Args:
        case_dir: an OpenFOAM case directory with ``constant/polyMesh``.
        cfg: the case config; loaded from ``case.yaml`` if omitted.
        tol: absolute tolerance for the sphere-radius check [m].

    Returns:
        Human-readable lines describing what was checked.

    Raises:
        MeshVerificationError: on the first check that fails.
    """
    case_dir = Path(case_dir)
    cfg = cfg or load_case_config(case_dir)
    R = cfg.geometry.sphere_radius_m
    report: list[str] = []

    def check(ok: bool, message: str) -> None:
        report.append(("  OK   " if ok else "  FAIL ") + message)
        if not ok:
            raise MeshVerificationError("\n".join(report))

    patches = read_boundary(case_dir)
    report.append(f"case: {case_dir}")
    report.append("patches: " + ", ".join(
        f"{p.name}({p.n_faces}, {p.type})" for p in patches.values()))

    inflow_name = cfg.geometry.patch
    check(inflow_name in patches, f"inflow patch {inflow_name!r} exists")

    # A warning, not a failure. `Unspecified` is a Pointwise V18.5R2 export
    # artifact and is not a registered OpenFOAM patch type (HA-07), but it is a
    # property of a mesh that already exists -- failing here would make this tool
    # useless on the committed mesh, which is exactly where someone would first
    # run it. Generated meshes never emit it, and a unit test enforces that.
    unspecified = [p.name for p in patches.values() if p.type == "Unspecified"]
    if unspecified:
        report.append(
            f"  WARN patches {unspecified} have type 'Unspecified' -- a Pointwise "
            f"export artifact, not a registered OpenFOAM patch type (HA-07). "
            f"Whether v1706 loads it is unresolved; a generated mesh emits 'patch'.")

    labels, faces, points = read_patch_geometry(case_dir, inflow_name)
    sizes = sorted({len(f) for f in faces})

    if cfg.mesh.type == "block_mesh_ogrid":
        nt = int(cfg.mesh.n_tangential)
        expected = 5 * nt * nt
        check(len(faces) == expected,
              f"inflow face count is 5*n_tangential^2 = {expected} (got {len(faces)})")
        check(sizes == [4], "inflow faces are quadrilaterals (blockMesh gives hexahedra)")
    else:
        # snappyHexMesh castellates and snaps, so faces are polygons of mixed
        # vertex count and the count is set by refinement, not by a formula.
        report.append(f"  ..   inflow face vertex counts: {sizes} (snapped surface)")

    used = sorted({p for f in faces for p in f})
    radii = np.linalg.norm(points[used], axis=1)
    deviation = np.abs(radii - R).max()

    # Exact for a projected or arc O-grid -- blockMesh places vertices analytically.
    # Approximate for snappyHexMesh, where points land within the snapping tolerance.
    radius_tol = tol if cfg.mesh.type == "block_mesh_ogrid" else 0.02 * R
    check(deviation < radius_tol,
          f"inflow vertices lie on the sphere r = {R} m to {radius_tol:.3g} m "
          f"(worst deviation {deviation:.3g}, span {radii.min():.6f} .. {radii.max():.6f})")

    # The arc/projection distinction shows up in the face interiors, not the
    # vertices, so it cannot be caught by the radius check above. Measure the
    # centroid deficit instead: an exact hemisphere keeps face centroids close to
    # the sphere, a ruled surface pulls the cap face centre ~16% inside.
    if cfg.mesh.type == "block_mesh_ogrid":
        c_dev = (R - np.linalg.norm(geometry.centroids(faces, points), axis=1)).max() / R
        if cfg.mesh.projection == "searchable_sphere":
            check(c_dev < 0.02,
                  f"projected faces follow the sphere (worst centroid deficit "
                  f"{100 * c_dev:.2f}% of R)")
        else:
            report.append(
                f"  WARN inflow patch is NOT a true hemisphere: worst face-centroid "
                f"deficit {100 * c_dev:.2f}% of R. mesh.projection is "
                f"{cfg.mesh.projection!r}; use 'searchable_sphere' for an exact surface.")

    cents = geometry.centroids(faces, points)
    check(bool(np.all(cents[:, 0] >= -radius_tol)),
          f"inflow is a hemisphere opening toward +x (min centroid x = {cents[:, 0].min():.6f})")

    c_r = np.linalg.norm(cents, axis=1)
    check(bool(np.all(c_r <= R + radius_tol)),
          f"face centroids sit inside the sphere (faceting): "
          f"{c_r.min():.6f} .. {c_r.max():.6f} <= {R}")

    # OpenFOAM winds boundary faces so the normal points OUT of the fluid domain.
    # The fluid is outside the cavity, so an inflow normal points toward the origin.
    norms = geometry.normals(faces, points)
    outward = np.einsum("ij,ij->i", norms, cents / c_r[:, None])
    check(bool(np.all(outward < 0) or np.all(outward > 0)),
          "inflow face normals are consistently oriented")

    report.append(f"  ..   {len(faces)} inflow faces, {len(points)} mesh points")
    return report


def main(argv: list[str]) -> int:
    case_dir = Path(argv[1]) if len(argv) > 1 else Path.cwd()
    try:
        for line in verify(case_dir):
            print(line)
    except MeshVerificationError as exc:
        print(str(exc), file=sys.stderr)
        print("\nmesh verification FAILED", file=sys.stderr)
        return 1
    print("\nmesh verification passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
