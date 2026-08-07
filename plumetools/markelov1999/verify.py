r"""Verify a generated mesh against the geometry ``case.yaml`` declares.

Run after ``snappyHexMesh``::

    python -m plumetools.markelov1999.verify cases/markelov1999/Cases/gap06in/p005psi

``checkMesh`` answers "is this a valid mesh?". This answers a question
``checkMesh`` cannot: **is this the geometry the model is about to be evaluated
on?** A mesh can be flawless and still have the cylinder in the wrong place, the
plate at the wrong gap, or an inflow surface of the wrong radius -- and every one
of those produces plausible output.

Everything is measured from the mesh itself, by fitting the patch's own vertices,
rather than read back from the dictionary that generated it. Reading the
dictionary would only prove the generator is self-consistent.

Checks
------

* required patch names present, with the right geometric types;
* no unexpected duplicate or empty physical patch;
* inflow radius and centre, from the patch vertices;
* inflow entirely on the ``+x`` side, and inside the plume cone;
* cylinder radius, axis direction, axial extent and centre;
* plate bounding box, hence its dimensions and position;
* the 6-inch cylinder-to-plate gap, measured between the two patches;
* the symmetry plane location.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

from plumetools.config import CaseConfig, load_case_config
from plumetools.inflow import read_patch_geometry
from plumetools.markelov1999 import sourceflow as sf
from plumetools.markelov1999.geometry import MarkelovGeometry, from_config
from plumetools.mesh.boundary import read_boundary

#: Relative tolerance on a snapped dimension.
#:
#: snappyHexMesh lands points on an analytic primitive to within its snapping
#: tolerance, and the *measured extent* of a patch is then bounded by the cell
#: size at that surface rather than by machine precision -- a cylinder's measured
#: radius, for instance, is the radius of the snapped polygon through its
#: vertices. 2% comfortably covers that at the shipped refinement levels while
#: still catching a wrong dimension, which would be wrong by tens of percent.
SNAP_TOLERANCE = 0.02


class MeshVerificationError(AssertionError):
    """A generated mesh does not match the geometry the config declares."""


def _patch_points(case_dir: Path, patch: str) -> np.ndarray:
    """The unique vertex coordinates of one patch [m]."""
    _labels, faces, points = read_patch_geometry(case_dir, patch)
    used = sorted({p for f in faces for p in f})
    return points[used]


def verify(case_dir: Path, cfg: CaseConfig | None = None) -> list[str]:
    """Check a case's mesh against its config.

    Args:
        case_dir: an OpenFOAM case directory with ``constant/polyMesh``.
        cfg: the case config; loaded from ``case.yaml`` if omitted.

    Returns:
        Human-readable lines describing what was checked and measured.

    Raises:
        MeshVerificationError: on the first check that fails, with the report so
            far attached so the failure has context.
    """
    case_dir = Path(case_dir)
    cfg = cfg or load_case_config(case_dir)
    geom = from_config(cfg)
    names = cfg.mesh.patch_names

    report: list[str] = [f"case: {case_dir}"]

    def check(ok: bool, message: str) -> None:
        report.append(("  OK   " if ok else "  FAIL ") + message)
        if not ok:
            raise MeshVerificationError("\n".join(report))

    def close(measured: float, expected: float, tol: float = SNAP_TOLERANCE) -> bool:
        return abs(measured - expected) <= tol * max(abs(expected), 1e-12)

    patches = read_boundary(case_dir)
    report.append("patches: " + ", ".join(
        f"{p.name}({p.n_faces}, {p.type})" for p in patches.values()))

    _verify_patches(cfg, patches, check)
    _verify_inflow(case_dir, cfg, geom, check, close, report)
    _verify_cylinder(case_dir, names, geom, check, close, report)
    plate_box = _verify_plate(case_dir, names, geom, check, close, report)
    _verify_gap(case_dir, names, geom, plate_box, check, close, report)
    _verify_symmetry(case_dir, names, geom, check, report)

    return report


def _verify_patches(cfg, patches, check) -> None:
    """Required patches, their types, and no empty physical patch."""
    names = cfg.mesh.patch_names
    required = {
        names["inflow"]: "patch",
        names["cylinder"]: "wall",
        names["plate"]: "wall",
        names["outer"]: "patch",
        names["symmetry"]: "symmetry",
        names["upstream_vacuum"]: "patch",
    }

    missing = sorted(set(required) - set(patches))
    check(not missing, f"all required patches exist (missing: {missing or 'none'})")

    for name, expected_type in required.items():
        actual = patches[name].type
        check(actual == expected_type,
              f"patch {name!r} has geometric type {expected_type!r} (got {actual!r})"
              + (
                  "  -- a body must be 'wall' or hitWallPatch never runs and it "
                  "exerts no force" if expected_type == "wall" else
                  "  -- an open boundary must be 'patch' or particles reflect "
                  "instead of leaving" if expected_type == "patch" else ""))

    empty = sorted(n for n in required if patches[n].n_faces == 0)
    check(not empty,
          f"no required patch is empty (empty: {empty or 'none'})"
          + ("  -- snappyHexMesh creates the patch even when it carved nothing, so "
             "an empty body patch means the body is absent from the mesh"
             if empty else ""))

    unexpected = sorted(set(patches) - set(required))
    if unexpected:
        check(False,
              f"the mesh has unexpected patch(es) {unexpected}; only "
              f"{sorted(required)} are configured, and an extra patch means either "
              f"a stale mesh or a snappyHexMesh surface nothing accounts for")


def _fit_sphere_centre(points: np.ndarray) -> tuple[np.ndarray, float]:
    """Least-squares sphere centre and radius through a set of points.

    For any point on a sphere, ``|p|**2 = 2 p.c + (R**2 - |c|**2)``, which is
    linear in the centre ``c`` and in the scalar ``k = R**2 - |c|**2``. Solving
    that system measures the centre from the mesh instead of assuming it.

    A vertex *centroid* cannot serve here: the meshed surface is only the
    ``x >= 0``, ``y >= 0`` quarter of the sphere, so its centroid sits well away
    from the centre by construction -- roughly ``(0.42R, 0.5R, 0)`` -- and
    comparing it against the origin would fail on a perfectly correct mesh.
    """
    p = np.asarray(points, dtype=np.float64)
    A = np.column_stack([2.0 * p, np.ones(len(p))])
    b = np.einsum("ij,ij->i", p, p)
    solution, *_ = np.linalg.lstsq(A, b, rcond=None)
    centre = solution[:3]
    radius = math.sqrt(max(solution[3] + float(centre @ centre), 0.0))
    return centre, radius


def _verify_inflow(case_dir, cfg, geom, check, close, report) -> None:
    """Inflow radius, centre, hemisphere orientation, and the plume cone."""
    name = cfg.mesh.patch_names["inflow"]
    points = _patch_points(case_dir, name)

    radii = np.linalg.norm(points, axis=1)
    measured = float(radii.mean())
    check(close(measured, geom.inflow_radius_m),
          f"inflow radius is {geom.inflow_radius_m:.6f} m "
          f"(measured mean {measured:.6f}, span {radii.min():.6f} .. {radii.max():.6f})")

    # The centre is the orifice, at the origin -- fitted from the vertices, not
    # assumed. A displaced centre would leave the radii above still tightly
    # clustered about their own mean, so the radius check alone cannot catch it.
    centre, fitted_radius = _fit_sphere_centre(points)
    offset = float(np.linalg.norm(centre))
    check(offset <= SNAP_TOLERANCE * geom.inflow_radius_m,
          f"the inflow surface is centred on the orifice at the origin "
          f"(fitted centre {centre.round(9).tolist()}, offset {offset:.3g} m)")
    check(close(fitted_radius, geom.inflow_radius_m),
          f"the fitted sphere radius is {geom.inflow_radius_m:.6f} m "
          f"(fitted {fitted_radius:.6f})")

    tolerance = SNAP_TOLERANCE * geom.inflow_radius_m
    check(bool(np.all(points[:, 0] >= -tolerance)),
          f"inflow is a hemisphere opening toward +x "
          f"(min vertex x = {points[:, 0].min():.6f})")

    theta_l = sf.limiting_angle(cfg.gas.gamma)
    interior = radii > 0
    theta = sf.off_axis_angle(points[interior])
    check(bool(np.all(theta < theta_l)),
          f"every inflow vertex is inside the {math.degrees(theta_l):.2f} deg plume "
          f"cone (worst {math.degrees(theta.max()):.2f} deg)")

    report.append(f"  ..   inflow: {len(points)} vertices")


def _verify_cylinder(case_dir, names, geom, check, close, report) -> None:
    """Cylinder radius, axis direction, axial extent and centre."""
    points = _patch_points(case_dir, names["cylinder"])

    # Radius in the plane normal to z, about the configured axis.
    radial = np.hypot(points[:, 0] - geom.cylinder_centre_x_m, points[:, 1])
    # End-cap vertices sit anywhere from 0 to R, so the radius is the maximum,
    # not the mean -- the mean would be dragged down by the caps.
    measured_r = float(radial.max())
    check(close(measured_r, geom.cylinder_radius_m),
          f"cylinder radius is {geom.cylinder_radius_m:.6f} m "
          f"(measured {measured_r:.6f})")

    z_min, z_max = float(points[:, 2].min()), float(points[:, 2].max())
    measured_length = z_max - z_min
    check(close(measured_length, geom.cylinder_length_m),
          f"cylinder length is {geom.cylinder_length_m:.6f} m "
          f"(measured {measured_length:.6f}, z {z_min:.6f} .. {z_max:.6f})")

    measured_centre_z = 0.5 * (z_min + z_max)
    check(abs(measured_centre_z - geom.cylinder_centre_z_m)
          <= SNAP_TOLERANCE * geom.cylinder_length_m,
          f"cylinder is centred at z = {geom.cylinder_centre_z_m:.6f} m "
          f"(measured {measured_centre_z:.6f})")

    x_min, x_max = float(points[:, 0].min()), float(points[:, 0].max())
    measured_centre_x = 0.5 * (x_min + x_max)
    check(close(measured_centre_x, geom.cylinder_centre_x_m),
          f"cylinder is centred at x = {geom.cylinder_centre_x_m:.6f} m "
          f"(measured {measured_centre_x:.6f})")

    # Axis direction: a z-axis cylinder has its x and y extents bounded by the
    # diameter while z spans the length. If the axis were along x or y the
    # extents would swap, which no radius check alone would catch.
    check(measured_length > 2.0 * measured_r,
          f"the cylinder axis is z: its z extent ({measured_length:.6f} m) exceeds "
          f"its diameter ({2 * measured_r:.6f} m)")
    check((x_max - x_min) <= 2.0 * geom.cylinder_radius_m * (1 + SNAP_TOLERANCE),
          f"the cylinder does not extend along x beyond its diameter "
          f"(x extent {x_max - x_min:.6f} m)")

    report.append(f"  ..   cylinder: {len(points)} vertices, "
                  f"upstream x = {x_min:.6f}, downstream x = {x_max:.6f}")


def _verify_plate(case_dir, names, geom, check, close, report):
    """Plate dimensions and position, from its bounding box."""
    points = _patch_points(case_dir, names["plate"])

    lo = points.min(axis=0)
    hi = points.max(axis=0)

    thickness = float(hi[0] - lo[0])
    check(close(thickness, geom.plate_thickness_m),
          f"plate thickness is {geom.plate_thickness_m:.6f} m "
          f"(measured {thickness:.6f}; an ASSUMPTION, not a paper value)")

    check(close(float(lo[0]), geom.plate_upstream_x_m),
          f"plate upstream face is at x = {geom.plate_upstream_x_m:.6f} m "
          f"(measured {lo[0]:.6f}) -- derived as cylinder centre + radius + gap")

    check(close(float(hi[0]), geom.plate_downstream_x_m),
          f"plate downstream face is at x = {geom.plate_downstream_x_m:.6f} m "
          f"(measured {hi[0]:.6f})")

    height = float(hi[2] - lo[2])
    check(close(height, geom.plate_height_z_m),
          f"plate height is {geom.plate_height_z_m:.6f} m (measured {height:.6f})")

    # Only the modelled half is present, so the y extent is half the plate width.
    width_half = float(hi[1] - lo[1])
    check(close(width_half, geom.plate_y_max_m),
          f"plate half-width is {geom.plate_y_max_m:.6f} m (measured "
          f"{width_half:.6f}); the full {geom.plate_width_y_m:.6f} m is recovered "
          f"by the y = 0 symmetry plane")

    check(abs(0.5 * (lo[2] + hi[2])) <= SNAP_TOLERANCE * geom.plate_height_z_m,
          f"plate centreline is on the plume axis "
          f"(measured z centre {0.5 * (lo[2] + hi[2]):.6f})")

    report.append(f"  ..   plate: {len(points)} vertices, box "
                  f"x [{lo[0]:.6f}, {hi[0]:.6f}]  y [{lo[1]:.6f}, {hi[1]:.6f}]  "
                  f"z [{lo[2]:.6f}, {hi[2]:.6f}]")
    return lo, hi


def _verify_gap(case_dir, names, geom, plate_box, check, close, report) -> None:
    """The 6-inch gap, measured between the two meshed patches.

    Measured, not recomputed: the cylinder's downstream extent comes from the
    cylinder patch's own vertices and the plate's upstream face from the plate
    patch's, so this checks the two bodies really are that far apart in the mesh
    rather than that the arithmetic in the config is self-consistent.
    """
    cylinder = _patch_points(case_dir, names["cylinder"])
    lo, _hi = plate_box

    cylinder_base_x = float(cylinder[:, 0].max())
    plate_face_x = float(lo[0])
    measured = plate_face_x - cylinder_base_x

    check(close(measured, geom.gap_m),
          f"cylinder-to-plate gap is {geom.gap_m:.6f} m "
          f"({geom.gap_m / 0.0254:.3f} in), measured {measured:.6f} m "
          f"({measured / 0.0254:.3f} in) between the cylinder base at "
          f"x = {cylinder_base_x:.6f} and the plate face at x = {plate_face_x:.6f}")


def _verify_symmetry(case_dir, names, geom, check, report) -> None:
    """The symmetry patch is a plane at the configured y, and the only one."""
    points = _patch_points(case_dir, names["symmetry"])

    y_span = float(points[:, 1].max() - points[:, 1].min())
    check(y_span < 1e-9,
          f"the symmetry patch is a plane of constant y (span {y_span:.3g} m)")

    y = float(points[:, 1].mean())
    check(abs(y - geom.symmetry_plane_y_m) < 1e-9,
          f"the symmetry plane is at y = {geom.symmetry_plane_y_m:.6f} "
          f"(measured {y:.9f})")

    # The full z extent must survive: halving z would delete the cylinder end
    # caps and the plate edges, which is the 3D wake this case exists to study.
    z_span = float(points[:, 2].max() - points[:, 2].min())
    check(points[:, 2].min() < 0.0 < points[:, 2].max(),
          f"the domain keeps the full z extent (symmetry patch spans "
          f"{points[:, 2].min():.4f} .. {points[:, 2].max():.4f} m, {z_span:.4f} m), "
          f"so the cylinder end caps and plate edges are inside it")

    report.append(f"  ..   symmetry: {len(points)} vertices at y = {y:.9f}")


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
