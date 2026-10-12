r"""Verify a built mesh against the geometry ``case.yaml`` declares.

Run after the mesh pipeline::

    python -m plumetools.nozzleplume verify <case>

``checkMesh`` answers "is this a valid mesh?". This answers a question
``checkMesh`` cannot: **is this the mesh that was planned?** A mesh can be
flawless and still have a nozzle of the wrong area, a plate at the wrong
stand-off, or a boundary of the wrong type -- and every one of those produces a
plausible result.

Everything is measured from the mesh itself, from the patches' own faces, rather
than read back from the dictionaries that generated it. Reading the dictionaries
would only prove the generator is self-consistent.

Checks
------

* required patches present, with the right geometric types, and no others;
* the cell count, against the estimate and against ``mesh.max_cells``;
* the nozzle: face count and area against the prediction, every face on
  ``x = 0`` and inside ``R0``, every normal along ``-x``;
* the plate: its bounding box against the configured stand-off, thickness and
  extents, and the area of its upstream face;
* each symmetry plane: flat, and where it should be.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from plumetools.mesh.boundary import read_boundary
from plumetools.nozzleplume import checks
from plumetools.nozzleplume.meshread import read_cell_count, read_patches

#: Relative tolerance on a measured area.
#:
#: The nozzle patch is made of whole lattice faces and the plate is snapped to
#: its box with feature snapping, so both are exact to round-off on a correct
#: mesh. 0.5% is loose enough for that never to trip and tight enough to catch a
#: chamfered plate edge -- which costs 1.4% of the upstream face when feature
#: snapping is off.
AREA_TOLERANCE = 0.005

#: Tolerance on a measured plate coordinate, as a fraction of the plate cell.
POSITION_TOLERANCE = 0.05


class MeshVerificationError(AssertionError):
    """A built mesh does not match the geometry the config declares."""


def verify(case_dir: Path, cfg, geom, plan, nozzle) -> list[str]:
    """Check a case's mesh against its config.

    Args:
        case_dir: an OpenFOAM case directory with ``constant/polyMesh``.
        cfg: the case config.
        geom: the :class:`~plumetools.nozzleplume.geometry.PlumeGeometry`.
        plan: the :class:`~plumetools.nozzleplume.mesh.MeshPlan`.
        nozzle: the :class:`~plumetools.nozzleplume.mesh.NozzlePrediction`.

    Returns:
        Human-readable lines describing what was checked and measured.

    Raises:
        MeshVerificationError: on the first check that fails, with the report so
            far attached so the failure has context.
    """
    case_dir = Path(case_dir)
    names = cfg.mesh.patch_names
    report: list[str] = ["Mesh verification (measured from constant/polyMesh)"]

    def check(ok: bool, message: str) -> None:
        report.append(("  OK   " if ok else "  FAIL ") + message)
        if not ok:
            raise MeshVerificationError("\n".join(report))

    boundary = read_boundary(case_dir)
    expected = checks.expected_patch_types(cfg, geom)

    # --- patches --------------------------------------------------------------
    missing = sorted(set(expected) - set(boundary))
    check(not missing,
          f"all required patches exist (missing: {missing or 'none'})")
    for name, kind in expected.items():
        check(boundary[name].type == kind,
              f"patch {name!r} is geometric type {kind!r} "
              f"(got {boundary[name].type!r})")
    empty = sorted(name for name in expected if boundary[name].n_faces == 0)
    check(not empty,
          f"no required patch is empty (empty: {empty or 'none'})"
          + ("  -- snappyHexMesh creates a body's patch even when it carved "
             "nothing, so an empty one means the body is absent" if empty else ""))
    unexpected = sorted(set(boundary) - set(expected))
    check(not unexpected,
          f"the mesh has no patch case.yaml does not describe "
          f"(unexpected: {unexpected or 'none'})"
          + ("  -- a stale mesh from an earlier configuration; re-run "
             "./Allmesh --force" if unexpected else ""))

    # --- cells ----------------------------------------------------------------
    n_cells = read_cell_count(case_dir)
    if n_cells is not None:
        check(n_cells <= plan.max_cells,
              f"{n_cells:,} cells, within mesh.max_cells = {plan.max_cells:,}"
              + ("" if n_cells <= plan.max_cells else
                 "  -- snappyHexMesh stops refining at twice that, so a mesh "
                 "this far over may also be quietly under-refined"))
        report.append(
            f"  ..   the estimate was {plan.estimated_cells:,} "
            f"({100.0 * (plan.estimated_cells / n_cells - 1.0):+.1f}%)")

    wanted = [names["nozzle"]]
    if geom.plate is not None:
        wanted.append(names["plate"])
    if geom.symmetry_planes:
        wanted.append(names["symmetry"])
    patches = read_patches(case_dir, wanted)

    _verify_nozzle(patches[names["nozzle"]], geom, plan, nozzle, check, report)
    if geom.plate is not None:
        _verify_plate(patches[names["plate"]], geom, plan, check, report)
    if geom.symmetry_planes:
        _verify_symmetry(patches[names["symmetry"]], geom, plan, check, report)
    return report


def _verify_nozzle(patch, geom, plan, nozzle, check, report) -> None:
    """Face count, area, position and orientation of the nozzle patch."""
    check(patch.n_faces == nozzle.n_faces,
          f"the nozzle has the {nozzle.n_faces} faces the plan predicted "
          f"(meshed {patch.n_faces})"
          + ("" if patch.n_faces == nozzle.n_faces else
             "  -- the exit is not at one refinement level, or the lattice is "
             "not where the plan put it"))

    area = patch.area_m2
    check(abs(area - nozzle.area_m2) <= AREA_TOLERANCE * nozzle.area_m2,
          f"the nozzle area is the predicted {nozzle.area_m2:.6e} m^2 "
          f"(measured {area:.6e})")

    tolerance = 1.0e-6 * plan.exit_cell_size_m
    off_plane = float(np.abs(patch.vertices[:, 0]).max())
    check(off_plane <= tolerance,
          f"every nozzle vertex is on the exit plane x = 0 "
          f"(worst |x| = {off_plane:.3g} m)")

    radii = np.hypot(patch.centres[:, 1], patch.centres[:, 2])
    check(float(radii.max()) <= geom.radius_m * (1.0 + 1.0e-9),
          f"every nozzle face centre is inside R0 = {geom.radius_m:.6g} m "
          f"(furthest {radii.max():.6g} m)")

    check(bool(np.all(patch.normals[:, 0] < -0.999)),
          "every nozzle face points out of the domain along -x, so the "
          "injected flux is along +x")

    error = (area - geom.exact_nozzle_area_m2) / geom.exact_nozzle_area_m2
    report.append(
        f"  ..   nozzle: {patch.n_faces} faces, {area:.6e} m^2, "
        f"{100.0 * error:+.2f}% from the exact {geom.exact_nozzle_area_m2:.6e} m^2")


def _verify_plate(patch, geom, plan, check, report) -> None:
    """The plate's box, from its vertices, and the area of its upstream face."""
    x0, x1, y0, y1, z0, z1 = geom.modelled_plate_bounds
    lo = patch.vertices.min(axis=0)
    hi = patch.vertices.max(axis=0)
    tolerance = POSITION_TOLERANCE * plan.plate_cell_size_m

    for label, measured, wanted in (
            ("upstream face x", lo[0], x0), ("lee face x", hi[0], x1),
            ("y min", lo[1], y0), ("y max", hi[1], y1),
            ("z min", lo[2], z0), ("z max", hi[2], z1)):
        check(abs(float(measured) - wanted) <= tolerance,
              f"plate {label} is {wanted:.6g} m (measured {float(measured):.6g})")

    # The upstream face: outward normals along +x, i.e. into the plate.
    normals = patch.normals
    front = normals[:, 0] > 0.999
    front_area = float(patch.areas[front].sum())
    wanted_area = (y1 - y0) * (z1 - z0)
    check(abs(front_area - wanted_area) <= AREA_TOLERANCE * wanted_area,
          f"the plate's upstream face is {wanted_area:.6g} m^2 "
          f"(measured {front_area:.6g} over {int(front.sum())} faces)"
          + ("" if abs(front_area - wanted_area) <= AREA_TOLERANCE * wanted_area
             else "  -- a shortfall with oblique faces present means the "
                  "plate's edges were chamfered, not snapped"))

    front_x = patch.centres[front, 0]
    check(float(np.abs(front_x - x0).max()) <= tolerance,
          f"the upstream face is flat at x = {x0:.6g} m "
          f"(face centres {front_x.min():.6g} .. {front_x.max():.6g})")

    oblique = int(np.count_nonzero(np.abs(normals).max(axis=1) < 0.999))
    report.append(
        f"  ..   plate: {patch.n_faces} faces, {patch.area_m2:.6g} m^2, "
        f"{oblique} not axis-aligned; box x [{lo[0]:.6g}, {hi[0]:.6g}] "
        f"y [{lo[1]:.6g}, {hi[1]:.6g}] z [{lo[2]:.6g}, {hi[2]:.6g}]")


def _verify_symmetry(patch, geom, plan, check, report) -> None:
    """Each symmetry plane is flat and at zero.

    To a thousandth of the smallest cell, not to round-off. Where the plate
    crosses a symmetry plane snappyHexMesh moves the vertices it snaps, and does
    not put them back exactly: measured, face centres up to 0.6 micrometres off
    the plane on the R-4D quarter model. That is not a displaced plane -- a
    plane in the wrong place is out by a cell -- and the patch is type
    ``symmetry``, which unlike ``symmetryPlane`` does not require flatness.
    """
    tolerance = 1.0e-3 * plan.min_cell_size_m
    for axis in geom.symmetry_planes:
        index = "xyz".index(axis)
        # A face belongs to this plane if its outward normal is along -axis.
        mine = patch.normals[:, index] < -0.999
        check(bool(mine.any()),
              f"the symmetry patch has faces on the plane {axis} = 0")
        worst = float(np.abs(patch.centres[mine, index]).max())
        check(worst <= tolerance,
              f"the symmetry plane {axis} = 0 is flat and at zero "
              f"(worst face centre {worst:.3g} m)")
    report.append(f"  ..   symmetry: {patch.n_faces} faces on "
                  + " and ".join(f"{axis} = 0" for axis in geom.symmetry_planes))


def main(case_dir: Path, cfg, geom, plan, nozzle) -> int:
    """Print the verification report; return a process exit code."""
    try:
        for line in verify(case_dir, cfg, geom, plan, nozzle):
            print(line)
    except MeshVerificationError as exc:
        print(str(exc), file=sys.stderr)
        print("\nmesh verification FAILED", file=sys.stderr)
        return 1
    print("\nmesh verification passed")
    return 0
