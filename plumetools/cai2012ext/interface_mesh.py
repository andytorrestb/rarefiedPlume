r"""The truncated free-plume mesh: the Cai box with the dense core carved out.

Pipeline (``baseCaseInterface/Allmesh``)::

    blockMesh                 the parent family's graded box, re-planned so its
                              uniform core encloses the interface sphere
    snappyHexMesh -overwrite  remove r < R_interface around the exit centre and
                              snap the cut to the analytic sphere; the snapped
                              surface becomes the 'nozzle' patch
    checkMesh                 quality, and one connected region
    interface_case verify-mesh   the patch is the surface the model was chosen for

Composition only. The lattice is :func:`plumetools.cai2012.mesh.plan` -- the
same segments, the same grading rule, the same blockMeshDict writer -- with the
core extents and the core cell re-derived from the interface instead of the
exit; the carve is :mod:`plumetools.foamio.primitives`' ``SearchableSphere`` and
``snappyHexMeshDict``, as ``plumetools.markelov1999.mesh`` carves its inflow
hemisphere out of ``x = 0``. No cylinder, plate or STL.

The patch keeps the name ``nozzle``
-----------------------------------
so ``constant/dsmcProperties`` (``plumeFieldInflow``'s ``patches ( nozzle )``)
and every field writer work unchanged. Its physical role is different: it is
the downstream continuum interface, not the exit disk, and every report says
so. ``x = 0`` outside the hemisphere stays ``upstreamVacuum``, an open
boundary; there is no wall at the angular cutoff and no symmetry plane.

The core cell
-------------
``target_cell_over_mfp`` times the closure's mean free path at the interface on
the axis -- the densest retained point -- capped by ``max_core_cell_over_D``
like the parent. The budget counts the cells that SURVIVE the carve (lattice
centres outside the sphere), not the lattice, and coarsens by the parent's rule
if they exceed ``mesh.max_cells``; the plan says so when it does.
"""

from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from plumetools.cai2012 import mesh as cai_mesh
from plumetools.foamio.primitives import (
    SearchableSphere,
    check_surfaces_are_resolvable,
    fmt,
    render_mesh_quality_dict,
    render_snappy_dict,
    write_lines,
)

GENERATOR = "plumetools.cai2012ext.interface_mesh"

#: Fine cells kept beyond the sphere before the expanding segments start, as a
#: fraction of ``interface.domain_clearance_over_D`` (and at most 0.25 D).
CORE_MARGIN_FRACTION = 0.5
CORE_MARGIN_MAX_OVER_D = 0.25

#: A budget the lattice may not exceed even before the carve, so that a
#: blockMesh of absurd size is refused rather than attempted.
LATTICE_HEADROOM = 4.0


@dataclass(frozen=True)
class InterfaceMeshPlan:
    """The truncated mesh, before a line of it is written.

    Attributes:
        lattice: the parent family's :class:`~plumetools.cai2012.mesh.MeshPlan`
            for the blockMesh box (cells inside the sphere included).
        radius_m: ``R_interface``.
        centre_m: the sphere centre, the source point.
        mean_free_path_m: the closure's ``lambda`` at ``(R, 0)``.
        requested_cell_size_m, coarsened: as the parent plan's.
        retained_cells: lattice cells whose centre lies outside the sphere --
            the mesh snappyHexMesh leaves, to within the snapped layer.
        max_cells: the budget ``retained_cells`` is held to.
        core_x_m, core_half_m: the uniform core's extent.
    """

    lattice: cai_mesh.MeshPlan
    radius_m: float
    centre_m: tuple
    mean_free_path_m: float
    requested_cell_size_m: float
    coarsened: bool
    retained_cells: int
    max_cells: int
    core_x_m: float
    core_half_m: float

    @property
    def core_cell_size_m(self) -> float:
        return float(self.lattice.core_cell_size_m)

    @property
    def min_cell_size_m(self) -> float:
        return float(self.lattice.min_cell_size_m)

    @property
    def core_cell_sizes_m(self) -> tuple:
        """``(dx, dy, dz)`` of a uniform-core cell: the segments holding the axis.

        Not :attr:`MeshPlan.core_cell_sizes_m`, which takes the first uniform
        lateral segment -- an outer segment, when one comes out uniform.
        """
        def holding_zero(segments):
            return next(s for s in segments if s.lo <= 0.0 < s.hi).first_cell_m
        return (self.lattice.x_segments[0].first_cell_m,
                holding_zero(self.lattice.y_segments),
                holding_zero(self.lattice.z_segments))

    @property
    def core_cell_volume_m3(self) -> float:
        dx, dy, dz = self.core_cell_sizes_m
        return float(dx * dy * dz)

    @property
    def cell_over_mfp(self) -> float:
        return self.core_cell_size_m / self.mean_free_path_m

    @property
    def predicted_inlet_faces(self) -> int:
        """``2 pi R^2`` over a core cell face: the hemisphere in core cells."""
        dx, dy, dz = self.core_cell_sizes_m
        return int(round(2.0 * math.pi * self.radius_m ** 2 / (dy * dz)))

    def as_dict(self) -> dict:
        out = self.lattice.as_dict()
        out.update({
            "lattice_cells": int(self.lattice.n_cells),
            "n_cells": int(self.retained_cells),
            "retained_cells_basis": "lattice cell centres outside the sphere",
            "interface_radius_m": float(self.radius_m),
            "interface_centre_m": list(self.centre_m),
            "interface_mean_free_path_m": float(self.mean_free_path_m),
            "cell_over_mean_free_path": float(self.cell_over_mfp),
            "requested_cell_size_m": float(self.requested_cell_size_m),
            "coarsened_for_budget": bool(self.coarsened),
            "max_cells_budget": int(self.max_cells),
            "core_x_m": float(self.core_x_m),
            "core_half_m": float(self.core_half_m),
            "predicted_inlet_faces": self.predicted_inlet_faces,
            "inlet_patch": "nozzle (role: downstream continuum interface)",
        })
        return out

    def describe(self) -> list[str]:
        lines = [
            "Interface mesh plan (graded Cartesian box, sphere carved by "
            "snappyHexMesh)",
            f"  interface                sphere R = {fmt(self.radius_m)} m at "
            f"{tuple(self.centre_m)}, patch 'nozzle'",
            f"  lambda at (R, 0)         {self.mean_free_path_m:.4e} m",
            f"  core cell                {fmt(self.core_cell_size_m)} m, "
            f"cell/lambda {self.cell_over_mfp:.3f}",
            f"  lattice / retained       {self.lattice.n_cells:,} / "
            f"{self.retained_cells:,} cells (budget {self.max_cells:,})",
            f"  inlet faces (predicted)  ~{self.predicted_inlet_faces:,}",
        ]
        if self.coarsened:
            lines.append(
                f"  COARSENED: lambda asked for {fmt(self.requested_cell_size_m)} "
                f"m, which would not fit in {self.max_cells:,} retained cells.")
        return lines


def _retained(lattice: cai_mesh.MeshPlan, radius_m: float) -> int:
    """Lattice cells whose centre lies outside the sphere at the origin."""
    xs = np.concatenate([s.cell_centres() for s in lattice.x_segments])
    ys = np.concatenate([s.cell_centres() for s in lattice.y_segments])
    zs = np.sort(np.concatenate([s.cell_centres() for s in lattice.z_segments]))
    r2 = radius_m * radius_m
    removed = 0
    for x in xs:
        disk = r2 - x * x
        if disk <= 0.0:
            continue
        room = disk - ys * ys
        half = np.sqrt(np.clip(room, 0.0, None))
        inside = (np.searchsorted(zs, half, side="left")
                  - np.searchsorted(zs, -half, side="right"))
        removed += int(np.sum(np.where(room > 0.0, inside, 0)))
    return int(len(xs) * len(ys) * len(zs) - removed)


def plan(cfg, geom, settings, *, radius_m: float, mean_free_path_m: float
         ) -> InterfaceMeshPlan:
    """Plan the truncated mesh for a sphere of ``radius_m`` [m].

    Raises:
        ValueError: if the sphere and its fine margin do not fit in the box, or
            the retained-cell budget cannot be met.
    """
    d = geom.diameter_m
    margin = min(CORE_MARGIN_MAX_OVER_D * d,
                 CORE_MARGIN_FRACTION * settings.domain_clearance_over_D * d)
    core_x = radius_m + margin
    core_half = radius_m + margin
    lateral = min(geom.y_max_m, -geom.y_min_m, geom.z_max_m, -geom.z_min_m)
    if not (core_x < geom.x_max_m and core_half < lateral):
        raise ValueError(
            f"the interface sphere (R = {radius_m / d:.3f} D) plus its "
            f"{margin / d:.3f} D fine margin does not fit inside the box "
            f"(x_max {geom.x_max_m / d:g} D, half-width {lateral / d:g} D); "
            f"enlarge the domain through geometry.* explicitly")

    if cfg.mesh.core_cell_size_m is not None:
        # Pinned, as the parent allows: a development mesh says so itself.
        requested = float(cfg.mesh.core_cell_size_m)
    else:
        requested = min(float(cfg.mesh.target_cell_over_mfp) * mean_free_path_m,
                        float(cfg.mesh.max_core_cell_over_D) * d)
    if cfg.mesh.min_core_cell_size_m is not None:
        requested = max(requested, float(cfg.mesh.min_core_cell_size_m))
    budget = int(cfg.mesh.max_cells)

    def lattice_for(cell: float) -> cai_mesh.MeshPlan:
        # The parent's outer cell counts assume a thick far field. Beyond a
        # sphere near the box edge the shell is thin, and 34 or 40 cells there
        # would be SMALLER than the core cell -- shrinking the time step and
        # inverting the grading. Never more outer cells than core cells fit.
        room_x, room_lateral = geom.x_max_m - core_x, lateral - core_half
        if min(room_x, room_lateral) < cell:
            raise ValueError(
                f"the {min(room_x, room_lateral) / d:.3f} D left between the fine "
                f"core and the box is less than one {cell:.4g} m core cell; "
                f"enlarge the domain through geometry.* explicitly")
        mesh_cfg = dataclasses.replace(
            cfg.mesh, core_cell_size_m=float(cell), min_core_cell_size_m=None,
            core_x_over_D=core_x / d, core_half_over_D=core_half / d,
            outer_x_cells=max(1, min(int(cfg.mesh.outer_x_cells), int(room_x // cell))),
            outer_lateral_cells=max(1, min(int(cfg.mesh.outer_lateral_cells),
                                           int(room_lateral // cell))),
            max_cells=2 ** 62)
        return cai_mesh.plan(dataclasses.replace(cfg, mesh=mesh_cfg), geom,
                             _MeanFreePath(mean_free_path_m))

    cell = requested
    coarsened = False
    # A cheap guard first: the lattice count from the segment counts alone, so
    # that a hopeless request is refused before any cell centre is enumerated.
    estimate = _lattice_cells_estimate(geom, core_x, core_half, cell, cfg)
    if estimate > LATTICE_HEADROOM * budget:
        coarsened = True
        cell = requested * (estimate / (LATTICE_HEADROOM * budget)) ** (1.0 / 3.0)
    lattice = lattice_for(cell)
    retained = _retained(lattice, radius_m)
    while retained > budget:
        coarsened = True
        cell *= max(1.02, (retained / budget) ** (1.0 / 3.0))
        if cell > radius_m:
            raise ValueError(
                f"meeting mesh.max_cells = {budget:,} retained cells needs a "
                f"core cell coarser than the interface radius; raise "
                f"mesh.max_cells or shrink the domain")
        lattice = lattice_for(cell)
        retained = _retained(lattice, radius_m)

    return InterfaceMeshPlan(
        lattice=lattice, radius_m=float(radius_m), centre_m=(0.0, 0.0, 0.0),
        mean_free_path_m=float(mean_free_path_m),
        requested_cell_size_m=float(requested), coarsened=coarsened,
        retained_cells=int(retained), max_cells=budget,
        core_x_m=float(core_x), core_half_m=float(core_half))


@dataclass(frozen=True)
class _MeanFreePath:
    """The one attribute :func:`plumetools.cai2012.mesh.plan` reads of a state."""

    mean_free_path_m: float


def _lattice_cells_estimate(geom, core_x, core_half, cell, cfg) -> float:
    """Core cells plus the fixed outer counts, without building the lattice."""
    nx = math.ceil(core_x / cell) + int(cfg.mesh.outer_x_cells)
    ny = math.ceil(2.0 * core_half / cell) + 2 * int(cfg.mesh.outer_lateral_cells)
    return float(nx) * ny * ny


def resolution_estimate(source, mesh_plan: InterfaceMeshPlan, geom, theta_c: float,
                        n_samples: int = 400) -> dict:
    """Largest planned ``cell/lambda`` along the axis and the core edge.

    Pre-run only: the closure's ``lambda`` is a continuum estimate, and beyond
    breakdown the gas is not that state. ``Alldiagnose.py`` measures the real
    ratio cell by cell after the run.
    """
    axes = [np.unique(np.concatenate([s.cell_edges() for s in segments]))
            for segments in (mesh_plan.lattice.x_segments,
                             mesh_plan.lattice.y_segments,
                             mesh_plan.lattice.z_segments)]

    def width(axis_edges, value):
        i = np.clip(np.searchsorted(axis_edges, value) - 1, 0, len(axis_edges) - 2)
        return axis_edges[i + 1] - axis_edges[i]

    out = {}
    for name, theta in (("axis", 0.0), ("core_edge", float(theta_c))):
        direction = np.array([math.cos(theta), math.sin(theta), 0.0])
        limits = [geom.x_max_m / max(direction[0], 1e-300),
                  geom.y_max_m / max(direction[1], 1e-300)]
        r_end = min(limits) * (1.0 - 1e-9)
        r = np.linspace(mesh_plan.radius_m * (1.0 + 1e-9), r_end, n_samples)
        points = r[:, None] * direction[None, :]
        cell = np.maximum.reduce([width(axes[k], points[:, k]) for k in range(3)])
        st = source.state(r, np.full(r.shape, theta))
        ratio = cell / st["lambda"]
        worst = int(np.nanargmax(ratio))
        out[name] = {"theta_deg": math.degrees(theta),
                     "max_cell_over_lambda": float(ratio[worst]),
                     "at_r_over_D": float(r[worst] / geom.diameter_m),
                     "cell_over_lambda_at_interface": float(ratio[0])}
    out["basis"] = ("planned cell edge (largest of dx, dy, dz) over the "
                    "closure's continuum mean free path; an estimate, measured "
                    "after the run by Alldiagnose.py")
    return out


def location_in_mesh(mesh_plan: InterfaceMeshPlan) -> tuple:
    """A cell centre in the fluid: the last x cell, just off both centre planes."""
    xs = np.concatenate([s.cell_centres() for s in mesh_plan.lattice.x_segments])
    ys = np.concatenate([s.cell_centres() for s in mesh_plan.lattice.y_segments])
    zs = np.concatenate([s.cell_centres() for s in mesh_plan.lattice.z_segments])
    point = (float(xs[-1]), float(ys[np.searchsorted(ys, 0.0)]),
             float(zs[np.searchsorted(zs, 0.0)]))
    if math.sqrt(sum(v * v for v in point)) <= mesh_plan.radius_m:
        raise ValueError(f"locationInMesh {point} is inside the interface sphere")
    return point


def sphere(cfg, mesh_plan: InterfaceMeshPlan) -> SearchableSphere:
    """The carved surface, named after the inflow patch, geometric type patch."""
    return SearchableSphere(
        name=cfg.mesh.patch_names["nozzle"], level=0, patch_type="patch",
        centre=tuple(mesh_plan.centre_m), radius=float(mesh_plan.radius_m))


def render_snappy_hex_mesh_dict(cfg, mesh_plan: InterfaceMeshPlan) -> str:
    surface = sphere(cfg, mesh_plan)
    check_surfaces_are_resolvable([surface], mesh_plan.core_cell_size_m)
    return render_snappy_dict(
        [surface],
        location_in_mesh=location_in_mesh(mesh_plan),
        n_cells_between_levels=1,
        generator=GENERATOR,
        note=(f"Cai 2012 interface mode: remove r < {fmt(mesh_plan.radius_m)} m "
              f"around the exit centre; the snapped sphere is the 'nozzle' "
              f"patch, injecting the continuum state, NOT the exit disk."),
        preamble=[
            "// Level 0: no refinement. The graded blockMesh core already has the",
            "// interface cell size; snappyHexMesh only removes the cells inside",
            "// the sphere and snaps the cut onto it. 'patch', never 'wall':",
            "// DSMC particles crossing back through the interface are deleted.",
        ])


def render_block_mesh_dict(cfg, geom, mesh_plan: InterfaceMeshPlan) -> str:
    """The parent family's blockMeshDict for the re-planned lattice, annotated."""
    mesh_cfg = dataclasses.replace(
        cfg.mesh, core_x_over_D=mesh_plan.core_x_m / geom.diameter_m,
        core_half_over_D=mesh_plan.core_half_m / geom.diameter_m)
    text = cai_mesh.render_block_mesh_dict(
        dataclasses.replace(cfg, mesh=mesh_cfg), geom, mesh_plan.lattice)
    note = [
        "// INTERFACE MODE (plumetools.cai2012ext.interface_mesh). The nozzle",
        "// disk is NOT carved by topoSet/createPatch in this case:",
        "// snappyHexMesh removes r < R_interface and the snapped sphere becomes",
        "// the 'nozzle' patch. The comments below are the parent family's.",
        "",
    ]
    lines = text.splitlines()
    at = next((i for i, line in enumerate(lines) if line.startswith("scale")),
              0)
    return "\n".join(lines[:at] + note + lines[at:]) + "\n"


def write_mesh_setup(case_dir: Path, cfg, geom, mesh_plan: InterfaceMeshPlan
                     ) -> list[Path]:
    """``blockMeshDict``, ``snappyHexMeshDict`` and ``meshQualityDict``.

    No ``topoSetDict``/``createPatchDict``: a stale pair from a direct-exit
    mesh would carve a disk that no longer exists, so they are removed.
    """
    system = Path(case_dir) / "system"
    for stale in ("topoSetDict", "createPatchDict"):
        (system / stale).unlink(missing_ok=True)
    return [
        write_lines(system / "blockMeshDict",
                    render_block_mesh_dict(cfg, geom, mesh_plan).splitlines()),
        write_lines(system / "snappyHexMeshDict",
                    render_snappy_hex_mesh_dict(cfg, mesh_plan).splitlines()),
        write_lines(system / "meshQualityDict",
                    render_mesh_quality_dict(GENERATOR).splitlines()),
    ]


# --------------------------------------------------------------------------- #
# the meshed patch
# --------------------------------------------------------------------------- #

def patch_report(*, centres, outward_normals, areas, points, radius_m: float,
                 diameter_m: float, theta_core: float, theta_limit: float
                 ) -> dict:
    """Measure the meshed inlet patch: radius, centre, normals, area, angles.

    Everything from the patch's own geometry, not from the dictionary that
    asked for it.
    """
    c = np.asarray(centres, dtype=np.float64)
    n = np.asarray(outward_normals, dtype=np.float64)
    a = np.asarray(areas, dtype=np.float64)
    p = np.asarray(points, dtype=np.float64)
    vertex_r = np.linalg.norm(p, axis=1)
    centroid_r = np.linalg.norm(c, axis=1)
    r_hat = c / centroid_r[:, None]
    facing = -np.einsum("ij,ij->i", n, r_hat)
    theta = np.arccos(np.clip(r_hat[:, 0], -1.0, 1.0))
    # Least-squares sphere through the vertices: |p|^2 = 2 p.c + (R^2 - |c|^2).
    lhs = np.hstack([2.0 * p, np.ones((len(p), 1))])
    solution, *_ = np.linalg.lstsq(lhs, np.einsum("ij,ij->i", p, p), rcond=None)
    centre = solution[:3]
    fitted = math.sqrt(max(solution[3] + centre @ centre, 0.0))
    hemisphere = 2.0 * math.pi * radius_m ** 2
    return {
        "faces": int(len(a)),
        "area_m2": float(a.sum()),
        "hemisphere_area_m2": hemisphere,
        "area_error": float(a.sum() / hemisphere - 1.0),
        "vertex_radius_over_R": [float(vertex_r.min() / radius_m),
                                 float(vertex_r.max() / radius_m)],
        "centroid_radius_over_R": [float(centroid_r.min() / radius_m),
                                   float(centroid_r.max() / radius_m)],
        "fitted_radius_over_R": float(fitted / radius_m),
        "fitted_centre_over_D": [float(v / diameter_m) for v in centre],
        "min_x_over_D": float(c[:, 0].min() / diameter_m),
        "theta_max_deg": float(np.degrees(theta.max())),
        "faces_in_core": int(np.count_nonzero(theta <= theta_core)),
        "faces_beyond_limit_angle": int(np.count_nonzero(theta >= theta_limit)),
        "min_inward_normal_dot_radial": float(facing.min()),
        "faces_facing_the_source": int(np.count_nonzero(facing <= 0.0)),
    }


def patch_problems(report: dict, *, snap_tolerance: float, min_core_faces: int,
                   area_tolerance: float) -> list[str]:
    """Why the meshed patch is not the surface the model was chosen for."""
    problems = []
    lo, hi = report["vertex_radius_over_R"]
    if lo < 1.0 - snap_tolerance or hi > 1.0 + snap_tolerance:
        problems.append(f"patch vertices lie at {lo:.4f}..{hi:.4f} R, outside "
                        f"+/-{snap_tolerance:g} of the sphere")
    if abs(report["fitted_radius_over_R"] - 1.0) > snap_tolerance:
        problems.append(f"fitted radius {report['fitted_radius_over_R']:.4f} R")
    if max(abs(v) for v in report["fitted_centre_over_D"]) > snap_tolerance:
        problems.append(f"fitted centre {report['fitted_centre_over_D']} D is not "
                        f"the source point")
    if report["min_x_over_D"] < -1.0e-9:
        problems.append("patch faces upstream of x = 0")
    if report["faces_facing_the_source"]:
        problems.append(f"{report['faces_facing_the_source']} face(s) oriented "
                        f"back toward the source")
    if report["faces_in_core"] < min_core_faces:
        problems.append(f"only {report['faces_in_core']} faces inside the core "
                        f"(need {min_core_faces})")
    if abs(report["area_error"]) > area_tolerance:
        problems.append(f"meshed area differs from 2 pi R^2 by "
                        f"{100 * report['area_error']:+.2f}%")
    return problems
