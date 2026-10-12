r"""The mesh: a cubic background, refined and carved by OpenFOAM's own tools.

Pipeline::

    blockMesh      a uniform cubic box; x = 0 is ONE patch
    snappyHexMesh  refine the plume volumes, carve the plate out as a 'wall'
    topoSet        select the x = 0 faces inside r <= R0 into a faceSet
    createPatch    move that set into its own 'nozzle' patch
    checkMesh      mesh quality

Nothing but OpenFOAM utilities and analytic primitives: no STL, no external
mesher, nothing to version beside ``case.yaml``.

Why not the graded ``blockMesh`` the Cai family uses
----------------------------------------------------
:mod:`plumetools.cai2012.mesh` grades a tensor-product block lattice away from a
fine core. That works for a nozzle alone in a box. It cannot put a body in the
plume: a tensor-product grading that is fine at the nozzle is fine along a whole
slab of the domain, and a plate tens of diameters downstream then sits in cells
with a 40:1 aspect ratio.

An octree does not have that problem. The background is uniform and **cubic**,
every refinement halves a cell in all three directions, and resolution goes
where it is asked for: in volumes around the plume, and at the plate surface.

Resolution
----------
One length sizes everything::

    cell size at level N  =  mesh.background_cell_size_m / 2**N

``mesh.refinement_regions`` are volumes held at a level -- cylinders coaxial with
the plume, or boxes. The plate surface is held at ``mesh.plate_refinement_level``.
Between two levels snappyHexMesh leaves ``mesh.n_cells_between_levels`` cells of
each intermediate level, so the mesh never jumps by more than a factor of two.

The nozzle patch is a staircase
-------------------------------
A circle has no exact representation on a Cartesian grid. The ``nozzle`` patch
is the set of ``x = 0`` faces whose centres satisfy :math:`y^2+z^2 \le R_0^2`, so
its area approaches :math:`\pi R_0^2` as the exit cell shrinks but never equals
it -- and does not approach it monotonically. Counting lattice squares in a
circle is the Gauss circle problem: the error wanders in sign and size, inside a
bound that falls with the cell. For the R-4D exit it is -1.98% at 12 cells
across the diameter, +0.64% at 24 and -0.02% at 48. The injected mass flow is
off by exactly that fraction, so it is computed and reported, never assumed.

For the set to be predictable at all the whole exit disk has to sit at **one**
refinement level, which :func:`plan` enforces: a region that covers part of the
disk at a finer level than the rest is a configuration error, not a mesh.

:func:`predicted_nozzle_faces` then computes the face count and area the mesh
will have, and :mod:`plumetools.nozzleplume.verify` measures them back out of
the built mesh. When the two disagree the mesh is not the mesh that was planned.

The plate is exact
------------------
The plate is a ``searchableBox`` snapped with ``implicitFeatureSnap``. Measured
on the R-4D plate (34.14 x 11.88 x 0.25 m at 0.125 m cells): every face of the
patch is axis-aligned and its upstream face is 405.5832 m^2, the product of the
two dimensions to the last digit. Without feature snapping the same mesh has
1096 chamfered faces along the plate's edges and a front face 1.4% short --
see :func:`plumetools.foamio.primitives.render_snap_controls`.

The cell count is an estimate
-----------------------------
An octree's cell count is not known until it is built. :func:`estimate_cells`
integrates the requested level over the domain, buffer layers included, and
lands within a few percent of what snappyHexMesh produces -- close enough to
refuse a mistyped level before a half-hour mesh, which is all it is for. The
real count is read back after meshing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from plumetools.foamio.primitives import (
    BoxPatch,
    RefinementRegion,
    SearchableBox,
    SearchableCylinder,
    check_surfaces_are_resolvable,
    fmt,
    header,
    render_box_block_mesh_dict,
    render_mesh_quality_dict,
    render_snappy_dict,
    write_lines,
)
from plumetools.nozzleplume.geometry import PlumeGeometry

#: Dotted path recorded in every dictionary this module writes.
GENERATOR = "plumetools.nozzleplume.mesh"

#: Half-thickness of the ``cylinderToFace`` selector, as a fraction of the exit
#: cell. Thin enough to catch nothing but the ``x = 0`` boundary faces, thick
#: enough that a face centre at exactly ``x = 0`` is unambiguously inside it.
SELECTOR_HALF_THICKNESS_FRACTION = 0.25

#: Where in a finest-level cell the ``locationInMesh`` seed sits, per axis.
#:
#: Not 0.5. A cell centre is interior to *that* cell, but it is a **vertex** of
#: the eight children the cell is split into if it is refined -- and a seed on a
#: vertex is the classic way to make snappyHexMesh keep the wrong region. These
#: fractions are not dyadic, so the seed is strictly inside a cell at every
#: refinement level.
SEED_FRACTIONS = (0.37, 0.41, 0.43)

#: Lattice points per axis used to integrate one refinement layer's volume.
ESTIMATE_SAMPLES = 64


# --------------------------------------------------------------------------- #
# refinement regions
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Region:
    """One refinement volume, in metres.

    Attributes:
        name: the ``geometry`` entry name.
        shape: ``cylinder`` (coaxial with the plume) or ``box``.
        level: the octree level held inside it.
        x_min_m, x_max_m: axial extent [m].
        radius_m: cylinder radius [m].
        y_min_m, y_max_m, z_min_m, z_max_m: box bounds [m].
    """

    name: str
    shape: str
    level: int
    x_min_m: float
    x_max_m: float
    radius_m: float | None = None
    y_min_m: float | None = None
    y_max_m: float | None = None
    z_min_m: float | None = None
    z_max_m: float | None = None

    def bounds(self, grow_m: float = 0.0) -> tuple:
        """Bounding box, grown by ``grow_m`` all round [m]."""
        g = float(grow_m)
        if self.shape == "cylinder":
            r = float(self.radius_m) + g
            return (self.x_min_m - g, self.x_max_m + g, -r, r, -r, r)
        return (self.x_min_m - g, self.x_max_m + g, self.y_min_m - g,
                self.y_max_m + g, self.z_min_m - g, self.z_max_m + g)

    def contains(self, x, y, z, grow_m: float = 0.0):
        """Whether points are inside the volume, grown by ``grow_m`` all round."""
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        z = np.asarray(z, dtype=np.float64)
        x0, x1, y0, y1, z0, z1 = self.bounds(grow_m)
        axial = (x >= x0) & (x <= x1)
        if self.shape == "cylinder":
            return axial & (y * y + z * z <= (float(self.radius_m) + grow_m) ** 2)
        return axial & (y >= y0) & (y <= y1) & (z >= z0) & (z <= z1)

    def searchable(self):
        """The analytic primitive snappyHexMesh tests cell centres against."""
        if self.shape == "cylinder":
            return SearchableCylinder(
                name=self.name, point1=(self.x_min_m, 0.0, 0.0),
                point2=(self.x_max_m, 0.0, 0.0), radius=float(self.radius_m))
        return SearchableBox(
            name=self.name,
            min=(self.x_min_m, self.y_min_m, self.z_min_m),
            max=(self.x_max_m, self.y_max_m, self.z_max_m))

    def describe(self) -> str:
        if self.shape == "cylinder":
            return (f"x {fmt(self.x_min_m)} .. {fmt(self.x_max_m)} m, "
                    f"r <= {fmt(self.radius_m)} m")
        return (f"x {fmt(self.x_min_m)} .. {fmt(self.x_max_m)}, "
                f"y {fmt(self.y_min_m)} .. {fmt(self.y_max_m)}, "
                f"z {fmt(self.z_min_m)} .. {fmt(self.z_max_m)} m")


def regions_from_config(cfg) -> tuple:
    """``mesh.refinement_regions`` as :class:`Region` objects, in file order."""
    return tuple(
        Region(
            name=str(r.name), shape=str(r.shape), level=int(r.level),
            x_min_m=float(r.x_min_m), x_max_m=float(r.x_max_m),
            radius_m=None if r.radius_m is None else float(r.radius_m),
            y_min_m=None if r.y_min_m is None else float(r.y_min_m),
            y_max_m=None if r.y_max_m is None else float(r.y_max_m),
            z_min_m=None if r.z_min_m is None else float(r.z_min_m),
            z_max_m=None if r.z_max_m is None else float(r.z_max_m),
        )
        for r in cfg.mesh.refinement_regions)


# --------------------------------------------------------------------------- #
# the plan
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class MeshPlan:
    """Everything about the mesh, computed before a line of it is written.

    Attributes:
        background_cell_size_m: the level-0 cell [m].
        divisions: ``(nx, ny, nz)`` of the background box.
        regions: the refinement volumes.
        plate_level: octree level at the plate surface; ``None`` with no plate.
        exit_level: the level the whole exit disk sits at.
        n_cells_between_levels: buffer cells between successive levels.
        max_cells: the configured ceiling on the estimate.
        cells_by_level: estimated cell count at each level.
        mean_free_path_m: ``lambda`` at the exit, for the cell/mfp ratio.
        seed_point: ``locationInMesh`` [m].
    """

    background_cell_size_m: float
    divisions: tuple
    regions: tuple
    plate_level: int | None
    exit_level: int
    n_cells_between_levels: int
    max_cells: int
    cells_by_level: dict
    mean_free_path_m: float
    seed_point: tuple

    def cell_size_m(self, level: int) -> float:
        """Cell size at an octree level [m]: the background halved per level."""
        return self.background_cell_size_m / 2 ** int(level)

    @property
    def n_background_cells(self) -> int:
        nx, ny, nz = self.divisions
        return int(nx) * int(ny) * int(nz)

    @property
    def max_level(self) -> int:
        """The finest level anywhere in the mesh."""
        levels = [r.level for r in self.regions]
        if self.plate_level is not None:
            levels.append(self.plate_level)
        return max(levels, default=0)

    @property
    def min_cell_size_m(self) -> float:
        """The smallest cell edge in the mesh [m]. Sets the Courant number."""
        return self.cell_size_m(self.max_level)

    @property
    def exit_cell_size_m(self) -> float:
        """Edge of a cell on the nozzle exit [m]."""
        return self.cell_size_m(self.exit_level)

    @property
    def exit_cell_volume_m3(self) -> float:
        """Volume of a cell on the nozzle exit [m^3]. Exact: the cell is a cube."""
        return self.exit_cell_size_m ** 3

    @property
    def plate_cell_size_m(self) -> float | None:
        """Edge of a cell on the plate surface [m], or ``None`` with no plate."""
        return None if self.plate_level is None else self.cell_size_m(self.plate_level)

    @property
    def estimated_cells(self) -> int:
        """The estimated cell count. See :func:`estimate_cells`."""
        return int(round(sum(self.cells_by_level.values())))

    @property
    def exit_cell_over_mfp(self) -> float:
        """Exit cell size over the exit mean free path. DSMC wants at most 1."""
        return self.exit_cell_size_m / self.mean_free_path_m

    @property
    def needs_snappy(self) -> bool:
        """Whether there is anything for snappyHexMesh to do."""
        return bool(self.regions) or self.plate_level is not None

    def as_dict(self) -> dict:
        nx, ny, nz = self.divisions
        return {
            "background_cell_size_m": self.background_cell_size_m,
            "background_divisions": [int(nx), int(ny), int(nz)],
            "n_background_cells": self.n_background_cells,
            "estimated_cells": self.estimated_cells,
            "estimated_cells_by_level": {
                int(level): int(round(count))
                for level, count in sorted(self.cells_by_level.items())},
            "max_cells": int(self.max_cells),
            "max_level": self.max_level,
            "min_cell_size_m": self.min_cell_size_m,
            "exit_level": self.exit_level,
            "exit_cell_size_m": self.exit_cell_size_m,
            "exit_cell_over_mean_free_path": self.exit_cell_over_mfp,
            "plate_level": self.plate_level,
            "plate_cell_size_m": self.plate_cell_size_m,
            "refinement_regions": [
                {"name": r.name, "shape": r.shape, "level": r.level,
                 "cell_size_m": self.cell_size_m(r.level)} for r in self.regions],
            "location_in_mesh": [float(v) for v in self.seed_point],
        }

    def describe(self, geom: PlumeGeometry) -> list[str]:
        nx, ny, nz = self.divisions
        lines = [
            "Mesh plan (cubic background, octree refinement by snappyHexMesh)",
            f"  background               {nx} x {ny} x {nz} = "
            f"{self.n_background_cells:,} cells of "
            f"{fmt(self.background_cell_size_m)} m",
        ]
        if self.regions:
            lines.append("  refinement regions")
            width = max(len(r.name) for r in self.regions)
            for r in self.regions:
                lines.append(
                    f"    {r.name:<{width}}  level {r.level}  "
                    f"{fmt(self.cell_size_m(r.level)):>9} m   {r.shape}: "
                    f"{r.describe()}")
        else:
            lines.append("  refinement regions       none")
        if self.plate_level is not None:
            across = geom.plate.thickness_m / self.plate_cell_size_m
            lines.append(
                f"  plate surface            level {self.plate_level}  "
                f"{fmt(self.plate_cell_size_m)} m   ({across:.1f} cells through "
                f"its {fmt(geom.plate.thickness_m)} m thickness)")
        lines += [
            f"  exit cells               level {self.exit_level}  "
            f"{fmt(self.exit_cell_size_m)} m   "
            f"= {geom.diameter_m / self.exit_cell_size_m:.1f} across D, "
            f"{self.exit_cell_over_mfp:.3g} exit mean free paths",
            f"  smallest cell            {fmt(self.min_cell_size_m)} m"
            f"   (level {self.max_level})",
            f"  estimated cells          {self.estimated_cells:,}"
            f"   (ceiling {self.max_cells:,}); by level: "
            + ", ".join(f"{level}: {int(round(count)):,}"
                        for level, count in sorted(self.cells_by_level.items())),
            f"  locationInMesh           ({fmt(self.seed_point[0])} "
            f"{fmt(self.seed_point[1])} {fmt(self.seed_point[2])})",
        ]
        return lines


def refinement_layers(geom: PlumeGeometry, background_cell_size_m: float,
                      regions, plate_level: int | None,
                      n_cells_between_levels: int) -> list:
    """Every requested refinement as a stack of nested volumes.

    Returns:
        ``[(level, contains(x, y, z), bounds, source)]``: a volume, the level
        cells inside it are held at, its bounding box, and the name of the
        region (or ``"the plate"``) that asked for it.

    Each feature contributes itself at its own level, then one shell per coarser
    level around it, ``n_cells_between_levels`` cells of that level thick -- the
    buffer snappyHexMesh adds so that levels change one at a time.

    A region and a surface differ in what "the feature at its own level" is. A
    region is exactly the cells whose centres are inside it. A **surface** is
    the cells it cuts *plus* ``n_cells_between_levels`` more at the same level
    on each side -- a skin ``n + 1`` cells thick, measured on the R-4D plate.
    """
    bg = float(background_cell_size_m)
    n_between = int(n_cells_between_levels)

    def size(level: int) -> float:
        return bg / 2 ** int(level)

    layers = []

    def add(level: int, contains, bounds, source: str) -> None:
        grow = 0.0
        layers.append((int(level), lambda x, y, z, g=grow: contains(x, y, z, g),
                       bounds(grow), source))
        for lower in range(int(level) - 1, 0, -1):
            grow += n_between * size(lower)
            layers.append((lower, lambda x, y, z, g=grow: contains(x, y, z, g),
                           bounds(grow), source))

    for region in regions:
        add(region.level, region.contains, region.bounds, f"region {region.name!r}")

    plate = geom.plate
    if plate is not None and plate_level is not None:
        skin = (n_between + 1) * size(plate_level)

        def plate_bounds(grow: float) -> tuple:
            x0, x1, y0, y1, z0, z1 = plate.bounds
            g = skin + grow
            return (x0 - g, x1 + g, y0 - g, y1 + g, z0 - g, z1 + g)

        add(plate_level,
            lambda x, y, z, g: plate.contains(x, y, z, skin + g), plate_bounds,
            "the plate")
    return layers


def exit_plane_cells(regions, geom: PlumeGeometry,
                     background_cell_size_m: float) -> list:
    """The first layer of cells over the exit disk, by snappyHexMesh's own rule.

    Returns:
        ``[(level, (x, y, z))]`` -- the level and centre of every cell that
        touches the exit plane and overlaps the modelled part of the disk.

    snappyHexMesh refines a cell when its **centre** is inside a region whose
    level it has not reached, and repeats until nothing changes. This does
    exactly that, starting from the background cells over the disk and following
    only the children that still overlap it.

    It is the rule itself rather than a margin standing in for it, because the
    rule has a consequence no margin captures: a region **smaller than the cells
    around it is never seen**. A 0.3 m cylinder in a 0.5 m background contains
    no cell centre, so nothing is ever refined, however fine the level it asks
    for.
    """
    bg = float(background_cell_size_m)
    radius_sq = geom.radius_m ** 2

    def overlaps(level: int, j: int, k: int) -> bool:
        """Whether the cell's footprint reaches inside the disk."""
        h = bg / 2 ** level
        y0, z0 = geom.y_min_m + j * h, geom.z_min_m + k * h
        # The point of the square nearest the axis.
        y = min(max(0.0, y0), y0 + h)
        z = min(max(0.0, z0), z0 + h)
        return y * y + z * z < radius_sq

    n_y = round((geom.y_max_m - geom.y_min_m) / bg)
    n_z = round((geom.z_max_m - geom.z_min_m) / bg)
    stack = [(0, j, k) for j in range(n_y) for k in range(n_z) if overlaps(0, j, k)]

    leaves = []
    while stack:
        level, j, k = stack.pop()
        h = bg / 2 ** level
        centre = (0.5 * h, geom.y_min_m + (j + 0.5) * h, geom.z_min_m + (k + 0.5) * h)
        wanted = max((r.level for r in regions if bool(r.contains(*centre))),
                     default=0)
        if level < wanted:
            stack += [(level + 1, 2 * j + dj, 2 * k + dk)
                      for dj in (0, 1) for dk in (0, 1)
                      if overlaps(level + 1, 2 * j + dj, 2 * k + dk)]
        else:
            leaves.append((level, centre))
    return leaves


def exit_level(regions, geom: PlumeGeometry, background_cell_size_m: float,
               layers=()) -> int:
    """The single refinement level the nozzle exit disk sits at.

    Args:
        regions: the refinement volumes.
        geom: the geometry.
        background_cell_size_m: the level-0 cell [m].
        layers: :func:`refinement_layers`, for the buffer check. Empty skips it.

    Returns:
        The level -- 0 if no region refines the exit.

    Raises:
        ValueError: if the exit would be meshed at more than one cell size. The
            staircase nozzle would then not be the one predicted, and the
            particle occupancy -- which is quoted for "the exit cell" -- would
            be right for only part of the inlet.

    Two things can split the exit. A region can reach only part of it, which
    :func:`exit_plane_cells` finds by applying snappyHexMesh's rule. Or a finer
    feature nearby -- another region, or the plate -- can spill its **buffer
    layers** onto the exit; those refine cells no region asked for, so they are
    checked separately against the cells the first pass produced.
    """
    cells = exit_plane_cells(regions, geom, background_cell_size_m)
    levels = sorted({level for level, _ in cells})
    if len(levels) > 1:
        finest = [r.name for r in regions if r.level == levels[-1]]
        raise ValueError(
            f"mesh.refinement_regions: the nozzle exit would be meshed at "
            f"refinement levels {levels}, not at one. snappyHexMesh refines a "
            f"cell when its CENTRE is inside a region, so the region that "
            f"resolves the exit ({finest}, level {levels[-1]}) has to contain "
            f"the centre of every cell over the disk (radius {geom.radius_m:g} m)"
            f" at each level on the way down. Give it a radius a cell or two "
            f"beyond the disk and x_min_m <= 0, or move it clear of the exit.")
    level = levels[0] if levels else 0
    if not cells or not layers:
        return level

    x, y, z = (np.array(axis) for axis in zip(*(centre for _, centre in cells)))
    for layer_level, contains, _, source in layers:
        if layer_level > level and bool(np.any(contains(x, y, z))):
            raise ValueError(
                f"mesh: the buffer layers of {source} reach the nozzle exit at "
                f"level {layer_level}, above the level {level} the exit is "
                f"otherwise meshed at, so the exit would be at two cell sizes. "
                f"A region smaller than the cells around it is never seen -- "
                f"snappyHexMesh tests cell CENTRES -- and only its buffer "
                f"arrives. Cover the whole exit disk (radius {geom.radius_m:g} "
                f"m) with a region at level {layer_level} or finer whose radius "
                f"reaches the surrounding cell centres, or move {source} away "
                f"from the exit.")
    return level


def estimate_cells(geom: PlumeGeometry, background_cell_size_m: float, regions,
                   plate_level: int | None, n_cells_between_levels: int) -> dict:
    """Estimate the cell count snappyHexMesh will produce, level by level.

    Args:
        geom: the geometry.
        background_cell_size_m: the level-0 cell [m].
        regions: the refinement volumes.
        plate_level: level at the plate surface, or ``None``.
        n_cells_between_levels: buffer cells between successive levels.

    Returns:
        ``{level: cells}``, as floats.

    The volumes of :func:`refinement_layers` -- each feature at its own level,
    wrapped in its buffer shells -- are claimed finest-first, so an overlap is
    counted once at its finest level, and what nobody claims is background.

    The plate's skin is where an estimate goes wrong. snappyHexMesh keeps
    ``n_cells_between_levels + 1`` cells at the surface level on each side of a
    surface, not one, and getting that wrong is most of the error there is to
    have: on the R-4D plate (811 m^2 of surface at 0.125 m) a one-cell skin
    undercounts the mesh by a third, 360 000 cells against the 529 126
    snappyHexMesh built, while this model gives 542 000 -- 2.5% over, and within
    8% at every level.

    Each volume is integrated on a lattice of :data:`ESTIMATE_SAMPLES` points
    per axis over its own bounding box, so a small fine region costs the same to
    integrate as a large coarse one and neither is under-sampled.
    """
    bg = float(background_cell_size_m)
    plate = geom.plate

    def size(level: int) -> float:
        return bg / 2 ** int(level)

    layers = refinement_layers(geom, bg, regions, plate_level,
                               n_cells_between_levels)

    # Finest first; the sort is stable, so equal levels keep file order.
    layers.sort(key=lambda layer: -layer[0])

    counts: dict = {}
    claimed = 0.0
    domain = geom.bounds
    for index, (level, contains, bounds, _) in enumerate(layers):
        lo = [max(bounds[2 * a], domain[2 * a]) for a in range(3)]
        hi = [min(bounds[2 * a + 1], domain[2 * a + 1]) for a in range(3)]
        if any(h <= l for l, h in zip(lo, hi)):
            continue
        axes = [l + (np.arange(ESTIMATE_SAMPLES) + 0.5) * (h - l) / ESTIMATE_SAMPLES
                for l, h in zip(lo, hi)]
        x, y, z = (a.ravel() for a in np.meshgrid(*axes, indexing="ij"))
        mine = contains(x, y, z)
        for _, earlier, _, _ in layers[:index]:
            if not mine.any():
                break
            mine &= ~earlier(x, y, z)
        volume = float(mine.mean()) * math.prod(h - l for l, h in zip(lo, hi))
        counts[level] = counts.get(level, 0.0) + volume / size(level) ** 3
        claimed += volume

    counts[0] = counts.get(0, 0.0) + max(0.0, geom.domain_volume_m3 - claimed) / bg ** 3

    # The plate itself is not fluid: its cells were counted above, at the plate
    # level, and snappyHexMesh removes them.
    if plate is not None and plate_level is not None:
        x0, x1, y0, y1, z0, z1 = geom.modelled_plate_bounds
        solid = (x1 - x0) * (y1 - y0) * (z1 - z0) / size(plate_level) ** 3
        counts[plate_level] = max(0.0, counts.get(plate_level, 0.0) - solid)

    return dict(sorted(counts.items()))


def location_in_mesh(geom: PlumeGeometry, finest_cell_size_m: float) -> tuple:
    """A seed point in the fluid: inside the box, outside the plate, on no face.

    snappyHexMesh keeps the region connected to this point and discards the
    rest, so a seed inside the plate, or on a cell face, edge or vertex, gives
    the wrong region or an outright failure.

    The point is placed between the exit and the plate, off the plume axis, and
    then moved to a fixed non-dyadic position inside the **finest-level** cell
    containing it (:data:`SEED_FRACTIONS`). Lattice planes at every level are
    whole multiples of the finest cell, so that position is on none of them
    whatever refinement the surrounding cells end up with.

    Raises:
        ValueError: if the point lands in the plate or on the boundary -- which
            takes a plate closer to the exit than one finest cell.
    """
    h = float(finest_cell_size_m)
    candidate = (
        0.5 * geom.flow_length_m,
        geom.y_min_m + 0.31 * (geom.y_max_m - geom.y_min_m),
        geom.z_min_m + 0.29 * (geom.z_max_m - geom.z_min_m),
    )
    lows = (geom.x_min_m, geom.y_min_m, geom.z_min_m)
    point = tuple(
        lo + (math.floor((value - lo) / h) + fraction) * h
        for value, lo, fraction in zip(candidate, lows, SEED_FRACTIONS))

    x, y, z = point
    if not (geom.x_min_m < x < geom.x_max_m and geom.y_min_m < y < geom.y_max_m
            and geom.z_min_m < z < geom.z_max_m):
        raise ValueError(f"locationInMesh {point} is on or outside the domain boundary")
    if geom.plate is not None and bool(geom.plate.contains(x, y, z, margin_m=h)):
        raise ValueError(
            f"locationInMesh {point} is inside the plate; snappyHexMesh would "
            f"keep the plate and discard the fluid. The plate is too close to "
            f"the exit plane to seed the gap between them.")
    return point


def plan(cfg, geom: PlumeGeometry, exit_state) -> MeshPlan:
    """Work out the whole mesh from the config, the geometry and the exit state.

    Raises:
        ValueError: for a refinement region outside the domain, an exit disk
            split across two levels, a plate too thin for its refinement level,
            or an estimated cell count over ``mesh.max_cells``. Each names the
            key to change.
    """
    bg = float(cfg.mesh.background_cell_size_m)
    x0, x1, y0, y1, z0, z1 = geom.bounds
    # Exact: geometry.from_config rounded every extent to the lattice.
    divisions = (round((x1 - x0) / bg), round((y1 - y0) / bg), round((z1 - z0) / bg))

    regions = regions_from_config(cfg)
    for region in regions:
        rx0, rx1, ry0, ry1, rz0, rz1 = region.bounds()
        if rx1 <= x0 or rx0 >= x1 or ry1 <= y0 or ry0 >= y1 or rz1 <= z0 or rz0 >= z1:
            raise ValueError(
                f"mesh.refinement_regions: {region.name!r} lies entirely outside "
                f"the domain and would refine nothing")

    plate_level = None
    if geom.plate is not None:
        plate_level = int(cfg.mesh.plate_refinement_level)
        try:
            check_surfaces_are_resolvable([_plate_surface(cfg, geom, plate_level)], bg)
        except ValueError as exc:
            raise ValueError(f"mesh.plate_refinement_level: {exc}") from None

    n_between = int(cfg.mesh.n_cells_between_levels)
    level = exit_level(
        regions, geom, bg,
        refinement_layers(geom, bg, regions, plate_level, n_between))
    cells = estimate_cells(geom, bg, regions, plate_level, n_between)

    total = int(round(sum(cells.values())))
    ceiling = int(cfg.mesh.max_cells)
    if total > ceiling:
        breakdown = ", ".join(
            f"level {lvl} ({bg / 2 ** lvl:g} m): {int(round(count)):,}"
            for lvl, count in cells.items())
        raise ValueError(
            f"the mesh is estimated at {total:,} cells, over mesh.max_cells = "
            f"{ceiling:,}. By level -- {breakdown}. Each level costs eight times "
            f"the one below per unit volume, so shrink or drop the finest "
            f"regions first; or raise mesh.max_cells if this is intended.")

    max_level = max([r.level for r in regions]
                    + ([plate_level] if plate_level is not None else []), default=0)
    return MeshPlan(
        background_cell_size_m=bg,
        divisions=tuple(int(n) for n in divisions),
        regions=regions,
        plate_level=plate_level,
        exit_level=level,
        n_cells_between_levels=n_between,
        max_cells=ceiling,
        cells_by_level=cells,
        mean_free_path_m=float(exit_state.mean_free_path_m),
        seed_point=location_in_mesh(geom, bg / 2 ** max_level),
    )


# --------------------------------------------------------------------------- #
# the nozzle patch
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class NozzlePrediction:
    """The staircase nozzle patch the mesh will produce.

    ``area_m2`` and ``exact_area_m2`` are both of the **modelled** part of the
    exit -- half or a quarter of it under symmetry.
    """

    n_faces: int
    area_m2: float
    exact_area_m2: float
    cells_across_diameter: float

    @property
    def area_error(self) -> float:
        """Signed relative area error, ``(meshed - exact) / exact``."""
        return (self.area_m2 - self.exact_area_m2) / self.exact_area_m2

    def as_dict(self) -> dict:
        return {
            "n_faces": self.n_faces,
            "area_m2": self.area_m2,
            "exact_area_m2": self.exact_area_m2,
            "area_error": self.area_error,
            "cells_across_diameter": self.cells_across_diameter,
        }

    def describe(self) -> list[str]:
        return [
            "Nozzle patch (a staircase on a Cartesian grid)",
            f"  faces                    {self.n_faces}",
            f"  cells across D           {self.cells_across_diameter:.1f}",
            f"  area                     {self.area_m2:.6e} m^2  against the "
            f"exact {self.exact_area_m2:.6e} m^2",
            f"  area error               {100.0 * self.area_error:+.2f}%"
            f"   -- the injected mass flow is off by the same amount",
        ]


def predicted_nozzle_faces(plan_: MeshPlan, geom: PlumeGeometry) -> NozzlePrediction:
    """Which ``x = 0`` faces will land in the nozzle patch, and their area.

    The exit sits at one refinement level, so the ``x = 0`` faces across it are
    squares of the exit cell size on a lattice that starts at the domain's
    ``y_min`` and ``z_min``. This applies ``topoSet``'s ``cylinderToFace`` rule
    -- face centre inside ``r <= R0`` -- to that lattice.

    The mesh has not been built yet, so this is a *prediction*;
    :mod:`plumetools.nozzleplume.verify` measures the built mesh and compares.
    """
    h = plan_.exit_cell_size_m
    reach = geom.radius_m + h

    def centres(lo: float, hi: float) -> np.ndarray:
        first = max(0, math.floor((-reach - lo) / h))
        last = min(round((hi - lo) / h), math.ceil((reach - lo) / h))
        return lo + (np.arange(first, last) + 0.5) * h

    yy, zz = np.meshgrid(centres(geom.y_min_m, geom.y_max_m),
                         centres(geom.z_min_m, geom.z_max_m), indexing="ij")
    count = int(np.count_nonzero(geom.is_on_nozzle(yy, zz)))
    return NozzlePrediction(
        n_faces=count,
        area_m2=count * h * h,
        exact_area_m2=geom.exact_nozzle_area_m2,
        cells_across_diameter=geom.diameter_m / h,
    )


# --------------------------------------------------------------------------- #
# dictionaries
# --------------------------------------------------------------------------- #

def _plate_surface(cfg, geom: PlumeGeometry, level: int) -> SearchableBox:
    """The plate as the ``searchableBox`` snappyHexMesh carves."""
    plate = geom.plate
    return SearchableBox(
        name=cfg.mesh.patch_names["plate"],
        level=int(level),
        # 'wall', never 'patch'. A patch DELETES the molecules that reach it; a
        # wall runs hitWallPatch, which records fD and q and applies the
        # gas-surface model. A plate declared 'patch' would absorb the plume and
        # report no load.
        patch_type="wall",
        # The WHOLE plate, also under symmetry: the domain boundary cuts it at
        # the symmetry plane, which leaves no spurious face there. A half-width
        # box would end exactly on that plane and leave one.
        min=(plate.x_upstream_m, plate.y_min_m, plate.z_min_m),
        max=(plate.x_downstream_m, plate.y_max_m, plate.z_max_m),
    )


def render_block_mesh_dict(cfg, geom: PlumeGeometry, plan_: MeshPlan) -> str:
    """The uniform cubic background box.

    ``x = 0`` is written as a **single** patch, ``upstreamVacuum``. The nozzle
    is carved out of it after snappyHexMesh has refined it, because a disk is
    not a block face and the faces it is made of do not exist until then.
    """
    names = cfg.mesh.patch_names
    nx, ny, nz = plan_.divisions
    lows = {"y": "y_min", "z": "z_min"}
    symmetric = tuple(lows[axis] for axis in ("y", "z")
                      if axis in geom.symmetry_planes)
    outer = tuple(face for face in ("x_max", "y_min", "y_max", "z_min", "z_max")
                  if face not in symmetric)

    patches = [BoxPatch(
        names["upstream_vacuum"], "patch", ("x_min",),
        "x = 0, outside the nozzle disk. An OPEN boundary: a molecule that "
        "comes back through the exit plane has left the domain. NOT a symmetry "
        "plane.")]
    if symmetric:
        patches.append(BoxPatch(
            names["symmetry"], "symmetry", symmetric,
            "the modelled plane(s) of symmetry; molecules reflect specularly"))
    patches.append(BoxPatch(
        names["outer"], "patch", outer,
        "the vacuum far field. 'patch', never 'wall': a molecule reaching it is "
        "DELETED, which is what expanding into vacuum means."))

    return render_box_block_mesh_dict(
        bounds=geom.bounds,
        divisions=(nx, ny, nz),
        patches=patches,
        generator=GENERATOR,
        note=(f"Background box: {nx}x{ny}x{nz} = {plan_.n_background_cells} cubic "
              f"cells of {fmt(plan_.background_cell_size_m)} m. Refinement and "
              f"the plate come from snappyHexMesh; the nozzle from topoSet."),
    )


def render_snappy_hex_mesh_dict(cfg, geom: PlumeGeometry, plan_: MeshPlan) -> str:
    """The refinement volumes and the plate, with castellate and snap controls."""
    surfaces = []
    carving = "refinement only -- nothing is carved"
    if geom.plate is not None:
        surfaces.append(_plate_surface(cfg, geom, plan_.plate_level))
        carving = (f"carve the plate ({fmt(geom.plate.size_y_m)} x "
                   f"{fmt(geom.plate.size_z_m)} x {fmt(geom.plate.thickness_m)} m "
                   f"at x = {fmt(geom.plate.x_upstream_m)} m) at "
                   f"{fmt(plan_.plate_cell_size_m)} m")

    return render_snappy_dict(
        surfaces,
        location_in_mesh=plan_.seed_point,
        n_cells_between_levels=plan_.n_cells_between_levels,
        generator=GENERATOR,
        # A backstop, not a budget: an octree that reaches maxGlobalCells stops
        # refining and yields a quietly under-resolved mesh. The estimate has
        # already been held to mesh.max_cells, and verify re-checks the real
        # count, so this only has to stop a runaway.
        max_global_cells=2 * plan_.max_cells,
        refinement_regions=[RefinementRegion(surface=r.searchable(), level=r.level)
                            for r in plan_.regions],
        # The plate is a box, so it is nothing but feature edges.
        implicit_feature_snap=geom.plate is not None,
        note=(f"{len(plan_.regions)} refinement region(s) down to "
              f"{fmt(plan_.min_cell_size_m)} m; {carving}."),
        preamble=[
            "// Analytic primitives, no STL: snappyHexMesh snaps to the true surface",
            "// rather than to a triangulation of it, and nothing outside the",
            "// repository has to be generated or versioned.",
            "//",
            "// Refinement regions are volumes held at an octree level; they create",
            "// no patch. The plate creates one, named after its geometry entry, and",
            "// its patchInfo type is load-bearing: 'wall' runs hitWallPatch, which",
            "// records fD/q and applies the WallInteractionModel. A plate declared",
            "// 'patch' would absorb the plume and feel no load.",
        ],
    )


def render_topo_set_dict(cfg, geom: PlumeGeometry, plan_: MeshPlan) -> str:
    """``system/topoSetDict``: the ``x = 0`` faces inside ``r <= R0``.

    Two actions, and the order is load-bearing. ``patchToFace`` first restricts
    the set to the ``x = 0`` boundary; only then does ``cylinderToFace`` subset
    it by radius. Run the other way round, ``cylinderToFace`` would also select
    every **internal** face inside the cylinder, and ``createPatch`` would be
    asked to turn interior faces into a boundary.
    """
    names = cfg.mesh.patch_names
    half = SELECTOR_HALF_THICKNESS_FRACTION * plan_.exit_cell_size_m

    lines = header(
        "topoSetDict", GENERATOR,
        f"Select the exit disk r <= {fmt(geom.radius_m)} m out of "
        f"{names['upstream_vacuum']}.")
    lines += [
        "// Order matters: patchToFace restricts the set to the x = 0 boundary,",
        "// and only then does cylinderToFace subset it by radius. Reversed, the",
        "// cylinder would also pick up every INTERNAL face inside it.",
        "//",
        "// The selector is a thin disc straddling x = 0 rather than a long",
        "// cylinder, so it cannot reach a boundary face anywhere else even if",
        "// the first action is ever changed.",
        "",
        "actions",
        "(",
        "    {",
        "        name    nozzleFaces;",
        "        type    faceSet;",
        "        action  new;",
        "        source  patchToFace;",
        f"        patch   {names['upstream_vacuum']};",
        "    }",
        "",
        "    {",
        "        name    nozzleFaces;",
        "        type    faceSet;",
        "        action  subset;",
        "        source  cylinderToFace;",
        f"        point1  ({fmt(-half)} 0 0);",
        f"        point2  ({fmt(half)} 0 0);",
        f"        radius  {fmt(geom.radius_m)};",
        "    }",
        ");",
        "",
        "// ************************************************************************* //",
    ]
    return "\n".join(lines) + "\n"


def render_create_patch_dict(cfg) -> str:
    """``system/createPatchDict``: move ``nozzleFaces`` into the nozzle patch.

    ``type patch``, never ``wall``. ``plumeFieldInflow`` refuses outright to
    inject across a ``wallPolyPatch``, so a mis-typed nozzle is a hard error
    rather than a silent physical change.
    """
    names = cfg.mesh.patch_names
    lines = header(
        "createPatchDict", GENERATOR,
        f"Promote the exit disk to its own patch, {names['nozzle']}.")
    lines += [
        "// pointSync false: the faces are already boundary faces of one patch,",
        "// so no point addressing changes and syncing would only cost a pass.",
        "pointSync   false;",
        "",
        "patches",
        "(",
        "    {",
        f"        name            {names['nozzle']};",
        "        patchInfo",
        "        {",
        "            // 'patch', never 'wall'. plumeFieldInflow aborts if an",
        "            // inflow patch is a wallPolyPatch, and rightly: a wall is a",
        "            // gas-surface interaction site, not a source.",
        "            type        patch;",
        "        }",
        "        constructFrom   set;",
        "        set             nozzleFaces;",
        "    }",
        ");",
        "",
        "// ************************************************************************* //",
    ]
    return "\n".join(lines) + "\n"


def write_mesh_setup(case_dir: Path, cfg, geom: PlumeGeometry,
                     plan_: MeshPlan) -> list[Path]:
    """Write every mesh dictionary, in the order the pipeline consumes them.

    ``snappyHexMeshDict`` is written only when there is something for
    snappyHexMesh to do, and a stale one is **removed** otherwise: ``./Allmesh``
    decides whether to run snappyHexMesh by whether the file is there, so a
    leftover from an earlier configuration would carve a plate ``case.yaml`` no
    longer has.
    """
    system = Path(case_dir) / "system"
    written = [write_lines(system / "blockMeshDict",
                           render_block_mesh_dict(cfg, geom, plan_).splitlines())]

    snappy = system / "snappyHexMeshDict"
    if plan_.needs_snappy:
        written.append(write_lines(
            snappy, render_snappy_hex_mesh_dict(cfg, geom, plan_).splitlines()))
    elif snappy.exists():
        snappy.unlink()

    written += [
        write_lines(system / "meshQualityDict",
                    render_mesh_quality_dict(GENERATOR).splitlines()),
        write_lines(system / "topoSetDict",
                    render_topo_set_dict(cfg, geom, plan_).splitlines()),
        write_lines(system / "createPatchDict",
                    render_create_patch_dict(cfg).splitlines()),
    ]
    return written
