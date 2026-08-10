r"""The graded multi-block Cartesian mesh, and the cost of the one Cai used.

Pipeline::

    blockMesh     the graded multi-block box; x = 0 is ONE patch
    topoSet       select the x = 0 faces inside r <= R0 into a faceSet
    createPatch   move that set into its own 'nozzle' patch
    checkMesh     mesh quality

Why not a literal translation of Cai's grid
-------------------------------------------
Cai's grid is uniform with :math:`\Delta x = \Delta y = \lambda_0`, referred to
the Kn = 0.01 exit properties -- 2 mm for D = 0.2 m. He could afford that because
his solver is **axisymmetric**: the mesh is a 2-D ``(x, r)`` sheet.

Translated literally into the 3-D Cartesian box this case asks for
(:math:`0 \le X/D \le 10`, :math:`|Y|,|Z|/D \le 10`), 2 mm cells give

.. code-block:: text

    1000 x 2000 x 2000 = 4.0e9 cells

and at least as many parcels. :func:`uniform_cost` computes that number for
whatever configuration is loaded and ``./Allmesh`` prints it *before* writing
anything, which is the point: the impracticality is a reported result, not a
thing discovered halfway through a mesh.

What is built instead
---------------------
A tensor-product multi-block box: a **uniform fine core** wrapped around the
nozzle and the near plume, and one expanding segment on each outward side.

.. code-block:: text

      z
      ^   +--------+---------------------------+
      |   |        |                           |   expanding
      |   +--------+---------------------------+
      |   | core   |  expanding in x           |
    --+---#========+===========================+--> x
      |   | (fine, |                           |
      |   +--------+---------------------------+
      |   |  uniform)                          |   expanding
          +--------+---------------------------+
         x=0    core_x                       x_max

Two properties are enforced rather than hoped for:

* **the exit disk lies entirely inside the uniform core**, so every injecting
  face is the same size and the staircase nozzle is as round as the core
  resolution allows (``mesh.core_half_over_D > 0.5`` is a config error);
* **the first cell of each expanding segment matches the core cell**, because
  the expansion ratio is *solved for* rather than typed in --
  :func:`geometric_segment`. A jump in cell size at the core boundary would
  appear in the sampled density as a step that looks like physics.

The nozzle patch is a staircase
-------------------------------
A circle has no exact representation on a Cartesian grid. The ``nozzle`` patch is
the set of ``x = 0`` faces whose centres satisfy :math:`y^2+z^2 \le R_0^2`, so its
area approaches :math:`\pi R_0^2` as the core cell shrinks but never equals it.
:func:`predicted_nozzle_faces` computes the count and the area the mesh will
have, and :mod:`plumetools.cai2012.checks` measures them back out of the built
mesh and fails the case if they disagree with the disk by more than
``checks.max_nozzle_area_error``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from plumetools.cai2012.geometry import CaiGeometry
from plumetools.foamio.primitives import (
    BOX_FACES,
    fmt,
    fmt_point,
    header,
    write_lines,
)

#: Dotted path recorded in every dictionary this module writes.
GENERATOR = "plumetools.cai2012.mesh"

#: Half-thickness of the ``cylinderToFace`` selector, as a fraction of the
#: smallest cell. The selector must be thin enough to catch nothing but the
#: ``x = 0`` boundary faces -- it is applied to a set already restricted to that
#: patch, so this only has to be non-degenerate -- and thick enough that a face
#: centre at exactly ``x = 0`` is unambiguously inside it.
SELECTOR_HALF_THICKNESS_FRACTION = 0.25


# --------------------------------------------------------------------------- #
# one graded segment
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Segment:
    """One axis interval of the block lattice.

    Attributes:
        lo, hi: the interval [m].
        n_cells: cells across it.
        expansion: ``simpleGrading`` ratio, last cell over first, **in the local
            +axis direction**. ``1.0`` is uniform.
        first_cell_m, last_cell_m: the resulting cell sizes [m].
    """

    lo: float
    hi: float
    n_cells: int
    expansion: float
    first_cell_m: float
    last_cell_m: float

    @property
    def length_m(self) -> float:
        return self.hi - self.lo

    @property
    def min_cell_m(self) -> float:
        return min(self.first_cell_m, self.last_cell_m)

    @property
    def max_cell_m(self) -> float:
        return max(self.first_cell_m, self.last_cell_m)

    def cell_edges(self) -> np.ndarray:
        """The ``n_cells + 1`` cell boundaries [m], including both ends."""
        n = int(self.n_cells)
        if n == 1:
            return np.array([self.lo, self.hi])
        q = float(self.expansion) ** (1.0 / (n - 1))
        if abs(q - 1.0) < 1e-12:
            return np.linspace(self.lo, self.hi, n + 1)
        widths = self.first_cell_m * q ** np.arange(n)
        edges = np.empty(n + 1)
        edges[0] = self.lo
        edges[1:] = self.lo + np.cumsum(widths)
        edges[-1] = self.hi  # exact, against accumulated round-off
        return edges

    def cell_centres(self) -> np.ndarray:
        edges = self.cell_edges()
        return 0.5 * (edges[:-1] + edges[1:])


def uniform_segment(lo: float, hi: float, cell_size_m: float) -> Segment:
    """A uniform segment of at least one cell, no coarser than ``cell_size_m``.

    The count is rounded **up**, so the realised cell is never larger than
    asked for -- a mesh that is silently 3% coarser than the DSMC criterion
    demanded is the kind of thing that is never noticed.
    """
    length = float(hi) - float(lo)
    if not length > 0.0:
        raise ValueError(f"segment {lo} .. {hi} is empty or inverted")
    n = max(1, int(math.ceil(length / float(cell_size_m) - 1e-9)))
    width = length / n
    return Segment(float(lo), float(hi), n, 1.0, width, width)


def geometric_segment(lo: float, hi: float, n_cells: int, first_cell_m: float,
                      max_expansion: float) -> Segment:
    r"""An expanding segment whose first cell matches ``first_cell_m`` if it can.

    Args:
        lo, hi: the interval [m]. Expansion runs from ``lo`` to ``hi``.
        n_cells: cells across it.
        first_cell_m: the cell size to match at the ``lo`` end -- the core cell.
        max_expansion: cap on the total ratio, last cell over first.

    Returns:
        A :class:`Segment` whose ``expansion`` is the ``simpleGrading`` ratio.

    The geometric progression with common ratio :math:`q` and first cell
    :math:`\delta` spans :math:`\delta(q^n-1)/(q-1)`. Setting that equal to the
    segment length and solving for :math:`q` is what makes the mesh continuous
    across the core boundary. The equation is monotone in :math:`q`, so a
    bisection is both sufficient and certain.

    If even ``max_expansion`` cannot span the segment, the cap wins: the ratio is
    clamped and the first cell comes out **larger** than ``first_cell_m``. That
    is a visible mesh jump, so :meth:`MeshPlan.describe` reports the ratio and
    :mod:`plumetools.cai2012.checks` warns on it, rather than either failing
    (the mesh is still usable) or silently ignoring the cap.
    """
    length = float(hi) - float(lo)
    if not length > 0.0:
        raise ValueError(f"segment {lo} .. {hi} is empty or inverted")
    n = int(n_cells)
    if n < 1:
        raise ValueError(f"n_cells must be >= 1, got {n_cells}")
    delta = float(first_cell_m)
    cap = max(1.0, float(max_expansion))

    if n == 1:
        return Segment(float(lo), float(hi), 1, 1.0, length, length)

    def span(q: float) -> float:
        if abs(q - 1.0) < 1e-14:
            return delta * n
        return delta * (q ** n - 1.0) / (q - 1.0)

    q_cap = cap ** (1.0 / (n - 1))

    if span(1.0) >= length:
        # A uniform segment already overshoots: the core cell is bigger than the
        # segment needs, so just divide it evenly.
        width = length / n
        return Segment(float(lo), float(hi), n, 1.0, width, width)

    if span(q_cap) < length:
        # The cap cannot reach; clamp and let the first cell grow.
        first = length * (q_cap - 1.0) / (q_cap ** n - 1.0)
        return Segment(float(lo), float(hi), n, q_cap ** (n - 1),
                       first, first * q_cap ** (n - 1))

    low, high = 1.0, q_cap
    for _ in range(200):
        mid = 0.5 * (low + high)
        if span(mid) < length:
            low = mid
        else:
            high = mid
    q = 0.5 * (low + high)
    # Renormalise so the widths sum to the length exactly, absorbing the
    # bisection residual into the first cell rather than into the last edge.
    first = length * (q - 1.0) / (q ** n - 1.0) if abs(q - 1.0) > 1e-14 else length / n
    return Segment(float(lo), float(hi), n, q ** (n - 1), first, first * q ** (n - 1))


# --------------------------------------------------------------------------- #
# the plan
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class MeshPlan:
    """Everything about the mesh, computed before a line of it is written.

    Attributes:
        x_segments, y_segments, z_segments: the block lattice.
        core_cell_size_m: the uniform cell size in the core [m].
        requested_cell_size_m: what the DSMC criterion asked for [m], before any
            budget coarsening.
        coarsened: whether ``mesh.max_cells`` forced ``core_cell_size_m`` up.
        mean_free_path_m: ``lambda0`` at the exit, for the cell/mfp ratio.
        max_cells: the configured budget.
    """

    x_segments: tuple
    y_segments: tuple
    z_segments: tuple
    core_cell_size_m: float
    requested_cell_size_m: float
    coarsened: bool
    mean_free_path_m: float
    max_cells: int

    # --- counts and sizes ----------------------------------------------------
    @property
    def divisions(self) -> tuple:
        """``(nx, ny, nz)`` summed over the segments of each axis."""
        return (sum(s.n_cells for s in self.x_segments),
                sum(s.n_cells for s in self.y_segments),
                sum(s.n_cells for s in self.z_segments))

    @property
    def n_cells(self) -> int:
        nx, ny, nz = self.divisions
        return int(nx) * int(ny) * int(nz)

    @property
    def n_blocks(self) -> int:
        return len(self.x_segments) * len(self.y_segments) * len(self.z_segments)

    @property
    def min_cell_size_m(self) -> float:
        """The smallest cell edge in the mesh [m]. Sets the Courant number."""
        return min(s.min_cell_m for s in
                   self.x_segments + self.y_segments + self.z_segments)

    @property
    def max_cell_size_m(self) -> float:
        """The largest cell edge in the mesh [m], out at the vacuum boundary."""
        return max(s.max_cell_m for s in
                   self.x_segments + self.y_segments + self.z_segments)

    @property
    def core_cell_sizes_m(self) -> tuple:
        """``(dx, dy, dz)`` of an actual core cell [m].

        Not ``core_cell_size_m`` three times. That is the size the criterion
        *asked* for; each axis then rounds its cell count **up** to fit the core
        extent, so the realised cells are a percent or two smaller and differ
        between axes. The particle weight and the occupancy are proportional to
        the cell volume, so using the requested size instead of the realised one
        would put a few percent of error into every derived quantity for no
        reason.
        """
        def core_of(segments) -> float:
            uniform = [s for s in segments if s.expansion == 1.0]
            return (uniform[0] if uniform else segments[0]).first_cell_m

        return (self.x_segments[0].first_cell_m,
                core_of(self.y_segments), core_of(self.z_segments))

    @property
    def exit_cell_volume_m3(self) -> float:
        """Volume of an actual cell at the exit [m^3]."""
        dx, dy, dz = self.core_cell_sizes_m
        return dx * dy * dz

    @property
    def cell_over_mfp(self) -> float:
        """Core cell size over the exit mean free path. Cai asks for 1."""
        return self.core_cell_size_m / self.mean_free_path_m

    @property
    def max_grading_jump(self) -> float:
        """Worst ratio of an expanding segment's core-side cell to the core cell.

        1.0 when every expansion was solved for continuity. Above 1 means a cap
        bound somewhere and the mesh steps at the core boundary.

        The core-side cell is the segment's ``min_cell_m``, not its
        ``first_cell_m``: the low-side lateral segments expand towards
        **decreasing** ``y`` or ``z``, so their small cell is the last one.
        """
        outer = [s for s in self.x_segments + self.y_segments + self.z_segments
                 if s.expansion != 1.0]
        if not outer:
            return 1.0
        return max(s.min_cell_m for s in outer) / self.core_cell_size_m

    def as_dict(self) -> dict:
        nx, ny, nz = self.divisions
        return {
            "n_cells": self.n_cells,
            "n_blocks": self.n_blocks,
            "divisions": [int(nx), int(ny), int(nz)],
            "core_cell_size_m": self.core_cell_size_m,
            "requested_cell_size_m": self.requested_cell_size_m,
            "coarsened_for_budget": bool(self.coarsened),
            "min_cell_size_m": self.min_cell_size_m,
            "max_cell_size_m": self.max_cell_size_m,
            "cell_over_mean_free_path": self.cell_over_mfp,
            "max_grading_jump": self.max_grading_jump,
            "max_cells_budget": int(self.max_cells),
        }

    def describe(self) -> list[str]:
        nx, ny, nz = self.divisions
        lines = [
            "Mesh plan (graded Cartesian, full 3-D)",
            f"  blocks                   {self.n_blocks}"
            f"  ({len(self.x_segments)} x {len(self.y_segments)} x "
            f"{len(self.z_segments)} segments)",
            f"  cells                    {nx} x {ny} x {nz} = {self.n_cells:,}"
            f"   (budget {self.max_cells:,})",
            f"  core cell                {fmt(self.core_cell_size_m)} m",
            f"  smallest / largest cell  {fmt(self.min_cell_size_m)} m  /  "
            f"{fmt(self.max_cell_size_m)} m",
            f"  cell / lambda0           {self.cell_over_mfp:.3f}"
            f"   (Cai asks for 1)",
        ]
        if self.coarsened:
            lines.append(
                f"  COARSENED: the criterion asked for "
                f"{fmt(self.requested_cell_size_m)} m, which would not fit in "
                f"{self.max_cells:,} cells. Raise mesh.max_cells to honour it.")
        if self.max_grading_jump > 1.001:
            lines.append(
                f"  grading jump             {self.max_grading_jump:.2f}x at the "
                f"core boundary -- mesh.outer_expansion capped the ratio; raise "
                f"it or add outer cells for a continuous mesh")
        for axis, segments in (("x", self.x_segments), ("y", self.y_segments),
                               ("z", self.z_segments)):
            for s in segments:
                kind = "uniform" if s.expansion == 1.0 else f"grading {s.expansion:.3g}"
                lines.append(
                    f"    {axis}  {fmt(s.lo):>12} .. {fmt(s.hi):<12} "
                    f"{s.n_cells:>4} cells  {kind:<16} "
                    f"{fmt(s.first_cell_m)} -> {fmt(s.last_cell_m)} m")
        return lines


def plan(cfg, geom: CaiGeometry, exit_state) -> MeshPlan:
    """Work out the whole mesh from the config, the geometry and the exit state.

    Args:
        cfg: the case config.
        geom: the geometry.
        exit_state: the :class:`~plumetools.cai2012.inflow.ExitState`, for
            ``lambda0``.

    Returns:
        A :class:`MeshPlan`.

    Raises:
        ValueError: if the cell budget cannot be met even at one cell per core
            extent, which means the domain or the budget is nonsense rather than
            the resolution being ambitious.

    The core cell size is the smaller of the collision criterion
    ``mesh.target_cell_over_mfp * lambda0`` and the geometric cap
    ``mesh.max_core_cell_over_D * D``, unless ``mesh.core_cell_size_m`` pins it.
    The cap is what makes ``Kn = 100`` meshable at all: there
    ``lambda0 = 20 m``, a hundred nozzle diameters, and the collision criterion
    alone would allow a cell larger than the nozzle. If the result would exceed
    ``mesh.max_cells``
    the cell is coarsened -- by a factor computed from the cube-root of the
    overshoot, then bisected up until it fits -- and the fact is recorded on the
    plan so every report says so. §5 of the case specification asks for exactly
    this: do not blindly generate a multi-billion-cell mesh.
    """
    d = geom.diameter_m
    lam = float(exit_state.mean_free_path_m)

    if cfg.mesh.core_cell_size_m is not None:
        requested = float(cfg.mesh.core_cell_size_m)
    else:
        requested = min(float(cfg.mesh.target_cell_over_mfp) * lam,
                        float(cfg.mesh.max_core_cell_over_D) * d)
    if cfg.mesh.min_core_cell_size_m is not None:
        requested = max(requested, float(cfg.mesh.min_core_cell_size_m))

    # Never so coarse the exit disk is one cell across, whatever was pinned.
    requested = min(requested, geom.radius_m)

    core_x = float(cfg.mesh.core_x_over_D) * d
    core_half = float(cfg.mesh.core_half_over_D) * d

    def build(cell: float):
        x_segments = [uniform_segment(0.0, core_x, cell)]
        x_segments.append(geometric_segment(
            core_x, geom.x_max_m, int(cfg.mesh.outer_x_cells),
            x_segments[0].last_cell_m, cfg.mesh.outer_expansion))

        def lateral(lo: float, hi: float) -> list:
            """Segments along a lateral axis: [-outer, core, +outer]."""
            segments = []
            if lo < -core_half:
                # Expansion runs outward, i.e. towards -y. In the block's local
                # +y direction the cells therefore CONTRACT, so the grading ratio
                # is the reciprocal and the segment is built from the core end.
                mirrored = geometric_segment(
                    lo, -core_half, int(cfg.mesh.outer_lateral_cells),
                    cell, cfg.mesh.outer_expansion)
                segments.append(Segment(
                    lo, -core_half, mirrored.n_cells, 1.0 / mirrored.expansion,
                    mirrored.last_cell_m, mirrored.first_cell_m))
            inner_lo = max(lo, -core_half)
            segments.append(uniform_segment(inner_lo, min(hi, core_half), cell))
            if hi > core_half:
                segments.append(geometric_segment(
                    core_half, hi, int(cfg.mesh.outer_lateral_cells),
                    segments[-1].last_cell_m, cfg.mesh.outer_expansion))
            return segments

        y_segments = lateral(geom.y_min_m, geom.y_max_m)
        z_segments = lateral(geom.z_min_m, geom.z_max_m)
        return x_segments, y_segments, z_segments

    def count(cell: float) -> int:
        xs, ys, zs = build(cell)
        return (sum(s.n_cells for s in xs) * sum(s.n_cells for s in ys)
                * sum(s.n_cells for s in zs))

    cell = requested
    coarsened = False
    budget = int(cfg.mesh.max_cells)
    if count(cell) > budget:
        coarsened = True
        # One cube-root jump gets very close; the loop then walks the rest.
        cell = requested * (count(requested) / budget) ** (1.0 / 3.0)
        while count(cell) > budget:
            cell *= 1.05
            if cell > geom.radius_m:
                # Stop here rather than coarsening to absurdity: past this point
                # the exit disk is under two cells across and the mesh has
                # stopped being a discretisation of this problem. The outer
                # segments also have FIXED cell counts, so a budget below their
                # product can never be met at any cell size -- both cases land
                # on the same, more useful message.
                raise ValueError(
                    f"meeting mesh.max_cells = {budget:,} needs a core cell "
                    f"coarser than the nozzle radius ({fmt(geom.radius_m)} m), "
                    f"so the exit disk would not be resolved at all. The "
                    f"smallest mesh this configuration can produce is "
                    f"{count(geom.radius_m):,} cells "
                    f"(mesh.outer_x_cells and mesh.outer_lateral_cells are "
                    f"fixed counts). Raise mesh.max_cells or shrink the domain.")

    x_segments, y_segments, z_segments = build(cell)
    return MeshPlan(
        x_segments=tuple(x_segments),
        y_segments=tuple(y_segments),
        z_segments=tuple(z_segments),
        core_cell_size_m=cell,
        requested_cell_size_m=requested,
        coarsened=coarsened,
        mean_free_path_m=lam,
        max_cells=budget,
    )


# --------------------------------------------------------------------------- #
# cost of the literal Cai translation
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class UniformCost:
    """What a literal 3-D translation of Cai's uniform grid would cost."""

    cell_size_m: float
    divisions: tuple
    n_cells: float
    particles_at_target: float
    target_particles_per_cell: float

    def describe(self) -> list[str]:
        nx, ny, nz = self.divisions
        return [
            "Literal 3-D translation of Cai's uniform grid (NOT built)",
            f"  dx = dy = dz             {fmt(self.cell_size_m)} m"
            f"   (= lambda0 at Kn = 0.01)",
            f"  cells                    {nx:,} x {ny:,} x {nz:,} = "
            f"{self.n_cells:.3e}",
            f"  parcels at {self.target_particles_per_cell:g}/cell     "
            f"{self.particles_at_target:.3e}",
            "  Cai's grid is affordable because his solver is AXISYMMETRIC: it is",
            "  a 2-D (x, r) sheet. In 3-D Cartesian the same criterion is a",
            "  multi-billion-cell mesh, so this case grades instead -- see",
            "  docs/cai2012-case.md.",
        ]


def uniform_cost(cfg, geom: CaiGeometry, reference_cell_size_m: float) -> UniformCost:
    """Cell and parcel count of the uniform mesh Cai's criterion asks for.

    Args:
        cfg: the case config, for the occupancy target.
        geom: the geometry.
        reference_cell_size_m: Cai's ``dx``, i.e.
            ``target_cell_over_mfp * lambda_ref`` at ``Kn = 0.01``.

    Returns:
        A :class:`UniformCost`. Nothing is written; this exists to be printed.

    The parcel figure is deliberately the *optimistic* one -- the target
    occupancy everywhere, when in reality the outer domain would hold far fewer
    and the core far more. It is a lower bound on an already impossible number.
    """
    cell = float(reference_cell_size_m)
    if not cell > 0.0:
        raise ValueError(f"reference cell size must be positive, got {cell}")
    nx = math.ceil((geom.x_max_m - geom.x_min_m) / cell)
    ny = math.ceil((geom.y_max_m - geom.y_min_m) / cell)
    nz = math.ceil((geom.z_max_m - geom.z_min_m) / cell)
    n_cells = float(nx) * float(ny) * float(nz)
    target = float(cfg.resolution.target_particles_per_cell)
    return UniformCost(
        cell_size_m=cell,
        divisions=(nx, ny, nz),
        n_cells=n_cells,
        particles_at_target=n_cells * target,
        target_particles_per_cell=target,
    )


# --------------------------------------------------------------------------- #
# the nozzle patch
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class NozzlePrediction:
    """The staircase nozzle patch the mesh will produce."""

    n_faces: int
    area_m2: float
    exact_area_m2: float
    cells_across_diameter: float

    @property
    def area_error(self) -> float:
        """Signed relative area error, ``(meshed - exact) / exact``."""
        return (self.area_m2 - self.exact_area_m2) / self.exact_area_m2

    def describe(self) -> list[str]:
        return [
            "Nozzle patch (a staircase on a Cartesian grid)",
            f"  faces                    {self.n_faces}",
            f"  cells across D           {self.cells_across_diameter:.1f}",
            f"  area                     {self.area_m2:.6e} m^2  against "
            f"pi R0^2 = {self.exact_area_m2:.6e} m^2",
            f"  area error               {100.0 * self.area_error:+.2f}%",
        ]


def predicted_nozzle_faces(plan_: MeshPlan, geom: CaiGeometry) -> NozzlePrediction:
    """Which ``x = 0`` faces will land in the nozzle patch, and their area.

    Applies the same rule ``topoSet``'s ``cylinderToFace`` applies -- face centre
    inside ``r <= R0`` -- to the ``(y, z)`` cell centres of the block lattice.
    The mesh has not been built yet, so this is a *prediction*;
    :func:`plumetools.cai2012.checks.check_nozzle_area` measures the built mesh
    and compares. When the two disagree, the mesh is not the mesh that was
    planned, which is worth finding out before a solver runs on it.
    """
    y_edges = np.concatenate([s.cell_edges()[:-1] for s in plan_.y_segments]
                             + [plan_.y_segments[-1].cell_edges()[-1:]])
    z_edges = np.concatenate([s.cell_edges()[:-1] for s in plan_.z_segments]
                             + [plan_.z_segments[-1].cell_edges()[-1:]])
    yc = 0.5 * (y_edges[:-1] + y_edges[1:])
    zc = 0.5 * (z_edges[:-1] + z_edges[1:])
    dy = np.diff(y_edges)
    dz = np.diff(z_edges)

    yy, zz = np.meshgrid(yc, zc, indexing="ij")
    areas = np.outer(dy, dz)
    inside = geom.is_on_nozzle(yy, zz)

    _, dy, dz = plan_.core_cell_sizes_m
    return NozzlePrediction(
        n_faces=int(np.count_nonzero(inside)),
        area_m2=float(areas[inside].sum()),
        exact_area_m2=geom.exact_nozzle_area_m2,
        cells_across_diameter=geom.diameter_m / max(dy, dz),
    )


# --------------------------------------------------------------------------- #
# particle-count estimate
# --------------------------------------------------------------------------- #

def estimate_particles(plan_: MeshPlan, geom: CaiGeometry, exit_state,
                       n_equivalent_particles: float, *,
                       samples_per_axis: int = 20) -> dict:
    """Estimate the steady-state parcel count from the analytical solution.

    Args:
        plan_: the mesh plan (for the report only).
        geom: the geometry.
        exit_state: the exit state.
        n_equivalent_particles: the particle weight.
        samples_per_axis: coarse grid used to integrate the analytical density
            over the domain.

    Returns:
        ``{"total_molecules", "total_particles", "mean_particles_per_cell"}``.

    The **collisionless** density is integrated over the box, which is exact for
    ``Kn = 100`` and an underestimate for the denser cases, where collisions push
    molecules out of the plume core into the surrounding vacuum without changing
    the total. Cheap and honest to within a factor that matters far less than
    knowing whether the answer is 1e6 or 1e11 parcels before submitting the run.
    """
    from plumetools.cai2012 import analytical

    n = max(4, int(samples_per_axis))
    x_edges = np.linspace(geom.x_min_m, geom.x_max_m, n + 1)
    y_edges = np.linspace(geom.y_min_m, geom.y_max_m, n + 1)
    z_edges = np.linspace(geom.z_min_m, geom.z_max_m, n + 1)
    xc = 0.5 * (x_edges[:-1] + x_edges[1:])
    yc = 0.5 * (y_edges[:-1] + y_edges[1:])
    zc = 0.5 * (z_edges[:-1] + z_edges[1:])
    cell_volume = ((x_edges[1] - x_edges[0]) * (y_edges[1] - y_edges[0])
                   * (z_edges[1] - z_edges[0]))

    xx, yy, zz = np.meshgrid(xc, yc, zc, indexing="ij")
    ratio = analytical.density_ratio(
        xx, yy, zz, geom.radius_m, exit_state.speed_ratio,
        n_radial=16, n_azimuthal=32)

    molecules = float(ratio.sum()) * cell_volume * exit_state.number_density_per_m3
    particles = molecules / float(n_equivalent_particles)
    return {
        "total_molecules": molecules,
        "total_particles": particles,
        "mean_particles_per_cell": particles / max(1, plan_.n_cells),
    }


# --------------------------------------------------------------------------- #
# dictionaries
# --------------------------------------------------------------------------- #

def _lattice(plan_: MeshPlan):
    """Vertex coordinates of the block lattice, and an index function."""
    xs = [plan_.x_segments[0].lo] + [s.hi for s in plan_.x_segments]
    ys = [plan_.y_segments[0].lo] + [s.hi for s in plan_.y_segments]
    zs = [plan_.z_segments[0].lo] + [s.hi for s in plan_.z_segments]

    def index(i: int, j: int, k: int) -> int:
        return i + j * len(xs) + k * len(xs) * len(ys)

    return xs, ys, zs, index


def render_block_mesh_dict(cfg, geom: CaiGeometry, plan_: MeshPlan) -> str:
    """The graded multi-block ``blockMeshDict``.

    ``x = 0`` is written as a **single** patch, ``upstreamVacuum``. The nozzle is
    carved out of it afterwards by ``topoSet`` + ``createPatch``, because a disk
    is not a block face and blockMesh can only name block faces.
    """
    xs, ys, zs, index = _lattice(plan_)
    names = cfg.mesh.patch_names
    nx, ny, nz = plan_.divisions

    lines = header(
        "blockMeshDict", GENERATOR,
        f"Cai 2012: {plan_.n_blocks} blocks, {nx}x{ny}x{nz} = {plan_.n_cells} "
        f"cells, core {fmt(plan_.core_cell_size_m)} m.")
    lines += [
        "// The exit plane x = 0 is ONE patch here. topoSet + createPatch carve",
        "// the circular 'nozzle' out of it afterwards: a disk is not a block",
        "// face, and blockMesh can only name block faces.",
        "//",
        "// Every outer boundary is geometric type 'patch', never 'wall'.",
        "// particle::hitBoundaryFace finds no handler for a plain patch and sets",
        "// keepParticle = false -- the molecule is DELETED, which is what",
        "// expanding into vacuum means. This is only possible because the case",
        "// runs plumeFieldInflow, which injects across a named patch list; stock",
        "// FreeStream injects on every patch-type boundary and would force these",
        "// to be walls, i.e. a closed reflecting box.",
        "",
        "scale   1;",
        "",
        "vertices",
        "(",
    ]
    for k, z in enumerate(zs):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                lines.append(f"    {fmt_point((x, y, z))}   // {index(i, j, k)}")
    lines += [");", "", "blocks", "("]

    for k, zseg in enumerate(plan_.z_segments):
        for j, yseg in enumerate(plan_.y_segments):
            for i, xseg in enumerate(plan_.x_segments):
                corners = (
                    index(i, j, k), index(i + 1, j, k),
                    index(i + 1, j + 1, k), index(i, j + 1, k),
                    index(i, j, k + 1), index(i + 1, j, k + 1),
                    index(i + 1, j + 1, k + 1), index(i, j + 1, k + 1),
                )
                lines.append(
                    "    hex (" + " ".join(str(c) for c in corners) + ") "
                    f"({xseg.n_cells} {yseg.n_cells} {zseg.n_cells}) "
                    f"simpleGrading ({fmt(xseg.expansion)} "
                    f"{fmt(yseg.expansion)} {fmt(zseg.expansion)})")
    lines += [
        ");",
        "",
        "// No curved edges: the only curved feature is the nozzle rim, and that",
        "// is a face SELECTION, not a geometric edge.",
        "edges",
        "(",
        ");",
        "",
        "boundary",
        "(",
    ]

    # Collect the exterior faces of the lattice, per patch role.
    faces: dict[str, list] = {}
    comments = {}

    def add(role: str, quad, comment: str) -> None:
        name = names[role]
        faces.setdefault(name, []).append(quad)
        comments.setdefault(name, comment)

    ni, nj, nk = (len(plan_.x_segments), len(plan_.y_segments),
                  len(plan_.z_segments))
    symmetry_role = ("symmetry" if geom.symmetry_mode == "half_y"
                     and "symmetry" in names else None)

    for k in range(nk):
        for j in range(nj):
            for i in range(ni):
                local = {
                    0: index(i, j, k), 1: index(i + 1, j, k),
                    2: index(i + 1, j + 1, k), 3: index(i, j + 1, k),
                    4: index(i, j, k + 1), 5: index(i + 1, j, k + 1),
                    6: index(i + 1, j + 1, k + 1), 7: index(i, j + 1, k + 1),
                }

                def quad(face_key: str):
                    return tuple(local[v] for v in BOX_FACES[face_key])

                if i == 0:
                    add("upstream_vacuum", quad("x_min"),
                        "x = 0, outside the nozzle disk. An OPEN boundary: a "
                        "molecule that scatters back through the exit plane has "
                        "left the domain. NOT a symmetry plane.")
                if i == ni - 1:
                    add("outer", quad("x_max"), "the vacuum far field")
                if j == 0:
                    if symmetry_role:
                        add(symmetry_role, quad("y_min"),
                            "y = 0, the modelled plane of symmetry")
                    else:
                        add("outer", quad("y_min"), "the vacuum far field")
                if j == nj - 1:
                    add("outer", quad("y_max"), "the vacuum far field")
                if k == 0:
                    add("outer", quad("z_min"), "the vacuum far field")
                if k == nk - 1:
                    add("outer", quad("z_max"), "the vacuum far field")

    order = [names["upstream_vacuum"]]
    if symmetry_role:
        order.append(names[symmetry_role])
    order.append(names["outer"])

    for name in order:
        patch_type = "symmetry" if (symmetry_role and name == names[symmetry_role]) \
            else "patch"
        lines += [f"    {name}", "    {", f"        type {patch_type};",
                  f"        // {comments[name]}", "        faces", "        ("]
        lines += [f"            ({a} {b} {c} {d})" for a, b, c, d in faces[name]]
        lines += ["        );", "    }", ""]

    lines += [");", "", "mergePatchPairs", "(", ");", "",
              "// ************************************************************************* //"]
    return "\n".join(lines) + "\n"


def render_topo_set_dict(cfg, geom: CaiGeometry, plan_: MeshPlan) -> str:
    """``system/topoSetDict``: the ``x = 0`` faces inside ``r <= R0``.

    Two actions, and the order is load-bearing. ``patchToFace`` first restricts
    the set to the ``x = 0`` boundary; only then does ``cylinderToFace`` subset it
    by radius. Run the other way round, ``cylinderToFace`` would also select every
    **internal** face inside the cylinder, and ``createPatch`` would be asked to
    turn interior faces into a boundary.
    """
    names = cfg.mesh.patch_names
    half = SELECTOR_HALF_THICKNESS_FRACTION * plan_.min_cell_size_m

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


def render_create_patch_dict(cfg, geom: CaiGeometry) -> str:
    """``system/createPatchDict``: move ``nozzleFaces`` into the nozzle patch.

    ``type patch``, never ``wall``. ``plumeFieldInflow`` refuses outright to
    inject across a ``wallPolyPatch`` -- injecting through a wall would create
    molecules inside a solid and suppress the surface fluxes ``hitWallPatch``
    records -- so a mis-typed nozzle is a hard error rather than a silent
    physical change.
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


def render_mesh_quality_dict() -> str:
    """``system/meshQualityDict``. Required by ``createPatch``'s mesh handling."""
    lines = header("meshQualityDict", GENERATOR,
                   "OpenFOAM's own defaults; a blockMesh box passes all of them.")
    lines += [
        '#includeEtc "caseDicts/meshQualityDict"',
        "",
        "// ************************************************************************* //",
    ]
    return "\n".join(lines) + "\n"


def write_mesh_setup(case_dir: Path, cfg, geom: CaiGeometry,
                     plan_: MeshPlan) -> list[Path]:
    """Write ``blockMeshDict``, ``topoSetDict``, ``createPatchDict``, quality.

    Returns:
        The paths written, in the order the mesh pipeline consumes them.
    """
    system = Path(case_dir) / "system"
    return [
        write_lines(system / "blockMeshDict",
                    render_block_mesh_dict(cfg, geom, plan_).splitlines()),
        write_lines(system / "topoSetDict",
                    render_topo_set_dict(cfg, geom, plan_).splitlines()),
        write_lines(system / "createPatchDict",
                    render_create_patch_dict(cfg, geom).splitlines()),
        write_lines(system / "meshQualityDict",
                    render_mesh_quality_dict().splitlines()),
    ]
