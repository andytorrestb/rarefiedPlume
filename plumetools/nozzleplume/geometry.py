r"""The geometry: an exit disk at ``x = 0``, a vacuum box, and perhaps a plate.

::

                          vacuum  (molecules deleted)
      y,z  ^
           |   +------------------------------------------+
           |   |                              ##          |
           |   |          expanding plume     ##  plate   |
    -------+---#====#-------------------------##--------->+  x
           |   |    ^                         ##          |
           |   |    nozzle: the disk          ##          |
           |   |    r <= R0 at x = 0                      |
           |   +------------------------------------------+
               x = 0                    standoff       x_max

               <---- upstreamVacuum: the REST of the x = 0 plane

Four facts about this layout are load-bearing. The first three are
:mod:`plumetools.cai2012.geometry`'s, unchanged; the fourth is new.

**The nozzle is the inlet.** The DSMC injection happens on the physical exit
disk, at the physical exit state.

**``upstreamVacuum`` is not a symmetry plane.** It is the part of the ``x = 0``
boundary outside the disk, and it is *open*: a molecule that comes back through
it -- and with a plate downstream, reflected molecules do -- has left the domain.
There is no spacecraft behind the nozzle in this model for it to hit.

**The outer boundaries are geometric type** ``patch``. ``particle::hitBoundaryFace``
finds no handler for a plain ``patch`` and deletes the molecule, which is what
"expands into vacuum" means.

**The plate is geometric type** ``wall``. ``DSMCParcel::hitWallPatch`` then runs:
it records the momentum and energy the gas delivers (``fD``, ``q``) and applies
the gas-surface model. A plate declared ``patch`` would silently absorb the
plume and feel no load at all.

The plate is a rectangular box **normal to the plume axis**: its thickness is
along ``x``, and it is positioned by the distance from the exit plane to its
upstream face. A tilted plate is not expressible here.

The box is built on the background lattice
------------------------------------------
Every domain extent is rounded **up** to a whole number of background cells.
That keeps the background mesh exactly cubic -- snappyHexMesh's octree halves a
cell in all three directions at once, so a cell that starts 3% off cubic stays
3% off at every level -- and it puts ``y = 0`` and ``z = 0`` on lattice planes,
which is what makes the staircase nozzle four-fold symmetric and a symmetry
plane a matter of deleting half the box.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


class GeometryError(ValueError):
    """The configured geometry is inconsistent or cannot be meshed."""


@dataclass(frozen=True)
class Plate:
    """The plate as a box in metres. Always the **whole** plate.

    Attributes:
        x_upstream_m: the face the plume strikes [m].
        thickness_m: extent along ``x`` [m].
        y_min_m, y_max_m, z_min_m, z_max_m: lateral extents [m].
        temperature_K: wall temperature [K].

    Under a symmetry plane only part of this box is in the domain;
    :meth:`PlumeGeometry.modelled_plate_bounds` is the part that is.
    """

    x_upstream_m: float
    thickness_m: float
    y_min_m: float
    y_max_m: float
    z_min_m: float
    z_max_m: float
    temperature_K: float = 300.0

    @property
    def x_downstream_m(self) -> float:
        """The lee face [m]."""
        return self.x_upstream_m + self.thickness_m

    @property
    def size_y_m(self) -> float:
        return self.y_max_m - self.y_min_m

    @property
    def size_z_m(self) -> float:
        return self.z_max_m - self.z_min_m

    @property
    def face_area_m2(self) -> float:
        """Area of the upstream face of the whole plate [m^2]."""
        return self.size_y_m * self.size_z_m

    @property
    def bounds(self) -> tuple:
        """``(x_min, x_max, y_min, y_max, z_min, z_max)`` [m]."""
        return (self.x_upstream_m, self.x_downstream_m, self.y_min_m,
                self.y_max_m, self.z_min_m, self.z_max_m)

    def contains(self, x, y, z, margin_m: float = 0.0):
        """Whether a point is inside the plate, grown by ``margin_m`` all round."""
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        z = np.asarray(z, dtype=np.float64)
        m = float(margin_m)
        return ((x >= self.x_upstream_m - m) & (x <= self.x_downstream_m + m)
                & (y >= self.y_min_m - m) & (y <= self.y_max_m + m)
                & (z >= self.z_min_m - m) & (z <= self.z_max_m + m))


@dataclass(frozen=True)
class PlumeGeometry:
    """The physical geometry, in metres, as it will be meshed.

    Attributes:
        diameter_m: nozzle exit diameter ``D`` [m].
        x_max_m: downstream boundary [m]. The exit plane is ``x = 0``.
        y_min_m, y_max_m, z_min_m, z_max_m: lateral extents of the **modelled**
            box [m]; a lower bound of 0 is a symmetry plane.
        symmetry_planes: the axes, of ``y`` and ``z``, with a symmetry plane.
        plate: the plate, or ``None`` for a free plume.
        requested_extents_m: ``(x_max, y_half, z_half)`` as ``case.yaml`` gave
            them, before rounding up to the background lattice.
    """

    diameter_m: float
    x_max_m: float
    y_min_m: float
    y_max_m: float
    z_min_m: float
    z_max_m: float
    symmetry_planes: tuple = ()
    plate: Plate | None = None
    requested_extents_m: tuple = ()

    # --- the exit disk -------------------------------------------------------
    @property
    def radius_m(self) -> float:
        """``R0 = D/2`` [m]."""
        return 0.5 * self.diameter_m

    @property
    def symmetry_multiplicity(self) -> int:
        """How many copies of the modelled domain make the whole: 1, 2 or 4.

        Every extensive quantity measured in the modelled domain -- a mass flow,
        a force on the plate, a parcel count -- is this fraction of the real
        one.
        """
        return 2 ** len(self.symmetry_planes)

    @property
    def full_nozzle_area_m2(self) -> float:
        """``pi R0^2`` [m^2], the whole exit."""
        return math.pi * self.radius_m ** 2

    @property
    def exact_nozzle_area_m2(self) -> float:
        """The part of the exit disk inside the modelled domain [m^2].

        The *exact* area. What a Cartesian mesh produces is a staircase
        approximation to it, predicted by
        :func:`plumetools.nozzleplume.mesh.predicted_nozzle_faces` and measured
        back out of the mesh by :mod:`plumetools.nozzleplume.verify`.
        """
        return self.full_nozzle_area_m2 / self.symmetry_multiplicity

    # --- the box -------------------------------------------------------------
    @property
    def x_min_m(self) -> float:
        """The exit plane. Always 0."""
        return 0.0

    @property
    def bounds(self) -> tuple:
        """``(x_min, x_max, y_min, y_max, z_min, z_max)`` [m]."""
        return (self.x_min_m, self.x_max_m, self.y_min_m, self.y_max_m,
                self.z_min_m, self.z_max_m)

    @property
    def domain_volume_m3(self) -> float:
        """Volume of the modelled box [m^3]."""
        return ((self.x_max_m - self.x_min_m) * (self.y_max_m - self.y_min_m)
                * (self.z_max_m - self.z_min_m))

    @property
    def modelled_plate_bounds(self) -> tuple | None:
        """The plate clipped to the modelled box, or ``None`` with no plate.

        Identical to :attr:`Plate.bounds` in a full 3-D domain. Under a symmetry
        plane the plate straddles it, and what is meshed is the half or quarter
        on the modelled side.
        """
        if self.plate is None:
            return None
        x0, x1, y0, y1, z0, z1 = self.plate.bounds
        return (x0, x1, max(y0, self.y_min_m), min(y1, self.y_max_m),
                max(z0, self.z_min_m), min(z1, self.z_max_m))

    @property
    def flow_length_m(self) -> float:
        """How far the jet travels before something stops it [m].

        The stand-off with a plate, the length of the box without one.
        """
        return self.plate.x_upstream_m if self.plate is not None else self.x_max_m

    # --- membership ----------------------------------------------------------
    def is_on_nozzle(self, y, z):
        """Whether ``(y, z)`` lies on the modelled part of the exit disk.

        The test a face centre is put through to decide whether it belongs to
        the ``nozzle`` patch. ``topoSet``'s ``cylinderToFace`` applies exactly
        this rule in OpenFOAM, so the Python answer and the meshed patch agree by
        construction -- which is what makes the predicted face count checkable
        against the real one.
        """
        y = np.asarray(y, dtype=np.float64)
        z = np.asarray(z, dtype=np.float64)
        return ((y * y + z * z <= self.radius_m ** 2)
                & (y >= self.y_min_m) & (z >= self.z_min_m))

    def contains(self, x, y, z):
        """Whether a point is inside the modelled box (boundary included)."""
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        z = np.asarray(z, dtype=np.float64)
        return ((x >= self.x_min_m) & (x <= self.x_max_m)
                & (y >= self.y_min_m) & (y <= self.y_max_m)
                & (z >= self.z_min_m) & (z <= self.z_max_m))

    # --- reporting -----------------------------------------------------------
    def describe(self) -> list[str]:
        """Report lines, printed by ``./Allmesh``."""
        d = self.diameter_m
        lines = [
            "Geometry",
            f"  nozzle exit              the disk r <= {self.radius_m:.6g} m at"
            f" x = 0; plume axis +x",
            f"  domain x                 0 .. {self.x_max_m:g} m"
            f"   = {self.x_max_m / d:.1f} D",
            f"  domain y                 {self.y_min_m:g} .. {self.y_max_m:g} m",
            f"  domain z                 {self.z_min_m:g} .. {self.z_max_m:g} m",
        ]
        if self.requested_extents_m:
            asked = tuple(float(v) for v in self.requested_extents_m)
            built = (self.x_max_m, self.y_max_m, self.z_max_m)
            if any(abs(a - b) > 1e-9 for a, b in zip(asked, built)):
                lines.append(
                    f"                           rounded up from x_max {asked[0]:g},"
                    f" y_half {asked[1]:g}, z_half {asked[2]:g} m to whole"
                    f" background cells")
        lines.append(f"  volume                   {self.domain_volume_m3:.6g} m^3")
        if self.symmetry_planes:
            planes = ", ".join(f"{axis} = 0" for axis in self.symmetry_planes)
            lines.append(
                f"  symmetry                 {planes}: 1/"
                f"{self.symmetry_multiplicity} of the flow is modelled")
        else:
            lines.append("  symmetry                 none; full 3-D")

        plate = self.plate
        if plate is None:
            lines.append("  plate                    none; a free plume into vacuum")
            return lines
        lines += [
            f"  plate                    {plate.size_y_m:g} (y) x {plate.size_z_m:g}"
            f" (z) m, {plate.thickness_m:g} m thick, normal to the plume axis",
            f"    upstream face          x = {plate.x_upstream_m:g} m"
            f"   = {plate.x_upstream_m / d:.1f} D from the exit",
            f"    spans                  y {plate.y_min_m:g} .. {plate.y_max_m:g} m,"
            f"  z {plate.z_min_m:g} .. {plate.z_max_m:g} m",
            f"    clear of the boundary  {self.x_max_m - plate.x_downstream_m:g} m"
            f" behind, {self.y_max_m - plate.y_max_m:g} m in y,"
            f" {self.z_max_m - plate.z_max_m:g} m in z",
            f"    wall temperature       {plate.temperature_K:g} K",
        ]
        return lines


def _round_up(value: float, step: float) -> float:
    """``value`` rounded up to a whole number of ``step``.

    The 1e-9 keeps a value that already is a whole number of steps from being
    pushed to the next one by floating-point division.
    """
    return max(1, math.ceil(float(value) / float(step) - 1.0e-9)) * float(step)


def from_config(cfg) -> PlumeGeometry:
    """Build and validate the geometry from a loaded ``case.yaml``.

    Raises:
        GeometryError: if the exit disk does not fit in the box, or the plate is
            not strictly inside it. A plate flush with an open boundary would be
            meshed as though it continued past it, which removes the edge the
            gas escapes around -- so it must clear every outer boundary by at
            least one background cell.
    """
    cell = float(cfg.mesh.background_cell_size_m)
    domain = cfg.domain
    planes = tuple(domain.symmetry_planes)

    x_max = _round_up(domain.x_max_m, cell)
    y_half = _round_up(domain.y_half_m, cell)
    z_half = _round_up(domain.z_half_m, cell)

    plate = None
    if cfg.plate is not None:
        p = cfg.plate
        plate = Plate(
            x_upstream_m=float(p.standoff_m),
            thickness_m=float(p.thickness_m),
            y_min_m=float(p.centre_y_m) - 0.5 * float(p.size_y_m),
            y_max_m=float(p.centre_y_m) + 0.5 * float(p.size_y_m),
            z_min_m=float(p.centre_z_m) - 0.5 * float(p.size_z_m),
            z_max_m=float(p.centre_z_m) + 0.5 * float(p.size_z_m),
            temperature_K=float(p.temperature_K),
        )

    geom = PlumeGeometry(
        diameter_m=float(cfg.nozzle.exit_diameter_m),
        x_max_m=x_max,
        y_min_m=0.0 if "y" in planes else -y_half,
        y_max_m=y_half,
        z_min_m=0.0 if "z" in planes else -z_half,
        z_max_m=z_half,
        symmetry_planes=planes,
        plate=plate,
        requested_extents_m=(float(domain.x_max_m), float(domain.y_half_m),
                             float(domain.z_half_m)),
    )

    for axis, half in (("y", y_half), ("z", z_half)):
        if not half > geom.radius_m:
            raise GeometryError(
                f"domain.{axis}_half_m is {half:g} m, which does not contain the "
                f"exit disk of radius {geom.radius_m:g} m")

    if plate is not None:
        slack = cell * (1.0 - 1.0e-9)
        too_close = []
        if geom.x_max_m - plate.x_downstream_m < slack:
            too_close.append(
                f"its lee face at x = {plate.x_downstream_m:g} m is "
                f"{geom.x_max_m - plate.x_downstream_m:g} m from x_max = "
                f"{geom.x_max_m:g} m (raise domain.x_max_m)")
        for axis in ("y", "z"):
            lo = getattr(plate, f"{axis}_min_m")
            hi = getattr(plate, f"{axis}_max_m")
            half = getattr(geom, f"{axis}_max_m")
            if half - hi < slack:
                too_close.append(
                    f"its {axis} edge at {hi:g} m is {half - hi:g} m from the "
                    f"boundary at {half:g} m (raise domain.{axis}_half_m)")
            # The low side is a symmetry plane under symmetry, which the plate
            # is meant to cross.
            if axis not in planes and lo + half < slack:
                too_close.append(
                    f"its {axis} edge at {lo:g} m is {lo + half:g} m from the "
                    f"boundary at {-half:g} m (raise domain.{axis}_half_m)")
        if too_close:
            raise GeometryError(
                f"the plate must clear every open boundary by at least one "
                f"background cell ({cell:g} m), or it is meshed as though it "
                f"ran on past the boundary: " + "; ".join(too_close))

    return geom
