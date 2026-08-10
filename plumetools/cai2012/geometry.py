r"""The Cai 2012 geometry: a circular exit disk at ``x = 0`` and a vacuum box.

::

                              vacuum  (particles deleted)
      y,z  ^ +y_half
           |   +--------------------------------------------+
           |   |                                            |
           |   |            expanding plume                 |
    -------+---#====#-------------------------------------->+  x
           |   |    ^                                       |
           |   |    nozzle: the disk r <= R0 at x = 0        |
           |   +--------------------------------------------+
           v -y_half
               x=0                                     x = x_max

               <---- upstreamVacuum: the REST of the x = 0 plane

Three facts about this layout are load-bearing.

**The nozzle is the inlet.** Not an analytical source surface standing in for
one: the DSMC injection happens on the physical exit disk, at the physical exit
state. There is no hemisphere anywhere in this case family.

**``upstreamVacuum`` is not a symmetry plane.** It is the part of the ``x = 0``
boundary outside the disk, and it is an *open* boundary. A molecule that scatters
back through it has left the domain. Calling it a symmetry plane would reflect
those molecules back in and quietly turn the half-space into a mirror box.

**The outer boundaries must be geometric type** ``patch``. In OpenFOAM
``particle::hitBoundaryFace`` finds no handler for a plain ``patch`` and sets
``keepParticle = false`` -- the molecule is deleted, which is what "expands into
vacuum" means. A ``wall`` there would reflect. This is only possible because the
inflow model takes an explicit patch list; stock ``FreeStream`` would inject on
every ``patch`` in the mesh and force the walls.

Coordinates are normalised by the **diameter** throughout -- ``X/D``, ``Z/D`` --
because that is what Cai plots.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CaiGeometry:
    """The physical geometry, in metres, derived from ``case.yaml``.

    Attributes:
        diameter_m: nozzle diameter ``D`` [m].
        x_max_m: downstream boundary [m]. The exit plane is ``x = 0``.
        y_min_m, y_max_m, z_min_m, z_max_m: lateral extents [m].
        symmetry_mode: ``none`` or ``half_y``; see
            :data:`plumetools.cai2012.config.SYMMETRY_MODES`.
    """

    diameter_m: float
    x_max_m: float
    y_min_m: float
    y_max_m: float
    z_min_m: float
    z_max_m: float
    symmetry_mode: str = "none"

    # --- the exit disk -------------------------------------------------------
    @property
    def radius_m(self) -> float:
        """``R0 = D/2`` [m]."""
        return 0.5 * self.diameter_m

    @property
    def exact_nozzle_area_m2(self) -> float:
        """``pi R0^2`` [m^2] -- or half of it under ``half_y``.

        The *exact* disk area. What a Cartesian mesh actually produces is a
        staircase approximation to it, measured out of the mesh by
        :func:`plumetools.cai2012.checks.check_nozzle_area` and compared with
        this. They are not the same number and the difference is reported.
        """
        area = math.pi * self.radius_m ** 2
        return 0.5 * area if self.symmetry_mode == "half_y" else area

    @property
    def x_min_m(self) -> float:
        """The exit plane. Always 0: it is the origin of ``X/D``."""
        return 0.0

    @property
    def bounds(self) -> tuple:
        """``(x_min, x_max, y_min, y_max, z_min, z_max)`` [m]."""
        return (self.x_min_m, self.x_max_m, self.y_min_m, self.y_max_m,
                self.z_min_m, self.z_max_m)

    @property
    def domain_volume_m3(self) -> float:
        """Volume of the box [m^3]."""
        return ((self.x_max_m - self.x_min_m)
                * (self.y_max_m - self.y_min_m)
                * (self.z_max_m - self.z_min_m))

    # --- normalisation -------------------------------------------------------
    def over_diameter(self, length_m):
        """A length in nozzle diameters. Cai's ``X/D``, ``Z/D``."""
        return np.asarray(length_m, dtype=np.float64) / self.diameter_m

    def from_diameter(self, length_over_d):
        """The inverse: nozzle diameters back to metres."""
        return np.asarray(length_over_d, dtype=np.float64) * self.diameter_m

    # --- membership ----------------------------------------------------------
    def is_on_nozzle(self, y, z):
        """Whether ``(y, z)`` lies on the exit disk, ``y^2 + z^2 <= R0^2``.

        The test a face centre is put through to decide whether it belongs to the
        ``nozzle`` patch. ``topoSet``'s ``cylinderToFace`` applies exactly this
        rule in OpenFOAM, so the Python answer and the meshed patch agree by
        construction rather than by coincidence -- which is what makes the
        predicted face count checkable against the real one.
        """
        y = np.asarray(y, dtype=np.float64)
        z = np.asarray(z, dtype=np.float64)
        inside = y * y + z * z <= self.radius_m ** 2
        if self.symmetry_mode == "half_y":
            inside = inside & (y >= self.y_min_m)
        return inside

    def contains(self, x, y, z):
        """Whether a point is inside the box (inclusive of the boundary)."""
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        z = np.asarray(z, dtype=np.float64)
        return ((x >= self.x_min_m) & (x <= self.x_max_m)
                & (y >= self.y_min_m) & (y <= self.y_max_m)
                & (z >= self.z_min_m) & (z <= self.z_max_m))

    def describe(self) -> list[str]:
        """Report lines, printed by ``./Allmesh``."""
        d = self.diameter_m
        return [
            "Cai 2012 geometry",
            f"  nozzle diameter D        {d:g} m   (R0 = {self.radius_m:g} m)",
            f"  exact exit area          {self.exact_nozzle_area_m2:.6e} m^2"
            f"  = pi R0^2",
            f"  domain x                 0 .. {self.x_max_m:g} m"
            f"   = 0 .. {self.x_max_m / d:g} D",
            f"  domain y                 {self.y_min_m:g} .. {self.y_max_m:g} m"
            f"   = {self.y_min_m / d:g} .. {self.y_max_m / d:g} D",
            f"  domain z                 {self.z_min_m:g} .. {self.z_max_m:g} m"
            f"   = {self.z_min_m / d:g} .. {self.z_max_m / d:g} D",
            f"  volume                   {self.domain_volume_m3:.6g} m^3",
            f"  symmetry                 {self.symmetry_mode}",
        ]


def from_config(cfg) -> CaiGeometry:
    """Build the geometry from a :class:`~plumetools.cai2012.config.CaiCaseConfig`.

    Every extent in ``case.yaml`` is expressed in nozzle diameters, so the
    dimensional geometry exists only here. Changing ``D`` rescales the whole
    domain, which is the intent: the validation quantities are normalised.
    """
    d = float(cfg.nozzle.diameter_m)
    g = cfg.geometry
    y_max = float(g.y_half_over_D) * d
    y_min = 0.0 if g.symmetry_mode == "half_y" else -y_max
    return CaiGeometry(
        diameter_m=d,
        x_max_m=float(g.x_max_over_D) * d,
        y_min_m=y_min,
        y_max_m=y_max,
        z_min_m=-float(g.z_half_over_D) * d,
        z_max_m=float(g.z_half_over_D) * d,
        symmetry_mode=str(g.symmetry_mode),
    )
