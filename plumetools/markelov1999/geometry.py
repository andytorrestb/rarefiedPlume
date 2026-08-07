r"""Geometry of the AIAA 99-3455 flat-plate-in-a-cylinder-wake configuration.

Everything is SI internally. The original inch dimensions are kept alongside, in
:data:`PAPER_DIMENSIONS_IN`, so a reader can check the conversion without
reaching for the paper.

Coordinate convention
---------------------

======================  ====================================================
plume axis              ``+x``
cylinder axis           ``z``
lateral direction       ``y``
origin                  the orifice (source) centre
symmetry plane          ``y = 0``
======================  ====================================================

The full cylinder length and the full plate height are retained in ``z``; only
``y`` is halved. That is the one plane the configuration actually possesses --
``z = 0`` is *not* a symmetry plane once the finite cylinder ends and the finite
plate top and bottom edges are present, and treating it as one would remove the
three-dimensional wake this case exists to study.

Derived, never hard-coded
-------------------------

The plate's ``x`` location is *computed*::

    plate_upstream_x  = cylinder_centre_x + cylinder_radius + gap
    plate_downstream_x = plate_upstream_x + plate_thickness

so changing the gap or the cylinder moves the plate, and no coordinate appears in
the code as an unexplained literal. :meth:`MarkelovGeometry.describe` prints the
arithmetic.

What the paper does not give
----------------------------

Two numbers are **assumptions**, flagged as such everywhere they appear:

``plate_thickness_m``
    The paper does not clearly establish it. The baseline here is 0.5 in
    (0.0127 m) -- thin enough to be a "flat plate" against a 6 in plate width, and
    thick enough that a mesh can resolve it as a solid body rather than a
    zero-thickness baffle. Configurable.

``inflow_radius_m``
    The radius of the hemispherical DSMC source boundary. This is a
    *computational-domain* choice, not a paper value: the paper's physical source
    is the 0.8255 mm orifice, and the hemisphere is simply where the analytical
    source-flow solution is handed to the particle solver. The baseline 0.1524 m
    (6 in) matches the inflow arc radius the archived ``wake-cylinder`` cases
    intended (recorded in ``cases/ARCHIVE.md``, finding AR-01) -- the *value*, not
    that lineage's implementation, which is defective.
    :meth:`MarkelovGeometry.inflow_to_cylinder_clearance_m` proves it does not
    reach the cylinder.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from plumetools.markelov1999.constants import INCH_TO_M, inch_to_m
from plumetools.markelov1999.sourceflow import ORIFICE_DIAMETER_M, ORIFICE_RADIUS_M

#: The dimensions as printed in AIAA 99-3455, in the paper's own units.
#:
#: Kept verbatim so the SI values below can be checked against the source without
#: the paper to hand. Every entry is a *paper value*; assumptions are listed
#: separately in :data:`ASSUMED_DIMENSIONS`.
PAPER_DIMENSIONS_IN = {
    "orifice_diameter_mm": 0.8255,
    "cylinder_diameter_in": 6.0,
    "cylinder_length_in": 18.0,
    "cylinder_centre_x_in": 11.75,
    "plate_width_y_in": 6.0,
    "plate_height_z_in": 15.0,
    "cylinder_to_plate_gap_in": 6.0,
}

#: Values the paper does not provide. See the module docstring.
ASSUMED_DIMENSIONS = {
    "plate_thickness_in": 0.5,
    "inflow_radius_in": 6.0,
}


class GeometryError(ValueError):
    """The configured geometry is inconsistent or physically impossible."""


@dataclass(frozen=True)
class MarkelovGeometry:
    """The full configuration in metres, with every derived coordinate computed.

    Attributes:
        orifice_radius_m: the **physical** orifice radius ``r_e`` [m]. Used by the
            source-flow equations. Paper value: 4.1275e-4 m.
        inflow_radius_m: radius of the hemispherical DSMC source boundary [m].
            An assumption -- see the module docstring. This is the surface the
            model is *evaluated on*; it is never substituted for
            ``orifice_radius_m`` in the equations.
        cylinder_radius_m: [m]. Paper: 3 in = 0.0762 m.
        cylinder_length_m: extent along ``z`` [m]. Paper: 18 in = 0.4572 m.
        cylinder_centre_x_m: [m]. Paper: 11.75 in = 0.29845 m.
        cylinder_centre_z_m: [m]. Zero -- the cylinder is centred on the plume
            axis, which is what "plate centreline aligned with the plume axis"
            plus an axisymmetric source implies.
        gap_m: cylinder downstream surface to plate upstream face [m].
            Paper: 6 in = 0.1524 m.
        plate_width_y_m: full width in ``y`` [m]. Paper: 6 in = 0.1524 m.
        plate_height_z_m: full height in ``z`` [m]. Paper: 15 in = 0.381 m.
        plate_thickness_m: extent in ``x`` [m]. An assumption; baseline 0.0127 m.
        symmetry_plane_y_m: the modelled half is ``y >= this``. Zero.
    """

    orifice_radius_m: float = ORIFICE_RADIUS_M
    inflow_radius_m: float = 6.0 * INCH_TO_M
    cylinder_radius_m: float = 3.0 * INCH_TO_M
    cylinder_length_m: float = 18.0 * INCH_TO_M
    cylinder_centre_x_m: float = 11.75 * INCH_TO_M
    cylinder_centre_z_m: float = 0.0
    gap_m: float = 6.0 * INCH_TO_M
    plate_width_y_m: float = 6.0 * INCH_TO_M
    plate_height_z_m: float = 15.0 * INCH_TO_M
    plate_thickness_m: float = 0.5 * INCH_TO_M
    symmetry_plane_y_m: float = 0.0

    # -- derived x coordinates ---------------------------------------------- #

    @property
    def cylinder_upstream_x_m(self) -> float:
        """Windward (source-facing) surface of the cylinder [m].

        ``cylinder_centre_x - cylinder_radius``. 0.22225 m at paper values.
        """
        return self.cylinder_centre_x_m - self.cylinder_radius_m

    @property
    def cylinder_downstream_x_m(self) -> float:
        """Leeward (base) surface of the cylinder [m].

        ``cylinder_centre_x + cylinder_radius``. 0.37465 m at paper values.
        """
        return self.cylinder_centre_x_m + self.cylinder_radius_m

    @property
    def plate_upstream_x_m(self) -> float:
        """Upstream face of the plate [m] -- **derived, not a literal**.

        ``cylinder_centre_x + cylinder_radius + gap``
        ``= 0.29845 + 0.0762 + 0.1524 = 0.52705 m`` at paper values.
        """
        return self.cylinder_downstream_x_m + self.gap_m

    @property
    def plate_downstream_x_m(self) -> float:
        """Downstream face of the plate [m].

        ``plate_upstream_x + plate_thickness``. Depends on the thickness
        assumption; 0.53975 m at the 0.5 in baseline.
        """
        return self.plate_upstream_x_m + self.plate_thickness_m

    @property
    def plate_centre_x_m(self) -> float:
        """Mid-thickness of the plate [m], for a ``searchableBox`` centre."""
        return 0.5 * (self.plate_upstream_x_m + self.plate_downstream_x_m)

    # -- extents ------------------------------------------------------------- #

    @property
    def cylinder_z_min_m(self) -> float:
        """Lower end cap of the cylinder [m]."""
        return self.cylinder_centre_z_m - 0.5 * self.cylinder_length_m

    @property
    def cylinder_z_max_m(self) -> float:
        """Upper end cap of the cylinder [m]."""
        return self.cylinder_centre_z_m + 0.5 * self.cylinder_length_m

    @property
    def plate_z_min_m(self) -> float:
        """Bottom edge of the plate [m]."""
        return -0.5 * self.plate_height_z_m

    @property
    def plate_z_max_m(self) -> float:
        """Top edge of the plate [m]."""
        return 0.5 * self.plate_height_z_m

    @property
    def plate_y_max_m(self) -> float:
        """Outboard edge of the plate on the modelled half [m].

        The plate straddles ``y = 0``, so the modelled half spans
        ``[0, plate_width_y/2]``.
        """
        return 0.5 * self.plate_width_y_m

    @property
    def cylinder_axis_point(self) -> tuple[float, float, float]:
        """A point on the cylinder axis [m] -- its lower end-cap centre."""
        return (self.cylinder_centre_x_m, 0.0, self.cylinder_z_min_m)

    @property
    def cylinder_axis_end_point(self) -> tuple[float, float, float]:
        """The cylinder's upper end-cap centre [m]."""
        return (self.cylinder_centre_x_m, 0.0, self.cylinder_z_max_m)

    # -- clearances and checks ----------------------------------------------- #

    def inflow_to_cylinder_clearance_m(self) -> float:
        """Minimum distance from the inflow surface to the cylinder surface [m].

        Positive means they are separated; zero or negative means the source
        boundary reaches into the body, which would make the inflow condition
        meaningless there.

        The nearest point of the cylinder to the origin is on its lateral surface
        at ``z`` within the cylinder's span, at distance
        ``cylinder_centre_x - cylinder_radius`` (the origin lies on the plane
        ``z = cylinder_centre_z`` when the cylinder is centred, and the end caps
        are strictly further away). Subtracting the inflow radius gives the gap.

        At baseline values: ``0.22225 - 0.1524 = 0.06985 m`` (2.75 in).
        """
        # Distance from the origin to the cylinder axis, in the plane normal to z.
        axial_offset = abs(self.cylinder_centre_z_m)
        radial = math.hypot(self.cylinder_centre_x_m, 0.0)

        if axial_offset <= 0.5 * self.cylinder_length_m:
            # The origin projects onto the cylinder's lateral surface.
            surface_distance = radial - self.cylinder_radius_m
        else:
            # Past an end cap: nearest point is on the cap rim or face.
            dz = axial_offset - 0.5 * self.cylinder_length_m
            dr = max(0.0, radial - self.cylinder_radius_m)
            surface_distance = math.hypot(dr, dz)

        return surface_distance - self.inflow_radius_m

    def inflow_to_plate_clearance_m(self) -> float:
        """Minimum distance from the inflow surface to the plate [m].

        The plate is far downstream, so this is comfortably positive; it is
        checked anyway because a large ``inflow_radius_m`` is a configurable
        mistake waiting to happen.

        Distance from the origin to the axis-aligned plate box, per axis: the box
        straddles ``y = 0`` and ``z = 0`` by construction (the plate is centred on
        the plume axis), so only the ``x`` offset contributes unless the plate has
        been moved off-centre.
        """
        dx = max(0.0, self.plate_upstream_x_m)
        dy = max(0.0, -self.plate_y_max_m, self.symmetry_plane_y_m - self.plate_y_max_m)
        dz = max(0.0, self.plate_z_min_m, -self.plate_z_max_m)
        return math.sqrt(dx * dx + dy * dy + dz * dz) - self.inflow_radius_m

    def validate(self) -> None:
        """Raise :class:`GeometryError` if the configuration is not physical.

        Checks, in order: every length positive; the orifice inside the inflow
        surface; the bodies placed downstream in the right order; then the
        clearances. Placement is checked before clearance deliberately -- a
        cylinder that swallows the origin makes the clearance number meaningless,
        and reporting "reduce the inflow radius below -0.026 m" would send the
        reader after the wrong parameter.

        These are hard errors, not warnings. A source boundary intersecting a body
        has no meaningful inflow condition, and requirement 10 of the case
        specification lists it among the cases that must fail rather than warn.
        """
        positives = {
            "orifice_radius_m": self.orifice_radius_m,
            "inflow_radius_m": self.inflow_radius_m,
            "cylinder_radius_m": self.cylinder_radius_m,
            "cylinder_length_m": self.cylinder_length_m,
            "gap_m": self.gap_m,
            "plate_width_y_m": self.plate_width_y_m,
            "plate_height_z_m": self.plate_height_z_m,
            "plate_thickness_m": self.plate_thickness_m,
        }
        for name, value in positives.items():
            if not (value > 0.0) or not math.isfinite(value):
                raise GeometryError(f"{name} must be finite and positive, got {value}")

        if self.orifice_radius_m >= self.inflow_radius_m:
            raise GeometryError(
                f"orifice_radius_m ({self.orifice_radius_m}) must be smaller than "
                f"inflow_radius_m ({self.inflow_radius_m}); the source-flow solution "
                f"is handed over on a surface outside the orifice, not inside it"
            )

        if self.cylinder_upstream_x_m <= 0.0:
            raise GeometryError(
                f"the cylinder straddles the source plane: its upstream surface is at "
                f"x = {self.cylinder_upstream_x_m:.6g} m, so the orifice at the origin "
                f"would be inside the body. Move the cylinder downstream "
                f"(bodies.cylinder.centre_x_m > {self.cylinder_radius_m:.6g})."
            )

        if self.plate_upstream_x_m <= self.cylinder_downstream_x_m:
            raise GeometryError(
                f"the plate is not downstream of the cylinder: plate face at "
                f"x = {self.plate_upstream_x_m:.6g} m, cylinder base at "
                f"x = {self.cylinder_downstream_x_m:.6g} m"
            )

        clearance = self.inflow_to_cylinder_clearance_m()
        if clearance <= 0.0:
            raise GeometryError(
                f"the inflow hemisphere (radius {self.inflow_radius_m} m) reaches the "
                f"cylinder: clearance is {clearance:.6g} m. Reduce "
                f"geometry.inflow_radius_m below "
                f"{self.cylinder_upstream_x_m:.6g} m, or move the cylinder downstream."
            )

        plate_clearance = self.inflow_to_plate_clearance_m()
        if plate_clearance <= 0.0:
            raise GeometryError(
                f"the inflow hemisphere reaches the plate: clearance is "
                f"{plate_clearance:.6g} m"
            )

    # -- reporting ------------------------------------------------------------ #

    def describe(self) -> list[str]:
        """Human-readable lines showing every derived coordinate and its arithmetic.

        Printed by ``./Allmesh`` so the plate's ``x`` never appears as a bare
        number anyone has to take on trust.
        """
        to_in = 1.0 / INCH_TO_M
        return [
            "AIAA 99-3455 geometry (SI; inches in brackets)",
            f"  orifice radius r_e      {self.orifice_radius_m:.6e} m "
            f"(diameter {ORIFICE_DIAMETER_M * 1e3:.4f} mm)          [PAPER]",
            f"  inflow hemisphere R     {self.inflow_radius_m:.6f} m "
            f"[{self.inflow_radius_m * to_in:.3f} in]              [ASSUMPTION]",
            f"  cylinder radius         {self.cylinder_radius_m:.6f} m "
            f"[{self.cylinder_radius_m * to_in:.3f} in]              [PAPER]",
            f"  cylinder length (z)     {self.cylinder_length_m:.6f} m "
            f"[{self.cylinder_length_m * to_in:.3f} in]             [PAPER]",
            f"  cylinder centre x       {self.cylinder_centre_x_m:.6f} m "
            f"[{self.cylinder_centre_x_m * to_in:.3f} in]             [PAPER]",
            f"  cylinder upstream  x    {self.cylinder_upstream_x_m:.6f} m  "
            f"= {self.cylinder_centre_x_m:.6f} - {self.cylinder_radius_m:.6f}",
            f"  cylinder downstream x   {self.cylinder_downstream_x_m:.6f} m  "
            f"= {self.cylinder_centre_x_m:.6f} + {self.cylinder_radius_m:.6f}",
            f"  gap                     {self.gap_m:.6f} m "
            f"[{self.gap_m * to_in:.3f} in]              [PAPER]",
            f"  plate upstream  x       {self.plate_upstream_x_m:.6f} m  "
            f"= {self.cylinder_downstream_x_m:.6f} + {self.gap_m:.6f}",
            f"  plate thickness         {self.plate_thickness_m:.6f} m "
            f"[{self.plate_thickness_m * to_in:.3f} in]              [ASSUMPTION]",
            f"  plate downstream x      {self.plate_downstream_x_m:.6f} m  "
            f"= {self.plate_upstream_x_m:.6f} + {self.plate_thickness_m:.6f}",
            f"  plate width  (y)        {self.plate_width_y_m:.6f} m "
            f"[{self.plate_width_y_m * to_in:.3f} in]              [PAPER]",
            f"  plate height (z)        {self.plate_height_z_m:.6f} m "
            f"[{self.plate_height_z_m * to_in:.3f} in]             [PAPER]",
            f"  symmetry plane          y = {self.symmetry_plane_y_m:.6f} m "
            f"(full length in z retained)",
            f"  inflow -> cylinder clearance {self.inflow_to_cylinder_clearance_m():.6f} m "
            f"[{self.inflow_to_cylinder_clearance_m() * to_in:.3f} in]",
            f"  inflow -> plate clearance    {self.inflow_to_plate_clearance_m():.6f} m "
            f"[{self.inflow_to_plate_clearance_m() * to_in:.3f} in]",
        ]

    def to_metadata(self) -> dict:
        """A flat, JSON-serialisable record for manifests and case summaries.

        Carries the SI values, the original inch dimensions, and which of them the
        paper actually prints -- so a downstream reader never has to guess.
        """
        return {
            "si": {
                "orifice_radius_m": self.orifice_radius_m,
                "inflow_radius_m": self.inflow_radius_m,
                "cylinder_radius_m": self.cylinder_radius_m,
                "cylinder_length_m": self.cylinder_length_m,
                "cylinder_centre_x_m": self.cylinder_centre_x_m,
                "cylinder_centre_z_m": self.cylinder_centre_z_m,
                "cylinder_upstream_x_m": self.cylinder_upstream_x_m,
                "cylinder_downstream_x_m": self.cylinder_downstream_x_m,
                "gap_m": self.gap_m,
                "plate_upstream_x_m": self.plate_upstream_x_m,
                "plate_downstream_x_m": self.plate_downstream_x_m,
                "plate_width_y_m": self.plate_width_y_m,
                "plate_height_z_m": self.plate_height_z_m,
                "plate_thickness_m": self.plate_thickness_m,
                "symmetry_plane_y_m": self.symmetry_plane_y_m,
                "inflow_to_cylinder_clearance_m": self.inflow_to_cylinder_clearance_m(),
            },
            "paper_dimensions": dict(PAPER_DIMENSIONS_IN),
            "assumed_dimensions": dict(ASSUMED_DIMENSIONS),
        }


def from_config(cfg) -> MarkelovGeometry:
    """Build a :class:`MarkelovGeometry` from a loaded ``case.yaml``.

    Args:
        cfg: a :class:`plumetools.config.CaseConfig` whose ``bodies`` section
            carries the cylinder and plate, and whose ``geometry`` section carries
            the inflow hemisphere and orifice.

    Returns:
        A validated geometry.

    Raises:
        GeometryError: if the configuration is not physical.

    The inflow hemisphere radius is read from ``geometry.sphere_radius_m`` -- the
    same key the mesh generator uses -- so the surface the model is evaluated on
    and the surface that gets meshed cannot drift apart. ``config.py`` separately
    enforces ``mesh.sphere_radius_m == geometry.sphere_radius_m``.
    """
    geom = MarkelovGeometry(
        orifice_radius_m=cfg.stagnation.throat_radius_m,
        inflow_radius_m=cfg.geometry.sphere_radius_m,
        cylinder_radius_m=cfg.bodies.cylinder.radius_m,
        cylinder_length_m=cfg.bodies.cylinder.length_m,
        cylinder_centre_x_m=cfg.bodies.cylinder.centre_x_m,
        cylinder_centre_z_m=cfg.bodies.cylinder.centre_z_m,
        gap_m=cfg.bodies.gap_m,
        plate_width_y_m=cfg.bodies.plate.width_y_m,
        plate_height_z_m=cfg.bodies.plate.height_z_m,
        plate_thickness_m=cfg.bodies.plate.thickness_m,
        symmetry_plane_y_m=cfg.bodies.symmetry_plane_y_m,
    )
    geom.validate()
    return geom


def cylinder_surface_angles(points, geom: MarkelovGeometry) -> np.ndarray:
    """Azimuth about the cylinder axis [rad], for each point.

    Args:
        points: ``(n, 3)`` positions [m], typically cylinder face centres.
        geom: the configuration.

    Returns:
        ``(n,)`` angles in ``(-pi, pi]``, measured in the ``x-y`` plane from the
        ``+x`` direction about the cylinder axis:

        * ``0``      -- the **leeward** (base, downstream) generator;
        * ``+-pi``   -- the **windward** (stagnation, source-facing) generator.

        The plume travels ``+x``, so the surface facing the source is the one at
        ``-x`` relative to the cylinder centre, i.e. ``|psi| = pi``. Getting this
        backwards silently swaps the windward and leeward pressures, so
        :func:`plumetools.markelov1999.postprocess.pressure_windows` states it
        again at the point of use.
    """
    p = np.asarray(points, dtype=np.float64)
    if p.ndim != 2 or p.shape[1] != 3:
        raise ValueError(f"expected (n, 3) points, got {p.shape}")
    return np.arctan2(p[:, 1], p[:, 0] - geom.cylinder_centre_x_m)


def report_paper_conversions() -> list[str]:
    """Lines showing each paper dimension converted to SI, for documentation."""
    lines = ["AIAA 99-3455 dimensions as printed, converted to SI:"]
    lines.append(f"  orifice diameter   {PAPER_DIMENSIONS_IN['orifice_diameter_mm']} mm"
                 f"   -> {ORIFICE_DIAMETER_M:.6e} m "
                 f"(radius {ORIFICE_RADIUS_M:.6e} m)")
    for key, value in PAPER_DIMENSIONS_IN.items():
        if key.endswith("_in"):
            lines.append(f"  {key[:-3]:<22} {value:>6} in -> {inch_to_m(value):.6f} m")
    for key, value in ASSUMED_DIMENSIONS.items():
        lines.append(f"  {key[:-3]:<22} {value:>6} in -> {inch_to_m(value):.6f} m"
                     f"   [ASSUMPTION -- not a paper value]")
    return lines
