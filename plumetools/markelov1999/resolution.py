r"""Particle-resolution estimate, and the post-run audit that checks it.

The target is ~20 simulated particles per collision cell near the cylinder, in
its immediate wake, and at the plate.

The uniform-weight constraint
-----------------------------
Standard ``dsmcFoam`` has **one** particle weight for the whole domain:
``DSMCCloud::nParticle_`` is a single scalar read from
``nEquivalentParticles``, and ``calculateFields`` multiplies by it globally. There
is no radial weighting, no per-cell weighting, and no adaptive scheme to
configure. (The MNF fork's ``maxRadialWeightingFactor`` has no counterpart here.)

So for a flow whose density falls as ``1/r**2`` on a mesh whose cell size also
varies, the occupancy target can be met in exactly **one** region; every other
region follows from the ratio

    occupancy(A) / occupancy(B) = [n(A) V_cell(A)] / [n(B) V_cell(B)]

which is a property of the flow and the mesh, not of anything configurable.

This module therefore:

* sizes the weight from one named region (``resolution.sizing_region``);
* reports the occupancy that implies in **every** region of interest;
* reports the expected total particle count;
* states the limitation explicitly rather than implying the target is met
  everywhere.

Estimate versus audit
---------------------
:func:`estimate` is analytical -- it evaluates the source-flow density at a
representative point in each region and multiplies by a cell volume. That is
enough to *choose* a weight before meshing, and it is what
``generate_cases.py`` uses. It is **not** evidence that the target was met.

:func:`audit` reads the actual sampled ``dsmcRhoN`` field, which
``DSMCCloud::calculateFields`` fills with the literal parcel count per cell, and
reports the occupancy that really occurred. Only the audit can say whether the
target was achieved.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from plumetools.markelov1999 import sourceflow as sf
from plumetools.markelov1999.geometry import MarkelovGeometry


@dataclass(frozen=True)
class Region:
    """A named place where occupancy is estimated.

    Attributes:
        name: ``"cylinder"``, ``"wake"`` or ``"plate"``.
        point: a representative position [m], where the density is evaluated.
        cell_size_m: the mesh cell size there [m].
        description: what the point represents, for the report.
    """

    name: str
    point: tuple
    cell_size_m: float
    description: str

    @property
    def cell_volume_m3(self) -> float:
        return self.cell_size_m ** 3

    @property
    def radius_m(self) -> float:
        return math.sqrt(sum(c * c for c in self.point))

    @property
    def theta_rad(self) -> float:
        return math.acos(self.point[0] / self.radius_m) if self.radius_m > 0 else 0.0


@dataclass(frozen=True)
class RegionEstimate:
    """Estimated occupancy in one region."""

    region: Region
    number_density_per_m3: float
    molecules_per_cell: float
    particles_per_cell: float
    mean_free_path_m: float

    @property
    def cell_over_mfp(self) -> float:
        """Cell size divided by the local mean free path. DSMC wants this < 1."""
        return self.region.cell_size_m / self.mean_free_path_m


@dataclass(frozen=True)
class ResolutionEstimate:
    """The full pre-run resolution report."""

    n_equivalent_particles: float
    sizing_region: str
    target_particles_per_cell: float
    regions: tuple
    total_particles: float
    total_molecules: float
    domain_volume_m3: float

    def by_name(self, name: str) -> RegionEstimate:
        for r in self.regions:
            if r.region.name == name:
                return r
        raise KeyError(f"no region {name!r}; have {[r.region.name for r in self.regions]}")

    @property
    def min_particles_per_cell(self) -> float:
        return min(r.particles_per_cell for r in self.regions)

    def report(self) -> list[str]:
        """Human-readable lines, including the limitation."""
        lines = [
            "Particle resolution (ANALYTICAL ESTIMATE -- not a measurement)",
            f"  particle weight          {self.n_equivalent_particles:.6e} "
            f"molecules per simulated particle",
            f"  sized on                 '{self.sizing_region}' at a target of "
            f"{self.target_particles_per_cell:g} particles/cell",
            "",
            f"  {'region':<10} {'cell [m]':>10} {'n [1/m^3]':>12} "
            f"{'part/cell':>10} {'cell/mfp':>9}",
        ]
        for r in self.regions:
            lines.append(
                f"  {r.region.name:<10} {r.region.cell_size_m:>10.6f} "
                f"{r.number_density_per_m3:>12.4e} {r.particles_per_cell:>10.1f} "
                f"{r.cell_over_mfp:>9.2f}")
        lines += [
            "",
            f"  expected total particles {self.total_particles:.4e} "
            f"(domain volume {self.domain_volume_m3:.4f} m^3)",
            f"  minimum occupancy across the regions of interest: "
            f"{self.min_particles_per_cell:.1f} particles/cell",
        ]

        shortfall = [r.region.name for r in self.regions
                     if r.particles_per_cell < self.target_particles_per_cell]
        if shortfall:
            lines += [
                "",
                f"  NOTE: {shortfall} fall below the "
                f"{self.target_particles_per_cell:g}-particle target.",
                "  Standard dsmcFoam has ONE particle weight for the whole domain",
                "  (DSMCCloud::nParticle_ is a single scalar), so the target can be",
                "  met in exactly one region and the rest follow from the density and",
                "  cell-size ratios. Set resolution.sizing_region to the sparsest",
                "  region to meet it everywhere, at a proportionally higher particle",
                "  count. No adaptive weighting is available to avoid the trade.",
            ]
        lines += [
            "",
            "  This is an ESTIMATE. Run ./Allpost after solving for the audit against",
            "  the sampled dsmcRhoN field, which is the only evidence of what the",
            "  occupancy actually was.",
        ]
        return lines


def mean_free_path(number_density_per_m3: float, diameter_m: float) -> float:
    """VHS mean free path at the reference temperature [m].

    ``lambda = 1 / (sqrt(2) * pi * d**2 * n)`` -- the hard-sphere expression, used
    here with the VHS reference diameter. That is an approximation: the true VHS
    mean free path carries a ``(Tref/T)**(omega-1/2)`` factor. It is not corrected
    because the temperature this expansion reaches varies by orders of magnitude
    across the domain and the check that consumes this number is a warning
    threshold, not a physical result. Documented rather than silently applied.
    """
    if number_density_per_m3 <= 0.0:
        return float("inf")
    return 1.0 / (math.sqrt(2.0) * math.pi * diameter_m ** 2 * number_density_per_m3)


def regions_of_interest(cfg, geom: MarkelovGeometry, cell_sizes: dict) -> list[Region]:
    """The three regions the resolution target applies to.

    Args:
        cfg: the case config.
        geom: the geometry.
        cell_sizes: ``{"inflow": ..., "cylinder": ..., "plate": ...}`` [m], from
            :func:`plumetools.markelov1999.mesh.surface_cell_sizes`.

    Returns:
        Cylinder stagnation, immediate wake, and plate stagnation.

    The wake point sits ``resolution.wake_offset_m`` behind the cylinder base on
    the axis. Its cell size is the cylinder's, since it is inside that surface's
    refinement halo.
    """
    offset = float(cfg.resolution.wake_offset_m)
    return [
        Region(
            name="cylinder",
            point=(geom.cylinder_upstream_x_m, 0.0, 0.0),
            cell_size_m=cell_sizes["cylinder"],
            description="cylinder windward stagnation point",
        ),
        Region(
            name="wake",
            point=(geom.cylinder_downstream_x_m + offset, 0.0, 0.0),
            cell_size_m=cell_sizes["cylinder"],
            description=f"{offset:.4f} m behind the cylinder base, on the axis",
        ),
        Region(
            name="plate",
            point=(geom.plate_upstream_x_m, 0.0, 0.0),
            cell_size_m=cell_sizes["plate"],
            description="plate upstream face, on the axis",
        ),
    ]


def number_density_at(point, cfg) -> float:
    """Source-flow number density at a point [1/m^3].

    The **undisturbed** source-flow value: it ignores the bodies entirely, so the
    density it reports at the cylinder stagnation point is what would be there if
    the cylinder were not. That is the right basis for sizing -- the compressed
    stagnation region is denser, so the estimate is conservative -- but it is
    another reason the estimate is not a measurement.
    """
    from plumetools.markelov1999.inflow import molecular_mass_from_config

    r = math.sqrt(sum(c * c for c in point))
    theta = math.acos(point[0] / r) if r > 0 else 0.0
    return float(sf.number_density(
        theta, r,
        gamma=cfg.gas.gamma,
        T0_K=cfg.stagnation.T0_K,
        mass_kg=molecular_mass_from_config(cfg),
        p0_pa=cfg.stagnation.p0_pa,
        orifice_radius_m=cfg.stagnation.throat_radius_m,
        exponent_offset=cfg.angular.exponent_offset,
        quadrature_points=cfg.angular.quadrature_points,
    ))


def total_molecules_in_domain(cfg, geom: MarkelovGeometry) -> float:
    """Molecules in the domain at steady state, by radial integration.

    The source flow is radial with ``n ~ 1/r**2``, so the molecule count in a cone
    element between ``r1`` and ``r2`` is ``n(r1) r1**2 (r2 - r1)`` per steradian.
    Integrating the angular profile over the modelled half domain gives

        N = n_axis(R) R**2 * Omega_half * (r_max - R)

    with ``Omega_half = pi * integral(sin(t) f(t) dt)`` -- half of the full
    ``2*pi`` azimuth because ``y >= 0``.

    This is a **collisionless** estimate. It ignores the bodies, which block part
    of the cone, and it ignores compression at the stagnation regions. It is used
    only to report an expected order of magnitude for the particle count, so that
    a configuration asking for 10^9 particles is caught before it is submitted.
    """
    from plumetools.markelov1999.inflow import molecular_mass_from_config

    R = geom.inflow_radius_m
    n_axis = float(sf.number_density(
        0.0, R,
        gamma=cfg.gas.gamma,
        T0_K=cfg.stagnation.T0_K,
        mass_kg=molecular_mass_from_config(cfg),
        p0_pa=cfg.stagnation.p0_pa,
        orifice_radius_m=cfg.stagnation.throat_radius_m,
        exponent_offset=cfg.angular.exponent_offset,
        quadrature_points=cfg.angular.quadrature_points,
    ))
    omega_half = math.pi * sf.normalization_integral(
        cfg.gas.gamma, cfg.angular.exponent_offset, cfg.angular.quadrature_points)

    # Representative outer radius: the far corner of the domain understates the
    # typical path, the near face overstates it, so use the downstream extent.
    r_max = float(cfg.mesh.x_max_m)
    return n_axis * R * R * omega_half * max(0.0, r_max - R)


def estimate(cfg, geom: MarkelovGeometry, cell_sizes: dict,
             n_equivalent_particles: float | None = None) -> ResolutionEstimate:
    """Estimate occupancy in every region, and pick a weight if none is given.

    Args:
        cfg: the case config.
        geom: the geometry.
        cell_sizes: cell size at each refined surface [m].
        n_equivalent_particles: an explicit weight. ``None`` derives one so that
            ``resolution.sizing_region`` reaches the target.

    Returns:
        A :class:`ResolutionEstimate`.

    Raises:
        ValueError: if the sizing region has zero density, which would make the
            derived weight infinite.
    """
    regions = regions_of_interest(cfg, geom, cell_sizes)
    diameter = cfg.dsmc.species.diameter_m
    target = float(cfg.resolution.target_particles_per_cell)

    densities = {r.name: number_density_at(r.point, cfg) for r in regions}

    if n_equivalent_particles is None:
        sizing = cfg.resolution.sizing_region
        chosen = next(r for r in regions if r.name == sizing)
        n_sizing = densities[sizing]
        if n_sizing <= 0.0:
            raise ValueError(
                f"the source-flow density at the '{sizing}' region is {n_sizing}; "
                f"a particle weight cannot be derived from it. Is the region "
                f"outside the plume cone?")
        n_equivalent_particles = n_sizing * chosen.cell_volume_m3 / target

    estimates = tuple(
        RegionEstimate(
            region=r,
            number_density_per_m3=densities[r.name],
            molecules_per_cell=densities[r.name] * r.cell_volume_m3,
            particles_per_cell=densities[r.name] * r.cell_volume_m3
            / n_equivalent_particles,
            mean_free_path_m=mean_free_path(densities[r.name], diameter),
        )
        for r in regions
    )

    total_molecules = total_molecules_in_domain(cfg, geom)
    volume = ((cfg.mesh.x_max_m - cfg.mesh.x_min_m)
              * (cfg.mesh.y_max_m - cfg.mesh.y_min_m)
              * (cfg.mesh.z_max_m - cfg.mesh.z_min_m))

    return ResolutionEstimate(
        n_equivalent_particles=float(n_equivalent_particles),
        sizing_region=cfg.resolution.sizing_region,
        target_particles_per_cell=target,
        regions=estimates,
        total_particles=total_molecules / n_equivalent_particles,
        total_molecules=total_molecules,
        domain_volume_m3=float(volume),
    )


# --------------------------------------------------------------------------- #
# post-run audit
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class AuditResult:
    """Measured occupancy, from the sampled ``dsmcRhoN`` field."""

    time: str
    field: str
    n_cells: int
    mean_particles_per_cell: float
    median_particles_per_cell: float
    max_particles_per_cell: float
    occupied_cells: int
    total_particles: float
    region_occupancy: dict
    target: float

    def meets_target(self) -> bool:
        """True only if every region of interest reached the target."""
        return bool(self.region_occupancy) and all(
            v >= self.target for v in self.region_occupancy.values())

    def report(self) -> list[str]:
        lines = [
            f"Particle resolution AUDIT (measured, from {self.field} at t = {self.time})",
            f"  cells                    {self.n_cells}",
            f"  occupied cells           {self.occupied_cells} "
            f"({100.0 * self.occupied_cells / max(self.n_cells, 1):.1f}%)",
            f"  total particles          {self.total_particles:.4e}",
            f"  particles/cell  mean {self.mean_particles_per_cell:.2f}  "
            f"median {self.median_particles_per_cell:.2f}  "
            f"max {self.max_particles_per_cell:.0f}",
        ]
        if self.region_occupancy:
            lines.append("  regions of interest:")
            for name, value in sorted(self.region_occupancy.items()):
                verdict = "OK " if value >= self.target else "LOW"
                lines.append(f"    {verdict} {name:<10} {value:>8.1f} particles/cell "
                             f"(target {self.target:g})")
        verdict = "MET" if self.meets_target() else "NOT met"
        lines.append(f"  target of {self.target:g} particles/cell: {verdict}")
        return lines


def audit(cell_particle_counts: np.ndarray, cell_centres: np.ndarray,
          cfg, geom: MarkelovGeometry, *, time: str, field: str,
          region_radius_m: float = 0.02) -> AuditResult:
    """Measure occupancy from a sampled per-cell particle count.

    Args:
        cell_particle_counts: ``(n_cells,)`` parcels per cell. This is exactly
            what ``dsmcRhoN`` holds: ``DSMCCloud::calculateFields`` increments it
            once per parcel and, unlike every other measurement field, never
            multiplies by ``nParticle`` or divides by the cell volume.
        cell_centres: ``(n_cells, 3)`` cell centres [m].
        cfg: the case config.
        geom: the geometry.
        time: the time directory the field came from, for the report.
        field: the field name, for the report.
        region_radius_m: how far around each region's representative point to
            average.

    Returns:
        An :class:`AuditResult`.

    Raises:
        ValueError: on mismatched array lengths.
    """
    counts = np.asarray(cell_particle_counts, dtype=np.float64)
    centres = np.asarray(cell_centres, dtype=np.float64)
    if centres.shape != (len(counts), 3):
        raise ValueError(
            f"got {len(counts)} cell values but {centres.shape} cell centres")

    from plumetools.markelov1999.mesh import surface_cell_sizes
    cell_sizes = surface_cell_sizes(cfg)
    regions = regions_of_interest(cfg, geom, cell_sizes)

    occupancy: dict[str, float] = {}
    for region in regions:
        distance = np.linalg.norm(centres - np.asarray(region.point), axis=1)
        near = distance <= region_radius_m
        if np.any(near):
            occupancy[region.name] = float(counts[near].mean())

    return AuditResult(
        time=time,
        field=field,
        n_cells=len(counts),
        mean_particles_per_cell=float(counts.mean()) if len(counts) else 0.0,
        median_particles_per_cell=float(np.median(counts)) if len(counts) else 0.0,
        max_particles_per_cell=float(counts.max()) if len(counts) else 0.0,
        occupied_cells=int(np.count_nonzero(counts > 0)),
        total_particles=float(counts.sum()),
        region_occupancy=occupancy,
        target=float(cfg.resolution.target_particles_per_cell),
    )
