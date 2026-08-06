r"""Mesh generation for the AIAA 99-3455 configuration.

Pipeline::

    blockMesh                 uniform Cartesian background box
    snappyHexMesh -overwrite  carve the inflow hemisphere, the cylinder and the
                              plate out of it, snapping to analytic primitives
    checkMesh                 mesh quality
    python -m plumetools.markelov1999.verify   the geometry the model assumes

Composed from :mod:`plumetools.foamio.primitives`, shared with the single-sphere
generator, so this module is the *geometry* of the case and not a second copy of
the dictionary syntax.

Patches
-------

============================  ======================  =======================
patch                         geometric type          source
============================  ======================  =======================
``inflow``                    ``patch``               snappyHexMesh
``cylinder``                  ``wall``                snappyHexMesh
``plate``                     ``wall``                snappyHexMesh
``symmetry``                  ``symmetry``            blockMesh, ``y = y_min``
``upstreamVacuum``            ``patch``               blockMesh, ``x = x_min``
``vacuum``                    ``patch``               blockMesh, the other four
============================  ======================  =======================

The geometric types are load-bearing, not labels:

* ``patch`` -- ``particle::hitBoundaryFace`` finds no special handler and sets
  ``td.keepParticle = false``, so the particle is **deleted**. That is what an
  open vacuum boundary must do, and it is why no custom deleting boundary model
  is needed once inflow is restricted to the inflow patch.
* ``wall`` -- ``DSMCParcel::hitWallPatch`` runs, recording the surface fluxes
  (``fD``, ``q``, ``rhoN``, ...) and applying the ``WallInteractionModel``. A body
  declared as ``patch`` would absorb every particle and exert no force at all.
* ``symmetry`` -- ``particle::hitSymmetryPatch`` reflects, which is what a plane
  of symmetry means.

``upstreamVacuum`` is deliberately **not** called a symmetry plane. The
configuration is symmetric about ``y = 0`` and about nothing else; ``x = 0`` is
simply where the source sits, and a particle that scatters back through it has
left the domain. Naming it ``symmetry`` would reflect those particles back in.
"""

from __future__ import annotations

from pathlib import Path

from plumetools.foamio.primitives import (
    BoxPatch,
    SearchableBox,
    SearchableCylinder,
    SearchableSphere,
    check_surfaces_are_resolvable,
    fmt,
    render_box_block_mesh_dict,
    render_mesh_quality_dict,
    render_snappy_dict,
    uniform_divisions,
    write_lines,
)
from plumetools.markelov1999.geometry import MarkelovGeometry, from_config

#: Dotted path recorded in every dictionary this module writes.
GENERATOR = "plumetools.markelov1999.mesh"


def domain_bounds(cfg) -> tuple[float, float, float, float, float, float]:
    """``(x_min, x_max, y_min, y_max, z_min, z_max)`` from ``mesh`` [m]."""
    m = cfg.mesh
    bounds = (m.x_min_m, m.x_max_m, m.y_min_m, m.y_max_m, m.z_min_m, m.z_max_m)
    if any(v is None for v in bounds):
        raise ValueError(
            "mesh.{x,y,z}_{min,max}_m are all required for mesh.type "
            "'snappy_markelov'; load the config through "
            "plumetools.config.load_case_config, which checks this"
        )
    return tuple(float(v) for v in bounds)


def check_domain_contains_geometry(cfg, geom: MarkelovGeometry) -> list[str]:
    """Verify the box encloses the source, both bodies and a downstream buffer.

    Returns:
        Report lines describing the clearance on each side.

    Raises:
        ValueError: if any body pokes out of the box, or the source is not inside
            it. Either makes snappyHexMesh's result depend on where the surface
            happens to be clipped, which is not a geometry anyone chose.

    The requirement that a *body* be strictly inside matters more than it looks: a
    cylinder end cap flush with the ``z`` boundary would be meshed as though the
    cylinder continued forever, removing the end effects this case exists to
    study.
    """
    x0, x1, y0, y1, z0, z1 = domain_bounds(cfg)
    report: list[str] = []

    def require(ok: bool, message: str) -> None:
        if not ok:
            raise ValueError(message)

    require(x0 <= 0.0 < x1,
            f"the source at x = 0 is outside the domain x in [{x0}, {x1}]")
    require(y0 <= geom.symmetry_plane_y_m,
            f"the symmetry plane y = {geom.symmetry_plane_y_m} is outside "
            f"y in [{y0}, {y1}]")

    require(geom.plate_downstream_x_m < x1,
            f"the plate (to x = {geom.plate_downstream_x_m:.6g}) reaches or exceeds "
            f"the downstream boundary x_max = {x1}; increase mesh.x_max_m")
    require(geom.cylinder_z_min_m > z0 and geom.cylinder_z_max_m < z1,
            f"the cylinder spans z {geom.cylinder_z_min_m:.6g} .. "
            f"{geom.cylinder_z_max_m:.6g}, which is not strictly inside "
            f"z in [{z0}, {z1}]; its end caps must be in the domain or the wake "
            f"around them cannot form")
    require(geom.plate_z_min_m > z0 and geom.plate_z_max_m < z1,
            f"the plate spans z {geom.plate_z_min_m:.6g} .. {geom.plate_z_max_m:.6g}, "
            f"which is not strictly inside z in [{z0}, {z1}]")
    require(geom.plate_y_max_m < y1,
            f"the plate reaches y = {geom.plate_y_max_m:.6g}, outside y_max = {y1}")
    require(geom.cylinder_centre_x_m + geom.cylinder_radius_m < x1,
            "the cylinder reaches the downstream boundary")
    require(geom.inflow_radius_m < min(x1, y1 - y0, z1 - z0),
            "the inflow hemisphere does not fit inside the domain")

    report.append(
        f"  domain x [{fmt(x0)}, {fmt(x1)}]  y [{fmt(y0)}, {fmt(y1)}]  "
        f"z [{fmt(z0)}, {fmt(z1)}] m")
    report.append(
        f"  downstream buffer past the plate  "
        f"{x1 - geom.plate_downstream_x_m:.4f} m "
        f"({(x1 - geom.plate_downstream_x_m) / (2 * geom.cylinder_radius_m):.2f} "
        f"cylinder diameters)")
    report.append(
        f"  lateral buffer past the plate edge {y1 - geom.plate_y_max_m:.4f} m")
    report.append(
        f"  axial buffer past the cylinder cap {z1 - geom.cylinder_z_max_m:.4f} m")
    return report


def refinement_levels(cfg) -> dict[str, int]:
    """Per-surface octree levels, each falling back to ``mesh.refinement_level``."""
    m = cfg.mesh
    return {
        "inflow": int(m.inflow_refinement_level
                      if m.inflow_refinement_level is not None else m.refinement_level),
        "cylinder": int(m.cylinder_refinement_level
                        if m.cylinder_refinement_level is not None else m.refinement_level),
        "plate": int(m.plate_refinement_level
                     if m.plate_refinement_level is not None else m.refinement_level),
    }


def surface_cell_sizes(cfg) -> dict[str, float]:
    """Cell size at each carved surface [m].

    ``background_cell_size_m / 2**level``, using the *actual* background cell size
    the divisions produce rather than the requested one -- rounding to an integer
    cell count can move it by a few percent, and every downstream estimate
    (particles per cell, Courant number, cell-to-mean-free-path ratio) should be
    based on the mesh that will exist.
    """
    actual = background_cell_size(cfg)
    return {name: actual / 2 ** level for name, level in refinement_levels(cfg).items()}


def background_cell_size(cfg) -> float:
    """The background cell size the integer divisions actually give [m].

    The largest of the three axis spacings: it is the one that decides whether a
    surface can be found, and reporting the smallest would flatter the mesh.
    """
    x0, x1, y0, y1, z0, z1 = domain_bounds(cfg)
    nx, ny, nz = uniform_divisions(domain_bounds(cfg), cfg.mesh.background_cell_size_m)
    return max((x1 - x0) / nx, (y1 - y0) / ny, (z1 - z0) / nz)


def searchable_surfaces(cfg, geom: MarkelovGeometry | None = None):
    """The three carved surfaces, in a deterministic order.

    Args:
        cfg: the case config.
        geom: the geometry; derived from ``cfg`` if omitted.

    Returns:
        ``[inflow sphere, cylinder, plate]``.

    Order is fixed rather than incidental: it decides the order of the
    ``geometry`` entries, hence the order snappyHexMesh creates the patches, hence
    the patch order in ``constant/polyMesh/boundary``. A mesh regenerated from the
    same config must produce the same patch indices, or every decomposed run and
    every stored field ordering shifts underneath.
    """
    geom = geom or from_config(cfg)
    names = cfg.mesh.patch_names
    levels = refinement_levels(cfg)

    return [
        SearchableSphere(
            name=names["inflow"],
            level=levels["inflow"],
            # 'patch', never 'wall': the inflow is where particles are created,
            # and a wall patch would make it a gas-surface interaction site.
            patch_type="patch",
            centre=(0.0, 0.0, 0.0),
            radius=geom.inflow_radius_m,
        ),
        SearchableCylinder(
            name=names["cylinder"],
            level=levels["cylinder"],
            patch_type="wall",
            point1=geom.cylinder_axis_point,
            point2=geom.cylinder_axis_end_point,
            radius=geom.cylinder_radius_m,
        ),
        SearchableBox(
            name=names["plate"],
            level=levels["plate"],
            patch_type="wall",
            # The full-width box is intersected by the domain at y = 0, so the
            # modelled half comes out of the symmetry plane rather than out of a
            # half-width box that would leave a spurious internal face at y = 0.
            min=(geom.plate_upstream_x_m, -geom.plate_y_max_m, geom.plate_z_min_m),
            max=(geom.plate_downstream_x_m, geom.plate_y_max_m, geom.plate_z_max_m),
        ),
    ]


def location_in_mesh(cfg, geom: MarkelovGeometry | None = None) -> tuple[float, float, float]:
    """A seed point in the fluid: inside the box, outside every carved body.

    snappyHexMesh keeps the region connected to this point and discards the rest,
    so a seed inside a body, on a face, on an edge or on a symmetry plane produces
    either the wrong region or an outright failure.

    The point is snapped to the **centre of the background cell containing it**. A
    cell centre is interior by construction whatever the cell size, so the seed
    tracks the grid instead of drifting onto it -- a fixed fraction of the domain
    lands exactly on a cell face whenever it happens to be a multiple of the cell
    spacing, which is common.

    The candidate is placed in the lateral far field, off the plume axis and away
    from every body: outboard of the plate edge in ``y`` and upstream of the
    cylinder in ``x``, where nothing is carved.
    """
    geom = geom or from_config(cfg)
    x0, x1, y0, y1, z0, z1 = domain_bounds(cfg)
    nx, ny, nz = uniform_divisions(domain_bounds(cfg), cfg.mesh.background_cell_size_m)

    def cell_centre(value: float, lo: float, hi: float, n: int) -> float:
        h = (hi - lo) / n
        index = min(n - 1, max(0, int((value - lo) / h)))
        return lo + (index + 0.5) * h

    # Lateral far field: outboard of the plate and the cylinder, level with the
    # cylinder mid-span but well clear of it in y.
    candidate = (
        0.5 * (geom.cylinder_upstream_x_m + x1),
        0.5 * (max(geom.plate_y_max_m, geom.cylinder_radius_m) + y1),
        0.5 * (geom.cylinder_z_max_m + z1),
    )
    point = (
        cell_centre(candidate[0], x0, x1, nx),
        cell_centre(candidate[1], y0, y1, ny),
        cell_centre(candidate[2], z0, z1, nz),
    )

    _assert_point_is_in_the_fluid(point, geom, cfg)
    return point


def _assert_point_is_in_the_fluid(point, geom: MarkelovGeometry, cfg) -> None:
    """Raise if the seed point is inside a body, or on a boundary of the domain."""
    x, y, z = point
    x0, x1, y0, y1, z0, z1 = domain_bounds(cfg)

    if not (x0 < x < x1 and y0 < y < y1 and z0 < z < z1):
        raise ValueError(f"locationInMesh {point} is on or outside the domain boundary")

    if (x * x + y * y + z * z) ** 0.5 <= geom.inflow_radius_m:
        raise ValueError(
            f"locationInMesh {point} is inside the inflow hemisphere "
            f"(radius {geom.inflow_radius_m}); snappyHexMesh would keep the cavity "
            f"and discard the fluid")

    in_cylinder_span = geom.cylinder_z_min_m <= z <= geom.cylinder_z_max_m
    radial = ((x - geom.cylinder_centre_x_m) ** 2 + y * y) ** 0.5
    if in_cylinder_span and radial <= geom.cylinder_radius_m:
        raise ValueError(f"locationInMesh {point} is inside the cylinder")

    if (geom.plate_upstream_x_m <= x <= geom.plate_downstream_x_m
            and abs(y) <= geom.plate_y_max_m
            and geom.plate_z_min_m <= z <= geom.plate_z_max_m):
        raise ValueError(f"locationInMesh {point} is inside the plate")

    if y == geom.symmetry_plane_y_m:
        raise ValueError(
            f"locationInMesh {point} lies on the symmetry plane "
            f"y = {geom.symmetry_plane_y_m}")


# --------------------------------------------------------------------------- #
# dictionaries
# --------------------------------------------------------------------------- #

def render_background_block_mesh_dict(cfg, geom: MarkelovGeometry | None = None) -> str:
    """The uniform Cartesian background box, with its three blockMesh patches."""
    geom = geom or from_config(cfg)
    bounds = domain_bounds(cfg)
    nx, ny, nz = uniform_divisions(bounds, cfg.mesh.background_cell_size_m)
    names = cfg.mesh.patch_names
    sizes = surface_cell_sizes(cfg)

    return render_box_block_mesh_dict(
        bounds=bounds,
        divisions=(nx, ny, nz),
        patches=[
            BoxPatch(
                names["symmetry"], "symmetry", ("y_min",),
                f"y = {fmt(bounds[2])}. The ONE plane of symmetry this "
                f"configuration has. z is NOT halved: the cylinder end caps and "
                f"the plate edges are inside the domain."),
            BoxPatch(
                names["upstream_vacuum"], "patch", ("x_min",),
                f"x = {fmt(bounds[0])}, outside the inflow cavity snappyHexMesh "
                f"carves from this face. An OPEN boundary, not a symmetry plane: "
                f"particles that scatter upstream have left the domain."),
            BoxPatch(
                names["outer"], "patch", ("x_max", "y_max", "z_min", "z_max"),
                "the far field. 'patch', so particles are deleted rather than "
                "reflected -- see the module docstring."),
        ],
        generator=GENERATOR,
        note=(f"Background box for snappyHexMesh: {nx}x{ny}x{nz} = {nx * ny * nz} "
              f"cells of ~{fmt(background_cell_size(cfg))} m. snappyHexMesh refines "
              f"to {fmt(sizes['inflow'])} m at the inflow, "
              f"{fmt(sizes['cylinder'])} m at the cylinder and "
              f"{fmt(sizes['plate'])} m at the plate."),
    )


def render_snappy_hex_mesh_dict(cfg, geom: MarkelovGeometry | None = None) -> str:
    """The three searchable primitives plus castellate and snap controls."""
    geom = geom or from_config(cfg)
    surfaces = searchable_surfaces(cfg, geom)
    check_surfaces_are_resolvable(surfaces, background_cell_size(cfg))
    sizes = surface_cell_sizes(cfg)

    across_plate = geom.plate_thickness_m / sizes["plate"]
    around_cylinder = (2 * 3.141592653589793 * geom.cylinder_radius_m) / sizes["cylinder"]

    return render_snappy_dict(
        surfaces,
        location_in_mesh=location_in_mesh(cfg, geom),
        n_cells_between_levels=int(cfg.mesh.n_cells_between_levels),
        generator=GENERATOR,
        note=(f"AIAA 99-3455: carve the inflow hemisphere (R = "
              f"{fmt(geom.inflow_radius_m)} m), the cylinder (R = "
              f"{fmt(geom.cylinder_radius_m)} m, L = {fmt(geom.cylinder_length_m)} m) "
              f"and the plate (x = {fmt(geom.plate_upstream_x_m)} .. "
              f"{fmt(geom.plate_downstream_x_m)} m). "
              f"~{around_cylinder:.0f} cells around the cylinder, "
              f"~{across_plate:.1f} across the plate thickness."),
        preamble=[
            "// Analytic primitives, no STL: snappyHexMesh snaps to the true surface",
            "// rather than to a triangulation of it, and nothing outside the",
            "// repository has to be generated or versioned. The geometry entry name",
            "// becomes the patch name.",
            "//",
            "// patchInfo types are load-bearing: 'patch' DELETES particles that reach",
            "// it (correct for the inflow surface), 'wall' runs hitWallPatch, which",
            "// records fD/q and applies the WallInteractionModel (required for a",
            "// solid body -- a body declared 'patch' would absorb the plume and",
            "// exert no force).",
        ],
    )


def write_mesh_setup(case_dir: Path, cfg) -> list[Path]:
    """Write ``blockMeshDict``, ``snappyHexMeshDict`` and ``meshQualityDict``.

    Args:
        case_dir: an OpenFOAM case directory.
        cfg: a config with ``mesh.type == "snappy_markelov"``.

    Returns:
        The paths written, in the order blockMesh / snappyHexMesh consume them.

    Raises:
        ValueError: for the wrong ``mesh.type``, a geometry outside the domain, or
            a background too coarse to find one of the surfaces.
    """
    if cfg.mesh.type != "snappy_markelov":
        raise ValueError(
            f"write_mesh_setup needs mesh.type 'snappy_markelov', got {cfg.mesh.type!r}")

    geom = from_config(cfg)
    check_domain_contains_geometry(cfg, geom)

    system = Path(case_dir) / "system"
    return [
        write_lines(system / "blockMeshDict",
                    render_background_block_mesh_dict(cfg, geom).splitlines()),
        write_lines(system / "snappyHexMeshDict",
                    render_snappy_hex_mesh_dict(cfg, geom).splitlines()),
        write_lines(system / "meshQualityDict",
                    render_mesh_quality_dict(GENERATOR).splitlines()),
    ]


def describe(cfg, geom: MarkelovGeometry | None = None) -> list[str]:
    """Report lines about the mesh this config will produce, for ``./Allmesh``."""
    geom = geom or from_config(cfg)
    nx, ny, nz = uniform_divisions(domain_bounds(cfg), cfg.mesh.background_cell_size_m)
    sizes = surface_cell_sizes(cfg)
    levels = refinement_levels(cfg)

    lines = [
        f"  background  {nx} x {ny} x {nz} = {nx * ny * nz} cells "
        f"of ~{fmt(background_cell_size(cfg))} m",
    ]
    for name in ("inflow", "cylinder", "plate"):
        lines.append(
            f"  {name:<10}  level {levels[name]} -> {fmt(sizes[name])} m")
    lines.append(
        f"  plate thickness spans {geom.plate_thickness_m / sizes['plate']:.1f} cells")
    lines.append(
        f"  cylinder circumference spans "
        f"{2 * 3.141592653589793 * geom.cylinder_radius_m / sizes['cylinder']:.0f} cells")
    lines += check_domain_contains_geometry(cfg, geom)
    lines.append(f"  locationInMesh {location_in_mesh(cfg, geom)}")
    return lines
