r"""Generate a snappyHexMesh setup for the box-with-hemispherical-cavity geometry.

The generator ``mesh.type: snappy_hex_sphere`` selects, and what
``cases/3d-inflow`` uses. The alternative is the projected O-grid in
:mod:`plumetools.foamio.blockmesh`.

The sphere is declared as a ``searchableSphere`` **primitive**, so no STL is
needed and the surface is analytic: snappyHexMesh snaps to the true sphere rather
than to a triangulation of it.

Three files are written:

    system/blockMeshDict      a plain uniform box -- the background mesh, and the
                              ONLY thing blockMesh builds: no sphere, no inflow
                              patch, no curved edges
    system/snappyHexMeshDict  the searchableSphere and the castellate/snap controls
    system/meshQualityDict    quality limits, #included by the above

and the pipeline is::

    blockMesh                 # background box
    snappyHexMesh -overwrite  # carve the cavity, snap to the sphere

How the cavity appears
----------------------
The sphere is centred at the origin, which lies **on** the ``x = 0`` face of the
background box (``x in [0, L]``). Only the ``+x`` half of the sphere therefore
intersects the mesh, and that is what gets carved out -- the hemisphere falls out
of the geometry rather than needing to be constructed. ``locationInMesh`` is
placed in the fluid, outside the sphere, so snappyHexMesh keeps the region
connected to it and discards the cavity interior.

Resolution
----------
Only ``background_cell_size_m`` sizes the box. Surface resolution comes from the
octree::

    surface cell size = background_cell_size_m / 2**refinement_level

so a coarse background with more levels concentrates cells at the source and
leaves the near-vacuum far field coarse -- which is what the O-grid's
``radial_grading`` was doing by hand. The background only has to *find* the
sphere, so a few background cells across the diameter is enough; it does not have
to resolve it.

Trade-off against the projected O-grid
--------------------------------------
The projected O-grid (``mesh.type: block_mesh_ogrid`` with
``mesh.projection: searchable_sphere``) gives an exact sphere, a pure hexahedral
mesh, and a face count of exactly ``5*n^2``.

What it does not give is **cell quality**. Five blocks meet at the cube corners
projected onto the sphere, and ``radial_grading`` shears every cell in the
graded direction: at the committed 20/24/10 settings ``checkMesh`` reports a max
non-orthogonality of **66.1** (against a 70 limit) and a mean of **31.2**.
Nothing about that improves with refinement -- it is the block topology.

snappyHexMesh keeps the background hexes axis-aligned and orthogonal everywhere
except the one cell layer it snaps. Measured on the same geometry at 0.25 m /
level 4: max non-orthogonality **36.9**, mean **11.9**, max skewness **0.76**
against 1.79, in **32 272** cells against 48 000. The price:

* **Polyhedral cells and polygonal faces near the surface.** Castellation splits
  hexes and snapping produces faces with more than four vertices.
  :func:`plumetools.geometry.centroid` divides by ``len(vertices)`` so it handles
  those correctly -- but the legacy ``/3.0`` divisor would be wrong on every one
  of them (finding AR-01).
* **An unpredictable inflow face count**, set by refinement level and snapping
  rather than by a formula -- 5452 at the case's settings. Comparisons against the
  2044-face Pointwise mesh stop being like-for-like.
* **Approximate radii, in principle.** Snapped points land on the sphere to within
  the snapping tolerance, not exactly, so :mod:`plumetools.verify_mesh` loosens
  its radius check for this mesh type. In practice a *primitive* costs nothing:
  the worst measured vertex deviation is 8.3e-16 m, i.e. round-off. An STL would
  not behave this way.
* **A stair-step risk at the rim.** The sphere is bisected by ``x = 0``, so it
  crosses that boundary at right angles -- the well-conditioned case, not a
  tangency -- but the rim where ``inflow`` meets ``sym`` is still a surface /
  domain-boundary intersection, which is where snapping is least predictable.
  Inspect that circle before trusting the mesh.
"""

from __future__ import annotations

import math
from pathlib import Path

BANNER = """\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1706                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/"""

SEPARATOR = "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"


def _fmt(v: float) -> str:
    return f"{v:.10g}"


def _header(object_name: str, note: str) -> list[str]:
    return [
        BANNER,
        "// GENERATED by plumetools.foamio.snappy from case.yaml -- do not edit.",
        "// Regenerate with ./Allmesh. See docs/mesh.md.",
        f"// {note}",
        "FoamFile",
        "{",
        "    version     2.0;",
        "    format      ascii;",
        "    class       dictionary;",
        f"    object      {object_name};",
        "}",
        SEPARATOR,
        "",
    ]


def _write(path: Path, lines: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return path


def _geometry(cfg) -> tuple[float, float, float]:
    m = cfg.mesh
    R, H, L = float(m.sphere_radius_m), float(m.box_half_width_m), float(m.box_length_m)
    if not 0 < R < H:
        raise ValueError(f"need 0 < sphere_radius_m ({R}) < box_half_width_m ({H})")
    if L <= R:
        raise ValueError(f"need box_length_m ({L}) > sphere_radius_m ({R})")
    return R, H, L


def background_divisions(cfg) -> tuple[int, int, int]:
    """Background cell counts from ``mesh.background_cell_size_m``.

    Returns ``(nx, ny, nz)``, each at least 1. Cells are kept as close to cubic as
    rounding allows, which is what snappyHexMesh's octree refinement assumes.

    The background does not have to *resolve* the sphere -- ``refinement_level``
    does that -- but it must **find** it, so a background cell may not exceed the
    radius. At ``size == R`` the diameter still spans two cells.
    """
    R, H, L = _geometry(cfg)
    size = float(cfg.mesh.background_cell_size_m)
    if size <= 0:
        raise ValueError(f"background_cell_size_m must be positive, got {size}")
    if size > R:
        raise ValueError(
            f"background_cell_size_m ({size}) exceeds the sphere radius ({R}); the "
            f"cavity would fall between cells. Use at most R = {R}, and set the "
            f"surface resolution with refinement_level instead."
        )
    return (max(1, round(L / size)),
            max(1, round(2 * H / size)),
            max(1, round(2 * H / size)))


def surface_cell_size(cfg) -> float:
    """Cell size at the sphere: the background size halved once per octree level."""
    _R, _H, L = _geometry(cfg)
    nx, _ny, _nz = background_divisions(cfg)
    return (L / nx) / (2 ** int(cfg.mesh.refinement_level))


def location_in_mesh(cfg) -> tuple[float, float, float]:
    """A point in the fluid: inside the box, outside the sphere.

    A ``locationInMesh`` sitting exactly on a face, an edge, or a plane of
    symmetry is a classic way to make snappyHexMesh keep the wrong region or fail
    outright -- and "off-axis" is not sufficient on its own. An off-axis fraction
    of the domain still lands on a cell face whenever it happens to be a multiple
    of the background cell size, which ``0.5 * L`` is for every even ``nx``.

    So the point is snapped to the **centre** of the background cell containing
    it. A cell centre is interior by construction, whatever the cell size, and the
    seed then moves with the grid instead of drifting onto it.
    """
    R, H, L = _geometry(cfg)
    nx, ny, nz = background_divisions(cfg)

    def index_of(value: float, lo: float, span: float, n: int) -> int:
        return min(n - 1, max(0, int((value - lo) / (span / n))))

    def centre(index: int, lo: float, span: float, n: int) -> float:
        return lo + (index + 0.5) * (span / n)

    # Fractions chosen to sit well inside the box and off both symmetry planes.
    iy = index_of(0.3313 * H, -H, 2 * H, ny)
    iz = index_of(0.2371 * H, -H, 2 * H, nz)
    # On a coarse background the two fractions can round into the same cell, which
    # would put the seed on the y = z diagonal -- also a plane of symmetry here,
    # since the cross-section is square. One cell of separation is enough.
    if ny == nz and iy == iz:
        iz = iz - 1 if iz > 0 else iz + 1

    point = (centre(index_of(0.4137 * L, 0.0, L, nx), 0.0, L, nx),
             centre(iy, -H, 2 * H, ny),
             centre(iz, -H, 2 * H, nz))

    if math.dist(point, (0.0, 0.0, 0.0)) <= R:
        raise ValueError(
            f"locationInMesh {point} is inside the sphere (radius {R}); "
            f"snappyHexMesh would keep the cavity and discard the fluid"
        )
    return point


def render_background_block_mesh_dict(cfg) -> str:
    """A plain uniform box for snappyHexMesh to refine.

    Patches: ``sym`` on ``x = 0``, ``vacuum`` on the other five faces. The
    ``inflow`` patch does not exist yet -- snappyHexMesh creates it when it carves
    the cavity.

    The outer patch takes ``mesh.outer_patch_type``, exactly as the O-grid does.
    It is not cosmetic: standard dsmcFoam's FreeStream injects on every
    ``patch``-type boundary, so leaving this as ``patch`` while the case runs the
    standard solver turns the vacuum boundary into a second inflow and aborts the
    run. See docs/solver-compatibility.md.
    """
    R, H, L = _geometry(cfg)
    nx, ny, nz = background_divisions(cfg)
    names = cfg.mesh.patch_names
    p_outer = names.get("outer", "vacuum")
    p_sym = names.get("symmetry", "sym")
    t_outer = cfg.mesh.outer_patch_type

    # Box corners: 0-3 at x = 0, 4-7 at x = L, both wound (-y-z, +y-z, +y+z, -y+z).
    verts = [
        (0.0, -H, -H), (0.0, H, -H), (0.0, H, H), (0.0, -H, H),
        (L, -H, -H), (L, H, -H), (L, H, H), (L, -H, H),
    ]

    lines = _header(
        "blockMeshDict",
        f"Background box for snappyHexMesh: {nx}x{ny}x{nz} = {nx * ny * nz} cells, "
        f"cell size ~{_fmt(L / nx)} m. snappyHexMesh refines this down to "
        f"~{_fmt(surface_cell_size(cfg))} m at the sphere.",
    )
    lines += ["scale   1;", "", "vertices", "("]
    lines += [f"    ({_fmt(x)} {_fmt(y)} {_fmt(z)})" for x, y, z in verts]
    lines += [
        ");",
        "",
        "blocks",
        "(",
        f"    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1)",
        ");",
        "",
        "// No curved edges: the sphere is introduced by snappyHexMesh, not here.",
        "edges",
        "(",
        ");",
        "",
        "boundary",
        "(",
        f"    {p_sym}",
        "    {",
        "        type symmetry;",
        "        // x = 0. snappyHexMesh will cut the hemisphere out of this face.",
        "        faces",
        "        (",
        "            (0 3 2 1)",
        "        );",
        "    }",
        "",
        f"    {p_outer}",
        "    {",
        f"        type {t_outer};",
        "        // the remaining five box faces"
        + ("  --  WALL: reflects, does not absorb (see docs/solver-compatibility.md)"
           if t_outer == "wall" else ""),
        "        faces",
        "        (",
        "            (4 5 6 7)   // x = L",
        "            (0 1 5 4)   // z = -H",
        "            (3 7 6 2)   // z = +H",
        "            (0 4 7 3)   // y = -H",
        "            (1 2 6 5)   // y = +H",
        "        );",
        "    }",
        "",
        ");",
        "",
        "mergePatchPairs",
        "(",
        ");",
        "",
        "// ************************************************************************* //",
    ]
    return "\n".join(lines) + "\n"


def render_snappy_hex_mesh_dict(cfg) -> str:
    """The searchableSphere primitive plus castellate and snap controls."""
    R, _H, _L = _geometry(cfg)
    m = cfg.mesh
    level = int(m.refinement_level)
    if level < 0:
        raise ValueError(f"refinement_level must be >= 0, got {level}")
    inflow = m.patch_names.get("inflow", "inflow")
    lx, ly, lz = location_in_mesh(cfg)
    finest = surface_cell_size(cfg)

    lines = _header(
        "snappyHexMeshDict",
        f"Carve a hemispherical cavity of radius {_fmt(R)} m out of the background box. "
        f"Refinement level {level} -> ~{_fmt(finest)} m at the surface "
        f"({_fmt(R / finest)} cells across the radius).",
    )
    lines += [
        "castellatedMesh true;",
        "snap            true;",
        "addLayers       false;   // DSMC: no boundary-layer stack needed",
        "",
        "// The inflow surface as an analytic primitive -- no STL, so snapping targets",
        "// the true sphere rather than a triangulation of it. The geometry entry name",
        f"// becomes the patch name, which is why it is '{inflow}'.",
        "geometry",
        "{",
        f"    {inflow}",
        "    {",
        "        type    searchableSphere;",
        f"        centre  ({_fmt(0.0)} {_fmt(0.0)} {_fmt(0.0)});",
        f"        radius  {_fmt(R)};",
        "        // v1706 spells this 'centre'; OpenFOAM v2006+ also accepts 'origin'.",
        "    }",
        "}",
        "",
        "castellatedMeshControls",
        "{",
        "    maxLocalCells       2000000;",
        "    maxGlobalCells      4000000;",
        "    minRefinementCells  10;",
        "    maxLoadUnbalance    0.10;",
        f"    nCellsBetweenLevels {int(m.n_cells_between_levels)};",
        "",
        "    // No feature edges: a sphere has none.",
        "    features ();",
        "",
        "    refinementSurfaces",
        "    {",
        f"        {inflow}",
        "        {",
        f"            level ({level} {level});",
        "            patchInfo",
        "            {",
        "                type patch;   // never 'wall': this is an inflow, and it must",
        "                              // not become a no-slip surface by default",
        "            }",
        "        }",
        "    }",
        "",
        "    resolveFeatureAngle 30;",
        "",
        "    refinementRegions",
        "    {",
        "    }",
        "",
        "    // In the fluid: inside the box, outside the sphere, and deliberately off",
        "    // every symmetry plane. snappyHexMesh keeps the region connected to this",
        "    // point and discards the cavity interior.",
        f"    locationInMesh ({_fmt(lx)} {_fmt(ly)} {_fmt(lz)});",
        "",
        "    allowFreeStandingZoneFaces true;",
        "}",
        "",
        "snapControls",
        "{",
        "    nSmoothPatch    3;",
        "    tolerance       2.0;",
        "    nSolveIter      50;",
        "    nRelaxIter      5;",
        "    // Feature snapping is off: a sphere has no feature edges to capture. The",
        "    // one edge-like curve here is the rim where the sphere crosses x = 0, and",
        "    // that is a surface / domain-boundary intersection, which feature snapping",
        "    // tends to pull off the surface rather than onto it.",
        "    nFeatureSnapIter        10;",
        "    implicitFeatureSnap     false;",
        "    explicitFeatureSnap     false;",
        "    multiRegionFeatureSnap  false;",
        "}",
        "",
        "addLayersControls",
        "{",
        "    relativeSizes       true;",
        "    layers {}",
        "    expansionRatio      1.0;",
        "    finalLayerThickness 0.3;",
        "    minThickness        0.1;",
        "    nGrow               0;",
        "    featureAngle        60;",
        "    nRelaxIter          3;",
        "    nSmoothSurfaceNormals 1;",
        "    nSmoothNormals      3;",
        "    nSmoothThickness    10;",
        "    maxFaceThicknessRatio 0.5;",
        "    maxThicknessToMedialRatio 0.3;",
        "    minMedialAxisAngle  90;",
        "    nBufferCellsNoExtrude 0;",
        "    nLayerIter          50;",
        "}",
        "",
        "meshQualityControls",
        "{",
        '    #include "meshQualityDict"',
        "}",
        "",
        "debug 0;",
        "mergeTolerance 1e-6;",
        "",
        "// ************************************************************************* //",
    ]
    return "\n".join(lines) + "\n"


def render_mesh_quality_dict(cfg) -> str:
    """Quality limits, ``#include``d by snappyHexMeshDict.

    Kept as a local file rather than including OpenFOAM's
    ``etc/caseDicts/meshQualityDict`` so the case does not depend on an
    installation path.
    """
    lines = _header("meshQualityDict", "Quality limits for snappyHexMesh.")
    lines += [
        "maxNonOrtho         65;",
        "maxBoundarySkewness 20;",
        "maxInternalSkewness  4;",
        "maxConcave          80;",
        "minVol              1e-13;",
        "minTetQuality       1e-30;",
        "minArea             -1;",
        "minTwist            0.02;",
        "minDeterminant      0.001;",
        "minFaceWeight       0.02;",
        "minVolRatio         0.01;",
        "minTriangleTwist    -1;",
        "",
        "nSmoothScale 4;",
        "errorReduction 0.75;",
        "",
        "relaxed",
        "{",
        "    maxNonOrtho 75;",
        "}",
        "",
        "// ************************************************************************* //",
    ]
    return "\n".join(lines) + "\n"


def write_snappy_setup(case_dir: Path, cfg) -> list[Path]:
    """Write all three dictionaries for the snappyHexMesh path.

    Args:
        case_dir: an OpenFOAM case directory.
        cfg: a :class:`plumetools.config.CaseConfig` with
            ``mesh.type == "snappy_hex_sphere"``.

    Returns:
        The paths written, in the order blockMesh / snappyHexMesh consume them.
    """
    if cfg.mesh.type != "snappy_hex_sphere":
        raise ValueError(
            f"write_snappy_setup needs mesh.type 'snappy_hex_sphere', "
            f"got {cfg.mesh.type!r}"
        )
    case_dir = Path(case_dir)
    system = case_dir / "system"
    return [
        _write(system / "blockMeshDict", render_background_block_mesh_dict(cfg).splitlines()),
        _write(system / "snappyHexMeshDict", render_snappy_hex_mesh_dict(cfg).splitlines()),
        _write(system / "meshQualityDict", render_mesh_quality_dict(cfg).splitlines()),
    ]
