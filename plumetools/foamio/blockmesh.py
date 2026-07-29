r"""Generate a ``blockMeshDict`` for a box with a hemispherical inflow cavity.

Replaces the binary Pointwise dependency for ``cases/3d-inflow`` (finding RP-16).
The geometry is a pure function of ``(R, H, L, n_tangential, n_radial, grading)``,
so the dictionary is text, diffable, and regenerable from ``case.yaml``.

Geometry
--------
A box ``x in [0, L]``, ``y, z in [-H, H]``, with a hemisphere of radius ``R``
removed at the origin on the ``x = 0`` face::

        z
        ^          vacuum (5 outer box faces)
        |   +---------------------------+
        |   |                           |
     ---+---|-- ) inflow                |----> x
        |   |  (hemisphere, radius R)   |
        |   +---------------------------+
            ^ sym (the x = 0 plane)

Topology -- a 5-block O-grid
----------------------------
Take the ``+x`` half of a cube inscribed in the sphere and connect each of its
five outward faces to the corresponding box face. With ``a = R/sqrt(3)`` (a cube
corner projected onto the sphere) and ``b = R/sqrt(2)`` (an equator point at 45
degrees), the blocks are: a *cap* block reaching the ``x = L`` face, and four
*side* blocks reaching ``y = +-H`` and ``z = +-H``. Each side block contributes
one face to ``sym``; the cap block does not touch ``x = 0``.

The five inner faces tile the hemisphere exactly: the ``+x`` half of a cube's six
faces is the ``+x`` face whole, half of each of ``+-y`` and ``+-z``, and no ``-x``
face.

Curvature comes from ``arc`` edges on the twelve inner edges, with interpolation
points computed analytically as ``R * (P1 + P2) / |P1 + P2|``. ``arc`` is used
rather than ``project``/``searchableSphere`` because the projection directive
syntax could not be verified against OpenFOAM v1706 -- see ``docs/mesh.md``. The
consequence is that block *faces* on the sphere are ruled surfaces, so interior
surface points bulge very slightly inside the true sphere.

Block and boundary-face orientation are **derived**, not hand-written: blocks are
flipped if their inner-to-outer normal points the wrong way, and each boundary
face is oriented outward from the domain. Getting either backwards yields
negative-volume cells or inverted patches, and neither is obvious by inspection.

Consequence for the model
-------------------------
``blockMesh`` produces **hexahedra**, so inflow faces are **quadrilaterals**,
where the committed Pointwise mesh is tetrahedral with triangular faces. The
legacy centroid divided the vertex sum by the literal ``3.0``, which on a quad
gives exactly 4/3 of the true centroid -- the defect that silently corrupted the
archived wake-cylinder cases (AR-01). :func:`plumetools.geometry.centroid`
divides by ``len(vertices)``, so this mesh is safe. Do not reintroduce the
``3.0``.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

_SQRT3 = math.sqrt(3.0)
_SQRT2 = math.sqrt(2.0)

#: Sign pattern for (y, z) within each group of four vertices.
_SIGNS = ((+1, +1), (+1, -1), (-1, -1), (-1, +1))

#: (name, inner four vertex indices, outer four vertex indices), before orientation.
BLOCKS = (
    ("cap", (0, 1, 2, 3), (8, 9, 10, 11)),
    ("+y", (1, 0, 4, 5), (9, 8, 12, 13)),
    ("-y", (3, 2, 6, 7), (11, 10, 14, 15)),
    ("+z", (0, 3, 7, 4), (8, 11, 15, 12)),
    ("-z", (2, 1, 5, 6), (10, 9, 13, 14)),
)

#: The twelve edges lying on the sphere: 4 cap, 4 meridian, 4 equator.
SPHERE_EDGES = (
    (0, 1), (1, 2), (2, 3), (3, 0),
    (0, 4), (1, 5), (2, 6), (3, 7),
    (4, 5), (5, 6), (6, 7), (7, 4),
)

#: The six faces of a hex, as indices into its local 0..7 vertex ordering.
_HEX_FACES = (
    (0, 1, 2, 3), (4, 5, 6, 7),
    (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7),
)

_TOL = 1e-9


def build_vertices(R: float, H: float, L: float) -> np.ndarray:
    """The 16 block-corner vertices.

    Args:
        R: sphere radius [m].
        H: box half-width in y and z [m].
        L: box length in x [m].

    Returns:
        ``(16, 3)`` array. Indices 0-3 are the inner cap corners (all at
        ``|r| = R``), 4-7 the inner equator ring (``x = 0``, ``|r| = R``), 8-11
        the outer ``x = L`` face, 12-15 the outer ``x = 0`` ring.
    """
    a = R / _SQRT3
    b = R / _SQRT2
    return np.array(
        [*[(a, sy * a, sz * a) for sy, sz in _SIGNS],
         *[(0.0, sy * b, sz * b) for sy, sz in _SIGNS],
         *[(L, sy * H, sz * H) for sy, sz in _SIGNS],
         *[(0.0, sy * H, sz * H) for sy, sz in _SIGNS]],
        dtype=np.float64,
    )


def arc_point(p1, p2, R: float) -> np.ndarray:
    """Midpoint of the great-circle arc between two points on a sphere.

    ``R * (p1 + p2) / |p1 + p2|`` -- exact, so the dictionary stays a pure
    function of the geometry parameters.
    """
    mid = np.asarray(p1, dtype=np.float64) + np.asarray(p2, dtype=np.float64)
    norm = np.linalg.norm(mid)
    if norm == 0.0:
        raise ValueError("antipodal points have no unique arc midpoint")
    return R * mid / norm


def _quad_normal(quad, verts: np.ndarray) -> np.ndarray:
    """Right-hand-rule normal of a quad given in winding order."""
    p = verts[list(quad)]
    return np.cross(p[1] - p[0], p[2] - p[1])


def orient_block(inner, outer, verts: np.ndarray):
    """Order a hex so its volume is positive.

    OpenFOAM requires the first four vertices to be wound so the right-hand
    normal points toward the second four. Reverses both quads together if not.
    """
    inner_centre = verts[list(inner)].mean(axis=0)
    outer_centre = verts[list(outer)].mean(axis=0)
    if np.dot(_quad_normal(inner, verts), outer_centre - inner_centre) < 0:
        return tuple(reversed(inner)), tuple(reversed(outer))
    return tuple(inner), tuple(outer)


def orient_face_outward(quad, verts: np.ndarray, outward: np.ndarray):
    """Wind a boundary face so its normal points out of the domain."""
    if np.dot(_quad_normal(quad, verts), outward) < 0:
        return tuple(reversed(quad))
    return tuple(quad)


def classify_boundary_faces(verts: np.ndarray, R: float, H: float, L: float):
    """Sort every block face into inflow / sym / outer / internal, geometrically.

    Classifying by position rather than by hard-coded index lists means a change
    to the block table cannot silently mis-assign a patch.

    Returns:
        ``(inflow, sym, outer)`` lists of outward-wound global vertex quads.
    """
    inflow, sym, outer = [], [], []

    for _, inner_raw, outer_raw in BLOCKS:
        inner, outer_q = orient_block(inner_raw, outer_raw, verts)
        hexa = (*inner, *outer_q)

        for local in _HEX_FACES:
            quad = tuple(hexa[i] for i in local)
            p = verts[list(quad)]
            centre = p.mean(axis=0)

            if np.all(np.abs(np.linalg.norm(p, axis=1) - R) < _TOL):
                # On the sphere: outward from the domain means into the cavity.
                inflow.append(orient_face_outward(quad, verts, -centre))
            elif np.all(np.abs(p[:, 0]) < _TOL):
                sym.append(orient_face_outward(quad, verts, np.array([-1.0, 0.0, 0.0])))
            elif np.all(np.abs(p[:, 0] - L) < _TOL):
                outer.append(orient_face_outward(quad, verts, np.array([1.0, 0.0, 0.0])))
            else:
                for axis in (1, 2):
                    for sign in (+1, -1):
                        if np.all(np.abs(p[:, axis] - sign * H) < _TOL):
                            hint = np.zeros(3)
                            hint[axis] = sign
                            outer.append(orient_face_outward(quad, verts, hint))

    return inflow, sym, outer


def _fmt(v: float) -> str:
    return f"{v:.10g}"


def _vec(p) -> str:
    return f"({_fmt(p[0])} {_fmt(p[1])} {_fmt(p[2])})"


def _patch_block(name: str, patch_type: str, faces, comment: str) -> list[str]:
    out = [f"    {name}", "    {", f"        type {patch_type};",
           f"        // {comment}", "        faces", "        ("]
    out += [f"            ({' '.join(str(i) for i in f)})" for f in faces]
    out += ["        );", "    }", ""]
    return out


def render_block_mesh_dict(cfg) -> str:
    """Render the dictionary text.

    Args:
        cfg: a :class:`plumetools.config.CaseConfig`.

    Returns:
        The complete ``blockMeshDict`` as a string.

    Raises:
        NotImplementedError: for a ``mesh.type`` other than ``block_mesh_ogrid``
            or a ``mesh.projection`` other than ``none``.
        ValueError: if the geometry parameters are not physically ordered.
    """
    mesh = cfg.mesh
    if mesh.type != "block_mesh_ogrid":
        raise NotImplementedError(f"mesh.type {mesh.type!r} is not implemented")
    if mesh.projection != "none":
        raise NotImplementedError(
            f"mesh.projection {mesh.projection!r} is not implemented. Only 'none' "
            f"(arc edges) is supported: the project/searchableSphere directive "
            f"syntax has not been verified against OpenFOAM v1706. See docs/mesh.md."
        )

    R, H, L = float(mesh.sphere_radius_m), float(mesh.box_half_width_m), float(mesh.box_length_m)
    if not 0 < R < H:
        raise ValueError(f"need 0 < sphere_radius_m ({R}) < box_half_width_m ({H})")
    if L <= R:
        raise ValueError(f"need box_length_m ({L}) > sphere_radius_m ({R})")
    nt, nr = int(mesh.n_tangential), int(mesh.n_radial)
    if nt < 1 or nr < 1:
        raise ValueError(f"n_tangential ({nt}) and n_radial ({nr}) must be >= 1")

    names = mesh.patch_names
    p_inflow = names.get("inflow", "inflow")
    p_outer = names.get("outer", "vacuum")
    p_sym = names.get("symmetry", "sym")

    verts = build_vertices(R, H, L)
    inflow, sym, outer = classify_boundary_faces(verts, R, H, L)

    lines = [
        "/*--------------------------------*- C++ -*----------------------------------*\\",
        "| =========                 |                                                 |",
        "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |",
        "|  \\\\    /   O peration     | Version:  v1706                                 |",
        "|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |",
        "|    \\\\/     M anipulation  |                                                 |",
        "\\*---------------------------------------------------------------------------*/",
        "// GENERATED by plumetools.foamio.blockmesh from case.yaml -- do not edit.",
        "// Regenerate with ./Allmesh. See docs/mesh.md.",
        "//",
        f"// Box x in [0, {_fmt(L)}], y,z in [{_fmt(-H)}, {_fmt(H)}], with a hemispherical",
        f"// inflow cavity of radius {_fmt(R)} at the origin on the x = 0 symmetry plane.",
        f"// 5-block O-grid -> {5 * nt * nt} inflow quads, {5 * nt * nt * nr} hexahedral cells.",
        "FoamFile",
        "{",
        "    version     2.0;",
        "    format      ascii;",
        "    class       dictionary;",
        "    object      blockMeshDict;",
        "}",
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //",
        "",
        "scale   1;",
        "",
        "vertices",
        "(",
    ]

    labels = (["inner cap"] * 4 + ["inner equator"] * 4
              + [f"outer x = {_fmt(L)}"] * 4 + ["outer x = 0"] * 4)
    for i, (v, label) in enumerate(zip(verts, labels)):
        lines.append(f"    {_vec(v):<46} // {i:2d}  {label}")
    lines += [");", "", "blocks", "("]

    grading = _fmt(mesh.radial_grading)
    for name, inner_raw, outer_raw in BLOCKS:
        inner, outer_q = orient_block(inner_raw, outer_raw, verts)
        order = " ".join(str(i) for i in (*inner, *outer_q))
        lines.append(
            f"    hex ({order}) ({nt} {nt} {nr}) simpleGrading (1 1 {grading})"
            f"   // {name}"
        )
    lines += [");", "",
              "// Sphere curvature. Without these the inflow surface would be flat.",
              "edges", "("]
    for v0, v1 in SPHERE_EDGES:
        lines.append(f"    arc {v0} {v1} {_vec(arc_point(verts[v0], verts[v1], R))}")
    lines += [");", "", "boundary", "("]

    lines += _patch_block(p_inflow, "patch", inflow,
                          "hemispherical plume source surface")
    lines += _patch_block(p_outer, "patch", outer,
                          "the five outer box faces")
    lines += _patch_block(p_sym, "symmetry", sym,
                          "the x = 0 plane")

    lines += [");", "", "mergePatchPairs", "(", ");", "",
              "// ************************************************************************* //"]
    return "\n".join(lines) + "\n"


def write_block_mesh_dict(case_dir: Path, cfg) -> Path:
    """Write ``<case_dir>/system/blockMeshDict``.

    Args:
        case_dir: an OpenFOAM case directory.
        cfg: a :class:`plumetools.config.CaseConfig`.

    Returns:
        The path written.

    Writing the dictionary is harmless. Running ``blockMesh`` afterwards is not --
    it overwrites ``constant/polyMesh``, which for ``cases/3d-inflow`` is the only
    copy of the Pointwise mesh and the substrate of the regression golden. See
    ``Allmesh`` and ``docs/mesh.md``.
    """
    case_dir = Path(case_dir)
    out = case_dir / "system" / "blockMeshDict"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write(render_block_mesh_dict(cfg))
    return out
