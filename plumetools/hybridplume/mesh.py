"""The full 3-D Cartesian mesh shared by the CFD precursor and the DSMC region.

One ``blockMesh`` box ``0 <= x <= X``, ``|y|, |z| <= W``: a uniform core around
the source and the near plume, and geometrically expanding segments outside it
(the expansion ratio is solved so the first outer cell matches the core cell;
:mod:`plumetools.cai2012.mesh` owns that solver and is reused unchanged). No
wedge, no symmetry plane, no empty patch.

Patches
-------
=============  =========  ===================================================
name           type       role
=============  =========  ===================================================
``source``     patch      the circular exit, ``x = 0, r <= Re``: CFD inflow;
                          DSMC injection (and absorption of returning molecules)
``exitPlane``  patch      ``x = 0, r > Re``: open in both solvers. CFD outflow
                          into vacuum; DSMC deletes molecules crossing it
``farField``   patch      the other five faces: open in both solvers
=============  =========  ===================================================

The same physical boundary condition outside the disk -- gas crossing ``x = 0``
leaves the problem -- is therefore imposed in both solvers, which is what makes
the two solutions comparable near the lip.

The source disk on a Cartesian mesh
-----------------------------------
``source`` is the set of ``x = 0`` faces whose centres satisfy ``y^2 + z^2 <=
Re^2`` (a staircase), selected by ``topoSet`` and promoted by ``createPatch``.
:func:`predicted_source` applies the same rule in Python before the mesh exists
so the area error -- a direct multiplier on the injected flux -- is known and
reported. It is never rescaled away: a flux rescaling would hide a geometric
error behind a physical-looking density change.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from plumetools.cai2012.mesh import Segment, geometric_segment, uniform_segment

from .config import CaseConfig
from .physics import SourceState

GENERATOR = "plumetools.hybridplume.mesh"
PATCHES = {"source": "source", "exit_plane": "exitPlane", "far_field": "farField"}
SELECTOR_HALF_THICKNESS_FRACTION = 0.25

BOX_FACES = {
    "x_min": (0, 4, 7, 3), "x_max": (1, 2, 6, 5),
    "y_min": (0, 1, 5, 4), "y_max": (3, 7, 6, 2),
    "z_min": (0, 3, 2, 1), "z_max": (4, 5, 6, 7),
}


@dataclass(frozen=True)
class MeshPlan:
    x_segments: tuple
    y_segments: tuple
    z_segments: tuple
    core_cell_m: float

    @property
    def divisions(self) -> tuple:
        return tuple(sum(s.n_cells for s in segs)
                     for segs in (self.x_segments, self.y_segments, self.z_segments))

    @property
    def n_cells(self) -> int:
        nx, ny, nz = self.divisions
        return nx * ny * nz

    @property
    def n_blocks(self) -> int:
        return len(self.x_segments) * len(self.y_segments) * len(self.z_segments)

    def edges(self, axis: int) -> np.ndarray:
        segs = (self.x_segments, self.y_segments, self.z_segments)[axis]
        parts = [segs[0].cell_edges()] + [s.cell_edges()[1:] for s in segs[1:]]
        return np.concatenate(parts)

    @property
    def max_expansion_ratio(self) -> float:
        ratios = []
        for segs in (self.x_segments, self.y_segments, self.z_segments):
            widths = np.diff(np.concatenate(
                [segs[0].cell_edges()] + [s.cell_edges()[1:] for s in segs[1:]]))
            ratios.append(float(np.max(widths[1:] / widths[:-1])))
            ratios.append(float(np.max(widths[:-1] / widths[1:])))
        return max(ratios)

    def describe(self) -> dict:
        return {
            "divisions": list(self.divisions), "n_cells": self.n_cells,
            "n_blocks": self.n_blocks, "core_cell_m": self.core_cell_m,
            "max_adjacent_cell_ratio": self.max_expansion_ratio,
            "x_range_m": [self.x_segments[0].lo, self.x_segments[-1].hi],
            "y_range_m": [self.y_segments[0].lo, self.y_segments[-1].hi],
            "z_range_m": [self.z_segments[0].lo, self.z_segments[-1].hi],
        }


def _outer(lo: float, hi: float, n: int, first: float, cap: float, reverse: bool) -> Segment:
    seg = geometric_segment(0.0, hi - lo, n, first, cap)
    if not reverse:
        return Segment(lo, hi, seg.n_cells, seg.expansion, seg.first_cell_m, seg.last_cell_m)
    # Expansion runs outward, i.e. towards lo: flip the grading.
    return Segment(lo, hi, seg.n_cells, 1.0 / seg.expansion, seg.last_cell_m, seg.first_cell_m)


def plan(cfg: CaseConfig, source: SourceState) -> MeshPlan:
    d = source.diameter_m
    h = float(cfg.require("mesh.core_cell_over_D")) * d
    x_max = float(cfg.require("domain.x_max_over_D")) * d
    w = float(cfg.require("domain.half_width_over_D")) * d
    core_x = float(cfg.get("mesh.core_x_over_D", cfg.require("domain.x_max_over_D"))) * d
    core_w = float(cfg.get("mesh.core_half_over_D",
                           cfg.require("domain.half_width_over_D"))) * d
    n_out = int(cfg.get("mesh.outer_cells", 0))
    cap = float(cfg.get("mesh.max_outer_expansion", 8.0))
    if core_w <= source.radius_m:
        raise ValueError("the uniform core must contain the whole source disk")

    # Snap the core so a whole number of cells spans it and the disk sits on a
    # cell-centred lattice symmetric about the axis.
    nc_half = max(1, int(math.ceil(core_w / h - 1e-9)))
    core_w = nc_half * h
    xs_core = uniform_segment(0.0, min(core_x, x_max), h)
    xs = [xs_core]
    if x_max > xs_core.hi + 1e-12:
        if n_out < 1:
            raise ValueError("domain.x_max_over_D exceeds mesh.core_x_over_D but "
                             "mesh.outer_cells is 0")
        xs.append(_outer(xs_core.hi, x_max, n_out, xs_core.last_cell_m, cap, False))
    ys = [uniform_segment(-core_w, core_w, h)]
    if w > core_w + 1e-12:
        if n_out < 1:
            raise ValueError("domain.half_width_over_D exceeds the core but "
                             "mesh.outer_cells is 0")
        ys = [_outer(-w, -core_w, n_out, h, cap, True), ys[0],
              _outer(core_w, w, n_out, h, cap, False)]
    return MeshPlan(tuple(xs), tuple(ys), tuple(ys), h)


@dataclass(frozen=True)
class SourcePrediction:
    n_faces: int
    area_m2: float
    exact_area_m2: float
    face_centres: np.ndarray       # (n, 3)

    @property
    def area_error(self) -> float:
        return self.area_m2 / self.exact_area_m2 - 1.0

    def describe(self) -> dict:
        return {"n_faces": self.n_faces, "area_m2": self.area_m2,
                "exact_area_m2": self.exact_area_m2,
                "relative_area_error": self.area_error,
                "note": "the injected and continuum source fluxes scale with this "
                        "area; it is reported, never rescaled"}


def predicted_source(p: MeshPlan, source: SourceState) -> SourcePrediction:
    """The staircase source faces the mesh pipeline will select."""
    ye = p.edges(1)
    ze = p.edges(2)
    yc = 0.5 * (ye[:-1] + ye[1:])
    zc = 0.5 * (ze[:-1] + ze[1:])
    yy, zz = np.meshgrid(yc, zc, indexing="ij")
    dy = np.diff(ye)[:, None] * np.ones_like(zz)
    dz = np.ones_like(yy) * np.diff(ze)[None, :]
    inside = yy ** 2 + zz ** 2 <= source.radius_m ** 2
    centres = np.stack([np.zeros(inside.sum()), yy[inside], zz[inside]], axis=1)
    return SourcePrediction(int(inside.sum()), float((dy * dz)[inside].sum()),
                            source.area_m2, centres)


def _fmt(v: float) -> str:
    return f"{float(v):.10g}"


def _header(obj: str, note: str) -> list:
    return [
        "/*--------------------------------*- C++ -*----------------------------------*\\",
        f"  Generated by {GENERATOR}. Do not edit: regenerate from hybrid.yaml.",
        f"  {note}",
        "\\*---------------------------------------------------------------------------*/",
        "FoamFile", "{", "    version     2.0;", "    format      ascii;",
        "    class       dictionary;", f"    object      {obj};", "}", "",
    ]


def render_block_mesh_dict(p: MeshPlan) -> str:
    xs = [p.x_segments[0].lo] + [s.hi for s in p.x_segments]
    ys = [p.y_segments[0].lo] + [s.hi for s in p.y_segments]
    zs = [p.z_segments[0].lo] + [s.hi for s in p.z_segments]

    def index(i, j, k):
        return i + j * len(xs) + k * len(xs) * len(ys)

    nx, ny, nz = p.divisions
    lines = _header("blockMeshDict",
                    f"{p.n_blocks} blocks, {nx}x{ny}x{nz} = {p.n_cells} cells, "
                    f"core cell {_fmt(p.core_cell_m)} m. Full 3-D: no wedge, no symmetry.")
    lines += ["scale 1;", "", "vertices", "("]
    for z in zs:
        for y in ys:
            for x in xs:
                lines.append(f"    ({_fmt(x)} {_fmt(y)} {_fmt(z)})")
    lines += [");", "", "blocks", "("]
    for k, zseg in enumerate(p.z_segments):
        for j, yseg in enumerate(p.y_segments):
            for i, xseg in enumerate(p.x_segments):
                c = (index(i, j, k), index(i + 1, j, k), index(i + 1, j + 1, k),
                     index(i, j + 1, k), index(i, j, k + 1), index(i + 1, j, k + 1),
                     index(i + 1, j + 1, k + 1), index(i, j + 1, k + 1))
                lines.append(
                    "    hex (" + " ".join(map(str, c)) + f") ({xseg.n_cells} "
                    f"{yseg.n_cells} {zseg.n_cells}) simpleGrading ({_fmt(xseg.expansion)} "
                    f"{_fmt(yseg.expansion)} {_fmt(zseg.expansion)})")
    lines += [");", "", "edges ();", "", "boundary", "("]
    faces = {"exit_plane": [], "far_field": []}
    ni, nj, nk = len(p.x_segments), len(p.y_segments), len(p.z_segments)
    for k in range(nk):
        for j in range(nj):
            for i in range(ni):
                local = {0: index(i, j, k), 1: index(i + 1, j, k),
                         2: index(i + 1, j + 1, k), 3: index(i, j + 1, k),
                         4: index(i, j, k + 1), 5: index(i + 1, j, k + 1),
                         6: index(i + 1, j + 1, k + 1), 7: index(i, j + 1, k + 1)}

                def quad(key):
                    return tuple(local[v] for v in BOX_FACES[key])
                if i == 0:
                    faces["exit_plane"].append(quad("x_min"))
                if i == ni - 1:
                    faces["far_field"].append(quad("x_max"))
                if j == 0:
                    faces["far_field"].append(quad("y_min"))
                if j == nj - 1:
                    faces["far_field"].append(quad("y_max"))
                if k == 0:
                    faces["far_field"].append(quad("z_min"))
                if k == nk - 1:
                    faces["far_field"].append(quad("z_max"))
    notes = {"exit_plane": "x = 0; the source disk is carved out of it by topoSet",
             "far_field": "open vacuum boundary"}
    for role in ("exit_plane", "far_field"):
        lines += [f"    {PATCHES[role]}", "    {", "        type patch;",
                  f"        // {notes[role]}", "        faces", "        ("]
        lines += [f"            ({a} {b} {c} {d})" for a, b, c, d in faces[role]]
        lines += ["        );", "    }"]
    lines += [");", "", "mergePatchPairs ();", ""]
    return "\n".join(lines) + "\n"


def render_topo_set_dict(p: MeshPlan, source: SourceState) -> str:
    half = SELECTOR_HALF_THICKNESS_FRACTION * p.core_cell_m
    lines = _header("topoSetDict", f"x = 0 faces with centre radius <= {_fmt(source.radius_m)} m")
    lines += [
        "// patchToFace first, then cylinderToFace subset: reversed, the cylinder",
        "// would also select interior faces.",
        "actions", "(",
        "    { name sourceFaces; type faceSet; action new; source patchToFace;"
        f" patch {PATCHES['exit_plane']}; }}",
        "    { name sourceFaces; type faceSet; action subset; source cylinderToFace;",
        f"      point1 ({_fmt(-half)} 0 0); point2 ({_fmt(half)} 0 0);"
        f" radius {_fmt(source.radius_m)}; }}",
        ");", ""]
    return "\n".join(lines) + "\n"


def render_create_patch_dict() -> str:
    lines = _header("createPatchDict", "promote the source disk to its own patch")
    lines += [
        "pointSync false;", "", "patches", "(", "    {",
        f"        name {PATCHES['source']};",
        "        // 'patch', never 'wall': plumeFieldInflow refuses a wall inflow.",
        "        patchInfo { type patch; }",
        "        constructFrom set;", "        set sourceFaces;", "    }", ");", ""]
    return "\n".join(lines) + "\n"


def write_mesh_setup(case_dir: Path, p: MeshPlan, source: SourceState) -> list:
    system = Path(case_dir) / "system"
    system.mkdir(parents=True, exist_ok=True)
    out = []
    for name, text in (("blockMeshDict", render_block_mesh_dict(p)),
                       ("topoSetDict", render_topo_set_dict(p, source)),
                       ("createPatchDict", render_create_patch_dict()),
                       ("meshQualityDict",
                        '#includeEtc "caseDicts/meshQualityDict"\n')):
        path = system / name
        if name == "meshQualityDict":
            text = "\n".join(_header("meshQualityDict", "OpenFOAM defaults")) + "\n" + text
        path.write_text(text)
        out.append(path)
    return out


MESH_COMMANDS = (
    ("blockMesh", ()),
    ("topoSet", ()),
    ("createPatch", ("-overwrite",)),
    ("checkMesh", ()),
)
