"""Capsule interface: watertight STL, feature extraction and snappyHexMesh inputs."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import math

import numpy as np

from plumetools.foamio.primitives import SearchableSphere, render_snappy_dict
from .ellipsoid_mesh import background, source_regions, finish_dicts

_FEATURES_PLACEHOLDER = "    features ();"


def _mesh(shape, n_phi, n_tip, back):
    """Vertices and triangles for given azimuthal and tip-ring counts."""
    R, L, a = shape.R, shape.L, shape.a
    theta = np.linspace(math.pi/2, 0.0, n_tip+1)[:-1]      # joint ring .. last ring before apex
    rings = [(-back, R), (L, R)] + [(L+a*math.cos(t), R*math.sin(t)) for t in theta[1:]]
    phi = np.arange(n_phi)*2*math.pi/n_phi
    verts = [(-back, 0.0, 0.0)]
    for x, r in rings:
        verts += [(x, r*math.cos(p), r*math.sin(p)) for p in phi]
    apex = len(verts)
    verts.append((L+a, 0.0, 0.0))
    ring = lambda i, j: 1+i*n_phi+j % n_phi
    tris = [(0, ring(0, j+1), ring(0, j)) for j in range(n_phi)]
    for i in range(len(rings)-1):
        for j in range(n_phi):
            a_, b_, c_, d_ = ring(i, j), ring(i, j+1), ring(i+1, j+1), ring(i+1, j)
            tris += [(a_, b_, c_), (a_, c_, d_)]
    last = len(rings)-1
    tris += [(ring(last, j), ring(last, j+1), apex) for j in range(n_phi)]
    return np.array(verts, float), np.array(tris, int)


def _bary(n=6):
    """Barycentric sample weights over a triangle, including vertices and edge midpoints."""
    w = [(i/n, j/n, (n-i-j)/n) for i in range(n+1) for j in range(n+1-i)]
    return np.array(w)


def _deviation(shape, vertices, triangles):
    """Max |sqrt(q) - 1| over barycentric samples of every triangle with x >= 0."""
    worst = 0.0
    weights = _bary()
    for start in range(0, len(triangles), 4096):
        pts = np.einsum("sk,tkd->tsd", weights, vertices[triangles[start:start+4096]]).reshape(-1, 3)
        pts = pts[pts[:, 0] >= 0]
        if len(pts):
            worst = max(worst, float(np.max(np.abs(np.sqrt(shape.q(pts))-1))))
    return worst


def capsule_stl(shape, chord_tol: float, back: float):
    """Closed outward triangulation: cylinder from x = -back to L, spheroidal tip, back disc.

    The ring counts are refined until the sagitta relative to the local radius
    (max |sqrt(q) - 1| over a barycentric sample grid of each triangle, x >= 0) is at most ``chord_tol``.
    """
    if not chord_tol > 0:
        raise ValueError("chord_tol must be positive")
    n_phi = max(8, math.ceil(math.pi/math.acos(1-min(chord_tol, .5))))
    n_tip = max(4, math.ceil(math.pi/2/math.sqrt(8*chord_tol)))
    for _ in range(40):
        vertices, triangles = _mesh(shape, n_phi, n_tip, float(back))
        if _deviation(shape, vertices, triangles) <= chord_tol:
            return vertices, triangles
        n_phi, n_tip = math.ceil(n_phi*1.2), math.ceil(n_tip*1.2)
    raise RuntimeError("capsule STL did not reach the chord tolerance")


def write_stl(path, vertices, triangles):
    """ASCII STL, solid ``interface``."""
    v = np.asarray(vertices, float)
    lines = ["solid interface"]
    for tri in np.asarray(triangles, int):
        p = v[tri]
        n = np.cross(p[1]-p[0], p[2]-p[0])
        norm = np.linalg.norm(n)
        n = n/norm if norm > 0 else n
        lines.append(f"  facet normal {n[0]:.9g} {n[1]:.9g} {n[2]:.9g}")
        lines.append("    outer loop")
        lines += [f"      vertex {q[0]:.12g} {q[1]:.12g} {q[2]:.12g}" for q in p]
        lines += ["    endloop", "  endfacet"]
    lines.append("endsolid interface")
    Path(path).write_text("\n".join(lines)+"\n", encoding="utf-8", newline="\n")


@dataclass(frozen=True)
class _TriSurface(SearchableSphere):
    def geometry_lines(self):
        return ["        type triSurfaceMesh;", f"        name {self.name};"]


def write_mesh_inputs(case_dir, cfg, source, shape):
    """Mesher inputs for a fitted capsule; returns the written dictionary paths."""
    case_dir = Path(case_dir)
    system = case_dir/"mesh-build/system"
    tri = case_dir/"mesh-build/constant/triSurface"
    system.mkdir(parents=True, exist_ok=True)
    tri.mkdir(parents=True, exist_ok=True)
    box = background(cfg, source)
    if shape.q(np.array([box["seed"]]))[0] <= 1:
        raise ValueError("background seed is inside the fitted capsule")
    level = int(cfg.get("dsmc.mesh.surface_level", 3))
    tol = .1*float(cfg.get("interface.fit.geometry_relative_tolerance", .02))
    vertices, triangles = capsule_stl(shape, tol, box["size"])
    write_stl(tri/"interface.stl", vertices, triangles)
    surface = _TriSurface(name="interface", level=level, radius=1.0, centre=(0, 0, 0),
                          patch_type="patch")
    regions = source_regions(cfg, source, box["size"], shape.R)
    snappy = render_snappy_dict([surface], location_in_mesh=box["seed"], n_cells_between_levels=3,
                                generator=__name__, note="fitted 3D capsule and source refinement",
                                max_global_cells=int(cfg.get("dsmc.mesh.max_global_cells", 12000000)),
                                refinement_regions=regions)
    # The geometry entry is keyed by the STL file; refinement is keyed by `name`.
    snappy = snappy.replace("    interface\n    {\n        type triSurfaceMesh;",
                            "    interface.stl\n    {\n        type triSurfaceMesh;", 1)
    if _FEATURES_PLACEHOLDER not in snappy:
        raise RuntimeError("snappyHexMeshDict features placeholder not found")
    snappy = snappy.replace(_FEATURES_PLACEHOLDER, "    features\n    (\n        {\n"
                            f'            file "interface.eMesh";\n            level {level};\n'
                            "        }\n    );", 1)
    paths = finish_dicts(system, cfg, source, box["block"], snappy)
    extract = system/"surfaceFeatureExtractDict"
    extract.write_text("FoamFile { version 2.0; format ascii; class dictionary; "
                       "object surfaceFeatureExtractDict; }\n"
                       "interface.stl\n{\n    extractionMethod extractFromSurface;\n"
                       "    includedAngle 150;\n    writeObj no;\n}\n", encoding="utf-8", newline="\n")
    return paths+[extract]
