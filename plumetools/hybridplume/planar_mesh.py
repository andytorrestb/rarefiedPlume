"""Small native blockMesh construction outside a clipped elliptical core."""
from __future__ import annotations

import math
import numpy as np
from .mesh import _header


def render_dsmc_mesh(cfg, source, fit):
    d = source.diameter_m
    depth = source.depth_m
    h = float(cfg.require("mesh.core_cell_over_D")) * d
    xmax = float(cfg.require("domain.x_max_over_D")) * d
    width = float(cfg.require("domain.half_width_over_D")) * d
    xc, a, b = (float(fit[k]) for k in ("xc", "a", "b"))
    if depth is None or not (0 <= xc < a and b > 0 and xc+a < xmax and b < width):
        raise ValueError("planar ellipse must cross the source and stay inside the domain")
    theta0 = math.acos(-xc/a)
    theta = np.linspace(theta0, -theta0, 7)
    inner = np.column_stack((xc + a*np.cos(theta), b*np.sin(theta)))
    inner[[0, -1], 0] = 0
    slit_half = float(inner[0, 1])
    if slit_half > source.radius_m + 1e-9*d:
        raise ValueError("ellipse source intersection extends outside the slit")
    outer = np.array([[0, width], [xmax/2, width], [xmax, width], [xmax, 0],
                      [xmax, -width], [xmax/2, -width], [0, -width]])
    fraction = (source.radius_m-slit_half)/(width-slit_half)
    split = fraction > 1e-8
    fractions = [0, float(fraction), 1] if split else [0, 1]
    layers = [(1-f)*inner + f*outer for f in fractions]
    n = len(inner)
    nl = len(layers)

    def index(k, layer, i):
        return k*nl*n + layer*n + i

    def fmt(p):
        return "(" + " ".join(f"{v:.14g}" for v in p) + ")"

    lines = _header("blockMeshDict", "planar exterior of conservative fitted ellipse")
    # Plane sections of this compensated ellipsoid have precisely axes a and b.
    compensation = math.sqrt(1 + (depth/(2*b))**2)
    radii = (a*compensation, b*compensation, b*compensation)
    lines += ["scale 1;", "geometry", "{", "    fittedEllipse",
              "    {", "        type searchableSphere;",
              f"        centre ({xc:.14g} 0 0);", f"        radius {fmt(radii)};", "    }"]
    for name, z in (("frontPlane", -depth/2), ("backPlane", depth/2)):
        lines += [f"    {name}", "    {", "        type searchablePlane;",
                  "        planeType pointAndNormal;", "        pointAndNormalDict",
                  "        {", f"            point (0 0 {z:.14g});",
                  "            normal (0 0 1);", "        }", "    }"]
    lines += ["}", "vertices", "("]
    for z in (-depth/2, depth/2):
        for layer in layers:
            lines += ["    " + fmt((x, y, z)) for x, y in layer]
    lines += [");", "blocks", "("]
    patches = {name: [] for name in ("interface", "source", "exitPlane", "farField", "front", "back")}
    for layer in range(nl-1):
        nr = max(2, int(math.ceil(np.linalg.norm(layers[layer+1]-layers[layer], axis=1).max()/h)))
        for i in range(n-1):
            c = (index(0, layer, i), index(0, layer, i+1),
                 index(0, layer+1, i+1), index(0, layer+1, i),
                 index(1, layer, i), index(1, layer, i+1),
                 index(1, layer+1, i+1), index(1, layer+1, i))
            ns = max(4, int(math.ceil(max(np.linalg.norm(outer[i+1]-outer[i]),
                                         a*abs(theta[i+1]-theta[i]))/h)))
            lines += ["    hex (" + " ".join(map(str, c)) + f") ({ns} {nr} 1) simpleGrading (1 1 1)"]
            quad = lambda q: tuple(c[j] for j in q)
            if layer == 0:
                patches["interface"].append(quad((0, 1, 5, 4)))
            if layer == nl-2:
                patches["farField"].append(quad((3, 7, 6, 2)))
            if i == 0:
                patches["source" if split and layer == 0 else "exitPlane"].append(quad((0, 4, 7, 3)))
            if i == n-2:
                patches["source" if split and layer == 0 else "exitPlane"].append(quad((1, 2, 6, 5)))
            patches["front"].append(quad((0, 3, 2, 1)))
            patches["back"].append(quad((4, 5, 6, 7)))
    lines += [");", "edges", "("]
    for k, plane in ((0, "frontPlane"), (1, "backPlane")):
        z = (-0.5 if k == 0 else 0.5)*depth
        for i in range(n-1):
            lines += [f"    project {index(k, 0, i)} {index(k, 0, i+1)} (fittedEllipse {plane})"]
            if split:
                f = fractions[1]
                angle = np.linspace(theta[i], theta[i+1], 17)[1:-1]
                arc = np.column_stack((xc+a*np.cos(angle), b*np.sin(angle)))
                t = np.linspace(0, 1, 17)[1:-1, None]
                outer_points = (1-t)*outer[i] + t*outer[i+1]
                controls = (1-f)*arc + f*outer_points
                lines += [f"    spline {index(k, 1, i)} {index(k, 1, i+1)}",
                          "    ("]
                lines += ["        "+fmt((x, y, z)) for x, y in controls]
                lines += ["    )"]
    lines += [");", "boundary", "("]
    for name, faces in patches.items():
        kind = "empty" if name in ("front", "back") else "patch"
        lines += [f"    {name}", "    {", f"        type {kind};", "        faces", "        ("]
        lines += ["            ("+" ".join(map(str, face))+")" for face in faces]
        lines += ["        );", "    }"]
    lines += [");", "mergePatchPairs ();", ""]
    return "\n".join(lines)
