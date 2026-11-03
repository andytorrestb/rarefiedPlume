"""A small, valid 3-D hex box ``polyMesh`` with the hybrid family's patches.

Written by Python so the hybridplume geometry, breakdown, interface and mapping
code can be tested without OpenFOAM. Faces follow OpenFOAM's conventions:
internal faces first, owner < neighbour, area vectors from owner to neighbour;
boundary faces grouped by patch, area vectors out of the domain.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

HEADER = ("FoamFile\n{{\n    version 2.0;\n    format ascii;\n    class {cls};\n"
          "    object {obj};\n}}\n\n")


def box_mesh(case: Path, edges_x, edges_y, edges_z, radius: float) -> dict:
    ex, ey, ez = (np.asarray(e, float) for e in (edges_x, edges_y, edges_z))
    nx, ny, nz = ex.size - 1, ey.size - 1, ez.size - 1

    def pid(i, j, k):
        return i + (nx + 1) * (j + (ny + 1) * k)

    def cid(i, j, k):
        return i + nx * (j + ny * k)

    gx, gy, gz = np.meshgrid(ex, ey, ez, indexing="ij")
    points = np.zeros(((nx + 1) * (ny + 1) * (nz + 1), 3))
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                points[pid(i, j, k)] = (ex[i], ey[j], ez[k])

    def quad(axis, i, j, k):
        """The face normal to ``axis`` at lattice corner (i, j, k), +axis oriented."""
        if axis == 0:
            q = [pid(i, j, k), pid(i, j + 1, k), pid(i, j + 1, k + 1), pid(i, j, k + 1)]
        elif axis == 1:
            q = [pid(i, j, k), pid(i, j, k + 1), pid(i + 1, j, k + 1), pid(i + 1, j, k)]
        else:
            q = [pid(i, j, k), pid(i + 1, j, k), pid(i + 1, j + 1, k), pid(i, j + 1, k)]
        p = points[q]
        n = np.cross(p[1] - p[0], p[2] - p[0])
        return q if n[axis] > 0 else q[::-1]

    faces, owner, neigh = [], [], []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c = cid(i, j, k)
                for axis, (di, dj, dk) in enumerate(((1, 0, 0), (0, 1, 0), (0, 0, 1))):
                    a, b, cc = i + di, j + dj, k + dk
                    if a < nx and b < ny and cc < nz:
                        faces.append(quad(axis, a, b, cc))
                        owner.append(c)
                        neigh.append(cid(a, b, cc))
    patches = {"source": [], "exitPlane": [], "farField": []}
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c = cid(i, j, k)
                if i == 0:
                    yc, zc = 0.5 * (ey[j] + ey[j + 1]), 0.5 * (ez[k] + ez[k + 1])
                    name = "source" if yc * yc + zc * zc <= radius ** 2 else "exitPlane"
                    patches[name].append((quad(0, 0, j, k)[::-1], c))
                if i == nx - 1:
                    patches["farField"].append((quad(0, nx, j, k), c))
                if j == 0:
                    patches["farField"].append((quad(1, i, 0, k)[::-1], c))
                if j == ny - 1:
                    patches["farField"].append((quad(1, i, ny, k), c))
                if k == 0:
                    patches["farField"].append((quad(2, i, j, 0)[::-1], c))
                if k == nz - 1:
                    patches["farField"].append((quad(2, i, j, nz), c))
    starts = {}
    for name, entries in patches.items():
        starts[name] = (len(faces), len(entries))
        for q, c in entries:
            faces.append(q)
            owner.append(c)

    poly = Path(case) / "constant" / "polyMesh"
    poly.mkdir(parents=True, exist_ok=True)
    (poly / "points").write_text(
        HEADER.format(cls="vectorField", obj="points") + f"{len(points)}\n(\n"
        + "\n".join(f"({x:.12g} {y:.12g} {z:.12g})" for x, y, z in points) + "\n)\n")
    (poly / "faces").write_text(
        HEADER.format(cls="faceList", obj="faces") + f"{len(faces)}\n(\n"
        + "\n".join(f"4({' '.join(map(str, q))})" for q in faces) + "\n)\n")
    (poly / "owner").write_text(
        HEADER.format(cls="labelList", obj="owner") + f"{len(owner)}\n(\n"
        + "\n".join(map(str, owner)) + "\n)\n")
    (poly / "neighbour").write_text(
        HEADER.format(cls="labelList", obj="neighbour") + f"{len(neigh)}\n(\n"
        + "\n".join(map(str, neigh)) + "\n)\n")
    body = "".join(f"    {name}\n    {{\n        type            patch;\n"
                   f"        nFaces          {n};\n        startFace       {s};\n    }}\n"
                   for name, (s, n) in starts.items())
    (poly / "boundary").write_text(
        HEADER.format(cls="polyBoundaryMesh", obj="boundary") + f"{len(starts)}\n(\n"
        + body + ")\n")
    return {"n_cells": nx * ny * nz, "shape": (nx, ny, nz), "patches": starts,
            "cid": cid}
