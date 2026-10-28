r"""Select a fixed CFD -> DSMC coupling surface on the 3-D mesh, and assess it.

The surface is a set of actual mesh faces: every internal face with exactly one
cell in the retained continuum region ``R``. Nothing is reduced to a centreline
radius or an azimuthal average. ``R`` comes from one of three policies:

* ``automatic`` -- the source-connected continuum component from
  :mod:`.breakdown` (thresholds never loosened to manufacture a region);
* ``manual`` -- a sphere, cylinder or box, reduced to the part connected to the
  source and then evaluated by **the same** diagnostics;
* ``none`` -- empty: full DSMC.

Orientation
-----------
``n`` points from the CFD region into the DSMC region. For an internal face
``f`` with area vector ``S_f`` (owner -> neighbour) that is ``+S_f/|S_f|`` when
the owner is in ``R`` and ``-S_f/|S_f|`` otherwise.

What is reported per face and integrated
----------------------------------------
state (``rho``, ``n``, ``T``, ``U``), ``Un``, ``Mn = Un/c``, the indicators on
both sides, the continuum NSF flux, the Maxwellian forward (injected) and
reverse (returning) half-range fluxes, and the local DSMC resolution. Integrated
budgets close against the region's whole boundary: interface + covered source +
any outer-boundary contact, whose outward area vectors must sum to zero.

Suitability
-----------
A fixed one-way surface is ``suitable`` only if every criterion in
``interface.criteria`` holds; otherwise the verdict is ``unsuitable`` with the
reasons. The criteria are a policy, not physics, and are all reported with
their measured values:

* returning mass fraction (reverse/forward Maxwellian mass flux) -- one-way
  coupling ignores what comes back;
* area fraction with subsonic normal Mach number -- an oblique surface can have
  ``|U| > c`` and ``Un < c``; total Mach alone proves nothing;
* outer-boundary contact -- a retained region that reaches the far field is
  not enclosed by the surface;
* overlap layers -- DSMC-side cell layers in which the CFD is still within
  ``overlap_ceiling_factor`` of the thresholds;
* every interface state valid, every retained cell valid.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from . import halfrange
from .breakdown import continuum_mask
from .fluxes import continuum_flux
from .meshgeom import MeshGeometry
from .mesh import PATCHES

STATUS_SUITABLE = "suitable"
STATUS_UNSUITABLE = "unsuitable"
STATUS_NONE = "no_continuum_region"
STATUS_FULL_DSMC = "full_dsmc"


@dataclass
class Surface:
    """The faces of a retained region and everything sampled on them."""

    retained: np.ndarray            # (n_cells,) bool
    faces: np.ndarray               # interface face labels
    cfd_cell: np.ndarray
    dsmc_cell: np.ndarray
    normal: np.ndarray              # (k, 3) CFD -> DSMC
    area: np.ndarray
    centre: np.ndarray
    source_covered: np.ndarray      # source face labels whose owner is in R
    source_uncovered: np.ndarray
    outer_contact: np.ndarray       # exitPlane/farField faces owned by R


def manual_region(mesh: MeshGeometry, spec: dict, diameter_m: float) -> np.ndarray:
    c = mesh.cell_centres
    kind = spec["type"]
    if kind == "sphere":
        centre = np.asarray(spec.get("centre_over_D", [0, 0, 0]), float) * diameter_m
        inside = np.linalg.norm(c - centre, axis=1) < float(spec["radius_over_D"]) * diameter_m
    elif kind == "cylinder":
        r = np.hypot(c[:, 1], c[:, 2])
        inside = (c[:, 0] < float(spec["length_over_D"]) * diameter_m) & (
            r < float(spec["radius_over_D"]) * diameter_m)
    elif kind == "box":
        lo = np.asarray(spec["lo_over_D"], float) * diameter_m
        hi = np.asarray(spec["hi_over_D"], float) * diameter_m
        inside = np.all((c > lo) & (c < hi), axis=1)
    else:
        raise ValueError(f"unknown manual surface {kind!r}")
    return inside


def connected_to_source(mesh: MeshGeometry, region: np.ndarray) -> np.ndarray:
    o, nb = mesh.owner[:mesh.n_internal], mesh.neighbour
    both = region[o] & region[nb]
    g = coo_matrix((np.ones(both.sum()), (o[both], nb[both])),
                   shape=(mesh.n_cells, mesh.n_cells))
    _, labels = connected_components(g, directed=False)
    src = mesh.owner[mesh.patch_slice(PATCHES["source"])]
    seeds = np.unique(labels[src[region[src]]])
    return region & np.isin(labels, seeds)


def extract(mesh: MeshGeometry, retained: np.ndarray) -> Surface:
    o, nb = mesh.owner[:mesh.n_internal], mesh.neighbour
    cut = retained[o] != retained[nb]
    faces = np.nonzero(cut)[0]
    owner_in = retained[o[faces]]
    s = mesh.face_areas[faces]
    a = np.linalg.norm(s, axis=1)
    normal = np.where(owner_in[:, None], s, -s) / a[:, None]
    cfd_cell = np.where(owner_in, o[faces], nb[faces])
    dsmc_cell = np.where(owner_in, nb[faces], o[faces])
    src = np.arange(mesh.patch_slice(PATCHES["source"]).start,
                    mesh.patch_slice(PATCHES["source"]).stop)
    covered = src[retained[mesh.owner[src]]]
    uncovered = src[~retained[mesh.owner[src]]]
    outer = []
    for key in ("exit_plane", "far_field"):
        sl = mesh.patch_slice(PATCHES[key])
        f = np.arange(sl.start, sl.stop)
        outer.append(f[retained[mesh.owner[f]]])
    return Surface(retained, faces, cfd_cell, dsmc_cell, normal, a,
                   mesh.face_centres[faces], covered, uncovered,
                   np.concatenate(outer) if outer else np.empty(0, int))


def face_weights(mesh: MeshGeometry, faces: np.ndarray) -> np.ndarray:
    """OpenFOAM linear interpolation weight of the owner, per internal face."""
    s = mesh.face_areas[faces]
    d_o = np.abs(np.einsum("ij,ij->i", s, mesh.face_centres[faces]
                           - mesh.cell_centres[mesh.owner[faces]]))
    d_n = np.abs(np.einsum("ij,ij->i", s, mesh.cell_centres[mesh.neighbour[faces]]
                           - mesh.face_centres[faces]))
    return d_n / (d_o + d_n)


def face_values(mesh: MeshGeometry, surf: Surface, cell_values: np.ndarray,
                method: str = "linear") -> np.ndarray:
    """Values on the interface faces: linear between the two cells, or the CFD cell's."""
    v = np.asarray(cell_values, float)
    if method == "continuum_cell":
        return v[surf.cfd_cell]
    if method != "linear":
        raise ValueError(f"unknown face-state method {method!r}")
    w = face_weights(mesh, surf.faces)
    o = mesh.owner[surf.faces]
    nb = mesh.neighbour[surf.faces]
    shape = (-1,) + (1,) * (v.ndim - 1)
    return w.reshape(shape) * v[o] + (1.0 - w).reshape(shape) * v[nb]


def closure(mesh: MeshGeometry, surf: Surface) -> dict:
    """Sum of outward area vectors over the region's boundary; zero if closed."""
    total = (surf.normal * surf.area[:, None]).sum(axis=0)
    total += mesh.face_areas[surf.source_covered].sum(axis=0)
    total += mesh.face_areas[surf.outer_contact].sum(axis=0)
    scale = surf.area.sum() + np.linalg.norm(mesh.face_areas[surf.source_covered], axis=1).sum()
    return {"sum_outward_area_vector_m2": total.tolist(),
            "relative_to_boundary_area": float(np.linalg.norm(total) / scale) if scale > 0 else None}


def overlap_layers(mesh: MeshGeometry, surf: Surface, ind: dict, thresholds: dict, criterion,
                   ceiling: float, n_layers: int = 4, required_fraction: float = 0.9) -> dict:
    ok_ceiling = continuum_mask(ind, thresholds, criterion, factor=ceiling)
    o, nb = mesh.owner[:mesh.n_internal], mesh.neighbour
    seen = surf.retained.copy()
    layer = np.zeros(mesh.n_cells, bool)
    layer[surf.dsmc_cell] = True
    rows = []
    consecutive = 0
    still = True
    for k in range(1, n_layers + 1):
        vol = mesh.cell_volumes[layer]
        frac = float(vol[ok_ceiling[layer]].sum() / vol.sum()) if vol.sum() > 0 else 0.0
        rows.append({"layer": k, "cells": int(layer.sum()),
                     "volume_fraction_within_ceiling": frac})
        if still and frac >= required_fraction:
            consecutive += 1
        else:
            still = False
        seen |= layer
        nxt = np.zeros(mesh.n_cells, bool)
        a = layer[o] & ~seen[nb]
        b = layer[nb] & ~seen[o]
        nxt[nb[a]] = True
        nxt[o[b]] = True
        layer = nxt
        if not layer.any():
            break
    return {"ceiling_factor": ceiling, "required_volume_fraction": required_fraction,
            "layers": rows, "overlap_layers": consecutive}


def assess(mesh: MeshGeometry, surf: Surface, state, grads, ind: dict, cfg, source,
           *, mode: str) -> dict:
    gas = source.gas
    crit = cfg.get("interface.criteria", {}) or {}
    thresholds = cfg.get("breakdown.thresholds")
    criterion = tuple(cfg.get("breakdown.criterion", ["Kn_GLL", "B_exp"]))
    method = cfg.get("interface.face_state", "linear")
    out = {"mode": mode, "face_state_method": method}
    reasons = []

    retained = surf.retained
    out["retained_cells"] = int(retained.sum())
    out["retained_volume_m3"] = float(mesh.cell_volumes[retained].sum())
    if mode == "none":
        out["status"] = STATUS_FULL_DSMC
        return out
    if not retained.any():
        out["status"] = STATUS_NONE
        out["reasons"] = ["no continuum cell is connected to the source"]
        return out

    # -- coverage ------------------------------------------------------------
    src_area = np.linalg.norm(mesh.face_areas[mesh.patch_slice(PATCHES["source"])], axis=1)
    cov_area = np.linalg.norm(mesh.face_areas[surf.source_covered], axis=1).sum()
    outer_area = np.linalg.norm(mesh.face_areas[surf.outer_contact], axis=1).sum()
    out["source"] = {"faces": int(src_area.size), "covered_faces": int(surf.source_covered.size),
                     "covered_area_fraction": float(cov_area / src_area.sum()),
                     "uncovered_faces_injected_directly_by_dsmc": int(surf.source_uncovered.size)}
    out["interface_faces"] = int(surf.faces.size)
    out["interface_area_m2"] = float(surf.area.sum())
    out["outer_boundary_contact_area_m2"] = float(outer_area)
    contact_fraction = float(outer_area / (outer_area + surf.area.sum()))
    out["outer_boundary_contact_fraction"] = contact_fraction
    out["closure"] = closure(mesh, surf)
    bad_retained = retained & ~ind["valid"]
    out["invalid_retained_cells"] = int(bad_retained.sum())
    noncont = retained & ~continuum_mask(ind, thresholds, criterion)
    out["noncontinuum_retained_volume_fraction"] = float(
        mesh.cell_volumes[noncont].sum() / mesh.cell_volumes[retained].sum())

    # -- face states ---------------------------------------------------------
    rho = face_values(mesh, surf, state.rho, method)
    t = face_values(mesh, surf, state.T, method)
    u = face_values(mesh, surf, state.U, method)
    gu = face_values(mesh, surf, np.stack([g.vectors for g in grads["U"]], axis=2), method)
    gt = face_values(mesh, surf, grads["T"].vectors, method)
    valid = (rho > 0) & (t > 0) & np.isfinite(rho) & np.isfinite(t) & np.all(np.isfinite(u), 1)
    grads_ok = np.all(np.isfinite(gu.reshape(-1, 9)), 1) & np.all(np.isfinite(gt), 1)
    out["invalid_face_states"] = int((~valid).sum())
    out["faces_without_gradients"] = int((~grads_ok).sum())
    un = np.einsum("ij,ij->i", u, surf.normal)
    c = gas.sound_speed(np.where(t > 0, t, np.nan))
    mn = un / c
    area = surf.area
    a_tot = area.sum()
    out["normal_mach"] = {
        "min": float(np.nanmin(mn)), "max": float(np.nanmax(mn)),
        "area_weighted_mean": float(np.nansum(mn * area) / a_tot),
        "subsonic_area_fraction": float(area[mn < 1.0].sum() / a_tot),
        "inflow_to_cfd_area_fraction": float(area[un < 0.0].sum() / a_tot),
        "total_mach_area_weighted": float(np.nansum(np.linalg.norm(u, axis=1) / c * area) / a_tot),
    }

    # -- fluxes --------------------------------------------------------------
    use = valid
    fwd = halfrange.incoming_flux(rho[use], u[use], t[use], surf.normal[use], gas.gas_constant)
    rev = halfrange.incoming_flux(rho[use], u[use], t[use], -surf.normal[use], gas.gas_constant)
    nsf = continuum_flux(rho[use], u[use], t[use], surf.normal[use], gas,
                         grad_u=np.where(grads_ok[use, None, None], gu[use], 0.0),
                         grad_t=np.where(grads_ok[use, None], gt[use], 0.0))
    euler = continuum_flux(rho[use], u[use], t[use], surf.normal[use], gas, inviscid=True)
    A = area[use]
    integ = {k: halfrange.integrate(v, A) for k, v in (("maxwellian_forward", fwd),
                                                        ("maxwellian_reverse", rev))}
    integ["maxwellian_net"] = {
        "mass": integ["maxwellian_forward"]["mass"] - integ["maxwellian_reverse"]["mass"],
        # momentum crossing towards -n carries -(J2 n) ...: net = fwd - rev
        "momentum": (np.asarray(integ["maxwellian_forward"]["momentum"])
                     - np.asarray(integ["maxwellian_reverse"]["momentum"])).tolist(),
        "energy": integ["maxwellian_forward"]["energy"] - integ["maxwellian_reverse"]["energy"],
    }
    integ["continuum_nsf"] = halfrange.integrate(nsf, A)
    integ["continuum_euler"] = halfrange.integrate(euler, A)
    integ["continuum_viscous_parts"] = {
        "viscous_momentum": (nsf["viscous_momentum"] * A[:, None]).sum(0).tolist(),
        "viscous_energy": float(np.dot(nsf["viscous_energy"], A)),
        "conductive_energy": float(np.dot(nsf["conductive_energy"], A)),
    }
    m = gas.mass_kg
    integ["number_rates_per_s"] = {
        "maxwellian_forward": integ["maxwellian_forward"]["mass"] / m,
        "maxwellian_reverse": integ["maxwellian_reverse"]["mass"] / m,
        "continuum": integ["continuum_nsf"]["mass"] / m}
    out["fluxes"] = integ
    ret = integ["maxwellian_reverse"]["mass"] / integ["maxwellian_forward"]["mass"]
    out["return_mass_fraction"] = float(ret)
    out["maxwellian_vs_nsf"] = {
        "mass_rel": float(integ["maxwellian_net"]["mass"] / integ["continuum_nsf"]["mass"] - 1),
        "energy_rel": float(integ["maxwellian_net"]["energy"] / integ["continuum_nsf"]["energy"] - 1),
        "momentum_abs_diff": (np.asarray(integ["maxwellian_net"]["momentum"])
                              - np.asarray(integ["continuum_nsf"]["momentum"])).tolist(),
        "note": "a Maxwellian carries no viscous stress or heat flux; the "
                "difference is the NSF terms the inlet cannot represent"}

    # -- indicators on the surface ------------------------------------------
    for side, cells in (("cfd_side", surf.cfd_cell), ("dsmc_side", surf.dsmc_cell)):
        out[f"indicators_{side}"] = {
            name: {"max": float(np.nanmax(ind[name][cells])) if np.isfinite(ind[name][cells]).any() else None,
                   "area_weighted_mean": float(np.nansum(ind[name][cells] * area)
                                               / area[np.isfinite(ind[name][cells])].sum())
                   if np.isfinite(ind[name][cells]).any() else None}
            for name in ("Kn_rho", "Kn_T", "Kn_U", "Kn_GLL", "B_exp")}

    ov = overlap_layers(mesh, surf, ind, thresholds, criterion,
                        float(cfg.get("breakdown.overlap_ceiling_factor", 2.0)))
    out["overlap"] = ov

    # -- verdict -------------------------------------------------------------
    def check(name, value, limit, ok):
        out.setdefault("criteria", {})[name] = {"value": value, "limit": limit, "pass": bool(ok)}
        if not ok:
            reasons.append(f"{name} = {value:.4g} (limit {limit})")

    check("return_mass_fraction", float(ret),
          crit.get("max_return_mass_fraction", 0.05),
          ret <= crit.get("max_return_mass_fraction", 0.05))
    sub = out["normal_mach"]["subsonic_area_fraction"]
    check("subsonic_normal_area_fraction", sub,
          crit.get("max_subsonic_normal_area_fraction", 0.25),
          sub <= crit.get("max_subsonic_normal_area_fraction", 0.25))
    check("outer_boundary_contact_fraction", contact_fraction,
          crit.get("max_outer_boundary_contact_fraction", 0.0),
          contact_fraction <= crit.get("max_outer_boundary_contact_fraction", 0.0))
    check("overlap_layers", float(ov["overlap_layers"]),
          crit.get("min_overlap_layers", 1),
          ov["overlap_layers"] >= crit.get("min_overlap_layers", 1))
    check("invalid_face_states", float(out["invalid_face_states"]), 0,
          out["invalid_face_states"] == 0)
    check("invalid_retained_cells", float(out["invalid_retained_cells"]), 0,
          out["invalid_retained_cells"] == 0)
    if mode == "manual":
        v = out["noncontinuum_retained_volume_fraction"]
        lim = crit.get("max_noncontinuum_retained_volume_fraction", 0.0)
        check("noncontinuum_retained_volume_fraction", v, lim, v <= lim)
    out["status"] = STATUS_SUITABLE if not reasons else STATUS_UNSUITABLE
    out["reasons"] = reasons
    out["_face_arrays"] = {"rho": rho, "T": t, "U": u, "valid": valid, "Mn": mn, "Un": un,
                           "grad_U": gu, "grad_T": gt}
    return out


def write_vtk(path: Path, mesh: MeshGeometry, faces: np.ndarray, cell_data: dict) -> Path:
    """Legacy ASCII VTK polydata of the given mesh faces, with per-face data."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    polys = [mesh.labels[mesh.offsets[f]:mesh.offsets[f + 1]] for f in faces]
    used = np.unique(np.concatenate(polys)) if polys else np.empty(0, int)
    remap = {int(p): i for i, p in enumerate(used)}
    lines = ["# vtk DataFile Version 3.0", "hybridplume interface", "ASCII",
             "DATASET POLYDATA", f"POINTS {used.size} double"]
    lines += [f"{x:.9g} {y:.9g} {z:.9g}" for x, y, z in mesh.points[used]]
    size = sum(len(p) + 1 for p in polys)
    lines.append(f"POLYGONS {len(polys)} {size}")
    lines += [f"{len(p)} " + " ".join(str(remap[int(v)]) for v in p) for p in polys]
    lines.append(f"CELL_DATA {len(polys)}")
    for name, values in cell_data.items():
        values = np.asarray(values, float)
        values = np.where(np.isfinite(values), values, -1.0)
        if values.ndim == 2:
            lines.append(f"VECTORS {name} double")
            lines += [f"{a:.6g} {b:.6g} {c:.6g}" for a, b, c in values]
        else:
            lines += [f"SCALARS {name} double 1", "LOOKUP_TABLE default"]
            lines += [f"{v:.6g}" for v in values]
    path.write_text("\n".join(lines) + "\n")
    return path
