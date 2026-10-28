r"""Continuum-breakdown indicators on the actual 3-D CFD cells.

Definitions
-----------
With ``lambda`` the diagnostic mean free path of :mod:`.physics` (convention
``breakdown.mean_free_path_convention``, default ``vhs``) at the local ``n``
and ``T``:

=================  ========================================================
``Kn_rho``         ``lambda |grad rho| / rho``
``Kn_T``           ``lambda |grad T| / T``
``Kn_U``           ``lambda |grad |U|| / max(|U|, c)``, ``c = sqrt(gamma R T)``
                   -- the speed normalisation keeps it finite at stagnation
                   and in the slow background
``Kn_GLL``         ``max(Kn_rho, Kn_T, Kn_U)`` (Boyd, Chen & Candler 1995;
                   Wang & Boyd 2003's ``Kn_max``)
``B_exp``          ``-tau_c U . grad(ln rho)``, ``tau_c = lambda / c_bar``,
                   ``c_bar = sqrt(8RT/pi)``. **Signed**: positive where the
                   density falls along the flow (an expansion). A rate
                   parameter in the spirit of Bird (1970); not interchangeable
                   with other "Bird parameter" definitions, and its collision
                   time is nominal (no DSMC collision rate exists yet)
=================  ========================================================

Gradients are weighted least squares over face neighbours
(:func:`plumetools.rarefaction.core.least_squares_gradient`, OpenFOAM's
``leastSquares``), of the **unaveraged** CFD cell values, one-sided at
boundaries. They are exact for linear fields; at the source lip the field is
not smooth, so the indicators there are mesh-dependent by construction:
cells within ``breakdown.lip_exclusion_over_D`` of the lip circle are counted
separately (``near_lip``) and the report shows their share of every region.

Validity
--------
A cell with ``rho <= 0``, ``T <= 0``, a non-finite value or an unsolvable
gradient stencil is **invalid**: its indicators are NaN and it is never
counted as continuum. No denominator is regularised to make a bad cell look
benign.

Continuum region
----------------
A cell is continuum-valid when every indicator named in
``breakdown.criterion`` (default ``[Kn_GLL, B_exp]``; ``B_exp`` by magnitude)
is below its threshold. The retained region is the connected component(s) of
continuum cells that touch the source patch -- 3-D face adjacency,
``scipy.sparse.csgraph``. Disconnected continuum pockets (e.g. a recompression
region) are reported but handed to DSMC, as in Petersen et al. (2025).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from plumetools.foamio.volfield import read_internal_field
from plumetools.rarefaction.core import least_squares_gradient

from .meshgeom import MeshGeometry
from .physics import Gas

INDICATORS = ("Kn_rho", "Kn_T", "Kn_U", "Kn_GLL", "B_exp")


@dataclass(frozen=True)
class CfdState:
    rho: np.ndarray
    T: np.ndarray
    U: np.ndarray
    p: np.ndarray
    time: str

    @property
    def speed(self) -> np.ndarray:
        return np.linalg.norm(self.U, axis=1)


def read_state(case_dir: Path, time: str, gas: Gas) -> CfdState:
    d = Path(case_dir) / time
    t = read_internal_field(d / "T")
    u = read_internal_field(d / "U").reshape(-1, 3)
    p = read_internal_field(d / "p")
    rho = read_internal_field(d / "rho") if (d / "rho").is_file() else p / (gas.gas_constant * t)
    return CfdState(rho, t, u, p, time)


def gradients(mesh: MeshGeometry, state: CfdState, graph=None) -> dict:
    graph = graph or mesh.cell_graph()
    valid = (state.rho > 0) & (state.T > 0) & np.isfinite(state.rho) & np.isfinite(state.T) \
        & np.all(np.isfinite(state.U), axis=1)
    out = {"valid": valid}
    out["rho"] = least_squares_gradient(state.rho, valid, graph)
    out["T"] = least_squares_gradient(state.T, valid, graph)
    out["speed"] = least_squares_gradient(state.speed, valid, graph)
    out["U"] = [least_squares_gradient(state.U[:, k], valid, graph) for k in range(3)]
    return out


def grad_u_tensor(grads: dict) -> np.ndarray:
    """``(n, 3, 3)`` with ``[i, j] = dU_j/dx_i``."""
    return np.stack([g.vectors for g in grads["U"]], axis=2)


def indicators(state: CfdState, grads: dict, gas: Gas, convention: str = "vhs") -> dict:
    n = state.rho / gas.mass_kg
    lam = gas.mean_free_path(n, state.T, convention)
    c = gas.sound_speed(np.where(state.T > 0, state.T, np.nan))
    cbar = gas.mean_speed(np.where(state.T > 0, state.T, np.nan))
    with np.errstate(divide="ignore", invalid="ignore"):
        kn_rho = lam * grads["rho"].magnitude / state.rho
        kn_t = lam * grads["T"].magnitude / state.T
        kn_u = lam * grads["speed"].magnitude / np.maximum(state.speed, c)
        b = -(lam / cbar) * np.einsum("ij,ij->i", state.U, grads["rho"].vectors) / state.rho
    valid = (grads["valid"] & grads["rho"].valid & grads["T"].valid & grads["speed"].valid
             & np.isfinite(lam))
    out = {"lambda": np.where(valid, lam, np.nan)}
    for name, v in (("Kn_rho", kn_rho), ("Kn_T", kn_t), ("Kn_U", kn_u), ("B_exp", b)):
        out[name] = np.where(valid, v, np.nan)
    out["Kn_GLL"] = np.where(valid, np.fmax(np.fmax(kn_rho, kn_t), kn_u), np.nan)
    out["valid"] = valid
    return out


def continuum_mask(ind: dict, thresholds: dict, criterion=("Kn_GLL", "B_exp"),
                   factor: float = 1.0) -> np.ndarray:
    mask = ind["valid"].copy()
    for name in criterion:
        thr = float(thresholds[name]) * factor
        values = np.abs(ind[name]) if name == "B_exp" else ind[name]
        mask &= np.isfinite(values) & (values < thr)
    return mask


def components(mesh: MeshGeometry, mask: np.ndarray) -> tuple:
    """Connected components of the masked cells over internal faces."""
    o, nb = mesh.owner[:mesh.n_internal], mesh.neighbour
    both = mask[o] & mask[nb]
    g = coo_matrix((np.ones(both.sum()), (o[both], nb[both])),
                   shape=(mesh.n_cells, mesh.n_cells))
    n_comp, labels = connected_components(g, directed=False)
    labels = np.where(mask, labels, -1)
    return n_comp, labels


def source_connected(mesh: MeshGeometry, mask: np.ndarray, source_patch: str) -> dict:
    """The retained region: continuum components touching the source patch."""
    _, labels = components(mesh, mask)
    src_cells = mesh.owner[mesh.patch_slice(source_patch)]
    seeds = np.unique(labels[src_cells])
    seeds = seeds[seeds >= 0]
    retained = np.isin(labels, seeds) & mask
    pockets = np.unique(labels[mask & ~retained])
    return {"retained": retained, "labels": labels, "seed_components": seeds.tolist(),
            "n_disconnected_continuum_components": int(pockets.size),
            "disconnected_continuum_volume_m3": float(mesh.cell_volumes[mask & ~retained].sum())}


def lip_distance(mesh: MeshGeometry, radius_m: float) -> np.ndarray:
    c = mesh.cell_centres
    r = np.hypot(c[:, 1], c[:, 2])
    return np.hypot(c[:, 0], r - radius_m)


def summarise(mesh: MeshGeometry, ind: dict, retained: np.ndarray, near_lip: np.ndarray) -> dict:
    v = mesh.cell_volumes
    out = {"cells": mesh.n_cells, "invalid_cells": int((~ind["valid"]).sum()),
           "retained_cells": int(retained.sum()),
           "retained_volume_m3": float(v[retained].sum()),
           "retained_cells_near_lip": int((retained & near_lip).sum()),
           "near_lip_cells": int(near_lip.sum())}
    for name in INDICATORS:
        vals = ind[name]
        ok = np.isfinite(vals)
        out[name] = {
            "max": float(np.nanmax(vals)) if ok.any() else None,
            "max_away_from_lip": float(np.nanmax(np.where(near_lip, np.nan, vals)))
            if (ok & ~near_lip).any() else None,
            "max_in_retained": float(np.nanmax(np.where(retained, vals, np.nan)))
            if (ok & retained).any() else None,
        }
    return out


def write_fields(case_dir: Path, time: str, mesh: MeshGeometry, fields: dict) -> list:
    """Write diagnostic volScalarFields into the CFD time directory (for ParaView).

    Invalid cells are written as -1 and listed in the header: a masked cell
    must not be drawn in a colour that reads as "continuum".
    """
    out = []
    d = Path(case_dir) / time
    for name, values in fields.items():
        vals = np.where(np.isfinite(values), values, -1.0)
        lines = ["FoamFile", "{", "    version 2.0;", "    format ascii;",
                 "    class volScalarField;", f'    location "{time}";', f"    object {name};",
                 "}", f"// plumetools.hybridplume.breakdown: {name}; -1 marks an invalid cell",
                 "dimensions [0 0 0 0 0 0 0];", "",
                 f"internalField nonuniform List<scalar> {vals.size}", "("]
        lines += [f"{x:.8g}" for x in vals]
        lines += [")", ";", "", "boundaryField", "{"]
        for patch, info in mesh.patches.items():
            lines += [f"    {patch} {{ type zeroGradient; }}" if info.type not in (
                "empty", "symmetry", "wedge", "cyclic") else f"    {patch} {{ type {info.type}; }}"]
        lines += ["}", ""]
        path = d / name
        path.write_text("\n".join(lines))
        out.append(path)
    return out


def threshold_sensitivity(mesh: MeshGeometry, ind: dict, thresholds: dict, criterion,
                          factors, source_patch: str) -> list:
    rows = []
    for f in factors:
        mask = continuum_mask(ind, thresholds, criterion, factor=float(f))
        reg = source_connected(mesh, mask, source_patch)
        r = reg["retained"]
        rows.append({"factor": float(f),
                     "thresholds": {k: float(thresholds[k]) * float(f) for k in criterion},
                     "retained_cells": int(r.sum()),
                     "retained_volume_m3": float(mesh.cell_volumes[r].sum()),
                     "max_x_retained_m": float(mesh.cell_centres[r, 0].max()) if r.any() else None,
                     "n_disconnected_components": reg["n_disconnected_continuum_components"]})
    return rows
