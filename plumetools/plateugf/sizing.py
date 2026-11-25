r"""Cells sized by the local mean free path, as nozzleplume refinement cones.

The design (§3) sizes a cell at a point from Cai's collisionless density::

    lambda_local = lambda0 * n0 / n          VHS; temperature dependence ignored
    s = clamp(min(lambda_local, geo_fraction * r), floor, ceiling)

with ``r`` the distance from the exit centre. The ``r/20`` term resolves the
plume's expansion where the mean free path is long -- it sizes every cell at
Kn_e >= 0.1 -- and the floor keeps the densest cases buildable; USP handles the
cells left larger than lambda there.

nozzleplume meshes an octree: level ``L`` cells are ``background / 2**L``. A
point needs level ``L`` or finer when ``s < background / 2**(L-1)``, so the
region of level ``L`` is that set. Each region is enclosed in the coaxial cone
nozzleplume accepts (``mesh.refinement_regions``, ``shape: cone``), fitted on a
dense sample grid with a margin of one level-``L`` cell plus the grid spacing.
The plume is axisymmetric about +x, so a point is ``(x, rho)``.

Gas reflected off the plate is not in this field (spec risk R4).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from plumetools.cai2012 import analytical

#: The smallest x sampled: half the finest (10 mm) cell, the first cell centre.
#: The collisionless integrand is singular in the exit plane itself, and
#: analytical.density_ratio raises its quadrature order as 1/x for the WHOLE
#: call -- hence one call per x row below.
X_MIN_M = 5.0e-3
_RHO_SHIFT = 1.0e-4
_TINY = 1.0e-300


@dataclass(frozen=True)
class DensityTable:
    """Collisionless n/n0 on a log-spaced (x, rho) grid, log-interpolated."""

    x_m: np.ndarray
    rho_m: np.ndarray
    log_n: np.ndarray

    def __post_init__(self):
        object.__setattr__(self, "_interp", RegularGridInterpolator(
            (np.log(self.x_m), np.log(self.rho_m + _RHO_SHIFT)), self.log_n,
            bounds_error=False, fill_value=None))

    def n_ratio(self, x, rho) -> np.ndarray:
        """n/n0 at ``(x, rho)``; 0 at and behind the exit plane."""
        x, rho = np.broadcast_arrays(np.asarray(x, float), np.asarray(rho, float))
        out = np.zeros(x.shape)
        ahead = x > 0.0
        xs = np.clip(x[ahead], self.x_m[0], self.x_m[-1])
        rs = np.clip(rho[ahead], 0.0, self.rho_m[-1])
        out[ahead] = np.exp(self._interp(
            np.column_stack([np.log(xs), np.log(rs + _RHO_SHIFT)])))
        return out


def density_table(x_max_m: float, rho_max_m: float, radius_m: float,
                  speed_ratio: float, n_points: int = 80) -> DensityTable:
    """Tabulate Cai's collisionless density out to ``x_max_m`` and ``rho_max_m``."""
    x = np.geomspace(X_MIN_M, float(x_max_m), n_points)
    rho = np.concatenate([[0.0], np.geomspace(_RHO_SHIFT, float(rho_max_m), n_points)])
    n = np.array([
        analytical.density_ratio(np.full_like(rho, xi), rho, np.zeros_like(rho),
                                 radius_m, speed_ratio, n_radial=32, n_azimuthal=64)
        for xi in x])
    return DensityTable(x_m=x, rho_m=rho, log_n=np.log(np.maximum(n, _TINY)))


def target_cell_size(table: DensityTable, x, rho, *, knudsen: float,
                     diameter_m: float, floor_m: float = 0.01,
                     ceiling_m: float = 0.64, geo_fraction: float = 0.05) -> np.ndarray:
    """The design's cell size at ``(x, rho)`` [m]."""
    x, rho = np.broadcast_arrays(np.asarray(x, float), np.asarray(rho, float))
    n = table.n_ratio(x, rho)
    lam0 = float(knudsen) * float(diameter_m)
    with np.errstate(divide="ignore"):
        lam = np.where(n > 0.0, lam0 / np.maximum(n, _TINY), np.inf)
    r = np.hypot(x, rho)
    s = np.minimum(lam, geo_fraction * r)
    return np.clip(s, floor_m, ceiling_m)


def _up(value: float) -> float:
    """Round up to a whole millimetre."""
    return math.ceil(round(value * 1000.0, 6)) / 1000.0


def _fit_grid(x_limit_m: float, rho_limit_m: float):
    x = np.unique(np.concatenate([np.geomspace(X_MIN_M, x_limit_m, 200),
                                  np.linspace(X_MIN_M, x_limit_m, 400)]))
    rho = np.unique(np.concatenate([[0.0], np.geomspace(_RHO_SHIFT, rho_limit_m, 200),
                                    np.linspace(0.0, rho_limit_m, 500)]))
    return x, rho, x_limit_m / 399.0, rho_limit_m / 499.0


def refinement_regions(table: DensityTable, *, knudsen: float, diameter_m: float,
                       x_limit_m: float, rho_limit_m: float,
                       background_m: float = 0.64, max_level: int = 6,
                       floor_m: float = 0.01, geo_fraction: float = 0.05) -> list[dict]:
    """One cone per octree level 1..``max_level``, as ``RegionConfig`` mappings.

    Empty levels are left out. Every cone is clipped to ``x_limit_m`` and
    ``rho_limit_m`` (the domain), never refused.
    """
    x, rho, dx, drho = _fit_grid(float(x_limit_m), float(rho_limit_m))
    X, R = np.meshgrid(x, rho, indexing="ij")
    s = target_cell_size(table, X, R, knudsen=knudsen, diameter_m=diameter_m,
                         floor_m=floor_m, ceiling_m=background_m,
                         geo_fraction=geo_fraction)
    radius = 0.5 * float(diameter_m)
    x_limit, rho_limit = float(x_limit_m), float(rho_limit_m)
    out = []
    for level in range(1, int(max_level) + 1):
        member = s < background_m / 2 ** (level - 1)
        rows = np.flatnonzero(member.any(axis=1))
        if rows.size == 0:
            continue
        cell = background_m / 2 ** level
        margin = cell + drho
        envelope = np.where(member, R, -np.inf).max(axis=1)
        x_end = min(float(x[rows].max()) + cell + dx, x_limit)
        for k, (xa, xb) in enumerate(_segments(x_end, cell)):
            inside = (x >= xa - dx) & (x <= xb + dx) & np.isfinite(envelope)
            if not inside.any():
                continue
            env = envelope[inside] + margin
            xs = x[inside]
            r_a = float(env[0]) if xa > 0.0 else max(float(env[0]), radius + margin)
            r_a = max(r_a, float(env[xs <= xa + dx].max()) if (xs <= xa + dx).any() else r_a)
            later = xs > xa + 1e-12
            slope = float(np.max((env[later] - r_a) / (xs[later] - xa))) if later.any() else 0.0
            r_b = r_a + slope * (xb - xa)
            if xa == 0.0:
                # The exit disk sits wholly in the first segment of every level:
                # nozzleplume wants it at one level, and snappy refines by centre.
                r_b = max(r_b, radius + margin) if r_a >= radius + margin else r_b
            out.append({
                "name": f"mfpL{level}s{k}",
                "level": level,
                "shape": "cone",
                "x_min_m": round(xa, 6),
                "x_max_m": min(_up(xb), x_limit),
                "radius1_m": min(_up(r_a), rho_limit),
                "radius2_m": min(_up(max(r_b, 0.0)), rho_limit),
            })
    return out


def _segments(x_end: float, cell: float, n: int = 10) -> list[tuple[float, float]]:
    """Contiguous x intervals [0, x_end]: geometric, none shorter than 4 cells."""
    if x_end <= 8.0 * cell:
        return [(0.0, x_end)]
    first = max(4.0 * cell, x_end / 2 ** n)
    edges = np.unique(np.concatenate([[0.0], np.geomspace(first, x_end, n)]))
    return [(float(a), float(b)) for a, b in zip(edges[:-1], edges[1:])]
