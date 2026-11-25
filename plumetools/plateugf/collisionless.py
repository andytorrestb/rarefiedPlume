r"""Collisionless number density around the plate: Cai's plume plus the plate's own gas.

Two parts, both collisionless:

* **the plume** -- Cai's free-plume density
  (:func:`plumetools.plateugf.sizing.density_table`, held to
  :mod:`plumetools.cai2012.analytical`), with nothing in the plate's shadow:
  behind it, wherever the straight line from the exit centre crosses it;
* **the plate gas** -- the upstream face re-emits Cai's incident number flux
  ``Gamma`` (:func:`plumetools.plateugf.caiplate.moments`) diffusely: cosine
  law, half-range Maxwellian at the wall temperature. From a patch of uniform
  ``Gamma`` the density at ``P`` is

  .. math::

      n(P) = \frac{\sqrt{\pi}}{2 c_w}\,\frac{\Gamma\,\Omega(P)}{\pi},
      \qquad c_w = \sqrt{2 k T_w / m},

  with ``Omega`` the patch's solid angle seen from ``P`` (closed form for a
  rectangle). Summed over patches it is exact for a piecewise-constant
  ``Gamma`` at any distance, and tends to ``Gamma sqrt(pi) / c_w`` on the wall.
  The WHOLE plate is summed: the mirrored quarters emit into this one too.

A flat plate's re-emitted gas never returns to it, and the domain's boundaries
only remove gas, so nothing else is missing -- except collisions, by design:
the simulation's departures from this field are the collisional effects.

On the 0.3M-cell meshes the plate gas is summed at every point; past
``GRID_MAX_POINTS`` (the 24-34M-cell meshes) it is tabulated on a (y, z, log h)
grid in front of the plate and interpolated, which smooths the step at the
plate's edge over one grid step (``GRID_STEP_M``) within a few centimetres of
the plate.

Used for the decomposition weights (:mod:`plumetools.plateugf.weights`) and
the density comparison (:mod:`plumetools.plateugf.viz`).
"""
from __future__ import annotations

import math

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from plumetools.nozzleplume import exitstate
from plumetools.nozzleplume.constants import BOLTZMANN_J_PER_K
from plumetools.plateugf import caiplate, sizing

#: Plate patch edge for the flux quadrature [m]: Gamma varies over metres.
PATCH_M = 0.25
#: Above this many points in front of the plate, tabulate instead of summing.
GRID_MAX_POINTS = 500_000
GRID_STEP_M = 0.125
GRID_H_POINTS = 48
#: Nearest tabulated distance from the plate [m]; closer points use it.
H_MIN_M = 1.0e-3


class Model:
    """The collisionless density of one plate case (``cfg``: nozzleplume's case config)."""

    def __init__(self, cfg, patch_m: float = PATCH_M):
        es = exitstate.from_config(cfg)
        plate = cfg.plate
        self.n0 = float(es.number_density_per_m3)
        self.radius_m = 0.5 * float(es.diameter_m)
        self.speed_ratio = float(es.speed_ratio)
        self.standoff_m = float(plate.standoff_m)
        self.half_y_m, self.half_z_m = 0.5 * float(plate.size_y_m), 0.5 * float(plate.size_z_m)
        ny = max(1, int(round(2.0 * self.half_y_m / patch_m)))
        nz = max(1, int(round(2.0 * self.half_z_m / patch_m)))
        ye = np.linspace(-self.half_y_m, self.half_y_m, ny + 1)
        ze = np.linspace(-self.half_z_m, self.half_z_m, nz + 1)
        yc, zc = np.meshgrid(0.5 * (ye[1:] + ye[:-1]), 0.5 * (ze[1:] + ze[:-1]), indexing="ij")
        centres = np.column_stack([np.full(yc.size, self.standoff_m), yc.ravel(), zc.ravel()])
        self.gamma = caiplate.moments(
            centres, radius_m=self.radius_m, number_density=self.n0, speed=es.speed_m_per_s,
            temperature_K=es.temperature_K, mass_kg=es.species.mass_kg)["flux"].reshape(ny, nz)
        # sum_p Gamma_p Omega_p = sum over patch corners of K_c F(corner - P), with
        # F(a, b) = atan(ab / (h sqrt(a^2 + b^2 + h^2))) and K_c the signed sum of
        # the (up to four) patches meeting at the corner.
        padded = np.zeros((ny + 2, nz + 2))
        padded[1:-1, 1:-1] = self.gamma
        k = padded[:-1, :-1] - padded[1:, :-1] - padded[:-1, 1:] + padded[1:, 1:]
        keep = k != 0.0
        self._cy = np.broadcast_to(ye[:, None], k.shape)[keep]
        self._cz = np.broadcast_to(ze[None, :], k.shape)[keep]
        self._k = k[keep]
        c_w = math.sqrt(2.0 * BOLTZMANN_J_PER_K * float(plate.temperature_K) / es.species.mass_kg)
        self.wall_density_factor = math.sqrt(math.pi) / c_w     # n on the wall = Gamma * this
        self._prefactor = self.wall_density_factor / (2.0 * math.pi)
        self._table = None

    # -- the plume ---------------------------------------------------------- #

    def plume(self, points: np.ndarray) -> np.ndarray:
        """Cai's plume density [1/m^3], zero in the plate's shadow and behind the exit plane."""
        points = np.atleast_2d(np.asarray(points, float))
        x, y, z = points[:, 0], np.abs(points[:, 1]), np.abs(points[:, 2])
        rho = np.hypot(y, z)
        x_need, rho_need = float(x.max(initial=1.0)), float(rho.max(initial=1.0))
        t = self._table
        if t is None or t.x_m[-1] < x_need or t.rho_m[-1] < rho_need:
            self._table = t = sizing.density_table(1.05 * x_need + 1.0, 1.05 * rho_need + 1.0,
                                                   self.radius_m, self.speed_ratio)
        n = self.n0 * t.n_ratio(x, rho)
        s = self.standoff_m
        with np.errstate(divide="ignore", invalid="ignore"):
            shadow = (x > s) & (y * s / x < self.half_y_m) & (z * s / x < self.half_z_m)
        n[shadow] = 0.0
        return n

    # -- the plate gas -------------------------------------------------------- #

    def _sum(self, h: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        out = np.empty(h.size)
        chunk = max(1, 20_000_000 // self._k.size)
        for i in range(0, h.size, chunk):
            hh = h[i:i + chunk, None]
            a = self._cy[None, :] - y[i:i + chunk, None]
            b = self._cz[None, :] - z[i:i + chunk, None]
            f = np.arctan(a * b / (hh * np.sqrt(a * a + b * b + hh * hh)))
            out[i:i + chunk] = f @ self._k
        return self._prefactor * out

    def plate_gas(self, points: np.ndarray, *, grid: bool | None = None) -> np.ndarray:
        """Density [1/m^3] of the gas the plate re-emits; zero at and behind its face.

        ``grid``: tabulate and interpolate (default: when more than
        ``GRID_MAX_POINTS`` points are in front of the plate).
        """
        points = np.atleast_2d(np.asarray(points, float))
        h = self.standoff_m - points[:, 0]
        front = np.flatnonzero(h > 0.0)
        out = np.zeros(len(points))
        if front.size == 0:
            return out
        y, z = np.abs(points[front, 1]), np.abs(points[front, 2])
        hf = np.maximum(h[front], H_MIN_M)
        if grid is None:
            grid = front.size > GRID_MAX_POINTS
        if not grid:
            out[front] = self._sum(hf, y, z)
            return out
        gy = np.arange(0.0, y.max() + 2.0 * GRID_STEP_M, GRID_STEP_M)
        gz = np.arange(0.0, z.max() + 2.0 * GRID_STEP_M, GRID_STEP_M)
        gh = np.geomspace(H_MIN_M, max(hf.max(), 2.0 * H_MIN_M), GRID_H_POINTS)
        Y, Z, H = np.meshgrid(gy, gz, gh, indexing="ij")
        values = self._sum(H.ravel(), Y.ravel(), Z.ravel()).reshape(Y.shape)
        table = RegularGridInterpolator((gy, gz, np.log(gh)), values)
        out[front] = table(np.column_stack([y, z, np.log(np.clip(hf, gh[0], gh[-1]))]))
        return out

    def density(self, points: np.ndarray, *, grid: bool | None = None) -> np.ndarray:
        """Plume plus plate gas [1/m^3]."""
        return self.plume(points) + self.plate_gas(points, grid=grid)
