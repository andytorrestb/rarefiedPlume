"""Uniform disk flux source and geometric (non-Maxwellian) ballistic f0."""
from __future__ import annotations

from dataclasses import dataclass
import math
from functools import lru_cache
import numpy as np
from scipy.special import erf, ndtri
from scipy.stats import qmc

from plumetools.cai2012.gas import ARGON, VhsSpecies, maxwellian_number_flux


@lru_cache(maxsize=32)
def _legendre(order):
    return np.polynomial.legendre.leggauss(order)


@dataclass(frozen=True)
class DiskSource:
    radius_m: float
    number_density_m3: float
    temperature_K: float
    speed_ratio: float
    species: VhsSpecies = ARGON

    def __post_init__(self):
        for name in ('radius_m', 'number_density_m3', 'temperature_K'):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(name + ' must be finite and positive')
        if not np.isfinite(self.speed_ratio) or self.speed_ratio < 0:
            raise ValueError('only finite nonnegative exit drift is supported')
        if self.species.internal_degrees_of_freedom != 0:
            raise ValueError('only single species monatomic elastic sources supported')

    @property
    def thermal_speed(self):
        return self.species.most_probable_speed(self.temperature_K)

    @property
    def emission_rate(self):
        return math.pi * self.radius_m**2 * maxwellian_number_flux(
            self.number_density_m3, self.speed_ratio, self.species, self.temperature_K)

    def sample_flux(self, count, seed):
        """Scrambled Sobol flux samples, not density-Maxwellian samples.

        Uniform disk area; vx density proportional to vx exp[-(vx/a-S)^2].
        Transverse components are independent N(0,a^2/2). The inverse axial CDF
        uses a fixed 32769-node interpolation, with a >=9 thermal-speed tail.
        """
        if count < 1:
            raise ValueError('count must be positive')
        u = qmc.Sobol(5, scramble=True, seed=int(seed)).random_base2(
            int(math.ceil(math.log2(count))))[:count]
        radius = self.radius_m * np.sqrt(u[:, 0])
        phi = 2*math.pi*u[:, 1]
        positions = np.column_stack((np.zeros(count), radius*np.cos(phi), radius*np.sin(phi)))
        s = self.speed_ratio
        grid = np.linspace(0, max(10., s+9.), 32769)
        cdf = (s*.5*(erf(grid-s)+erf(s))
               +(np.exp(-s*s)-np.exp(-(grid-s)**2))/(2*math.sqrt(math.pi)))
        norm = (math.exp(-s*s)+math.sqrt(math.pi)*s*(1+erf(s)))/(2*math.sqrt(math.pi))
        cdf /= norm
        vx = self.thermal_speed * np.interp(u[:, 2], cdf, grid)
        vt = self.thermal_speed / math.sqrt(2) * ndtri(u[:, 3:5])
        return positions, np.column_stack((vx, vt))

    def local_quadrature(self, point, options):
        """Velocities and density weights for actual visible-disk f0.

        Integrates solid angle directly: d^3v=c^2 dc dmu dphi. Disk intersection
        lengths determine mu bounds, removing the near-exit area-Jacobian peak.
        This is a directional/speed quadrature, never a local moment closure.
        """
        x, y, z = np.asarray(point, dtype=float)
        if x <= 0:
            return np.empty((0, 3)), np.empty(0)
        rho = math.hypot(y, z)
        er = np.array([y, z])/rho if rho > 0 else np.array([1., 0.])
        et = np.array([-er[1], er[0]])
        pm, pw = _legendre(options.partner_azimuthal)
        if rho < self.radius_m:
            # Full azimuth is periodic: the midpoint rule converges much faster
            # here than a nonperiodic Gauss rule at the same node count.
            phi = 2*math.pi*(np.arange(options.partner_azimuthal)+.5)/options.partner_azimuthal
            wphi = np.full(options.partner_azimuthal, 2*math.pi/options.partner_azimuthal)
            disc = np.sqrt(self.radius_m**2-rho**2*np.sin(phi)**2)
            low_l = np.zeros_like(phi)
            high_l = -rho*np.cos(phi)+disc
            qdir = np.cos(phi)[:, None]*er+np.sin(phi)[:, None]*et
        else:
            limit = math.asin(self.radius_m/rho)
            phi, wphi = limit*pm, limit*pw
            disc = np.sqrt(self.radius_m**2-rho**2*np.sin(phi)**2)
            low_l = rho*np.cos(phi)-disc
            high_l = rho*np.cos(phi)+disc
            qdir = -np.cos(phi)[:, None]*er+np.sin(phi)[:, None]*et
        mu_low = x/np.hypot(x, high_l)
        mu_high = x/np.hypot(x, low_l)
        gm, gw = _legendre(options.partner_radial)
        mu = mu_low[:, None]+.5*(gm+1)*(mu_high-mu_low)[:, None]
        wmu = .5*gw*(mu_high-mu_low)[:, None]
        tangent = np.sqrt(1-mu**2)
        directions = np.concatenate((mu[..., None], -tangent[..., None]*qdir[:, None, :]), axis=2)
        cn, cw = _legendre(options.partner_speed)
        cmax = self.thermal_speed*(self.speed_ratio+options.speed_tail_thermal)
        speed, wspeed = .5*cmax*(cn+1), .5*cmax*cw
        velocities = directions[:, :, None, :]*speed[None, None, :, None]
        drift = np.array([self.speed_ratio*self.thermal_speed, 0, 0])
        exponential = np.exp(-np.sum((velocities-drift)**2, axis=-1)/self.thermal_speed**2)
        weights = (self.number_density_m3/(math.pi**1.5*self.thermal_speed**3)
                   * exponential * speed[None, None, :]**2 * wspeed[None, None, :]
                   * wmu[:, :, None]*wphi[:, None, None])
        return velocities.reshape(-1, 3), weights.ravel()
