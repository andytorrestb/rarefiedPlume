"""Reduced, mutually thinned representative pair events in a finite plume.

There is no collision time stepping or feedback into f0. See the mathematical
contract in docs/cai2012-virtual-collision.md before interpreting predictions.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import hashlib
import json
import math
import time
import numpy as np
from scipy.special import erf, erfc

from plumetools.cai2012 import analytical
from .moments import derive

MODEL_VERSION = 'mutual-first-event-1'


@dataclass(frozen=True)
class NumericalOptions:
    samples: int = 128
    hazard_steps: int = 12
    partner_radial: int = 4
    partner_azimuthal: int = 12
    partner_speed: int = 16
    seed: int = 2026
    domain_length_D: float = 20.
    domain_radius_D: float = 12.
    batches: int = 4
    chunk_size: int = 128
    bin_order: int = 2
    residence_order: int = 4
    cai_radial: int = 48
    cai_azimuthal: int = 96
    speed_tail_thermal: float = 8.

    def __post_init__(self):
        for key in ('samples', 'hazard_steps', 'partner_radial', 'partner_azimuthal',
                    'partner_speed', 'batches', 'chunk_size', 'bin_order', 'residence_order',
                    'cai_radial', 'cai_azimuthal'):
            value = getattr(self, key)
            if int(value) != value or value < 2:
                raise ValueError(key+' must be an integer >=2')
            object.__setattr__(self, key, int(value))
        if self.cai_azimuthal < 4:
            raise ValueError('cai_azimuthal must be >=4, as in the existing Cai evaluator')
        if not np.isfinite(self.seed) or int(self.seed) != self.seed or self.seed < 0:
            raise ValueError('seed must be a nonnegative integer')
        object.__setattr__(self, 'seed', int(self.seed))
        if self.samples < self.batches*2:
            raise ValueError('samples must be >=2*batches')
        if self.samples % self.batches:
            raise ValueError('samples must be divisible by batches for equal-rate batch accounting')
        for key in ('domain_length_D', 'domain_radius_D', 'speed_tail_thermal'):
            if not np.isfinite(getattr(self, key)) or getattr(self, key) <= 0:
                raise ValueError(key+' must be finite and positive')


def ray_exit(position, velocity, length, radius):
    """First boundary crossing of an interior half-cylinder ray."""
    p, v = np.asarray(position), np.asarray(velocity)
    cap = (length-p[0])/v[0] if v[0] > 0 else (-p[0]/v[0] if v[0] < 0 else np.inf)
    a, b = np.dot(v[1:], v[1:]), np.dot(p[1:], v[1:])
    side = (-b+math.sqrt(max(0., b*b+a*(radius**2-np.dot(p[1:], p[1:])))))/a if a > 0 else np.inf
    if cap <= side:
        return max(0., cap), ('downstream' if v[0] > 0 else 'return')
    return max(0., side), 'side'


def _cai_raw(points, source, options):
    """Extend the existing disk moment integral to its second tensor.

    Existing scalar Cai fields are called directly, preserving disabled mode.
    Basis is (x,r,phi), represented at meridian (x,r,0).
    """
    p = np.asarray(points)
    p = np.column_stack((p[:, 0], np.hypot(p[:, 1], p[:, 2]), np.zeros(len(p))))
    n0, a, s0 = source.number_density_m3, source.thermal_speed, source.speed_ratio
    legacy = analytical.moments(*p.T, source.radius_m, s0,
                                n_radial=options.cai_radial, n_azimuthal=options.cai_azimuthal)
    M0 = legacy['n_over_n0']*n0
    U = legacy['U_sqrt_beta0']*a
    M1, M2 = M0[:, None]*U, np.empty((len(p), 3, 3))
    nr, na = analytical._refined_order(source.radius_m, float(p[:, 0].min()),
                                      options.cai_radial, options.cai_azimuthal)
    rho, wr, phi, wp = analytical._disk_quadrature(source.radius_m, nr, na)
    q = np.column_stack((np.zeros(nr*na), (rho[:, None]*np.cos(phi)).ravel(),
                         (rho[:, None]*np.sin(phi)).ravel()))
    weight = (wr[:, None]*wp).ravel()
    for start in range(0, len(p), max(1, min(options.chunk_size, 500000//len(q)))):
        stop = min(len(p), start+max(1, min(options.chunk_size, 500000//len(q))))
        delta = p[start:stop, None, :]-q
        radius = np.linalg.norm(delta, axis=-1)
        e = delta/radius[..., None]
        s = s0*e[..., 0]
        decay = np.exp(-s0*s0*(1-e[..., 0]**2))*(1+erf(s))
        fourth = ((.5*s**3+1.25*s)*math.exp(-s0*s0)
                  +math.sqrt(math.pi)*decay*(.375+1.5*s*s+.5*s**4))
        weights = (n0*a*a*math.pi**-1.5*e[..., 0]/radius**2*weight*fourth)
        M2[start:stop] = np.einsum('ij,ijk,ijl->ikl', weights, e, e)
    out = derive(M0, M1, M2, source.species.mass_kg)
    # These are the EXACT existing evaluator's arrays, not sampled limits.
    out['U'] = U
    out['T_trans'] = legacy['T_over_T0']*source.temperature_K
    return out


class VirtualCollisionModel:
    def __init__(self, source, kernel, options=None):
        from .source import DiskSource
        from .kernel import VHSKernel, ZeroKernel, UnsupportedCollisionLaw
        if type(source) is not DiskSource:
            raise ValueError('unsupported geometry/source: uniform DiskSource required')
        if type(kernel) not in (VHSKernel, ZeroKernel):
            raise UnsupportedCollisionLaw('only the verified VHSKernel or explicit ZeroKernel is supported')
        self.source, self.kernel = source, kernel
        self.options = options or NumericalOptions()
        if hasattr(kernel, 'species') and kernel.species != source.species:
            raise ValueError('source and kernel species differ')
        self.length = 2*source.radius_m*self.options.domain_length_D
        self.radius = 2*source.radius_m*self.options.domain_radius_D
        if self.radius <= source.radius_m:
            raise ValueError('control cylinder must contain emitting disk')
        metadata = dict(version=MODEL_VERSION, source=asdict(source), options=asdict(self.options),
                        kernel=asdict(kernel) if hasattr(kernel, '__dataclass_fields__') else repr(kernel))
        self.cache_key = hashlib.sha256(json.dumps(metadata, sort_keys=True).encode()).hexdigest()
        self.events = None
        self.diagnostics = dict(model_version=MODEL_VERSION, cache_key=self.cache_key,
                                coordinate_basis=['axial', 'cylindrical_radial', 'azimuthal'],
                                domain_length_m=self.length, domain_radius_m=self.radius,
                                domain_tail_error='unavailable; enlarge domain to estimate',
                                quadrature_error='unavailable; refine quadrature to estimate',
                                dsmc_uncertainty='unavailable', cache_state='cold')
        self.diagnostics['validity_scope'] = ('valid arrays check moment realizability only; '
            'quadrature, domain, sampling convergence and single-event adequacy require separate checks')
        b = self.options.speed_tail_thermal
        self.diagnostics['speed_cutoff_m_s'] = source.thermal_speed*(source.speed_ratio+b)
        self.diagnostics['omitted_ballistic_density_upper_bound_over_exit_density'] = float(
            erfc(b)+2*b/math.sqrt(math.pi)*math.exp(-b*b))
        self.diagnostics['speed_tail_scope'] = ('absolute f0 density tail bound from centered Maxwellian; '
            'not a bound on collision-frequency or corrected-field quadrature error')

    def _frequency(self, point, velocity):
        partners, density = self.source.local_quadrature(point, self.options)
        rates = density*self.kernel.rate(np.linalg.norm(partners-velocity, axis=-1))
        return float(rates.sum()), partners, rates

    def _hazard(self, position, velocity, end_time):
        # Quadratically spaced time panels resolve the dense near-exit plume
        # without adding a fitted collision distance. Refine hazard_steps to
        # assess this midpoint rule independently of sampling uncertainty.
        edges = end_time*np.linspace(0., 1., self.options.hazard_steps+1)**2
        widths = np.diff(edges)
        times = .5*(edges[1:]+edges[:-1])
        increments = np.array([self._frequency(position+t*velocity, velocity)[0]*dt
                               for t, dt in zip(times, widths)])
        return increments, edges

    def preprocess(self):
        if self.events is not None:
            self.diagnostics['cache_state'] = 'warm'
            return self.diagnostics
        started = time.perf_counter()
        self.events = []
        exposures, probabilities, remaining, weights = [], [], [], []
        boundary = dict(downstream=0., side=0., return_=0.)
        pair_residual = np.zeros(5)
        count_by_batch = np.full(self.options.batches, self.options.samples//self.options.batches)
        count_by_batch[:self.options.samples % self.options.batches] += 1
        from .kernel import elastic_pair
        for batch, count in enumerate(count_by_batch):
            seed = self.options.seed+batch*104729
            positions, velocities = self.source.sample_flux(int(count), seed)
            rng = np.random.default_rng(seed+17)
            rays = []
            for p, v in zip(positions, velocities):
                end, crossing = ray_exit(p, v, self.length, self.radius)
                boundary[crossing if crossing != 'return' else 'return_'] += 1/self.options.samples
                inc, edges = self._hazard(p, v, end)
                tau = float(inc.sum())
                probability = -math.expm1(-tau)
                exposures.append(tau)
                probabilities.append(probability)
                if probability == 0:
                    continue
                target = -math.log1p(-rng.random()*probability)
                cumul = np.r_[0., np.cumsum(inc)]
                index = min(np.searchsorted(cumul, target, side='right')-1, len(inc)-1)
                event_time = edges[index]+(target-cumul[index])/inc[index]*(edges[index+1]-edges[index])
                location = p+event_time*v
                nu, partners, rates = self._frequency(location, v)
                if nu <= 0 or not np.isfinite(nu):
                    raise FloatingPointError('invalid event collision frequency')
                w = partners[rng.choice(len(rates), p=rates/nu)]
                backward_time = location[0]/w[0]
                q = location-backward_time*w
                partner_hazard = self._hazard(q, w, backward_time)[0].sum()
                weight = .5*probability*math.exp(-partner_hazard)
                direction = rng.normal(size=3)
                direction /= np.linalg.norm(direction)
                vp, wp = elastic_pair(v, w, direction)
                this = [(location, vp, weight), (location, wp, weight),
                        (location, v, -weight), (location, w, -weight)]
                rays.extend(this)
                weights.append(weight)
                residual = np.r_[0., vp+wp-v-w,
                                  np.dot(vp, vp)+np.dot(wp, wp)-np.dot(v, v)-np.dot(w, w)]
                pair_residual += weight/self.options.samples*residual
                for vv in (vp, wp):
                    after, _ = ray_exit(location, vv, self.length, self.radius)
                    remaining.append(float(self._hazard(location, vv, after)[0].sum()))
                for loc, vv, ww in this:
                    _, crossing = ray_exit(loc, vv, self.length, self.radius)
                    boundary[crossing if crossing != 'return' else 'return_'] += ww/self.options.samples
            # A batch estimates the full rate independently; average batches below.
            self.events.append([(p, v, w*self.source.emission_rate/count) for p, v, w in rays])
        exp = np.asarray(exposures)
        rem = np.asarray(remaining)
        self.diagnostics.update(preprocessing_wall_s=time.perf_counter()-started,
            samples=self.options.samples, batches=self.options.batches,
            exposure=exp.tolist(), event_probability=probabilities,
            uncollided_weight=(np.exp(-exp)).tolist(), symmetric_pair_weights=weights,
            probability_scope='event_probability and uncollided_weight are frozen first-event proposal weights, before mutual thinning',
            ensemble_event_rate_over_emission=float(np.sum(weights)/self.options.samples),
            ensemble_uncollided_source_number_fraction=1.-2*float(np.sum(weights)/self.options.samples),
            remaining_exposure=rem.tolist(),
            exposure_mean=float(exp.mean()), remaining_exposure_mean=float(rem.mean()) if rem.size else 0.,
            remaining_exposure_weighted_mean=(float(np.average(rem, weights=np.repeat(weights, 2)))
                if rem.size and np.sum(weights) > 0 else 0.),
            multiple_event_warning=bool(np.any(rem > 1.) or np.any(exp > 1.)),
            closed_control_volume_flux=dict(number_relative_residual=float(sum(boundary.values())-1),
                momentum_residual_m_s=pair_residual[1:4].tolist(),
                twice_energy_residual_m2_s2=float(pair_residual[4]),
                number_residual_per_s=float(sum(boundary.values())-1)*self.source.emission_rate,
                momentum_residual_N=(pair_residual[1:4]*self.source.emission_rate*self.source.species.mass_kg).tolist(),
                energy_residual_W=float(.5*pair_residual[4]*self.source.emission_rate*self.source.species.mass_kg),
                outward_number_by_boundary_over_emission=boundary,
                scope='sampled source plus all signed event rays; all boundary exits included'))
        return self.diagnostics

    def _ray_moments(self, p, v, rate, xb, rb):
        """Exact ring residence intervals, Gauss integration of rotating basis."""
        n = len(xb)
        end, _ = ray_exit(p, v, self.length, self.radius)
        cuts = [np.zeros(n), np.full(n, end)]
        if v[0] != 0:
            cuts.extend([(xb[:, 0]-p[0])/v[0], (xb[:, 1]-p[0])/v[0]])
        a, b, c = np.dot(v[1:], v[1:]), np.dot(p[1:], v[1:]), np.dot(p[1:], p[1:])
        if a > 0:
            for radius in (rb[:, 0], rb[:, 1]):
                disc = b*b-a*(c-radius**2)
                root = np.sqrt(np.maximum(disc, 0.))
                # No real intersection adds repeated endpoints, never fake time.
                cuts.extend([np.where(disc >= 0, (-b-root)/a, 0.),
                             np.where(disc >= 0, (-b+root)/a, 0.)])
        cuts = np.sort(np.clip(np.stack(cuts, axis=1), 0, end), axis=1)
        out0, out1, out2 = np.zeros(n), np.zeros((n, 3)), np.zeros((n, 3, 3))
        gn, gw = np.polynomial.legendre.leggauss(self.options.residence_order)
        for left, right in zip(cuts[:, :-1].T, cuts[:, 1:].T):
            midpoint = p+(left+right)[:, None]*.5*v
            r = np.linalg.norm(midpoint[:, 1:], axis=1)
            inside = ((midpoint[:, 0] >= xb[:, 0]) & (midpoint[:, 0] <= xb[:, 1])
                      & (r >= rb[:, 0]) & (r <= rb[:, 1]) & (right > left))
            active = np.flatnonzero(inside)
            if not len(active):
                continue
            dt = (right-left)[active]
            out0[active] += rate*dt
            for node, weight in zip(gn, gw):
                pos = p+(left[active]+.5*(node+1)*(right-left)[active])[:, None]*v
                radius = np.linalg.norm(pos[:, 1:], axis=1)
                safe = np.where(radius > 1e-30, radius, 1.)
                components = np.column_stack((np.full(len(active), v[0]),
                    (pos[:, 1]*v[1]+pos[:, 2]*v[2])/safe,
                    (pos[:, 1]*v[2]-pos[:, 2]*v[1])/safe))
                out1[active] += rate*.5*weight*dt[:, None]*components
                out2[active] += rate*.5*weight*dt[:, None, None]*components[:, :, None]*components[:, None, :]
        return out0, out1, out2

    def _candidate_bins(self, p, v, xb, rb):
        """Conservative axis/radius envelope of a finite straight ray.

        Radial norm is convex along the segment; its maximum is at an endpoint
        and its minimum is at the clipped quadratic minimum. This is only a
        search prefilter. Exact residence intervals remain in _ray_moments.
        """
        end, _ = ray_exit(p, v, self.length, self.radius)
        finish = p+end*v
        xmin, xmax = min(p[0], finish[0]), max(p[0], finish[0])
        a = np.dot(v[1:], v[1:])
        closest = min(end, max(0., -np.dot(p[1:], v[1:])/a)) if a > 0 else 0.
        rmin = np.linalg.norm(p[1:]+closest*v[1:])
        rmax = max(np.linalg.norm(p[1:]), np.linalg.norm(finish[1:]))
        # Expand the envelope by roundoff only, so a tangential/bin-boundary
        # trajectory cannot be rejected through last-bit coordinate arithmetic.
        tol = 32*np.finfo(float).eps*max(self.length, self.radius)
        return np.flatnonzero((xb[:, 0] <= xmax+tol) & (xb[:, 1] >= xmin-tol)
                             & (rb[:, 0] <= rmax+tol) & (rb[:, 1] >= rmin-tol))

    def evaluate(self, points, *, half_width_m=None, collisions=True):
        points = np.asarray(points, dtype=float).reshape(-1, 3)
        if not len(points) or not np.isfinite(points).all() or np.any(points[:, 0] <= 0):
            raise ValueError('finite points strictly downstream required')
        if not collisions or np.all(self.kernel.rate(np.array([1., 100., 1000.])) == 0):
            out = _cai_raw(points, self.source, self.options)
            # Point API uses the original Cartesian coordinates, exactly as Cai.
            # Ring-bin API below instead returns the local cylindrical basis.
            legacy = analytical.moments(*points.T, self.source.radius_m, self.source.speed_ratio,
                n_radial=self.options.cai_radial, n_azimuthal=self.options.cai_azimuthal)
            rho = np.hypot(points[:, 1], points[:, 2])
            ey = np.divide(points[:, 1], rho, out=np.ones(len(rho)), where=rho > 0)
            ez = np.divide(points[:, 2], rho, out=np.zeros(len(rho)), where=rho > 0)
            rotation = np.zeros((len(points), 3, 3))
            rotation[:, 0, 0] = 1
            rotation[:, 1, 1], rotation[:, 1, 2] = ey, -ez
            rotation[:, 2, 1], rotation[:, 2, 2] = ez, ey
            U = legacy['U_sqrt_beta0']*self.source.thermal_speed
            M0 = legacy['n_over_n0']*self.source.number_density_m3
            M2 = np.einsum('nij,njk,nlk->nil', rotation, out['M2'], rotation)
            out = derive(M0, M0[:, None]*U, M2, self.source.species.mass_kg)
            out['U'], out['T_trans'] = U, legacy['T_over_T0']*self.source.temperature_K
            out.update(standard_error={}, diagnostics=dict(self.diagnostics, correction_disabled=True))
            out['diagnostics']['coordinate_basis'] = ['Cartesian_x', 'Cartesian_y', 'Cartesian_z']
            return out
        if half_width_m is None:
            raise ValueError('enabled point evaluation requires explicit cylindrical bin half_width_m')
        hx, hr = half_width_m
        rho = np.hypot(points[:, 1], points[:, 2])
        return self.evaluate_bins(np.column_stack((points[:, 0]-hx, points[:, 0]+hx)),
                                  np.column_stack((np.maximum(0., rho-hr), rho+hr)))

    def raw_correction_bins(self, x_bounds, r_bounds):
        """Signed transport correction only, without evaluating a Cai grid.

        These are delta raw moments, not a standalone distribution. This adapter
        permits interpolation of the correction followed by addition to the
        exact pointwise Cai raw moments. Independent batch buffers retain
        spatial covariance. No central moments are formed from a signed delta.
        """
        started = time.perf_counter()
        xb, rb = np.asarray(x_bounds, dtype=float), np.asarray(r_bounds, dtype=float)
        if (xb.shape != rb.shape or xb.ndim != 2 or xb.shape[1] != 2 or not len(xb)
            or not np.isfinite(xb).all() or not np.isfinite(rb).all()
            or np.any(xb[:, 0] < 0) or np.any(xb[:, 1] > self.length)
            or np.any(rb[:, 0] < 0) or np.any(rb[:, 1] > self.radius)
            or np.any(np.diff(xb) <= 0) or np.any(np.diff(rb) <= 0)):
            raise ValueError('positive-volume (N,2) bins must lie inside the control cylinder')
        self.preprocess()
        volume = math.pi*np.diff(xb)[:, 0]*(rb[:, 1]**2-rb[:, 0]**2)
        estimates = []
        for rays in self.events:
            delta = [np.zeros(len(xb)), np.zeros((len(xb), 3)), np.zeros((len(xb), 3, 3))]
            for p, v, rate in rays:
                candidates = self._candidate_bins(p, v, xb, rb)
                for start in range(0, len(candidates), self.options.chunk_size):
                    ids = candidates[start:start+self.options.chunk_size]
                    raw = self._ray_moments(p, v, rate, xb[ids], rb[ids])
                    for dest, values in zip(delta, raw):
                        dest[ids] += values
            estimates.append([values/volume.reshape((-1,)+(1,)*(values.ndim-1)) for values in delta])
        batches = dict(zip(('M0', 'M1', 'M2'),
                          [np.stack([batch[j] for batch in estimates]) for j in range(3)]))
        return {**{key: values.mean(axis=0) for key, values in batches.items()},
                'batch_raw_moments': batches, 'bin_volumes_m3': volume,
                'diagnostics': dict(self.diagnostics, correction_evaluation_wall_s=time.perf_counter()-started,
                                    bins_evaluated=len(xb), estimator='signed cylindrical-ring residence correction')}

    def evaluate_bins(self, x_bounds, r_bounds, *, collisions=True):
        started = time.perf_counter()
        xb, rb = np.asarray(x_bounds, dtype=float), np.asarray(r_bounds, dtype=float)
        if xb.shape != rb.shape or xb.ndim != 2 or xb.shape[1] != 2 or not len(xb):
            raise ValueError('bounds must have matching nonempty (N,2) shapes')
        if (not np.isfinite(xb).all() or not np.isfinite(rb).all() or np.any(xb[:, 0] < 0)
            or np.any(xb[:, 1] > self.length) or np.any(rb[:, 0] < 0)
            or np.any(rb[:, 1] > self.radius) or np.any(np.diff(xb) <= 0) or np.any(np.diff(rb) <= 0)):
            raise ValueError('positive-volume bins must lie inside the control cylinder')
        volume = math.pi*np.diff(xb)[:, 0]*(rb[:, 1]**2-rb[:, 0]**2)
        base = [np.zeros(len(xb)), np.zeros((len(xb), 3)), np.zeros((len(xb), 3, 3))]
        nodes, weights = np.polynomial.legendre.leggauss(self.options.bin_order)
        for nx, wx in zip(nodes, weights):
            for nr, wr in zip(nodes, weights):
                x = xb[:, 0]+.5*(nx+1)*np.diff(xb)[:, 0]
                r = np.sqrt(rb[:, 0]**2+.5*(nr+1)*(rb[:, 1]**2-rb[:, 0]**2))
                for start in range(0, len(xb), self.options.chunk_size):
                    stop = start+self.options.chunk_size
                    raw = _cai_raw(np.column_stack((x[start:stop], r[start:stop], np.zeros(len(x[start:stop])))), self.source, self.options)
                    for dest, key in zip(base, ('M0', 'M1', 'M2')):
                        dest[start:stop] += .25*wx*wr*raw[key]
        if collisions:
            self.preprocess()
        estimates = []
        for rays in self.events if collisions else [[]]*self.options.batches:
            delta = [np.zeros_like(a) for a in base]
            for p, v, rate in rays:
                candidates = self._candidate_bins(p, v, xb, rb)
                for start in range(0, len(candidates), self.options.chunk_size):
                    stop = start+self.options.chunk_size
                    indices = candidates[start:stop]
                    raw = self._ray_moments(p, v, rate, xb[indices], rb[indices])
                    for dest, val in zip(delta, raw):
                        dest[indices] += val
            estimates.append([b+d/volume.reshape((-1,)+(1,)*(d.ndim-1)) for b, d in zip(base, delta)])
        arrays = [np.stack([batch[j] for batch in estimates]) for j in range(3)]
        means = [a.mean(axis=0) for a in arrays]
        out = derive(*means, self.source.species.mass_kg)
        batch_derived = [derive(*batch, self.source.species.mass_kg) for batch in estimates]
        errors = {key: np.std(np.stack([b[key] for b in batch_derived]), axis=0, ddof=1)/math.sqrt(len(estimates))
                  for key in ('n', 'U', 'T_tensor', 'T_trans')}
        errors.update({key: a.std(axis=0, ddof=1)/math.sqrt(len(estimates)) for key, a in zip(('M0', 'M1', 'M2'), arrays)})
        # Retain the already-computed independent batch arrays so adapters can
        # propagate spatial covariance through interpolation/aggregation before
        # deriving central moments. These are existing buffers, not copies.
        out.update(standard_error=errors,
                   batch_raw_moments=dict(zip(('M0', 'M1', 'M2'), arrays)),
                   ballistic_raw_moments=dict(zip(('M0', 'M1', 'M2'), base)),
                   bin_volumes_m3=volume,
                   diagnostics=dict(self.diagnostics, evaluation_wall_s=time.perf_counter()-started,
                       bins_evaluated=len(xb), numerical_uncertainty='independent scrambled batch standard error only',
                       invalid_bin_count=int((~out['valid']).sum())))
        return out
