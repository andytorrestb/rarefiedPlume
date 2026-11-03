"""Coordinate-only, empirical deformations of the unchanged Cai field.

``W`` maps physical comparison coordinates to analytical sampling coordinates;
``density = Cai(W)``. Coordinates are nozzle diameters from the exit centre.
No Jacobian or fitted amplitude multiplies density. These maps do not define a
conservative flow solution. The radial/angular family preserves angle, whereas
the two affine families can change it.
"""

from __future__ import annotations

from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import math
import multiprocessing
import os
from pathlib import Path

import numpy as np

from plumetools.cai2012 import analytical


FAMILIES = ("baseline", "axial", "axial_transverse", "radial_angular")
FEATURES = ("inverse_kn", "inverse_kn_squared")
# Cells within 1e-7 diameters share one Cai value: well below the cell size, and
# axisymmetry makes most of a 3D mesh repeat the same (X, P).
COORDINATE_DECIMALS = 7
_UNIQUE_CACHE_SIZE = 64
_unique_cache = OrderedDict()
_workers = 1
_executor = None


def configure_workers(cores):
    """Set the Cai worker processes, never more than this job's allocated CPUs."""
    global _workers, _executor
    if _executor is not None:
        _executor.shutdown()
        _executor = None
    _workers = max(1, min(int(cores), len(os.sched_getaffinity(0))))
    return _workers


def _map_chunks(tasks):
    """A worker lost to the OS raises BrokenProcessPool rather than hanging."""
    global _executor
    if _workers == 1 or len(tasks) == 1:
        return [_analytical_chunk(task) for task in tasks]
    if _executor is None:
        _executor = ProcessPoolExecutor(_workers, mp_context=multiprocessing.get_context("fork"))
    return list(_executor.map(_analytical_chunk, tasks,
                              chunksize=max(1, len(tasks)//(16*_workers))))


def _analytical_chunk(task):
    kind, x, rho, speed_ratio, options = task
    if kind == "density":
        return analytical.density_ratio(x, rho, 0., .5, speed_ratio, **options)
    values = analytical.moments(x, rho, 0., .5, speed_ratio, **options)
    return values["T_over_T0" if kind == "temperature" else "U_sqrt_beta0"]


def unique_points(points):
    """Representative points and the inverse restoring ``points``; cached per array."""
    key = (points.shape, hashlib.blake2b(np.ascontiguousarray(points).data, digest_size=16).digest())
    cached = _unique_cache.get(key)
    if cached is None:
        _, first, inverse = np.unique(np.round(points, COORDINATE_DECIMALS), axis=0,
                                      return_index=True, return_inverse=True)
        cached = first, inverse.ravel()
        _unique_cache[key] = cached
        if len(_unique_cache) > _UNIQUE_CACHE_SIZE:
            _unique_cache.popitem(last=False)
    else:
        _unique_cache.move_to_end(key)
    return points[cached[0]], cached[1]


def _evaluate_sorted(points, kind, speed_ratio, n_radial, n_azimuthal, chunk_size):
    """Evaluate in X order: each chunk's nearest point sets its adaptive disk order."""
    order = np.argsort(points[:, 0], kind="stable")
    options = {"n_radial": n_radial, "n_azimuthal": n_azimuthal}
    tasks = [(kind, *points[order[start:start+chunk_size]].T, speed_ratio, options)
             for start in range(0, len(order), chunk_size)]
    values = np.concatenate(_map_chunks(tasks)) if tasks else np.empty((0,))
    result = np.empty_like(values)
    result[order] = values
    return result


def cai_field(points, kind, speed_ratio, n_radial=48, n_azimuthal=96, chunk_size=128):
    """Unmapped Cai ``density``, ``temperature`` or ``velocity`` at ``[X, P]`` points."""
    unique, inverse = unique_points(np.asarray(points, dtype=float).reshape(-1, 2))
    return _evaluate_sorted(unique, kind, speed_ratio, n_radial, n_azimuthal, chunk_size)[inverse]


def kn_feature(kn, feature="inverse_kn"):
    """Empirical feature, not a collision probability.

    The default is ``1/(1+Kn)``. ``inverse_kn_squared`` is an explicitly
    configured sensitivity basis, ``1/(1+Kn)**2``. Both saturate at low Kn;
    neither can distinguish arbitrarily small Kn. Zero is allowed as the
    limiting endpoint used in bound diagnostics; physical fits use Kn > 0.
    """
    k = np.asarray(kn, dtype=float)
    if np.any(np.isnan(k)) or np.any(k < 0):
        raise ValueError("Kn must be nonnegative (infinity is the identity limit)")
    if feature not in FEATURES:
        raise ValueError(f"unknown Kn feature {feature!r}; choose {FEATURES}")
    value = 1.0 / (1.0 + k)
    return value if feature == "inverse_kn" else value * value


def _points(points, columns=2):
    points = np.asarray(points, dtype=float)
    if points.ndim < 1 or points.shape[-1] != columns:
        raise ValueError(f"points must have a trailing coordinate dimension {columns}")
    if not np.all(np.isfinite(points)):
        raise ValueError("coordinates must be finite")
    if columns == 2 and np.any(points < 0):
        raise ValueError("supported downstream coordinates require X >= 0 and P >= 0")
    return points


def _angular_extrema(coefficients):
    """Exact extrema of c0+c1*cos(2theta)+c2*cos(4theta)."""
    c = np.pad(np.asarray(coefficients, dtype=float), (0, 3-len(coefficients)))
    candidates = [-1.0, 1.0]
    if c[2] != 0:
        vertex = -c[1] / (4*c[2])
        if -1 < vertex < 1:
            candidates.append(vertex)
    u = np.asarray(candidates)
    values = c[0] + c[1]*u + c[2]*(2*u*u-1)
    return float(np.min(values)), float(np.max(values))


class CoordinateMap:
    """Small bounded mapping; coefficients are flattened in serialized order.

    Affine ordering is ``[x_h, x_h2, rho_h, rho_h2]`` with the second-degree
    entries present only for degree 2 and the rho entries only for
    ``axial_transverse``. Radial ordering is ``[c_00,c_01,c_02,c_10,c_11,c_12]``
    truncated to ``angular_order + 1`` entries in each radial row.

    ``log_bound`` bounds each affine log scale, or the radial log derivative,
    globally over all Kn >= 0, distance >= 0, and downstream angles. Invalid
    candidates raise rather than clipping coefficients or dropping samples.
    """

    def __init__(self, family="baseline", coefficients=None, degree=1,
                 angular_order=2, quadrature_order=24, log_bound=2.0,
                 feature="inverse_kn"):
        if family not in FAMILIES:
            raise ValueError(f"unknown map family {family!r}; choose {FAMILIES}")
        if degree not in (1, 2):
            raise ValueError("degree must be 1 or 2")
        if angular_order not in (0, 1, 2):
            raise ValueError("angular_order must be 0, 1 or 2")
        if int(quadrature_order) != quadrature_order or quadrature_order < 2:
            raise ValueError("quadrature_order must be an integer >= 2")
        if not np.isfinite(log_bound) or not 0 < log_bound <= 20:
            raise ValueError("log_bound must be finite and in (0,20]")
        kn_feature(1.0, feature)
        self.family, self.degree = family, int(degree)
        self.angular_order, self.quadrature_order = int(angular_order), int(quadrature_order)
        self.log_bound, self.feature = float(log_bound), feature
        self.coefficients = np.zeros(self.parameter_count) if coefficients is None else np.asarray(coefficients, dtype=float).ravel().copy()
        if len(self.coefficients) != self.parameter_count or not np.all(np.isfinite(self.coefficients)):
            raise ValueError(f"{family} requires {self.parameter_count} finite coefficients")
        if self.bound_diagnostics()["max_abs_log"] > self.log_bound + 1e-12:
            raise ValueError("map exceeds the configured global log_bound")
        self.coefficients.setflags(write=False)

    @property
    def parameter_count(self):
        return {"baseline": 0, "axial": self.degree,
                "axial_transverse": 2*self.degree,
                "radial_angular": 2*(self.angular_order+1)}[self.family]

    def _kn(self, kn):
        h = kn_feature(kn, self.feature)
        if h.ndim != 0:
            raise ValueError("one scalar Kn is required per map evaluation")
        return float(h)

    def _affine_logs(self, kn):
        h = self._kn(kn)
        powers = np.asarray([h**(i+1) for i in range(self.degree)])
        lx = float(self.coefficients[:self.degree] @ powers)
        lr = float(self.coefficients[self.degree:] @ powers) if self.family == "axial_transverse" else 0.0
        return lx, lr

    def scale_factors(self, kn):
        if self.family == "radial_angular":
            raise ValueError("radial/angular map has no single pair of affine scales")
        logs = (0.0, 0.0) if self.family == "baseline" else self._affine_logs(kn)
        return {"axial": math.exp(logs[0]), "transverse": math.exp(logs[1])}

    def _angular(self, theta):
        return np.stack([np.cos(2*l*theta) for l in range(self.angular_order+1)], axis=-1)

    def radial_coordinates(self, xi, theta, kn):
        """Return ``(eta, d_eta/d_xi)``; deterministic Gauss-Legendre integral.

        Integrate in y=log(1+t), where B1=1-exp(-y). This resolves the source
        and distant plume using one small rule across many radial decades.
        """
        if self.family != "radial_angular":
            raise ValueError("radial_coordinates applies only to radial_angular")
        xi, theta = np.broadcast_arrays(np.asarray(xi, float), np.asarray(theta, float))
        if np.any(~np.isfinite(xi)) or np.any(xi < 0) or np.any(~np.isfinite(theta)) or np.any((theta < 0) | (theta > np.pi/2)):
            raise ValueError("finite xi >= 0 and theta in [0,pi/2] required")
        h = self._kn(kn)
        c = self.coefficients.reshape(2, self.angular_order+1)
        a = self._angular(theta)
        intercept, distance = h*(a @ c[0]), h*(a @ c[1])
        derivative = np.exp(intercept + distance*xi/(1+xi))
        if h == 0 or not np.any(self.coefficients):
            return xi.copy(), np.ones_like(xi)
        if not np.any(c[1]):
            return xi*np.exp(intercept), derivative
        nodes, weights = np.polynomial.legendre.leggauss(self.quadrature_order)
        upper = np.log1p(xi).reshape(-1)
        intercept_flat, distance_flat = intercept.reshape(-1), distance.reshape(-1)
        eta = np.empty_like(upper)
        # Bound temporary arrays independently of the number of mesh cells.
        # Three main work arrays at 1e6 float64 pairs use about 24 MiB.
        chunk = max(1, 1_000_000 // self.quadrature_order)
        for start in range(0, upper.size, chunk):
            stop = min(start+chunk, upper.size)
            y = 0.5*upper[start:stop, None]*(nodes+1)
            integrand = np.exp(intercept_flat[start:stop, None] +
                               distance_flat[start:stop, None]*(-np.expm1(-y)) + y)
            eta[start:stop] = 0.5*upper[start:stop]*np.sum(weights*integrand, axis=-1)
        return eta.reshape(xi.shape), derivative

    def map_points(self, points, kn):
        points = _points(points)
        h = self._kn(kn)
        if self.family == "baseline" or h == 0 or not np.any(self.coefficients):
            return points.copy()
        if self.family in ("axial", "axial_transverse"):
            lx, lr = self._affine_logs(kn)
            return points*np.exp(-np.asarray([lx, lr]))
        xi = np.hypot(points[..., 0], points[..., 1])
        eta, _ = self.radial_coordinates(xi, np.arctan2(points[..., 1], points[..., 0]), kn)
        scale = np.divide(eta, xi, out=np.ones_like(xi), where=xi > 0)
        return points*scale[..., None]

    def map_cartesian(self, points, kn, *, axis=0, origin=(0., 0., 0.), diameter=1.0):
        """Map physical 3-D coordinates preserving azimuth about selected axis.

        ``axis`` is the already adapted solver downstream axis. The caller
        obtains it from the case's existing geometry adapter.
        """
        p = _points(points, 3)
        if axis not in (0, 1, 2) or not np.isfinite(diameter) or diameter <= 0:
            raise ValueError("axis must be 0,1,2 and diameter positive")
        origin = np.asarray(origin, dtype=float)
        if origin.shape != (3,) or not np.all(np.isfinite(origin)):
            raise ValueError("origin must contain three finite coordinates")
        q = (p-origin)/diameter
        transverse = [i for i in range(3) if i != axis]
        rho = np.hypot(q[..., transverse[0]], q[..., transverse[1]])
        meridional = np.stack((q[..., axis], rho), axis=-1)
        mapped = self.map_points(meridional, kn)
        if np.array_equal(mapped, meridional):
            return p.copy()
        ratio = np.divide(mapped[..., 1], rho, out=np.ones_like(rho), where=rho > 0)
        result = q.copy()
        result[..., axis] = mapped[..., 0]
        result[..., transverse] *= ratio[..., None]
        return result*diameter+origin

    def jacobian(self, points, kn):
        """3-D volume Jacobian diagnostic, never a density multiplier.

        At a radial origin with nonzero B0 angular terms no unique derivative
        exists: return NaN explicitly rather than choosing a directional limit.
        """
        points = _points(points)
        if self.family == "baseline" or self._kn(kn) == 0 or not np.any(self.coefficients):
            return np.ones(points.shape[:-1])
        if self.family in ("axial", "axial_transverse"):
            lx, lr = self._affine_logs(kn)
            return np.full(points.shape[:-1], math.exp(-lx-2*lr))
        xi = np.hypot(points[..., 0], points[..., 1])
        eta, derivative = self.radial_coordinates(xi, np.arctan2(points[..., 1], points[..., 0]), kn)
        ratio = np.divide(eta, xi, out=np.ones_like(xi), where=xi > 0)
        result = ratio*ratio*derivative
        c = self.coefficients.reshape(2, self.angular_order+1)
        origin_value = math.nan if np.any(c[0, 1:] != 0) else math.exp(3*self._kn(kn)*c[0, 0])
        return np.where(xi == 0, origin_value, result)

    def bound_diagnostics(self, kn=0.0, tolerance=1e-5):
        h = self._kn(kn)
        values = [0.0]
        if self.family in ("axial", "axial_transverse"):
            for coefficients in self.coefficients.reshape(-1, self.degree):
                candidates = [0., h]
                if self.degree == 2 and coefficients[1] != 0:
                    vertex = -coefficients[0]/(2*coefficients[1])
                    if 0 < vertex < h:
                        candidates.append(vertex)
                values.extend(sum(c*t**(i+1) for i, c in enumerate(coefficients)) for t in candidates)
        elif self.family == "radial_angular":
            c = self.coefficients.reshape(2, self.angular_order+1)
            for row in (c[0], c[0]+c[1]):
                values.extend(h*v for v in _angular_extrema(row))
        maximum = float(np.max(np.abs(values)))
        return {"max_abs_log": maximum, "minimum_log": float(min(values)),
                "maximum_log": float(max(values)), "log_bound": self.log_bound,
                "bound_active": maximum >= self.log_bound-tolerance}

    def source_diagnostics(self, kn, near_x=1e-3):
        """Source disk/lip geometry; passing the exit plane is insufficient."""
        points = np.asarray([[0., .5], [near_x, .5], [near_x, 0.]])
        mapped = self.map_points(points, kn)
        angular_origin = self.family == "radial_angular" and self._kn(kn) != 0 and np.any(self.coefficients.reshape(2, self.angular_order+1)[0, 1:])
        return {"exit_lip_physical": points[0].tolist(), "exit_lip_mapped": mapped[0].tolist(),
                "lip_radial_displacement_over_D": float(mapped[0, 1]-.5),
                "emitting_disk_boundary_preserved": bool(abs(mapped[0, 1]-.5) < 1e-10),
                "near_exit_points": points[1:].tolist(), "near_exit_mapped": mapped[1:].tolist(),
                "origin_differentiable": not bool(angular_origin),
                "source_plane_analytical_boundary": "undefined by the off-axis disk quadrature; no source-plane density claim"}

    def quadrature_convergence(self, points, kn):
        """Compare q, 2q and 4q integration rules at exactly the same points."""
        original = self.map_points(points, kn)
        results = []
        for factor in (2, 4):
            configuration = self.to_dict()
            configuration["quadrature_order"] = factor*self.quadrature_order
            results.append(CoordinateMap.from_dict(configuration).map_points(points, kn))
        scale = np.maximum(np.linalg.norm(results[-1], axis=-1), 1e-15)
        return {"orders": [self.quadrature_order, 2*self.quadrature_order, 4*self.quadrature_order],
                "q_vs_4q_max_relative": float(np.max(np.linalg.norm(original-results[-1], axis=-1)/scale, initial=0)),
                "2q_vs_4q_max_relative": float(np.max(np.linalg.norm(results[0]-results[-1], axis=-1)/scale, initial=0))}

    def to_dict(self):
        return {"family": self.family, "coefficients": self.coefficients.tolist(),
                "degree": self.degree, "angular_order": self.angular_order,
                "quadrature_order": self.quadrature_order, "log_bound": self.log_bound,
                "feature": self.feature}

    @classmethod
    def from_dict(cls, configuration):
        return cls(**configuration)


class CaiPredictor:
    """Portable direct evaluator without DSMC inputs, resampling, or amplitude.

    Input ``points[...,2]`` are ``[x/D_e,rho/D_e]``. ``temperature`` transfers
    the same density-derived map, without refitting or changing its coefficients.
    Source-plane points have no defined off-axis disk integral and raise; they
    must be counted by the fixed preprocessing mask, never model-dependently
    discarded. Disk quadrature is the repository's existing adaptive rule.
    """

    def __init__(self, mapping, speed_ratio=2.0, n_radial=48, n_azimuthal=96,
                 chunk_size=128):
        if not isinstance(mapping, CoordinateMap):
            raise TypeError("mapping must be a CoordinateMap")
        if not np.isfinite(speed_ratio) or speed_ratio < 0:
            raise ValueError("speed_ratio must be finite and nonnegative")
        if int(n_radial) != n_radial or n_radial < 2 or int(n_azimuthal) != n_azimuthal or n_azimuthal < 4:
            raise ValueError("disk quadrature needs n_radial >=2, n_azimuthal >=4")
        if int(chunk_size) != chunk_size or chunk_size < 1:
            raise ValueError("chunk_size must be a positive integer")
        self.mapping, self.speed_ratio = mapping, float(speed_ratio)
        self.n_radial, self.n_azimuthal, self.chunk_size = int(n_radial), int(n_azimuthal), int(chunk_size)

    def _evaluate(self, points, kn, temperature):
        points = _points(points)
        shape = points.shape[:-1]
        unique, inverse = unique_points(points.reshape(-1, 2))
        flat = self.mapping.map_points(unique, kn).reshape(-1, 2)
        if np.any(flat[:, 0] <= 0):
            raise ValueError("Cai direct evaluation requires X > 0; source-plane boundary is undefined")
        result = _evaluate_sorted(flat, "temperature" if temperature else "density",
                                  self.speed_ratio, self.n_radial, self.n_azimuthal, self.chunk_size)
        if not np.all(np.isfinite(result)):
            raise ValueError("mapped Cai field is unavailable/nonfinite; candidate must be penalized")
        return result[inverse].reshape(shape)

    def density(self, points, kn):
        return self._evaluate(points, kn, False)

    def temperature(self, points, kn):
        return self._evaluate(points, kn, True)

    def evaluation_diagnostics(self, points, kn):
        """Report the unchanged evaluator's near-plane quadrature limitations.

        A finite numerical value at its maximum disk rule is not proof that a
        near-lip peak is resolved. Flag points whose requested adaptive order
        would exceed the existing evaluator's caps; retain them in field fits.
        """
        mapped = self.mapping.map_points(points, kn).reshape(-1, 2)
        x = mapped[:, 0]
        positive = x > 0
        radial_need = np.full(x.shape, np.inf)
        azimuthal_need = np.full(x.shape, np.inf)
        radial_need[positive] = np.ceil(math.pi/x[positive])
        azimuthal_need[positive] = np.ceil(2*math.pi/x[positive])
        capped = positive & ((radial_need > analytical.MAX_RADIAL) |
                             (azimuthal_need > analytical.MAX_AZIMUTHAL))
        return {"n_points": len(x), "n_undefined_source_plane": int(np.count_nonzero(~positive)),
                "n_adaptive_disk_order_capped": int(np.count_nonzero(capped)),
                "adaptive_disk_order_caps": [analytical.MAX_RADIAL, analytical.MAX_AZIMUTHAL],
                "minimum_mapped_X": float(np.min(x)) if len(x) else None,
                "direct_evaluation": True, "interpolation_cache": False,
                "cap_interpretation": "finite values are retained; convergence near the source/lip is not established by hitting the cap"}

    def quadrature_convergence(self, points, kn, *, temperature=False):
        """Disk-rule sensitivity, separate from the map integration check."""
        results = []
        for factor in (1, 2, 4):
            predictor = CaiPredictor(self.mapping, self.speed_ratio,
                                     min(factor*self.n_radial, analytical.MAX_RADIAL),
                                     min(factor*self.n_azimuthal, analytical.MAX_AZIMUTHAL), self.chunk_size)
            results.append(predictor._evaluate(points, kn, temperature))
        scale = np.maximum(np.abs(results[-1]), 1e-15)
        return {"quantity": "temperature" if temperature else "density",
                "requested_order_factors": [1, 2, 4],
                "base_disk_order": [self.n_radial, self.n_azimuthal],
                "q_vs_4q_max_relative": float(np.max(np.abs(results[0]-results[-1])/scale, initial=0)),
                "2q_vs_4q_max_relative": float(np.max(np.abs(results[1]-results[-1])/scale, initial=0)),
                "evaluation": self.evaluation_diagnostics(points, kn)}

    def to_dict(self):
        return {"schema": "cai-coordinate-predictor-v1", "mapping": self.mapping.to_dict(),
                "speed_ratio": self.speed_ratio, "n_radial": self.n_radial,
                "n_azimuthal": self.n_azimuthal, "chunk_size": self.chunk_size,
                "coordinates": "X=x/D_e,P=rho/D_e; W maps physical to Cai sampling point",
                "density_definition": "Cai(W); no amplitude and no determinant multiplier"}

    def save(self, path):
        Path(path).write_text(json.dumps(self.to_dict(), indent=2)+"\n", encoding="utf-8")

    @classmethod
    def from_dict(cls, data):
        if data.get("schema") != "cai-coordinate-predictor-v1":
            raise ValueError("unsupported predictor schema")
        return cls(CoordinateMap.from_dict(data["mapping"]), **{key: data[key] for key in
                   ("speed_ratio", "n_radial", "n_azimuthal", "chunk_size")})

    @classmethod
    def load(cls, path):
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))
