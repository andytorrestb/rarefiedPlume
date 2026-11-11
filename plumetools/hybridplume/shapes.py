"""Fitted interface shapes: spheroid and capsule (cylinder plus spheroidal tip).

Both shapes are axisymmetric about the x-axis and start at the exit plane x=0.
`q` is a scaled quadratic form that equals 1 on the surface and is < 1 inside.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .ellipsoid import clipped_volume, verify_ellipsoid


def _rho2(points):
    p = np.asarray(points, float)
    return p[:, 1]**2+p[:, 2]**2


@dataclass(frozen=True)
class Spheroid:
    xc: float
    a: float
    b: float

    def q(self, points):
        p = np.asarray(points, float)
        return ((p[:, 0]-self.xc)/self.a)**2+_rho2(p)/self.b**2

    def deviation(self, vertices):
        return float(np.max(np.abs(np.sqrt(self.q(vertices))-1)))

    def outward_gradient(self, points):
        p = np.asarray(points, float)
        return np.column_stack((2*(p[:, 0]-self.xc)/self.a**2,
                                2*p[:, 1]/self.b**2, 2*p[:, 2]/self.b**2))

    def verify(self, lower, upper, retained, clearance):
        return verify_ellipsoid({"xc": self.xc, "a": self.a, "b": self.b,
                                 "clearance_m": clearance}, lower, upper, retained)

    def grown(self, factor):
        return Spheroid(self.xc, self.a*factor, self.b*factor)

    def section(self, n):
        theta0 = math.acos(min(1.0, max(-1.0, -self.xc/self.a)))
        t = np.linspace(theta0, 0.0, n)
        x, y = self.xc+self.a*np.cos(t), self.b*np.sin(t)
        return np.vstack((np.column_stack((x, y)), np.column_stack((x[::-1], -y[::-1]))))

    def volume(self):
        return clipped_volume([self.xc, self.a, self.b])

    def base_radius(self):
        return float(self.b*math.sqrt(max(0.0, 1-(self.xc/self.a)**2)))

    def params(self):
        return [self.xc, self.a, self.b]

    def to_json(self):
        return {"shape": "spheroid", "xc": self.xc, "a": self.a, "b": self.b}


@dataclass(frozen=True)
class Capsule:
    R: float
    L: float
    a: float

    def q(self, points):
        p = np.asarray(points, float)
        tip = np.maximum(p[:, 0]-self.L, 0)/self.a
        return tip**2+_rho2(p)/self.R**2

    def deviation(self, vertices):
        return float(np.max(np.abs(np.sqrt(self.q(vertices))-1)))

    def outward_gradient(self, points):
        p = np.asarray(points, float)
        return np.column_stack((2*np.maximum(p[:, 0]-self.L, 0)/self.a**2,
                                2*p[:, 1]/self.R**2, 2*p[:, 2]/self.R**2))

    def verify(self, lower, upper, retained, clearance):
        margin = float(clearance)
        bad, checked, minimum = 0, 0, float("inf")
        tip = {"xc": self.L, "a": self.a, "b": self.R, "clearance_m": margin}
        for start in range(0, len(lower), 65536):
            excluded = ~np.asarray(retained[start:start+65536], bool)
            lo = np.asarray(lower[start:start+65536])[excluded]-margin
            hi = np.asarray(upper[start:start+65536])[excluded]+margin
            checked += len(lo)
            if not len(lo):
                continue
            rho = np.hypot(*np.maximum(np.maximum(lo[:, 1:], -hi[:, 1:]), 0).T)
            overlap = (lo[:, 0] <= self.L) & (hi[:, 0] >= 0)
            scaled = (rho/self.R)**2
            bad += int(np.count_nonzero(overlap & (scaled < 1-1e-8)))
            if overlap.any():
                minimum = min(minimum, float(scaled[overlap].min()))
            beyond = hi[:, 0] > self.L
            if beyond.any():
                sub = np.asarray(lower[start:start+65536])[excluded][beyond]
                sup = np.asarray(upper[start:start+65536])[excluded][beyond]
                result = verify_ellipsoid(tip, sub, sup, np.zeros(len(sub), bool))
                bad += result["intersected_excluded_voxels"]
                if result["min_scaled_distance_squared"] is not None:
                    minimum = min(minimum, result["min_scaled_distance_squared"])
        return {"feasible": bad == 0, "intersected_excluded_voxels": bad,
                "excluded_voxels_checked": checked,
                "min_scaled_distance_squared": minimum if np.isfinite(minimum) else None}

    def grown(self, factor):
        return Capsule(self.R*factor, self.L, self.a*factor)

    def section(self, n):
        t = np.linspace(math.pi/2, 0.0, n)
        x, y = self.L+self.a*np.cos(t), self.R*np.sin(t)
        if self.L > 0:
            x, y = np.r_[0.0, x], np.r_[self.R, y]
        return np.vstack((np.column_stack((x, y)), np.column_stack((x[::-1], -y[::-1]))))

    def volume(self):
        return math.pi*self.R**2*self.L+2/3*math.pi*self.R**2*self.a

    def base_radius(self):
        return self.R

    def params(self):
        return [self.R, self.L, self.a]

    def to_json(self):
        return {"shape": "capsule", "R": self.R, "L": self.L, "a": self.a}


def from_json(doc):
    kind = doc.get("shape", "spheroid")
    if kind == "spheroid":
        return Spheroid(float(doc["xc"]), float(doc["a"]), float(doc["b"]))
    if kind == "capsule":
        return Capsule(float(doc["R"]), float(doc["L"]), float(doc["a"]))
    raise ValueError(f"unknown fitted shape {kind!r}")
