r"""Moments and boundary crossings measured from the written DSMC parcels.

Standard ``dsmcFoam`` writes, at every write time, each parcel's position and
cell (``positions``: ``(x y z) celli``), velocity ``U``, ``typeId``, and its
identity ``(origProcId, origId)``. Two measurements are built on that.

Tensor moments
--------------
For each cell of volume ``V`` and a set of snapshots, with parcel weight ``W``::

    N  = <sum W>,   P1 = <sum W v>,   P2 = <sum W v v>
    n  = N / V,     U  = P1 / N,      C  = P2 / N - U U
    T  = tr(C) / (3R),   T_ii = C_ii / R,
    T_par = (u.C.u)/R,   T_perp = (tr C - u.C.u) / (2R),   u = U/|U|

``<>`` is the mean over snapshots (instantaneous samples at write times; they
are fewer than ``fieldAverage``'s per-step samples, but they carry the full
tensor, which ``fieldAverage`` of the standard fields does not). Uncertainty is
the standard error over independent blocks of snapshots. Occupancy -- parcels
per cell -- is kept separately from ``n``. A cell whose thermal spread is not
resolved is flagged, never clamped.

Exact crossings (the "crossing window")
---------------------------------------
``DSMCCloud::evolve`` (v2512) runs ``inflow -> move -> collisions``. A parcel
written at the end of step ``k`` at ``x_k`` with velocity ``v_k`` (post-
collision) moves along ``x_k + v_k t``, ``0 <= t <= dt``, during step ``k+1``.
If it is absent at the end of step ``k+1`` it was deleted at a boundary on that
ray. Writing the cloud every step for a window therefore gives every exit
**exactly**: which boundary (the first one the ray meets), and the physical
mass ``W m``, momentum ``W m v_k`` and energy ``W m |v_k|^2 / 2`` it carried.
Newly injected parcels are counted from the inflow model's log line (exact
number); their momentum and energy are not reconstructable from the post-
collision written velocity, so the injection's momentum/energy is reported as
the expected half-range value and labelled so.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_HEADER = re.compile(rb"^\s*(\d+)\s*$", re.M)


def _list_body(path: Path) -> tuple:
    data = Path(path).read_bytes()
    brace = data.find(b"}")
    m = _HEADER.search(data, brace if brace >= 0 else 0)
    if m is None:
        # OpenFOAM writes "0()" for an empty cloud on one line
        if re.search(rb"\n\s*0\s*\(\s*\)", data):
            return 0, b""
        raise ValueError(f"{path}: no list header")
    n = int(m.group(1))
    start = data.find(b"(", m.end())
    end = data.rfind(b")")
    return n, data[start + 1:end]


def _numbers(body: bytes, n: int, per: int, path: Path) -> np.ndarray:
    if n == 0:
        return np.empty((0, per))
    clean = body.replace(b"(", b" ").replace(b")", b" ")
    arr = np.array(clean.split(), dtype=np.float64)
    if arr.size != n * per:
        raise ValueError(f"{path}: expected {n} x {per} numbers, found {arr.size} "
                         "(binary format or barycentric 'coordinates' are not supported)")
    return arr.reshape(n, per)


@dataclass(frozen=True)
class Cloud:
    position: np.ndarray     # (p, 3)
    cell: np.ndarray         # (p,)
    velocity: np.ndarray     # (p, 3)
    type_id: np.ndarray
    ident: np.ndarray        # (p,) int64: origProc << 32 | origId

    @property
    def size(self) -> int:
        return int(self.cell.size)


_UNIFORM = re.compile(rb"^\s*(\d+)\s*\{\s*([^}]*?)\s*\}", re.M)


def _read_list(path: Path, per: int, expect: int | None = None) -> np.ndarray:
    data = Path(path).read_bytes()
    brace = data.find(b"}")
    m = _UNIFORM.search(data, brace + 1 if brace >= 0 else 0)
    if m is not None:     # OpenFOAM's compact all-equal form: N{value}
        n = int(m.group(1))
        value = np.array(m.group(2).replace(b"(", b" ").replace(b")", b" ").split(),
                         dtype=np.float64)
        if expect is not None and n != expect:
            raise ValueError(f"{path}: {n} entries, expected {expect}")
        return np.tile(value.reshape(1, -1), (n, 1)).reshape(n, per)
    n, body = _list_body(path)
    if expect is not None and n != expect:
        raise ValueError(f"{path}: {n} entries, expected {expect}")
    return _numbers(body, n, per, path)


def read_cloud(time_dir: Path, cloud: str = "dsmc") -> Cloud:
    d = Path(time_dir) / "lagrangian" / cloud
    if not (d / "positions").is_file():
        return Cloud(np.empty((0, 3)), np.empty(0, np.int64), np.empty((0, 3)),
                     np.empty(0, np.int64), np.empty(0, np.int64))
    pos = _read_list(d / "positions", 4)
    n = pos.shape[0]
    vel = _read_list(d / "U", 3, n)
    tid = _read_list(d / "typeId", 1, n).ravel() if (d / "typeId").is_file() else np.zeros(n)
    proc = (_read_list(d / "origProcId", 1, n).ravel() if (d / "origProcId").is_file()
            else np.zeros(n))
    oid = _read_list(d / "origId", 1, n).ravel() if (d / "origId").is_file() else np.arange(n)
    ident = (proc.astype(np.int64) << 32) | oid.astype(np.int64)
    return Cloud(pos[:, :3], pos[:, 3].astype(np.int64), vel, tid.astype(np.int64), ident)


@dataclass
class MomentAccumulator:
    n_cells: int

    def __post_init__(self):
        self.snapshots = 0
        self.count = np.zeros(self.n_cells)
        self.w = np.zeros(self.n_cells)
        self.p1 = np.zeros((self.n_cells, 3))
        self.p2 = np.zeros((self.n_cells, 3, 3))

    def add(self, cloud: Cloud, weight: float) -> None:
        c = cloud.cell
        n = self.n_cells
        self.snapshots += 1
        self.count += np.bincount(c, minlength=n)
        self.w += weight * np.bincount(c, minlength=n)
        for i in range(3):
            self.p1[:, i] += weight * np.bincount(c, cloud.velocity[:, i], minlength=n)
            for j in range(i, 3):
                v = weight * np.bincount(c, cloud.velocity[:, i] * cloud.velocity[:, j],
                                         minlength=n)
                self.p2[:, i, j] += v
                if j != i:
                    self.p2[:, j, i] += v

    def fields(self, volumes: np.ndarray, gas, min_parcels: float = 10.0) -> dict:
        k = max(self.snapshots, 1)
        N = self.w / k
        with np.errstate(divide="ignore", invalid="ignore"):
            n = N / volumes
            u = self.p1 / k / N[:, None]
            cov = self.p2 / k / N[:, None, None] - np.einsum("ci,cj->cij", u, u)
            r = gas.gas_constant
            tr = np.trace(cov, axis1=1, axis2=2)
            speed = np.linalg.norm(u, axis=1)
            uhat = u / speed[:, None]
            t_par = np.einsum("ci,cij,cj->c", uhat, cov, uhat) / r
            t_perp = (tr - t_par * r) / (2 * r)
            occupancy = self.count / k
            resolved = (occupancy >= min_parcels) & (tr > 0)
        return {"n": n, "rho": gas.mass_kg * n, "U": u, "C": cov, "T": tr / (3 * r),
                "T_dir": np.diagonal(cov, axis1=1, axis2=2) / r, "T_parallel": t_par,
                "T_perpendicular": t_perp, "parcels_per_snapshot": occupancy,
                "resolved": resolved, "snapshots": self.snapshots}


def block_statistics(per_block: list) -> dict:
    """Mean and standard error over independent blocks of a list of arrays."""
    arr = np.stack(per_block)
    b = arr.shape[0]
    mean = np.nanmean(arr, axis=0)
    se = np.nanstd(arr, axis=0, ddof=1) / np.sqrt(b) if b > 1 else np.full_like(mean, np.nan)
    return {"mean": mean, "stderr": se, "blocks": b}


# --------------------------------------------------------------------------- #
# exact crossings
# --------------------------------------------------------------------------- #

def classify_exits(x0: np.ndarray, v: np.ndarray, dt: float, box: tuple,
                   in_retained, radius_m: float, n_samples: int = 8) -> np.ndarray:
    """Which boundary each straight path ``x0 + v t`` (``0 <= t <= dt``) leaves through.

    Returns labels: 0 interface (enters the retained continuum region),
    1 source (``x = 0``, ``r <= Re``), 2 exitPlane (``x = 0``, ``r > Re``),
    3 farField, -1 still inside (should not happen for a deleted parcel).
    The path is sampled at ``n_samples`` points; the first point outside the
    box or inside the retained region decides. ``in_retained(points) -> bool``.
    """
    lo, hi = np.asarray(box[0], float), np.asarray(box[1], float)
    labels = np.full(x0.shape[0], -1, dtype=np.int64)
    undecided = np.ones(x0.shape[0], bool)
    for s in np.linspace(1.0 / n_samples, 1.0, n_samples):
        p = x0 + v * (s * dt)
        out = np.any((p < lo) | (p > hi), axis=1)
        idx = np.nonzero(undecided & out)[0]
        if idx.size:
            # time to each box face along the ray; the earliest exit wins
            t_face = np.full((idx.size, 6), np.inf)
            for a in range(3):
                va = v[idx, a]
                with np.errstate(divide="ignore", invalid="ignore"):
                    t_lo = (lo[a] - x0[idx, a]) / va
                    t_hi = (hi[a] - x0[idx, a]) / va
                t_face[:, 2 * a] = np.where(va < 0, t_lo, np.inf)
                t_face[:, 2 * a + 1] = np.where(va > 0, t_hi, np.inf)
            which = np.argmin(t_face, axis=1)
            t_hit = t_face[np.arange(idx.size), which]
            hit = x0[idx] + v[idx] * t_hit[:, None]
            r = np.hypot(hit[:, 1], hit[:, 2])
            lab = np.where(which == 0, np.where(r <= radius_m, 1, 2), 3)
            labels[idx] = lab
            undecided[idx] = False
        idx = np.nonzero(undecided & ~out)[0]
        if idx.size:
            inside = in_retained(p[idx])
            labels[idx[inside]] = 0
            undecided[idx[inside]] = False
        if not undecided.any():
            break
    return labels


EXIT_NAMES = {0: "interface", 1: "source", 2: "exitPlane", 3: "farField", -1: "unclassified"}


def crossing_window(snapshots: list, dt: float, weight: float, gas, box: tuple,
                    in_retained, radius_m: float) -> dict:
    """Exact exit fluxes over consecutive per-step snapshots ``[(time, Cloud), ...]``."""
    m = gas.mass_kg
    totals = {name: {"parcels": 0, "mass": 0.0, "momentum": np.zeros(3), "energy": 0.0}
              for name in EXIT_NAMES.values()}
    new_parcels = 0
    steps = 0
    for (t0, a), (t1, b) in zip(snapshots[:-1], snapshots[1:]):
        if not np.isclose(t1 - t0, dt, rtol=1e-6):
            raise ValueError(f"snapshots {t0} -> {t1} are not one step ({dt}) apart")
        steps += 1
        gone = ~np.isin(a.ident, b.ident)
        new_parcels += int((~np.isin(b.ident, a.ident)).sum())
        if gone.any():
            lab = classify_exits(a.position[gone], a.velocity[gone], dt, box,
                                 in_retained, radius_m)
            v = a.velocity[gone]
            for code, name in EXIT_NAMES.items():
                k = lab == code
                if not k.any():
                    continue
                totals[name]["parcels"] += int(k.sum())
                totals[name]["mass"] += weight * m * k.sum()
                totals[name]["momentum"] += weight * m * v[k].sum(axis=0)
                totals[name]["energy"] += weight * m * 0.5 * (v[k] ** 2).sum()
    duration = steps * dt
    rates = {name: {"parcels_per_s": d["parcels"] / duration, "mass_kg_s": d["mass"] / duration,
                    "momentum_N": (d["momentum"] / duration).tolist(),
                    "energy_W": d["energy"] / duration, "parcels": d["parcels"]}
             for name, d in totals.items()}
    return {"steps": steps, "duration_s": duration, "exit_rates": rates,
            "new_parcels_surviving_their_first_step": new_parcels,
            "method": "exact: deletions between consecutive per-step snapshots, ray-traced"}


def lattice_lookup(cell_centres: np.ndarray, flag: np.ndarray, edges: tuple):
    """``in_region(points)`` for a structured (blockMesh box) mesh, from cell centres."""
    ex, ey, ez = (np.asarray(e) for e in edges)
    shape = (ex.size - 1, ey.size - 1, ez.size - 1)
    i = np.searchsorted(ex, cell_centres[:, 0]) - 1
    j = np.searchsorted(ey, cell_centres[:, 1]) - 1
    k = np.searchsorted(ez, cell_centres[:, 2]) - 1
    grid = np.zeros(shape, bool)
    grid[i, j, k] = flag

    def in_region(points: np.ndarray) -> np.ndarray:
        a = np.clip(np.searchsorted(ex, points[:, 0]) - 1, 0, shape[0] - 1)
        b = np.clip(np.searchsorted(ey, points[:, 1]) - 1, 0, shape[1] - 1)
        c = np.clip(np.searchsorted(ez, points[:, 2]) - 1, 0, shape[2] - 1)
        return grid[a, b, c]
    return in_region


def inserted_from_log(log_path: Path) -> np.ndarray:
    """Per-step 'Particles inserted' counts from a dsmcFoam log."""
    text = Path(log_path).read_text(errors="replace")
    return np.array([int(v) for v in re.findall(r"Particles inserted\s*=\s*(\d+)", text)])


def cell_half_range_estimate(cloud: Cloud, face_cells: np.ndarray, face_normals_out: np.ndarray,
                             face_areas: np.ndarray, cell_volumes: np.ndarray, weight: float,
                             mass: float) -> dict:
    """Snapshot estimate of the flux leaving through each face from its owner cell.

    ``sum_p W m (v.n)+ A / V`` over parcels in the cell, which is the exact
    expected crossing rate for a distribution uniform over the cell. Labelled
    an ESTIMATE; :func:`crossing_window` is the exact measurement.
    """
    out_mass = np.zeros(face_cells.size)
    order = np.argsort(cloud.cell, kind="stable")
    cells_sorted = cloud.cell[order]
    for f, (c, nrm) in enumerate(zip(face_cells, face_normals_out)):
        lo = np.searchsorted(cells_sorted, c)
        hi = np.searchsorted(cells_sorted, c, side="right")
        v = cloud.velocity[order[lo:hi]]
        vn = v @ nrm
        out_mass[f] = weight * mass * np.clip(vn, 0, None).sum() * face_areas[f] / cell_volumes[c]
    return {"mass_kg_s_per_face": out_mass, "mass_kg_s": float(out_mass.sum())}
