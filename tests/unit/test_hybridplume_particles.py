"""hybridplume: parcel files, moment estimates, exact crossing accounting.

The moment test samples a drifting Maxwellian with a fixed seed; with 200 000
parcels the standard error of n per cell is ~1/sqrt(2e5/8) = 0.6 % and of T
~0.9 %, so the 3 % tolerance is about 5 standard errors.
"""

from __future__ import annotations

import numpy as np
import pytest

from plumetools.hybridplume import particles
from plumetools.hybridplume.physics import BIRD_ARGON

HEAD = "FoamFile\n{{\n    version 2.0;\n    format ascii;\n    class {cls};\n    object {obj};\n}}\n\n"


def _write_cloud(time_dir, pos, cells, vel, ids, compact_type=True):
    d = time_dir / "lagrangian" / "dsmc"
    d.mkdir(parents=True)
    n = len(cells)
    (d / "positions").write_text(HEAD.format(cls="Cloud<passiveParticle>", obj="positions")
                                 + f"{n}\n(\n" + "\n".join(
        f"({p[0]:.12g} {p[1]:.12g} {p[2]:.12g}) {c}" for p, c in zip(pos, cells)) + "\n)\n")
    (d / "U").write_text(HEAD.format(cls="vectorField", obj="U") + f"{n}\n(\n" + "\n".join(
        f"({v[0]:.12g} {v[1]:.12g} {v[2]:.12g})" for v in vel) + "\n)\n")
    if compact_type:
        (d / "typeId").write_text(HEAD.format(cls="labelField", obj="typeId") + f"{n}{{0}}\n")
    else:
        (d / "typeId").write_text(HEAD.format(cls="labelField", obj="typeId")
                                  + f"{n}\n(\n" + "\n".join("0" for _ in range(n)) + "\n)\n")
    (d / "origProcId").write_text(HEAD.format(cls="labelField", obj="origProcId") + f"{n}{{0}}\n")
    (d / "origId").write_text(HEAD.format(cls="labelField", obj="origId") + f"{n}\n(\n"
                              + "\n".join(map(str, ids)) + "\n)\n")


def test_reads_both_list_forms(tmp_path):
    pos = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
    _write_cloud(tmp_path / "1", pos, [3, 7], [[1, 2, 3], [4, 5, 6]], [10, 11])
    _write_cloud(tmp_path / "2", pos, [3, 7], [[1, 2, 3], [4, 5, 6]], [10, 11], False)
    for t in ("1", "2"):
        c = particles.read_cloud(tmp_path / t)
        assert c.size == 2 and np.array_equal(c.cell, [3, 7])
        assert np.allclose(c.velocity[1], [4, 5, 6]) and np.array_equal(c.type_id, [0, 0])
        assert np.array_equal(c.ident & 0xFFFFFFFF, [10, 11])


def test_moments_recover_a_drifting_maxwellian():
    rng = np.random.default_rng(7)
    gas = BIRD_ARGON
    t, u = 250.0, np.array([500.0, -100.0, 30.0])
    sigma = np.sqrt(gas.gas_constant * t)
    n_p, n_cells, w = 200_000, 8, 1e12
    cells = rng.integers(0, n_cells, n_p)
    vel = u + sigma * rng.normal(size=(n_p, 3))
    cloud = particles.Cloud(np.zeros((n_p, 3)), cells, vel, np.zeros(n_p, int), np.arange(n_p))
    acc = particles.MomentAccumulator(n_cells)
    acc.add(cloud, w)
    vol = np.full(n_cells, 1e-3)
    f = acc.fields(vol, gas)
    assert np.allclose(f["n"], w * n_p / n_cells / 1e-3, rtol=0.03)
    assert np.allclose(f["T"], t, rtol=0.03)
    assert np.allclose(f["U"], u, atol=0.03 * sigma)
    assert np.allclose(f["T_dir"], t, rtol=0.05)
    assert np.allclose(f["T_parallel"], t, rtol=0.05)


def test_block_statistics_standard_error():
    blocks = [np.array([1.0, 2.0]), np.array([3.0, 2.0]), np.array([2.0, 2.0])]
    s = particles.block_statistics(blocks)
    assert np.allclose(s["mean"], [2.0, 2.0])
    assert s["stderr"][0] == pytest.approx(np.std([1, 3, 2], ddof=1) / np.sqrt(3))
    assert s["stderr"][1] == 0.0


def _cloud(pos, vel, ids):
    pos = np.asarray(pos, float)
    return particles.Cloud(pos, np.zeros(len(ids), int), np.asarray(vel, float),
                           np.zeros(len(ids), int), np.asarray(ids, np.int64))


def test_crossings_are_exact_and_classified():
    gas = BIRD_ARGON
    dt, w, m = 1e-6, 1e10, gas.mass_kg
    box = (np.array([0.0, -1.0, -1.0]), np.array([2.0, 1.0, 1.0]))
    edges = (np.linspace(0, 2, 21), np.linspace(-1, 1, 21), np.linspace(-1, 1, 21))
    centres = np.stack(np.meshgrid(*(0.5 * (e[:-1] + e[1:]) for e in edges), indexing="ij"),
                       -1).reshape(-1, 3)
    # the "CFD" region: a ball away from the source, as no DSMC parcel starts inside it
    retained = np.linalg.norm(centres - [1.0, -0.6, 0.0], axis=1) < 0.3
    in_region = particles.lattice_lookup(centres, retained, edges)
    # 0: into the retained region (interface); 1: out through x=0 inside r<=0.1
    # (source); 2: x=0 outside (exitPlane); 3: far field; 4: stays
    pos = [[1.0, -0.2, 0.0], [0.05, 0.0, 0.02], [0.05, 0.6, 0.0], [1.95, 0.0, 0.0],
           [1.0, 0.5, 0.5]]
    vel = [[0, -2e5, 0], [-1e5, 0, 0], [-1e5, 0, 0], [1e5, 0, 0], [1.0, 0, 0]]
    a = _cloud(pos, vel, [1, 2, 3, 4, 5])
    b = _cloud([[1.0, 0.5, 0.5], [0.1, 0.1, 0.1]], [[1.0, 0, 0], [5.0, 0, 0]], [5, 6])
    out = particles.crossing_window([(0.0, a), (dt, b)], dt, w, gas, box, in_region, 0.1)
    rates = out["exit_rates"]
    for name in ("interface", "source", "exitPlane", "farField"):
        assert rates[name]["parcels"] == 1
    assert rates["unclassified"]["parcels"] == 0
    assert rates["interface"]["mass_kg_s"] == pytest.approx(w * m / dt)
    assert rates["interface"]["momentum_N"][1] == pytest.approx(w * m * -2e5 / dt)
    assert rates["farField"]["energy_W"] == pytest.approx(w * m * 0.5 * 1e10 / dt)
    assert out["new_parcels_surviving_their_first_step"] == 1


def test_crossings_refuse_non_consecutive_snapshots():
    a = _cloud([[0.5, 0, 0]], [[0, 0, 0]], [1])
    with pytest.raises(ValueError):
        particles.crossing_window([(0.0, a), (3e-6, a)], 1e-6, 1.0, BIRD_ARGON,
                                  (np.zeros(3), np.ones(3)), lambda p: np.zeros(len(p), bool),
                                  0.1)


def test_log_insertion_counts(tmp_path):
    log = tmp_path / "log.dsmcFoam"
    log.write_text("Time = 1\n    Particles inserted              = 590\n"
                   "Time = 2\n    Particles inserted              = 611\n")
    assert np.array_equal(particles.inserted_from_log(log), [590, 611])
