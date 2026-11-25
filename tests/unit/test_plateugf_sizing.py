"""Mean-free-path sizing for the plate-impingement family (spec §3)."""
from __future__ import annotations

import numpy as np
import pytest

from plumetools.plateugf import sizing

D = 0.2
S = 2.0
X_LIMIT = 7.0          # z/D = 25: plate at 5 m, 0.25 m thick, 1.75 m behind
RHO_LIMIT = 21.0


@pytest.fixture(scope="module")
def table():
    return sizing.density_table(X_LIMIT, RHO_LIMIT, 0.5 * D, S)


def regions(table, kn, x_limit=X_LIMIT):
    return sizing.refinement_regions(
        table, knudsen=kn, diameter_m=D, x_limit_m=x_limit, rho_limit_m=RHO_LIMIT)


def _inside(region, x, rho):
    x0, x1 = region["x_min_m"], region["x_max_m"]
    r = region["radius1_m"] + (region["radius2_m"] - region["radius1_m"]) * (x - x0) / (x1 - x0)
    return (x >= x0) & (x <= x1) & (rho <= r)


def _by_level(rs):
    out = {}
    for r in rs:
        out.setdefault(r["level"], []).append(r)
    return out


def test_rarefied_regions_are_geometric(table):
    assert regions(table, 100.0) == regions(table, 0.1)


@pytest.mark.parametrize("kn", [100.0, 1e-3, 1e-4])
def test_every_region_point_is_inside_its_cone(table, kn):
    rng = np.random.default_rng(1)
    x = np.exp(rng.uniform(np.log(5e-3), np.log(X_LIMIT), 20000))
    rho = np.exp(rng.uniform(np.log(1e-4), np.log(RHO_LIMIT), 20000))
    s = sizing.target_cell_size(table, x, rho, knudsen=kn, diameter_m=D)
    levels = _by_level(regions(table, kn))
    for level in levels:
        member = s < 0.64 / 2 ** (level - 1)
        covered = np.zeros(member.sum(), dtype=bool)
        for finer in (lv for lv in levels if lv >= level):
            for region in levels[finer]:
                covered |= _inside(region, x[member], rho[member])
        assert covered.all(), (level, x[member][~covered][:5], rho[member][~covered][:5])


def test_regions_nest_and_finest_holds_exit_disk(table):
    rs = regions(table, 1e-3)
    levels = _by_level(rs)
    assert sorted(levels) == [1, 2, 3, 4, 5, 6]
    reach = {lv: max(r["x_max_m"] for r in levels[lv]) for lv in levels}
    assert all(reach[lv + 1] <= reach[lv] for lv in range(1, 6))
    first = min(levels[6], key=lambda r: r["x_min_m"])
    assert first["x_min_m"] == 0.0
    assert min(first["radius1_m"], first["radius2_m"]) >= 0.5 * D + 0.01
    assert all(r["shape"] == "cone" for r in rs)
    assert len({r["name"] for r in rs}) == len(rs)
    assert all(r["name"].startswith(f"mfpL{r['level']}") for r in rs)


def test_dense_cones_hug_the_plume(table):
    # The level-6 (10 mm) region at Kn 1e-4 holds ~40 m^3; one enclosing cone
    # took ~190 m^3. The fitted cones must stay within 2x the true volume.
    x = np.linspace(5e-3, X_LIMIT, 700)
    rho = np.linspace(0.0, RHO_LIMIT, 2100)
    X, R = np.meshgrid(x, rho, indexing="ij")
    s = sizing.target_cell_size(table, X, R, knudsen=1e-4, diameter_m=D)
    true = (2 * np.pi * R * (s < 0.02)).sum() * (x[1] - x[0]) * (rho[1] - rho[0])
    fitted = sum(np.pi / 3 * (r["x_max_m"] - r["x_min_m"])
                 * (r["radius1_m"] ** 2 + r["radius1_m"] * r["radius2_m"] + r["radius2_m"] ** 2)
                 for r in regions(table, 1e-4) if r["level"] == 6)
    assert true <= fitted < 2.0 * true


def test_dense_case_reaches_further_than_rarefied(table):
    dense = _by_level(regions(table, 1e-4))
    rare = _by_level(regions(table, 100.0))
    assert max(r["x_max_m"] for r in dense[3]) > max(r["x_max_m"] for r in rare[3])


def test_cones_clip_to_domain(table):
    for r in regions(table, 1e-4, x_limit=3.0):
        assert r["x_max_m"] <= 3.0
        assert r["radius2_m"] <= RHO_LIMIT


def test_floor_and_ceiling(table):
    x = np.geomspace(5e-3, X_LIMIT, 50)
    rho = np.geomspace(1e-4, RHO_LIMIT, 50)
    X, R = np.meshgrid(x, rho)
    for kn in (100.0, 1e-4):
        s = sizing.target_cell_size(table, X, R, knudsen=kn, diameter_m=D)
        assert s.min() >= 0.01 - 1e-12
        assert s.max() <= 0.64 + 1e-12
