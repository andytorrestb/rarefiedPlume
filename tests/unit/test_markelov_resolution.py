"""Particle-resolution estimation and the post-run audit.

The distinction the tests enforce: the estimate *chooses* a weight, the audit
*measures* what happened. Only the audit can say whether the target was met.
"""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
import pytest

from plumetools.markelov1999 import resolution as res
from plumetools.markelov1999.constants import PSI_TO_PA
from plumetools.markelov1999.geometry import MarkelovGeometry

from test_markelov_mesh import make_cfg  # noqa: E402  (shared fixture builder)

CELL_SIZES = {"inflow": 0.00625, "cylinder": 0.00625, "plate": 0.00625}


@pytest.fixture
def cfg():
    from dataclasses import replace as dc_replace
    from plumetools.config import StagnationConfig
    base = make_cfg()
    return dc_replace(base, stagnation=StagnationConfig(
        p0_pa=5 * PSI_TO_PA, T0_K=300.0, throat_radius_m=0.00041275))


@pytest.fixture
def geom():
    return MarkelovGeometry()


# --------------------------------------------------------------------------- #
# mean free path
# --------------------------------------------------------------------------- #

def test_mean_free_path_matches_the_hard_sphere_formula():
    """lambda = 1/(sqrt(2) pi d^2 n)."""
    n, d = 2.9674e18, 4.17e-10
    assert res.mean_free_path(n, d) == pytest.approx(
        1.0 / (math.sqrt(2.0) * math.pi * d * d * n), rel=1e-15)


def test_mean_free_path_is_infinite_in_a_vacuum():
    assert res.mean_free_path(0.0, 4.17e-10) == float("inf")


def test_mean_free_path_scales_inversely_with_density():
    a = res.mean_free_path(1e18, 4.17e-10)
    b = res.mean_free_path(1e19, 4.17e-10)
    assert a / b == pytest.approx(10.0, rel=1e-12)


# --------------------------------------------------------------------------- #
# regions
# --------------------------------------------------------------------------- #

def test_the_three_regions_of_interest(cfg, geom):
    regions = res.regions_of_interest(cfg, geom, CELL_SIZES)
    assert [r.name for r in regions] == ["cylinder", "wake", "plate"]


def test_regions_sit_where_the_requirement_says(cfg, geom):
    regions = {r.name: r for r in res.regions_of_interest(cfg, geom, CELL_SIZES)}

    assert regions["cylinder"].point[0] == pytest.approx(geom.cylinder_upstream_x_m)
    assert regions["wake"].point[0] == pytest.approx(
        geom.cylinder_downstream_x_m + cfg.resolution.wake_offset_m)
    assert regions["plate"].point[0] == pytest.approx(geom.plate_upstream_x_m)
    for region in regions.values():
        assert region.point[1] == 0.0, "regions sit on the plume axis"


def test_the_wake_region_uses_the_cylinder_cell_size(cfg, geom):
    """It is inside the cylinder surface's refinement halo, not the plate's."""
    regions = {r.name: r for r in res.regions_of_interest(cfg, geom, CELL_SIZES)}
    assert regions["wake"].cell_size_m == CELL_SIZES["cylinder"]


def test_region_geometry_helpers(cfg, geom):
    region = res.regions_of_interest(cfg, geom, CELL_SIZES)[0]
    assert region.cell_volume_m3 == pytest.approx(0.00625 ** 3)
    assert region.radius_m == pytest.approx(geom.cylinder_upstream_x_m)
    assert region.theta_rad == pytest.approx(0.0, abs=1e-12)


# --------------------------------------------------------------------------- #
# the estimate
# --------------------------------------------------------------------------- #

def test_the_sizing_region_hits_the_target_exactly(cfg, geom):
    """That is what "sized on" means: the weight is chosen so it does."""
    estimate = res.estimate(cfg, geom, CELL_SIZES)
    assert estimate.by_name("cylinder").particles_per_cell == pytest.approx(20.0, rel=1e-9)


def test_the_other_regions_follow_from_the_flow_and_the_mesh(cfg, geom):
    """With one weight for the whole domain, the occupancy ratio between two
    regions is n*V_cell there over n*V_cell here -- nothing configurable."""
    estimate = res.estimate(cfg, geom, CELL_SIZES)
    cylinder = estimate.by_name("cylinder")
    plate = estimate.by_name("plate")

    expected = ((plate.number_density_per_m3 * plate.region.cell_volume_m3)
                / (cylinder.number_density_per_m3 * cylinder.region.cell_volume_m3))
    assert plate.particles_per_cell / cylinder.particles_per_cell == pytest.approx(
        expected, rel=1e-12)


def test_sizing_on_the_plate_meets_the_target_there_instead(cfg, geom):
    on_plate = replace(cfg, resolution=replace(cfg.resolution, sizing_region="plate"))
    estimate = res.estimate(on_plate, geom, CELL_SIZES)

    assert estimate.by_name("plate").particles_per_cell == pytest.approx(20.0, rel=1e-9)
    assert estimate.by_name("cylinder").particles_per_cell > 20.0


def test_sizing_on_the_sparsest_region_costs_more_particles(cfg, geom):
    """The trade the report describes, measured rather than asserted."""
    on_cylinder = res.estimate(cfg, geom, CELL_SIZES)
    on_plate = res.estimate(
        replace(cfg, resolution=replace(cfg.resolution, sizing_region="plate")),
        geom, CELL_SIZES)

    assert on_plate.total_particles > on_cylinder.total_particles
    assert on_plate.n_equivalent_particles < on_cylinder.n_equivalent_particles


def test_the_weight_scales_with_reservoir_pressure(cfg, geom):
    """Density is linear in p0, so at a fixed occupancy target the weight is too
    -- and the total particle count is therefore pressure-independent."""
    from plumetools.config import StagnationConfig
    weights, totals = [], []
    for psi in (5, 25, 100, 475):
        case = replace(cfg, stagnation=StagnationConfig(
            p0_pa=psi * PSI_TO_PA, T0_K=300.0, throat_radius_m=0.00041275))
        estimate = res.estimate(case, geom, CELL_SIZES)
        weights.append(estimate.n_equivalent_particles)
        totals.append(estimate.total_particles)

    for i, psi in enumerate((5, 25, 100, 475)):
        assert weights[i] / weights[0] == pytest.approx(psi / 5.0, rel=1e-12)
    for total in totals:
        assert total == pytest.approx(totals[0], rel=1e-12)


def test_an_explicit_weight_is_honoured(cfg, geom):
    estimate = res.estimate(cfg, geom, CELL_SIZES, n_equivalent_particles=1.0e11)
    assert estimate.n_equivalent_particles == 1.0e11
    assert estimate.by_name("cylinder").particles_per_cell != pytest.approx(20.0)


def test_halving_the_cell_divides_occupancy_by_eight(cfg, geom):
    """A volume ratio, which is why the plate refinement level matters so much
    to the particle budget."""
    coarse = res.estimate(cfg, geom, CELL_SIZES, n_equivalent_particles=1.0e11)
    fine_sizes = dict(CELL_SIZES, plate=CELL_SIZES["plate"] / 2)
    fine = res.estimate(cfg, geom, fine_sizes, n_equivalent_particles=1.0e11)

    assert (coarse.by_name("plate").particles_per_cell
            / fine.by_name("plate").particles_per_cell) == pytest.approx(8.0, rel=1e-9)


def test_the_report_states_the_uniform_weight_limitation(cfg, geom):
    """Requirement 9: if variable weighting is unsupported, say so."""
    text = "\n".join(res.estimate(cfg, geom, CELL_SIZES).report())
    assert "ONE particle weight for the whole domain" in text
    assert "No adaptive weighting" in text
    assert "ESTIMATE" in text
    assert "not a measurement" in text


def test_the_report_names_the_regions_that_fall_short(cfg, geom):
    text = "\n".join(res.estimate(cfg, geom, CELL_SIZES).report())
    assert "'wake'" in text or "wake" in text
    assert "fall below" in text


def test_the_report_lists_every_region_whichever_is_sized_on(cfg, geom):
    for region in ("cylinder", "plate", "wake"):
        sized = replace(cfg, resolution=replace(cfg.resolution, sizing_region=region))
        text = "\n".join(res.estimate(sized, geom, CELL_SIZES).report())
        for name in ("cylinder", "wake", "plate"):
            assert name in text


def test_cell_over_mean_free_path_is_reported_per_region(cfg, geom):
    estimate = res.estimate(cfg, geom, CELL_SIZES)
    for region in estimate.regions:
        assert region.cell_over_mfp == pytest.approx(
            region.region.cell_size_m / region.mean_free_path_m, rel=1e-12)


def test_density_is_zero_outside_the_plume_cone(cfg):
    """The condition the sizing guard exists for. All three regions sit on the
    axis, so none can currently trigger it -- but the density function it calls
    genuinely returns zero out there, so the guard is not guessing."""
    on_axis = res.number_density_at((0.2, 0.0, 0.0), cfg)
    behind = res.number_density_at((-0.2, 0.0, 0.0), cfg)
    assert on_axis > 0.0
    assert behind == 0.0


def test_a_region_with_no_density_cannot_size_the_weight(cfg, geom, monkeypatch):
    """Better an error naming the region than an infinite weight reaching
    dsmcProperties, where it would silently produce an empty simulation."""
    monkeypatch.setattr(res, "number_density_at", lambda point, config: 0.0)
    with pytest.raises(ValueError, match="cannot be derived"):
        res.estimate(cfg, geom, CELL_SIZES)


def test_total_particle_count_is_the_molecule_count_over_the_weight(cfg, geom):
    estimate = res.estimate(cfg, geom, CELL_SIZES)
    assert estimate.total_particles == pytest.approx(
        estimate.total_molecules / estimate.n_equivalent_particles, rel=1e-12)


# --------------------------------------------------------------------------- #
# the audit
# --------------------------------------------------------------------------- #

def test_audit_measures_occupancy_near_each_region(cfg, geom):
    """dsmcRhoN holds the literal parcel count per cell: DSMCCloud::calculateFields
    increments it once per parcel and never scales it."""
    points = np.array([
        [geom.cylinder_upstream_x_m, 0.0, 0.0],
        [geom.cylinder_upstream_x_m + 0.005, 0.0, 0.0],
        [geom.plate_upstream_x_m, 0.0, 0.0],
        [geom.cylinder_downstream_x_m + cfg.resolution.wake_offset_m, 0.0, 0.0],
        [0.8, 0.25, 0.3],   # far field, outside every region
    ])
    counts = np.array([20.0, 22.0, 4.0, 7.0, 0.0])

    result = res.audit(counts, points, cfg, geom, time="0.004", field="dsmcRhoN")

    assert result.region_occupancy["cylinder"] == pytest.approx(21.0)
    assert result.region_occupancy["plate"] == pytest.approx(4.0)
    assert result.region_occupancy["wake"] == pytest.approx(7.0)
    assert result.n_cells == 5
    assert result.total_particles == pytest.approx(53.0)
    assert result.occupied_cells == 4


def test_audit_reports_the_target_as_not_met_when_any_region_falls_short(cfg, geom):
    points = np.array([
        [geom.cylinder_upstream_x_m, 0.0, 0.0],
        [geom.plate_upstream_x_m, 0.0, 0.0],
        [geom.cylinder_downstream_x_m + cfg.resolution.wake_offset_m, 0.0, 0.0],
    ])
    short = res.audit(np.array([25.0, 3.0, 30.0]), points, cfg, geom,
                      time="t", field="dsmcRhoN")
    assert not short.meets_target()
    assert "NOT met" in "\n".join(short.report())

    met = res.audit(np.array([25.0, 21.0, 30.0]), points, cfg, geom,
                    time="t", field="dsmcRhoN")
    assert met.meets_target()
    assert "MET" in "\n".join(met.report())


def test_audit_with_no_region_data_does_not_claim_the_target_is_met(cfg, geom):
    """An empty region set must not read as success."""
    far = np.array([[0.85, 0.28, 0.33]])
    result = res.audit(np.array([0.0]), far, cfg, geom, time="t", field="dsmcRhoN")
    assert result.region_occupancy == {}
    assert not result.meets_target()


def test_audit_rejects_mismatched_arrays(cfg, geom):
    with pytest.raises(ValueError, match="cell values but"):
        res.audit(np.array([1.0, 2.0]), np.zeros((3, 3)), cfg, geom,
                  time="t", field="dsmcRhoN")


def test_audit_report_marks_low_regions(cfg, geom):
    points = np.array([[geom.cylinder_upstream_x_m, 0.0, 0.0],
                       [geom.plate_upstream_x_m, 0.0, 0.0]])
    text = "\n".join(res.audit(np.array([25.0, 1.0]), points, cfg, geom,
                               time="t", field="dsmcRhoN").report())
    assert "LOW plate" in text
    assert "OK  cylinder" in text
    assert "measured" in text
