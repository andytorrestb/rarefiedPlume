"""DSMC quality checks: which conditions warn, and which are hard errors.

The split is the point of the module, so most tests are about *severity* rather
than about the numbers.
"""

from __future__ import annotations

import json
import math
from dataclasses import replace

import numpy as np
import pytest

from plumetools.config import StagnationConfig
from plumetools.markelov1999 import checks as mc
from plumetools.markelov1999 import resolution as res
from plumetools.markelov1999.constants import PSI_TO_PA
from plumetools.markelov1999.geometry import MarkelovGeometry
from plumetools.mesh.boundary import PatchInfo

from test_markelov_mesh import make_cfg  # noqa: E402

CELL_SIZES = {"inflow": 0.00625, "cylinder": 0.00625, "plate": 0.00625}


@pytest.fixture
def cfg():
    base = make_cfg()
    return replace(base, stagnation=StagnationConfig(
        p0_pa=5 * PSI_TO_PA, T0_K=300.0, throat_radius_m=0.00041275))


@pytest.fixture
def geom():
    return MarkelovGeometry()


@pytest.fixture
def patches():
    return {
        "symmetry": PatchInfo("symmetry", "symmetry", 2102, 0),
        "upstreamVacuum": PatchInfo("upstreamVacuum", "patch", 530, 2102),
        "vacuum": PatchInfo("vacuum", "patch", 2208, 2632),
        "inflow": PatchInfo("inflow", "patch", 1596, 4840),
        "cylinder": PatchInfo("cylinder", "wall", 3054, 6436),
        "plate": PatchInfo("plate", "wall", 1608, 9490),
    }


class FakeInflow:
    """The subset of InflowResult the checks read."""

    def __init__(self, n=1596, rhoN=None, T=None, theta=None):
        self.n_faces = n
        self.rhoN = np.full(n, 1.0e18) if rhoN is None else np.asarray(rhoN)
        self.T = np.full(n, 300.0) if T is None else np.asarray(T)
        self.theta = np.linspace(0.03, 1.55, n) if theta is None else np.asarray(theta)
        self.U = np.tile([789.0, 0.0, 0.0], (n, 1))
        self.areas = np.full(n, 4.6e-5)


def status(report, name):
    return next(r.status for r in report.results if r.name == name)


# --------------------------------------------------------------------------- #
# a healthy case
# --------------------------------------------------------------------------- #

def test_the_baseline_case_passes_everything_except_occupancy(cfg, geom, patches):
    estimate = res.estimate(cfg, geom, CELL_SIZES)
    report = mc.run_checks(cfg, geom, inflow=FakeInflow(), patches=patches,
                           estimate=estimate, cell_sizes=CELL_SIZES)

    assert report.failures == []
    report.raise_if_failed()   # must not raise
    # The plate is under target at the shipped weighting; that is a warning.
    assert status(report, "particles_per_cell") == "warn"


def test_the_report_is_machine_readable(cfg, geom, patches):
    estimate = res.estimate(cfg, geom, CELL_SIZES)
    report = mc.run_checks(cfg, geom, inflow=FakeInflow(), patches=patches,
                           estimate=estimate, cell_sizes=CELL_SIZES)
    document = report.as_dict()

    json.dumps(document)   # must not raise
    assert document["passed"] is True
    assert document["n_failures"] == 0
    assert all({"name", "status", "message"} <= set(r) for r in document["results"])


# --------------------------------------------------------------------------- #
# hard errors
# --------------------------------------------------------------------------- #

def test_a_nonpositive_pressure_or_temperature_fails(cfg, geom):
    for bad in (replace(cfg, stagnation=replace(cfg.stagnation, p0_pa=0.0)),
                replace(cfg, stagnation=replace(cfg.stagnation, T0_K=-1.0))):
        report = mc.run_checks(bad, geom)
        assert status(report, "positive_inputs") == "fail"
        with pytest.raises(mc.DsmcCheckError):
            report.raise_if_failed()


def test_a_nan_input_fails(cfg, geom):
    bad = replace(cfg, stagnation=replace(cfg.stagnation, p0_pa=float("nan")))
    assert status(mc.run_checks(bad, geom), "positive_inputs") == "fail"


def test_an_inflow_surface_reaching_the_cylinder_fails(cfg):
    """Requirement 5 and 10 both single this out."""
    intersecting = MarkelovGeometry(inflow_radius_m=0.25)
    report = mc.run_checks(cfg, intersecting)
    assert status(report, "inflow_cylinder_clearance") == "fail"


def test_the_clearance_is_reported_when_it_passes(cfg, geom):
    """Requirement 5 asks for the minimum distance, not just a verdict."""
    report = mc.run_checks(cfg, geom)
    message = next(r.message for r in report.results
                   if r.name == "inflow_cylinder_clearance")
    assert "0.069850" in message
    assert "2.75 in" in message


def test_a_missing_patch_fails(cfg, geom, patches):
    del patches["plate"]
    report = mc.run_checks(cfg, geom, patches=patches)
    assert status(report, "required_patches") == "fail"


def test_an_empty_physical_patch_fails(cfg, geom, patches):
    """snappyHexMesh creates a patch even when it carved nothing, so an empty
    body patch means the body is simply not in the mesh."""
    patches["plate"] = PatchInfo("plate", "wall", 0, 9490)
    report = mc.run_checks(cfg, geom, patches=patches)
    assert status(report, "empty_patches") == "fail"
    assert "not in the mesh at all" in next(
        r.message for r in report.results if r.name == "empty_patches")


def test_a_time_step_crossing_more_than_one_cell_fails(cfg, geom):
    """Requirement 10: this is one of the conditions that must fail, not warn."""
    coarse = replace(cfg, dsmc=replace(cfg.dsmc, delta_t_s=1.0e-5))
    report = mc.run_checks(coarse, geom, cell_sizes=CELL_SIZES)
    assert status(report, "time_step_courant") == "fail"
    message = next(r.message for r in report.results if r.name == "time_step_courant")
    assert "Reduce dsmc.delta_t_s" in message


def test_a_marginal_time_step_only_warns(cfg, geom):
    marginal = replace(cfg, dsmc=replace(cfg.dsmc, delta_t_s=3.0e-6))
    assert status(mc.run_checks(marginal, geom, cell_sizes=CELL_SIZES),
                  "time_step_courant") == "warn"


def test_the_shipped_time_step_passes(cfg, geom):
    assert status(mc.run_checks(cfg, geom, cell_sizes=CELL_SIZES),
                  "time_step_courant") == "pass"


def test_non_finite_inflow_values_fail(cfg, geom):
    rhoN = np.full(10, 1.0e18)
    rhoN[3] = float("nan")
    report = mc.run_checks(cfg, geom, inflow=FakeInflow(n=10, rhoN=rhoN))
    assert status(report, "inflow_finite") == "fail"


def test_a_nonpositive_inflow_density_fails(cfg, geom):
    rhoN = np.full(10, 1.0e18)
    rhoN[7] = 0.0
    report = mc.run_checks(cfg, geom, inflow=FakeInflow(n=10, rhoN=rhoN))
    assert status(report, "inflow_density_positive") == "fail"


def test_a_zero_inflow_temperature_fails(cfg, geom):
    T = np.full(10, 300.0)
    T[2] = 0.0
    report = mc.run_checks(cfg, geom, inflow=FakeInflow(n=10, T=T))
    assert status(report, "inflow_temperature_positive") == "fail"
    assert "Zero boundary temperature detected" in next(
        r.message for r in report.results if r.name == "inflow_temperature_positive")


def test_an_inflow_face_outside_the_plume_cone_fails(cfg, geom):
    """Requirement 10 lists it. On a hemispherical source with a diatomic gas it
    cannot happen, so if it does the mesh is not the described geometry."""
    theta = np.linspace(0.03, 1.55, 10)
    theta[5] = 2.5   # beyond the 2.277 rad limiting angle
    report = mc.run_checks(cfg, geom, inflow=FakeInflow(n=10, theta=theta))
    assert status(report, "inflow_inside_plume_cone") == "fail"


def test_the_hemisphere_is_entirely_inside_the_cone_for_nitrogen(cfg, geom):
    report = mc.run_checks(cfg, geom, inflow=FakeInflow())
    assert status(report, "inflow_inside_plume_cone") == "pass"
    assert "130.45" in next(r.message for r in report.results
                            if r.name == "inflow_inside_plume_cone")


def test_every_failure_is_listed_not_just_the_first(cfg, geom, patches):
    """They usually share a cause, and fixing them one run at a time is slow."""
    broken = replace(cfg,
                     stagnation=replace(cfg.stagnation, p0_pa=0.0),
                     dsmc=replace(cfg.dsmc, delta_t_s=1.0e-4))
    report = mc.run_checks(broken, geom, cell_sizes=CELL_SIZES)
    assert len(report.failures) >= 2
    with pytest.raises(mc.DsmcCheckError) as excinfo:
        report.raise_if_failed()
    assert "positive_inputs" in str(excinfo.value)
    assert "time_step_courant" in str(excinfo.value)


# --------------------------------------------------------------------------- #
# warnings
# --------------------------------------------------------------------------- #

def test_a_cell_coarser_than_the_mean_free_path_warns_not_fails(cfg, geom):
    """It biases the collision rate; it does not make the case meaningless."""
    dense = replace(cfg, stagnation=replace(cfg.stagnation, p0_pa=475 * PSI_TO_PA))
    coarse = {k: 0.05 for k in CELL_SIZES}
    estimate = res.estimate(dense, geom, coarse)
    report = mc.run_checks(dense, geom, estimate=estimate, cell_sizes=coarse)

    assert status(report, "cell_over_mean_free_path") == "warn"
    assert report.failures == []


def test_low_occupancy_warns_and_names_the_trade(cfg, geom):
    estimate = res.estimate(cfg, geom, CELL_SIZES)
    report = mc.run_checks(cfg, geom, estimate=estimate, cell_sizes=CELL_SIZES)
    message = next(r.message for r in report.results if r.name == "particles_per_cell")
    assert "one weight for the whole domain" in message
    assert "sizing_region" in message


def test_occupancy_is_never_a_hard_error(cfg, geom):
    """It is an estimate; only the post-run audit is evidence."""
    tiny = {k: 0.0005 for k in CELL_SIZES}
    estimate = res.estimate(cfg, geom, tiny, n_equivalent_particles=1.0e15)
    report = mc.run_checks(cfg, geom, estimate=estimate, cell_sizes=tiny)
    assert status(report, "particles_per_cell") == "warn"
    assert report.failures == []


# --------------------------------------------------------------------------- #
# disabling
# --------------------------------------------------------------------------- #

def test_disabling_the_checks_is_itself_reported(cfg, geom):
    """Silently running none of them would be worse than running none loudly."""
    off = replace(cfg, checks=replace(cfg.checks, enabled=False))
    report = mc.run_checks(off, geom)
    assert status(report, "checks_enabled") == "warn"
    assert len(report.results) == 1


# --------------------------------------------------------------------------- #
# the summary file
# --------------------------------------------------------------------------- #

def test_the_summary_is_written_as_json(tmp_path):
    path = mc.write_summary(tmp_path, {"case": {"name": "p005psi"}, "checks": {}})
    assert path.name == "case-summary.json"
    assert json.loads(path.read_text(encoding="utf-8"))["case"]["name"] == "p005psi"


def test_the_summary_is_stable_across_writes(tmp_path):
    """Sorted keys, so a re-run produces the same bytes for the same data."""
    summary = {"b": 2, "a": 1}
    first = mc.write_summary(tmp_path, summary).read_bytes()
    second = mc.write_summary(tmp_path, summary).read_bytes()
    assert first == second
