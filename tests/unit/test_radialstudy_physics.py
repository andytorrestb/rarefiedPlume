from dataclasses import replace

import numpy as np
import pytest

from test_radialstudy_config import STUDY


def model(kn):
    from plumetools.radialstudy.config import load_study
    from plumetools.radialstudy.physics import source_for_kn
    cfg = load_study(STUDY)
    # Existing mechanics tests exercise the original fixed-duration configuration.
    cfg = replace(cfg,domain=replace(cfg.domain,outer_over_D=25,angle_deg=None),
                  mesh=replace(cfg.mesh,initial_cells=None),controller=replace(cfg.controller,enabled=False))
    return source_for_kn(cfg, kn)


def test_mapping_and_density_scaling_reuse_existing_source():
    from plumetools.cai2012ext.interface import equivalent_source
    a, b = model(.001), model(.01)
    assert a.source.r_a / a.source.diameter_m == pytest.approx(.7172871)
    assert a.source.T_a == pytest.approx(299.845, abs=.001)
    assert a.state(a.r_in)['n'] / b.state(b.r_in)['n'] == pytest.approx(10)
    expected = equivalent_source(a.reference)
    for key, value in a.state(a.r_in).items():
        np.testing.assert_allclose(value, expected.state(a.r_in, 0)[key])


def test_breakdown_outcomes_and_definitions():
    from plumetools.radialstudy.physics import predict_breakdown
    assert predict_breakdown(model(.0001))['r_over_D'] == pytest.approx(14.992, abs=.002)
    assert predict_breakdown(model(.001))['r_over_D'] == pytest.approx(2.788, abs=.002)
    assert predict_breakdown(model(.01))['status'] == 'ALREADY_NONCONTINUUM_AT_INLET'
    m = model(.0001)
    assert predict_breakdown(replace(m, r_out=2.))['status'] == 'NOT_REACHED_IN_DOMAIN'
    assert predict_breakdown(m, definition='radial_approx')['r_over_D'] > predict_breakdown(m)['r_over_D']
    assert predict_breakdown(replace(m, r_in=float('nan')))['status'] == 'INVALID_STATE'


def test_inlet_override_must_keep_upstream_continuum():
    from plumetools.radialstudy.physics import source_for_kn
    m = model(.001)
    with pytest.raises(ValueError, match='upstream'):
        source_for_kn(replace(m.config, domain=replace(m.config.domain, inlet_over_D=3.)), .001)


def test_revision_cannot_change_domain_or_inlet():
    from plumetools.radialstudy.physics import require_same_physics
    m = model(.001)
    require_same_physics(m.identity(angle_rad=.001), m.identity(angle_rad=.001))
    changed = {**m.identity(angle_rad=.001), 'r_in_m': m.r_in * 1.1}
    with pytest.raises(ValueError, match='physical'):
        require_same_physics(m.identity(angle_rad=.001), changed)
