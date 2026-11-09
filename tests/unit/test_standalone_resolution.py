"""The standalone extraction must reproduce the production physics."""
from dataclasses import replace
from pathlib import Path
import csv
import os
import runpy
import shutil
import subprocess
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT/'cases/cai2012-radial-resolution/standalone_resolution.py'


def standalone():
    assert SCRIPT.is_file(), 'Missing standalone resolution script'
    return runpy.run_path(str(SCRIPT))


@pytest.mark.parametrize('definition', ['general', 'radial_approx'])
def test_extraction_matches_existing_states_and_breakdown(definition):
    from plumetools.radialstudy.config import load_study
    from plumetools.radialstudy.physics import predict_breakdown, source_for_kn
    script = standalone()
    cfg = load_study(ROOT/'cases/cai2012-radial-resolution/study.yaml')
    cfg = replace(cfg, domain=replace(cfg.domain, outer_over_D=25))
    for kn in script['Kn_e']:
        model = source_for_kn(cfg, float(kn))
        source = script['make_source'](kn)
        radii = np.geomspace(model.r_in, 25*model.source.diameter_m, 12)
        actual = script['state'](radii, source)
        expected = model.state(radii)
        for key in ('M', 'n', 'T', 'U', 'lambda', 'c_bar'):
            np.testing.assert_allclose(actual[key], expected[key], rtol=1e-10)
        np.testing.assert_allclose(script['bird'](radii, source, definition),
                                   model.source.bird(radii, 0, definition), rtol=1e-10)
        expected_b = predict_breakdown(model, definition=definition)
        actual_b = script['breakdown'](kn, definition=definition)
        assert actual_b['status'] == expected_b['status']
        if 'r_over_D' in expected_b:
            assert actual_b['r_over_D'] == pytest.approx(expected_b['r_over_D'], rel=1e-9)
            expected_lambda = float(model.state(expected_b['r_m'])['lambda'])
            assert actual_b['lambda_breakdown_m'] == pytest.approx(expected_lambda, rel=1e-9)
            assert actual_b['Kn_breakdown_D_e'] == pytest.approx(expected_lambda/model.source.diameter_m, rel=1e-9)
        else:
            assert actual_b['r_over_D'] is None
            assert actual_b['lambda_breakdown_m'] is None
            assert actual_b['Kn_breakdown_D_e'] is None


def test_missing_crossing_and_invalid_inputs_are_not_fake_radii():
    s = standalone()
    result = s['breakdown'](.0001, max_radius_over_D=10)
    assert result['status'] == 'NOT_REACHED_IN_DOMAIN'
    assert result['r_over_D'] is None
    assert result['lambda_breakdown_m'] is None
    assert result['Kn_breakdown_D_e'] is None
    for kn in (0, -1, float('nan')):
        with pytest.raises(ValueError):
            s['breakdown'](kn)


def test_copied_script_runs_outside_repository(tmp_path):
    standalone()
    script = tmp_path/'standalone_resolution.py'
    shutil.copyfile(SCRIPT, script)
    environment = {**os.environ, 'PYTHONPATH': '', 'MPLBACKEND': 'Agg',
                   'MPLCONFIGDIR': str(tmp_path/'matplotlib')}
    result = subprocess.run([sys.executable, '-B', str(script)], cwd=tmp_path,
                            env=environment, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    output = tmp_path/'standalone-resolution-results'
    assert (output/'cell_sizes.png').stat().st_size > 1000
    assert (output/'breakdown_radius.png').stat().st_size > 1000
    assert len((output/'resolution.csv').read_text().splitlines()) == 8
    with (output/'resolution.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        if row['status'] == 'CROSSING_IN_DOMAIN':
            assert row['sizing_basis'] == 'breakdown'
            local_lambda = float(row['lambda_breakdown_m'])
        elif row['status'] == 'ALREADY_NONCONTINUUM_AT_INLET':
            assert row['sizing_basis'] == 'exit'
            local_lambda = float(row['lambda_e_m'])
        else:
            assert row['dx_min_m'] == row['dx_max_m'] == ''
            continue
        assert float(row['lambda_inlet_m']) == pytest.approx(local_lambda)
        assert float(row['dx_min_m']) == pytest.approx(.33*local_lambda)
        assert float(row['dx_max_m']) == pytest.approx(local_lambda)
        assert float(row['dx_exit_min_m']) == pytest.approx(.33*float(row['lambda_e_m']))
        assert float(row['dx_exit_max_m']) == pytest.approx(float(row['lambda_e_m']))
