import csv
import sys

import pytest

from test_radialstudy_run import prepared


def test_inventory_and_csv_work_without_matplotlib(tmp_path,monkeypatch):
    from plumetools.radialstudy.post import analyse_study
    from plumetools.radialstudy.report import write_study_report
    p,r = prepared(tmp_path)
    analysis = analyse_study(tmp_path)
    assert len(analysis['cases']) == 1
    assert analysis['cases'][0]['status'] == 'missing_measurements'
    assert analysis['comparisons'] == []
    monkeypatch.setitem(sys.modules,'matplotlib',None)
    files = write_study_report(tmp_path/'report',analysis,plots=False)
    assert any(p.suffix == '.csv' for p in files)
    rows = list(csv.DictReader((tmp_path/'report/summary.csv').open()))
    assert rows[0]['status'] == 'missing_measurements'
    assert 'missing_measurements' in (tmp_path/'report/report.md').read_text()


def test_optional_metric_requires_complete_independent_evidence():
    from plumetools.radialstudy.post import comparison_metrics
    import numpy as np
    profile = {'boundaries_m':[1.,2.,3.], 'rho':[1.,1.], 'momentum':[[2.,0,0],[2.,0,0]],
               'energy':[500000.,500000.]}
    from plumetools.cai2012.gas import ARGON
    out = comparison_metrics(profile,profile,ARGON,targets={'n':.1},trusted=False)
    assert out['n']['verdict'] == 'inconclusive'
    assert out['n']['E'] == pytest.approx(0)
