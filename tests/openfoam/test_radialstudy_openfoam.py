"""Built native mesh geometry, later extended with a bounded solver mechanics run."""
from dataclasses import replace
from pathlib import Path
import shutil
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'unit'))
from test_radialstudy_physics import model

pytestmark = pytest.mark.needs_openfoam


def test_actual_radial_mesh_geometry(tmp_path):
    if not shutil.which('blockMesh'):
        pytest.skip('requires OpenFOAM v2512')
    from plumetools.radialstudy.mesh import initial_mesh, materialize
    from plumetools.refinement.geometry import read_geometry
    m = model(100)
    m = replace(m, r_out=.164)
    mesh = initial_mesh(m)
    # Deliberately exercise transverse subdivisions too.
    mesh = replace(mesh, transverse=3)
    destination = tmp_path/'ancestor with spaces'/'mesh'
    materialize(mesh, destination)
    actual = read_geometry(destination, time='0')
    np.testing.assert_allclose(actual.centres_m, mesh.geometry.centres_m, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(actual.volumes_m3, mesh.geometry.volumes_m3, rtol=1e-9)
    np.testing.assert_allclose(actual.max_edge_m, mesh.geometry.max_edge_m, rtol=1e-9)


def test_bounded_plus_mechanics_run(tmp_path):
    if not shutil.which('dsmcFoam+'):
        pytest.skip('requires the v2512 hyStrath plus port')
    from plumetools.radialstudy.case import write_case
    from plumetools.radialstudy.physics import source_for_kn
    from plumetools.radialstudy.mesh import initial_mesh
    from plumetools.radialstudy.run import mesh_case,run_case,collect_case
    m = model(100)
    cfg = replace(m.config,domain=replace(m.config.domain,outer_over_D=.82,angle_deg=.5),
                  mesh=replace(m.config.mesh,geometric_cap_over_D=.04))
    m = source_for_kn(cfg,100)
    mesh = initial_mesh(m)
    record = write_case(tmp_path/'mechanics only',m,mesh)
    assert mesh.n_cells <= 100
    assert record['numerics']['planned_steps'] <= 2000
    assert record['numerics']['predicted_total_parcels'] <= 5000
    case = tmp_path/'mechanics only'
    assert mesh_case(case)['n_cells'] == mesh.n_cells
    result = run_case(case)
    assert result['status'] == 'completed'
    windows = collect_case(case)
    assert len(windows) == 12
    for w in windows:
        assert all((case/w['time']/name).is_file() for name in ('rhoN_Ar','Ttra_Ar','U_Ar','dsmcNMean_Ar','mfp_Ar','mct_Ar'))
        assert (case/w['time']/'lagrangian/dsmc/positions').is_file()
    assert run_case(case,resume=True)['status'] == 'already_finished'
    from plumetools.radialstudy.post import case_plan
    from plumetools.refinement.plan import verify_sources
    plan,_ = case_plan(case)
    for relative in ('0/C','0/V','run-state.json','sampling-windows.json'):
        assert relative in plan.provenance['source_files']
        path = case/relative
        original = path.read_bytes()
        path.write_bytes(original+b'\n')
        with pytest.raises(ValueError,match='changed|stale'):
            verify_sources(case,plan)
        path.write_bytes(original)
    verify_sources(case,plan)
