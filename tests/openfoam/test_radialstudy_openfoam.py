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
    materialize(mesh, tmp_path/'mesh')
    actual = read_geometry(tmp_path/'mesh', time='0')
    np.testing.assert_allclose(actual.centres_m, mesh.geometry.centres_m, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(actual.volumes_m3, mesh.geometry.volumes_m3, rtol=1e-9)
    np.testing.assert_allclose(actual.max_edge_m, mesh.geometry.max_edge_m, rtol=1e-9)
