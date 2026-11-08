from __future__ import annotations

import importlib
import math
import numpy as np
import pytest
from _refinement import write_cartesian_case, write_field


def test_rotated_stretched_cell_uses_actual_edges(tmp_path):
    from plumetools.mesh.polymesh import read_points_array
    from _helpers import HEADER
    case = write_cartesian_case(tmp_path / "case")
    points = read_points_array(case)
    points[:, 0] *= 2
    rotation = np.array([[1, -1, 0], [1, 1, 0], [0, 0, math.sqrt(2)]]) / math.sqrt(2)
    points = points @ rotation.T
    path = case / "constant/polyMesh/points"
    path.write_text(HEADER.format(cls="vectorField", obj="points") + f"\n{len(points)}\n(\n"
                    + "\n".join("(" + " ".join(format(v, ".17g") for v in row) + ")" for row in points) + "\n)\n")
    write_field(case / "1/C", np.array([[1, .5, .5], [3, .5, .5]]) @ rotation.T, [0, 1, 0, 0, 0, 0, 0])
    write_field(case / "1/V", [2, 2], [0, 3, 0, 0, 0, 0, 0])
    geometry = importlib.import_module("plumetools.refinement.geometry").read_geometry(case, time="1")
    assert geometry.max_edge_m == pytest.approx([2, 2])
    assert geometry.volumes_m3 == pytest.approx([2, 2])


def test_mesh_identity_changes_with_points_not_absolute_directory(tmp_path):
    import shutil
    reader = importlib.import_module("plumetools.refinement.geometry").read_geometry
    case = write_cartesian_case(tmp_path / "first")
    copied = tmp_path / "second"; shutil.copytree(case, copied)
    a = reader(case, time="1").mesh_sha256
    assert reader(copied, time="1").mesh_sha256 == a
    path = copied / "constant/polyMesh/points"
    path.write_text(path.read_text().replace("(0 0 0)", "(-0.1 0 0)"))
    assert reader(copied, time="1").mesh_sha256 != a
