"""Flow slices on the symmetry planes (phase B, B4)."""
from __future__ import annotations

import numpy as np

from plumetools.plateugf import viz


def test_slice_cells_touch_the_plane():
    centres = np.array([[1.0, 0.005, 0.3],     # 10 mm cell touching y = 0
                        [1.0, 0.015, 0.3],     # 10 mm cell one layer off
                        [2.0, 0.32, 0.9],      # 0.64 m cell touching y = 0
                        [2.0, 0.96, 0.9]])     # 0.64 m cell one layer off
    volumes = np.array([1e-6, 1e-6, 0.64 ** 3, 0.64 ** 3])
    assert viz.touching(centres, volumes, axis=1).tolist() == [True, False, True, False]


def test_slices_write_one_image_per_field_and_plane(tmp_path):
    n = 60
    rng = np.random.default_rng(0)
    centres = np.column_stack([rng.uniform(0.01, 7, n), rng.uniform(0, 0.004, n),
                               rng.uniform(0.01, 8, n)])
    centres[n // 2:, [1, 2]] = centres[n // 2:, [2, 1]]
    data = {"C": centres, "V": np.full(n, 1e-6), "n_over_n0": rng.uniform(1e-6, 1, n),
            "U": rng.uniform(0, 900, n), "T": rng.uniform(50, 400, n)}
    written = viz.draw(data, tmp_path, label="synthetic", plate=(5.0, 5.25, 17.07, 5.94))
    names = sorted(p.name for p in written)
    assert names == sorted(f"slice-{f}-{p}.png" for f in ("n", "U", "T") for p in ("y0", "z0"))
