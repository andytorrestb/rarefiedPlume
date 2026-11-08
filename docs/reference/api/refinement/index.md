# `plumetools.refinement`

Case-independent measured DSMC sizing and immutable target plans. See the
[workflow](../../../mesh-refinement.md). NumPy/PyYAML and the existing shared
VHS/OpenFOAM readers are reused; native construction requires OpenFOAM v2512.

::: plumetools.refinement
    options:
      members: false

| Module | Purpose | Principal interfaces |
|---|---|---|
| `plumetools.refinement.model` | normalized SI contracts | `Geometry`, `Measurements`, `Targets`, `Assessment`, `RefinementPlan`, `RefinementResult`, `MeshRefiner` |
| `plumetools.refinement.assess` | measured stopping decisions | `assess`, `usable_window` |
| `plumetools.refinement.config` | standalone gas/reader configuration | `RefinementConfig`, `load_settings` |
| `plumetools.refinement.geometry` | actual edges, C/V, content hashes | `read_geometry`, `mesh_fingerprint`, `file_sha256` |
| `plumetools.refinement.foam` | standard, hyStrath, custom adapters | `Snapshot`, `read_snapshot` |
| `plumetools.refinement.plan` | immutable portable artifacts | `plan_refinement`, `write_plan`, `read_plan`, `verify_sources`, `plan_identity` |
| `plumetools.refinement.openfoam` | explicit local hex materialization | `HexRefiner` |
| `plumetools.refinement.report` | JSON/CSV export | `write_assessment`, `write_cells_csv` |
| `plumetools.refinement.__main__` | assess, plan, apply | `main` |

<!-- docs-example: portable-refinement -->
```python
import numpy as np
from plumetools.nozzleplume.gas import VhsSpecies
from plumetools.refinement.model import Geometry, Measurements
from plumetools.refinement.assess import assess

gas = VhsSpecies(name="Ar", mass_kg=6.63e-26, diameter_m=4.17e-10, omega=.74, t_ref_K=273)
geometry = Geometry(np.array([[0., 0., 0.]]), np.array([1e-9]), np.array([1e-3]), "mesh-sha")
measured = Measurements(np.array([1e18]), np.array([300.]), np.array([20.]),
    {"start_s": 1., "end_s": 2., "source": "declared", "basis": "mean_per_cell"}, "fields-sha")
assert assess(geometry, measured, gas).summary["action"] == "criteria_satisfied"
```
