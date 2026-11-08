# `plumetools.radialstudy`

The [fixed-inlet radial study](../../../cai2012-radial-resolution.md) adapts Cai
exit physics to the portable measured-target API. See the
[tutorial](../../../tutorials/cai2012-radial-resolution.md).

::: plumetools.radialstudy
    options:
      members: false

| Module | Purpose | Principal interfaces |
|---|---|---|
| `plumetools.radialstudy.config` | strict study controls | `RadialStudyConfig`, `load_study` |
| `plumetools.radialstudy.physics` | existing source and Bird wrappers | `RadialModel`, `source_for_kn`, `predict_breakdown` |
| `plumetools.radialstudy.mesh` | actual streamtube geometry and targets | `RadialMesh`, `initial_mesh`, `RadialRefiner`, `render_block_mesh`, `read_radial_mesh` |
| `plumetools.radialstudy.case` | native dictionaries and sizing | `write_case`, `load_case`, `sizing` |
| `plumetools.radialstudy.generation` | seven initial cases and revisions | `generate_initial`, `revise_case` |
| `plumetools.radialstudy.run` | serial lifecycle and sampling provenance | `mesh_case`, `run_case`, `read_windows`, `collect_case` |
| `plumetools.radialstudy.post` | moments, diagnostics and comparisons | `assess_case`, `analyse_study`, `case_plan`, `post_case`, `project_moments` |
| `plumetools.radialstudy.report` | machine-readable and plotted outputs | `write_study_report`, `write_case_report` |
| `plumetools.radialstudy.__main__` | study/case CLI | `main` |

<!-- docs-example: radial-prediction -->
```python
from pathlib import Path
from plumetools.radialstudy.config import load_study
from plumetools.radialstudy.physics import source_for_kn, predict_breakdown

cfg = load_study(Path("cases/cai2012-radial-resolution/study.yaml"))
assert len(cfg.knudsen_numbers) == 7
prediction = predict_breakdown(source_for_kn(cfg, .0001))
assert prediction["status"] == "CROSSING_IN_DOMAIN"
assert abs(prediction["r_over_D"] - 14.992) < .002
```
