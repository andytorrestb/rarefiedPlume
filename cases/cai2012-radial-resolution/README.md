# Cai exit radial resolution

Seven Kn values, a fixed mapped upstream inlet, and measured corrections only
where needed. Native OpenFOAM v2512 hyStrath `dsmcFoam+`/`dsmcInitialise+`, serial.

Edit [study.yaml](study.yaml), then:

```bash
python generate_cases.py --dry-run
python generate_cases.py
bash ./Allrun --dry-run
bash ./Allrun
# After an interruption, or after generating more revisions:
bash ./Allrestart
```

The study-level `Allrun` builds/verifies each mesh and runs all generated
`Cases/Kn*/rN_ppc*` cases one at a time, in decreasing Kn and increasing revision
order. It checks the whole queue first and stops at the first solver failure.
`Allrestart` skips finished cases, resumes complete checkpoints, and starts
untouched cases. Neither script deletes results. Both accept `--cases OTHER_TREE`
and `--dry-run`. The default densest case exceeds the step budget: inspect the
dry run, then pass `--allow-budget-exceedance` explicitly if that cost is intended.

For a single case, use its own `Allmesh`, `Allrun` (`--resume` to continue), and
`Allpost` scripts. Postprocessing and refinement remain explicit steps.

Ensure this checkout is on PYTHONPATH and use a Python with numpy/scipy/PyYAML.
The [tutorial](../../docs/tutorials/cai2012-radial-resolution.md) includes the
owner's WSL setup example, correction commands, restart and output locations.
The [method](../../docs/cai2012-radial-resolution.md) explains assumptions and
evidence. `baseCase/` contains wrapper templates; generated cases contain their
own complete dictionaries and embedded physical reference.

Initial predicted target: `min(lambda_predicted, 0.05D)`. Final measured target:
largest actual edge `<= 0.8 lambda_measured`, mean parcels/cell `>= 20`, plus
separate timestep and sampling checks. Passing cases need no further mesh sweep.
Revisions retain parents and fix inlet/domain. Full scientific runs remain unrun.
