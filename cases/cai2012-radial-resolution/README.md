# Cai exit radial resolution

Seven Kn values, a fixed mapped upstream inlet, and measured corrections only
where needed. Native OpenFOAM v2512 hyStrath `dsmcFoam+`/`dsmcInitialise+`, serial.

Edit [study.yaml](study.yaml), then:

```bash
python generate_cases.py --dry-run
python generate_cases.py
cd Cases/Kn100/r0_ppc20
bash ./Allmesh
bash ./Allrun
bash ./Allpost
```

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
