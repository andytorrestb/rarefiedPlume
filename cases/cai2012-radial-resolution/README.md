# Cai exit radial resolution

Seven argon exit Kn values, a fixed mapped upstream inlet, a 2.5-degree
half-angle, and an outlet at 10 exit diameters. Each case starts with 100 cells.
Native OpenFOAM v2512 hyStrath dsmcFoam+ runs one case at a time in WSL.

With this checkout on PYTHONPATH, Python with numpy/scipy/PyYAML, and OpenFOAM
sourced, edit [study.yaml](study.yaml), then:

```bash
python generate_cases.py --dry-run
python generate_cases.py
bash ./Allrun --dry-run
bash ./Allrun
# Continue after an interruption:
bash ./Allrestart
# Follow the persistent summary in another terminal:
tail -f Cases-adaptive/adaptive-study.log
```

The new defaults use **Cases-adaptive**, preserving the earlier **Cases** tree.
Both runners accept `--cases OTHER_TREE`. Legacy trees keep their fixed-duration
behavior; use `bash ./Allrestart --cases "$PWD/Cases"` to continue one.

After startup, the controller checks complete native write windows. Four stable
windows allow a numerical sizing decision. Passing cases stop early. Trusted
failures prepare a corrected revision, record the results/settings, clean the
superseded attempt's native outputs, and restart from initial conditions.
Unsettled measurements continue sampling on the same cloud.

`Cases-adaptive/adaptive-study.log` is append-only across retries and command
restarts. It records Kn, attempt, time, cells, particle weight, timestep, measured
requirements and actions. Small per-attempt records and plans remain; accepted
and crashed runs are retained. Retry/resource/sampling limits stop with an honest
non-passing outcome. `--allow-budget-exceedance` explicitly overrides step/parcel
budgets, but never the retry limit or mesh cell cap.

The 100-cell mesh is exploratory. Final requirements remain largest actual edge
`<= 0.8 lambda_measured`, mean parcels/cell `>= 20`, and separate timestep checks.
Numerical sizing acceptance does not certify statistical independence or flow
convergence. Dense cases and the wider angle can require large final meshes.
The geometry is a square expanding tube, with 2.5 degrees in each coordinate
plane, rather than a circular axisymmetric wedge.

The [tutorial](../../docs/tutorials/cai2012-radial-resolution.md) gives WSL setup,
recovery and report commands. The [method](../../docs/cai2012-radial-resolution.md)
explains physics and evidence. Dictionary generation never starts a solver.

The native hyStrath measurement code requires `writeControl runTime` and a
`writeInterval` in seconds. A step-count interval can overflow its 32-bit
averaging-step counter during field construction. The generator uses the frozen
sampling cadence in seconds so output and measurement windows stay aligned.
