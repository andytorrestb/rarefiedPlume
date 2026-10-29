# cases/cai2012 — Cai & Wang 2012 collisionless circular plume

Cai, C., and Wang, L., *Numerical Validations for a Set of Collisionless Rocket
Plume Solutions*, J. Spacecraft and Rockets 49(1), 2012.

Argon expands from a circular nozzle exit (D = 0.2 m, S₀ = 2) **directly into
vacuum**; the exit is the DSMC inlet, through the `plumeFieldInflow` library.
Three Knudsen numbers (100, 0.1, 0.01) × four numerical-particle levels
(1×, 2×, 5×, 10×) = 12 cases, compared with Cai's collisionless analytical
solution. This is **not** a reproduction of Cai's DSMC numbers, and no result is
validated.

```bash
./generate_cases.py --dry-run     # what would be generated; writes nothing
./generate_cases.py               # Cases/<Kn>_<level>, manifest.yaml
./AllmeshCases                    # or --dict-only (no OpenFOAM)
./validate_cases.py               # the tree is the study (reads generated dictionaries)
./AllrunCases                     # --list / --index N / --case NAME; resumes
./AllpostCases                    # validation outputs, study table, imagery
```

Only `exit.knudsen` varies per case; everything else is derived from
[`baseCase/case.yaml`](baseCase/case.yaml), the one authoritative input.
Generated dictionaries, `Cases/`, `manifest.yaml` and `results/` are gitignored.

- Tutorial (inputs, cheap checks, production, resume, outputs, destructive flags):
  [`docs/tutorials/cai2012.md`](../../docs/tutorials/cai2012.md)
- Paper values, assumptions, mesh, time step, measured results:
  [`docs/cai2012-case.md`](../../docs/cai2012-case.md)
- Statistical health of the Kn = 100 case: [`../cai2012-health/`](../cai2012-health/README.md)
