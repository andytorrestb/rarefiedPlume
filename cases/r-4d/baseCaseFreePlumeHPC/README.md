# R-4D free-plume HPC case

This case uses the R-4D exit state and production mesh/run settings from
`baseCaseHPC`, but has no plate (`plate: null`). The plume expands into a
50 × 80 × 32 m³ vacuum box; particles reaching any outer boundary are removed.
The refinement regions follow the plume out to 30 m.

The configuration estimates about 7.35 million cells and uses 12 MPI
subdomains. It is an HPC production configuration, not a laptop smoke test.
Its results are not evidence of numerical convergence or physical accuracy.

From this directory on an HPC node with OpenFOAM v2512, the repository's
`plumeFieldInflow` library built, and `plumetools` installed:

```bash
./Allmesh --describe
./Allmesh
./Allrun
./Allpost
```

`Allmesh --describe` validates the configuration and reports estimated resource
use without writing files. Regenerate dictionaries and mesh from `case.yaml` in
the target environment; do not rely on a laptop-generated mesh being portable.
`Allrun` resumes an interrupted simulation. `Allpost` writes the centreline,
rarefaction diagnostics, and the configured flow-field plots. Plate-load files
are not produced for this free-plume case.

The OpenFOAM mesh/solver run and ParaView visualization have not been validated
for this case. The config loader, derived mesh/run plan, and visualization
specification should be validated locally before HPC execution.
