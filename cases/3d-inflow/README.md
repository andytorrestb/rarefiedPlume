# cases/3d-inflow — legacy source-flow reference case

A rarefied plume expands from a hemispherical inflow surface (radius 0.5 m,
opening toward +x) into a 5 m vacuum box. The per-face inflow comes from the
**legacy** analytical source-flow model (`plumetools.sourceflow`), with every
known defect frozen as a `LEGACY` value in [`case.yaml`](case.yaml); standard
OpenFOAM v2512 `dsmcFoam` solves the field. It is the case the regression golden
in `tests/regression` is bound to.

**Not physically validated.** The model has no recoverable citation.

```bash
./Allmesh --dict-only     # mesh dictionaries from case.yaml; no OpenFOAM
./Allmesh --yes           # OVERWRITES constant/polyMesh: blockMesh + snappyHexMesh + verify
./Allrun --no-solve       # 0/ inflow fields only (needs a mesh)
./Allrun                  # full run; --serial, -np N
./Allclean --dry-run      # preview the reset
```

`constant/polyMesh` is gitignored here, so a fresh checkout has **no mesh**.
The historical Pointwise mesh is recoverable with
`git checkout f83233f^ -- cases/3d-inflow/constant/polyMesh`
(the `pointwise-mesh-3d-inflow` tag named in older text does not exist).

- Tutorial (inputs, checks, run, outputs, HA-06/HA-07, frozen defects, files):
  [`docs/tutorials/3d-inflow.md`](../../docs/tutorials/3d-inflow.md)
- The model and its defects: [`docs/source-flow-model.md`](../../docs/source-flow-model.md)
- The mesh pipeline: [`docs/mesh.md`](../../docs/mesh.md)
- Why the outer boundary reflects: [`docs/solver-compatibility.md`](../../docs/solver-compatibility.md)
