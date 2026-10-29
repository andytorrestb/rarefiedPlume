# cases/r-4d — R-4D exit-plane plume, with or without a plate

The exit plane of an R-4D nozzle (D = 0.377 m, 3187 m/s, 498 K, 1.3 × 10²² m⁻³;
Kn ≈ 4 × 10⁻⁴) is the DSMC inlet; the plume expands along +x into vacuum onto a
34.14 × 11.88 × 0.25 m plate 20 m downstream. It is the first case of the
general nozzle-plume family, `plumetools/nozzleplume`. At this exit density the
case is **under-resolved by necessity** (exit cell ≈ 220 mean free paths) and
says so; **no result is validated**.

| Directory | What |
|---|---|
| `baseCase/` | the case as described, 0.5 m background, ≈ 459 000 cells |
| `baseCaseHPC/` | the same plate on a finer, larger production mesh |
| `baseCaseFreePlumeHPC/` | no plate (`plate: null`), production mesh ([README](baseCaseFreePlumeHPC/README.md)) |

```bash
cd baseCase
./Allmesh --describe     # what case.yaml implies and costs; writes nothing
./Allmesh && ./Allrun && ./Allpost
./Allviz                 # slices (pvpython)
python -m plumetools.iteration baseCase <new-dir> --note "why"   # from cases/r-4d: a fresh mesh iteration
```

Change a mesh in a **new iteration**, never by re-meshing a case that holds results.

- Tutorial: [`docs/tutorials/r-4d.md`](../../docs/tutorials/r-4d.md)
- What was given and assumed, mesh, cost, outputs, what has been checked
  (moved from this README): [`docs/r-4d-case.md`](../../docs/r-4d-case.md)
- Diagnostics and the remesh workflow:
  [`docs/rarefaction-diagnostics.md`](../../docs/rarefaction-diagnostics.md)
