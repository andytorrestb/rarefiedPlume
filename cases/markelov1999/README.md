# cases/markelov1999 — AIAA 99-3455, a flat plate in a cylinder wake

Lumpkin, F. E., Stewart, B. D., and Markelov, G. N., *Study of 3D Rarefied Flow
on a Flat Plate in the Wake of a Cylinder*, AIAA 99-3455, 1999.

A rarefied N₂ plume from a small orifice passes a finite cylinder and impinges
on a flat plate 6 in behind it, at reservoir pressures 5, 25, 100 and 475 psi.
The inflow is a paper-faithful source-flow model (`markelov1999_axisymmetric`),
injected per face by `plumeFieldInflow`. **No result here is validated** against
DAC or experiment.

```bash
./generate_cases.py --dry-run     # needs CaseFoam: pip install -e ".[cases]"
./generate_cases.py               # Cases/gap06in/p005psi ... p475psi
./AllmeshCases                    # blockMesh + snappyHexMesh + verification
./AllrunCases                     # sequentially -- not CaseFoam's concurrent Allrun
./AllpostCases                    # case summaries, study table, pressure plots
```

[`baseCase/case.yaml`](baseCase/case.yaml) is the authoritative input; every
OpenFOAM dictionary and `0/` field is generated from it and gitignored.

- Tutorial: [`docs/tutorials/markelov1999.md`](../../docs/tutorials/markelov1999.md)
- Paper values, model, assumptions, flux verification, surface pressure:
  [`docs/markelov1999-case.md`](../../docs/markelov1999-case.md)
- The custom inflow model: [`docs/plume-field-inflow.md`](../../docs/plume-field-inflow.md)
