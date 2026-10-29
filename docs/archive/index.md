# Archived cases

!!! danger "Historical record only — not supported quickstarts"
    These directories are **frozen**: nothing is added to, edited in or removed
    from them, and nothing here is fixed. Several produced results under labels
    that do not describe what was run. Do not treat any committed output in them
    as validated, and do not copy their scripts into new work. The authoritative
    record of what is wrong in each is [`cases/ARCHIVE.md`](../../cases/ARCHIVE.md);
    this page summarises it.

All five predate the `plumetools` package (2022–2023), name the MNF fork's
`dsmcFoam+` in their `controlDict`s, and would need that fork — whose build is
not recoverable — to run as committed ([Environment](../environment.md)).

## Two unrelated inflow lineages

| Lineage | Cases | Inflow |
|---|---|---|
| **Uniform radial** | `1d`, `2d-planar`, `2d-wedge`, `caseFoamEx` | `processInflowData.py` takes `U_mag T rhoN` on the command line and writes constant values around a 90° arc. This is the lineage that reproduces the Stewart & Lumpkin (2012) DAC isentropic-expansion case — the only V&V plots in the repository, and the only literature citation (the header of `2d-wedge/mesh/generateBlockMeshDict.m`, "Stuart" being a typo for Stewart) |
| **Analytical source flow** | `wake-cylinder/{5,25,100,475}psi` | the source-flow model that became `plumetools.sourceflow`, run on quadrilateral 2D meshes where its triangle assumption breaks |

## The cases

| Case | Period | What it was for | Known problems (IDs from `ARCHIVE.md`) |
|---|---|---|---|
| `1d` | 2022-02 – 2022-05 | the first isentropic-expansion validation against the digitised DAC curves (`n_nstar_radius.dat`, `T_Tstar_radius.dat`) | duplicate reference data (AR-15) |
| `2d-planar` | 2022-04 | a planar variant of the same expansion | `Allrun` is the stock OpenFOAM tutorial script: no `dsmcInitialise+`, no inflow generation (AR-12) |
| `2d-wedge` | 2022-04 – 2022-06 | the axisymmetric wedge V&V case of Stewart & Lumpkin (2012) | wrong isentropic exponent in `calcCaseParams.py` (AR-06); a digit-transposed sonic normaliser among three (AR-07); a 45 mm radius offset on its validation plots (AR-08); no `reconstructPar` before `postProcess`, and a `monitor` reading a log nothing writes (AR-11); a stale point count copied from another mesh (AR-13) |
| `caseFoamEx` | 2022-05 | a CaseFoam sweep of the wedge case over `nEquivalentParticles` (6 levels) × `maxRadialWeightingFactor` (4 levels), plus `baseCase` | the validation plot reuses a leaked loop variable, so all 24 curves are the same data under 24 labels (AR-05); `Allrun`/`Allclean` contain only a comment (AR-12); the repository's only git symlink (AR-14) |
| `wake-cylinder` | 2022-03 – 2023-02 | the first attempt at the AIAA 99-3455 cylinder-wake configuration, at four reservoir pressures, plus `3d` and `5psi_old` | quad-face centroid divided by 3 (AR-01); **the pressure sweep varied nothing** — `5psi` ran 100 psi, the other three 475 psi (AR-02); blocking `input()` calls (AR-03); three incompatible angular models (AR-04); byte-identical stale inputs (AR-09); hard-coded absolute paths and tutorial residue (AR-10); an undefined shell alias (AR-11) |

The active successor of the wake-cylinder idea is
[`markelov1999`](../tutorials/markelov1999.md), which is new work and does not
inherit these defects. The active successor of the source-flow lineage is
[`3d-inflow`](../tutorials/3d-inflow.md), which freezes the model's defects
deliberately so its extraction could be proven behaviour-preserving.

## Provenance of the historical results

Recoverable from committed artefacts: OpenFOAM build v1706
(`wake-cylinder/5psi/sample_out.txt`), Pointwise V18.5R2 (`.pw` headers),
ParaView 5.10.0 (`util/paraView/trace_lineplot.py`), the solver `dsmcFoam+`
(version not recorded). Python and library versions are unknown, so bit-level
reproduction is not achievable. The findings in `ARCHIVE.md` were verified by
reading the code and measuring the committed meshes at commit `7997ccc` (the
tag `pre-refactor-baseline` that `ARCHIVE.md` names is not present in this
repository).

## On Windows

`cases/caseFoamEx/2d-wedge/system/blockMeshDict` is a git symlink; without
symlink privileges it materialises as a text file and the tree stays
"modified". `git config --global core.symlinks true` before cloning avoids it.
It is not repaired, because that would modify an archived case.
