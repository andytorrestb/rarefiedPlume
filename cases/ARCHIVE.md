# Archived cases — historical record, frozen

This file records which cases are **frozen** and what is known to be wrong inside them.
It lives at the `cases/` level deliberately: nothing is added to, edited in, or removed
from an archived case directory.

## Scope

| Status | Cases |
|---|---|
| **ARCHIVED** — frozen, do not modify | `1d` · `2d-planar` · `2d-wedge` · `caseFoamEx` (27 dirs) · `wake-cylinder` (6 dirs) |
| **ACTIVE** — under refactor | `3d-inflow` (reference) · `iss_solar_panels{,_2,_3}` · `solar_panel_particle_resolution_study/Fnum_*` |

Archived cases are one-off studies kept for historical record. They are **not** deleted,
**not** refactored, and **not** fixed. Their results stand as produced.

Two different inflow models live in this repository, and they are unrelated:

* **Uniform radial** (`1d`, `2d-planar`, `2d-wedge`, `caseFoamEx`) — `processInflowData.py`
  takes `U_mag T rhoN` on the command line and writes constant values around a 90° arc.
  This is the **validated** lineage: it reproduces the Stewart & Lumpkin (2012) DAC case
  and is the only one with real V&V plots.
* **Analytical plume source flow** (`3d-inflow`, `wake-cylinder/{5,25,100,475}psi`,
  `iss_solar_panels*`, `Fnum_*`) — parses the mesh and evaluates the source-flow model.

`wake-cylinder/{5,25,100,475}psi` are archived even though they use the source-flow
model, because they run it on quadrilateral 2D meshes where its triangle assumption
breaks (AR-01).

## Why this matters

**Do not treat archived output as validated output.** Several archived cases produced
results under labels that do not describe what was actually run. The findings below are
recorded so that nobody — including future automated comparisons — mistakes a committed
PNG or `.dat` for a verified result.

## Known findings in archived cases

Verified by reading the code and measuring the committed meshes at tag
`pre-refactor-baseline`. None of these are fixed.

| ID | Sev | Case | Finding |
|---|---|---|---|
| **AR-01** | S1 | `wake-cylinder/*` | **Quad-face centroid defect.** All inflow faces are `4(...)` but `calculateCentroid` divides the vertex sum by the literal `3.0`. Measured: code centroid mean \|r\| **0.3946 m** vs true **0.2960 m** — ratio 1.3331 = exactly 4/3. Every θ, φ, U and n on that patch is computed at the wrong location. Compounded by a hard-coded `r = 0.5` where the true inflow arc radius varies 0.2223–0.4033 m (an arc centred at (−0.37465, 0), radius 0.1524 m, per `5psi_old/processInflowData.py:20-22`). |
| **AR-02** | S1 | `wake-cylinder/{5,25,100,475}psi` | **The pressure sweep varied nothing.** All four hard-code `Po = 475*6894.75729` except `5psi`, which reads `physical_props['Po'] = 100` × 6894.76 — so the directory named **5psi ran 100 psi** and the other three all ran 475 psi. The `{100,25,475}psi` scripts differ from each other **only** in the value on line 558, which is then ignored. |
| **AR-03** | S2 | `wake-cylinder/{25,100,475}psi` | **Blocking `input()`** at line 387, inside `angularDependence`, called once per face (99 faces), preceded by `print(theta, phi)`. Entered `main` in commit `7f6ec54`. These scripts cannot run non-interactively. A second `input()` sits in `wake-cylinder/digitized_data/plot_digitized_data.py`. |
| **AR-04** | S2 | `wake-cylinder` vs `3d-inflow` | **Three incompatible angular models** in nominally-equivalent cases: `475psi` computes `f_theta * f_phi` then `return f_theta` (discarding f_phi); `5psi` uses a θ-only function with `asin(x/r)`; `3d-inflow` returns the product with `acos(z/r)`. Cross-case comparison within the sweep is not meaningful. |
| **AR-05** | S2 | `caseFoamEx/graphValidation.py:36-64` | **Leaked loop variable.** `oneGraph()`'s sampling loop leaves `curr_case` bound to the last case; the two plotting loops below reuse it instead of rebuilding the path. All 24 curves plot the same data under 24 different labels. Conclusions about nEquivalentParticles / maxRadialWeightingFactor sensitivity drawn from `digitized_vs_analytical_{n,T}.png` are not supported by those figures. |
| **AR-06** | S2 | `2d-wedge/calcCaseParams.py:44` | `nnt_Ma` uses exponent `1/(ga+1)`; the isentropic relation is `(T/T*)^(1/(γ−1))`, as used correctly in `graphCaseValidation.py:64`. At γ=5/3 that is 0.375 instead of 1.5. |
| **AR-07** | S2 | 3 validation scripts | **Three different sonic normalizers**: `rhoN_Ar / 8.377e20` (`2d-wedge`), `/ 8.773e20` (`util/graphValidation.py`), `/ 6.02e20` (`caseFoamEx`); temperature normalizers `1000.0` vs `800.0`. `8.377` vs `8.773` is a digit transposition — at least one is wrong. All undocumented. |
| **AR-08** | S3 | `2d-wedge` | **Radius inconsistency.** `mesh/generateBlockMeshDict.m:17` sets `rstar = 1.05` (measured true inflow \|r\| = 1.0499), but `mkSampleDict.py:9` uses `r_inlet = 1.005` and `graphCaseValidation.py:45` adds `+ 1.005`. Radial coordinates on the 2d-wedge validation plots are offset by 45 mm. |
| **AR-09** | S3 | `wake-cylinder/*` | **Stale identical inputs.** All five `inflowData.dat` are byte-identical (md5 `1bd45e7c54c7bc0973390b335494ae18`), inherited from the `5psi_old` lineage. Four of the five cases no longer read them, but each still ships `displayInflowData.py`, which does — so running it in `475psi` plots 5psi_old data under that case's name. |
| **AR-10** | S3 | `wake-cylinder/*/contourPlot.py` | Hard-coded absolute path `/home/andy/OpenFOAM/andy-v1706/run/rarefiedPlume/cases/wake-cylinder`, a hard-coded `(240, 920)` structured grid matching no mesh in the repository, and a `Circle((0,0), 5.0)` labelled "the turbine" — copy-paste residue from an unrelated tutorial. |
| **AR-11** | S3 | `2d-wedge/run.sh`, `monitor` | `decomposePar` + `mpirun` with **no `reconstructPar`** before `postProcess`; no output redirection, while `monitor` greps `log.dsmcFoam+` — a file nothing writes. `wake-cylinder/*/run.sh` calls an undefined shell alias `of1706`. |
| **AR-12** | S3 | `caseFoamEx/{Allrun,Allclean}`, `2d-planar/Allrun` | The two `caseFoamEx` scripts contain only a comment line. `2d-planar/Allrun` is the stock OpenFOAM tutorial script (`restore0Dir`, `runApplication blockMesh`, `postProcess -func streamFunction`) — no `dsmcInitialise+`, no inflow generation. |
| **AR-13** | S4 | `2d-wedge/displayInflowData.py:12` | Uses 33 points and radius 1.05, while the same directory's `processInflowData.py:10` uses `r = 132` and the mesh has 132 inflow faces. **33** is the `caseFoamEx/2d-wedge` mesh's inflow count, left behind by a directory copy. |
| **AR-14** | S3 | `caseFoamEx/2d-wedge/system/blockMeshDict` | The repository's **only** git symlink (mode `120000` → `../mesh/blockMeshDict`). On Windows without symlink privileges it materialises as a text file containing the path, leaving the working tree permanently dirty, and `git diff` on it fails with `Function not implemented`. Workaround: `git config core.symlinks true` before cloning. Not fixed, because fixing it means modifying an archived case. |
| **AR-15** | S4 | 60 files | `n_nstar_radius.dat` (30 copies, one hash) and `T_Tstar_radius.dat` (30 copies, one hash) — perfect duplicates of the digitised Stewart & Lumpkin (2012) reference data. A deduplicated copy is provided under `reference/` for **active** cases; the archived copies stay where they are. |

## Provenance of the reference data

The only literature citation anywhere in this repository is the header of
`2d-wedge/mesh/generateBlockMeshDict.m`:

> This program creates a blockMeshDict file to create an axisymmetric wedge mesh for
> reproducing the isentropic expansion V&V case in Stuart and Lumpkin (2012).
>
> Jonathan S. Pitt, Aegis Aerospace, Inc., for LENS
> Created: 21 Sep 2021 · Updated: 10 Nov 2021

Plot legends elsewhere read "DAC (Stewart, Lumpkin)", so "Stuart" is a typo for
**Stewart**. This reference covers the `1d`/`2d-wedge` isentropic-expansion validation
case and the digitised `n_nstar_radius.dat` / `T_Tstar_radius.dat` data.

**It does not cover the 3D plume source-flow model.** That model — the limiting angle,
the `(γ+0.41)/(γ−1)` angular exponent, the normalisation integral — has no citation in
any comment, commit message, or document in this repository. See
`docs/source-flow-model.md` for what could and could not be reconstructed from the code.

## Environment (partial — see `docs/environment.md`)

Recoverable from committed artefacts:

* OpenFOAM build **v1706** (`wake-cylinder/5psi/sample_out.txt`, line 8)
* Pointwise **V18.5R2** (`.pw` file headers)
* ParaView **5.10.0** (`util/paraView/trace_lineplot.py`, line 1)
* Solver: `dsmcFoam+` / `dsmcInitialise+` — the MNF fork; version not recorded

Python, numpy, scipy, pandas and matplotlib versions are **unknown and unrecoverable**.
Bit-level reproduction of archived results is therefore not achievable.
