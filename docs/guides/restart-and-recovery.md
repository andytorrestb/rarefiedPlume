# Restart and recovery

Re-running a case's `Allrun` **continues** it; it never starts over and never
deletes a time directory. This page explains how the decision is made, what to
do when a case is refused, how to run further, and how to change a mesh
without damaging results. It applies to the `Allrun` of `cai2012`,
`cai2012-health`, `cai2012-kn-extension`, `markelov1999`, `r-4d` (all three
variants) and every solver stage of the hybrid `Allrun`. `3d-inflow`'s `Allrun`
has no resume logic.

## How a resume is decided

The generated `controlDict` says `startFrom latestTime`, so `dsmcFoam` itself
resumes from the newest time directory each rank has. What `Allrun` decides is
what the **supporting steps** must do — two of them would destroy the state
being resumed from: `dsmcInitialise` rebuilds the parcel cloud, and
`decomposePar -force` deletes `processor*/`. `Allrun` asks

```bash
python -m plumetools.foamio.checkpoints plan . --np N
```

which inspects every rank and prints shell assignments:

| Variable | Meaning | `Allrun` does |
|---|---|---|
| `RESUME` | the newest checkpoint complete everywhere it is needed (empty: fresh start) | skips `runInflow.py`/`0/` generation and `dsmcInitialise` when set |
| `DECOMPOSE=none` | the processor directories hold `RESUME` and their count matches `N` | no `decomposePar` |
| `DECOMPOSE=latest` | resuming from the case directory, or the rank count changed | `decomposePar -force -latestTime` |
| `DECOMPOSE=all` | fresh start | `decomposePar -force` |
| `RECONSTRUCT_FIRST=1` | the processors hold a newer complete state than the case directory, and re-decomposition would delete it | `reconstructPar -newTimes` first |
| `REACHED_END=1` | `RESUME` is at `endTime` (with the half-step slack `Time::run()` uses) | reports "already reached endTime" and exits 0 |

After a parallel solve `Allrun` runs `reconstructPar -newTimes`, so earlier
writes are not reconstructed again. `fieldAverage` keeps its accumulators in
each time's `uniform/`, so averages span the whole sampled period across legs.

**Complete** means, on every rank: `T/uniform/time` exists, every field file in
`T/` and `T/lagrangian/<cloud>/` ends with OpenFOAM's end divider, `T/` holds
every file any other rank holds at `T`, and a rank that had a parcel cloud at
the previous time still has one. A time that fails any of these is
*incomplete* — what a job killed mid-write leaves. Inspect without running:

```bash
python -m plumetools.foamio.checkpoints report <case>
```

## When a run is refused

If the newest time is incomplete on some rank, `plan` exits **3**, prints why on
stderr, and `Allrun` stops before changing anything. The tool's own recovery
instruction is: move the incomplete time aside on every rank and in the case
directory, then run again:

```bash
mv processor0/<t> processor0/<t>.incomplete      # ... and every other processorN that has it
mv <t> <t>.incomplete                             # if the case directory has it too
./Allrun                                          # resumes from the previous complete time
```

A renamed directory is no longer a time directory, so nothing reads it; delete
it once the resumed run has passed it. A serial run is also refused when the
processor directories hold a newer time than the case directory: reconstruct
first (`reconstructPar -newTimes`) or run with the same `-np`.

## Running further

```bash
# raise the end time in case.yaml (e.g. dsmc.sampling_domain_transits, or end_time_s), then
./Allmesh --dict-only      # rewrites controlDict; does NOT re-mesh
./Allrun                   # continues from the latest time
```

`0/` is **not** regenerated on a resume: the solver reads the boundary fields
and parcels from the resume time, so a physics edit in `case.yaml` would appear
to take hold while changing nothing. To change physics — or the particle weight
— start a fresh run (clean, or better, a fresh iteration).

## Changing the rank count

`./Allrun -np M` on a case decomposed for `N`: `plan` reports `DECOMPOSE=latest`
and, if the processors hold newer state, `RECONSTRUCT_FIRST=1`; `Allrun`
reconstructs, then decomposes only the resume time for `M` ranks.

## Changing the mesh

`dsmcFoam` does not refine a mesh while it runs (the owner checked v2512: it
builds a static mesh and ignores `dynamicMeshDict`), and
re-meshing a case in place would leave results addressing cells that no longer
exist. Every family's `Allmesh` refuses to re-mesh a case with results unless
`--force` is given. The supported path is a **new iteration**:

```bash
python -m plumetools.iteration <case> <new-case> --note "why"
# edit <new-case>/case.yaml, then ./Allmesh && ./Allrun from t = 0
```

It copies inputs only and checks that no run state came with them
([r-4d tutorial](../tutorials/r-4d.md#change-the-mesh-a-fresh-iteration)).

## Family differences

| Family | Notes |
|---|---|
| `cai2012-health` | `AllrunCases` prunes reconstructed parcel clouds from all but the newest time after each case; `processor*/`, which a resume reads, is untouched. `--keep-clouds` disables it |
| `cai2012-kn-extension` | `runCase.sh` is resubmitted by hand after a walltime kill; `hpc finished` stops a finished case from being re-timed. Interface cases also check `interface-identity.yaml`: a restart on a mesh built for another coupling is refused |
| `cai2012-hybrid` | each stage resumes independently; a stage at `endTime` is skipped; nothing is deleted — to start over, remove the case directory. The `dsmcFoam` restart insertion burst recorded by the owner is unexplained |
| `3d-inflow` | no resume: rerun from scratch after `./Allclean` |
| Generated studies | a fix to a template's runner scripts does not reach generated cases; the Kn extension's `generate_cases.py --refresh` copies them into cases whose physics is unchanged |

## What never happens

Nothing in `Allrun`, `checkpoints` or the PBS wrappers deletes, moves or edits
a time directory. The commands that do delete are explicit and named:
`Allclean`, `AllcleanCases` (`--purge` for the whole tree), `generate_cases.py
--overwrite`, `Allmesh --force` (the mesh) and `3d-inflow/Allmesh --yes`.
Preview the cleaners with `--dry-run`.
