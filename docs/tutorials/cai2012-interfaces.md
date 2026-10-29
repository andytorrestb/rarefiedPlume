# Kn extension: the interface track

!!! abstract "At a glance"

    | | |
    |---|---|
    | **Purpose** | At Kn = 1e-3 and 1e-4, remove the dense core from the DSMC domain and inject a continuum continuation of Cai's exit on a downstream sphere |
    | **Model** | Cai's exit continued as a Markelov-type angular source to a sphere chosen by Bird's breakdown parameter; truncated Cartesian mesh with the core carved out (`cai2012ext.interface*`) |
    | **Evidence** | Implemented and unit-tested; one tiny development case ran end to end on 1 and 2 ranks (mechanics only, 2026-10-22). **Nothing validated** |
    | **Status** | The overlap validation **cannot run as declared** (no admissible sphere at Kn = 0.01 for B = 0.02, 0.03, 0.05); the two production cases are generated, meshable, numerically unsuitable, and **refused at the run step** without passing evidence |
    | **Runs on** | Python for generation, reports and `--dict-only`; OpenFOAM v2512 + `plumeFieldInflow` for meshing (blockMesh + snappyHexMesh) and running |
    | **Start in** | `cases/cai2012-kn-extension/` |

The method, the equations, the audit, the decisions this track does **not**
make for you, and why the angular law does not fit Cai's exit are in
[The interface track, explained](../cai2012-kn-interface.md). Read §4 and §5
there before choosing any value below.

## Three kinds of case

`interfaces.yaml` declares three groups, generated into `Interface/` with
`interface-manifest.yaml`. Their inputs differ, and so does what each is
allowed to do:

| Group | Cases | `purpose` / `run_mode` | What it is for | Gate |
|---|---|---|---|---|
| `production` | `Kn0p001_np1x_sphere`, `Kn0p0001_np1x_sphere` | `production` / `production` | the scientific cases | refused by `hpc label` (exit 3) until `validation_evidence` points at a passing overlap for the **same coupling** |
| `validation` | `val_Kn0p01_*`, one per `bird_thresholds` entry (0.02, 0.03, 0.05) | — / `production` | truncated runs at Kn = 0.01 compared with the full-domain `cai2012` case | currently `interface-infeasible`: no admissible sphere |
| `smoke` | `smoke_Kn0p001_B0p05` | `development-smoke` / `exploratory` | mechanics only, on a shrunken domain (3 D box, ≤ 50 000 cells) | runs, marked NON-PRODUCTION |

The production entries also leave two decisions open by design: the core
half-angle (`core_half_angle_deg`) and the radius below which the source is not
trusted (`source_validity_min_radius_over_D`). The committed file sets
45°/60° and 0.7173 (the closure's r_a/D) as initial values; the explanation page
§4 says what each choice implies.

Statuses written by the assessment (`cai2012ext.interface_assessment`):

| Status | Meaning | `runCase.sh` |
|---|---|---|
| `interface-infeasible` | no admissible sphere, or a decision missing | `hpc guard` refuses (exit 2) |
| `unsuitable-for-production` | a production case misses a numerical criterion | `hpc label` refuses (exit 3) unless `run_mode: exploratory` |
| `pending-coupling-validation` | numerically sound, no passed overlap at this threshold | `hpc label` refuses unless exploratory |
| `production-candidate` | everything met, and passed evidence for the same coupling | runs, unmarked |
| `interface-validation` | an overlap run | runs if numerically sound (or exploratory), marked NON-PRODUCTION |
| `development-smoke` | mechanics only | runs, marked NON-PRODUCTION |

## Inputs

| File | Role |
|---|---|
| `interfaces.yaml` | the track: `overlay` (`baseCaseInterface`), `tree` (`Interface`), `manifest_name`, `resources` (memory ceiling), the three groups, and the validation `overlap` region, targets and uncertainty settings |
| `study.yaml`, `baseCase/case.yaml` | read with the direct-exit loader and policy, unchanged; each interface case is `baseCase` + `baseCaseInterface` |
| `<case>/interface.yaml` | written per case from `interfaces.yaml`; an **input**, compared like `case.yaml` — except `validation_evidence` and `label`, which are provenance (a change is *stale*, not a conflict) |

Every key, with allowed values: [configuration reference](../reference/configuration.md#interfaceyaml-and-interfacesyaml).

## Cheap checks (no OpenFOAM)

From `cases/cai2012-kn-extension/`, in WSL or Linux (or native Windows Python
for the Python steps):

```bash
python generate_interfaces.py --dry-run     # plan and every verdict; stages in .staging-*, keeps nothing
python generate_interfaces.py --report      # Interface/, interface-manifest.yaml, and
                                            # results/interface-report.yaml per case with the
                                            # (core angle x threshold) decision table
(cd Interface/smoke_Kn0p001_B0p05 && bash ./Allmesh --dict-only)
```

`generate_interfaces.py` has the same `--refresh`, `--regenerate` and
`--overwrite` semantics as [`generate_cases.py`](cai2012-kn-extension.md#generation-exactly):
conflicts are reported and nothing is changed unless a flag says so;
`--overwrite` destroys results. It does not read or write the direct-exit trees.
Exit status is 0 even when cases are infeasible or pending validation — those
are recorded verdicts.

One case at a time, through the same code the scripts call:

```bash
python -m plumetools.cai2012ext.interface_case report Interface/<case> --decision-table   # exit 3 if not meshable
python -m plumetools.cai2012ext.interface_case mesh-dicts Interface/<case>                # REFUSED (exit 3) if infeasible
```

The unit tests for this track run in seconds and need no OpenFOAM:

```bash
python -m pytest tests/unit/test_cai2012ext_interface.py \
                 tests/unit/test_cai2012ext_interface_case.py \
                 tests/unit/test_cai2012ext_overlap.py -q -p no:cacheprovider
```

## The development smoke run (OpenFOAM, minutes)

Recorded on the owner's laptop as about 37 000 cells and 2 × 10⁴ parcels, under
two minutes per run. It checks mechanics, never physics.

```bash
cd cases/cai2012-kn-extension/Interface/smoke_Kn0p001_B0p05
timeout 600 bash ./Allmesh
timeout 900 bash ./Allrun --serial
bash ./Allclean && timeout 900 bash ./Allrun -np 2
python -m plumetools.cai2012ext.interface_case check-decomposed .
bash ./Allpost && python Alldiagnose.py && python Allcompare.py && python Allconverge.py
python -m plumetools.cai2012ext.interface_case collect . --output results/interface-evidence.yaml
```

!!! note "Owner's laptop: two OpenFOAM installs"
    On the machine these were recorded on, a user-bin `decomposePar` from
    another OpenFOAM build fails every case; putting the v2512 binaries first
    fixes it: `export PATH="$FOAM_APPBIN:$PATH"` after sourcing OpenFOAM's
    `bashrc`. See [Troubleshooting](../guides/troubleshooting.md#a-different-decomposepar-is-first-on-path).

## Cluster sequence (you submit; nothing here does)

1. On the cluster: update the checkout, rebuild the inflow library if needed
   (`applications/dsmcBoundaryModels/Allwmake`), `pip install -e ".[test,plots]"`
   in the job environment.
2. `cd cases/cai2012-kn-extension && python generate_interfaces.py --report`.
3. Smoke first, from `Interface/smoke_Kn0p001_B0p05`:
   `qsub -N knx-imesh-smoke runMesh.sh`, then `qsub -N knx-irun-smoke runCase.sh`
   (exploratory; `dsmc.n_subdomains = 2` in its overrides), then `qsub runPost.sh`.
4. **Overlap cases: blocked.** Once a threshold set that admits spheres is
   chosen and the two open decisions are set, regenerate; then for each
   `val_Kn0p01_*`: `runMesh.sh`, `runCase.sh`, `runPost.sh`. `runPost.sh` runs
   `Alloverlap.py` against `overlap.baseline_case`
   (`../../../cai2012/Cases/Kn0p01_np1x`), or
   `python Alloverlap.py --baseline <path>`.
5. Collect: each `runPost.sh` writes `results/interface-evidence.yaml`;
   `python -m plumetools.cai2012ext.interface_case collect <case>` regenerates
   it. Copy those small files and `results/overlap.yaml` back; no fields needed.
6. Production: set `validation_evidence` in `interfaces.yaml` to the passing
   `overlap.yaml` for the same coupling (Bird definition and threshold, core
   angle, closure, tails, source problem), `generate_interfaces.py --refresh`,
   then `runMesh.sh`, `runCase.sh`, `runPost.sh` for the two `_sphere` cases.
   Without passed evidence `hpc label` refuses them.

`runViz.sh` refuses interface cases: no imagery spec exists for them.

## What counts as overlap evidence

`Alloverlap.py` (`python -m plumetools.cai2012ext.overlap <case> [--baseline PATH]`)
compares the truncated run with the full-domain benchmark in the configured
region and writes `results/overlap.yaml` with a verdict per quantity (n, speed,
T) against the targets, with block-bootstrap intervals and empty-bin
exclusions. A production case may rest on a `pass` only if the evidence records,
and matches, the model, `purpose: interface-validation`, the Bird definition
and threshold, the core half-angle, the tail treatment, the angular exponent
offset and the source problem. Kn, R, mesh, time step and weight may differ.
Evidence for another coupling reports `incompatible`; evidence that does not
say which coupling it tested (anything written before 2026-10-27) reports
`insufficient-provenance` — re-run `Alloverlap.py` to rewrite it.

## Outputs

| Path | Contents |
|---|---|
| `interface-manifest.yaml` | every case's derived plan and verdict |
| `<case>/results/interface-report.yaml` | the report, with the decision table under `--report` |
| `<case>/interface-identity.yaml` | what the mesh and parcels were built for; a restart on another mesh is refused |
| `<case>/results/interface-inflow.yaml` | the per-face inflow written by `interface_case inflow` |
| `<case>/results/overlap.yaml` | overlap verdict (validation cases) |
| `<case>/results/interface-evidence.yaml` | what `collect` brings back from a cluster run |

**Judging the output.** The measured fact the explanation page reports is that
Markelov's sonic-orifice angular law, imposed on Cai's supersonic exit, puts
0.36–0.45 of Cai's collisionless density on the axis (§6 there). Until an
overlap passes, nothing from a production interface case is a result.

## Resources

`interfaces.yaml` sets `resources.memory_ceiling_bytes: 1.5e12`; the
assessment compares the estimated memory and parcel counts with it. The smoke
case's size is the only measured one in the repository (above). Production
sizes are estimates in `interface-manifest.yaml`.

## Figure placeholder

!!! figure "Figure placeholder — overlap region, truncated vs full-domain density"
    **Subject:** n/n₀ in the overlap region (4 ≤ r/D ≤ 9.5, θ ≤ 15°) for a `val_Kn0p01_*` case against `cases/cai2012/Cases/Kn0p01_np1x`.
    **Source run:** a validation case that is not `interface-infeasible` — none exists with the committed thresholds.
    **Plotting command:** none verified; `Alloverlap.py` writes the numbers, not a figure.
    **Destination:** `docs/assets/figures/kn-interface-overlap.png`
    **Suggested caption:** "Overlap validation of the interface coupling at Kn = 0.01 (state threshold, core angle and verdict)."

## Related

- [The interface track, explained](../cai2012-kn-interface.md)
- [Kn extension, direct cases](cai2012-kn-extension.md) — generation semantics shared by this track
- API: [`interface`](../reference/api/cai2012ext/interface.md), [`interface_case`](../reference/api/cai2012ext/interface_case.md), [`interface_assessment`](../reference/api/cai2012ext/interface_assessment.md), [`interface_evidence`](../reference/api/cai2012ext/interface_evidence.md), [`overlap`](../reference/api/cai2012ext/overlap.md)
