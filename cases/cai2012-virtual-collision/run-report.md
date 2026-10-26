# Bounded implementation assessment, 2026-10-26

Audited source revision: `ad4cdc636b094b3bf1e0b601e7967b8b927370cc`.
The implementation is delivered as new sibling files; legacy Cai formulas,
base-case scripts, generation/run policy gates and unrelated work are unchanged.

**Historical smoke report:** the later `compare_exports.py` workflow supersedes
the geometry limitation below for the four archived Cartesian result sets. Their
geometry is recoverable after checking the saved mesh dictionary and exported
coordinates/counts. Updated three-method comparisons are indexed in
[`cai2012-kn-extension/results/virtual_collision/README.md`](../cai2012-kn-extension/results/virtual_collision/README.md).
Missing original `C/V` files alone no longer block those comparisons. DSMC tensor
and independent time-block uncertainty limitations remain.

The scientific hypothesis remains **unassessed**. Local discovery found 17
records: two templates, four direct-nozzle compact-result cases, two direct
configurations without results, six interface cases and three missing direct
Kn cases. Direct Kn 100, 0.1 and the additional Cai-family 0.01 are unavailable.
The Kn 1/10 graded exports have inconclusive stationarity. Kn 0.001/0.0001
direct cases retain their non-production resolution/time-step exclusions.
Interface cases have a different source and are excluded. None provides the
unaggregated geometry/moments required by the new matched-bin comparison.
Consequently no matched three-way accuracy plots or accuracy improvements are
claimed. Exact reasons and available diagnostics are in
[inventory.json](results/smoke/inventory.json),
[skipped.json](results/smoke/skipped.json) and
[assessment_skips.csv](results/smoke/assessment_skips.csv).

The predictive Kn=10 smoke used 64 flux samples, four batches, eight quadratic
hazard panels, 4 x 8 x 12 partner quadrature, seed 2026, a 12D-long/12D-radius
cylinder, and three finite central ring bins. All three bins passed moment
realizability checks. Mean proposed exposure was 0.04012; event-weighted mean
remaining exposure was 0.02749. The modeled ensemble uncollided source-number
fraction was 0.96198. These finite-sample diagnostics establish no accuracy
regime or uncertainty bound.

The sampled closed-volume number residual was zero; momentum residual magnitude
was below 2e-22 N and energy residual was 3.1e-19 W. All side, downstream and
return exits were counted. No return crossing happened in this small sample;
the model permits and separately counts backscatter. This benchmark checks
event bookkeeping, not a discretized spatial-divergence theorem.

| Same three bins and raw/derived outputs | Wall seconds | Core-seconds |
|---|---:|---:|
| Cai evaluation | 0.0302 | 0.0218 |
| Virtual preprocessing, cold cache | 0.3478 | 0.3185 |
| Virtual evaluation, prepared events | 0.2997 | 0.2973 |

Total virtual cost was 21.4 times Cai for this tiny request. Peak process-lifetime
RSS was 98.7 MiB, not an incremental-memory measurement. The host was WSL2 Linux
x86_64 (24 reported logical CPUs), with BLAS/OpenMP thread counts explicitly one.
This bounded run is compatible with the requested laptop resource ceiling; it
is not a benchmark performed on a specified eight-core machine. Credible DSMC
production resource accounting was unavailable, so no speed advantage over
production DSMC is established.

Independent refinements changed density-vector norms by 0.336% (samples), 0.070%
(hazard), 0.502% (partner quadrature), 0.0017% (spatial quadrature) and 0.303%
(domain). These are sensitivity measurements, **not convergence acceptance**.
Sampling and quadrature remain material relative to the small correction.
Changing the finite domain also changes the hazard-panel placement and sampled
events, so the domain difference is not an isolated tail-error bound. The
absolute omitted ballistic density has a Maxwellian tail bound; collision-rate
and corrected-field truncation errors are not thereby bounded. Raw moments,
tensors, realizability flags, batch standard errors and refinement diagnostics
are retained in NPZ/JSON alongside [the model-only plot](results/smoke/Kn10_model_only.png).
DSMC confidence intervals and tensor validation remain unavailable.

Verification completed:

- New kernel/model/assessment tests: **33 passed** (32.10 s), including the
  off-axis exact Cai regression and synthetic three-way plotting path.
- Affected existing analytical/gas/comparison/policy/interface tests:
  **189 passed, 3 failed** (273.50 s).
- All three failures reproduced with identical assertions from a pristine
  `git archive HEAD` copy of the audited revision (30.92 s): policy's
  `target_particles_per_cell` override expectation, interface manifest's
  `CONFIGURATION_INCOMPLETE` expectation and the corresponding interface guard.
  No legacy policy changes were made to suppress them.
- CLI help, inventory, bounded smoke and PBS shell syntax checks passed.
- No production DSMC, case generation or scheduler submission was performed.

Reproduce from the repository root in the configured Python environment:

```bash
python cases/cai2012-virtual-collision/assess.py smoke --output cases/cai2012-virtual-collision/results/smoke
python cases/cai2012-virtual-collision/assess.py full --output cases/cai2012-virtual-collision/results-full
python -m pytest tests/unit/test_virtualcollision_kernel.py tests/unit/test_virtualcollision_model.py tests/unit/test_virtualcollision_assessment.py -q
```

For manual HPC submission: `cd cases/cai2012-virtual-collision` followed by
`qsub runAssess.sh`. The environment override and numerical controls are in
[README.md](README.md) and [assessment.yaml](assessment.yaml). Full assessment
was not run locally. Read the [mathematical contract](../../docs/cai2012-virtual-collision.md)
before interpreting this reduced, no-fit, single-event model.
