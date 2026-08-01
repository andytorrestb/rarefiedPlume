# Regression goldens for `cases/3d-inflow`

## What these are

`golden/3d_inflow_v0.npz` and `golden/boundary*.txt` are a frozen snapshot of the
inflow data produced by the **pre-refactor** code, captured at git tag
`pre-refactor-baseline` (commit `7997ccc`).

They exist so the extraction of `plumetools` can be proven not to have changed any
number.

> **These goldens prove the refactor did not change behaviour. They do NOT prove the
> behaviour is scientifically correct.** Findings SM-01…SM-10 — the ignored per-case
> stagnation pressure, the split stagnation temperature, the N₂/Ar species mismatch,
> the θ-from-+z convention, the un-normalised angular function — are all deliberately
> frozen *into* this snapshot. Agreement with it means "unchanged", never "right".

## How they were produced

```bash
git checkout pre-refactor-baseline
python tests/regression/capture_golden.py
```

`capture_golden.py` imports `cases/3d-inflow/processInflowData.py` **without editing
it**, working around two obstacles:

1. The legacy module calls `plumeSourceFlowModel()` at module scope (line 556), so the
   trailing call is stripped from the source before `exec` and invoked explicitly.
2. `readMeshStats()` shells out to `checkMesh`, which needs OpenFOAM. It is replaced by
   a stub. Only its `points` and `faces` keys are ever read by the legacy script, and
   both are taken **from the mesh files themselves** rather than hard-coded, so the stub
   cannot drift from the committed mesh.

`printInflowSurface` is wrapped to capture the in-memory `[rhoN, U, T]` dictionaries
before it writes, then delegated to so the `0/` files are produced as well.

## Contents

| File | Contents |
|---|---|
| `golden/3d_inflow_v0.npz` | `labels` (2044 face IDs, **in insertion order**), `rhoN` (2044,), `U` (2044, 3), `T` (2044,) |
| `golden/boundaryU.txt` | `0/boundaryU` as written by the legacy code |
| `golden/boundaryT.txt` | `0/boundaryT` |
| `golden/boundaryNumberDensity_Ar.txt` | `0/boundaryNumberDensity_Ar` |

Captured values, for orientation:

```
faces          : 2044
rhoN [min,max] : 2.702132e+19  5.580291e+21
|U|            : 788.164111 on every face   (the field is exactly radial: U = v_l * r_hat)
T              : 300.0 everywhere
all finite     : yes
```

**Face order is part of the contract.** `printInflow.py` writes values in
dict-iteration order and `dsmcFoam+` consumes them in patch-face order, so `labels` is
stored explicitly and must be reproduced, not just the value multiset.

## How the golden files are kept stable

These are byte-exact reference data, so nothing may rewrite them silently. Three things
enforce that:

* **`.gitattributes` marks `tests/regression/golden/** -text`.** Git therefore performs
  no line-ending conversion on checkout or commit. Without it, `core.autocrlf=true` —
  common on Windows — stored the files as LF but checked them out as CRLF, so the index
  and the working tree disagreed byte-for-byte while `git status` still looked clean.
  Any tool that rewrote a golden then surfaced it as a spurious modification.
* **`*.npz binary`**, rather than relying on git's content heuristic.
* **`capture_golden.py` writes LF explicitly.** The legacy code emits CRLF on Windows
  and LF on Linux for identical inputs, so capturing verbatim would have made the golden
  depend on the platform that captured it. Regeneration is now idempotent: re-running
  the capture on a clean tree produces byte-identical files.

## Two formatting quirks the comparison must tolerate

Both were discovered while capturing this snapshot. Neither is a behaviour change; both
make naive byte comparison non-portable.

1. **Line endings were platform-dependent.** Pinned as described above; the comparison
   still normalises them, so a golden captured before that fix does not fail.

2. **`boundaryT` values are emitted as `300`, not `300.0`.** `calculateT` returns the
   Python `int` literal `300`, and the writer interpolates `str(...)` of it. A config
   carrying `T0_K: 300.0` yields `300.0`. `U` and `rhoN` are unaffected — they are
   `numpy.float64`, whose `str()` matches `repr(float(x))`.

Consequently `test_3d_inflow.py` compares:

* **all structural lines byte-for-byte** after newline normalisation — header, dimensions,
  `internalField`, patch names, BC types, the face count, and the delimiters; and
* **the numeric payload numerically**, at `rtol=1e-12`, by parsing values back out.

That still catches every regression that matters (wrong count, wrong order, wrong patch
list, wrong dimensions, wrong values) without pinning the suite to one OS's newline
convention or to an int-vs-float formatting accident.

## Regenerating

Don't, unless you are deliberately re-baselining after an approved physics change. In
that case add a **new** file (`3d_inflow_v1.npz`) alongside this one and record which
`legacy:` flag flip it corresponds to — keep `v0` for provenance.
