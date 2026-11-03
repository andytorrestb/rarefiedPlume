# Shared workflow modules

Two top-level modules every study family uses.

| Module | Purpose | Principal interfaces |
|---|---|---|
| [`plumetools.casemanifest`](casemanifest.md) | the case list a study's runner scripts work through: read whole, checked, then listed; CLI | `load`, `case_paths`, `select`, `ManifestError`, `main` |
| [`plumetools.iteration`](iteration.md) | start the next mesh iteration of a case: its inputs and nothing from its run; CLI | `prepare`, `Iteration`, `IterationError`, `resumable_state`, `CASE_INPUTS`, `main` |

## `casemanifest`

The runners of `cai2012`, `cai2012-health` and `markelov1999` take their cases
from `manifest.yaml` — the record of what was generated — not from a glob. The
whole manifest is checked before a single path is printed, so a manifest broken
at its fourth entry no longer launches the first three (defect 3 of the
2026-10-27 refactor). The Kn extension lists through
`plumetools.cai2012ext.policy`, which applies the suite policy as well.

<!-- docs-example: casemanifest -->
```python
from plumetools.casemanifest import ManifestError, case_paths, select

paths = case_paths({"cases": [{"path": "Cases/Kn100_np1x"}, {"path": "Cases/Kn0p1_np1x"}]})
assert select(paths, ["Cases/Kn0p1_np1x"]) == ["Cases/Kn0p1_np1x"]
try:
    case_paths({"cases": [{"path": "Cases/A"}, {"path": "Cases/A"}]})   # a duplicate is refused
except ManifestError:
    pass
```

## `iteration`

`prepare(source, destination, *, note=None, include=())` makes `destination` a
fresh case holding `source`'s inputs: `case.yaml`, `viz.yaml`, `open.foam` and
the `All*` scripts (`CASE_INPUTS`), plus any `include` paths. No mesh,
dictionary, `0/`, time or processor directory, parcel or averaging state is
carried over, and the new directory is checked for all of them
(`resumable_state`). It writes `iteration.yaml` and `provenance/` (the parent's
`case.yaml`) and only reads the source. It refuses a destination that is not
new and empty. This is the supported way to change a mesh after a run;
re-meshing a result-bearing case in place leaves results addressing cells that
no longer exist. Usage: [r-4d tutorial](../../../tutorials/r-4d.md#change-the-mesh-a-fresh-iteration).
