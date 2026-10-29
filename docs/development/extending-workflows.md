# Extending the workflows

Where new work goes, and what has to change with it so the documentation and
the checks stay true. Module ownership and data flow are in
[Architecture](../plumetools-architecture.md); this page is the checklist.

## Rules the repository keeps

- **Generalise by adding, not editing.** The paper-specific packages
  (`cai2012`, `markelov1999`) are kept as they are; general code goes in a new
  package, as `nozzleplume`, `cai2012ext` and `hybridplume` did. A defect found
  in a frozen package is reported (see [Known gaps](known-gaps.md)), not fixed
  in passing.
- **Archived cases are never edited**, not even their READMEs.
- **One authoritative input per case.** Physics lives in `case.yaml` (or
  `hybrid.yaml`); every OpenFOAM dictionary and `0/` field is generated from it
  and gitignored. A new parameter is a new key in the family's loader, with a
  default, a unit in its name and a validation rule — not a literal in a script.
- **Outputs are pinned.** Generated dictionaries and fields of the paper
  families are compared byte for byte across changes; do not merge two
  families' writers to save code.
- **Independent checks stay independent.** Audits, `verify` modules and test
  oracles compute what they check without calling the code under check.

## Adding …

| Change | Where | Then also |
|---|---|---|
| a Knudsen number, particle level, pressure or sweep row | the family's `study.yaml` (data only) | regenerate; the derivation, manifest and runners follow |
| a derived quantity | the family's derivation (`cai2012.inflow`, `cai2012ext.direct_case.derive_config`, `interface_case.derive`, `nozzleplume.run`) — once | it appears in generation, the PBS steps and post-processing alike |
| a production criterion | `cai2012ext.feasibility` (direct) or `interface_assessment` (interface) | the manifest and `hpc label` pick it up |
| a mesh type | the family's mesh module: `plan(...)` → frozen `MeshPlan`, `render_*` per dictionary, `write_mesh_setup` as the only writer | read built meshes through `plumetools.mesh`, never a new parser |
| a rarefaction diagnostic | one registered function in `rarefaction/core.py`; new inputs read in `rarefaction/foam.py` | [Rarefaction diagnostics § Adding a diagnostic](../rarefaction-diagnostics.md#adding-a-diagnostic) |
| a family post-processing quantity | that family's `post.py`, reading through `foamio.volfield` and `mesh.patches` | translate `FoamDataError` into the family's error class |
| a new family | a new `plumetools/<family>/` package and `cases/<family>/` | see below |

## Keeping the documentation true

`python docs/_tools/check_api.py` and the strict site build fail when the code
and the reference disagree in the ways that can be checked mechanically:

- **A new module** must be listed, by dotted name, on its family's page
  (`docs/reference/api/<family>/index.md`); its own page is generated. A new
  **subpackage** or **top-level module** must also be assigned a family in
  `docs/_tools/hooks.py` (`FAMILIES` / `TOP_LEVEL`) and given an index page.
- **A renamed or removed public name** that a family page lists in its
  "Principal interfaces" column fails the check until the table is updated.
- **A removed module** with a hand-written page fails until the page goes.

What no check can catch, and the author of a change must update by hand:

| Changed | Update |
|---|---|
| a script's options or a `python -m` entry point | [CLI reference](../reference/cli.md) |
| a configuration key, default or rule | [Configuration reference](../reference/configuration.md) |
| an output file or column | [Outputs](../reference/outputs.md) |
| what a workflow does, needs or writes | the family's tutorial |
| a case's runnability, evidence or blockers | the [case catalog](../reference/case-catalog.md) |
| a new figure produced | replace its placeholder ([Figures](../assets/figures/README.md)) |

## A new case family

1. A package under `plumetools/` with its own loader (refusing other families'
   files by `model:`), derivation, dictionaries, `0/` writer and post-processing;
   shared I/O from `foamio` and `mesh`.
2. `cases/<family>/` with a `baseCase` (or equivalent template), a study file
   if it is a matrix, and runner scripts modelled on an existing family's
   (`Allmesh`, `Allrun` using `foamio.checkpoints plan`, `Allpost`, `Allclean`).
3. Narrow `.gitignore` entries for its generated dictionaries and trees,
   listed per family — never generalised to `**/constant/polyMesh/`.
4. Tests in `tests/unit` (and `tests/openfoam` for solver behaviour).
5. Documentation: a family page under `docs/reference/api/`, an entry in
   `hooks.FAMILIES`, a tutorial in `docs/tutorials/` and `mkdocs.yml`'s nav, a
   row in the case catalog, and sections in the CLI, configuration and outputs
   references. A short `cases/<family>/README.md` linking to the tutorial.
