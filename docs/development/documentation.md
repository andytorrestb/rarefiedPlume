# Writing and building the documentation

The site is MkDocs with the Material theme, built locally into `site/`
(gitignored). Nothing deploys it: there is no Pages configuration, no hosting
and no publication step, and CI builds and checks it without publishing.

## Environment

The documentation has its own pinned dependencies in `docs/requirements.txt`
and its own environment. It needs **no** OpenFOAM, ParaView, MPI, CaseFoam,
numpy or scipy: the API pages are read from source, not imported. It does need
**Python ≥ 3.11** (Material 9.8, mkdocstrings 1.x and mkdocstrings-python 2.x
require it), independently of the package's own ≥ 3.9 floor.

=== "Linux / WSL / macOS"

    ```bash
    python3 -m venv .venv-docs
    .venv-docs/bin/python -m pip install -r docs/requirements.txt
    ```

=== "Windows (PowerShell or cmd)"

    ```bash
    py -3 -m venv .venv-docs
    .venv-docs\Scripts\python -m pip install -r docs/requirements.txt
    ```

`.venv-docs/` is gitignored. `pyproject.toml` and `Allsetup` are unchanged, so
a simulation install never acquires these packages.

## Commands

From the repository root, with the documentation environment's Python:

```bash
python -m mkdocs serve                    # http://127.0.0.1:8000, rebuilds on save
python -m mkdocs build --strict           # into site/; any warning fails
python docs/_tools/check_api.py           # module inventory and documented names
python docs/_tools/check_links.py         # internal links and anchors, site and GitHub views
python docs/_tools/run_examples.py        # the allowlisted examples (needs plumetools' own deps)
```

`run_examples.py` imports `plumetools`, so it needs numpy, scipy and PyYAML —
install the package (`pip install -e .`) into a separate environment, or the
documentation environment if you do not mind the extra packages. CI runs it in
its own job. A full strict build takes about half a minute on a laptop.

## How the site is assembled

| Piece | File | Does |
|---|---|---|
| configuration | `mkdocs.yml` | theme, Markdown extensions, explicit navigation, mkdocstrings options |
| hooks | `docs/_tools/hooks.py` | generates one API page per module at build time and fills the "API reference" navigation from the source tree; rewrites repository-relative links to GitHub source URLs after checking the target exists |
| docstring bridge | `docs/_tools/griffe_rst.py` | a Griffe extension translating the docstrings' reStructuredText to Markdown in memory |
| checks | `docs/_tools/check_api.py`, `check_links.py`, `run_examples.py` | see below |
| styles, scripts | `docs/stylesheets/extra.css`, `docs/javascripts/mathjax.js` | figure placeholders, MathJax configuration |

`docs/_tools/` and `docs/requirements.txt` are excluded from the site
(`exclude_docs`).

### API pages

A module gets a generated page (`::: plumetools.<module>`) unless a
hand-written one exists at `docs/reference/api/<family>/<name>.md`. mkdocstrings
reads the package statically (`allow_inspection: false`): no import, no
optional runtime. Each family's `index.md` is hand-written and must list every
module of the family; it carries the companion notes — units, conventions,
examples — that docstrings leave implicit.

### The docstring bridge

The package's docstrings are Google-style sections written in
reStructuredText. `griffe_rst.py` translates, at render time and without
touching the source: `.. math::` blocks (to `$$` blocks, with several rows as
`aligned`/`gathered` and their `\tag`s shown as labels, since arithmatex wraps
every block in `\[…\]` where `align` cannot nest), `:math:` roles, `:func:`,
`:class:`, `:mod:`, `:meth:`, `:attr:`, `:data:`, `:exc:`, `:obj:` roles (to
cross-references when the target is rendered, otherwise to inline code),
`.. code-block::` and `text::` literal blocks, `.. warning::`/`.. note::`,
simple `=====` tables, the one grid table (as preformatted text), and section
underlines.

Known limits: other directives are left as written (none occur today);
cross-references to private helpers become inline code; a multi-parameter
entry such as `cell_width, cell_height: …` in an `Args:` section renders as one
row. Griffe's warnings about missing type annotations and such entries are
turned off in `mkdocs.yml` (`warn_missing_types`, `warn_unknown_params`)
because they are source style, not build errors — every other warning fails
the strict build.

### Links

Write ordinary relative links. Inside `docs/`, MkDocs validates them and their
anchors. A link out of `docs/` (to a case file, a script, a module's source) is
checked to name a path **git tracks** and rewritten to its GitHub URL on
`DOCS_SOURCE_REF` (default `main`); a link to a missing or untracked path — a
local result, a build product — fails the strict build, so a link that works
only in one checkout fails there rather than in CI. In the repository's other
Markdown, `check_links.py` lists links into gitignored outputs as notes (two
virtual-collision records link to local review results) rather than failures. Because the pages are also read on GitHub, anchors must work
under both slug rules: avoid em dashes and non-ASCII letters in headings that
are link targets. `check_links.py` checks both. Links to *generated* module
pages work only in the built site; on GitHub, link to the family page.

Equations use `$…$` and `$$…$$` (arithmatex, generic mode). MathJax is loaded
from a CDN (`cdn.jsdelivr.net`, pinned version), so equations need a network
connection when a page is viewed; everything else works offline.

## Checks

| Check | Command | Catches |
|---|---|---|
| strict build | `mkdocs build --strict` | broken page links and anchors, nav entries without pages, unresolved cross-references, mkdocstrings errors, repository links to missing paths |
| API inventory | `check_api.py` | modules added or removed without their family page, hand-written pages for removed modules, `:::` targets that do not resolve, listed interface names that no longer exist |
| links | `check_links.py` | broken links and fragments in the built HTML (docstrings, navigation included); repository Markdown links and anchors under GitHub's rules |
| examples | `run_examples.py` | allowlisted page examples and read-only preparation commands, each in a temporary copy, with a time limit, verifying they leave the copied tree as it was where they claim to |

External URLs are not checked (to keep the build deterministic and offline).

## Figures

Figures are **placeholders** until the owner adds real ones: an admonition
`!!! figure "Figure placeholder — <subject>"` stating the subject, the source
run, a verified plotting command if one exists, the destination under
`docs/assets/figures/` and a suggested caption. No placeholder contains an
image link. The list is in [Figures to add](../assets/figures/README.md). To
add a figure, save it at the destination and replace the admonition with
`![caption](../assets/figures/<file>)` and a caption that names the run.

## Tooling status

Pinned: MkDocs 1.6.1, Material for MkDocs 9.8.0, mkdocstrings 1.0.8,
mkdocstrings-python 2.0.10, pymdown-extensions 12.3, Black 26.10.1 (signature
formatting), with every transitive package in `docs/requirements.txt`. Two
facts to watch, as of 2026-10-28:

- Material for MkDocs 9.8.0 is a security release whose package metadata marks
  the project inactive, with an end-of-life date of **5 May 2027**; its
  maintainers point to Zensical as the successor.
- MkDocs 2.0 pre-releases exist and are not compatible with Material 9.x (which
  requires `mkdocs<2`); the pins keep 1.6.1. Material prints a notice about
  this at build time (`NO_MKDOCS_2_WARNING=1` silences it).

Neither affects building today. Moving to another generator would mean
re-checking the hooks, the Griffe extension and the admonition syntax.
