# Post-processing

What each family's post-processing reads and writes, the separate ParaView
interpreter imagery needs, and the paths that work without it.

## Per family

| Family | Command | Reads | Writes | Needs |
|---|---|---|---|---|
| cai2012 | `Cases/<c>/Allpost`; `AllpostCases` | time-averaged fields + `C`, `V` | `results/` centreline, plane, metrics, figures; study table | `postProcess`, Python; matplotlib for figures |
| cai2012-health | `AllpostCases` | as cai2012, every frame | per-case audit; matrix sheets, curves, overlap | + ffmpeg (sheets), `pvpython` (imagery) |
| Kn extension | `AllpostCases`, `runPost.sh` | as cai2012 | + diagnostics, comparisons, stationarity; study table | + `pvpython` for imagery |
| markelov1999 | `Allpost`; `AllpostCases` | `log.checkMesh`, `log.dsmcFoam`, `fDMean`, `dsmcRhoN` | `case-summary.json`; study table and pressure plots | `writeCellCentres` (optional), matplotlib |
| r-4d | `Allpost`, `Allviz` | averaged fields, `C`, `V` | plate loads, centreline, diagnostics fields and summary, images | `postProcess`; `pvpython` for images |
| hybrid | `Allrun <case> post-dsmc` | the DSMC and CFD results | `results/*.json`, `summary.json`, figures | Python, matplotlib |
| any dsmcFoam case | `python -m plumetools.rarefaction <case> --settings FILE` | a reconstructed sampled time | diagnostic fields, summary | Python |

Everything reads the **time-averaged** fields (`*Mean`); a single-step field is
one step's parcels and measures shot noise. Check the parcel count in
`log.dsmcFoam` has levelled off before reading any of it.

## Imagery needs ParaView's own interpreter

`plumetools/viz/render_slices.py` imports `paraview.simple`, which exists only
inside ParaView's `pvpython`. **`paraview.simple` is not pip-installable**, and
running the script with `python` fails at the import with a message saying so.

```bash
pvpython --force-offscreen-rendering plumetools/viz/render_slices.py <case> [<case> …] --spec <viz.yaml>
```

- `--force-offscreen-rendering` keeps windows off the screen; the owner's
  cluster environment also sets `VTK_DEFAULT_OPENGL_WINDOW=vtkOSOpenGLRenderWindow`
  for headless nodes.
- `plumetools` is found beside the script; it need not be installed into
  ParaView's Python. PyYAML and numpy must be importable by `pvpython`: the
  paraview.org binaries ship them; a distribution `pvpython` uses the system
  Python's; a ParaView built on site may have neither, in which case the script
  falls back to `$PV_VENV` (or the active virtual environment) of the same
  Python `major.minor`. Check with
  `pvpython -c "import sys, yaml, numpy; print(sys.version)"`.
- `PVPYTHON_PATH` names a `pvpython` that is not on `PATH`; the study scripts
  use it.
- On headless or WSLg displays `pvpython` can die at exit with `GLXBadContext`
  after every image is written; `render_slices.py` reports its own exit status,
  so that does not turn a successful render into a failure.

Verified ParaView versions, and what they draw differently:
[Slice imagery § ParaView versions](../viz-slices.md#paraview-versions).

## Without ParaView

Every family's numerical post-processing works without ParaView:
`AllpostCases --no-viz` (cai2012, health, Kn extension, markelov1999),
`Allpost --no-viz` (r-4d). Without `pvpython` on `PATH`, imagery steps are
skipped with a one-line note (r-4d's `Allpost --viz` makes that an error
instead). `open.foam` in each case opens it in a ParaView GUI elsewhere.

## Videos and dissection (ffmpeg, no ParaView)

```bash
python plumetools/viz/video.py <case>/results/viz                # one video per frame folder
python plumetools/viz/video.py <viz-dir> --spec cases/cai2012/viz.yaml --framerate 10
python plumetools/viz/dissect.py sheet --matrix <matrix.yaml> --output sheet.png
python plumetools/viz/dissect.py converge <viz-dir> … --metric psnr
```

Without ffmpeg the frames are still written and the video step prints one line
saying no videos were made. The image metrics of `dissect converge` are
perceptual proxies for "has this stopped changing", never physical errors;
report them beside a physical metric.

## Reading what was drawn

The pitfalls that produce plausible wrong pictures are documented in
[Slice imagery](../viz-slices.md): `prefer_mean` falling back to instantaneous
fields where no average exists (use a sampling window), colour ranges that
differ between cases unless pinned, constant fields skipped, and cells holding
the invalid value (−1) drawn in the **bottom** colour of the bar — read
`mfpResolved` and `diagValid` before `KnCellMax`. `dsmcRhoN` is a parcel count,
not a density.
