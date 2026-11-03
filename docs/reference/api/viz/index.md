# `plumetools.viz`

Configurable slice imagery from a `dsmcFoam` case through ParaView, and the
ParaView-free steps that follow it: videos from the frames and a dissection of
a rendered sweep (contact sheets, image-convergence curves). The YAML spec says
what to draw; `resolve` decides how, with no ParaView in it; `scene` and
`render` drive `paraview.simple`.

Two runtimes: `render`, `scene` and the script `render_slices.py` run **only
under ParaView's `pvpython`**; `spec`, `resolve`, `catalog`, `video` and
`dissect` import under plain Python (the last two need ffmpeg to do anything).
The documentation build reads all of them statically, so it needs neither.

Depends on: `PyYAML`, `numpy`; `paraview.simple` (render, scene); ffmpeg
(video, dissect). Workflows: every family's imagery
([Post-processing](../../../guides/postprocessing.md)). Guide and caveats:
[Slice imagery](../../../viz-slices.md).

::: plumetools.viz
    options:
      show_root_heading: false
      members: false

## Modules

| Module | Purpose | Principal interfaces | Runtime |
|---|---|---|---|
| [`plumetools.viz.spec`](spec.md) | the sampling configuration: schema, loading, validation | `load_spec`, `VizSpec`, `FieldSpec`, `PlaneSpec`, `SamplingSpec`, `VizSpecError`, `DEFAULT_SPEC_PATH` | Python |
| [`plumetools.viz.catalog`](catalog.md) | what the fields `dsmcFoam` writes actually are | `catalog_entry`, `CatalogEntry`, `FieldKind`, `mean_name`, `base_name` | Python |
| [`plumetools.viz.resolve`](resolve.md) | times, field names, colour ranges, view geometry | `resolve_times`, `resolve_field_name`, `colour_range`, `is_constant`, `slice_frame`, `RenderError`, `NothingToRender` | Python |
| [`plumetools.viz.scene`](scene.md) | ParaView driven directly: reader, layout, views, one PNG | `Case`, `Layout`, `SliceView`, `SurfaceView`, `ColourMap`, `SessionLog` | pvpython |
| [`plumetools.viz.render`](render.md) | a `VizSpec` into PNGs | `Renderer`, `render_case`, `discover_case`, `RenderedImage`, `SkippedImage` | pvpython |
| [`plumetools.viz.render_slices`](render_slices.md) | the command-line entry point | `main`, `build_parser` | pvpython |
| [`plumetools.viz.video`](video.md) | encode frame folders with ffmpeg | `encode_all`, `encode_series`, `ffmpeg_command`, `series_from_manifest`, `have_ffmpeg`, `main` | Python + ffmpeg |
| [`plumetools.viz.dissect`](dissect.md) | contact sheets and image-difference curves | `build_contact_sheet`, `measure_series`, `convergence_curve`, `convergence_report`, `main` | Python + ffmpeg |

`plumetools/viz/slices.yaml` is the package-default spec (shipped as package
data). Its schema is in the [configuration reference](../../configuration.md#slice-imagery-specs).

<!-- docs-example: viz-spec -->
```python
from plumetools.viz.spec import load_spec

spec = load_spec()                         # plumetools/viz/slices.yaml, validated
planes = [plane.name for plane in spec.planes]
assert spec.sampling.prefer_mean            # draw *Mean fields when they exist
```

A spec only says what to draw. Whether the picture means anything depends on
the four traps [Slice imagery](../../../viz-slices.md) lists — among them that
`prefer_mean` silently falls back to the instantaneous field where no average
exists, and that cells holding the invalid value are drawn in the bottom colour
of the bar.
