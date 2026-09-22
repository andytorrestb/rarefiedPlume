r"""Configurable slice imagery from a dsmcFoam case, via VifPara and ParaView.

What this is for
----------------
``plumetools.cai2012.post`` extracts the four *validation* outputs and nothing
else -- that is deliberate, and §9 of the case specification says so. This
package is the other half: **survey imagery**. Any field the solver wrote, on any
cutting plane, as a PNG, driven by a YAML file rather than by edited code.

The two do not overlap and do not share code. ``post`` produces numbers that get
scored against Cai; this produces pictures that get looked at.

The layering, and why it is split this way
------------------------------------------
``paraview.simple`` only exists inside ParaView's own interpreter, so any module
that imports it can only be imported by ``pvpython``. That would put the whole
package out of reach of the test suite, which the repository runs with no
OpenFOAM and no ParaView installed.

So the split is strict:

======================  ====================  ==============================
:mod:`~.spec`           no ParaView           schema, YAML, validation
:mod:`~.catalog`        no ParaView           what the dsmcFoam fields *are*
:mod:`~.resolve`        no ParaView           which time, which field, what scale
:mod:`~.video`          no ParaView           frames -> video, through ffmpeg
:mod:`~.dissect`        no ParaView           contact sheets and convergence
                                              curves, through ffmpeg
:mod:`~.render`         **imports ParaView**  turns a spec into PNGs
``render_slices.py``    **imports ParaView**  the command-line entry point
======================  ====================  ==============================

:mod:`~.video` and :mod:`~.dissect` are on the ParaView-free side deliberately,
not incidentally: re-encoding a study at another frame rate, or re-cutting a
case matrix into a contact sheet, is then seconds of work rather than twenty
minutes of re-rendering, and both run on a machine with no ParaView at all.

Neither is imported here. They need ``ffmpeg`` at *run* time but not at import
time, and keeping them out of the package namespace means importing
:mod:`plumetools.viz` stays as cheap as it is for the schema modules::

    from plumetools.viz import dissect, video

Importing :mod:`plumetools.viz` gives you the first three and never touches
ParaView. :mod:`plumetools.viz.render` is imported lazily, by name, only by the
entry script. That is what lets ``tests/unit/test_viz.py`` run in the ordinary
suite -- and the split is drawn where it is on purpose: the decisions most
likely to be quietly wrong (a stale colour range, the wrong time, the
instantaneous field instead of the averaged one) are all on the testable side.

Running it
----------
Not with ``python``. VifPara re-launches the script inside ``pvpython`` with the
virtual environment injected, so the launcher has to be used::

    vifpara plumetools/viz/render_slices.py cases/cai2012/Cases/Kn100

See ``docs/viz-slices.md`` for the configuration file and
``docs/environment.md`` for the install.
"""

from plumetools.viz.catalog import (
    FIELD_CATALOG,
    FieldKind,
    base_name,
    catalog_entry,
    mean_name,
)
from plumetools.viz.resolve import (
    NothingToRender,
    RenderError,
    colour_range,
    component_range,
    decade_labels,
    expression_fields,
    is_constant,
    nudge_plane_origin,
    resolve_field_name,
    resolve_times,
    view_extents,
    view_half_height,
    view_shape,
)
from plumetools.viz.spec import (
    DEFAULT_SPEC_PATH,
    FieldSpec,
    ImageSpec,
    OutputSpec,
    PlaneSpec,
    RenderTask,
    SamplingSpec,
    VideoSpec,
    VizSpec,
    VizSpecError,
    load_spec,
)
from plumetools.viz.video import (
    EncodedVideo,
    Series,
    VideoError,
    concat_list,
    encode_all,
    encode_series,
    ffmpeg_command,
    have_ffmpeg,
    series_from_manifest,
    series_from_records,
)

__all__ = [
    "DEFAULT_SPEC_PATH",
    "FIELD_CATALOG",
    "FieldKind",
    "FieldSpec",
    "ImageSpec",
    "EncodedVideo",
    "NothingToRender",
    "OutputSpec",
    "PlaneSpec",
    "RenderError",
    "RenderTask",
    "SamplingSpec",
    "Series",
    "VideoError",
    "VideoSpec",
    "VizSpec",
    "VizSpecError",
    "base_name",
    "catalog_entry",
    "colour_range",
    "component_range",
    "concat_list",
    "decade_labels",
    "encode_all",
    "encode_series",
    "expression_fields",
    "ffmpeg_command",
    "have_ffmpeg",
    "is_constant",
    "load_spec",
    "mean_name",
    "nudge_plane_origin",
    "resolve_field_name",
    "resolve_times",
    "series_from_manifest",
    "series_from_records",
    "view_extents",
    "view_half_height",
    "view_shape",
]
