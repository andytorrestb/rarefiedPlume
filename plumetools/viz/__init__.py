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
:mod:`~.render`         **imports ParaView**  turns a spec into PNGs
``render_slices.py``    **imports ParaView**  the command-line entry point
======================  ====================  ==============================

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
    RenderError,
    colour_range,
    component_range,
    decade_labels,
    expression_fields,
    is_constant,
    resolve_field_name,
    resolve_times,
)
from plumetools.viz.spec import (
    DEFAULT_SPEC_PATH,
    FieldSpec,
    ImageSpec,
    OutputSpec,
    PlaneSpec,
    RenderTask,
    SamplingSpec,
    VizSpec,
    VizSpecError,
    load_spec,
)

__all__ = [
    "DEFAULT_SPEC_PATH",
    "FIELD_CATALOG",
    "FieldKind",
    "FieldSpec",
    "ImageSpec",
    "OutputSpec",
    "PlaneSpec",
    "RenderError",
    "RenderTask",
    "SamplingSpec",
    "VizSpec",
    "VizSpecError",
    "base_name",
    "catalog_entry",
    "colour_range",
    "component_range",
    "decade_labels",
    "expression_fields",
    "is_constant",
    "load_spec",
    "mean_name",
    "resolve_field_name",
    "resolve_times",
]
