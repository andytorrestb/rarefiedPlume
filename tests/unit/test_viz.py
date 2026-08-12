"""Tier 0 for :mod:`plumetools.viz` -- the half that needs no ParaView.

Nothing here imports :mod:`plumetools.viz.render`. What is covered is the part
that decides *what* gets drawn and *how it is scaled*, which is where a wrong
answer is invisible: a mirrored plume, an instantaneous field passed off as an
averaged one, or a colour range that silently discards six decades all produce
an image that looks entirely plausible.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

from plumetools.viz import (
    DEFAULT_SPEC_PATH,
    FieldKind,
    VizSpecError,
    base_name,
    catalog_entry,
    colour_range,
    component_range,
    decade_labels,
    expression_fields,
    is_constant,
    load_spec,
    mean_name,
    resolve_field_name,
)
from plumetools.viz.resolve import (
    NothingToRender,
    RenderError,
    nudge_plane_origin,
    resolve_times,
    view_half_height,
)
from plumetools.viz.spec import DerivedSpec, FieldSpec, PlaneSpec, VizSpec

REPO = Path(__file__).resolve().parents[2]

#: The sampling settings each study's ``AllpostCases`` renders with. One per
#: family, in the family's own directory -- they are not shared.
STUDY_SPECS = {
    "cai2012": REPO / "cases" / "cai2012" / "viz.yaml",
    "markelov1999": REPO / "cases" / "markelov1999" / "viz.yaml",
}


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

def write_spec(tmp_path, document, name="spec.yaml"):
    path = tmp_path / name
    path.write_text(yaml.safe_dump(document), encoding="utf-8")
    return path


MINIMAL = {
    "planes": [{"name": "xz", "normal": [0, 1, 0], "camera_up": [0, 0, 1]}],
    "fields": [{"name": "rhoN"}],
}


# --------------------------------------------------------------------------- #
# the shipped default
# --------------------------------------------------------------------------- #

def test_the_shipped_default_spec_loads_and_validates():
    spec = load_spec()
    assert spec.source == DEFAULT_SPEC_PATH
    assert [plane.name for plane in spec.planes] == ["xz", "xy", "yz-5D"]
    assert spec.sampling.prefer_mean is True
    assert spec.sampling.field_type == "CELLS"


def test_the_default_planes_put_the_flow_to_the_right():
    """screen-right is ``normal x camera_up``; a sign slip mirrors the plume.

    This is the check that a mirrored image would otherwise pass silently:
    every default plane must put its in-plane axis coordinate to the right.
    """
    spec = load_spec()
    wanted = {"xz": (1.0, 0.0, 0.0),      # +X, the plume axis
              "xy": (1.0, 0.0, 0.0),      # +X again
              "yz-5D": (0.0, 1.0, 0.0)}   # +Y, looking upstream
    for plane in spec.planes:
        normal, up = plane.normal, plane.camera_up
        right = (normal[1] * up[2] - normal[2] * up[1],
                 normal[2] * up[0] - normal[0] * up[2],
                 normal[0] * up[1] - normal[1] * up[0])
        assert right == pytest.approx(wanted[plane.name]), (
            f"plane {plane.name} is mirrored: screen-right is {right}")


def test_the_default_spec_draws_q_as_a_surface_field():
    """q is boundary-only. Cutting a plane through it gives a blank image."""
    spec = load_spec()
    derived = spec.derived_by_name
    kinds = {field.name: field.resolve_kind(derived) for field in spec.fields}
    assert kinds["q"] is FieldKind.SURFACE
    assert kinds["fD"] is FieldKind.SURFACE
    assert kinds["rhoN"] is FieldKind.VOLUME
    assert kinds["dsmcRhoN"] is FieldKind.VOLUME


def test_the_default_spec_expands_to_one_task_per_image():
    spec = load_spec()
    tasks = spec.tasks()
    volume = [task for task in tasks
              if task.field.resolve_kind(spec.derived_by_name)
              is FieldKind.VOLUME]
    surface = [task for task in tasks if task not in volume]

    # seven enabled volume fields on three planes
    assert len(volume) == 7 * 3
    # a surface field is not cut by anything, so it gets one view, not three
    assert len(surface) == 2
    assert {task.field.name for task in surface} == {"q", "fD"}


# --------------------------------------------------------------------------- #
# the catalogue
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("field, expected", [
    ("rhoN", "rhoNMean"),
    ("q", "qMean"),
    ("rhoNMean", "rhoNMean"),     # idempotent: never rhoNMeanMean
])
def test_mean_name_applies_the_suffix_once(field, expected):
    assert mean_name(field) == expected


@pytest.mark.parametrize("field, expected", [
    ("rhoNMean", "rhoN"),
    ("rhoN", "rhoN"),
    ("Mean", "Mean"),             # not a suffix on an empty base
])
def test_base_name_strips_the_suffix(field, expected):
    assert base_name(field) == expected


def test_an_averaged_field_inherits_its_base_catalogue_entry():
    assert catalog_entry("rhoNMean").log is True
    assert catalog_entry("qMean").kind is FieldKind.SURFACE
    assert catalog_entry("momentumMean").component == "Magnitude"


def test_an_unknown_field_gets_permissive_defaults_rather_than_an_error():
    """A field the catalogue has not heard of is new, not wrong."""
    entry = catalog_entry("someNewSolverField")
    assert entry.kind is FieldKind.VOLUME
    assert entry.log is False
    assert entry.component == ""


# --------------------------------------------------------------------------- #
# choosing the field
# --------------------------------------------------------------------------- #

def test_prefer_mean_reaches_for_the_averaged_field_first():
    field = FieldSpec(name="rhoN")
    assert field.candidate_names(True, {}) == ("rhoNMean", "rhoN")


def test_prefer_mean_off_reaches_for_the_instantaneous_field_first():
    field = FieldSpec(name="rhoN")
    assert field.candidate_names(False, {}) == ("rhoN", "rhoNMean")


def test_naming_the_averaged_field_pins_it_regardless_of_prefer_mean():
    field = FieldSpec(name="rhoNMean")
    assert field.candidate_names(False, {}) == ("rhoNMean",)


def test_a_derived_field_has_no_averaged_counterpart():
    derived = {"U": DerivedSpec(name="U", expression="momentum/rhoM")}
    field = FieldSpec(name="U")
    assert field.candidate_names(True, derived) == ("U",)


def test_the_instantaneous_field_is_the_fallback_before_averaging_starts():
    """Before fieldAverage's timeStart there is no rhoNMean, only rhoN."""
    field = FieldSpec(name="rhoN")
    early = ["rhoN", "momentum"]
    assert resolve_field_name(field.candidate_names(True, {}), early) == "rhoN"


def test_a_field_the_case_never_wrote_resolves_to_nothing():
    assert resolve_field_name(("qMean", "q"), ["rhoN"]) is None


# --------------------------------------------------------------------------- #
# colour ranges
# --------------------------------------------------------------------------- #

def test_a_log_range_floors_an_exactly_zero_minimum_by_decades():
    """A DSMC far field is exactly 0; a log scale has no colour for it."""
    assert colour_range(0.0, 1e16, log=True, decades=8) == (1e8, 1e16)


def test_a_log_range_keeps_a_positive_minimum_that_is_inside_the_window():
    assert colour_range(1e12, 1e16, log=True, decades=8) == (1e12, 1e16)


def test_a_log_range_raises_a_minimum_that_is_too_many_decades_down():
    low, high = colour_range(1e2, 1e16, log=True, decades=6)
    assert (low, high) == (1e10, 1e16)


def test_a_linear_range_is_left_alone():
    assert colour_range(-5.0, 300.0, log=False, decades=6) == (-5.0, 300.0)


def test_a_constant_field_has_no_range_to_draw():
    assert colour_range(0.0, 0.0, log=False, decades=6) is None
    assert colour_range(1.93e-16, 1.93e-16, log=True, decades=6) is None


def test_a_wholly_non_positive_field_cannot_be_log_scaled():
    assert colour_range(-3.0, -1.0, log=True, decades=6) is None


def test_a_non_finite_range_is_refused():
    assert colour_range(0.0, math.inf, log=False, decades=6) is None
    assert colour_range(math.nan, 1.0, log=False, decades=6) is None


def test_is_constant_tolerates_floating_point_noise_but_not_real_variation():
    assert is_constant(1.0, 1.0 + 1e-15)
    assert not is_constant(1.0, 1.0001)


# --------------------------------------------------------------------------- #
# colour bar labels
# --------------------------------------------------------------------------- #

def test_decade_labels_thin_a_wide_range_to_fit_the_bar():
    labels = decade_labels(1e8, 1e16, limit=5)
    assert len(labels) <= 5
    assert all(low <= value <= high
               for value, low, high in ((v, 1e8, 1e16) for v in labels))


def test_decade_labels_are_powers_of_ten():
    for value in decade_labels(1e8, 1e16, limit=5):
        assert math.log10(value) == pytest.approx(round(math.log10(value)))


def test_decade_labels_give_up_on_a_range_spanning_no_whole_decade():
    assert decade_labels(1.1e5, 9.9e5) == []
    assert decade_labels(0.0, 1e6) == []
    assert decade_labels(1e6, 1e6) == []


# --------------------------------------------------------------------------- #
# components
# --------------------------------------------------------------------------- #

def test_component_range_reads_the_magnitude_for_a_scalar():
    assert component_range((0.0, 5.0), "") == (0.0, 5.0)


def test_component_range_picks_a_named_vector_component():
    flattened = (-1.0, 1.0, -2.0, 2.0, -3.0, 3.0, 0.0, 9.0)
    assert component_range(flattened, "X") == (-1.0, 1.0)
    assert component_range(flattened, "Y") == (-2.0, 2.0)
    assert component_range(flattened, "Z") == (-3.0, 3.0)
    assert component_range(flattened, "Magnitude") == (0.0, 9.0)


# --------------------------------------------------------------------------- #
# times
# --------------------------------------------------------------------------- #

TIMES = [0.002, 0.004, 0.006]


@pytest.mark.parametrize("wanted, expected", [
    ("latest", [0.006]),
    ("first", [0.002]),
    ("all", TIMES),
    (0.004, [0.004]),
])
def test_resolve_times_selects_what_was_asked_for(wanted, expected):
    assert resolve_times(TIMES, wanted) == expected


def test_an_unwritten_time_is_an_error_not_a_snap_to_the_nearest():
    with pytest.raises(RenderError, match="no time 0.005"):
        resolve_times(TIMES, 0.005)


def test_a_case_that_never_ran_says_so():
    with pytest.raises(RenderError, match="Allrun"):
        resolve_times([], "latest")


def test_a_case_that_never_ran_is_not_counted_as_a_failure():
    """A part-way study is normal, and both families' per-case Allpost tolerate
    it. Imagery must not be stricter than the validation it accompanies."""
    with pytest.raises(NothingToRender):
        resolve_times([], "latest")
    # ...but a genuinely wrong time still is an error
    with pytest.raises(RenderError) as caught:
        resolve_times([0.002], 0.005)
    assert not isinstance(caught.value, NothingToRender)


# --------------------------------------------------------------------------- #
# derived expressions
# --------------------------------------------------------------------------- #

def test_mean_expands_to_follow_prefer_mean():
    entry = DerivedSpec(name="U", expression="momentum{mean}/rhoM{mean}")
    assert entry.resolved_expression(True) == "momentumMean/rhoMMean"
    assert entry.resolved_expression(False) == "momentum/rhoM"


def test_the_gas_constant_is_substituted_from_the_case():
    entry = DerivedSpec(name="Ttra", expression="linearKE{mean}/(3*{R})")
    resolved = entry.resolved_expression(False, {"R": "208.24"})
    assert resolved == "linearKE/(3*208.24)"


def test_an_unfillable_placeholder_is_an_error_not_a_guess():
    """A temperature computed with the wrong R is wrong by an invisible factor."""
    entry = DerivedSpec(name="Ttra", expression="linearKE/(3*{R})")
    with pytest.raises(VizSpecError, match="R"):
        entry.resolved_expression(False, {})


def test_expression_fields_finds_the_arrays_and_ignores_the_functions():
    fields = expression_fields(
        "(2*linearKEMean/rhoMMean - mag(momentumMean/rhoMMean)^2)/(3*208.2)")
    assert fields == {"linearKEMean", "rhoMMean", "momentumMean"}


# --------------------------------------------------------------------------- #
# schema validation
# --------------------------------------------------------------------------- #

def test_an_unknown_key_is_rejected(tmp_path):
    """A silently ignored 'prefere_mean' would draw shot noise for months."""
    document = dict(MINIMAL, sampling={"prefere_mean": True})
    with pytest.raises(VizSpecError, match="prefere_mean"):
        load_spec(write_spec(tmp_path, document))


def test_an_unknown_field_key_is_rejected(tmp_path):
    document = dict(MINIMAL, fields=[{"name": "rhoN", "logarithmic": True}])
    with pytest.raises(VizSpecError, match="logarithmic"):
        load_spec(write_spec(tmp_path, document))


def test_a_bare_string_is_a_valid_field_entry(tmp_path):
    document = dict(MINIMAL, fields=["rhoN", "dsmcRhoN"])
    spec = load_spec(write_spec(tmp_path, document))
    assert [field.name for field in spec.fields] == ["rhoN", "dsmcRhoN"]


def test_duplicate_plane_names_are_rejected(tmp_path):
    document = dict(MINIMAL, planes=[{"name": "xz"}, {"name": "xz"}])
    with pytest.raises(VizSpecError, match="duplicate"):
        load_spec(write_spec(tmp_path, document))


def test_duplicate_field_names_are_rejected(tmp_path):
    document = dict(MINIMAL, fields=[{"name": "rhoN"}, {"name": "rhoN"}])
    with pytest.raises(VizSpecError, match="duplicate"):
        load_spec(write_spec(tmp_path, document))


def test_a_camera_up_parallel_to_the_normal_has_no_horizontal_axis(tmp_path):
    document = dict(MINIMAL, planes=[
        {"name": "bad", "normal": [0, 0, 1], "camera_up": [0, 0, 1]}])
    with pytest.raises(VizSpecError, match="parallel"):
        load_spec(write_spec(tmp_path, document))


def test_a_zero_normal_is_rejected(tmp_path):
    document = dict(MINIMAL, planes=[{"name": "bad", "normal": [0, 0, 0]}])
    with pytest.raises(VizSpecError, match="zero vector"):
        load_spec(write_spec(tmp_path, document))


def test_origin_and_origin_over_L_are_mutually_exclusive(tmp_path):
    document = dict(MINIMAL, planes=[
        {"name": "xz", "origin": [0, 0, 0], "origin_over_L": [5, 0, 0]}])
    with pytest.raises(VizSpecError, match="not both"):
        load_spec(write_spec(tmp_path, document))


def test_an_inverted_range_is_rejected(tmp_path):
    document = dict(MINIMAL, fields=[{"name": "rhoN", "range": [1e16, 1e8]}])
    with pytest.raises(VizSpecError, match="must exceed"):
        load_spec(write_spec(tmp_path, document))


def test_an_unknown_time_selector_is_rejected(tmp_path):
    document = dict(MINIMAL, sampling={"time": "lastest"})
    with pytest.raises(VizSpecError, match="lastest"):
        load_spec(write_spec(tmp_path, document))


def test_an_unknown_filename_placeholder_is_rejected(tmp_path):
    document = dict(MINIMAL, output={"filename": "{fieldname}_{plane}"})
    with pytest.raises(VizSpecError, match="placeholder"):
        load_spec(write_spec(tmp_path, document))


def test_fields_with_no_planes_would_draw_nothing(tmp_path):
    document = {"planes": [], "fields": [{"name": "rhoN"}]}
    with pytest.raises(VizSpecError, match="nothing would be drawn"):
        load_spec(write_spec(tmp_path, document))


def test_a_future_schema_version_is_refused(tmp_path):
    document = dict(MINIMAL, version=2)
    with pytest.raises(VizSpecError, match="version"):
        load_spec(write_spec(tmp_path, document))


def test_a_field_naming_an_undefined_plane_is_rejected(tmp_path):
    document = dict(MINIMAL, fields=[{"name": "rhoN", "planes": ["nope"]}])
    spec = load_spec(write_spec(tmp_path, document))
    with pytest.raises(VizSpecError, match="nope"):
        spec.tasks()


# --------------------------------------------------------------------------- #
# origins in diameters
# --------------------------------------------------------------------------- #

def test_origin_over_L_is_scaled_by_the_length_scale():
    plane = PlaneSpec(name="yz", origin_over_L=(5.0, 0.0, 0.0))
    assert plane.origin_m(0.2) == (1.0, 0.0, 0.0)


def test_origin_over_L_without_a_length_scale_is_an_error_not_a_guess():
    """Defaulting L to 1 m would put the plane somewhere arbitrary, and still
    draw a picture that looked fine."""
    plane = PlaneSpec(name="yz", origin_over_L=(5.0, 0.0, 0.0))
    with pytest.raises(VizSpecError, match="length scale"):
        plane.origin_m(None)


def test_an_origin_in_metres_ignores_the_length_scale():
    plane = PlaneSpec(name="xz", origin=(0.0, 0.0, 0.0))
    assert plane.origin_m(None) == (0.0, 0.0, 0.0)
    assert plane.origin_m(0.2) == (0.0, 0.0, 0.0)


# --------------------------------------------------------------------------- #
# inheritance
# --------------------------------------------------------------------------- #

def test_extends_default_inherits_the_shipped_planes(tmp_path):
    path = write_spec(tmp_path, {"extends": "default",
                                 "fields": [{"name": "rhoN"}]})
    spec = load_spec(path)
    assert [plane.name for plane in spec.planes] == ["xz", "xy", "yz-5D"]
    assert [field.name for field in spec.fields] == ["rhoN"]


def test_extends_merges_mappings_key_by_key(tmp_path):
    path = write_spec(tmp_path, {"extends": "default",
                                 "sampling": {"prefer_mean": False}})
    spec = load_spec(path)
    assert spec.sampling.prefer_mean is False
    # untouched keys survive the merge
    assert spec.sampling.field_type == "CELLS"
    assert spec.sampling.skip_constant_fields is True


def test_extends_replaces_lists_wholesale(tmp_path):
    """Three listed fields must draw three fields, not three plus the parent's."""
    path = write_spec(tmp_path, {"extends": "default",
                                 "fields": [{"name": "rhoN"}]})
    spec = load_spec(path)
    assert len(spec.fields) == 1


def test_extends_follows_a_relative_path(tmp_path):
    write_spec(tmp_path, MINIMAL, name="parent.yaml")
    child = write_spec(tmp_path, {"extends": "parent.yaml",
                                  "image": {"height": 1200}}, name="child.yaml")
    spec = load_spec(child)
    assert spec.image.height == 1200
    assert [field.name for field in spec.fields] == ["rhoN"]


def test_a_circular_extends_is_caught(tmp_path):
    write_spec(tmp_path, {"extends": "b.yaml"}, name="a.yaml")
    write_spec(tmp_path, {"extends": "a.yaml"}, name="b.yaml")
    with pytest.raises(VizSpecError, match="circular"):
        load_spec(tmp_path / "a.yaml")


def test_a_missing_spec_file_says_so(tmp_path):
    with pytest.raises(VizSpecError, match="no spec file"):
        load_spec(tmp_path / "absent.yaml")


# --------------------------------------------------------------------------- #
# task filtering
# --------------------------------------------------------------------------- #

def test_a_disabled_field_draws_nothing(tmp_path):
    document = dict(MINIMAL, fields=[{"name": "rhoN", "enabled": False},
                                     {"name": "dsmcRhoN"}])
    spec = load_spec(write_spec(tmp_path, document))
    assert [task.field.name for task in spec.tasks()] == ["dsmcRhoN"]


def test_filtering_by_field_and_plane(tmp_path):
    document = {
        "planes": [{"name": "xz"}, {"name": "xy", "normal": [0, 0, -1],
                                    "camera_up": [0, 1, 0]}],
        "fields": [{"name": "rhoN"}, {"name": "dsmcRhoN"}],
    }
    spec = load_spec(write_spec(tmp_path, document))
    assert len(spec.tasks()) == 4
    tasks = spec.tasks(only_fields=["rhoN"], only_planes=["xz"])
    assert [(task.field.name, task.plane_name) for task in tasks] == [("rhoN", "xz")]


def test_filtering_on_an_unknown_name_is_an_error(tmp_path):
    spec = load_spec(write_spec(tmp_path, MINIMAL))
    with pytest.raises(VizSpecError, match="no such field"):
        spec.tasks(only_fields=["rohN"])
    with pytest.raises(VizSpecError, match="no such plane"):
        spec.tasks(only_planes=["zx"])


def test_a_field_can_restrict_itself_to_one_plane(tmp_path):
    document = {
        "planes": [{"name": "xz"}, {"name": "xy", "normal": [0, 0, -1],
                                    "camera_up": [0, 1, 0]}],
        "fields": [{"name": "rhoN", "planes": ["xz"]}],
    }
    spec = load_spec(write_spec(tmp_path, document))
    assert [task.plane_name for task in spec.tasks()] == ["xz"]


def test_a_surface_field_can_be_viewed_from_several_planes(tmp_path):
    document = {
        "planes": [{"name": "xz"}, {"name": "xy", "normal": [0, 0, -1],
                                    "camera_up": [0, 1, 0]}],
        "fields": [{"name": "q", "planes": ["xz", "xy"]}],
    }
    spec = load_spec(write_spec(tmp_path, document))
    assert [task.plane_name for task in spec.tasks()] == ["xz", "xy"]


# --------------------------------------------------------------------------- #
# the shipped study specs
#
# One per family, in the family's own directory. These tests are what stops a
# spec that ./AllpostCases renders unattended from silently rotting -- an
# unknown key or a flipped normal would otherwise surface as a wrong picture
# on somebody's next study run.
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("study", sorted(STUDY_SPECS))
def test_a_study_ships_its_own_sampling_settings(study):
    assert STUDY_SPECS[study].is_file(), (
        f"cases/{study}/viz.yaml is missing; its AllpostCases defaults to it")


@pytest.mark.parametrize("study", sorted(STUDY_SPECS))
def test_a_study_spec_loads_and_validates(study):
    """Unknown keys are rejected, so this catches a typo in a shipped file."""
    spec = load_spec(STUDY_SPECS[study])
    assert spec.fields and spec.planes


@pytest.mark.parametrize("study", sorted(STUDY_SPECS))
def test_a_study_spec_draws_exactly_one_plane(study):
    """./AllpostCases renders these unattended; three planes would be an hour."""
    spec = load_spec(STUDY_SPECS[study])
    assert [plane.name for plane in spec.planes] == ["midplane"]


@pytest.mark.parametrize("study", sorted(STUDY_SPECS))
def test_a_study_cut_contains_the_plume_axis(study):
    """Both families put the plume on +x with their centre plane at y = 0."""
    plane = load_spec(STUDY_SPECS[study]).planes[0]
    assert plane.normal == (0.0, 1.0, 0.0)
    assert plane.origin == (0.0, 0.0, 0.0)
    # screen-right = normal x camera_up must be +x, or the plume is mirrored
    normal, up = plane.normal, plane.camera_up
    right = (normal[1] * up[2] - normal[2] * up[1],
             normal[2] * up[0] - normal[0] * up[2],
             normal[0] * up[1] - normal[1] * up[0])
    assert right == pytest.approx((1.0, 0.0, 0.0))


@pytest.mark.parametrize("study", sorted(STUDY_SPECS))
def test_a_study_spec_draws_one_image_per_field(study):
    spec = load_spec(STUDY_SPECS[study])
    tasks = spec.tasks()
    assert len(tasks) == len([f for f in spec.fields if f.enabled])
    assert {task.plane_name for task in tasks} == {"midplane"}


@pytest.mark.parametrize("study", sorted(STUDY_SPECS))
def test_a_study_spec_draws_the_plume_and_the_particle_statistics(study):
    """Whatever else differs, every study wants the plume and the parcel count
    that says whether the rest of the images are noise."""
    names = {field.name for field in load_spec(STUDY_SPECS[study]).fields}
    assert {"rhoN", "dsmcRhoN", "U", "Ttra"} <= names


def test_only_the_study_with_walls_asks_for_the_surface_fields():
    """This is the whole reason the two studies stopped sharing a spec.

    q and fD are accumulated where a molecule strikes a WALL. markelov1999 has
    a cylinder and a plate; cai2012's three patches are two vacuums and an
    inlet, so asking for them there would render nothing on every case of every
    run.
    """
    markelov = {f.name for f in load_spec(STUDY_SPECS["markelov1999"]).fields}
    cai = {f.name for f in load_spec(STUDY_SPECS["cai2012"]).fields}
    assert {"q", "fD"} <= markelov
    assert not ({"q", "fD"} & cai)


def test_the_surface_fields_are_classified_as_surfaces():
    spec = load_spec(STUDY_SPECS["markelov1999"])
    derived = spec.derived_by_name
    kinds = {field.name: field.resolve_kind(derived) for field in spec.fields}
    assert kinds["q"] is FieldKind.SURFACE
    assert kinds["fD"] is FieldKind.SURFACE
    assert kinds["rhoN"] is FieldKind.VOLUME


@pytest.mark.parametrize("study", sorted(STUDY_SPECS))
def test_a_study_spec_inherits_the_derived_fields(study):
    """U and Ttra come from `extends: default`, not from a copy in each file."""
    derived = load_spec(STUDY_SPECS[study]).derived_by_name
    assert "U" in derived and "Ttra" in derived
    assert "{R}" in derived["Ttra"].expression      # gas taken from the case


# --------------------------------------------------------------------------- #
# planes on a domain boundary
# --------------------------------------------------------------------------- #

#: A half domain: y runs from the symmetry plane at 0 up to 0.3.
HALF_DOMAIN = (0.0, 0.9, 0.0, 0.3, -0.35, 0.35)


def test_a_plane_on_the_symmetry_face_is_moved_inside():
    """cases/markelov1999 models y >= 0, so its centre plane is the mesh edge
    and VTK cuts nothing there."""
    origin, moved = nudge_plane_origin((0.0, 0.0, 0.0), (0.0, 1.0, 0.0),
                                       HALF_DOMAIN)
    assert moved
    assert 0.0 < origin[1] < 0.3
    assert origin[0] == 0.0 and origin[2] == 0.0


def test_the_nudge_stays_far_inside_the_first_cell():
    """Small enough that the picture is still the symmetry plane."""
    origin, _ = nudge_plane_origin((0.0, 0.0, 0.0), (0.0, 1.0, 0.0), HALF_DOMAIN)
    assert origin[1] == pytest.approx(3e-5)      # 1e-4 of the 0.3 m extent


def test_a_plane_at_the_far_face_is_moved_inside_too():
    origin, moved = nudge_plane_origin((0.0, 0.3, 0.0), (0.0, 1.0, 0.0),
                                       HALF_DOMAIN)
    assert moved
    assert origin[1] == pytest.approx(0.3 - 3e-5)


def test_a_plane_through_the_middle_is_left_alone():
    """cases/cai2012 is a full box, so y = 0 is genuinely interior."""
    full_box = (0.0, 2.0, -2.0, 2.0, -2.0, 2.0)
    origin, moved = nudge_plane_origin((0.0, 0.0, 0.0), (0.0, 1.0, 0.0), full_box)
    assert not moved
    assert origin == (0.0, 0.0, 0.0)


def test_a_tilted_plane_is_not_guessed_at():
    """Several directions would be 'inward' and none is obviously right."""
    origin, moved = nudge_plane_origin((0.0, 0.0, 0.0), (0.0, 1.0, 1.0),
                                       HALF_DOMAIN)
    assert not moved
    assert origin == (0.0, 0.0, 0.0)


def test_a_degenerate_extent_is_left_alone():
    flat = (0.0, 0.9, 0.0, 0.0, -0.35, 0.35)
    _, moved = nudge_plane_origin((0.0, 0.0, 0.0), (0.0, 1.0, 0.0), flat)
    assert not moved


# --------------------------------------------------------------------------- #
# framing a 3-D view
# --------------------------------------------------------------------------- #

def test_the_view_fits_a_tall_body_without_clipping():
    """The cylinder is 0.457 m along z and the camera looks along y."""
    bounds = (0.22, 0.53, 0.0, 0.0762, -0.2286, 0.2286)
    half = view_half_height(bounds, (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
                            aspect=4 / 3)
    assert half >= 0.2286                    # the half-extent along up
    assert half == pytest.approx(0.2286 * 1.06)


def test_a_wide_body_is_fitted_by_width_not_height():
    """Half the largest Cartesian extent would clip this one."""
    bounds = (0.0, 4.0, 0.0, 0.1, -0.2, 0.2)
    half = view_half_height(bounds, (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
                            aspect=1.0)
    assert half == pytest.approx(2.0 * 1.06)     # driven by the 4 m x-extent


def test_the_view_half_height_is_always_positive():
    """A degenerate body must not produce a zero parallel scale."""
    degenerate = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    assert view_half_height(degenerate, (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
                            aspect=4 / 3) > 0.0


# --------------------------------------------------------------------------- #
# legend titles
# --------------------------------------------------------------------------- #

def test_a_legend_title_carries_the_catalogue_units():
    field = FieldSpec(name="rhoN")
    assert field.resolve_legend_title("rhoNMean", {}) == "rhoNMean  [m^-3]"


def test_a_derived_field_uses_its_own_units():
    derived = {"U": DerivedSpec(name="U", expression="x", units="m/s")}
    assert FieldSpec(name="U").resolve_legend_title("U", derived) == "U  [m/s]"


def test_an_explicit_legend_title_wins():
    field = FieldSpec(name="rhoN", legend_title="n / n0")
    assert field.resolve_legend_title("rhoNMean", {}) == "n / n0"


def test_a_field_with_no_known_units_is_titled_by_name_alone():
    assert FieldSpec(name="mystery").resolve_legend_title("mystery", {}) == "mystery"
