r"""What the fields ``dsmcFoam`` writes actually are.

This is a lookup table, not a policy. Its whole job is to stop every YAML file in
the repository from having to repeat that ``rhoN`` is a number density in
``m^-3`` that wants a logarithmic colour scale, or -- the one that matters --
that ``q`` is a **surface** field and cannot be sliced.

The volume/surface distinction is not cosmetic
----------------------------------------------
``dsmcFoam`` writes two kinds of ``volField``. Most carry a real internal field
sampled from the parcels in each cell. Four do not:

.. code-block:: text

    q          internalField uniform 0     heat flux INTO a boundary   [W/m^2]
    fD         internalField uniform (0 0 0)  force density ON a boundary [Pa]
    boundaryT  internalField uniform 0     wall temperature            [K]
    boundaryU  internalField uniform (0 0 0)  wall velocity            [m/s]

These live entirely in their ``boundaryField``; ``DSMCCloud`` accumulates them
when a molecule *strikes a wall*. Cutting a plane through the interior of the
domain and colouring it by ``q`` produces a uniformly zero image every time --
not a bug in the plot, an unsampled quantity. So they are tagged
:data:`FieldKind.SURFACE` and rendered on the boundary patches instead.

For a free plume with no solid surface -- ``cases/cai2012`` is exactly that: its
three patches are two vacuums and an inlet -- ``q`` and ``fD`` are identically
zero *everywhere*, patches included, and there is nothing to draw. The renderer
detects that as a constant field and says so rather than writing a blank PNG.
See ``skip_constant_fields`` in :class:`~plumetools.viz.spec.SamplingSpec`.

Averaged fields
---------------
``fieldAverage`` writes ``<name>Mean`` beside every field it is given. Those are
the ones worth looking at: a single-timestep ``rhoN`` is one step's worth of
parcels, because ``DSMCCloud::resetFields`` zeroes it every step, so an image of
it is a picture of shot noise. :func:`mean_name` applies the suffix and
``sampling.prefer_mean`` (on by default) makes the renderer reach for it
whenever the solver wrote one.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

#: Suffix ``fieldAverage`` appends to every field it averages.
MEAN_SUFFIX = "Mean"


class FieldKind(str, Enum):
    """Where a field carries data, and therefore how it can be drawn.

    ``VOLUME`` fields have a meaningful internal field and can be cut by a
    plane. ``SURFACE`` fields exist only on boundary patches and cannot -- see
    the module docstring.
    """

    VOLUME = "volume"
    SURFACE = "surface"


@dataclass(frozen=True)
class CatalogEntry:
    """Defaults for one known ``dsmcFoam`` field.

    Anything a YAML ``fields:`` entry does not say is taken from here; anything
    it does say wins. A field absent from the catalogue is not an error -- it
    renders as a linear-scale volume scalar under its own name, which is the
    right guess for a field this table has not been told about yet.

    Attributes:
        kind: volume or surface; a surface field is never sliced.
        units: printed in the colour-bar title. Not parsed.
        component: for vectors, the ParaView component to colour by --
            ``"Magnitude"``, ``"X"``, ``"Y"``, ``"Z"``. Empty for scalars.
        log: whether a logarithmic colour scale is the sensible default. True
            for densities, which span decades across a plume, and false for
            anything that changes sign or spans one decade.
        description: what the quantity is, for ``--list-fields``.
    """

    kind: FieldKind = FieldKind.VOLUME
    units: str = ""
    component: str = ""
    log: bool = False
    description: str = ""

    @property
    def is_vector(self) -> bool:
        """True when this field needs a component to colour by."""
        return bool(self.component)


#: Every field standard ``dsmcFoam`` (OpenFOAM v2512) writes, plus the cell
#: geometry ``postProcess`` adds. Keys are the base names; the ``Mean``
#: counterparts resolve to the same entry through :func:`catalog_entry`.
FIELD_CATALOG: dict[str, CatalogEntry] = {
    # --- the sampled flow field ------------------------------------------- #
    "rhoN": CatalogEntry(
        units="m^-3", log=True,
        description="number density n -- the plume itself"),
    "rhoM": CatalogEntry(
        units="kg/m^3", log=True,
        description="mass density rho = m n"),
    "momentum": CatalogEntry(
        units="kg/(m^2 s)", component="Magnitude", log=True,
        description="momentum density rho u; divide by rhoM for velocity"),
    "linearKE": CatalogEntry(
        units="J/m^3", log=True,
        description="translational kinetic energy density (1/2) rho <v^2>"),
    "internalE": CatalogEntry(
        units="J/m^3", log=True,
        description="internal energy density; identically 0 for a monatomic gas"),
    "iDof": CatalogEntry(
        units="m^-3", log=True,
        description="internal degrees of freedom density; 0 for a monatomic gas"),

    # --- DSMC statistics, not physics -------------------------------------- #
    "dsmcRhoN": CatalogEntry(
        units="m^-3", log=True,
        description="PARCEL number density -- times cell volume, the parcels "
                    "per cell. A statistical-quality field, not a physical one"),
    "dsmcSigmaTcRMax": CatalogEntry(
        units="m^3/s", log=True,
        description="(sigma_T c_r)_max, the NTC collision-selection ceiling"),

    # --- boundary-only: see the module docstring --------------------------- #
    "q": CatalogEntry(
        kind=FieldKind.SURFACE, units="W/m^2",
        description="heat flux into a boundary. Zero unless the case has a WALL"),
    "fD": CatalogEntry(
        kind=FieldKind.SURFACE, units="Pa", component="Magnitude",
        description="force density on a boundary. Zero unless the case has a WALL"),
    "boundaryT": CatalogEntry(
        kind=FieldKind.SURFACE, units="K",
        description="boundary temperature"),
    "boundaryU": CatalogEntry(
        kind=FieldKind.SURFACE, units="m/s", component="Magnitude",
        description="boundary velocity"),

    # --- written by postProcess, not by the solver -------------------------- #
    "V": CatalogEntry(
        units="m^3", log=True,
        description="cell volume (postProcess -func writeCellVolumes). The mesh "
                    "is graded, so this is worth looking at once"),
    "C": CatalogEntry(
        component="Magnitude",
        description="cell centres (postProcess -func writeCellCentres)"),
}


def mean_name(field: str) -> str:
    """The name ``fieldAverage`` writes for the running average of *field*.

    Applying the suffix twice is a no-op, so a spec that already names
    ``rhoNMean`` is left alone rather than asking for ``rhoNMeanMean``.

    Args:
        field: a base field name, e.g. ``"rhoN"``.

    Returns:
        e.g. ``"rhoNMean"``.
    """
    if field.endswith(MEAN_SUFFIX):
        return field
    return f"{field}{MEAN_SUFFIX}"


def base_name(field: str) -> str:
    """Strip a trailing ``Mean``, so ``rhoNMean`` looks up ``rhoN``."""
    if field.endswith(MEAN_SUFFIX) and len(field) > len(MEAN_SUFFIX):
        return field[: -len(MEAN_SUFFIX)]
    return field


def catalog_entry(field: str) -> CatalogEntry:
    """Catalogue defaults for *field*, averaged or not.

    An unknown field returns the default :class:`CatalogEntry` -- a linear-scale
    volume scalar -- rather than raising. A field this table has not heard of is
    a field somebody added to the solver, not a mistake in the YAML, and
    refusing to draw it would be the wrong call.

    Args:
        field: ``"rhoN"``, ``"rhoNMean"``, or anything else.

    Returns:
        The catalogue entry, or a permissive default.
    """
    name = base_name(field)
    return FIELD_CATALOG.get(name, CatalogEntry())
