r"""The numerics: mean free path, gradients, validity, and the diagnostic registry.

Pure numpy. Nothing here opens a file, knows what a ``case.yaml`` is, or draws
anything -- :mod:`~plumetools.rarefaction.foam` reads and writes OpenFOAM,
:mod:`~plumetools.rarefaction.summary` condenses, and the case families decide
where the numbers come from.

Inputs
------
Every diagnostic is built from a handful of named **inputs**, one cell value
each, that the I/O layer reads from a time directory:

==========================  ============================  =====================
input                       from (default)                meaning
==========================  ============================  =====================
``number_density``          ``rhoNMean``                  ``n`` [1/m^3]
``mass_density``            ``rhoMMean``                  ``rho`` [kg/m^3]
``momentum_density``        ``momentumMean``              ``rho u`` [kg/(m^2 s)]
``kinetic_energy_density``  ``linearKEMean``              ``rho <c^2>/2`` [J/m^3]
``parcels_per_step``        ``dsmcRhoNMean``              parcels in the cell,
                                                          averaged over steps
``averaged_steps``          ``uniform/.../functionObject  steps in the window
                            Properties``
``cell_volume``             ``V``                         [m^3]
``cell_extent``             ``constant/polyMesh``         bounding box [m]
``cell_graph``              ``polyMesh`` + ``C``          face neighbours
==========================  ============================  =====================

``dsmcRhoN`` is a **count**. ``DSMCCloud::calculateFields`` adds one per parcel
and never divides by the cell volume, while ``rhoN`` is multiplied by
``nParticle / V`` -- yet both files declare ``[0 -3 0 0 0 0 0]``. So
``rhoNMean * V / dsmcRhoNMean`` is the particle weight in every occupied cell,
which :mod:`~plumetools.rarefaction.foam` checks on every run.

Temperature
-----------
The translational temperature is taken from the averaged **moments**, never
averaged itself::

    3 R T = 2 <(1/2) rho c^2> / <rho>  -  |<rho u> / <rho>|^2

which is the cumulative-moment estimate a DSMC sample is meant to be read with.
In a sparsely sampled cell the difference can come out at or below zero, or
at the round-off of the moments -- one parcel, sampled for many steps, has no
thermal spread at all. Either way the cell has no temperature and is flagged;
nothing is clamped.

Mean free path
--------------
:func:`mean_free_path` is equation (G1) of :mod:`plumetools.nozzleplume.gas`,
element by element:

.. math::

    \lambda = \frac{1}{\sqrt2\,\pi d_{ref}^2\, n}
              \left(\frac{T}{T_{ref}}\right)^{\omega - 1/2}.

At ``omega = 0.5`` -- a hard sphere, the R-4D case -- the correction is exactly
one and ``lambda = 1/(sqrt(2) pi d^2 n)``: ``T`` does not enter. For VHS it is
the **equilibrium** mean free path at the local translational temperature. It
assumes a single Maxwellian at that temperature: in an expanding plume
(where the parallel and perpendicular temperatures part company), in a shock
layer, and wherever the jet meets the gas re-emitted from a wall (two streams,
whose collision rate is set by their relative speed and not by ``T``), it is a
nominal local length, not the distance a molecule actually flies.

Gradients
---------
:func:`least_squares_gradient` is a weighted least-squares fit over each cell's
face neighbours, weights ``1/|d|^2`` -- OpenFOAM's ``leastSquares``. It is
exact for a linear field on any stencil, including the 2:1 transitions of an
octree mesh, where a Gauss gradient is not. It is taken of the **averaged**
field: the gradient of a mean, never a mean of gradients.

* A valid cell's fit uses only VALID neighbours: beside an empty or
  under-sampled one it is one-sided, and with fewer than three non-coplanar
  valid neighbours it is not trusted, and its mask says so. Every other
  computable cell gets a second fit, from its computable neighbours, that is
  written for display and never trusted or summarised.
* Symmetry and empty patches contribute a mirror image of the cell, which
  imposes the zero normal gradient the plane implies. Other boundaries are
  one-sided.
* A gradient amplifies sampling noise: a relative noise ``e`` in ``n`` puts
  roughly ``e * lambda / Delta`` into ``Kn_grad`` -- negligible where cells are
  many mean free paths wide, and dominant where they are smaller than one.
* At a level transition the stencil is lopsided, and a field that is not
  linear across it -- a plume's density is not -- picks up a first-order error.

Validity
--------
Each cell gets a status, a sum of these flags (0 is valid):

=====  ===================  ===================================================
flag   name                 meaning
=====  ===================  ===================================================
1      ``NO_PARCELS``       no parcel was ever sampled there
2      ``FEW_PARCELS``      sampled, but not MORE than ``min_sampled_particles``
4      ``BAD_DENSITY``      ``n <= 0`` or not finite
8      ``BAD_TEMPERATURE``  no resolvable ``T > 0`` (see Temperature above)
=====  ===================  ===================================================

Two kinds of invalid cell are handled differently:

* **uncomputable** -- flags 1, 4 or 8 (:data:`UNCOMPUTABLE`): there is nothing
  to compute from. NaN in memory, ``invalid_value`` on disk.
* **low-confidence** -- flag 2 alone: the cell was sampled and has n and T,
  just from too few parcels. Its values are computed and written, so the
  fields show the whole sampled flow; ``diagValid`` is 0 there.

Either way the cell is not ``valid``, and every summary leaves it out. Nothing
is ever replaced by the exit state.

Adding a diagnostic
-------------------
One decorated function. It names what it needs, and the registry works out
whether a case has it::

    @register("KnGradP", requires=("mean_free_path", "pressure", "grad_pressure"),
              dimensions=DIMENSIONLESS, units="-",
              description="lambda |grad p| / p")
    def _kn_grad_p(ctx):
        lam = ctx["mean_free_path"]
        ...
        return Field(values, valid)

A requirement that is neither an input nor a registered quantity -- or an input
the I/O layer could not read -- makes that diagnostic, and only that one, be
skipped with the reason. See :func:`evaluate`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

#: OpenFOAM dimension sets, ``[kg m s K mol A cd]``.
DIMENSIONLESS = "[0 0 0 0 0 0 0]"
LENGTH = "[0 1 0 0 0 0 0]"

NO_PARCELS = 1
FEW_PARCELS = 2
BAD_DENSITY = 4
BAD_TEMPERATURE = 8

#: Flags that leave a cell with nothing to compute from. FEW_PARCELS is not
#: one: such a cell has a value, written and flagged, just not trusted.
UNCOMPUTABLE = NO_PARCELS | BAD_DENSITY | BAD_TEMPERATURE

#: Status flag -> what it means, for reports.
STATUS_FLAGS = {
    NO_PARCELS: "no_parcels_sampled",
    FEW_PARCELS: "too_few_sampled_parcels",
    BAD_DENSITY: "nonpositive_number_density",
    BAD_TEMPERATURE: "no_resolved_temperature",
}

#: How far a stencil's moment matrix may be from singular, relative to its
#: own scale, before the gradient is refused. Three orthogonal neighbours give
#: 1; three coplanar ones give 0 to round-off.
GRADIENT_CONDITION = 1.0e-6

#: The smallest thermal part of ``<c^2>`` that is a temperature rather than
#: round-off. dsmcFoam writes 10 significant figures, so the moments carry
#: ~1e-10 relative error; 1e-8 leaves a factor of 100. No flow comes near it:
#: at 3000 m/s it is 7e-5 K. A cell sampled by one parcel sits at ~1e-12.
TEMPERATURE_RESOLUTION = 1.0e-8


@dataclass(frozen=True)
class Field:
    """One cell value per cell, and whether each can be trusted.

    ``values`` is NaN only where nothing could be computed. A value in a cell
    that is not ``valid`` -- computable, but from too few samples -- is kept, to
    be written and looked at; every statistic reads ``valid`` and leaves it out.
    """

    values: np.ndarray
    valid: np.ndarray

    @classmethod
    def masked(cls, values, valid) -> "Field":
        """Values only where ``valid``: NaN everywhere else."""
        values = np.asarray(values, dtype=np.float64)
        valid = np.asarray(valid, dtype=bool) & np.isfinite(values)
        return cls(np.where(valid, values, np.nan), valid)

    @classmethod
    def flagged(cls, values, valid) -> "Field":
        """Values wherever they are finite; ``valid`` says which to trust."""
        values = np.asarray(values, dtype=np.float64)
        return cls(values, np.asarray(valid, dtype=bool) & np.isfinite(values))


@dataclass(frozen=True)
class Gradient:
    """``(n_cells, 3)`` gradient vectors, NaN where ``valid`` is false."""

    vectors: np.ndarray
    valid: np.ndarray

    @property
    def magnitude(self) -> np.ndarray:
        return np.linalg.norm(self.vectors, axis=1)


@dataclass(frozen=True)
class CellGraph:
    """Which cells share a face, and where the cells are.

    Attributes:
        centres: ``(n_cells, 3)`` cell centroids [m].
        owner, neighbour: the two cells of each INTERNAL face.
        mirror_cells: cells with a face on a symmetry or empty patch, once per
            such face.
        mirror_offsets: ``(k, 3)``, from each of those cells' centres to its
            mirror image in the patch [m].
    """

    centres: np.ndarray
    owner: np.ndarray
    neighbour: np.ndarray
    mirror_cells: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int64))
    mirror_offsets: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))

    @property
    def n_cells(self) -> int:
        return int(self.centres.shape[0])


# --------------------------------------------------------------------------- #
# physics
# --------------------------------------------------------------------------- #

def mean_free_path(number_density, temperature, species, *,
                   convention: str = "vhs") -> np.ndarray:
    """Equation (G1), element by element [m]. NaN where it is undefined.

    Args:
        number_density: ``n`` [1/m^3].
        temperature: translational ``T`` [K]. Unused for a hard sphere.
        species: anything with ``diameter_m``, ``omega`` and ``t_ref_K`` --
            :class:`plumetools.nozzleplume.gas.VhsSpecies`, for one.
        convention: ``vhs`` (temperature-corrected) or ``hard_sphere``, as in
            :data:`plumetools.nozzleplume.gas.MEAN_FREE_PATH_CONVENTIONS`.

    ``n <= 0`` gives NaN, not the infinity the scalar version returns: here it
    is a cell with no data, and a statistic must not count it.
    """
    if convention not in ("vhs", "hard_sphere"):
        raise ValueError(f"unknown mean-free-path convention {convention!r}")
    n = np.asarray(number_density, dtype=np.float64)
    sigma = math.pi * float(species.diameter_m) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.where((n > 0.0) & np.isfinite(n),
                       1.0 / (math.sqrt(2.0) * sigma * n), np.nan)
        exponent = float(species.omega) - 0.5
        if convention == "vhs" and exponent != 0.0:
            t = np.asarray(temperature, dtype=np.float64)
            lam = np.where((t > 0.0) & np.isfinite(t),
                           lam * (t / float(species.t_ref_K)) ** exponent, np.nan)
    return lam


def translational_temperature(mass_density, momentum_density,
                              kinetic_energy_density,
                              gas_constant_j_per_kg_k: float, *,
                              resolution: float = TEMPERATURE_RESOLUTION
                              ) -> np.ndarray:
    """``T`` [K] from the averaged moments; NaN where there is none.

    ``3 R T`` is the difference of ``<c^2>`` and ``|u|^2``, and where the
    first is no more than ``1 + resolution`` times the second the difference
    is round-off in the written moments, not a temperature. One parcel sampled
    many times gives exactly that -- a single molecule has no thermal spread --
    and a lingering one can pass a sample-count threshold on its own. NaN
    there, as for ``rho <= 0`` or a negative difference; the status flags it.
    Never clamped.
    """
    rho = np.asarray(mass_density, dtype=np.float64)
    momentum = np.asarray(momentum_density, dtype=np.float64).reshape(-1, 3)
    kinetic = np.asarray(kinetic_energy_density, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        total = 2.0 * kinetic / rho                     # <c^2>
        thermal = total - (momentum ** 2).sum(axis=1) / rho ** 2
        t = thermal / (3.0 * float(gas_constant_j_per_kg_k))
        resolved = (rho > 0.0) & (thermal > float(resolution) * total)
    return np.where(resolved, t, np.nan)


def least_squares_gradient(values, valid, graph: CellGraph) -> Gradient:
    """Weighted least-squares gradient over valid face neighbours.

    For cell ``P`` with valid neighbours ``j`` at ``d_j = x_j - x_P``::

        (sum w d d^T) grad = sum w d (phi_j - phi_P),     w = 1 / |d|^2

    Accumulated with ``bincount`` -- six moment components and three
    right-hand sides, no per-face matrices -- and solved by cofactors, so the
    memory is a few arrays the size of the mesh. See the module docstring for
    what is and is not in the stencil.
    """
    phi = np.asarray(values, dtype=np.float64)
    ok = np.asarray(valid, dtype=bool) & np.isfinite(phi)
    n = graph.n_cells

    owner, neighbour = graph.owner, graph.neighbour
    use = ok[owner] & ok[neighbour]
    p, q = owner[use], neighbour[use]
    d = graph.centres[q] - graph.centres[p]
    w = 1.0 / np.einsum("ij,ij->i", d, d)
    dphi = phi[q] - phi[p]

    # Both cells of a face get the same w d d^T, and the same w d dphi: the
    # neighbour's view is (-d)(-dphi).
    cells = np.concatenate([p, q])
    mirror = graph.mirror_cells
    keep = ok[mirror] if mirror.size else np.zeros(0, dtype=bool)
    mirror, offsets = mirror[keep], graph.mirror_offsets[keep]
    w_m = 1.0 / np.einsum("ij,ij->i", offsets, offsets) if mirror.size else np.empty(0)

    def moment(i: int, j: int) -> np.ndarray:
        face = w * d[:, i] * d[:, j]
        total = np.bincount(cells, np.concatenate([face, face]), minlength=n)
        if mirror.size:
            total += np.bincount(mirror, w_m * offsets[:, i] * offsets[:, j],
                                 minlength=n)
        return total

    def rhs(i: int) -> np.ndarray:
        face = w * d[:, i] * dphi
        return np.bincount(cells, np.concatenate([face, face]), minlength=n)

    gxx, gyy, gzz = moment(0, 0), moment(1, 1), moment(2, 2)
    gxy, gxz, gyz = moment(0, 1), moment(0, 2), moment(1, 2)
    bx, by, bz = rhs(0), rhs(1), rhs(2)

    cxx = gyy * gzz - gyz * gyz
    cxy = gxz * gyz - gxy * gzz
    cxz = gxy * gyz - gxz * gyy
    cyy = gxx * gzz - gxz * gxz
    cyz = gxy * gxz - gxx * gyz
    czz = gxx * gyy - gxy * gxy
    det = gxx * cxx + gxy * cxy + gxz * cxz

    # Each neighbour adds a unit trace, so trace/3 is the scale a well-spread
    # stencil's eigenvalues have.
    scale = ((gxx + gyy + gzz) / 3.0) ** 3
    solvable = ok & (scale > 0.0) & (det > GRADIENT_CONDITION * scale)
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = np.where(solvable, 1.0 / det, np.nan)
        grad = np.stack([
            (cxx * bx + cxy * by + cxz * bz) * inv,
            (cxy * bx + cyy * by + cyz * bz) * inv,
            (cxz * bx + cyz * by + czz * bz) * inv,
        ], axis=1)
    return Gradient(np.where(solvable[:, None], grad, np.nan), solvable)


def cell_status(number_density, temperature, parcels_per_step, parcel_count,
                min_sampled_particles: float) -> np.ndarray:
    """The status flags of every cell; 0 is valid. See the module docstring."""
    n = np.asarray(number_density, dtype=np.float64)
    t = np.asarray(temperature, dtype=np.float64)
    per_step = np.asarray(parcels_per_step, dtype=np.float64)
    count = np.asarray(parcel_count, dtype=np.float64)

    status = np.zeros(n.shape, dtype=np.int64)
    empty = ~(per_step > 0.0)
    status[empty] |= NO_PARCELS
    status[~empty & ~(count > float(min_sampled_particles))] |= FEW_PARCELS
    status[~((n > 0.0) & np.isfinite(n))] |= BAD_DENSITY
    status[~((t > 0.0) & np.isfinite(t))] |= BAD_TEMPERATURE
    return status


# --------------------------------------------------------------------------- #
# the registry
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Output:
    """How a quantity is written as a field.

    ``kind`` is ``value`` -- invalid cells get the fill value -- or ``mask``,
    an integer code that is always defined and never filled.
    """

    dimensions: str
    units: str
    description: str
    kind: str = "value"


@dataclass(frozen=True)
class Quantity:
    """Something computable from inputs and other quantities.

    ``requires`` may be a callable taking the settings, for a quantity whose
    needs depend on them.
    """

    name: str
    requires: tuple | Callable
    compute: Callable
    output: Output | None = None

    def requirements(self, settings) -> tuple:
        return tuple(self.requires(settings) if callable(self.requires)
                     else self.requires)


#: Every quantity, in the order output fields are written.
REGISTRY: dict[str, Quantity] = {}


def register(name: str, requires, *, dimensions: str | None = None,
             units: str = "", description: str = "", kind: str = "value"):
    """Decorator: add a quantity, and make it an output if it has dimensions."""
    def decorate(compute):
        output = (Output(dimensions, units, description, kind)
                  if dimensions is not None else None)
        REGISTRY[name] = Quantity(name, requires, compute, output)
        return compute
    return decorate


class MissingInput(LookupError):
    """A quantity's inputs, followed all the way down, are not all there."""

    def __init__(self, name: str, missing: dict):
        self.name, self.missing = name, missing
        super().__init__(f"{name}: {explain(missing)}")


def explain(missing: dict) -> str:
    """``{input: reason}`` as one sentence."""
    return "needs " + "; ".join(f"{name} ({reason})"
                                for name, reason in missing.items())


class Context:
    """The inputs of one case at one time, and everything computed from them.

    Args:
        inputs: ``{input name: value}`` -- whatever the I/O layer could read.
        species: the gas; see :func:`mean_free_path`. Also needs
            ``gas_constant_j_per_kg_k``.
        settings: a :class:`~plumetools.rarefaction.settings.DiagnosticSettings`.
        convention: the mean-free-path convention.
        unavailable: ``{input name: why it could not be read}``, so a skipped
            diagnostic says what to do about it.
        extra: quantities beyond :data:`REGISTRY` -- the per-length Knudsen
            numbers, for one.

    Quantities are computed on first use and kept.
    """

    def __init__(self, inputs: dict, *, species, settings,
                 convention: str = "vhs", unavailable: dict | None = None,
                 extra: dict | None = None):
        self.inputs = dict(inputs)
        self.species = species
        self.settings = settings
        self.convention = convention
        self.unavailable = dict(unavailable or {})
        # The per-length Knudsen numbers are written beside the mean free path.
        self.registry: dict = {}
        for name, quantity in REGISTRY.items():
            self.registry[name] = quantity
            if name == "mfp":
                self.registry.update(extra or {})
        self.registry.update(extra or {})
        self._cache: dict = {}

    def missing(self, name: str) -> dict:
        """``{input: reason}`` for every input ``name`` needs and lacks."""
        if name in self.inputs:
            return {}
        quantity = self.registry.get(name)
        if quantity is None:
            return {name: self.unavailable.get(name, "not provided")}
        out: dict = {}
        for requirement in quantity.requirements(self.settings):
            out.update(self.missing(requirement))
        return out

    def __getitem__(self, name: str):
        if name in self.inputs:
            return self.inputs[name]
        if name not in self._cache:
            missing = self.missing(name)
            if missing:
                raise MissingInput(name, missing)
            self._cache[name] = self.registry[name].compute(self)
        return self._cache[name]


def knudsen_on(name: str, length_m: float) -> Quantity:
    """The quantity ``Kn<name> = lambda / L``, for a characteristic length."""
    length_m = float(length_m)

    def compute(ctx):
        lam = ctx["mean_free_path"]
        return Field.flagged(lam.values / length_m, lam.valid)

    return Quantity(
        f"Kn{name}", ("mean_free_path",), compute,
        Output(DIMENSIONLESS, "-",
               f"lambda / {name}, {name} = {length_m:g} m: rarefaction on that "
               f"length scale"))


def output_names(ctx: Context, select=None) -> list[str]:
    """The output fields to evaluate, in registry order.

    Raises:
        ValueError: for a ``select`` entry that is not an output.
    """
    names = [name for name, q in ctx.registry.items() if q.output is not None]
    if select is None:
        return names
    unknown = sorted(set(select) - set(names))
    if unknown:
        raise ValueError(f"diagnostics.select names {unknown}, which are not "
                         f"outputs; valid: {names}")
    return [name for name in names if name in set(select)]


def evaluate(ctx: Context, names) -> tuple[dict, dict]:
    """Compute each named output that the case can supply.

    Returns:
        ``(computed, skipped)``: ``{name: Field}`` and ``{name: reason}``. A
        diagnostic whose inputs are missing is skipped; the rest still run.
    """
    computed, skipped = {}, {}
    for name in names:
        missing = ctx.missing(name)
        if missing:
            skipped[name] = explain(missing)
            continue
        computed[name] = ctx[name]
    return computed, skipped


# --------------------------------------------------------------------------- #
# intermediate quantities
# --------------------------------------------------------------------------- #

@register("temperature", ("mass_density", "momentum_density",
                          "kinetic_energy_density"))
def _temperature(ctx):
    return translational_temperature(
        ctx["mass_density"], ctx["momentum_density"],
        ctx["kinetic_energy_density"], ctx.species.gas_constant_j_per_kg_k)


def _count_requirements(settings):
    if settings.particle_count_basis == "accumulated":
        return ("parcels_per_step", "averaged_steps")
    return ("parcels_per_step",)


@register("parcel_count", _count_requirements)
def _parcel_count(ctx):
    per_step = np.asarray(ctx["parcels_per_step"], dtype=np.float64)
    if ctx.settings.particle_count_basis == "accumulated":
        return per_step * float(ctx["averaged_steps"])
    return per_step


@register("status", ("number_density", "temperature", "parcels_per_step",
                     "parcel_count"))
def _status(ctx):
    return cell_status(ctx["number_density"], ctx["temperature"],
                       ctx["parcels_per_step"], ctx["parcel_count"],
                       ctx.settings.min_sampled_particles)


@register("valid", ("status",))
def _valid(ctx):
    return ctx["status"] == 0


@register("computable", ("status",))
def _computable(ctx):
    """Sampled, with n > 0 and a resolved T > 0: a value exists, trusted or not."""
    return (ctx["status"] & UNCOMPUTABLE) == 0


@register("mean_free_path", ("number_density", "temperature", "valid",
                             "computable"))
def _mean_free_path(ctx):
    lam = mean_free_path(ctx["number_density"], ctx["temperature"], ctx.species,
                         convention=ctx.convention)
    return Field.flagged(np.where(ctx["computable"], lam, np.nan), ctx["valid"])


def _two_pass_gradient(ctx, name: str) -> Gradient:
    """A valid cell's gradient from valid neighbours only -- what is trusted
    and summarised. Every other computable cell gets one from its computable
    neighbours, to be written and looked at, never trusted."""
    phi, graph = ctx[name], ctx["cell_graph"]
    trusted = least_squares_gradient(phi, ctx["valid"], graph)
    shown = least_squares_gradient(phi, ctx["computable"], graph)
    vectors = np.where(trusted.valid[:, None], trusted.vectors, shown.vectors)
    return Gradient(vectors, trusted.valid)


@register("grad_number_density", ("number_density", "valid", "computable",
                                  "cell_graph"))
def _grad_n(ctx):
    return _two_pass_gradient(ctx, "number_density")


@register("grad_temperature", ("temperature", "valid", "computable",
                               "cell_graph"))
def _grad_t(ctx):
    return _two_pass_gradient(ctx, "temperature")


# --------------------------------------------------------------------------- #
# outputs, in the order they are written
# --------------------------------------------------------------------------- #

@register("mfp", ("mean_free_path",), dimensions=LENGTH, units="m",
          description="local mean free path lambda, from the averaged n and "
                      "translational T (equation G1)")
def _mfp(ctx):
    return ctx["mean_free_path"]


@register("cellSizeV", ("cell_volume",), dimensions=LENGTH, units="m",
          description="volume-equivalent cell size V^(1/3)")
def _cell_size_v(ctx):
    volume = np.asarray(ctx["cell_volume"], dtype=np.float64)
    return Field.masked(np.cbrt(volume), volume > 0.0)


@register("cellSizeMax", ("cell_extent",), dimensions=LENGTH, units="m",
          description="largest bounding-box dimension of the cell")
def _cell_size_max(ctx):
    extent = np.asarray(ctx["cell_extent"], dtype=np.float64)
    size = extent.max(axis=1)
    return Field.masked(size, size > 0.0)


def _ratio(numerator: Field, denominator: Field) -> Field:
    with np.errstate(divide="ignore", invalid="ignore"):
        return Field.flagged(numerator.values / denominator.values,
                            numerator.valid & denominator.valid)


@register("KnCellV", ("mean_free_path", "cellSizeV"), dimensions=DIMENSIONLESS,
          units="-", description="lambda / V^(1/3): mesh resolution, not "
                                 "rarefaction")
def _kn_cell_v(ctx):
    return _ratio(ctx["mean_free_path"], ctx["cellSizeV"])


@register("KnCellMax", ("mean_free_path", "cellSizeMax"), dimensions=DIMENSIONLESS,
          units="-", description="lambda / Delta_max: mesh resolution, not "
                                 "rarefaction")
def _kn_cell_max(ctx):
    return _ratio(ctx["mean_free_path"], ctx["cellSizeMax"])


@register("cellOverMfpV", ("mean_free_path", "cellSizeV"), dimensions=DIMENSIONLESS,
          units="-", description="V^(1/3) / lambda; above 1/k the cell is "
                                 "coarser than the target")
def _cell_over_mfp_v(ctx):
    return _ratio(ctx["cellSizeV"], ctx["mean_free_path"])


@register("cellOverMfpMax", ("mean_free_path", "cellSizeMax"),
          dimensions=DIMENSIONLESS, units="-",
          description="Delta_max / lambda; above 1/k the cell is coarser than "
                      "the target")
def _cell_over_mfp_max(ctx):
    return _ratio(ctx["cellSizeMax"], ctx["mean_free_path"])


@register("mfpResolved", ("KnCellMax",), dimensions=DIMENSIONLESS,
          units="1 yes / 0 no / -1 no data",
          description="the mesh criterion: 1 where lambda / Delta_max >= k, the "
                      "cell resolves the local mean free path; 0 where it does "
                      "not; -1 where there is no value. Trusted only where "
                      "diagValid = 1", kind="mask")
def _mfp_resolved(ctx):
    kn = ctx["KnCellMax"].values
    k = float(ctx.settings.cells_per_mean_free_path)
    with np.errstate(invalid="ignore"):
        code = np.where(np.isfinite(kn), (kn >= k).astype(np.float64), -1.0)
    return Field(code, np.ones(code.shape, dtype=bool))


def _gradient_knudsen(ctx, field_name: str, gradient_name: str) -> Field:
    lam = ctx["mean_free_path"]
    phi = np.asarray(ctx[field_name], dtype=np.float64)
    grad = ctx[gradient_name]
    with np.errstate(divide="ignore", invalid="ignore"):
        values = lam.values * grad.magnitude / phi
    # A valid zero gradient is a valid zero Knudsen number.
    return Field.flagged(values, lam.valid & grad.valid)


@register("KnGradN", ("mean_free_path", "number_density", "grad_number_density"),
          dimensions=DIMENSIONLESS, units="-",
          description="lambda |grad n| / n: rarefaction on the density's own "
                      "length scale")
def _kn_grad_n(ctx):
    return _gradient_knudsen(ctx, "number_density", "grad_number_density")


@register("KnGradT", ("mean_free_path", "temperature", "grad_temperature"),
          dimensions=DIMENSIONLESS, units="-",
          description="lambda |grad T| / T: rarefaction on the temperature's own "
                      "length scale")
def _kn_grad_t(ctx):
    return _gradient_knudsen(ctx, "temperature", "grad_temperature")


@register("cellSizeTarget", ("mean_free_path",), dimensions=LENGTH, units="m",
          description="lambda / k, the cell size that resolves the local mean "
                      "free path. Never clipped to a practical minimum")
def _cell_size_target(ctx):
    lam = ctx["mean_free_path"]
    return Field.flagged(lam.values / float(ctx.settings.cells_per_mean_free_path),
                        lam.valid)


@register("parcelsPerStep", ("parcels_per_step",), dimensions=DIMENSIONLESS,
          units="parcels", description="simulated parcels in the cell, averaged "
                                       "over the steps sampled (dsmcRhoNMean)")
def _parcels_per_step(ctx):
    per_step = np.asarray(ctx["parcels_per_step"], dtype=np.float64)
    return Field.masked(per_step, np.ones(per_step.shape, dtype=bool))


@register("sampledParcels", ("parcels_per_step", "averaged_steps"),
          dimensions=DIMENSIONLESS, units="parcels",
          description="parcel samples the cell's means are built from: "
                      "parcels per step x steps averaged")
def _sampled_parcels(ctx):
    total = (np.asarray(ctx["parcels_per_step"], dtype=np.float64)
             * float(ctx["averaged_steps"]))
    return Field.masked(total, np.ones(total.shape, dtype=bool))


@register("diagStatus", ("status",), dimensions=DIMENSIONLESS, units="flags",
          description="0 valid; else a sum of 1 no parcels, 2 too few sampled "
                      "parcels, 4 n <= 0, 8 no resolved T > 0", kind="mask")
def _diag_status(ctx):
    status = ctx["status"].astype(np.float64)
    return Field(status, np.ones(status.shape, dtype=bool))


@register("diagValid", ("valid",), dimensions=DIMENSIONLESS, units="0/1",
          description="1 where the cell's diagnostics are valid", kind="mask")
def _diag_valid(ctx):
    valid = ctx["valid"]
    return Field(valid.astype(np.float64), np.ones(valid.shape, dtype=bool))


@register("diagValidGradN", ("grad_number_density",), dimensions=DIMENSIONLESS,
          units="0/1", description="1 where grad n could be fitted", kind="mask")
def _diag_valid_grad_n(ctx):
    ok = ctx["grad_number_density"].valid
    return Field(ok.astype(np.float64), np.ones(ok.shape, dtype=bool))


@register("diagValidGradT", ("grad_temperature",), dimensions=DIMENSIONLESS,
          units="0/1", description="1 where grad T could be fitted", kind="mask")
def _diag_valid_grad_t(ctx):
    ok = ctx["grad_temperature"].valid
    return Field(ok.astype(np.float64), np.ones(ok.shape, dtype=bool))
