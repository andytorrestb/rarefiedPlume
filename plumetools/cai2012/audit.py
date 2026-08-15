r"""What the sampled solution actually carried: occupancy, overlap, and the sweep.

Three measurements that ``cases/cai2012`` does not make, all of them on output
that already exists on disk.

**1. The occupancy audit.** :func:`plumetools.cai2012.checks._check_occupancy`
prints an *estimate* of the exit cell's parcel count and says outright that only
a post-run audit of the sampled ``dsmcRhoN`` field can say what it was.
:func:`occupancy_audit` is that audit.

``dsmcRhoN`` is **the parcel count in the cell**, dimensionless -- not a parcel
number density. ``DSMCCloud::calculateFields`` adds 1 per parcel and never
divides by the cell volume, where ``rhoN`` accumulates ``nParticle/V``. So with
one global particle weight:

.. code-block:: text

    dsmcRhoN            parcels in the cell     [-]
    rhoN * V / dsmcRhoN = nParticle             the global weight

and the time-averaged ``dsmcRhoNMean`` is the mean occupancy over the sampling
window directly, with no volume factor.

This is worth stating flatly because the repository had it wrong until this
family measured it: ``catalog.py``, ``slices.yaml`` and ``cases/cai2012/viz.yaml``
all described ``dsmcRhoN`` as a number density in ``m^-3`` whose product with
the cell volume was the occupancy. Measured on ``cases/cai2012/Cases/Kn100``,
``rhoN*V/dsmcRhoN`` is the particle weight to 1.5e-9 and ``dsmcRhoN`` peaks at
20.47 -- the configured 20 parcels per exit cell. The wrong reading is out by a
factor of ``1/V``, which on this mesh is 10^6, and on a logarithmic colour scale
it looks entirely plausible.

The identity is therefore **checked, not assumed** (:func:`weight_consistency`):
if ``rhoNMean * V / dsmcRhoNMean`` is not the particle weight everywhere, the
reading of these fields is wrong again and every occupancy number here is wrong
with it.

**2. The overlap check.** A 4.5-transit run is a 1.5-transit run that kept
going. Same mesh, same seed, same time step, same decomposition -- so at any
time both have reached, the two ``rhoNMean`` fields are the *same numbers*, and
:func:`compare_overlap` asserts it byte for byte. A disagreement is not a small
statistical wobble to be tolerated; it means something in the
generate-mesh-run-sample chain is not deterministic, and that finding outranks
everything else this family measures.

**3. The sweep.** :func:`sweep_rows` collects the nine cases' metrics,
occupancy and statistical budget into one table, so the two axes can be compared
against the one quantity that ought to govern both -- the total number of
samples the average was built from.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from plumetools.cai2012.post import PostError, read_internal_field

#: Fraction of ``n0`` below which a cell is outside the plume Cai draws.
#:
#: His lowest contour. Occupancy in the far vacuum is meaninglessly low and
#: says nothing about the quality of the answer -- no parcel was ever supposed
#: to be there. Quoting a domain-wide median instead of a plume median is the
#: easiest way to make a well-resolved run look terrible.
PLUME_FLOOR_OVER_N0 = 1.0e-3


# --------------------------------------------------------------------------- #
# 1. occupancy
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class RegionOccupancy:
    """Parcels per cell over one region of the mesh.

    Attributes:
        name: the region.
        description: what it covers, for the report.
        n_cells: cells in it.
        mean, median, p10, p90, minimum, maximum: parcels per cell.
        below_floor: fraction of cells under ``checks.min_particles_per_cell``.
        below_target: fraction under ``resolution.target_particles_per_cell``.
        samples_per_cell: median occupancy times the number of averaged steps
            -- the statistical budget the mean in a typical cell was built
            from, which is the quantity both sweep axes actually move.
    """

    name: str
    description: str
    n_cells: int
    mean: float
    median: float
    p10: float
    p90: float
    minimum: float
    maximum: float
    below_floor: float
    below_target: float
    samples_per_cell: float

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "description": self.description,
            "n_cells": self.n_cells,
            "mean": self.mean,
            "median": self.median,
            "p10": self.p10,
            "p90": self.p90,
            "min": self.minimum,
            "max": self.maximum,
            "fraction_below_floor": self.below_floor,
            "fraction_below_target": self.below_target,
            "samples_per_cell": self.samples_per_cell,
        }


def region_occupancy(name: str, description: str, occupancy: np.ndarray, *,
                     floor: float, target: float,
                     averaged_steps: float) -> RegionOccupancy:
    """Summarise parcels per cell over a selection of cells.

    Args:
        name: region name.
        description: what the selection covers.
        occupancy: parcels per cell, one entry per cell in the region.
        floor: ``checks.min_particles_per_cell``.
        target: ``resolution.target_particles_per_cell``.
        averaged_steps: time steps ``fieldAverage`` accumulated over.

    Returns:
        A :class:`RegionOccupancy`. An empty selection gives zeros rather than
        raising: a region can legitimately be empty (no cell of the mesh in
        that band), and a study that stopped on it would be reporting a mesh
        property as a failure.

    The **median** rather than the mean is the headline: occupancy spans orders
    of magnitude between the exit and the plume edge, and a mean over that is
    dominated by the few densest cells -- it would report the exit cell's
    occupancy as though it were typical of the region.
    """
    values = np.asarray(occupancy, dtype=np.float64).ravel()
    if values.size == 0:
        return RegionOccupancy(name, description, 0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0)
    median = float(np.median(values))
    return RegionOccupancy(
        name=name,
        description=description,
        n_cells=int(values.size),
        mean=float(values.mean()),
        median=median,
        p10=float(np.percentile(values, 10.0)),
        p90=float(np.percentile(values, 90.0)),
        minimum=float(values.min()),
        maximum=float(values.max()),
        below_floor=float(np.count_nonzero(values < floor) / values.size),
        below_target=float(np.count_nonzero(values < target) / values.size),
        samples_per_cell=median * float(averaged_steps),
    )


def weight_consistency(number_density: np.ndarray, parcel_count: np.ndarray,
                       volumes: np.ndarray,
                       n_equivalent_particles: float) -> dict:
    """Check ``rhoN * V / dsmcRhoN`` really is the global particle weight.

    Args:
        number_density: ``rhoNMean`` [1/m^3].
        parcel_count: ``dsmcRhoNMean`` -- parcels per cell, dimensionless.
        volumes: cell volumes [m^3].
        n_equivalent_particles: the weight from ``case.yaml``/the manifest.

    Returns:
        ``{"n_compared", "measured_weight", "max_rel_error", "consistent"}``.

    Every occupancy figure in this module rests on two things: ``dsmcRhoN``
    being a parcel **count**, and one weight covering the whole domain
    (``DSMCCloud::nParticle_`` is a single scalar). Both can fail silently. The
    repository already had the first one wrong -- it read ``dsmcRhoN`` as a
    number density, which is out by ``1/V`` -- and a solver with radial
    weighting would break the second. Either way the ratio stops being
    constant, so it is measured rather than trusted.

    Cells with no parcels are skipped: ``0/0`` says nothing about the weight.
    """
    number_density = np.asarray(number_density, dtype=np.float64).ravel()
    parcel_count = np.asarray(parcel_count, dtype=np.float64).ravel()
    volumes = np.asarray(volumes, dtype=np.float64).ravel()
    occupied = (parcel_count > 0.0) & (number_density > 0.0)
    if not np.any(occupied):
        return {"n_compared": 0, "measured_weight": float("nan"),
                "max_rel_error": float("nan"), "consistent": False}

    ratio = (number_density[occupied] * volumes[occupied]
             / parcel_count[occupied])
    expected = float(n_equivalent_particles)
    error = (float(np.max(np.abs(ratio - expected)) / expected)
             if expected else float("inf"))
    return {
        "n_compared": int(np.count_nonzero(occupied)),
        "measured_weight": float(np.median(ratio)),
        "max_rel_error": error,
        # Loose, because both fields are written at 10 significant figures and
        # the ratio of two rounded numbers is not exact -- measured at 1.5e-9
        # on cases/cai2012/Cases/Kn100. Tight enough that a second weight, or a
        # missing volume factor, could not pass.
        "consistent": error < 1.0e-6,
    }


def occupancy_audit(sampled, cfg, geom, exit_state, *,
                    n_equivalent_particles: float,
                    parcel_count: np.ndarray,
                    averaged_steps: float,
                    core_cell_size_m: float | None = None,
                    x_stations_over_D=(0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0)) -> dict:
    """Measure parcels per cell, per region, from the sampled fields.

    Args:
        sampled: a :class:`~plumetools.cai2012.post.SampledField`.
        cfg: the case config.
        geom: the geometry.
        exit_state: the exit state, for ``n0``.
        n_equivalent_particles: the global particle weight.
        parcel_count: ``dsmcRhoNMean`` -- the mean number of parcels in each
            cell over the sampling window. Dimensionless, and **not** to be
            multiplied by the cell volume; see the module docstring.
        averaged_steps: time steps ``fieldAverage`` accumulated over.
        core_cell_size_m: the uniform core cell, from
            :func:`plumetools.cai2012.mesh.plan`. ``None`` derives it from the
            cell centres, since ``mesh.core_cell_size_m`` is normally ``null``
            and the real value only exists after the mesh is planned.
        x_stations_over_D: upper edges of the centreline bands to report.

    Returns:
        A dict with ``weight_check``, ``regions`` and ``centerline_stations``.

    Raises:
        PostError: if ``dsmcRhoNMean`` does not match the mesh.

    The regions are chosen to answer one question each:

    ==================  ==========================================================
    ``exit``            what ``checks.py`` estimates -- the only place the target
                        is met by construction
    ``core``            the uniform fine block, where the mesh is one cell size
                        and occupancy differences are density differences
    ``plume``           everywhere inside Cai's lowest contour: the region whose
                        numbers are supposed to mean something
    ``centerline``      the tube ``post.centerline`` averages over, which is what
                        the published errors are computed from
    ``domain``          everything, including the vacuum no parcel ever reached
    ==================  ==========================================================
    """
    parcel_count = np.asarray(parcel_count, dtype=np.float64).ravel()
    if parcel_count.size != sampled.n_cells:
        raise PostError(
            f"dsmcRhoNMean has {parcel_count.size:,} values against "
            f"{sampled.n_cells:,} cells. These are not the same mesh.")

    # No volume factor. dsmcRhoN IS the count; see the module docstring.
    occupancy = parcel_count

    d = geom.diameter_m
    x = sampled.centres[:, 0]
    radial = np.hypot(sampled.centres[:, 1], sampled.centres[:, 2])
    ratio = sampled.number_density / exit_state.number_density_per_m3

    floor = float(cfg.checks.min_particles_per_cell)
    target = float(cfg.resolution.target_particles_per_cell)

    core_x = float(cfg.mesh.core_x_over_D) * d
    core_half = float(cfg.mesh.core_half_over_D) * d
    # The exit cells: the first layer of the uniform core inside the nozzle
    # disk. One core cell deep -- the same cells whose volume set the particle
    # weight, so this is the number checks.py was estimating.
    core_cell = float(core_cell_size_m or cfg.mesh.core_cell_size_m
                      or _core_cell(x))
    # One layer, not 1.5 of them. Cell centres in the uniform core sit at
    # 0.5, 1.5, 2.5 ... cells, so a 1.5-cell cut takes the second layer too --
    # and the second layer is not the cell whose volume set the particle
    # weight, so the "measured against estimated" comparison would be against
    # a different cell than checks.py estimated.
    first_layer = x <= core_cell
    tube = radial <= float(cfg.post.centerline_radius_over_D) * d
    in_range = (x >= 0.0) & (x <= float(cfg.post.centerline_x_over_D_max) * d)

    selections = [
        ("exit", "first core-cell layer inside the nozzle disk (r <= R0)",
         first_layer & (radial <= geom.radius_m)),
        ("core", "the uniform fine block around the nozzle and near plume",
         (x >= 0.0) & (x <= core_x) & (np.abs(sampled.centres[:, 1]) <= core_half)
         & (np.abs(sampled.centres[:, 2]) <= core_half)),
        ("plume", f"every cell with n/n0 >= {PLUME_FLOOR_OVER_N0:g} "
                  f"(Cai's lowest contour)",
         ratio >= PLUME_FLOOR_OVER_N0),
        ("centerline", f"the r <= {cfg.post.centerline_radius_over_D:g} D tube "
                       f"the published errors are averaged over",
         tube & in_range),
        ("domain", "every cell, vacuum included", np.ones_like(x, dtype=bool)),
    ]

    regions = [
        region_occupancy(name, description, occupancy[mask], floor=floor,
                         target=target, averaged_steps=averaged_steps)
        for name, description, mask in selections
    ]

    stations = []
    lower = 0.0
    for upper in x_stations_over_D:
        band = tube & (x > lower * d) & (x <= upper * d)
        summary = region_occupancy(
            f"X/D {lower:g}-{upper:g}",
            "centreline tube", occupancy[band], floor=floor, target=target,
            averaged_steps=averaged_steps)
        stations.append({"x_over_D_min": float(lower),
                         "x_over_D_max": float(upper),
                         **summary.as_dict()})
        lower = upper

    return {
        "weight_check": weight_consistency(sampled.number_density,
                                           parcel_count, sampled.volumes,
                                           n_equivalent_particles),
        "n_equivalent_particles": float(n_equivalent_particles),
        "averaged_steps": float(averaged_steps),
        "target_particles_per_cell": target,
        "min_particles_per_cell": floor,
        "regions": [r.as_dict() for r in regions],
        "centerline_stations": stations,
    }


def _core_cell(x: np.ndarray) -> float:
    """Fallback core cell size, from the mesh itself.

    ``mesh.core_cell_size_m`` is normally ``null`` -- derived at mesh time --
    so the audit takes the smallest positive cell-centre abscissa, which on a
    uniform core starting at ``x = 0`` is half a cell.
    """
    positive = x[x > 0.0]
    return 2.0 * float(positive.min()) if positive.size else 1.0


def occupancy_report(audit: dict) -> list[str]:
    """The audit as printed lines."""
    check = audit["weight_check"]
    lines = [
        "Measured parcel occupancy (dsmcRhoNMean -- the parcel COUNT per cell)",
        f"  particle weight          {audit['n_equivalent_particles']:.6e}"
        f"   (measured {check['measured_weight']:.6e}, "
        f"max error {check['max_rel_error']:.2e})",
    ]
    if not check["consistent"]:
        lines.append(
            "  WARNING: rhoNMean * V / dsmcRhoNMean is NOT the particle weight. "
            "Every occupancy\n"
            "           number below assumes dsmcRhoN is a parcel COUNT and "
            "that one weight\n"
            "           covers the domain. One of those is false here.")
    lines += [
        f"  averaged over            {audit['averaged_steps']:.0f} time steps",
        "",
        f"  {'region':<12} {'cells':>10} {'median':>9} {'p10':>9} {'p90':>9} "
        f"{'<floor':>7} {'samples':>10}",
    ]
    for region in audit["regions"]:
        lines.append(
            f"  {region['name']:<12} {region['n_cells']:>10,} "
            f"{region['median']:>9.2f} {region['p10']:>9.2f} "
            f"{region['p90']:>9.2f} "
            f"{100.0 * region['fraction_below_floor']:>6.1f}% "
            f"{region['samples_per_cell']:>10.3e}")
    lines += [
        "",
        f"  {'station':<12} {'cells':>10} {'median':>9} {'<floor':>7}",
    ]
    for station in audit["centerline_stations"]:
        lines.append(
            f"  {station['name']:<12} {station['n_cells']:>10,} "
            f"{station['median']:>9.2f} "
            f"{100.0 * station['fraction_below_floor']:>6.1f}%")
    return lines


# --------------------------------------------------------------------------- #
# 2. the overlap
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class FieldComparison:
    """Two runs' versions of one field at one time.

    Attributes:
        field: the field name.
        time: the time directory both were read from.
        identical: byte-for-byte equal.
        n_values: how many numbers were compared.
        max_abs_diff, max_rel_diff: zero when ``identical``.
        note: why a numeric comparison was impossible, if it was.
    """

    field: str
    time: str
    identical: bool
    n_values: int = 0
    max_abs_diff: float = 0.0
    max_rel_diff: float = 0.0
    note: str = ""

    def as_dict(self) -> dict:
        return {"field": self.field, "time": self.time,
                "identical": self.identical, "n_values": self.n_values,
                "max_abs_diff": self.max_abs_diff,
                "max_rel_diff": self.max_rel_diff, "note": self.note}


def file_digest(path: Path) -> str:
    """SHA-256 of a field file.

    The fast path for the overlap check. Two runs that took the same
    trajectory write the same ASCII at the same ten significant figures, so
    equal digests settle it without parsing 20 MB of numbers twice. Unequal
    digests do not settle anything -- a differing header would do it -- so the
    numbers are then read and compared.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def compare_field(long_dir: Path, short_dir: Path, field: str,
                  time: str) -> FieldComparison:
    """Compare one field between two cases at one time.

    Args:
        long_dir: the case that kept running.
        short_dir: the case that stopped earlier.
        field: e.g. ``rhoNMean``.
        time: a time directory both wrote.

    Returns:
        A :class:`FieldComparison`.

    Raises:
        PostError: if either file is missing -- the caller chose this time
            because both cases were supposed to have it, so an absent file is a
            bug in the selection, not a result.
    """
    left = Path(long_dir) / time / field
    right = Path(short_dir) / time / field
    for path in (left, right):
        if not path.is_file():
            raise PostError(
                f"no {field} at {path}; the overlap check was told both cases "
                f"wrote time {time}")

    if file_digest(left) == file_digest(right):
        return FieldComparison(field=field, time=time, identical=True,
                               n_values=0, note="byte for byte")

    a = read_internal_field(left).ravel()
    b = read_internal_field(right).ravel()
    if a.size != b.size:
        return FieldComparison(
            field=field, time=time, identical=False,
            note=f"{a.size:,} values against {b.size:,} -- different meshes, "
                 f"so this is not a determinism result")

    difference = np.abs(a - b)
    scale = np.maximum(np.abs(a), np.abs(b))
    nonzero = scale > 0.0
    relative = np.zeros_like(difference)
    relative[nonzero] = difference[nonzero] / scale[nonzero]
    return FieldComparison(
        field=field, time=time, identical=False, n_values=int(a.size),
        max_abs_diff=float(difference.max()),
        max_rel_diff=float(relative.max()),
    )


def common_times(long_dir: Path, short_dir: Path) -> list[str]:
    """Time directories both cases wrote, ascending.

    Matched as **strings**, which is what OpenFOAM named them: two runs at the
    same time step and the same write interval produce the same names, and a
    float comparison here would invent a tolerance where there is an exact
    answer.
    """
    def times(directory: Path) -> dict:
        found = {}
        for entry in Path(directory).iterdir():
            if not entry.is_dir():
                continue
            try:
                value = float(entry.name)
            except ValueError:
                continue
            if value > 0.0:
                found[entry.name] = value
        return found

    left, right = times(long_dir), times(short_dir)
    shared = set(left) & set(right)
    return sorted(shared, key=lambda name: left[name])


def compare_overlap(long_dir: Path, short_dir: Path, *,
                    fields=("rhoNMean", "dsmcRhoNMean", "momentumMean"),
                    times=None) -> dict:
    """Check that a long run and a short run agree where they overlap.

    Args:
        long_dir: the case that kept running.
        short_dir: the case that stopped earlier.
        fields: which averaged fields to compare.
        times: the times to check, or ``None`` for every one both wrote.

    Returns:
        ``{"long", "short", "n_times", "n_compared", "all_identical",
        "worst", "comparisons"}``.

    What this establishes
    ---------------------
    The two cases differ in ``endTime`` and in nothing else: same mesh, same
    particle weight, same time step, same decomposition, same RNG seeding
    (``dsmcFoam`` seeds per rank from a fixed constant). The short run is
    therefore a **prefix** of the long one, and ``fieldAverage`` -- which starts
    at the same ``timeStart`` in both and stores its accumulators per time
    directory -- must have accumulated exactly the same numbers by any time both
    reached.

    So an agreement is not evidence that the physics is right. It is evidence
    that the pipeline is **reproducible**, which is the precondition for every
    other comparison in this family: six of the nine cases are partly redundant
    by construction, and if the redundancy does not hold, the differences the
    sweep attributes to statistics are partly something else.
    """
    long_dir, short_dir = Path(long_dir), Path(short_dir)
    selected = list(times) if times is not None else common_times(long_dir,
                                                                  short_dir)
    comparisons = []
    for time in selected:
        for name in fields:
            comparisons.append(compare_field(long_dir, short_dir, name, time))

    differing = [c for c in comparisons if not c.identical]
    worst = max(differing, key=lambda c: c.max_rel_diff, default=None)
    return {
        "long": str(long_dir).replace("\\", "/"),
        "short": str(short_dir).replace("\\", "/"),
        "n_times": len(selected),
        "times": selected,
        "n_compared": len(comparisons),
        "all_identical": not differing,
        "worst": worst.as_dict() if worst is not None else None,
        "comparisons": [c.as_dict() for c in comparisons],
    }


def overlap_report(result: dict) -> list[str]:
    """The overlap check as printed lines, with the verdict spelled out."""
    lines = [
        f"Overlap: {result['short']} against {result['long']}",
        f"  {result['n_compared']} field(s) over {result['n_times']} shared "
        f"time(s)",
    ]
    if result["n_compared"] == 0:
        lines.append(
            "  NOTHING COMPARED. The two cases share no time directory, which "
            "means the\n"
            "  write schedules differ -- so they are not the same run "
            "truncated, and the\n"
            "  matrix's redundancy assumption does not hold.")
        return lines
    if result["all_identical"]:
        lines.append(
            "  IDENTICAL, byte for byte. The short run is a prefix of the long "
            "one, as it\n"
            "  must be: same mesh, same weight, same step, same seed. The "
            "pipeline is\n"
            "  reproducible.")
        return lines
    worst = result["worst"] or {}
    lines += [
        f"  DIFFERENT. Worst: {worst.get('field')} at t = {worst.get('time')}, "
        f"max relative {worst.get('max_rel_diff', float('nan')):.3e}",
        "  These two runs differ only in endTime, so agreeing was not optional. "
        "Something",
        "  in generate -> mesh -> run -> sample is not deterministic, and that "
        "is a more",
        "  important finding than anything else in this study: the sweep "
        "attributes its",
        "  differences to statistics, and this says part of them are something "
        "else.",
    ]
    if worst.get("note"):
        lines.append(f"  note: {worst['note']}")
    return lines


# --------------------------------------------------------------------------- #
# 3. the sweep table
# --------------------------------------------------------------------------- #

#: Columns of ``results/sweep-table.csv``, in order.
SWEEP_COLUMNS = (
    "case", "particles_per_cell", "sampling_transits",
    "exit_occupancy_estimated", "exit_occupancy_measured",
    "plume_median_occupancy", "plume_fraction_below_floor",
    "centerline_median_occupancy", "samples_per_cell", "statistical_budget",
    "density_max_rel_error_percent", "density_mean_rel_error_percent",
    "density_rms_rel_error_percent", "velocity_max_rel_error_percent",
    "temperature_max_rel_error_percent",
)


def statistical_budget(particles_per_cell: float, sampling_transits: float,
                       steps_per_transit: float) -> float:
    """Samples the average in a target-occupancy cell was built from.

    .. code-block:: text

        budget = parcels per cell  x  time steps averaged over

    The quantity the whole sweep is about. Both axes multiply it and neither is
    supposed to matter separately -- doubling the parcels and halving the
    averaging leaves it unchanged -- so plotting error against *this* rather
    than against either axis is the test of whether that is true. Where the two
    axes stop being interchangeable, something other than sampling noise is
    setting the error.

    Successive samples of a DSMC field are **not** independent -- a parcel
    survives many steps -- so this is a budget, not an effective sample size.
    It is proportional to the true one for a fixed flow, which is all the
    comparison needs.
    """
    return float(particles_per_cell) * float(sampling_transits) * float(
        steps_per_transit)


def sweep_row(entry: dict, metrics: dict, audit: dict | None,
              steps_per_transit: float) -> dict:
    """One row of the sweep table.

    Args:
        entry: the case's manifest entry.
        metrics: its ``results/metrics.yaml``.
        audit: its ``results/occupancy.yaml``, or ``None`` if not run.
        steps_per_transit: time steps in one domain transit.

    Returns:
        A dict keyed by :data:`SWEEP_COLUMNS`.
    """
    errors = (metrics or {}).get("centerline", {})
    regions = {r["name"]: r for r in (audit or {}).get("regions", [])}
    plume = regions.get("plume", {})
    centre = regions.get("centerline", {})
    exit_region = regions.get("exit", {})

    return {
        "case": entry.get("name"),
        "particles_per_cell": entry.get("target_particles_per_cell"),
        "sampling_transits": entry.get("sampling_domain_transits"),
        "exit_occupancy_estimated": entry.get("exit_particles_per_cell_estimated"),
        "exit_occupancy_measured": exit_region.get("median"),
        "plume_median_occupancy": plume.get("median"),
        "plume_fraction_below_floor": plume.get("fraction_below_floor"),
        "centerline_median_occupancy": centre.get("median"),
        "samples_per_cell": centre.get("samples_per_cell"),
        "statistical_budget": statistical_budget(
            entry.get("target_particles_per_cell") or 0.0,
            entry.get("sampling_domain_transits") or 0.0,
            steps_per_transit),
        "density_max_rel_error_percent": _percent(errors.get("density_max_rel_error")),
        "density_mean_rel_error_percent": _percent(errors.get("density_mean_rel_error")),
        "density_rms_rel_error_percent": _percent(errors.get("density_rms_rel_error")),
        "velocity_max_rel_error_percent": _percent(errors.get("velocity_max_rel_error")),
        "temperature_max_rel_error_percent": _percent(
            errors.get("temperature_max_rel_error")),
    }


def _percent(value):
    return None if value is None else 100.0 * float(value)


def write_sweep_table(path: Path, rows: list) -> Path:
    """Write ``results/sweep-table.csv``, one row per case."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(",".join(SWEEP_COLUMNS) + "\n")
        for row in rows:
            cells = []
            for column in SWEEP_COLUMNS:
                value = row.get(column)
                if value is None:
                    cells.append("")
                elif isinstance(value, float):
                    cells.append(f"{value:.6g}")
                else:
                    cells.append(str(value))
            handle.write(",".join(cells) + "\n")
    return path


def write_yaml(path: Path, document: dict, title: str) -> Path:
    """Write a generated YAML record with a do-not-edit header."""
    import yaml

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"# GENERATED by plumetools.cai2012.audit -- do not edit.\n"
        f"# {title}\n\n"
        + yaml.safe_dump(document, sort_keys=False, default_flow_style=False),
        encoding="utf-8", newline="\n")
    return path


__all__ = [
    "PLUME_FLOOR_OVER_N0",
    "SWEEP_COLUMNS",
    "FieldComparison",
    "RegionOccupancy",
    "common_times",
    "compare_field",
    "compare_overlap",
    "file_digest",
    "occupancy_audit",
    "occupancy_report",
    "overlap_report",
    "region_occupancy",
    "statistical_budget",
    "sweep_row",
    "weight_consistency",
    "write_sweep_table",
    "write_yaml",
]
