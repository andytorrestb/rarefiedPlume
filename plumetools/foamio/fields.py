"""Write the ``0/`` inflow boundary field files read by dsmcFoam+.

Three files are produced, all consumed by ``dsmcFreeStreamInflowFieldPatch``:

    0/boundaryU                     volVectorField   [0 1 -1 0 0 0 0]
    0/boundaryT                     volVectorField   [0 0 0 1 0 0 0]
    0/boundaryNumberDensity_<sp>    volScalarField   [0 -3 0 0 0 0 0]

``boundaryT`` is a *vector* field: dsmcFoam+ takes per-component translational
temperature, and the model writes ``(T 0 0)`` with the y and z components zero.

The layout below reproduces the original byte-for-byte, including its
irregularities -- five-space indents inside ``FoamFile``, nine-space indents
inside the trailing patch blocks, and the extra space in ``type  calculated;``.
It is not tidied because the regression golden depends on it.

Unlike the original this uses ``f.write`` rather than rebinding ``sys.stdout``
(finding AD-05: an exception mid-block left the interpreter's stdout pointing at
a closed file), and writes ``newline="\\n"`` explicitly so output is identical on
Windows and Linux (the original emitted CRLF on Windows via ``print``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

BANNER = r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1706                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/"""

SEPARATOR = "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"


def _foam_file_header(field_class: str, object_name: str) -> list[str]:
    return [
        BANNER,
        "FoamFile",
        "{",
        "     version     2.0;",
        "     format      ascii;",
        f"     class       {field_class};",
        '     location    "0";',
        f"     object      {object_name};",
        "}",
        SEPARATOR,
    ]


def _patch_block(name: str, spec: dict, zero_value: str) -> list[str]:
    """Render one non-inflow patch entry.

    Args:
        name: patch name.
        spec: ``{"type": ..., "value": ...}``. A ``value`` key of ``"uniform"``
            renders ``value uniform <zero_value>;``.
        zero_value: ``"(0 0 0)"`` for vector fields, ``"0"`` for scalars.

    The two indentation styles are the original's, not a choice: patches with a
    value use ``type  calculated;`` with two spaces, those without use a single.
    """
    patch_type = spec.get("type", "zeroGradient")
    lines = [f"    {name}", "    {"]
    if "value" in spec:
        lines.append(f"         type  {patch_type};")
        lines.append(f"         value uniform {zero_value};")
    else:
        lines.append(f"         type {patch_type};")
    lines.append("    }")
    return lines


def _boundary_field(entry_lines: list[str], patches: dict, zero_value: str) -> list[str]:
    """Assemble the ``boundaryField`` block: the inflow entry then the rest.

    Patch blocks after the first are separated by a blank line; there is no blank
    between the inflow entry and the first trailing patch. That asymmetry is the
    original's.
    """
    lines = ["boundaryField", "{", *entry_lines]
    for i, (name, spec) in enumerate(patches.items()):
        if i:
            lines.append("")
        lines.extend(_patch_block(name, spec, zero_value))
    lines.append("}")
    return lines


def _inflow_entry(list_type: str, values: Iterable[str]) -> list[str]:
    values = list(values)
    return [
        "    inflow",
        "    {",
        "        type            fixedValue;",
        f"        value           nonuniform List<{list_type}>",
        f"        {len(values)}",
        "        (",
        *values,
        "        );",
        "    }",
    ]


def _write(path: Path, lines: list[str]) -> None:
    """Write LF-terminated lines, with a trailing newline."""
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


def write_vol_vector_field(path: Path, object_name: str, dimensions: str,
                           rows: Iterable[str], patches: dict) -> None:
    """Write a ``volVectorField`` with a nonuniform inflow list.

    Args:
        path: destination file.
        object_name: the ``object`` entry, matching the file name.
        dimensions: e.g. ``"[0 1 -1 0 0 0 0]"``.
        rows: pre-formatted, pre-indented value lines.
        patches: non-inflow patch specs, in output order.
    """
    _write(Path(path), [
        *_foam_file_header("volVectorField", object_name),
        f"dimensions      {dimensions};",
        "internalField   uniform (0 0 0);",
        *_boundary_field(_inflow_entry("vector", rows), patches, "(0 0 0)"),
    ])


def write_vol_scalar_field(path: Path, object_name: str, dimensions: str,
                           rows: Iterable[str], patches: dict) -> None:
    """Write a ``volScalarField`` with a nonuniform inflow list. See above."""
    _write(Path(path), [
        *_foam_file_header("volScalarField", object_name),
        f"dimensions      {dimensions};",
        "internalField   uniform 0;",
        *_boundary_field(_inflow_entry("scalar", rows), patches, "0"),
    ])


def write_inflow_fields(out_dir: Path, inflow, cfg) -> None:
    """Write all three inflow field files for a case.

    Args:
        out_dir: the case's ``0/`` directory. Created if absent.
        inflow: an :class:`plumetools.inflow.InflowResult`.
        cfg: the case :class:`plumetools.config.CaseConfig`.

    Value formatting matches the original, which interpolated ``str()`` of each
    number: ``repr(float(x))`` for velocity and density (identical to
    ``str(numpy.float64)``), and the temperature rendered without a decimal point
    when it is integral -- ``calculateT`` returned the Python int ``300``, so the
    golden reads ``( 300 0.0 0.0 )``.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    patches = dict(cfg.output.patches)
    species = cfg.gas.species_name

    dialect = getattr(cfg.output, "dialect", "mnf")

    u_rows = [
        f"            ({repr(float(u))} {repr(float(v))} {repr(float(w))})"
        for u, v, w in inflow.U
    ]
    write_vol_vector_field(
        out_dir / "boundaryU", "boundaryU", "[0 1 -1 0 0 0 0]", u_rows, patches,
    )

    # boundaryT's TYPE differs between the two solvers, and getting it wrong is a
    # hard read failure. Verified against OpenFOAM v2512: FreeStream.C declares
    #     const volScalarField::Boundary& boundaryT = cloud.boundaryT()...
    # whereas the MNF fork takes a vector holding per-component translational
    # temperature, which is what the pre-refactor code wrote and what the
    # regression golden pins.
    if dialect == "standard":
        t_rows = [f"            {_fmt_temperature(t)}" for t in inflow.T]
        write_vol_scalar_field(
            out_dir / "boundaryT", "boundaryT", "[0 0 0 1 0 0 0]", t_rows, patches,
        )
    else:
        t_rows = [f"            ( {_fmt_temperature(t)} 0.0 0.0 )" for t in inflow.T]
        write_vol_vector_field(
            out_dir / "boundaryT", "boundaryT", "[0 0 0 1 0 0 0]", t_rows, patches,
        )

    # Per-face number density. The MNF fork's dsmcFreeStreamInflowFieldPatch reads
    # this; standard dsmcFoam has no equivalent and takes ONE scalar per species
    # from constant/dsmcProperties (FreeStream.C:93,
    # `numberDensities_[i] = numberDensitiesDict.get<scalar>(molecules[i])`).
    #
    # It is still written under the standard dialect: unused files in 0/ are
    # ignored, and it is the only record of what the model actually computed.
    # The uniform value standard dsmcFoam *will* use has to be put in
    # dsmcProperties by hand -- see area_weighted_number_density().
    n_rows = [f"              {repr(float(n))}" for n in inflow.rhoN]
    write_vol_scalar_field(
        out_dir / f"boundaryNumberDensity_{species}",
        f"boundaryNumberDensity_{species}",
        "[0 -3 0 0 0 0 0]", n_rows, patches,
    )


#: The measurement fields standard dsmcFoam requires in 0/ but dsmcInitialise
#: does not create. Specs read from the v2512 freeSpaceStream tutorial's 0.orig.
#: (name, class, dimensions, uniform internal value)
MEASUREMENT_FIELDS = (
    ("dsmcRhoN",  "volScalarField", "[0 -3 0 0 0 0 0]",  "0"),
    ("fD",        "volVectorField", "[1 -1 -2 0 0 0 0]", "(0 0 0)"),
    ("iDof",      "volScalarField", "[0 -3 0 0 0 0 0]",  "0"),
    ("internalE", "volScalarField", "[1 -1 -2 0 0 0 0]", "0"),
    ("linearKE",  "volScalarField", "[1 -1 -2 0 0 0 0]", "0"),
    ("momentum",  "volVectorField", "[1 -2 -1 0 0 0 0]", "(0 0 0)"),
    ("q",         "volScalarField", "[1 0 -3 0 0 0 0]",  "0"),
    ("rhoM",      "volScalarField", "[1 -3 0 0 0 0 0]",  "0"),
    ("rhoN",      "volScalarField", "[0 -3 0 0 0 0 0]",  "0"),
)


#: Geometric patch types whose patchField type must match exactly. OpenFOAM
#: rejects anything else with "inconsistent patch and patchField types".
CONSTRAINT_PATCH_TYPES = frozenset({
    "symmetry", "symmetryPlane", "empty", "wedge", "cyclic", "cyclicAMI",
    "processor", "processorCyclic", "nonuniformTransformCyclic",
})


def patch_field_type(geometric_type: str) -> str:
    """The patchField type a given geometric patch type requires.

    Constraint patches -- ``symmetry``, ``empty``, ``wedge``, ``cyclic`` and
    friends -- must carry a patchField of the same name; OpenFOAM fails with
    "inconsistent patch and patchField types" otherwise. Ordinary ``patch`` and
    ``wall`` boundaries take a normal condition, and for zeroed measurement
    fields that is ``zeroGradient``.
    """
    return geometric_type if geometric_type in CONSTRAINT_PATCH_TYPES else "zeroGradient"


def write_measurement_fields(out_dir: Path, patches) -> list[Path]:
    """Write the zeroed measurement fields standard dsmcFoam expects in ``0/``.

    Args:
        out_dir: the case's ``0/`` directory.
        patches: the mesh's patches -- a mapping of name to
            :class:`plumetools.mesh.boundary.PatchInfo`, as
            :func:`plumetools.mesh.read_boundary` returns. The geometric type of
            each is needed, not just its name.

    Returns:
        The paths written.

    ``dsmcFoam`` constructs its cloud from these and fails hard if any is absent
    -- ``cannot find file "0/q"`` and so on. ``dsmcInitialise`` does not create
    them: OpenFOAM's own tutorials ship them in ``0.orig`` and copy them in with
    ``restore0Dir``.

    Generated rather than shipped as a static ``0.orig``, because the boundary
    entries have to name the patches this mesh actually has. The pre-refactor code
    hardcoded a patch list and wrote ``cylinder`` and ``plate`` into every field on
    meshes that had neither (finding AD-03); a checked-in ``0.orig`` would repeat
    that mistake for a lineage whose cases have inflow/sym/vacuum,
    inflow/vacuum, and inflow/panel/vacuum respectively.

    All nine are ``uniform 0`` internally -- they are outputs the solver fills in,
    not inputs. Boundary conditions are ``zeroGradient`` except on constraint
    patches, which must repeat their own type; see :func:`patch_field_type`.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Accept either {name: PatchInfo} or a bare sequence of names, treating the
    # latter as ordinary patches.
    if hasattr(patches, "items"):
        entries = [(n, patch_field_type(getattr(p, "type", "patch")))
                   for n, p in patches.items()]
    else:
        entries = [(n, "zeroGradient") for n in patches]

    written = []
    for name, field_class, dimensions, zero in MEASUREMENT_FIELDS:
        lines = [
            *_foam_file_header(field_class, name),
            f"dimensions      {dimensions};",
            f"internalField   uniform {zero};",
            "boundaryField",
            "{",
        ]
        for patch_name, bc in entries:
            lines += [f"    {patch_name}", "    {", f"        type            {bc};", "    }"]
        lines.append("}")
        path = out_dir / name
        _write(path, lines)
        written.append(path)
    return written


def _fmt_temperature(value: float) -> str:
    """Render temperature the way the original did.

    ``calculateT`` returned a Python ``int``, so ``str()`` gave ``300``; a config
    carrying ``T0_K: 300.0`` would give ``300.0``. Integral values are rendered
    without the decimal point so the generated file matches the golden byte for
    byte. Non-integral temperatures fall through to the normal float repr.
    """
    value = float(value)
    return str(int(value)) if value.is_integer() else repr(value)
