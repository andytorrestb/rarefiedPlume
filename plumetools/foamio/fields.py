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

    u_rows = [
        f"            ({repr(float(u))} {repr(float(v))} {repr(float(w))})"
        for u, v, w in inflow.U
    ]
    write_vol_vector_field(
        out_dir / "boundaryU", "boundaryU", "[0 1 -1 0 0 0 0]", u_rows, patches,
    )

    t_rows = [f"            ( {_fmt_temperature(t)} 0.0 0.0 )" for t in inflow.T]
    write_vol_vector_field(
        out_dir / "boundaryT", "boundaryT", "[0 0 0 1 0 0 0]", t_rows, patches,
    )

    n_rows = [f"              {repr(float(n))}" for n in inflow.rhoN]
    write_vol_scalar_field(
        out_dir / f"boundaryNumberDensity_{species}",
        f"boundaryNumberDensity_{species}",
        "[0 -3 0 0 0 0 0]", n_rows, patches,
    )


def _fmt_temperature(value: float) -> str:
    """Render temperature the way the original did.

    ``calculateT`` returned a Python ``int``, so ``str()`` gave ``300``; a config
    carrying ``T0_K: 300.0`` would give ``300.0``. Integral values are rendered
    without the decimal point so the generated file matches the golden byte for
    byte. Non-integral temperatures fall through to the normal float repr.
    """
    value = float(value)
    return str(int(value)) if value.is_integer() else repr(value)
