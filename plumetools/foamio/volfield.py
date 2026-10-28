r"""Read values out of an OpenFOAM field file: the internal field, or one patch.

One parser for every family. It replaces three that had grown apart -- in
``nozzleplume.post``, ``cai2012.post`` and ``markelov1999.foamfields`` -- and
reads the union of what they read:

==============================================  ==============================
entry                                           :attr:`FieldValues.values`
==============================================  ==============================
``uniform 1.5;``                                ``(1,)``
``uniform (1 2 3);``                            ``(1, 3)``
``nonuniform List<scalar> N ( ... );``          ``(N,)``
``nonuniform List<vector> N ( (x y z) ... );``  ``(N, 3)``
``nonuniform List<symmTensor>`` / ``<tensor>``  ``(N, 6)`` / ``(N, 9)``
``nonuniform List<scalar> 3(1 2 3);``           the inline form OpenFOAM
                                                writes for ten entries or fewer
``nonuniform List<scalar> 0()`` or ``0;``       ``(0,)``: a patch with no faces
==============================================  ==============================

Shapes: a scalar field is ``(n,)``; every other type carries a component axis,
``(n, k)``. A ``uniform`` entry holds its one value, ``n = 1``, and says so in
:attr:`FieldValues.uniform` -- a caller that needs one value per cell expands it
with :meth:`FieldValues.over`, which is the only place a value is broadcast.

Finding the entry
-----------------
``internalField`` and ``boundaryField`` are found where a keyword can start --
the beginning of a line, or after a ``;`` -- so a header comment that mentions
either does not count. A patch is found by name inside ``boundaryField``, its
dictionary by balanced braces (so a nested sub-dictionary does not end it
early), and its ``value`` keyword at the top level of that dictionary. A list
payload holds no ``;``, so the first one after the entry closes it and the last
``)`` before that closes the list; neither search walks the payload in Python,
which on a field of a few million cells is the difference between milliseconds
and most of a minute.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from plumetools.foamio.errors import FoamDataError

#: Components per value, by OpenFOAM primitive type.
COMPONENTS = {"scalar": 1, "sphericalTensor": 1, "vector": 3, "symmTensor": 6,
              "tensor": 9}

#: The type a parenthesised ``uniform`` value of that many components has.
_UNIFORM_TYPES = {1: "sphericalTensor", 3: "vector", 6: "symmTensor", 9: "tensor"}

_HEAD = re.compile(r"\s*(uniform|nonuniform)\b")
_LIST_TYPE = re.compile(r"\s*List\s*<\s*(\w+)\s*>")
_COUNT = re.compile(r"\s*(\d+)")
_INTERNAL = re.compile(r"(?:^|;)[ \t]*internalField\b", re.M)
_BOUNDARY = re.compile(r"(?:^|[;{}])[ \t]*boundaryField\b", re.M)
_BRACES = re.compile(r"[{}]")
_VALUE = re.compile(r"^[ \t]*value\b", re.M)
_NO_PARENS = str.maketrans("()", "  ")


@dataclass(frozen=True)
class FieldValues:
    """One field entry -- an ``internalField`` or a patch ``value`` -- as written.

    Attributes:
        values: ``(n,)`` for a scalar, ``(n, k)`` for a type with components
            (vector 3, symmTensor 6, tensor 9, sphericalTensor 1). A
            ``uniform`` entry holds its one value, ``n = 1``.
        uniform: ``True`` for a ``uniform`` entry -- one value standing for
            every cell or face -- and ``False`` for a list.
        value_type: ``scalar``, ``vector``, ... -- from ``List<...>``, or from
            the number of components of a ``uniform`` value.
        declared_count: the length the list declares, or ``None`` for a
            ``uniform`` entry or a list written without one.
        source: the file and entry, for messages.
    """

    values: np.ndarray
    uniform: bool
    value_type: str
    declared_count: int | None
    source: str

    @property
    def n_values(self) -> int:
        """Values held: 1 for a ``uniform`` entry, the list length otherwise."""
        return int(self.values.shape[0])

    @property
    def n_components(self) -> int:
        return COMPONENTS[self.value_type]


def _uniform(token: str, source: str) -> FieldValues:
    try:
        if token.startswith("(") and token.endswith(")"):
            numbers = [float(v) for v in token[1:-1].split()]
            value_type = _UNIFORM_TYPES.get(len(numbers))
            if value_type is None:
                raise ValueError(token)
            return FieldValues(np.asarray([numbers], dtype=np.float64), True,
                               value_type, None, source)
        return FieldValues(np.asarray([float(token)], dtype=np.float64), True,
                           "scalar", None, source)
    except ValueError:
        raise FoamDataError(f"{source}: malformed uniform value {token!r}") from None


def parse_value(text: str, start: int, source: str) -> FieldValues:
    """Parse the ``uniform`` or ``nonuniform List<>`` entry that follows ``start``.

    Args:
        text: the whole field file.
        start: an index just past the keyword (``internalField`` or ``value``).
        source: the file and entry, for messages.

    Raises:
        FoamDataError: for anything that is not a uniform value or a typed list.
    """
    head = _HEAD.match(text, start)
    if not head:
        raise FoamDataError(f"{source}: expected 'uniform' or 'nonuniform'")
    end = text.find(";", head.end())
    if end < 0:
        raise FoamDataError(f"{source}: unterminated entry")
    if head.group(1) == "uniform":
        return _uniform(text[head.end():end].strip(), source)

    kind = _LIST_TYPE.match(text, head.end())
    if not kind:
        raise FoamDataError(f"{source}: nonuniform value with no List<> type")
    value_type = kind.group(1)
    if value_type not in COMPONENTS:
        raise FoamDataError(
            f"{source}: List<{value_type}> is not a field type this reader knows "
            f"({', '.join(COMPONENTS)})")
    count = _COUNT.match(text, kind.end())
    declared = int(count.group(1)) if count else None

    opening = text.find("(", kind.end(), end)
    if opening < 0:
        # 'nonuniform List<scalar> 0;' -- an empty list, written with no payload.
        numbers = np.empty(0, dtype=np.float64)
    else:
        closing = text.rfind(")", opening, end)
        payload = text[opening + 1:closing]
        if value_type != "scalar":
            payload = payload.translate(_NO_PARENS)
        numbers = (np.fromstring(payload, dtype=np.float64, sep=" ")
                   if payload.strip() else np.empty(0, dtype=np.float64))

    if value_type == "scalar":
        return FieldValues(numbers, False, value_type, declared, source)
    k = COMPONENTS[value_type]
    if numbers.size % k:
        raise FoamDataError(
            f"{source}: {numbers.size} numbers is not a whole number of "
            f"{k}-component {value_type} values")
    return FieldValues(numbers.reshape(-1, k), False, value_type, declared, source)


def _read(path: Path) -> str:
    path = Path(path)
    if not path.is_file():
        raise FoamDataError(f"no field file at {path}")
    return path.read_text(encoding="utf-8", errors="replace")


def read_internal(path: Path) -> FieldValues:
    """The ``internalField`` of a volume field, with its metadata.

    Raises:
        FoamDataError: if the file is missing or has no readable internal field.
    """
    text = _read(path)
    match = _INTERNAL.search(text)
    if not match:
        raise FoamDataError(f"{path}: no internalField entry")
    return parse_value(text, match.end(), f"{path}: internalField")


def _top_level(text: str, start: int, stop: int) -> bool:
    """Whether ``stop`` sits at the top level of the dictionary opened at ``start``."""
    segment = text[start:stop]
    return segment.count("{") == segment.count("}")


def read_patch(path: Path, patch: str) -> FieldValues:
    """One patch's ``value`` entry, with its metadata.

    Raises:
        FoamDataError: if the file or the patch is missing, or the patch has no
            ``value`` entry -- which is what a ``zeroGradient`` boundary
            condition looks like on disk, and means the data was never written
            rather than that it was zero.
    """
    text = _read(path)
    boundary = _BOUNDARY.search(text)
    if not boundary:
        raise FoamDataError(f"{path}: no boundaryField block")
    entry = re.compile(rf"^\s*{re.escape(patch)}(?:\s*$|\s*(?=\{{))", re.M).search(
        text, boundary.end())
    if not entry:
        raise FoamDataError(f"{path}: patch {patch!r} is not in the boundaryField")

    body_start = text.find("{", entry.end())
    if body_start < 0:
        raise FoamDataError(f"{path}: patch {patch!r} has no dictionary")
    depth, body_end = 0, -1
    for brace in _BRACES.finditer(text, body_start):
        depth += 1 if brace.group() == "{" else -1
        if depth == 0:
            body_end = brace.start()
            break
    if body_end < 0:
        raise FoamDataError(f"{path}: patch {patch!r}: unbalanced braces")

    for value in _VALUE.finditer(text, body_start + 1, body_end):
        if _top_level(text, body_start + 1, value.start()):
            return parse_value(text, value.end(), f"{path}: patch {patch!r}")
    raise FoamDataError(
        f"{path}: patch {patch!r} has no 'value' entry. A zeroGradient "
        f"boundary condition writes none, which means the solver's data for "
        f"this patch was discarded at write time rather than being zero -- "
        f"fD and q must be 'calculated' on a wall.")


def read_internal_field(path: Path) -> np.ndarray:
    """The ``internalField`` values: ``(n,)`` scalar, ``(n, k)`` otherwise.

    A ``uniform`` entry comes back as one row; :func:`read_internal` says which
    it was.
    """
    return read_internal(path).values


def read_patch_field(path: Path, patch: str) -> np.ndarray:
    """One patch's values: ``(n,)`` scalar, ``(n, k)`` otherwise.

    A ``uniform`` entry comes back as one row; :func:`read_patch` says which it
    was.
    """
    return read_patch(path, patch).values
