"""Parse ``constant/polyMesh/sets/<name>`` face sets."""

from __future__ import annotations

from pathlib import Path

from plumetools.mesh._parse import find_list_block, read_text, strip_comments


def read_face_set(case_dir: Path, name: str) -> list[int]:
    """Read the global face labels of a ``faceSet``.

    Args:
        case_dir: an OpenFOAM case directory.
        name: the set name, e.g. ``"inflow"`` -- created by ``topoSet`` from
            ``system/topoSetDict``.

    Returns:
        Global face labels **in file order**. Order is part of the contract with
        dsmcFoam+: the values written to ``0/boundary*`` are consumed in patch-face
        order, so this must not be sorted or de-duplicated.

    Raises:
        FileNotFoundError: if the set file is missing (run ``topoSet`` first).
    """
    path = Path(case_dir) / "constant" / "polyMesh" / "sets" / name
    if not path.is_file():
        raise FileNotFoundError(
            f"face set {name!r} not found at {path}; run `topoSet` to create it "
            f"from system/topoSetDict"
        )
    text = strip_comments(read_text(path))
    count, first = find_list_block(text, path)
    lines = text.split("\n")

    labels = []
    for i in range(count):
        raw = lines[first + i].strip()
        try:
            labels.append(int(raw))
        except ValueError as exc:
            raise ValueError(
                f"{path}: line {first + i + 1}: expected a face label, got {raw!r}"
            ) from exc
    return labels
