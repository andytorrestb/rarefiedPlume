"""Parse ``constant/polyMesh/boundary``."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from plumetools.mesh._parse import read_text, strip_comments

#: A patch entry: a name followed by a brace-delimited body.
_PATCH = re.compile(r"(\w+)\s*\{([^{}]*)\}", re.S)


@dataclass(frozen=True)
class PatchInfo:
    """One entry of ``constant/polyMesh/boundary``.

    Attributes:
        name: patch name as written in the mesh.
        type: OpenFOAM patch type. May be ``"Unspecified"`` -- Pointwise V18.5R2
            emits that for unassigned patches, and three cases in this repository
            still carry it (finding HA-07). Parsing tolerates it; whether OpenFOAM
            loads it is a separate, unresolved question.
        n_faces: number of faces in the patch.
        start_face: index of the patch's first face in ``constant/polyMesh/faces``.
    """

    name: str
    type: str
    n_faces: int
    start_face: int


def read_boundary(case_dir: Path) -> dict[str, PatchInfo]:
    """Read every patch from ``<case_dir>/constant/polyMesh/boundary``.

    Args:
        case_dir: an OpenFOAM case directory.

    Returns:
        Mapping of patch name to :class:`PatchInfo`, in file order.

    Raises:
        FileNotFoundError: if the boundary file is missing.
        ValueError: if no patch entries could be parsed.

    Parsing is brace-structured rather than line-offset based. The legacy reader
    searched for the patch name as a substring and then skipped exactly three
    lines, which assumed the key order ``{`` / ``type`` / ``nFaces`` and broke on
    any file containing ``inGroups`` -- which is exactly what ``blockMesh`` writes
    (findings RB-03, RB-04).
    """
    path = Path(case_dir) / "constant" / "polyMesh" / "boundary"
    text = strip_comments(read_text(path))

    patches: dict[str, PatchInfo] = {}
    for match in _PATCH.finditer(text):
        name, body = match.group(1), match.group(2)
        n_faces = re.search(r"\bnFaces\s+(\d+)\s*;", body)
        start_face = re.search(r"\bstartFace\s+(\d+)\s*;", body)
        if not (n_faces and start_face):
            continue  # the FoamFile header block, not a patch
        patch_type = re.search(r"\btype\s+(\w+)\s*;", body)
        patches[name] = PatchInfo(
            name=name,
            type=patch_type.group(1) if patch_type else "Unspecified",
            n_faces=int(n_faces.group(1)),
            start_face=int(start_face.group(1)),
        )

    if not patches:
        raise ValueError(f"{path}: no patch entries with nFaces and startFace found")
    return patches


def require_patch(patches: dict[str, PatchInfo], name: str) -> PatchInfo:
    """Look up a patch by exact name, or raise listing what is available.

    Exact lookup, not substring matching: the legacy ``while patch not in line``
    would match ``inflow`` against a patch named ``inflowOuter`` (finding RB-03).
    """
    try:
        return patches[name]
    except KeyError:
        raise KeyError(
            f"patch {name!r} not found; mesh has {sorted(patches)}"
        ) from None
