"""Shared helpers for reading OpenFOAM ASCII list files.

Every parser here is bounded and fails loudly. The legacy code searched with
``while token not in line: line = f.readline()``, which spins forever at EOF
because ``readline()`` keeps returning ``''`` -- a missing patch or a failed
``checkMesh`` hung the script with no output (finding RB-01).
"""

from __future__ import annotations

from pathlib import Path


def read_text(path: Path) -> str:
    """Read an OpenFOAM file, normalising line endings.

    Meshes exported on one platform and parsed on another otherwise differ by
    stray ``\\r`` characters, which broke the legacy fixed-width slicing (RB-02).
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"expected an OpenFOAM file at {path}")
    return path.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n").replace("\r", "\n")


def strip_comments(text: str) -> str:
    """Remove ``//`` line comments and ``/* */`` block comments."""
    out = []
    i, n = 0, len(text)
    while i < n:
        if text.startswith("//", i):
            j = text.find("\n", i)
            i = n if j < 0 else j
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            i = n if j < 0 else j + 2
        else:
            out.append(text[i])
            i += 1
    return "".join(out)


def find_list_block(text: str, path: Path) -> tuple[int, int]:
    """Locate an OpenFOAM ``<count> ( ... )`` list.

    Returns ``(count, index_of_first_entry_line)`` into ``text.splitlines()``.

    The count is taken from the last non-empty line before the opening ``(``,
    rather than by searching for a line that happens to contain the digits --
    the legacy ``while str(nFaces) not in line`` matched headers by accident
    (finding RB-03).
    """
    lines = text.split("\n")
    for i, line in enumerate(lines):
        if line.strip() == "(":
            for j in range(i - 1, -1, -1):
                stripped = lines[j].strip()
                if not stripped:
                    continue
                try:
                    return int(stripped), i + 1
                except ValueError as exc:
                    raise ValueError(
                        f"{path}: expected a list count before '(' on line {i + 1}, "
                        f"found {stripped!r}"
                    ) from exc
            break
    raise ValueError(f"{path}: no '<count> ( ... )' list block found")
