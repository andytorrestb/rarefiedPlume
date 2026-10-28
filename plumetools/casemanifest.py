"""The case list a study's runner scripts work through: read whole, checked, then listed.

``AllmeshCases``, ``AllrunCases``, ``AllpostCases`` and ``AllcleanCases`` of
the cai2012, cai2012-health and markelov1999 families take their cases from
``manifest.yaml`` -- the record of what was generated -- rather than from a
glob. They used to read it with a few lines of inline Python inside
``mapfile < <(...)``, which has two holes: the exit status of a process
substitution never reaches the ``|| die`` after it, and the paths were printed
as they were read, so a manifest that broke at its fourth entry launched the
first three cases and reported success.

Here the whole manifest is read and checked before a single path is printed,
and a problem is an exit status the script can see::

    LISTING=$($PYTHON -m plumetools.casemanifest manifest.yaml) || die "..."

.. code-block:: text

    python -m plumetools.casemanifest [manifest.yaml] [--only PATH ...]

Exit status 0 with one path per line, in manifest order; 1 with the reason on
stderr and nothing on stdout. (The Kn extension's scripts list through
``plumetools.cai2012ext.policy``, which applies the study's Kn policy too.)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml


class ManifestError(ValueError):
    """The manifest cannot be read, or does not list cases a runner can use."""


def case_paths(document, where: str = "manifest.yaml") -> list[str]:
    """Every case path in a parsed manifest, in its order, all checked.

    Raises:
        ManifestError: if the document is not a mapping with a ``cases`` list,
            an entry is not a mapping with a non-empty single-line ``path``, or
            a path appears twice.
    """
    if not isinstance(document, dict):
        raise ManifestError(f"{where}: not a mapping")
    cases = document.get("cases", [])
    if not isinstance(cases, list):
        raise ManifestError(f"{where}: 'cases' is not a list")
    paths = []
    for index, entry in enumerate(cases):
        path = entry.get("path") if isinstance(entry, dict) else None
        if not isinstance(path, str) or not path.strip() or "\n" in path:
            raise ManifestError(
                f"{where}: cases[{index}] has no usable 'path' (got {entry!r}); "
                f"no case was listed")
        paths.append(path)
    duplicates = sorted({p for p in paths if paths.count(p) > 1})
    if duplicates:
        raise ManifestError(f"{where}: case path(s) listed twice: {duplicates}")
    return paths


def select(paths: list[str], only) -> list[str]:
    """The paths in ``only``, in manifest order; every one must be listed.

    Raises:
        ManifestError: naming the ones that are not -- a mistyped ``--only``
            that silently ran nothing would look exactly like a study that
            finished instantly.
    """
    wanted = set(only)
    missing = sorted(wanted - set(paths))
    if missing:
        raise ManifestError(f"--only names no generated case: {missing}\n"
                            f"the manifest has: {paths}")
    return [p for p in paths if p in wanted]


def load(path: Path) -> list[str]:
    """:func:`case_paths` of the manifest file at ``path``."""
    path = Path(path)
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise ManifestError(f"{path}: {exc}") from None
    return case_paths(document, str(path))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m plumetools.casemanifest",
        description="List a study manifest's case paths, after checking all of them.")
    parser.add_argument("manifest", nargs="?", type=Path, default=Path("manifest.yaml"))
    parser.add_argument("--only", nargs="+", default=None, metavar="PATH",
                        help="list only these case paths (each must be in the manifest)")
    args = parser.parse_args(argv)
    try:
        paths = load(args.manifest)
        if args.only:
            paths = select(paths, args.only)
    except ManifestError as exc:
        print(f"casemanifest: {exc}", file=sys.stderr)
        return 1
    if paths:
        print("\n".join(paths))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
