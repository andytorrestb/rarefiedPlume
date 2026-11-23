#!/usr/bin/env python3
"""Check that the generated cases are the kn-extension's cases, on uniGasFoam.

    ./validate_cases.py

Read-only. Fails (exit 1) and names every difference when:

* baseCase/case.yaml differs from the kn-extension's outside meta;
* study.yaml's case matrix, overrides or particle level differ;
* a case's blockMeshDict, topoSetDict, createPatchDict, meshQualityDict or
  decomposeParDict is not byte-identical to the kn-extension's for that Kn;
* a case's particle weight differs from the kn-extension's (1e-12 relative);
* a case's unigas.yaml is invalid;
* a copied post-processing script differs from its kn-extension original
  beyond the one allowed import line.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    from plumetools.cai2012ugf import generation
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')


def main() -> int:
    problems = generation.validate(HERE)
    if problems:
        print(f"validate_cases: {len(problems)} problem(s):", file=sys.stderr)
        for problem in problems:
            print(f"    {problem}", file=sys.stderr)
        return 1
    print("validate_cases: the family matches the kn-extension case for case.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
