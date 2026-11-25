"""Decomposition weights, decomposeParDict, and the case scripts' calls into plateugf."""
from __future__ import annotations

import re
import shlex
from pathlib import Path

import numpy as np

from plumetools.foamio.checkpoints import END_DIVIDER
from plumetools.foamio.volfield import read_internal_field
from plumetools.plateugf import dictionaries, weights
from plumetools.plateugf.generation import SCRIPTS
from plumetools.plateugf.workflow import build_parser

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "cases" / "plate-unigas" / "baseCase"


def test_the_weight_field_reads_back_and_ends_like_an_openfoam_write(tmp_path):
    values = np.array([1.0, 2.5, 1234.5, 1.0])
    path = weights.write_field(tmp_path / weights.FIELD, values)
    assert np.allclose(np.asarray(read_internal_field(path)).ravel(), values)
    # the checkpoint scan calls a file without the divider "cut off"
    assert path.read_bytes().rstrip().endswith(END_DIVIDER)
    assert "weightField         decompositionWeight;" in dictionaries.render_decompose_par_dict(
        type("Cfg", (), {"dsmc": type("D", (), {"n_subdomains": 192,
                                                 "decomposition_method": "scotch"})()})())


def _plateugf_calls(text: str) -> list[list[str]]:
    """Every ``$PYTHON -m plumetools.plateugf ...`` command, continuation lines joined."""
    calls = []
    for line in text.replace("\\\n", " ").splitlines():
        if line.lstrip().startswith("#"):
            continue
        match = re.search(r"\$PYTHON -m plumetools\.plateugf (.*?)(\|\||&&|\)|$)", line)
        if match:
            tokens = shlex.split(match.group(1), comments=True)
            # "${ARRAY[@]}" carries optional flags (validated where it is built): drop it
            calls.append(["1" if t.startswith("$") else t for t in tokens if not t.endswith("[@]}")])
    return calls


def test_every_plateugf_call_in_the_case_scripts_parses():
    parser = build_parser()
    seen = 0
    for name in SCRIPTS:
        for tokens in _plateugf_calls((BASE / name).read_text(encoding="utf-8")):
            parser.parse_args(tokens)        # SystemExit on a stray or missing argument
            seen += 1
    assert seen >= 8
