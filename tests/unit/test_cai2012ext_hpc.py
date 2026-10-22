"""What a case's own PBS job needs from the extension's drivers.

``AllpostCases --no-analyse`` runs one case's post-processing chain without
the study-level analysis, which waits for every case. It is run here against a
stand-in ``analyse.py``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
FAMILY = REPO / "cases" / "cai2012-kn-extension"

needs_bash = pytest.mark.skipif(shutil.which("bash") is None, reason="no bash")


def python_env(tmp_path: Path) -> dict:
    """PATH with a ``python`` that is this interpreter, and plumetools on it."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    wrapper = bin_dir / "python"
    wrapper.write_text(f'#!/bin/sh\nexec "{sys.executable}" "$@"\n',
                       encoding="utf-8")
    wrapper.chmod(0o755)
    env = {key: value for key, value in os.environ.items()
           if not key.startswith("PBS_")}
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    env["PYTHONPATH"] = os.pathsep.join(
        filter(None, [str(REPO), env.get("PYTHONPATH")]))
    return env


# --------------------------------------------------------------------------- #
# AllpostCases --no-analyse
# --------------------------------------------------------------------------- #

@needs_bash
def test_allpostcases_no_analyse_skips_the_study_analysis(tmp_path):
    family = tmp_path / "fam"
    (family / "Cases" / "Kn1_np1x").mkdir(parents=True)
    shutil.copy2(FAMILY / "AllpostCases", family / "AllpostCases")
    (family / "manifest.yaml").write_text(yaml.safe_dump({"cases": [
        {"name": "Kn1_np1x", "path": "Cases/Kn1_np1x", "Kn": 1.0,
         "feasibility": {"status": "production-candidate"}}]}), encoding="utf-8")
    ran = family / "analysed"
    (family / "analyse.py").write_text(
        "from pathlib import Path\nPath('analysed').write_text('yes')\n",
        encoding="utf-8")
    env = python_env(tmp_path)
    run = lambda *flags: subprocess.run(  # noqa: E731
        ["bash", str(family / "AllpostCases"), *flags], capture_output=True,
        text=True, cwd=tmp_path, env=env, timeout=300)

    result = run("--case", "Kn1_np1x", "--no-viz", "--no-analyse")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "no solver output yet" in result.stdout
    assert not ran.exists()

    result = run("--case", "Kn1_np1x", "--no-viz")
    assert result.returncode == 0, result.stdout + result.stderr
    assert ran.is_file()

    result = run("--analyse-only", "--no-analyse")
    assert result.returncode == 1
    assert "opposite" in result.stderr
