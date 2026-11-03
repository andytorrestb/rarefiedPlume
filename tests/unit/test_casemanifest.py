"""The runner scripts' case list: plumetools.casemanifest, and the scripts using it.

The defect this pins: the cai2012, cai2012-health and markelov1999 runners read
``manifest.yaml`` with inline Python inside ``mapfile -t CASES < <(...) || die``.
A process substitution's exit status never reaches the ``||``, and the paths
were printed as they were read, so a manifest that broke at its third entry ran
the first two cases and the script reported success.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from plumetools import casemanifest

REPO = Path(__file__).resolve().parents[2]

RUNNERS = [REPO / "cases" / family / script
           for family in ("cai2012", "cai2012-health", "markelov1999")
           for script in ("AllmeshCases", "AllrunCases", "AllpostCases", "AllcleanCases")]

BROKEN = {"cases": [{"path": "Cases/a"}, {"path": "Cases/b"}, {"name": "no path"},
                    {"path": "Cases/d"}]}


def write(directory: Path, document) -> Path:
    path = directory / "manifest.yaml"
    path.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    return path


def test_every_path_is_listed_in_order(tmp_path):
    path = write(tmp_path, {"cases": [{"path": "Cases/b"}, {"path": "Cases/a"}]})
    assert casemanifest.load(path) == ["Cases/b", "Cases/a"]


@pytest.mark.parametrize("document,message", [
    (BROKEN, r"cases\[2\] has no usable 'path'"),
    ({"cases": [{"path": "Cases/a"}, {"path": "Cases/a"}]}, "listed twice"),
    ({"cases": {"path": "x"}}, "'cases' is not a list"),
    (["Cases/a"], "not a mapping"),
    ({"cases": [{"path": ""}]}, "no usable 'path'"),
])
def test_a_bad_manifest_lists_nothing(tmp_path, capsys, document, message):
    path = write(tmp_path, document)
    with pytest.raises(casemanifest.ManifestError, match=message):
        casemanifest.load(path)
    assert casemanifest.main([str(path)]) == 1
    captured = capsys.readouterr()
    assert captured.out == "" and "casemanifest:" in captured.err


def test_only_keeps_manifest_order_and_refuses_a_typo(tmp_path, capsys):
    path = write(tmp_path, {"cases": [{"path": "Cases/a"}, {"path": "Cases/b"},
                                      {"path": "Cases/c"}]})
    assert casemanifest.main([str(path), "--only", "Cases/c", "Cases/a"]) == 0
    assert capsys.readouterr().out.split() == ["Cases/a", "Cases/c"]
    with pytest.raises(casemanifest.ManifestError, match="--only names no generated case"):
        casemanifest.select(["Cases/a"], ["Cases/z"])


def test_an_empty_or_missing_manifest(tmp_path):
    assert casemanifest.load(write(tmp_path, {"cases": []})) == []
    with pytest.raises(casemanifest.ManifestError):
        casemanifest.load(tmp_path / "absent.yaml")


@pytest.mark.parametrize("script", RUNNERS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_no_runner_reads_the_manifest_inside_a_process_substitution(script):
    text = script.read_text(encoding="utf-8")
    assert "plumetools.casemanifest" in text
    assert "< <($PYTHON" not in text


# --------------------------------------------------------------------------- #
# the scripts themselves (bash)
# --------------------------------------------------------------------------- #

needs_bash = pytest.mark.skipif(
    sys.platform == "win32" or shutil.which("bash") is None,
    reason="bash on native Windows is WSL's and cannot open a Windows path; "
           "the suite runs under WSL")


def python_env(tmp_path: Path) -> dict:
    """PATH with a ``python`` that is this interpreter, and plumetools on it."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    wrapper = bin_dir / "python"
    wrapper.write_text(f'#!/bin/sh\nexec "{Path(sys.executable).as_posix()}" "$@"\n',
                       encoding="utf-8")
    wrapper.chmod(0o755)
    env = dict(os.environ)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [str(REPO), env.get("PYTHONPATH")]))
    return env


OLD_READER = '''
die() { echo "die: $*" >&2; exit 1; }
PYTHON=python
mapfile -t CASES < <($PYTHON - <<'PY'
import yaml
document = yaml.safe_load(open("manifest.yaml", encoding="utf-8"))
for entry in document.get("cases", []):
    print(entry["path"])
PY
) || die "could not read manifest.yaml"
echo "would run ${#CASES[@]} case(s): ${CASES[*]}"
'''


@needs_bash
def test_the_old_reader_ran_the_head_of_a_broken_manifest(tmp_path):
    """The trigger, reproduced with the snippet every runner carried."""
    write(tmp_path, BROKEN)
    result = subprocess.run(["bash", "-c", OLD_READER], cwd=tmp_path, capture_output=True,
                            text=True, env=python_env(tmp_path), timeout=120)
    assert result.returncode == 0
    assert "would run 2 case(s): Cases/a Cases/b" in result.stdout


@needs_bash
@pytest.mark.parametrize("family,args", [("cai2012", ["--list"]), ("markelov1999", [])])
def test_a_runner_refuses_a_broken_manifest_before_any_case(tmp_path, family, args):
    shutil.copy2(REPO / "cases" / family / "AllmeshCases", tmp_path / "AllmeshCases")
    write(tmp_path, BROKEN)
    result = subprocess.run(["bash", "./AllmeshCases", *args], cwd=tmp_path,
                            capture_output=True, text=True, env=python_env(tmp_path),
                            timeout=120)
    assert result.returncode != 0
    assert "Cases/a" not in result.stdout and "could not read manifest.yaml" in result.stderr


@needs_bash
def test_a_runner_lists_a_good_manifest(tmp_path):
    shutil.copy2(REPO / "cases" / "cai2012" / "AllmeshCases", tmp_path / "AllmeshCases")
    write(tmp_path, {"cases": [{"path": "Cases/a"}, {"path": "Cases/b"}]})
    result = subprocess.run(["bash", "./AllmeshCases", "--list"], cwd=tmp_path,
                            capture_output=True, text=True, env=python_env(tmp_path),
                            timeout=120)
    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == ["0", "Cases/a", "1", "Cases/b"]
