"""The post-processing chain of the cai2012-unigas family (spec §5, §7).

The analysis scripts are the kn-extension's, byte for byte, except that
Allcompare.py and Allconverge.py import this family's post context (the
post-run verdict) in place of the kn-extension's pre-run one -- one line.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

from plumetools.cai2012ugf.generation import COPIED_FROM_KN_EXTENSION, copy_differences

REPO = Path(__file__).resolve().parents[2]
UGF = REPO / "cases" / "cai2012-unigas"


def copy_family(tmp_path: Path) -> Path:
    target = tmp_path / "cai2012-unigas"
    shutil.copytree(UGF, target, ignore=shutil.ignore_patterns(
        "Cases", "Pilot", "results", "manifest-*.yaml", "__pycache__"))
    return target


def read_steps(path: Path) -> list[str]:
    """The commands Allpostugf runs, in order: script names and module names."""
    steps = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith("#"):
            continue
        match = re.search(r"plumetools\.cai2012ugf\.(\w+)|\b(Allpost|All\w+\.py)\b", line)
        if match and re.match(r"^(step|\$PYTHON|bash)\b", line):
            steps.append(match.group(1) or match.group(2))
    return steps


def test_copies_match_originals_except_the_context_import():
    assert copy_differences(UGF) == []


def test_every_copied_file_exists():
    for relative in COPIED_FROM_KN_EXTENSION:
        assert (UGF / relative).is_file(), relative


def test_a_stray_edit_in_a_copy_is_named(tmp_path):
    fam = copy_family(tmp_path)
    with open(fam / "baseCase/Allcompare.py", "a", encoding="utf-8") as handle:
        handle.write("# edit\n")
    assert any("Allcompare.py" in d for d in copy_differences(fam))


def test_a_byte_copy_edit_is_named(tmp_path):
    fam = copy_family(tmp_path)
    with open(fam / "baseCase/Allaudit.py", "a", encoding="utf-8") as handle:
        handle.write("\n")
    assert any("Allaudit.py" in d for d in copy_differences(fam))


def test_the_context_import_is_this_familys():
    text = (UGF / "baseCase/Allcompare.py").read_text(encoding="utf-8")
    assert "from plumetools.cai2012ugf import context as interface_case" in text
    assert "from plumetools.cai2012ext import figures, interface_case" not in text


def test_cases_listing_reads_the_manifest(tmp_path, capsys):
    from plumetools.cai2012ugf.generation import main
    (tmp_path / "manifest-unigas.yaml").write_text(
        "cases:\n- {name: Kn1_np1x, path: Cases/Kn1_np1x}\n- {name: Kn10_np1x, path: Cases/Kn10_np1x}\n")
    assert main(["cases", str(tmp_path)]) == 0
    assert capsys.readouterr().out.split() == ["Cases/Kn1_np1x", "Cases/Kn10_np1x"]


def test_cases_listing_without_a_manifest_fails(tmp_path):
    from plumetools.cai2012ugf.generation import main
    assert main(["cases", str(tmp_path)]) == 1


def test_post_order():
    steps = read_steps(UGF / "baseCase/Allpostugf")
    assert steps.index("convert") < steps.index("Allpost") < steps.index("Allconverge.py") \
        < steps.index("Allverdict.py") < steps.index("Allcompare.py")
    assert "Alldiagnose.py" in steps and "Allaudit.py" in steps
