"""Run the documentation's allowlisted examples and preparation checks, safely.

    python docs/_tools/run_examples.py [--list] [--only NAME ...] [--keep]

Every check runs in a **temporary copy of the tracked files** it needs (from
``git ls-files``), never in the checkout, with a time limit and, on Linux, a
memory limit. After each one the copy is compared with its state before the
check: a check documented as read-only fails if it created or changed any file,
and a preparation check fails if it wrote outside the paths it is documented to
write. Nothing here meshes, solves, submits a job or needs OpenFOAM, ParaView,
MPI or CaseFoam; checks that need something absent are skipped with the reason.

Two kinds of check:

* **page examples** -- the fenced ``python`` block that follows a
  ``<!-- docs-example: NAME -->`` marker in a page. Only the names in
  ``EXAMPLES`` run; a marker not listed there is an error, so nothing executes
  by accident. Fenced blocks without a marker are never executed.
* **commands** -- documented preparation steps (``--dry-run`` generation,
  ``--describe``, ``setup``, ``plan`` on a synthetic case, …) listed in
  ``COMMANDS``.

``plumetools`` must be importable with its own dependencies (numpy, scipy,
PyYAML); the documentation environment alone is not enough. Exit status 0 when
every check that ran passed.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
TIMEOUT_S = 600
MEMORY_LIMIT_BYTES = 4 * 1024**3      # per check, Linux only

#: Page examples allowed to run, and the tracked subtrees each needs.
EXAMPLES = {
    "cai2012-derivation": ["cases/cai2012"],
    "cai2012-analytical": [],
    "first-run-analytical": [],
    "cai2012ext-feasibility": ["cases/cai2012-kn-extension", "cases/cai2012-health"],
    "markelov-sourceflow": [],
    "nozzleplume-load": ["cases/r-4d"],
    "hybrid-config": ["cases/cai2012-hybrid"],
    "rarefaction-mfp": [],
    "viz-spec": [],
    "foamio-checkpoints": [],
    "casemanifest": [],
}


@dataclass
class Command:
    """One documented command, run in a copy of ``needs``."""
    name: str
    argv: list
    needs: list
    writes: tuple = ()               # path prefixes it may create or change
    expect_exit: int = 0
    expect_stdout: str | None = None  # a regular expression the output must match
    shell: bool = False               # needs bash
    before: list = field(default_factory=list)   # setup steps (argv lists) run first


PY = "{python}"
COMMANDS = [
    Command("cai2012 generate --dry-run",
            [PY, "cases/cai2012/generate_cases.py", "--dry-run"], ["cases/cai2012"],
            expect_stdout=r"Kn100_np1x"),
    Command("cai2012-health generate --dry-run (no CaseFoam)",
            [PY, "cases/cai2012-health/generate_cases.py", "--dry-run"], ["cases/cai2012-health"],
            expect_stdout=r"--dry-run: nothing written\."),
    Command("Kn extension generate --dry-run",
            [PY, "cases/cai2012-kn-extension/generate_cases.py", "--dry-run"],
            ["cases/cai2012-kn-extension", "cases/cai2012-health"],
            expect_stdout=r"--dry-run: nothing written\."),
    Command("Kn extension generate, then validate --inputs-only",
            [PY, "cases/cai2012-kn-extension/validate_cases.py", "--inputs-only"],
            ["cases/cai2012-kn-extension", "cases/cai2012-health"],
            before=[[PY, "cases/cai2012-kn-extension/generate_cases.py"]]),
    Command("Kn extension policy listing (--all)",
            [PY, "-m", "plumetools.cai2012ext.policy", "--study",
             "cases/cai2012-kn-extension/study.yaml", "--all"],
            ["cases/cai2012-kn-extension", "cases/cai2012-health"],
            before=[[PY, "cases/cai2012-kn-extension/generate_cases.py"]],
            expect_stdout=r"Kn0p0001_np1x"),
    Command("Kn extension smoke_cases.py (preparation: writes the tree and dictionaries)",
            [PY, "cases/cai2012-kn-extension/smoke_cases.py"],
            ["cases/cai2012-kn-extension", "cases/cai2012-health"],
            writes=("cases/cai2012-kn-extension/Cases_full_suite/",
                    "cases/cai2012-kn-extension/manifest-full-suite.yaml",
                    "cases/cai2012-kn-extension/results/full-suite-preparation/")),
    Command("Kn extension scaling smoke (synthetic, bounded)",
            [PY, "cases/cai2012-kn-extension/Allscale.py", "smoke", "--no-plots"],
            ["cases/cai2012-kn-extension", "cases/cai2012-health"],
            writes=("cases/cai2012-kn-extension/results/scaling/full-suite/",)),
    Command("interfaces generate --dry-run",
            [PY, "cases/cai2012-kn-extension/generate_interfaces.py", "--dry-run"],
            ["cases/cai2012-kn-extension", "cases/cai2012-health", "cases/cai2012"]),
    Command("nozzleplume describe (r-4d)",
            [PY, "-m", "plumetools.nozzleplume", "describe", "cases/r-4d/baseCase"],
            ["cases/r-4d"]),
    Command("hybrid setup into a new case directory",
            [PY, "-m", "plumetools.hybridplume", "setup", "cases/cai2012-hybrid/runs/doc-check",
             "--config", "cases/cai2012-hybrid/hybrid.yaml", "--profile", "smoke-kn0p01-coarse",
             "--np", "2"],
            ["cases/cai2012-hybrid"], writes=("cases/cai2012-hybrid/runs/doc-check/",),
            expect_stdout=r"smoke-kn0p01-coarse"),
    Command("iteration: a fresh r-4d iteration",
            [PY, "-m", "plumetools.iteration", "cases/r-4d/baseCase", "cases/r-4d/iter-doc-check",
             "--note", "documentation check"],
            ["cases/r-4d"], writes=("cases/r-4d/iter-doc-check/",)),
    Command("Kn extension scaling inventory",
            [PY, "cases/cai2012-kn-extension/Allscale.py", "inventory"],
            ["cases/cai2012-kn-extension", "cases/cai2012-health"],
            writes=("cases/cai2012-kn-extension/results/scaling/full-suite/",)),
    Command("virtual-collision inventory into a scratch directory",
            [PY, "cases/cai2012-virtual-collision/assess.py", "inventory",
             "--output", "scratch/vc-inventory"],
            ["cases/cai2012-virtual-collision", "cases/cai2012-kn-extension"],
            writes=("scratch/vc-inventory/",)),
    Command("3d-inflow Allmesh --dict-only (bash, no OpenFOAM)",
            ["bash", "cases/3d-inflow/Allmesh", "--dict-only"], ["cases/3d-inflow"],
            shell=True,
            writes=("cases/3d-inflow/system/blockMeshDict", "cases/3d-inflow/system/snappyHexMeshDict",
                    "cases/3d-inflow/system/meshQualityDict")),
    Command("r-4d Allmesh --dict-only (bash, no OpenFOAM)",
            ["bash", "cases/r-4d/baseCase/Allmesh", "--dict-only"], ["cases/r-4d"],
            shell=True, writes=("cases/r-4d/baseCase/system/", "cases/r-4d/baseCase/constant/")),
]

# Synthetic checkpoint fixtures for `checkpoints plan`: a two-rank case with one
# complete time, and the same with a newer time one rank wrote only partly.
DIVIDER = "// " + "*" * 73 + " //\n"


def _write_rank(root: Path, rank: int, time_name: str, complete: bool) -> None:
    d = root / f"processor{rank}" / time_name
    (d / "uniform").mkdir(parents=True, exist_ok=True)
    (d / "uniform" / "time").write_text("value 0;\n" + DIVIDER)
    (d / "rhoN").write_text("internalField uniform 0;\n" + (DIVIDER if complete else ""))


def _checkpoint_fixture(copy: Path, partial: bool) -> Path:
    case = copy / "scratch" / ("ckpt-partial" if partial else "ckpt-complete")
    (case / "system").mkdir(parents=True)
    (case / "system" / "controlDict").write_text("endTime 1;\ndeltaT 0.1;\n")
    for rank in (0, 1):
        _write_rank(case, rank, "0.1", True)
    _write_rank(case, 0, "0.2", True)
    _write_rank(case, 1, "0.2", not partial)
    return case


CHECKPOINT_CASES = [
    ("checkpoints plan: complete on every rank -> resume, reuse processors", False, 0,
     r"RESUME=0\.2\s+.*DECOMPOSE=none"),
    ("checkpoints plan: partial time on one rank -> refused (exit 3)", True, 3, None),
]


# --------------------------------------------------------------------------- #

def tracked(paths: list) -> list:
    out = subprocess.run(["git", "ls-files", "-z", "--", *paths], cwd=ROOT,
                         capture_output=True, check=True)
    return [p for p in out.stdout.decode("utf-8").split("\0") if p]


def make_copy(needs: list) -> Path:
    copy = Path(tempfile.mkdtemp(prefix="rp-docs-check-"))
    for rel in tracked(["plumetools", *needs]):
        source = ROOT / rel
        if source.is_symlink() or not source.is_file():
            continue
        target = copy / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    (copy / "scratch").mkdir()
    return copy


def snapshot(root: Path) -> dict:
    state = {}
    for path in root.rglob("*"):
        if path.is_file() and "__pycache__" not in path.parts:
            stat = path.stat()
            state[path.relative_to(root).as_posix()] = (stat.st_size, stat.st_mtime_ns)
    return state


def changes(before: dict, after: dict) -> list:
    return sorted({p for p in set(before) | set(after) if before.get(p) != after.get(p)})


def _limits():
    try:
        import resource
        resource.setrlimit(resource.RLIMIT_AS, (MEMORY_LIMIT_BYTES, MEMORY_LIMIT_BYTES))
    except (ImportError, ValueError, OSError):
        pass


def run(argv, cwd: Path, timeout=TIMEOUT_S):
    env = dict(os.environ, PYTHONPATH=str(cwd), MPLBACKEND="Agg", PYTHONDONTWRITEBYTECODE="1")
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[key] = "1"
    argv = [sys.executable if a == PY else a for a in argv]
    kwargs = {"preexec_fn": _limits} if os.name == "posix" else {}
    return subprocess.run(argv, cwd=cwd, env=env, capture_output=True, text=True,
                          timeout=timeout, **kwargs)


def page_examples() -> dict:
    marker = re.compile(r"<!--\s*docs-example:\s*([\w-]+)\s*-->")
    found = {}
    for page in sorted(DOCS.rglob("*.md")):
        if "_tools" in page.parts:
            continue
        lines = page.read_text(encoding="utf-8").splitlines()
        for i, line in enumerate(lines):
            match = marker.search(line)
            if not match:
                continue
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j >= len(lines) or not lines[j].strip().startswith("```python"):
                raise SystemExit(f"{page}: marker {match.group(1)} is not followed by a python block")
            k = j + 1
            while k < len(lines) and not lines[k].strip().startswith("```"):
                k += 1
            found[match.group(1)] = (page.relative_to(ROOT).as_posix(), "\n".join(lines[j + 1:k]))
    return found


def check_examples(selected) -> list:
    results = []
    found = page_examples()
    for name in sorted(set(found) - set(EXAMPLES)):
        results.append((f"example {name}", "FAIL", f"marked in {found[name][0]} but not allowlisted"))
    for name, needs in EXAMPLES.items():
        if selected and name not in selected:
            continue
        if name not in found:
            results.append((f"example {name}", "FAIL", "allowlisted but no page marks it"))
            continue
        page, code = found[name]
        copy = make_copy(needs)
        try:
            before = snapshot(copy)
            start = time.perf_counter()
            done = run([PY, "-c", code], copy)
            spent = time.perf_counter() - start
            changed = changes(before, snapshot(copy))
            if done.returncode != 0:
                results.append((f"example {name}", "FAIL", f"{page}: exit {done.returncode}\n{done.stderr[-1500:]}"))
            elif changed:
                results.append((f"example {name}", "FAIL", f"{page}: changed {changed[:5]}"))
            else:
                results.append((f"example {name}", "pass", f"{page}, {spent:.1f} s, read-only"))
        finally:
            shutil.rmtree(copy, ignore_errors=True)
    return results


def check_commands(selected, keep: bool) -> list:
    results = []
    have_bash = shutil.which("bash") is not None and os.name == "posix"
    for command in COMMANDS:
        if selected and command.name not in selected:
            continue
        if command.shell and not have_bash:
            results.append((command.name, "skip", "needs a POSIX bash (on Windows, bash is WSL's and cannot open Windows paths)"))
            continue
        copy = make_copy(command.needs)
        try:
            for step in command.before:
                prepared = run(step, copy)
                if prepared.returncode != 0:
                    raise RuntimeError(f"setup step {step[1:]} exited {prepared.returncode}\n{prepared.stderr[-1500:]}")
            before = snapshot(copy)
            start = time.perf_counter()
            done = run(command.argv, copy)
            spent = time.perf_counter() - start
            after = snapshot(copy)
            changed = changes(before, after)
            stray = [p for p in changed if not p.startswith(command.writes)] if command.writes else changed
            staging = [p for p in after if "/.staging-" in f"/{p}"]
            problems = []
            if done.returncode != command.expect_exit:
                problems.append(f"exit {done.returncode}, expected {command.expect_exit}\n{(done.stderr or done.stdout)[-1500:]}")
            if command.expect_stdout and not re.search(command.expect_stdout, done.stdout, re.S):
                problems.append(f"output does not match {command.expect_stdout!r}")
            if stray:
                problems.append(f"wrote outside its documented outputs: {stray[:5]}")
            if staging:
                problems.append(f"left a staging directory behind: {staging[:2]}")
            kind = "writes only " + ", ".join(command.writes) if command.writes else "read-only"
            results.append((command.name, "FAIL" if problems else "pass",
                            "; ".join(problems) if problems else f"{spent:.1f} s, {kind}"))
        except (RuntimeError, subprocess.TimeoutExpired) as exc:
            results.append((command.name, "FAIL", str(exc)[:1500]))
        finally:
            if keep:
                print(f"    kept {copy}")
            else:
                shutil.rmtree(copy, ignore_errors=True)
    for name, partial, expect, pattern in CHECKPOINT_CASES:
        if selected and name not in selected:
            continue
        copy = make_copy([])
        try:
            case = _checkpoint_fixture(copy, partial)
            before = snapshot(copy)
            done = run([PY, "-m", "plumetools.foamio.checkpoints", "plan", str(case), "--np", "2"], copy)
            changed = changes(before, snapshot(copy))
            ok = done.returncode == expect and not changed and (
                pattern is None or re.search(pattern, done.stdout, re.S))
            detail = f"exit {done.returncode}" + ("" if ok else f"\n{done.stdout}{done.stderr}")
            results.append((name, "pass" if ok else "FAIL", detail + ", read-only" if ok else detail))
        finally:
            shutil.rmtree(copy, ignore_errors=True)
    return results


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--list", action="store_true", help="list the checks and stop")
    parser.add_argument("--only", nargs="+", default=None, metavar="NAME")
    parser.add_argument("--keep", action="store_true", help="keep command copies for inspection")
    args = parser.parse_args(argv)
    if args.list:
        for name in EXAMPLES:
            print(f"example   {name}")
        for command in COMMANDS:
            print(f"command   {command.name}")
        for name, *_ in CHECKPOINT_CASES:
            print(f"command   {name}")
        return 0
    try:
        import numpy, scipy, yaml  # noqa: F401,E401 -- plumetools' own dependencies
    except ImportError as exc:
        print(f"run_examples: {exc}; install the package first (pip install -e .)", file=sys.stderr)
        return 1
    if shutil.which("git") is None:
        print("run_examples: git is needed to copy the tracked files", file=sys.stderr)
        return 1
    selected = set(args.only) if args.only else None
    results = check_examples(selected) + check_commands(selected, args.keep)
    width = max(len(name) for name, *_ in results) if results else 0
    for name, status, detail in results:
        print(f"{status:<4}  {name:<{width}}  {detail.splitlines()[0] if status != 'FAIL' else detail}")
    failed = [r for r in results if r[1] == "FAIL"]
    passed = [r for r in results if r[1] == "pass"]
    skipped = [r for r in results if r[1] == "skip"]
    print(f"\nrun_examples: {len(passed)} passed, {len(skipped)} skipped, {len(failed)} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
