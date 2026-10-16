"""The fresh-iteration helper: the previous run untouched, nothing to resume from.

The source here is a case as a finished run leaves it -- a mesh, a reconstructed
time with parcels and ``fieldAverage`` state in it, processor directories, logs
and results -- so that anything the helper carried over would show.
"""

from __future__ import annotations

import hashlib
import os

import pytest
import yaml

from _nozzleplume import R4D_CASE
from plumetools import iteration


def finished_run(case):
    """A case directory in the state ./Allrun and ./Allpost leave it."""
    case.mkdir()
    for name in ("case.yaml", "viz.yaml", "open.foam", "Allmesh", "Allrun",
                 "Allpost", "Allviz", "Allclean"):
        (case / name).write_bytes((R4D_CASE / name).read_bytes())
    (case / "Allrun").chmod(0o755)
    for path, text in {
            "constant/polyMesh/owner": "mesh", "constant/dsmcProperties": "dict",
            "system/controlDict": "startFrom latestTime;", "0/boundaryT": "field",
            "0.0517723/rhoNMean": "mean", "0.0517723/lagrangian/dsmc/positions": "p",
            "0.0517723/uniform/functionObjects/functionObjectProperties": "state",
            "processor0/0.0517723/rhoN": "field", "results/diagnostics_summary.yaml":
            "time: '0.0517723'\n", "log.dsmcFoam": "log", "case-summary.json": "{}"}.items():
        (case / path).parent.mkdir(parents=True, exist_ok=True)
        (case / path).write_text(text, encoding="utf-8")
    return case


def tree(root):
    """Every file under ``root`` and its content hash."""
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}


def test_a_fresh_iteration_copies_the_inputs_and_leaves_the_run_alone(tmp_path):
    source = finished_run(tmp_path / "baseCase")
    before = tree(source)
    made = iteration.prepare(source, tmp_path / "iter01", note="exit cone to level 5")

    assert tree(source) == before                         # read, never written
    destination = made.destination
    assert sorted(made.copied) == ["Allclean", "Allmesh", "Allpost", "Allrun",
                                   "Allviz", "case.yaml", "open.foam", "viz.yaml"]
    assert (destination / "case.yaml").read_bytes() == (source / "case.yaml").read_bytes()
    if os.name == "posix":
        assert os.access(destination / "Allrun", os.X_OK)

    # Nothing a solver, a mesher or post-processing could pick up.
    assert iteration.resumable_state(destination) == []
    assert iteration.resumable_state(source) != []
    assert iteration.times(destination) == []
    left = {p.name for p in destination.iterdir()}
    assert left == set(made.copied) | {"iteration.yaml", "provenance"}

    record = yaml.safe_load((destination / "iteration.yaml").read_text(encoding="utf-8"))
    assert record["parent"]["latest_time"] == "0.0517723"
    assert record["parent"]["case_yaml_sha256"] == before["case.yaml"]
    assert record["note"] == "exit cone to level 5"
    assert record["lineage"] == ["baseCase"]
    assert (destination / "provenance" / "parent-case.yaml").read_bytes() \
        == (source / "case.yaml").read_bytes()
    assert (destination / "provenance" / "parent-diagnostics_summary.yaml").is_file()

    # An iteration of an iteration extends the lineage, and does not inherit
    # its parent's record as its own.
    second = iteration.prepare(destination, tmp_path / "iter02")
    assert second.record["lineage"] == ["baseCase", "iter01"]
    assert "iteration.yaml" not in second.copied


@pytest.mark.parametrize("where, message", [
    ("existing", "not empty"), ("inside", "inside"), ("state", "solver state")])
def test_an_iteration_refuses_what_could_mix_runs(tmp_path, where, message):
    source = finished_run(tmp_path / "baseCase")
    before = tree(source)
    destination, include = tmp_path / "iter01", ()
    if where == "existing":
        (destination / "0.01").mkdir(parents=True)
    elif where == "inside":
        destination = source / "iter01"
    else:
        include = ("0.0517723/rhoNMean",)
    with pytest.raises(iteration.IterationError, match=message):
        iteration.prepare(source, destination, include=include)
    assert tree(source) == before
