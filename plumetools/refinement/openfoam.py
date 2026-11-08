"""Portable v2512 hex refinement into fresh artifacts, never source cases."""
from __future__ import annotations

from contextlib import ExitStack
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
from typing import Optional

import numpy as np

from plumetools.foamio.patchfields import patch_field_type
from plumetools.foamio.volfield import read_internal
from plumetools.mesh.boundary import read_boundary
from plumetools.mesh.polymesh import read_faces_flat
from .assess import EDGE_COMPARISON_RTOL
from .geometry import read_geometry
from .model import RefinementPlan, RefinementResult
from .plan import _validate, json_safe, plan_identity, verify_sources, _publish_stage


def write_scalar_field(path: Path, values, dimensions, patches):
    """Stream an ASCII planning field with matching constraint patch types."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    values = np.asarray(values, dtype=float)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(f'FoamFile\n{{ version 2.0; format ascii; class volScalarField; location "{path.parent.name}"; object {path.name}; }}\n')
        stream.write(f'dimensions [{" ".join(map(str, dimensions))}];\ninternalField nonuniform List<scalar>\n{values.size}\n(\n')
        for start in range(0, values.size, 65536):
            stream.write("\n".join(format(float(v), ".17g") for v in values[start:start+65536]) + "\n")
        stream.write(");\nboundaryField\n{\n")
        for name, patch in patches.items():
            stream.write(f" {name} {{ type {patch_field_type(patch.type)}; }}\n")
        stream.write("}\n// ************************************************************************* //\n")


def write_cell_set(case: Path, name: str, cells):
    """Write a cellSet belonging only to the currently materialized mesh."""
    path = Path(case) / "constant/polyMesh/sets" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(f'FoamFile\n{{version 2.0; format ascii; class cellSet; location "constant/polyMesh/sets"; object {name};}}\n')
        stream.write(f"{len(cells)}\n(\n")
        for cell in cells:
            stream.write(f"{int(cell)}\n")
        stream.write(")\n// ************************************************************************* //\n")


def _run(command, cwd, log, *, env=None, timeout=300):
    """Execute argv, using a logical case path when OpenFOAM rejects whitespace."""
    log = Path(log)
    cwd = Path(cwd).resolve()
    with log.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("# " + shlex.join(map(str, command)) + "\n"); stream.flush()
        try:
            with ExitStack() as cleanup:
                logical = cwd
                if any(c.isspace() for c in str(cwd)):
                    temporary = Path(cleanup.enter_context(tempfile.TemporaryDirectory(prefix="foam-case-")))
                    logical = temporary / "case"
                    logical.symlink_to(cwd,target_is_directory=True)
                argv = [str(logical) if str(arg)==str(cwd) else str(arg) for arg in command]
                environment = dict(os.environ if env is None else env)
                environment["PWD"] = str(logical)
                result = subprocess.run(argv, cwd=logical, stdout=stream,
                                        stderr=subprocess.STDOUT, timeout=timeout, env=environment)
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise RuntimeError(f"{command[0]} failed: {exc}; log: {log}") from exc
    if result.returncode:
        tail = "\n".join(log.read_text(errors="replace").splitlines()[-24:])
        raise RuntimeError(f"{command[0]} exited {result.returncode}; log: {log}\n{tail}")


def _stage_system(case):
    system = case / "system"; system.mkdir()
    files = {
        "controlDict": ('application refineHexMesh;\nstartFrom startTime;\nstartTime 0;\nstopAt endTime;\n'
                        'endTime 1;\ndeltaT 1;\nwriteControl timeStep;\nwriteInterval 1;\n'
                        'writeFormat ascii;\nwritePrecision 17;\nwriteCompression off;\ntimeFormat general;\n'
                        'timePrecision 12;\nrunTimeModifiable false;\n'),
        "fvSchemes": ('ddtSchemes {default Euler;}\ngradSchemes {default Gauss linear;}\n'
                      'divSchemes {default none;}\nlaplacianSchemes {default Gauss linear corrected;}\n'
                      'interpolationSchemes {default linear;}\nsnGradSchemes {default corrected;}\n'),
        "fvSolution": 'solvers {}\n',
    }
    for name, body in files.items():
        (system / name).write_text(f"FoamFile\n{{version 2.0; format ascii; class dictionary; object {name};}}\n" + body, encoding="utf-8")


def _require_mesh_ok(log):
    """v2512 checkMesh returns zero even when it reports failed mesh checks."""
    text = Path(log).read_text(encoding="utf-8", errors="replace")
    if "Mesh OK." not in text or " mesh checks." in text:
        raise ValueError(f"mesh quality checks did not pass; see {log}")


def _check_topology(source):
    patches = read_boundary(source)
    if any(p.type in ("empty", "wedge", "processor", "processorCyclic") for p in patches.values()):
        raise ValueError("hex backend requires a reconstructed three-dimensional (3D) mesh topology; export targets for another backend")
    offsets, _ = read_faces_flat(source)
    if np.any(np.diff(offsets) != 4) and not (source / "constant/polyMesh/cellLevel").exists():
        raise ValueError("hex backend needs hex topology or compatible refinement history; use target export for another mesher")
    return patches


class HexRefiner:
    """Refine only measured failing cells plus v2512 conformity extensions.

    Work is done in a new staging case. The output carries inherited targets
    and density estimates, not a restarted cloud or new measured flow evidence.
    """
    def __init__(self, max_stages: int = 20, max_cells: Optional[int] = None):
        if isinstance(max_stages, bool) or int(max_stages) != max_stages or max_stages < 1:
            raise ValueError("max_stages must be a positive integer")
        if max_cells is not None and (int(max_cells) != max_cells or max_cells < 1):
            raise ValueError("max_cells must be a positive integer")
        self.max_stages, self.max_cells = int(max_stages), max_cells

    def refine(self, source: Path, plan: RefinementPlan, destination: Path) -> RefinementResult:
        """Publish a fresh mesh only once its inherited geometric bounds pass."""
        source, destination = Path(source).resolve(), Path(destination).resolve()
        _validate(plan)
        verify_sources(source, plan)
        if plan.action == "criteria_satisfied":
            return RefinementResult(None, "openfoam-hex", plan.mesh_sha256,
                                    len(plan.source_geometry.volumes_m3), "criteria_satisfied", False, [])
        if plan.action != "needs_mesh":
            raise ValueError(f"{plan.action}: a complete measured mesh plan is required")
        if destination == source or source in destination.parents:
            raise ValueError("output must be outside the source case")
        if destination.exists():
            raise ValueError(f"output already exists: {destination}")
        patches = _check_topology(source)
        initial_cells = len(plan.source_geometry.volumes_m3)
        if self.max_cells is not None and initial_cells + 7*len(plan.selected_cells) > self.max_cells:
            raise ValueError("requested cell refinement exceeds the cell budget; targets were not relaxed")
        if os.environ.get("WM_PROJECT_VERSION") not in ("v2512", "2512"):
            raise ValueError("source the OpenFOAM v2512 environment before mesh materialization")
        for command in ("refineHexMesh", "postProcess", "checkMesh"):
            if not shutil.which(command):
                raise ValueError(f"missing {command}; source the OpenFOAM v2512 environment")
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".hex-refinement-", dir=destination.parent) as temporary:
            stage = Path(temporary) / "case"; stage.mkdir()
            logs = stage / "logs"; logs.mkdir()
            history, log_names = [], []
            try:
                (stage / "constant").mkdir()
                shutil.copytree(source / "constant/polyMesh", stage / "constant/polyMesh")
                _stage_system(stage)
                write_scalar_field(stage / "0/refinementTargetEdge", plan.target_edge_m,
                                   [0, 1, 0, 0, 0, 0, 0], patches)
                write_scalar_field(stage / "0/refinementDensityEstimate", plan.density_estimate,
                                   [0, -3, 0, 0, 0, 0, 0], patches)
                current = plan.source_geometry
                target = plan.target_edge_m
                for step in range(1, self.max_stages+1):
                    selected = np.flatnonzero(current.max_edge_m > target * (1 + EDGE_COMPARISON_RTOL))
                    if not selected.size:
                        break
                    before = len(current.volumes_m3)
                    if self.max_cells is not None and before + 7*len(selected) > self.max_cells:
                        raise ValueError("cell budget exceeded; no incomplete mesh was published")
                    write_cell_set(stage, "refinementCells", selected)
                    name = f"{step:02d}-refineHexMesh.log"; log_names.append("logs/" + name)
                    _run(["refineHexMesh", "refinementCells", "-case", stage, "-overwrite"], stage, logs / name)
                    for function in ("writeCellCentres", "writeCellVolumes"):
                        name = f"{step:02d}-{function}.log"; log_names.append("logs/" + name)
                        _run(["postProcess", "-case", stage, "-func", function, "-time", "0"], stage, logs / name)
                    current = read_geometry(stage, time="0")
                    after = len(current.volumes_m3)
                    if self.max_cells is not None and after > self.max_cells:
                        raise ValueError("conformity refinement exceeded the cell budget; output was not published")
                    if after <= before:
                        raise ValueError("refinement made no mesh progress; topology/backend needs review")
                    target = read_internal(stage / "0/refinementTargetEdge").over(after, "cells")
                    if not np.all(np.isfinite(target) & (target > 0)):
                        raise ValueError("refinement lost valid inherited target fields")
                    history.append({"stage": step, "selected_cells": len(selected), "cells_before": before,
                                    "cells_after": after, "conformity_cells": (after-before)//7 - len(selected)})
                if np.any(current.max_edge_m > target * (1 + EDGE_COMPARISON_RTOL)):
                    raise ValueError("refinement stage limit reached before inherited targets were met")
                name = "checkMesh.log"; log_names.append("logs/" + name)
                _run(["checkMesh", "-case", stage, "-constant"], stage, logs / name)
                _require_mesh_ok(logs / name)
                document = json_safe({"kind": "plumetools.refinement.mesh", "schema_version": 1,
                                      "backend": "openfoam-hex", "plan_id": plan_identity(plan),
                                      "source_mesh_sha256": plan.mesh_sha256, "source_sha256": plan.source_sha256,
                                      "mesh_sha256": current.mesh_sha256, "cells": len(current.volumes_m3),
                                      "target_status": "inherited_targets_met", "requires_flow_rerun": True,
                                      "field_kind": "inherited_parent_estimates", "stages": len(history),
                                      "history": history, "logs": log_names, "provenance": plan.provenance})
                (stage / "refinement-result.json").write_text(json.dumps(document, indent=2, allow_nan=False) + "\n", encoding="utf-8")
                if destination.exists():
                    raise ValueError(f"output already exists: {destination}")
                _publish_stage(stage,destination)
            except Exception as exc:
                # Keep failure evidence outside the unpublished mesh. The
                # source and any existing requested output remain untouched.
                try:
                    failure = Path(tempfile.mkdtemp(prefix=destination.name + ".failure-", dir=destination.parent))
                    shutil.copytree(logs, failure / "logs")
                    (failure / "failure.txt").write_text(str(exc) + "\n", encoding="utf-8")
                except OSError:
                    pass
                raise
        return RefinementResult(destination, "openfoam-hex", current.mesh_sha256,
                                len(current.volumes_m3), "inherited_targets_met", True,
                                [str(destination / name) for name in log_names])
