"""The steps a plate-unigas case's scripts call: ``python -m plumetools.plateugf``.

::

    describe      what case.yaml implies: nozzleplume's report, plus memory
                  and disk estimates. Writes nothing.
    setup         describe, the checks, then every mesh and uniGasFoam dictionary.
    resume-plan   KEY=value lines for Allrun: where to resume, which leg.
    mark-restart  make uniGasFoam restart, not re-initialise, from a time.
    finished      exit 0 when the sampling leg has reached endTime.
    check-build   the built uniGasFoam is the commit case.yaml expects.
    smoke-check   the laptop smoke run's pass criteria (spec §8) on one write.

``resume-plan`` and ``finished`` read the run window from the generated
``system/controlDict.startup`` / ``.sampling`` rather than re-deriving the
case: nozzleplume's cell estimate takes ~30 s for these meshes.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import yaml

from plumetools.cai2012ugf import hpc as ugf_hpc
from plumetools.cai2012ugf.settings import load_settings
from plumetools.foamio.checkpoints import plan_restart, scan
from plumetools.foamio.errors import FoamDataError
from plumetools.foamio.volfield import read_internal_field, read_patch_field

#: Bytes per cell and per particle, measured on the laptop runs (memory note:
#: 5.6 GB at 1.31M cells + ~240 B/particle).
BYTES_PER_CELL = 4.3e3
BYTES_PER_PARTICLE = 240.0
#: Rough ASCII bytes per cell per written time (two uniGasVolFields sets).
DISK_BYTES_PER_CELL_WRITE = 600.0

_END = re.compile(r"^\s*endTime\s+([^;\s]+)\s*;", re.MULTILINE)
PLATE_FIELDS = ("p", "surfaceShearStress", "fD", "surfaceHeatTransfer")


# --------------------------------------------------------------------------- #
# the case, cheaply
# --------------------------------------------------------------------------- #

def _case_yaml(case_dir: Path) -> dict:
    return yaml.safe_load((Path(case_dir) / "case.yaml").read_text(encoding="utf-8")) or {}


def run_window(case_dir: Path) -> tuple[float, float]:
    """``(average_start_s, end_time_s)`` from the generated controlDicts."""
    out = []
    for leg in ("startup", "sampling"):
        path = Path(case_dir) / "system" / f"controlDict.{leg}"
        match = _END.search(path.read_text(encoding="utf-8")) if path.is_file() else None
        if not match:
            raise FoamDataError(f"{path}: no endTime -- run ./Allmesh (or --dict-only) first")
        out.append(float(match.group(1)))
    return out[0], out[1]


def resume_plan(case_dir: Path, n_subdomains: int) -> ugf_hpc.UgfRestartPlan:
    """cai2012ugf's restart plan, with this case's run window."""
    case_dir = Path(case_dir)
    base = plan_restart(case_dir, n_subdomains)
    average_start, end_time = run_window(case_dir)
    resume = float(base.resume) if base.resume else None
    leg = ugf_hpc.leg_for(resume, average_start, end_time,
                          ugf_hpc._slack(case_dir, base.resume))
    reasons = [base.refusal] if base.refusal else []
    if base.resume:
        decomposed = n_subdomains > 1 and base.proc_time == base.resume
        lost = ugf_hpc.lost_clouds(case_dir, base.resume, decomposed=decomposed)
        if lost:
            reasons.append(
                f"t = {base.resume}: {', '.join(lost)} held gas but "
                f"{'have' if len(lost) > 1 else 'has'} no lagrangian/{ugf_hpc.CLOUD} cloud "
                f"there; resuming would silently re-fill them. Move every later time aside "
                f"on every rank ('mv processorN/<t> processorN/<t>.lost') back to a time "
                f"whose clouds are complete, and resubmit. Nothing was changed.")
    return ugf_hpc.UgfRestartPlan(base=base, leg=leg, refusal="\n".join(reasons) or None)


def finished(case_dir: Path) -> bool:
    """The newest complete time has reached the end of the sampling leg."""
    latest = scan(case_dir).latest
    average_start, end_time = run_window(case_dir)
    return ugf_hpc.leg_for(latest[0] if latest else None, average_start, end_time,
                           ugf_hpc._slack(case_dir, latest[1] if latest else None)) == "done"


def check_build(case_dir: Path) -> str | None:
    """Why the built uniGasFoam is not the one this case expects, or ``None``."""
    expected = str((_case_yaml(case_dir).get("meta") or {}).get("unigas_commit") or "")
    if not expected:
        return f"{case_dir}/case.yaml has no meta.unigas_commit"
    built = ugf_hpc._built_commit()
    if built is None:
        return (f"no $FOAM_USER_APPBIN/{ugf_hpc.COMMIT_FILE}: build uniGasFoam with "
                f"cases/cai2012-unigas/runBuild.sh, which records the commit it built")
    if not (built.startswith(expected) or expected.startswith(built)):
        return f"uniGasFoam was built from {built}, but this case expects {expected}"
    return None


# --------------------------------------------------------------------------- #
# the smoke check (spec §8)
# --------------------------------------------------------------------------- #

def _latest_time(case_dir: Path) -> str | None:
    times = []
    for entry in Path(case_dir).iterdir():
        try:
            value = float(entry.name)
        except ValueError:
            continue
        if entry.is_dir() and value > 0.0:
            times.append((value, entry.name))
    return max(times)[1] if times else None


def _symmetry_patches(case_dir: Path) -> list[str]:
    boundary = Path(case_dir) / "constant" / "polyMesh" / "boundary"
    if not boundary.is_file():
        return []
    from plumetools.mesh.boundary import read_boundary
    return [name for name, info in read_boundary(case_dir).items()
            if info.type in ("symmetry", "symmetryPlane")]


def smoke_check(case_dir: Path, time: str | None = None) -> list[tuple[str, bool, str]]:
    """``(criterion, passed, detail)`` for one reconstructed write."""
    case_dir = Path(case_dir)
    time = time or _latest_time(case_dir)
    if time is None:
        return [("write_exists", False, "no reconstructed time directory")]
    t = case_dir / time
    cfg = _case_yaml(case_dir)
    names = (cfg.get("mesh") or {}).get("patch_names") or {}
    plate = names.get("plate", "plate")
    field = load_settings(case_dir).fields.mean_field
    results = []

    plate_values = {}
    for name in PLATE_FIELDS:
        path = t / f"{name}_{field}"
        try:
            plate_values[name] = np.asarray(read_patch_field(path, plate), dtype=float)
            ok = bool(np.all(np.isfinite(plate_values[name])))
            results.append((f"plate_{name}", ok,
                            f"{path.name} on '{plate}': {plate_values[name].shape[0]} faces"
                            + ("" if ok else ", NOT all finite")))
        except (OSError, FoamDataError) as exc:
            results.append((f"plate_{name}", False, f"{path.name} on '{plate}': {exc}"))

    p = plate_values.get("p")
    loaded = p is not None and bool(np.any(p > 0.0))
    results.append(("plate_loaded", loaded,
                    f"max p on the plate = {float(np.max(p)) if p is not None and p.size else 'n/a'} Pa"))

    rho_path = t / f"rhoN_{field}"
    try:
        rho = read_internal_field(rho_path)
        results.append(("volume_fields", bool(np.any(rho > 0.0)),
                        f"{rho_path.name}: {int(np.count_nonzero(rho > 0))} of {rho.size} cells hold gas"))
    except (OSError, FoamDataError) as exc:
        results.append(("volume_fields", False, f"{rho_path.name}: {exc}"))

    # A symmetry patch writes no value, so look at the cells that own its faces:
    # an absorbing plane would leave them empty next to the plume axis.
    symmetry = _symmetry_patches(case_dir)
    if symmetry:
        from plumetools.mesh.boundary import read_boundary
        from plumetools.mesh.polymesh import read_label_list
        try:
            owner = read_label_list(case_dir / "constant" / "polyMesh" / "owner")
            rho = read_internal_field(rho_path)
            patches = read_boundary(case_dir)
        except (OSError, FoamDataError) as exc:
            owner = None
            for patch in symmetry:
                results.append((f"gas_at_{patch}", False, str(exc)))
        if owner is not None:
            for patch in symmetry:
                info = patches[patch]
                cells = owner[info.start_face:info.start_face + info.n_faces]
                values = np.asarray(rho, dtype=float)[cells]
                results.append((f"gas_at_{patch}", bool(np.any(values > 0.0)),
                                f"{int(np.count_nonzero(values > 0))} of {values.size} cells "
                                f"beside '{patch}' hold gas"))
    return results


# --------------------------------------------------------------------------- #
# the symmetry planes, after snappyHexMesh
# --------------------------------------------------------------------------- #

_COUNT = re.compile(r"^\s*(\d+)\s*\n\s*\(", re.MULTILINE)


def flatten_symmetry(case_dir: Path, tolerance_m: float = 1.0e-3) -> int:
    """Put every point of the symmetry patch within ``tolerance_m`` of y = 0 or
    z = 0 exactly on that plane; return how many points moved.

    snappyHexMesh nudges boundary points where the plate's edge meets a
    symmetry plane (54 um on the z/D = 100 prototype), and nozzleplume's mesh
    verification holds the planes to 1e-3 of the smallest cell. The moves are
    tens of micrometres in cells of 10 mm or more.
    """
    from plumetools.mesh.boundary import read_boundary
    from plumetools.mesh.polymesh import read_faces_flat, read_points_array
    case_dir = Path(case_dir)
    path = case_dir / "constant" / "polyMesh" / "points"
    points = read_points_array(case_dir)
    offsets, labels = read_faces_flat(case_dir)
    on_plane = np.zeros(len(points), dtype=bool)
    for info in read_boundary(case_dir).values():
        if info.type not in ("symmetry", "symmetryPlane"):
            continue
        start, stop = offsets[info.start_face], offsets[info.start_face + info.n_faces]
        on_plane[np.asarray(labels[start:stop], dtype=int)] = True
    before = points.copy()
    for axis in (1, 2):
        near = on_plane & (np.abs(points[:, axis]) < tolerance_m)
        points[near, axis] = 0.0
    moved = int(np.count_nonzero(np.any(points != before, axis=1)))
    if moved:
        text = path.read_text(encoding="utf-8")
        match = _COUNT.search(text, text.index("}") + 1)
        if match is None or int(match.group(1)) != len(points):
            raise FoamDataError(f"{path}: not an ASCII point list this can rewrite")
        body = "\n".join(f"({x:.17g} {y:.17g} {z:.17g})" for x, y, z in points)
        path.write_text(text[:match.start()] + f"{len(points)}\n(\n{body}\n)\n\n"
                        "// ************************************************************************* //\n",
                        encoding="utf-8", newline="\n")
    return moved


# --------------------------------------------------------------------------- #
# describe / setup
# --------------------------------------------------------------------------- #

def estimates(case) -> list[str]:
    """Memory and disk estimates for ``case`` (a nozzleplume Case)."""
    cells = float(sum(case.plan.cells_by_level.values()))
    parcels = float(case.run.estimated_parcels)
    memory = cells * BYTES_PER_CELL + parcels * BYTES_PER_PARTICLE
    startup_writes = max(1, round(case.run.average_start_s / case.run.transit_s))
    sampling_writes = 4 * max(1, round((case.run.end_time_s - case.run.average_start_s)
                                       / case.run.transit_s))
    disk = (startup_writes + sampling_writes) * cells * DISK_BYTES_PER_CELL_WRITE * 2
    return [
        "Resources (uniGasFoam, rough)",
        f"  cells (estimate)         {cells:,.0f}",
        f"  parcels (estimate)       {parcels:.3g}",
        f"  memory                   {memory / 1e9:.1f} GB  "
        f"({BYTES_PER_CELL / 1e3:.1f} kB/cell + {BYTES_PER_PARTICLE:.0f} B/parcel, laptop-measured)",
        f"  writes                   {startup_writes} startup + {sampling_writes} sampling",
        f"  disk                     ~{disk / 1e9:.0f} GB (per-core + reconstructed copies)",
    ]


def describe(case_dir: Path) -> int:
    from plumetools.nozzleplume import workflow as npw
    case = npw.load(case_dir)
    report = case.checks()
    print("\n".join(npw.report(case, report) + estimates(case)))
    return 1 if report.failures else 0


def setup(case_dir: Path) -> int:
    from plumetools.nozzleplume import checks
    from plumetools.nozzleplume import workflow as npw
    from plumetools.plateugf import dictionaries
    case = npw.load(case_dir)
    report = case.checks()
    print("\n".join(npw.report(case, report) + estimates(case)))
    try:
        report.raise_if_failed()
    except checks.DsmcCheckError as exc:
        print(str(exc), file=sys.stderr)
        print("\nNo dictionaries were written.", file=sys.stderr)
        return 1
    for path in dictionaries.write_case_dictionaries(case, load_settings(case_dir)):
        print(f"    wrote {path}")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m plumetools.plateugf")
    steps = parser.add_subparsers(dest="step", required=True)
    for name in ("describe", "setup", "finished", "check-build", "flatten-symmetry"):
        steps.add_parser(name).add_argument("case_dir", type=Path)
    plan = steps.add_parser("resume-plan")
    plan.add_argument("case_dir", type=Path)
    plan.add_argument("--np", type=int, required=True)
    mark = steps.add_parser("mark-restart")
    mark.add_argument("case_dir", type=Path)
    mark.add_argument("--time", required=True)
    mark.add_argument("--np", type=int, required=True)
    smoke = steps.add_parser("smoke-check")
    smoke.add_argument("case_dir", type=Path)
    smoke.add_argument("--time")
    post_step = steps.add_parser("post")
    post_step.add_argument("case_dir", type=Path)
    post_step.add_argument("--time")
    table_step = steps.add_parser("family-table")
    table_step.add_argument("root", type=Path)
    window = steps.add_parser("stop-at")
    window.add_argument("case_dir", type=Path)
    window.add_argument("transits", type=float)
    args = parser.parse_args(argv)

    if args.step == "describe":
        return describe(args.case_dir)
    if args.step == "setup":
        return setup(args.case_dir)
    if args.step == "post":
        from plumetools.plateugf import post as post_module
        for path in post_module.run_post(args.case_dir, args.time):
            print(f"    wrote {path}")
        return 0
    if args.step == "family-table":
        from plumetools.plateugf import post as post_module
        print("
".join(post_module.family_table(args.root)))
        return 0
    if args.step == "flatten-symmetry":
        print(f"    flatten-symmetry: {flatten_symmetry(args.case_dir)} point(s) put on y = 0 / z = 0")
        return 0
    if args.step == "finished":
        return 0 if finished(args.case_dir) else 1
    if args.step == "check-build":
        problem = check_build(args.case_dir)
        if problem:
            print(f"check-build: {problem}", file=sys.stderr)
            return 1
        return 0
    if args.step == "resume-plan":
        result = resume_plan(args.case_dir, args.np)
        sys.stdout.write(result.shell())
        if result.refusal:
            print(result.refusal, file=sys.stderr)
            return 3
        return 0
    if args.step == "mark-restart":
        made = ugf_hpc.mark_restart(args.case_dir, args.time, decomposed=args.np > 1)
        if made is not None:
            print(f"    created {made} so uniGasFoam restarts on every rank")
        return 0
    if args.step == "stop-at":
        average_start, end_time = run_window(args.case_dir)
        transits = (_case_yaml(args.case_dir).get("dsmc") or {}).get("transient_domain_transits", 6.0)
        print("%.12g" % (average_start / float(transits) * args.transits))
        return 0
    results = smoke_check(args.case_dir, args.time)
    for name, ok, detail in results:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}: {detail}")
    return 0 if all(ok for _, ok, _ in results) else 1
