"""The steps a case's ``Allmesh``, ``Allrun`` and ``Allpost`` call.

A case directory holds a ``case.yaml`` and four thin shell scripts. Everything
they do that is not an OpenFOAM utility is one of these::

    describe   what case.yaml implies -- the exit state, the mesh, the cost.
               Writes nothing, needs no OpenFOAM. The command to run after
               editing dsmc.n_equivalent_particles.
    setup      describe, then write system/ and constant/ dictionaries
    verify     measure the built mesh against case.yaml
    fields     write 0/ for the mesh that was built
    post       plate loads and centreline into results/
    diagnostics
               mean free path, Knudsen numbers and mesh resolution: volume
               fields into the time directory, a summary into results/

Kept here rather than as Python embedded in each case's scripts, so a second
case of this family is a ``case.yaml`` and a copy of four files with no logic in
them.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

from plumetools.mesh.boundary import read_boundary
from plumetools.nozzleplume import (
    checks,
    diagnostics as diagnostics_module,
    dictionaries,
    exitstate,
    foamfields,
    gas,
    geometry,
    mesh,
    post as post_module,
    run as run_module,
    verify as verify_module,
)
from plumetools.nozzleplume.config import CaseConfig, ConfigError, load_case_config


@dataclass(frozen=True)
class Case:
    """Everything derived from one ``case.yaml``, before any mesh exists.

    Attributes:
        case_dir: the case directory.
        cfg: the validated ``case.yaml``.
        geom: domain, nozzle and plate geometry [m], symmetry planes.
        exit_state: the uniform exit state and the gas.
        plan: the octree mesh, planned, not written.
        nozzle: the nozzle patch ``topoSet`` + ``createPatch`` will carve.
        run: time step, averaging window, particle weight.
    """

    case_dir: Path
    cfg: CaseConfig
    geom: geometry.PlumeGeometry
    exit_state: exitstate.ExitState
    plan: mesh.MeshPlan
    nozzle: mesh.NozzlePrediction
    run: run_module.RunSettings

    @property
    def label(self) -> str:
        """The case's name for reports and figure titles."""
        return str(self.cfg.meta.get("case_name") or self.case_dir.resolve().name)

    def checks(self, patches=None):
        """The quality checks, with the mesh's patches once it has some."""
        return checks.run_checks(
            self.cfg, self.geom, self.exit_state, self.plan, self.run,
            nozzle=self.nozzle, patches=patches)


def load(case_dir: Path) -> Case:
    """Load ``case.yaml`` and derive everything that follows from it.

    Raises:
        ConfigError, GeometryError, ValueError: each naming the ``case.yaml``
            key responsible.
    """
    case_dir = Path(case_dir)
    cfg = load_case_config(case_dir)
    geom = geometry.from_config(cfg)
    exit_state = exitstate.from_config(cfg)
    plan = mesh.plan(cfg, geom, exit_state)
    nozzle = mesh.predicted_nozzle_faces(plan, geom)
    run = run_module.derive(cfg, exit_state, geom, plan, nozzle)
    return Case(case_dir=case_dir, cfg=cfg, geom=geom, exit_state=exit_state,
                plan=plan, nozzle=nozzle, run=run)


def report(case: Case, check_report) -> list[str]:
    """The whole pre-run report: geometry, gas, exit, mesh, run, checks."""
    blocks = [
        [f"### {case.label}"],
        case.geom.describe(),
        gas.describe(case.exit_state.species),
        case.exit_state.describe(),
        case.plan.describe(case.geom),
        case.nozzle.describe(),
        case.run.describe(),
        check_report.report(),
    ]
    lines: list[str] = []
    for block in blocks:
        lines += block + [""]
    return lines


def describe(case_dir: Path) -> int:
    """Print what ``case.yaml`` implies. Writes nothing."""
    case = load(case_dir)
    check_report = case.checks()
    print("\n".join(report(case, check_report)))
    if check_report.failures:
        print(f"{len(check_report.failures)} check(s) FAILED; ./Allmesh would "
              f"refuse to write this case.", file=sys.stderr)
        return 1
    return 0


def setup(case_dir: Path) -> int:
    """Print the report, then write every mesh and case dictionary."""
    case = load(case_dir)
    check_report = case.checks()
    print("\n".join(report(case, check_report)))
    try:
        check_report.raise_if_failed()
    except checks.DsmcCheckError as exc:
        print(str(exc), file=sys.stderr)
        print("\nNo dictionaries were written. Fix case.yaml and re-run.",
              file=sys.stderr)
        return 1

    written = mesh.write_mesh_setup(case.case_dir, case.cfg, case.geom, case.plan)
    written += dictionaries.write_case_dictionaries(
        case.case_dir, case.cfg, case.geom, case.exit_state, case.run)
    for path in written:
        print(f"    wrote {path}")
    return 0


def _patches(case: Case):
    """The built mesh's patches, or ``None`` with a message if there is no mesh."""
    if not (case.case_dir / "constant" / "polyMesh" / "boundary").is_file():
        print(f"no constant/polyMesh in {case.case_dir}. Generate it first:\n"
              f"    ./Allmesh", file=sys.stderr)
        return None
    return read_boundary(case.case_dir)


def verify(case_dir: Path) -> int:
    """Measure the built mesh against ``case.yaml``."""
    case = load(case_dir)
    if _patches(case) is None:
        return 1
    return verify_module.main(case.case_dir, case.cfg, case.geom, case.plan,
                              case.nozzle)


def fields(case_dir: Path, run_checks: bool = True) -> int:
    """Write ``0/`` for the mesh that was built, and ``case-summary.json``."""
    case = load(case_dir)
    patches = _patches(case)
    if patches is None:
        return 1

    name = case.cfg.mesh.patch_names["nozzle"]
    meshed = patches[name].n_faces if name in patches else 0
    print("\n".join(case.exit_state.describe()))
    print()
    print(f"  nozzle patch             {name}: {meshed} faces "
          f"(predicted {case.nozzle.n_faces}), {case.nozzle.area_m2:.6e} m^2")
    print(f"  injection                "
          f"{case.run.injection_rate_per_s:.6e} molecules/s = "
          f"{case.run.parcels_per_step:.4g} parcels/step at "
          f"nEquivalentParticles {case.run.n_equivalent_particles:.4e}")
    print()

    check_report = None
    if run_checks:
        check_report = case.checks(patches)
        print("\n".join(check_report.report()))
        print()
        try:
            check_report.raise_if_failed()
        except checks.DsmcCheckError as exc:
            print(str(exc), file=sys.stderr)
            print("\nNo 0/ fields were written. Fix case.yaml and re-run "
                  "./Allmesh.", file=sys.stderr)
            return 1

    for path in foamfields.write_all(case.case_dir, case.cfg, case.geom,
                                     case.exit_state, patches):
        print(f"    wrote {path}")

    checks.write_summary(case.case_dir, {
        "case": dict(case.cfg.meta),
        "exit": case.exit_state.as_dict(),
        "mesh": case.plan.as_dict(),
        "nozzle": {**case.nozzle.as_dict(), "patch": name, "meshed_faces": meshed},
        "run": case.run.as_dict(),
        "patches": {n: {"type": p.type, "n_faces": p.n_faces}
                    for n, p in patches.items()},
        "checks": check_report.as_dict() if check_report is not None else None,
    })
    print("    wrote case-summary.json")
    return 0


def post(case_dir: Path, time: str | None = None,
         results: str | None = None) -> int:
    """Extract the plate loads and the centreline into ``results/``."""
    case = load(case_dir)
    out = Path(results) if results else case.case_dir / "results"
    print(f"### {case.label}")

    loads = None
    if case.geom.plate is not None:
        try:
            loads = post_module.plate_loads(case.case_dir, case.cfg, case.geom,
                                            time=time)
        except (post_module.PostError, KeyError) as exc:
            print(f"post: {exc}", file=sys.stderr)
            return 1
        print("\n".join(loads.report(case.exit_state)))
        print()
        for path in post_module.write_plate_results(out, loads, case.exit_state):
            print(f"    wrote {path}")

    # The centreline needs cell centres, which only ./Allpost writes. Without
    # them the plate loads above are still a complete result, so this is
    # reported and skipped rather than made a failure.
    profile = None
    try:
        profile = post_module.centreline(case.case_dir, case.cfg, case.geom,
                                         case.exit_state, time=time)
    except post_module.PostError as exc:
        print(f"    no centreline: {exc}")
    else:
        print(f"    centreline           {profile['x_m'].size} of "
              f"{case.cfg.post.centreline_points} bins hold cells")
        print(f"    wrote {post_module.write_centreline_results(out, profile)}")

    figures = post_module.plot_results(out, loads, profile, case.label)
    for figure in figures:
        print(f"    wrote {figure}")
    if not figures and (loads is not None or profile is not None):
        print('    no plots: matplotlib is not installed '
              '(pip install -e ".[plots]")')
    return 0


def diagnostics(case_dir: Path, time: str | None = None,
                results: str | None = None) -> int:
    """Rarefaction and resolution diagnostics of one sampled time.

    Needs the averaged fields and the cell centres and volumes ``./Allpost``
    writes. Changes no mesh and no run: it is what to read before editing the
    refinement in ``case.yaml``. See :mod:`plumetools.nozzleplume.diagnostics`.
    """
    case = load(case_dir)
    print(f"### {case.label}")
    if not case.cfg.diagnostics.enabled:
        print("    diagnostics.enabled is false in case.yaml; nothing to do")
        return 0
    try:
        fields, written = diagnostics_module.run(
            case, time=time, results=Path(results) if results else None)
    except post_module.PostError as exc:
        print(f"diagnostics: {exc}", file=sys.stderr)
        return 1
    if fields:
        print(f"    wrote {len(fields)} fields into {fields[0].parent}: "
              + " ".join(path.name for path in fields))
    for path in written:
        print(f"    wrote {path}")
    return 0


def main(argv=None) -> int:
    """``python -m plumetools.nozzleplume <command> [case]``."""
    parser = argparse.ArgumentParser(
        prog="python -m plumetools.nozzleplume",
        description="Nozzle-exit plume cases: everything OpenFOAM reads is "
                    "generated from case.yaml.")
    commands = parser.add_subparsers(dest="command", required=True)

    def add(name: str, summary: str):
        sub = commands.add_parser(name, help=summary, description=summary)
        sub.add_argument("case", nargs="?", default=".",
                         help="case directory (default: the current one)")
        return sub

    add("describe", "print what case.yaml implies; write nothing")
    add("setup", "write the system/ and constant/ dictionaries")
    add("verify", "measure the built mesh against case.yaml")
    sub = add("fields", "write 0/ for the mesh that was built")
    sub.add_argument("--no-checks", action="store_true",
                     help="write 0/ without running the DSMC quality checks")
    sub = add("post", "plate loads and centreline into results/")
    sub.add_argument("--time", default=None,
                     help="time directory to read (default: the latest)")
    sub.add_argument("--results", default=None, metavar="DIR",
                     help="where to write (default: <case>/results)")
    sub = add("diagnostics", "mean free path, Knudsen numbers and mesh "
                             "resolution of a sampled time")
    sub.add_argument("--time", default=None,
                     help="time directory to read (default: the latest)")
    sub.add_argument("--results", default=None, metavar="DIR",
                     help="where the summary goes (default: <case>/results); "
                          "the fields always go into the time directory")

    args = parser.parse_args(argv)
    case_dir = Path(args.case)
    try:
        if args.command == "describe":
            return describe(case_dir)
        if args.command == "setup":
            return setup(case_dir)
        if args.command == "verify":
            return verify(case_dir)
        if args.command == "fields":
            return fields(case_dir, run_checks=not args.no_checks)
        if args.command == "diagnostics":
            return diagnostics(case_dir, time=args.time, results=args.results)
        return post(case_dir, time=args.time, results=args.results)
    except (ConfigError, geometry.GeometryError) as exc:
        print(f"{args.command}: {exc}", file=sys.stderr)
        return 1
    except ValueError as exc:
        # Raised by the mesh plan and the run settings, each naming the
        # case.yaml key to change. A traceback would bury that.
        print(f"{args.command}: {case_dir / 'case.yaml'}: {exc}", file=sys.stderr)
        return 1
