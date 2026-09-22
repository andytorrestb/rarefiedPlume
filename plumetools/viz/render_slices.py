#!/usr/bin/env python3
r"""Export slice imagery from a dsmcFoam case.

    vifpara plumetools/viz/render_slices.py <case-dir> [options]

NOT ``python``. ``paraview.simple`` exists only inside ParaView's own
interpreter; the ``vifpara`` launcher re-runs this file under ``pvpython`` with
this virtual environment injected. Running it with ``python`` fails at the
import, and the error says so.

Examples::

    # everything in the default spec, latest time
    vifpara plumetools/viz/render_slices.py cases/cai2012/Cases/Kn100

    # a whole study, in ONE ParaView process -- what AllpostCases does, with
    # that study's own sampling settings
    vifpara plumetools/viz/render_slices.py cases/cai2012/Cases/* \
        --spec cases/cai2012/viz.yaml

    # one field, one plane
    vifpara plumetools/viz/render_slices.py cases/cai2012/Cases/Kn100 \
        --field rhoN --plane xz

    # a spec of your own, and every written time
    vifpara plumetools/viz/render_slices.py cases/cai2012/Cases/Kn100 \
        --spec my-slices.yaml --time all

    # what does this case actually have?
    vifpara plumetools/viz/render_slices.py cases/cai2012/Cases/Kn100 --list

Options that duplicate a spec key (``--time``, ``--output``, ``--preset``,
``--field-type``) override it for this run only. Everything else lives in the
YAML -- see ``plumetools/viz/slices.yaml`` and ``docs/viz-slices.md``.
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import replace
from pathlib import Path

try:
    from plumetools.viz.spec import (FIELD_TYPES, TIME_SELECTORS, VizSpecError,
                                     load_spec)
except ImportError:                                        # pragma: no cover
    sys.exit('plumetools is not importable. From the repository root:\n'
             '    pip install -e ".[test]"')


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="render_slices.py",
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Run through the 'vifpara' launcher, not through 'python'.")
    parser.add_argument("case", nargs="+",
                        help="one or more OpenFOAM case directories. Several "
                             "are rendered in a single ParaView process, which "
                             "is what a study loop should do -- starting one "
                             "per case pays ParaView's startup every time")
    parser.add_argument("--spec", default=None, metavar="FILE",
                        help="sampling configuration "
                             "(default: plumetools/viz/slices.yaml)")
    parser.add_argument("--field", action="append", default=None, metavar="NAME",
                        help="draw only this field; repeatable")
    parser.add_argument("--plane", action="append", default=None, metavar="NAME",
                        help="draw only on this plane; repeatable")
    parser.add_argument("--time", default=None, metavar="T",
                        help=f"override sampling.time: "
                             f"{', '.join(TIME_SELECTORS)}, or a value in seconds")
    parser.add_argument("--output", default=None, metavar="DIR",
                        help="override output.directory")
    parser.add_argument("--preset", default=None, metavar="NAME",
                        help="override image.preset, e.g. 'Turbo'")
    parser.add_argument("--field-type", default=None, choices=FIELD_TYPES,
                        help="override sampling.field_type")
    parser.add_argument("--instantaneous", action="store_true",
                        help="turn OFF prefer_mean and draw the single-timestep "
                             "fields. This is a picture of shot noise; it is "
                             "useful for judging particle statistics and for "
                             "nothing else")
    parser.add_argument("--stop-on-error", action="store_true",
                        help="give up at the first case that fails, instead of "
                             "rendering the rest and reporting at the end")
    parser.add_argument("--allow-constant", action="store_true",
                        help="draw fields that are constant over the geometry "
                             "instead of skipping them")
    parser.add_argument("--list", action="store_true",
                        help="list the arrays the case actually has at the "
                             "resolved time, then exit without drawing")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the images that would be drawn, then exit")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="list the available arrays at each time")
    return parser


def main(argv: list[str]) -> int:
    args = build_parser().parse_args(argv[1:])

    try:
        spec = load_spec(args.spec)
    except VizSpecError as exc:
        print(f"render_slices: {exc}", file=sys.stderr)
        return 2

    # -- command-line overrides ------------------------------------------- #
    sampling = spec.sampling
    if args.time is not None:
        time = args.time if args.time in TIME_SELECTORS else float(args.time)
        sampling = replace(sampling, time=time)
    if args.field_type is not None:
        sampling = replace(sampling, field_type=args.field_type)
    if args.instantaneous:
        sampling = replace(sampling, prefer_mean=False)

    image = spec.image
    if args.preset is not None:
        image = replace(image, preset=args.preset)

    output = spec.output
    if args.output is not None:
        output = replace(output, directory=args.output)

    spec = replace(spec, sampling=sampling, image=image, output=output)

    # -- dry run needs no ParaView ----------------------------------------- #
    if args.dry_run:
        try:
            tasks = spec.tasks(only_fields=args.field, only_planes=args.plane)
        except VizSpecError as exc:
            print(f"render_slices: {exc}", file=sys.stderr)
            return 2
        derived = spec.derived_by_name
        print(f"### {len(tasks)} image(s) per case from {spec.source}")
        for task in tasks:
            candidates = task.field.candidate_names(
                spec.prefer_mean_for(task.field), derived)
            print(f"    {task.field.name:<20} {task.plane_name:<10} "
                  f"{task.field.resolve_kind(derived).value:<8} "
                  f"<- {' or '.join(candidates)}")
        print(f"    x {len(args.case)} case(s) = "
              f"{len(tasks) * len(args.case)} image(s)")
        return 0

    # -- everything past here needs ParaView -------------------------------- #
    try:
        from plumetools.viz.render import (NothingToRender, RenderError,
                                           Renderer, discover_case)
    except ImportError as exc:
        print(f"render_slices: {exc}\n"
              f"\n"
              f"       This script must be run through the VifPara launcher, "
              f"which puts\n"
              f"       this virtual environment inside pvpython:\n"
              f"\n"
              f"           vifpara {Path(__file__).name} <case-dir>\n"
              f"\n"
              f"       'python {Path(__file__).name}' cannot work: "
              f"paraview.simple only\n"
              f"       exists inside ParaView's own interpreter. See "
              f"docs/environment.md.",
              file=sys.stderr)
        return 2

    written = 0
    skipped = 0
    not_run: list[str] = []
    failed: list[str] = []

    for index, case_dir in enumerate(args.case):
        if index:
            print()
            print("=" * 62)
        renderer = None
        try:
            facts = discover_case(case_dir)
            case_spec = spec
            if args.output is not None and len(args.case) > 1:
                # The default output directory is relative and so already lands
                # inside each case. One explicit --output shared by several
                # cases would have each of them overwrite the last, silently,
                # leaving a directory of images all labelled with the same
                # field names but only the final case's data.
                case_spec = replace(spec, output=replace(
                    spec.output,
                    directory=str(Path(args.output) / facts.name)))
            renderer = Renderer(case_spec, facts,
                                allow_constant=args.allow_constant,
                                verbose=args.verbose)

            if args.list:
                _list_arrays(renderer, case_spec)
                continue

            renderer.run(only_fields=args.field, only_planes=args.plane)
            written += len(renderer.written)
            skipped += len(renderer.skipped)
            print()
            print(f"### {len(renderer.written)} image(s) in "
                  f"{renderer.output_dir}")
            if renderer.skipped:
                print(f"    {len(renderer.skipped)} skipped; "
                      f"see the 'skip' lines above and manifest.yaml")
        except NothingToRender as exc:
            # Not a failure. A study is routinely part-way through, and the
            # per-case Allpost this step is bolted onto reports what it can and
            # carries on; exiting non-zero over a picture of a case nobody has
            # run yet would be stricter than the validation it accompanies.
            print(f"    nothing to draw: {exc}")
            not_run.append(case_dir)
        except (RenderError, VizSpecError) as exc:
            print(f"\nrender_slices: {case_dir}: {exc}", file=sys.stderr)
            failed.append(case_dir)
            if args.stop_on_error:
                return 1
        finally:
            # Between cases, not after the last one only: readers and views are
            # global session state and the next case would inherit them.
            if renderer is not None:
                renderer.close()

    if len(args.case) > 1 and not args.list:
        print()
        print("=" * 62)
        print(f"### {written} image(s) over "
              f"{len(args.case) - len(not_run)} of {len(args.case)} case(s), "
              f"{skipped} field(s) skipped")
        if not_run:
            print(f"    {len(not_run)} case(s) have not been run and had "
                  f"nothing to draw:")
            for case_dir in not_run:
                print(f"      {case_dir}")

    if failed:
        print(f"\nrender_slices: {len(failed)} case(s) failed:", file=sys.stderr)
        for case_dir in failed:
            print(f"    {case_dir}", file=sys.stderr)
        return 1
    return 0


def _list_arrays(renderer, spec) -> int:
    """Print what the case has at the resolved time, and what the spec wants."""
    from plumetools.viz.catalog import catalog_entry
    from plumetools.viz.render import array_ranges, component_range, resolve_times

    times = resolve_times(renderer.reader.TimestepValues, spec.sampling.time,
                          time_min=spec.sampling.time_min,
                          time_max=spec.sampling.time_max)
    time = times[-1]
    available = array_ranges(renderer.reader, time,
                             spec.sampling.field_type)

    print(f"### {renderer.facts.name} at t = {format(time, 'g')} "
          f"({spec.sampling.field_type.lower()} data)")
    print()
    print(f"    {'array':<26} {'min':>12} {'max':>12}  description")
    for name in sorted(available):
        low, high = component_range(available[name], "")
        entry = catalog_entry(name)
        note = entry.description
        if entry.kind.value == "surface":
            note = f"[surface field] {note}"
        print(f"    {name:<26} {low:>12.4g} {high:>12.4g}  {note}")

    print()
    print("    The spec asks for:")
    derived = spec.derived_by_name
    for field in spec.fields:
        if not field.enabled:
            continue
        candidates = field.candidate_names(spec.prefer_mean_for(field), derived)
        found = next((name for name in candidates if name in available), None)
        if field.name in derived:
            status = "derived at render time"
        elif found:
            status = f"-> {found}"
        else:
            status = f"MISSING (looked for {', '.join(candidates)})"
        print(f"      {field.name:<20} {status}")
    return 0


def _exit(status: int) -> None:
    """Disconnect from ParaView's server before the interpreter tears it down.

    On a headless/WSLg display, ``pvpython --force-offscreen-rendering`` dies at
    exit with::

        X Error of failed request:  GLXBadContext
        Minor opcode of failed request:  5 (X_GLXMakeCurrent)

    releasing its GL contexts during interpreter finalisation, and the process
    ends up with status 1. Every image and the manifest are already written and
    closed by then, so a successful render reports failure -- exactly backwards
    for anything that loops over cases and checks ``$?``, which is what
    ``AllpostCases`` does.

    Disconnecting first releases the contexts in a defined order and the crash
    does not happen: measured here, ``pv.Disconnect()`` gives exit 0 with the
    full log intact.

    (The obvious shortcut, ``os._exit``, is wrong. VifPara imports colorama,
    which replaces ``sys.stdout`` with a ``StreamWrapper``; flushing that does
    not push the underlying buffer out, so ``os._exit`` silently discards every
    line the run printed. That was measured too.)

    ``PLUMETOOLS_VIZ_RAW_EXIT=1`` skips the disconnect and lets the crash
    happen, which is what you want when debugging ParaView rather than using it.
    """
    if os.environ.get("PLUMETOOLS_VIZ_RAW_EXIT") != "1":
        paraview = sys.modules.get("paraview.simple")
        if paraview is not None:
            try:
                paraview.Disconnect()
            except Exception:
                # Teardown is a courtesy. Whatever happened here, the images are
                # already on disk and the status below is the real answer.
                pass
    sys.stdout.flush()
    sys.stderr.flush()
    raise SystemExit(status)


if __name__ == "__main__":
    _exit(main(sys.argv))
