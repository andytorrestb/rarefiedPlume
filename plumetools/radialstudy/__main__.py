"""Study entry points; case wrappers delegate to the same decisions."""
import argparse
import json
from pathlib import Path
import sys


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command',required=True)
    g = sub.add_parser('generate')
    g.add_argument('--study',type=Path,default=Path('study.yaml'))
    g.add_argument('--output',type=Path,default=Path('Cases'))
    g.add_argument('--kn',type=float,nargs='+')
    g.add_argument('--dry-run',action='store_true')
    r = sub.add_parser('refine')
    r.add_argument('case',type=Path)
    r.add_argument('--plan',type=Path)
    r.add_argument('--correction',choices=['needed','mesh','parcels','timestep'],default='needed')
    r.add_argument('--delta-t',type=float)
    r.add_argument('--dry-run',action='store_true')
    batch = sub.add_parser('run-all',help='run generated cases sequentially; stop on the first failure')
    batch.add_argument('--cases',type=Path,default=Path('Cases'),help='generated case tree')
    batch.add_argument('--resume',action='store_true',help='skip finished cases, resume checkpoints and start untouched cases')
    batch.add_argument('--dry-run',action='store_true',help='report the queue without writing or launching commands')
    batch.add_argument('--allow-budget-exceedance',action='store_true')
    for name in ('mesh','run','post','analyse'):
        a = sub.add_parser(name)
        a.add_argument('case',type=Path)
        if name == 'run':
            a.add_argument('--resume',action='store_true')
            a.add_argument('--allow-budget-exceedance',action='store_true')
    a = p.parse_args(argv)
    try:
        if a.command == 'generate':
            from .config import load_study
            from .generation import generate_initial
            result = generate_initial(load_study(a.study),a.output,kn_values=a.kn,dry_run=a.dry_run)
        elif a.command == 'refine':
            from .generation import revise_case
            if a.plan:
                from plumetools.refinement.plan import read_plan
                plan = read_plan(a.plan)
                dt = a.delta_t
            else:
                from .post import case_plan
                plan, recommended_dt = case_plan(a.case)
                dt = a.delta_t or recommended_dt
            result = revise_case(a.case,plan,correction=a.correction,delta_t_s=dt,dry_run=a.dry_run)
        elif a.command == 'run-all':
            from .run import run_study
            result = run_study(a.cases,resume=a.resume,dry_run=a.dry_run,
                               allow_budget_exceedance=a.allow_budget_exceedance,
                               progress=lambda message:print(message,file=sys.stderr,flush=True))
        elif a.command in ('mesh','run'):
            from . import run
            from .case import load_case
            if a.command=='run' and load_case(a.case)[0].config.controller.enabled:
                from .adaptive import run_adaptive
                result = run_adaptive(a.case,resume=a.resume,allow_budget_exceedance=a.allow_budget_exceedance,
                                      progress=lambda message:print(message,file=sys.stderr,flush=True))
            else:
                result = (run.mesh_case(a.case) if a.command == 'mesh' else
                          run.run_case(a.case,resume=a.resume,allow_budget_exceedance=a.allow_budget_exceedance))
        else:
            from . import post
            result = post.analyse(a.case) if a.command == 'analyse' else post.post_case(a.case)
        print(json.dumps(result,indent=2,default=str,allow_nan=False))
        return 0
    except (ValueError,OSError,RuntimeError) as exc:
        print(f'error: {exc}',file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
