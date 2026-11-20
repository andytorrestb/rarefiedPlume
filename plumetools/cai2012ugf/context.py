"""What this family's ``Allcompare.py`` / ``Allconverge.py`` need of a case.

The kn-extension scripts call ``interface_case.post_context``, whose verdict is
the dsmcFoam pre-run label -- the one spec §7 rules out as a gate. This
family's copies import this module under that name instead (their only
difference from the originals), so they show the post-run verdict, or
``post-run-verdict-pending`` before ``Allverdict.py`` has run.
"""

from __future__ import annotations

from pathlib import Path

from plumetools.cai2012ext.interface_case import PostContext
from plumetools.cai2012ugf.verdict import read_verdict


def post_context(case_dir: Path) -> PostContext:
    """The direct-exit derivation with the post-run verdict."""
    from plumetools.cai2012ext import direct_case

    case_dir = Path(case_dir)
    case = direct_case.derive(case_dir)
    full = read_verdict(case_dir)
    verdict = {"status": full["status"],
               "production_ready": bool(full["production_ready"]),
               "unresolved": list(full["unresolved"])}
    return PostContext("direct", case.cfg, case.geom, case.exit_state, case.run,
                       verdict, float(case.nozzle.area_error), None, None)
