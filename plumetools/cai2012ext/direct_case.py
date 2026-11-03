"""One direct-exit case of the extension: everything it runs at, derived once.

The counterpart of :func:`plumetools.cai2012ext.interface_case.derive` for a
case that injects at Cai's physical exit. Generation, the manifest, the PBS
steps, the post-processing context and the scaling inventory all derive a case
here, so the mesh they plan, the weight they write and the verdict they print
come from one sequence::

    case.yaml -> geometry -> exit state -> mesh plan -> run settings
              -> predicted nozzle -> parcel estimate -> DSMC checks

The calls are the parent family's (:mod:`plumetools.cai2012`), in the order
``cases/cai2012/generate_cases.py`` makes them. Nothing here writes a file.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from plumetools.cai2012 import checks, inflow, mesh
from plumetools.cai2012.checks import CheckReport
from plumetools.cai2012.config import CaiCaseConfig, load_case_config
from plumetools.cai2012.geometry import CaiGeometry, from_config
from plumetools.cai2012.inflow import ExitState, RunSettings
from plumetools.cai2012.mesh import MeshPlan, NozzlePrediction


@dataclass(frozen=True)
class DirectCase:
    """Everything a direct-exit case runs at.

    Attributes:
        cfg: the validated ``case.yaml``.
        geom: the domain and nozzle geometry [m].
        exit_state: the uniform exit state (``n0`` [1/m^3], ``U0`` [m/s],
            ``T0`` [K]).
        plan: the graded block mesh, planned, not written.
        run: time step, end time, averaging window, particle weight.
        nozzle: the staircase nozzle patch ``createPatch`` will carve.
        particles: the steady-state parcel estimate
            (:func:`plumetools.cai2012.mesh.estimate_particles`).
        report: the DSMC quality checks against all of the above.
        case_dir: where the case lives, or ``None`` when derived from a config
            alone.

    Iterating gives ``(cfg, geom, exit_state, plan, run, nozzle, particles,
    report)`` -- the tuple :func:`plumetools.cai2012ext.generation.derive`
    used to return -- for code written against it. Use the attributes.
    """

    cfg: CaiCaseConfig
    geom: CaiGeometry
    exit_state: ExitState
    plan: MeshPlan
    run: RunSettings
    nozzle: NozzlePrediction
    particles: dict
    report: CheckReport
    case_dir: Path | None = None

    def __iter__(self):
        return iter((self.cfg, self.geom, self.exit_state, self.plan, self.run,
                     self.nozzle, self.particles, self.report))


def derive_config(cfg: CaiCaseConfig, case_dir: Path | None = None) -> DirectCase:
    """Derive a direct-exit case from an already loaded configuration."""
    geom = from_config(cfg)
    exit_state = inflow.from_config(cfg)
    plan = mesh.plan(cfg, geom, exit_state)
    run = inflow.derive_run_settings(
        cfg, exit_state, geom,
        min_cell_size_m=plan.min_cell_size_m,
        exit_cell_volume_m3=plan.exit_cell_volume_m3)
    nozzle = mesh.predicted_nozzle_faces(plan, geom)
    particles = mesh.estimate_particles(
        plan, geom, exit_state, run.n_equivalent_particles)
    report = checks.run_checks(cfg, geom, exit_state, plan, run,
                               nozzle=nozzle, particles=particles)
    return DirectCase(cfg=cfg, geom=geom, exit_state=exit_state, plan=plan,
                      run=run, nozzle=nozzle, particles=particles, report=report,
                      case_dir=None if case_dir is None else Path(case_dir))


def derive(case_dir: Path) -> DirectCase:
    """Derive a direct-exit case from its own ``case.yaml``."""
    return derive_config(load_case_config(case_dir), case_dir)
