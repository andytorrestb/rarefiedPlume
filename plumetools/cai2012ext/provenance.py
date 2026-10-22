r"""What every generated case was built from, recorded so it can be re-derived.

:func:`plumetools.cai2012.study.manifest_entry` already records Kn, D, lambda0,
T0, S0, U0, n0, the weight, the mesh and the schedule. This adds what an
extension study has to be able to prove afterwards:

* that the density reproduces the target Kn -- **recomputed independently** of
  :mod:`plumetools.cai2012.gas`, from the VHS expression written out again here
  (:func:`independent_knudsen`), so a defect in the module that produced ``n0``
  cannot also pass the audit of it;
* the mass density, pressure, thermal speeds and collision time of the exit
  state;
* the gas, as the numbers that went into ``dsmcProperties``;
* the startup in three units -- Cai's reference collision times (which must be
  10^4), the case's OWN exit collision times, and domain transits -- because the
  extension keeps Cai's Kn = 0.01 reference for every case and a reader should
  be able to see what that means at Kn = 1e-4 and at Kn = 10;
* the averaging window, separately from the startup;
* the source files and the code version.
"""

from __future__ import annotations

import hashlib
import math
import subprocess
from pathlib import Path

#: Boltzmann's constant [J/K], CODATA 2018 exact. Stated here rather than
#: imported: the independent recomputation should share as little as possible
#: with the code it audits.
K_BOLTZMANN = 1.380649e-23


def sha256(path: Path) -> str | None:
    """SHA-256 of a file, or ``None`` if it does not exist."""
    path = Path(path)
    if not path.is_file():
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_provenance(repo_root: Path) -> dict:
    """The code version: commit, whether the tree was dirty, and a describe.

    Never raises. A machine without git, or an exported tree, records
    ``"unavailable"`` -- the manifest still says what was generated, only not
    from which commit.
    """
    def git(*args) -> str | None:
        try:
            result = subprocess.run(
                ["git", "-C", str(repo_root), *args], capture_output=True,
                text=True, timeout=30, check=False)
        except (OSError, subprocess.SubprocessError):
            return None
        return result.stdout.strip() if result.returncode == 0 else None

    commit = git("rev-parse", "HEAD")
    if commit is None:
        return {"commit": "unavailable", "dirty": None, "describe": None}
    status = git("status", "--porcelain", "--untracked-files=no")
    return {
        "commit": commit,
        "dirty": bool(status) if status is not None else None,
        "describe": git("describe", "--always", "--dirty"),
    }


def independent_knudsen(number_density_per_m3: float, cfg) -> float:
    r"""The Knudsen number ``n0`` implies, recomputed from the gas section.

    .. math::

        \lambda = \frac{(T_0/T_{ref})^{\omega - 1/2}}{\sqrt2\,\pi\,d_{ref}^2\,n_0}
        \quad (\text{VHS}),\qquad
        \lambda = \frac{1}{\sqrt2\,\pi\,d_{ref}^2\,n_0}\quad(\text{hard sphere}),
        \qquad \mathrm{Kn} = \lambda / L

    Written out here on purpose rather than calling
    :func:`plumetools.cai2012.gas.knudsen_from_number_density`: that module
    produced ``n0``, and an audit that called it back would agree with it
    whatever it got wrong.
    """
    gas = cfg.gas
    n = float(number_density_per_m3)
    hard_sphere = 1.0 / (math.sqrt(2.0) * math.pi * float(gas.diameter_m) ** 2 * n)
    if gas.mean_free_path_convention == "vhs":
        factor = (float(cfg.exit.T0_K) / float(gas.t_ref_K)) ** (
            float(gas.omega) - 0.5)
    else:
        factor = 1.0
    length = float(cfg.nozzle.diameter_m)
    if cfg.exit.characteristic_length == "radius":
        length *= 0.5
    return hard_sphere * factor / length


def exit_conditions(cfg, geom, exit_state) -> dict:
    """The exit state, every quantity the extension promises to record."""
    gas = cfg.gas
    t0 = float(exit_state.T0_K)
    n0 = float(exit_state.number_density_per_m3)
    mass = float(gas.mass_kg)
    gas_constant = K_BOLTZMANN / mass
    most_probable = math.sqrt(2.0 * gas_constant * t0)
    recomputed = independent_knudsen(n0, cfg)
    target = float(cfg.exit.knudsen)
    return {
        "Kn_target": target,
        "Kn_recomputed": recomputed,
        "Kn_relative_deviation": (recomputed - target) / target,
        "Kn_recomputed_by": "plumetools.cai2012ext.provenance.independent_knudsen",
        "characteristic_length": cfg.exit.characteristic_length,
        "characteristic_length_m": float(exit_state.characteristic_length_m),
        "D_m": float(geom.diameter_m),
        "R0_m": float(geom.radius_m),
        "S0": float(exit_state.speed_ratio),
        "S0_recomputed": float(exit_state.velocity_m_per_s) / most_probable,
        "T0_K": t0,
        "U0_m_per_s": float(exit_state.velocity_m_per_s),
        "most_probable_speed_m_per_s": most_probable,
        "mean_thermal_speed_m_per_s": float(exit_state.mean_thermal_speed_m_per_s),
        "beta0_s2_per_m2": float(exit_state.beta0),
        "n0_per_m3": n0,
        "rho0_kg_per_m3": mass * n0,
        "p0_Pa": n0 * K_BOLTZMANN * t0,
        "lambda0_m": float(exit_state.mean_free_path_m),
        "exit_collision_time_s": float(exit_state.collision_time_s()),
        "number_flux_per_m2_s": float(exit_state.number_flux_per_m2_s),
        "mass_flux_kg_per_m2_s": mass * float(exit_state.number_flux_per_m2_s),
        "mean_free_path_convention": exit_state.convention,
    }


def gas_parameters(cfg) -> dict:
    """The gas exactly as ``case.yaml`` states it, plus its gas constant."""
    gas = cfg.gas
    return {
        "species_name": gas.species_name,
        "mass_kg": float(gas.mass_kg),
        "diameter_m": float(gas.diameter_m),
        "omega": float(gas.omega),
        "t_ref_K": float(gas.t_ref_K),
        "internal_degrees_of_freedom": int(gas.internal_degrees_of_freedom),
        "mean_free_path_convention": gas.mean_free_path_convention,
        "gas_constant_J_per_kg_K": K_BOLTZMANN / float(gas.mass_kg),
        "binary_collision_model": cfg.dsmc.binary_collision_model,
        "collisions_enabled": bool(cfg.dsmc.collisions_enabled),
    }


def startup(cfg, run, exit_state) -> dict:
    """When sampling starts, and in which units that is a statement about what.

    ``average_start_s`` is ``transient_collision_times`` reference collision
    times, the reference being the Kn = ``resolution.reference_knudsen`` exit
    state. That is the same number of seconds in every case; what it is
    worth in each case's own collision times varies by five decades, and both
    are recorded so the choice is visible rather than implied.
    """
    reference_time = float(run.reference_collision_time_s)
    local_time = float(exit_state.collision_time_s())
    return {
        "basis": run.transient_basis,
        "average_start_s": float(run.average_start_s),
        "reference_knudsen": float(cfg.resolution.reference_knudsen),
        "reference_mean_free_path_m": float(run.reference_mean_free_path_m),
        "reference_collision_time_s": reference_time,
        "requested_reference_collision_times": float(
            cfg.dsmc.transient_collision_times),
        "achieved_reference_collision_times": float(run.average_start_s)
        / reference_time,
        "exit_collision_times": float(run.average_start_s) / local_time,
        "domain_transits": float(run.average_start_s) / float(run.domain_transit_s),
        "domain_transit_s": float(run.domain_transit_s),
        "note": (f"Cai's 10^4 t0 with t0 the Kn = "
                 f"{float(cfg.resolution.reference_knudsen):g} exit collision "
                 f"time, held for every case. No case is generated at the "
                 f"reference Kn."),
    }


def averaging(cfg, run) -> dict:
    """The sampling window, the run length and the write schedule."""
    requested_interval = cfg.dsmc.write_interval_s
    return {
        "average_start_s": float(run.average_start_s),
        "end_time_s": float(run.end_time_s),
        "sampling_time_s": float(run.sampling_time_s),
        "requested_sampling_domain_transits": float(
            cfg.dsmc.sampling_domain_transits),
        "achieved_sampling_domain_transits": float(run.sampling_time_s)
        / float(run.domain_transit_s),
        "run_time_multiplier": float(run.run_time_multiplier),
        "requested_write_interval_s": (None if requested_interval is None
                                       else float(requested_interval)),
        "achieved_write_interval_s": float(run.write_interval_s),
        "write_interval_steps": int(run.write_interval_steps),
        "n_steps": int(run.n_steps),
        "n_writes": int(run.n_writes),
        "n_sampled_writes": int(run.n_sampled_writes),
        "note": ("endTime is snapped up to a whole number of write intervals, "
                 "so the achieved averaging is at most one interval longer "
                 "than requested; writeControl timeStep rounds the requested "
                 "interval to whole steps."),
    }


def particle_weight(cfg, run) -> dict:
    """The particle level, applied exactly once."""
    multiplier = float(run.numerical_particle_multiplier)
    target = float(cfg.resolution.target_particles_per_cell)
    return {
        "numerical_particle_multiplier": multiplier,
        "baseline_n_equivalent_particles": float(
            run.baseline_n_equivalent_particles),
        "n_equivalent_particles": float(run.n_equivalent_particles),
        "baseline_exit_particles_per_cell": target,
        "exit_particles_per_cell": float(run.exit_particles_per_cell),
        "weight_source": run.weight_source,
        "rule": "nEquivalentParticles = (n0 * V_exit / "
                "target_particles_per_cell) / numerical_particle_multiplier",
    }


def sources(root: Path, case_dir: Path, template_dir: Path,
            study_path: Path) -> dict:
    """Digests of every file the case was built from."""
    root = Path(root)

    def rel(path: Path) -> str:
        try:
            return str(Path(path).resolve().relative_to(root.resolve())).replace(
                "\\", "/")
        except ValueError:
            return str(path)

    return {
        "study": {"path": rel(study_path), "sha256": sha256(study_path)},
        "template": {"path": rel(template_dir / "case.yaml"),
                     "sha256": sha256(Path(template_dir) / "case.yaml")},
        "case": {"path": rel(Path(case_dir) / "case.yaml"),
                 "sha256": sha256(Path(case_dir) / "case.yaml")},
    }
