r"""Check a generated extension tree against what the extension promises.

:func:`plumetools.cai2012.verify.verify_tree` already checks the part shared with
the parent study, from the generated OpenFOAM dictionaries: that each case's
weight, time step and schedule match the manifest (and, were there more than
one particle level, that the levels of a Kn differ in ``nEquivalentParticles``
alone). This module adds the extension's own promises, each **recomputed here**
from first principles rather than read back from the code that produced it:

* exactly the explicitly selected Kn suite at the baseline particle level;
  an unmarked historical manifest retains the original four-case policy;
* ``n0`` reproduces the target Kn (:func:`plumetools.cai2012ext.provenance.independent_knudsen`);
* D, S0, T0 and the gas identical in every case;
* ``nEquivalentParticles * multiplier * configured target = n0 V_exit`` -- the
  level applied exactly once, including declared per-Kn occupancy targets;
* sampling starts at 10^4 Kn = 0.01 collision times, recomputed from
  ``lambda_ref = 0.01 D`` and the mean thermal speed, in every case;
* averaging lasts 4.5 domain transits, to within the one write interval the end
  time is snapped by;
* nobody has edited a case's ``case.yaml`` since the manifest was written.
"""

from __future__ import annotations

import math
from pathlib import Path

from plumetools.cai2012 import inflow, mesh
from plumetools.cai2012.config import CaiConfigError, load_case_config
from plumetools.cai2012.geometry import from_config
from plumetools.cai2012.verify import VerificationReport, verify_tree
from plumetools.cai2012ext import (
    PARTICLE_MULTIPLIERS,
    REFERENCE_KNUDSEN,
    policy,
    provenance,
)

#: Agreement required of values recomputed here in double precision.
EXACT_RTOL = 1.0e-9


def _close(a, b, rtol=EXACT_RTOL) -> bool:
    return a is not None and b is not None and math.isclose(
        float(a), float(b), rel_tol=rtol, abs_tol=0.0)


def reference_startup_s(cfg) -> float:
    """``transient_collision_times * lambda_ref / c_bar``, recomputed.

    ``lambda_ref = Kn_ref * D`` by the definition of the Knudsen number, and
    ``c_bar = sqrt(8 k T0 / (pi m))``: neither needs the gas module.
    """
    mean_speed = math.sqrt(8.0 * provenance.K_BOLTZMANN * float(cfg.exit.T0_K)
                           / (math.pi * float(cfg.gas.mass_kg)))
    reference_path = REFERENCE_KNUDSEN * float(cfg.nozzle.diameter_m)
    return float(cfg.dsmc.transient_collision_times) * reference_path / mean_speed


def domain_transit_s(cfg) -> float:
    """``x_max / U0`` with ``U0 = S0 sqrt(2 k T0 / m)``, recomputed."""
    speed = float(cfg.exit.speed_ratio) * math.sqrt(
        2.0 * provenance.K_BOLTZMANN * float(cfg.exit.T0_K)
        / float(cfg.gas.mass_kg))
    return float(cfg.geometry.x_max_over_D) * float(cfg.nozzle.diameter_m) / speed


def verify_extension(root: Path, manifest: dict, *,
                     with_dictionaries: bool = True) -> VerificationReport:
    """Every extension check, plus the parent family's tree check.

    Args:
        root: the family directory, holding ``manifest.yaml``.
        manifest: the loaded manifest.
        with_dictionaries: also run :func:`plumetools.cai2012.verify.verify_tree`,
            which reads the generated dictionaries (``./AllmeshCases
            --dict-only`` writes them without OpenFOAM).

    Returns:
        A :class:`~plumetools.cai2012.verify.VerificationReport`.
    """
    root = Path(root)
    report = VerificationReport()
    entries = manifest.get("cases") or []

    # --- the matrix, and the exclusions ---------------------------------------
    for problem in policy.check_manifest(manifest):
        report.fail("manifest.yaml", problem)

    try:
        allowed_knudsen = policy.suite_knudsen(manifest)
    except policy.ExtensionError as exc:
        report.fail("manifest.yaml", str(exc))
        return report
    expected_count = policy.expected_case_count(manifest)
    if len(entries) != expected_count:
        report.fail("manifest.yaml", f"lists {len(entries)} cases, not "
                                     f"{expected_count}")
    pairs = sorted((float(e.get("Kn", "nan")),
                    float(e.get("numerical_particle_multiplier", "nan")))
                   for e in entries)
    wanted = sorted((k, m) for k in allowed_knudsen for m in PARTICLE_MULTIPLIERS)
    if len(pairs) == len(wanted) and all(
            policy.same_knudsen(a[0], b[0]) and a[1] == b[1]
            for a, b in zip(pairs, wanted)):
        report.ok("manifest.yaml", f"exactly {len(wanted)} cases: Kn "
                                   f"{list(allowed_knudsen)} x "
                                   f"{list(PARTICLE_MULTIPLIERS)}")
    else:
        report.fail("manifest.yaml", f"the (Kn, multiplier) pairs are {pairs}, "
                                     f"not the extension matrix {wanted}")

    cases_dir = root / str(manifest.get("cases_dir", "Cases"))
    if cases_dir.is_dir():
        listed = {Path(e["path"]).name for e in entries}
        for entry in sorted(cases_dir.iterdir()):
            if not entry.is_dir() or entry.name in listed:
                continue
            try:
                cfg = load_case_config(entry)
                kn = cfg.exit.knudsen
            except CaiConfigError:
                kn = None
            if kn is not None and not any(policy.same_knudsen(kn, k) for k in allowed_knudsen):
                report.fail(str(entry.relative_to(root)),
                            f"holds a case at excluded Kn = {kn:g} for this suite. It is not "
                            f"in the manifest and nothing schedules it, but the "
                            f"extension tree must not contain one: remove it.")
            else:
                report.ok(str(entry.relative_to(root)),
                          "not in the manifest; ignored (stray directory)")

    # --- per case, recomputed -------------------------------------------------
    template = None
    for entry in entries:
        name = entry.get("name", "?")
        case_dir = root / entry["path"]
        try:
            cfg = load_case_config(case_dir)
        except CaiConfigError as exc:
            report.fail(name, f"case.yaml does not load: {exc}")
            continue

        digest = provenance.sha256(case_dir / "case.yaml")
        recorded = ((entry.get("sources") or {}).get("case") or {}).get("sha256")
        if digest != recorded:
            report.fail(name, "case.yaml has changed since manifest.yaml was "
                              "written; regenerate (./generate_cases.py) rather "
                              "than editing a generated case")
        else:
            report.ok(name, "case.yaml is the file the manifest describes")

        try:
            case_knudsen = policy.suite_knudsen(cfg)
        except policy.ExtensionError as exc:
            report.fail(name, str(exc))
            continue
        if case_knudsen != allowed_knudsen:
            report.fail(name, "case.yaml and manifest select different Knudsen suites")
        if policy.suite_name(manifest) == policy.FULL_SUITE:
            source_limitations = policy.source_limitations(cfg)
            for where, recorded_limitations in (
                    ("case.yaml meta", cfg.meta.get("known_source_limitations")),
                    ("manifest row", entry.get("known_source_limitations"))):
                if recorded_limitations != source_limitations:
                    report.fail(name, f"{where} does not retain the declared source limitations")
        if not any(policy.same_knudsen(cfg.exit.knudsen, k) for k in allowed_knudsen):
            report.fail(name, f"case.yaml is at excluded Kn = "
                              f"{cfg.exit.knudsen:g}")
            continue
        if not policy.same_knudsen(cfg.exit.knudsen, entry.get("Kn", -1.0)):
            report.fail(name, f"case.yaml says Kn = {cfg.exit.knudsen:g}, the "
                              f"manifest {entry.get('Kn')}")

        geom = from_config(cfg)
        exit_state = inflow.from_config(cfg)
        recomputed = provenance.independent_knudsen(
            exit_state.number_density_per_m3, cfg)
        if _close(recomputed, cfg.exit.knudsen):
            report.ok(name, f"n0 = {exit_state.number_density_per_m3:.6e} "
                            f"reproduces Kn = {recomputed:.9g} (recomputed)")
        else:
            report.fail(name, f"n0 = {exit_state.number_density_per_m3:.6e} "
                              f"implies Kn = {recomputed:.9g}, not "
                              f"{cfg.exit.knudsen:g}")

        invariants = {
            "nozzle.diameter_m": cfg.nozzle.diameter_m,
            "exit.speed_ratio": cfg.exit.speed_ratio,
            "exit.T0_K": cfg.exit.T0_K,
            "exit.characteristic_length": cfg.exit.characteristic_length,
            "gas": cfg.gas,
            "resolution.reference_knudsen": cfg.resolution.reference_knudsen,
            "dsmc.transient_basis": cfg.dsmc.transient_basis,
            "dsmc.transient_collision_times": cfg.dsmc.transient_collision_times,
            "dsmc.sampling_domain_transits": cfg.dsmc.sampling_domain_transits,
            "dsmc.write_interval_s": cfg.dsmc.write_interval_s,
            "dsmc.run_time_multiplier": cfg.dsmc.run_time_multiplier,
            "dsmc.output_frequency_multiplier":
                cfg.dsmc.output_frequency_multiplier,
        }
        if template is None:
            template = (name, invariants)
        else:
            differing = [k for k in invariants if invariants[k] != template[1][k]]
            if differing:
                report.fail(name, f"differs from {template[0]} in {differing}; "
                                  f"these are fixed across the extension")
            else:
                report.ok(name, f"D, S0, T0, gas, startup "
                                f"and averaging as in {template[0]}")

        plan = mesh.plan(cfg, geom, exit_state)
        run = inflow.derive_run_settings(
            cfg, exit_state, geom, min_cell_size_m=plan.min_cell_size_m,
            exit_cell_volume_m3=plan.exit_cell_volume_m3)
        multiplier = float(cfg.dsmc.numerical_particle_multiplier)
        target = float(cfg.resolution.target_particles_per_cell)
        recorded_target = (entry.get("weight") or {}).get("baseline_exit_particles_per_cell")
        if _close(target, recorded_target):
            report.ok(name, f"configured target {target:g} parcels/cell matches its manifest row")
        else:
            report.fail(name, f"configured target {target:g} parcels/cell differs from "
                              f"its manifest row {recorded_target!r}")
        product = (run.n_equivalent_particles * multiplier
                   * float(cfg.resolution.target_particles_per_cell))
        molecules = exit_state.number_density_per_m3 * plan.exit_cell_volume_m3
        if _close(product, molecules) and _close(
                run.n_equivalent_particles,
                float(entry.get("n_equivalent_particles", -1.0))):
            report.ok(name, f"weight {run.n_equivalent_particles:.6e} x "
                            f"{multiplier:g} x "
                            f"{cfg.resolution.target_particles_per_cell:g} = "
                            f"n0 V_exit: the level is applied once")
        else:
            report.fail(name, f"weight {run.n_equivalent_particles:.6e} x "
                              f"{multiplier:g} x target is {product:.6e}, "
                              f"n0 V_exit is {molecules:.6e} (manifest "
                              f"{entry.get('n_equivalent_particles')})")

        start = reference_startup_s(cfg)
        if _close(run.average_start_s, start) and _close(
                entry.get("average_start_s"), start):
            report.ok(name, f"sampling starts at {start:.6e} s = "
                            f"{cfg.dsmc.transient_collision_times:g} Kn = "
                            f"{REFERENCE_KNUDSEN:g} collision times")
        else:
            report.fail(name, f"sampling starts at {run.average_start_s!r} "
                              f"(manifest {entry.get('average_start_s')!r}), "
                              f"but 10^4 reference collision times is {start!r}")

        transit = domain_transit_s(cfg)
        sampled = run.end_time_s - run.average_start_s
        wanted = float(cfg.dsmc.sampling_domain_transits) * transit
        if wanted * (1.0 - EXACT_RTOL) <= sampled <= wanted + run.write_interval_s:
            report.ok(name, f"averages for {sampled / transit:.4f} domain "
                            f"transits (requested "
                            f"{cfg.dsmc.sampling_domain_transits:g}; the end is "
                            f"snapped up to a write)")
        else:
            report.fail(name, f"averages for {sampled:.6e} s = "
                              f"{sampled / transit:.4f} transits, not "
                              f"{cfg.dsmc.sampling_domain_transits:g}")

        status = (entry.get("feasibility") or {}).get("status")
        if status not in ("production-candidate", "unsuitable-for-production"):
            report.fail(name, f"no feasibility verdict in the manifest "
                              f"({status!r})")

    if with_dictionaries:
        tree = verify_tree(root, manifest)
        for message in tree.passes:
            report.passes.append(message)
        for message in tree.failures:
            report.failures.append(message)
    return report
