"""uniGasFoam dictionaries for a nozzleplume plate case (spec §2, §4, §6).

``plumetools.cai2012ugf``'s renderers read a cai2012-schema configuration.
Neither package may be edited, so :func:`ugf_views` hands them light views of a
nozzleplume :class:`~plumetools.nozzleplume.workflow.Case` carrying exactly the
attributes they read. ``system/boundariesDict`` is rendered here: it adds the
plate, a ``uniGasDiffuseWallPatch``, which also makes ``uniGasVolFields`` write
the plate loads (``p_``, ``surfaceShearStress_``, ``fD_``,
``surfaceHeatTransfer_``). Symmetry patches get no entry: OpenFOAM's particle
tracking reflects specularly at them (``uniGasParcel::hitPatch`` returns false).
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from plumetools.cai2012ugf import dictionaries as ugf
from plumetools.foamio.primitives import write_lines
from plumetools.nozzleplume import dictionaries as npd
from plumetools.nozzleplume import mesh as npmesh

#: Sampling-leg writes per domain transit (spec §6).
SAMPLING_WRITES_PER_TRANSIT = 4


def ugf_views(case):
    """``(cfg, exit_state, run)`` views of ``case`` for the cai2012ugf renderers."""
    cfg, es, run = case.cfg, case.exit_state, case.run
    gas = cfg.gas
    cfg_view = SimpleNamespace(
        gas=SimpleNamespace(
            species_name=str(gas.species_name),
            t_ref_K=float(gas.t_ref_K),
            mass_kg=float(gas.molecular_mass_kg),
            diameter_m=float(gas.molecular_diameter_m),
            internal_degrees_of_freedom=int(gas.internal_degrees_of_freedom),
            omega=float(gas.omega),
        ),
        geometry=SimpleNamespace(plume_axis=(1.0, 0.0, 0.0)),
        mesh=SimpleNamespace(patch_names=dict(cfg.mesh.patch_names)),
        dsmc=SimpleNamespace(
            initial_number_density_per_m3=float(cfg.dsmc.initial_number_density_per_m3),
            initial_temperature_K=float(cfg.dsmc.initial_temperature_K),
            n_subdomains=int(cfg.dsmc.n_subdomains),
        ),
    )
    exit_view = SimpleNamespace(
        T0_K=float(es.temperature_K),
        number_density_per_m3=float(es.number_density_per_m3),
        velocity_m_per_s=float(es.speed_m_per_s),
    )
    run_view = SimpleNamespace(
        n_equivalent_particles=float(run.n_equivalent_particles),
        average_start_s=float(run.average_start_s),
        end_time_s=float(run.end_time_s),
        delta_t_s=float(run.delta_t_s),
        domain_transit_s=float(run.transit_s),
        write_interval_s=float(run.transit_s) / SAMPLING_WRITES_PER_TRANSIT,
    )
    return cfg_view, exit_view, run_view


def render_boundaries_dict(case, settings) -> str:
    """``system/boundariesDict``: inject on the nozzle, reflect diffusely at the plate,
    delete at every open patch; nothing on the symmetry planes."""
    cfg_view, exit_view, _ = ugf_views(case)
    names = case.cfg.mesh.patch_names
    sp = cfg_view.gas.species_name
    t0 = ugf._num(exit_view.T0_K)
    plate_t = ugf._num(case.cfg.plate.temperature_K)
    lines = ugf._header(
        "boundariesDict",
        "Cai's uniform drifting Maxwellian on the nozzle disk; the plate is a "
        "fully diffuse wall at its temperature; every open patch deletes the "
        "particles that cross it; the symmetry planes reflect specularly "
        "(OpenFOAM's own tracking, so they have no entry).",
        "system")
    lines += ["uniGasPatchBoundaries", "("]
    for patch in (names["nozzle"], names["upstream_vacuum"], names["outer"]):
        lines += ugf._patch_boundary(patch, "uniGasDeletionPatch", ["allSpecies        yes;"])
    lines += ugf._patch_boundary(names["plate"], "uniGasDiffuseWallPatch",
                                 ["velocity          (0 0 0);",
                                  f"temperature       {plate_t};"])
    lines += [");", "", "uniGasGeneralBoundaries", "("]
    lines += [
        "    boundary",
        "    {",
        "        generalBoundaryProperties",
        "        {",
        f"            patch       {names['nozzle']};",
        "        }",
        "        boundaryModel   uniGasFreeStreamInflowPatch;",
        "        uniGasFreeStreamInflowPatchProperties",
        "        {",
        f"            typeIds                     ({sp});",
        f"            numberDensities             {{{sp}    {ugf._num(exit_view.number_density_per_m3)};}}",
        f"            translationalTemperature    {t0};",
        f"            rotationalTemperature       {t0};",
        f"            vibrationalTemperature      {t0};",
        f"            electronicTemperature       {t0};",
        f"            velocity                    {ugf._axis_velocity(cfg_view, exit_view.velocity_m_per_s)};",
        "        }",
        "    }",
        ");",
        "",
        "uniGasCyclicBoundaries",
        "(",
        ");",
    ]
    return ugf._text(lines)


def render_decompose_par_dict(cfg) -> str:
    """nozzleplume's decomposeParDict, balancing ``decompositionWeight``
    (:mod:`plumetools.plateugf.weights`) instead of cells."""
    from plumetools.plateugf.weights import FIELD
    out = []
    for line in npd.render_decompose_par_dict(cfg).splitlines():
        out.append(line)
        if line.startswith("method"):
            out += ["",
                    "// Balance the expected particles (+1 per cell), not the cells: the",
                    "// collisionless plume-plus-plate-gas count per cell, written to 0/ by",
                    "// Allmesh (python -m plumetools.plateugf weights).",
                    f"weightField         {FIELD};"]
    return "\n".join(out) + "\n"


def write_case_dictionaries(case, settings) -> list[Path]:
    """Write the mesh setup and every uniGasFoam dictionary of ``case``."""
    case_dir = Path(case.case_dir)
    cfg_view, _, run_view = ugf_views(case)
    system, constant = case_dir / "system", case_dir / "constant"
    written = list(npmesh.write_mesh_setup(case_dir, case.cfg, case.geom, case.plan))
    startup = ugf.render_control_dict(cfg_view, run_view, settings, "startup")
    texts = {
        constant / "uniGasProperties": ugf.render_uni_gas_properties(cfg_view, run_view, settings),
        system / "boundariesDict": render_boundaries_dict(case, settings),
        system / "fieldPropertiesDict": ugf.render_field_properties_dict(cfg_view, run_view, settings),
        system / "uniGasInitialisationDict": ugf.render_initialisation_dict(cfg_view),
        system / "hybridDecompositionDict": ugf.render_hybrid_decomposition_dict(run_view, settings),
        system / "chemReactDict": ugf.render_chem_react_dict(),
        system / "controlDict.startup": startup,
        system / "controlDict.sampling": ugf.render_control_dict(cfg_view, run_view, settings, "sampling"),
        system / "controlDict": startup,
        system / "decomposeParDict": render_decompose_par_dict(case.cfg),
        system / "fvSchemes": ugf.render_fv_schemes(),
        system / "fvSolution": ugf.render_fv_solution(),
    }
    for path, text in texts.items():
        written.append(write_lines(path, text.splitlines()))
    (case_dir / "0").mkdir(parents=True, exist_ok=True)
    return written
