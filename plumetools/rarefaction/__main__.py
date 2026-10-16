"""``python -m plumetools.rarefaction CASE --settings FILE``: any dsmcFoam case.

A nozzle-plume case does not need this -- ``./Allpost`` runs the diagnostics
with everything taken from its ``case.yaml``. This is the door for every other
case: the gas, the reference lengths and the thresholds come from a small YAML
file instead. See :mod:`plumetools.rarefaction` for its shape.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

from plumetools.nozzleplume.gas import MEAN_FREE_PATH_CONVENTIONS, VhsSpecies
from plumetools.rarefaction import foam, summary
from plumetools.rarefaction.run import diagnose
from plumetools.rarefaction.settings import DiagnosticSettings, SettingsError

_GAS_KEYS = {"molecular_mass_kg", "molecular_diameter_m", "omega", "t_ref_K",
             "mean_free_path_convention", "species_name"}


def load_settings_file(path: Path) -> tuple:
    """``(species, convention, settings, background_cell_size_m)`` from a file.

    Raises:
        SettingsError: naming the key at fault.
    """
    path = Path(path)
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise SettingsError(f"{path}: {exc}") from None
    unknown = sorted(set(raw) - {"gas", "diagnostics", "background_cell_size_m"})
    if unknown:
        raise SettingsError(f"{path}: unknown key(s) {unknown}; valid: "
                            f"['background_cell_size_m', 'diagnostics', 'gas']")
    gas = raw.get("gas") or {}
    stray = sorted(set(gas) - _GAS_KEYS)
    if stray:
        raise SettingsError(f"{path}:gas: unknown key(s) {stray}; valid: "
                            f"{sorted(_GAS_KEYS)}")
    missing = sorted({"molecular_mass_kg", "molecular_diameter_m"} - set(gas))
    if missing:
        raise SettingsError(f"{path}:gas: missing {missing}")
    convention = str(gas.get("mean_free_path_convention", "vhs"))
    if convention not in MEAN_FREE_PATH_CONVENTIONS:
        raise SettingsError(f"{path}:gas: mean_free_path_convention must be one "
                            f"of {list(MEAN_FREE_PATH_CONVENTIONS)}")
    try:
        species = VhsSpecies(
            name=str(gas.get("species_name", "gas")),
            mass_kg=float(gas["molecular_mass_kg"]),
            diameter_m=float(gas["molecular_diameter_m"]),
            omega=float(gas.get("omega", 0.5)),
            t_ref_K=float(gas.get("t_ref_K", 273.0)))
    except (TypeError, ValueError) as exc:
        raise SettingsError(f"{path}:gas: {exc}") from None
    settings = DiagnosticSettings.from_mapping(raw.get("diagnostics"),
                                               f"{path}:diagnostics")
    background = raw.get("background_cell_size_m")
    return species, convention, settings, (None if background is None
                                           else float(background))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m plumetools.rarefaction",
        description="Local mean free path, Knudsen numbers and mesh resolution "
                    "of a sampled dsmcFoam case. Writes the fields into the time "
                    "directory and a summary YAML.")
    parser.add_argument("case", help="case directory (reconstructed)")
    parser.add_argument("--settings", required=True, metavar="FILE",
                        help="gas, reference lengths and thresholds (YAML)")
    parser.add_argument("--time", default=None,
                        help="time directory to read (default: the latest)")
    parser.add_argument("--summary", default=None, metavar="FILE",
                        help="where to write the summary (default: "
                             "<case>/results/diagnostics_summary.yaml)")
    args = parser.parse_args(argv)

    case_dir = Path(args.case)
    try:
        species, convention, settings, background = load_settings_file(args.settings)
        if not settings.enabled:
            print("diagnostics.enabled is false; nothing to do")
            return 0
        diagnosis = diagnose(case_dir, species, settings, convention=convention,
                             time=args.time, background_cell_size_m=background)
    except (SettingsError, foam.PostError, ValueError) as exc:
        print(f"rarefaction: {exc}", file=sys.stderr)
        return 1

    print("\n".join(summary.report(diagnosis.document)))
    path = summary.write_yaml(
        Path(args.summary) if args.summary
        else case_dir / "results" / "diagnostics_summary.yaml", diagnosis.document)
    print(f"    wrote {len(diagnosis.written)} fields into "
          f"{case_dir / diagnosis.data.time}")
    print(f"    wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
