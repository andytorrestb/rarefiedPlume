"""Writers for OpenFOAM dictionaries and field files."""

from __future__ import annotations

from pathlib import Path

from plumetools.foamio.fields import write_inflow_fields

__all__ = ["write_inflow_fields", "write_mesh_setup", "mesh_pipeline"]


def write_mesh_setup(case_dir: Path, cfg) -> list[Path]:
    """Write whichever mesh dictionaries ``mesh.type`` calls for.

    Args:
        case_dir: an OpenFOAM case directory.
        cfg: a :class:`plumetools.config.CaseConfig`.

    Returns:
        The paths written.

    Raises:
        NotImplementedError: for an unknown ``mesh.type``.

    Dispatching here keeps ``Allmesh`` from having to know which module renders
    which mesh type.
    """
    # Imported lazily so a case using one generator need not parse the other.
    if cfg.mesh.type == "block_mesh_ogrid":
        from plumetools.foamio.blockmesh import write_block_mesh_dict
        return [write_block_mesh_dict(case_dir, cfg)]
    if cfg.mesh.type == "snappy_hex_sphere":
        from plumetools.foamio.snappy import write_snappy_setup
        return write_snappy_setup(case_dir, cfg)
    if cfg.mesh.type == "snappy_markelov":
        from plumetools.markelov1999.mesh import write_mesh_setup as write_markelov
        return write_markelov(case_dir, cfg)
    raise NotImplementedError(f"mesh.type {cfg.mesh.type!r} is not implemented")


def mesh_pipeline(cfg) -> list[str]:
    """The mesh commands to run, in order, for this ``mesh.type``.

    Returned as bare command names so ``Allmesh`` can check them against PATH and
    run them without duplicating the knowledge of which pipeline goes with which
    generator.
    """
    if cfg.mesh.type == "block_mesh_ogrid":
        return ["blockMesh"]
    if cfg.mesh.type in ("snappy_hex_sphere", "snappy_markelov"):
        return ["blockMesh", "snappyHexMesh"]
    raise NotImplementedError(f"mesh.type {cfg.mesh.type!r} is not implemented")
