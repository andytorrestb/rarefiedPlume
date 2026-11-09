"""Native ellipsoid cavity dictionaries, isolated from DSMC checkpoints."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import math

from plumetools.foamio.primitives import (
    SearchableSphere, SearchableCone, RefinementRegion, BoxPatch, fmt_point,
    render_box_block_mesh_dict, render_snappy_dict, render_mesh_quality_dict)
from . import mesh as cartesian


@dataclass(frozen=True)
class _Spheroid(SearchableSphere):
    transverse: float = 1.0

    def geometry_lines(self):
        return ["        type searchableSphere;", f"        centre {fmt_point(self.centre)};",
                f"        radius {fmt_point((self.radius, self.transverse, self.transverse))};"]


def write_mesh_inputs(case_dir, cfg, source, fit):
    case_dir = Path(case_dir)
    system = case_dir/"mesh-build/system"
    system.mkdir(parents=True, exist_ok=True)
    d = source.diameter_m
    length = float(cfg.require("domain.x_max_over_D"))*d
    width = float(cfg.require("domain.half_width_over_D"))*d
    requested = float(cfg.get("dsmc.mesh.background_cell_over_D", .2))*d
    divisions = tuple(max(2, math.ceil(span/requested)) for span in (length, 2*width, 2*width))
    size = max(length/divisions[0], 2*width/divisions[1], 2*width/divisions[2])
    source_size = float(cfg.get("dsmc.mesh.source_cell_over_mfp", .8))*source.mean_free_path
    level = max(0, math.ceil(math.log2(size/source_size)))
    surface_level = int(cfg.get("dsmc.mesh.surface_level", 3))
    spheroid = _Spheroid(name="interface", level=surface_level, radius=fit["a"],
                        transverse=fit["b"], centre=(fit["xc"], 0, 0), patch_type="patch")
    regions = []
    # Nested thin annuli find the source at coarse levels without refining an
    # entire near-source volume to the exit mean free path.
    for k in range(1, level+1):
        parent = size/2**(k-1)
        inner = max(0, fit["source_radius_m"]-2*parent)
        region = SearchableCone(name=f"sourceRefinement{k}", point1=(-parent, 0, 0),
                                point2=(4*parent, 0, 0), radius1=source.radius_m+2*parent,
                                radius2=source.radius_m+2*parent,
                                inner_radius1=inner, inner_radius2=inner)
        regions.append(RefinementRegion(region, k))
    seed = (length-.5*length/divisions[0], width-.5*2*width/divisions[1],
            width-1.5*2*width/divisions[2])
    if ((seed[0]-fit["xc"])/fit["a"])**2+(seed[1]**2+seed[2]**2)/fit["b"]**2 <= 1:
        raise ValueError("background seed is inside the fitted ellipsoid")
    block = render_box_block_mesh_dict(
        (0, length, -width, width, -width, width), divisions,
        [BoxPatch("exitPlane", "symmetryPlane", ("x_min",), "meshing constraint; converted to open patch before DSMC"),
         BoxPatch("farField", "patch", ("x_max", "y_min", "y_max", "z_min", "z_max"), "open exterior")],
        generator=__name__, note="3D DSMC background; not a CFD subset")
    snappy = render_snappy_dict([spheroid], location_in_mesh=seed, n_cells_between_levels=3,
                                generator=__name__, note="fitted 3D ellipsoid and source refinement",
                                max_global_cells=int(cfg.get("dsmc.mesh.max_global_cells", 12000000)),
                                refinement_regions=regions)
    # The primitive renderer defaults to half the global ceiling per rank.
    import re
    snappy = re.sub(r"maxLocalCells\s+\d+;", "maxLocalCells "
                    +str(int(cfg.get("dsmc.mesh.max_local_cells", 1000000)))+";", snappy)
    control = ("FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }\n"
               "application snappyHexMesh; startFrom startTime; startTime 0; stopAt endTime; "
               "endTime 1; deltaT 1; writeControl timeStep; writeInterval 1; "
               "writeFormat ascii; writePrecision 12;\n")
    plan = cartesian.plan(cfg, source)
    values = {"blockMeshDict": block, "snappyHexMeshDict": snappy,
              "meshQualityDict": render_mesh_quality_dict(__name__),
              "topoSetDict": cartesian.render_topo_set_dict(plan, source),
              "createPatchDict": cartesian.render_create_patch_dict(), "controlDict": control,
              "fvSchemes": 'FoamFile {version 2.0; format ascii; class dictionary; object fvSchemes;}\n'
                           'ddtSchemes {default none;} gradSchemes {default none;} divSchemes {default none;} '
                           'laplacianSchemes {default none;} interpolationSchemes {default none;} snGradSchemes {default none;}\n',
              "fvSolution": 'FoamFile {version 2.0; format ascii; class dictionary; object fvSolution;}\n'}
    for name, text in values.items():
        (system/name).write_text(text, encoding="utf-8", newline="\n")
    return [system/name for name in values]
