"""Analytical plume source-flow inflow generation for dsmcFoam cases.

Extracted from the per-case ``processInflowData.py`` scripts so that the model
lives in one tested place instead of being copied between case directories.

Targets OpenFOAM's own ``dsmcFoam``, verified against v2512. The MNF fork's
``dsmcFoam+`` is supported through ``output.dialect: mnf`` -- the two read
``0/boundaryT`` as different types -- but is not required by any active case;
see docs/solver-compatibility.md.

Layout:
    constants   physical constants (legacy precision -- see the module docstring)
    species     gas property table
    config      case.yaml loading and validation
    mesh/       OpenFOAM polyMesh parsing
    geometry    centroids, normals, spherical coordinates
    sourceflow  the analytical model, equations E1-E8
    foamio/     OpenFOAM dictionary and field writers
    inflow      the orchestration that ties the above together
"""

__version__ = "0.1.0"
