"""Verified single-species elastic OpenFOAM v2512 VHS collision kernel.

Provenance: installed OpenFOAM v2512 source,
``src/lagrangian/DSMC/submodels/BinaryCollisionModel/VariableHardSphere/``
``VariableHardSphere.C``, SHA256
``d43f948867a526f0452a2aabe5ffdb610a2db318ae427c600b83dc1429246731``.
The sigmaTcR method includes Gamma(5/2-omega); collide samples an isotropic
relative direction. This is not inferred from the equilibrium mean free path.
The installed solver's ``etc/controlDict`` defines k=1.38065e-23 J/K.
That value is retained here to match sigmaTcR; legacy Cai source moments use
the exact SI value 1.380649e-23. Their relative difference is 7.24e-7.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from plumetools.cai2012.gas import ARGON, VhsSpecies

OPENFOAM_V2512_BOLTZMANN_J_PER_K = 1.38065e-23


class UnsupportedCollisionLaw(ValueError):
    """The requested law has not been verified for this reduced model."""


@dataclass(frozen=True)
class VHSKernel:
    """Equal-species VHS sigma(g)*g [m^3/s], with isotropic scattering.

    Only the verified standard ``VariableHardSphere`` model is accepted.
    No collision-strength fitting parameter is exposed.
    """

    species: VhsSpecies = ARGON
    law: str = "VariableHardSphere"
    boltzmann_J_K: float = OPENFOAM_V2512_BOLTZMANN_J_PER_K

    def __post_init__(self):
        if self.law != "VariableHardSphere":
            raise UnsupportedCollisionLaw(self.law)
        if self.species.internal_degrees_of_freedom != 0:
            raise UnsupportedCollisionLaw("internal energy exchange is unsupported")
        if not math.isfinite(self.boltzmann_J_K) or self.boltzmann_J_K <= 0:
            raise ValueError("Boltzmann constant must be finite and positive")

    def rate(self, relative_speed):
        """Return sigma(g)*g in m^3/s for nonnegative g [m/s].

        At exactly zero relative speed return zero, as the verified solver
        does, including omega=1's discontinuous zero-speed convention.
        """
        g = np.asarray(relative_speed, dtype=float)
        if np.any(~np.isfinite(g)) or np.any(g < 0):
            raise ValueError("relative speed must be finite and nonnegative")
        sp = self.species
        exponent = sp.omega - 0.5
        coefficient = (math.pi * sp.diameter_m**2
                       * (4 * self.boltzmann_J_K * sp.t_ref_K / sp.mass_kg)**exponent
                       / math.gamma(2.5 - sp.omega))
        result = coefficient * np.power(g, 2 - 2 * sp.omega)
        return np.where(g > 0, result, 0.0)


@dataclass(frozen=True)
class ZeroKernel:
    """Explicit collision-free numerical limit; not a fitted physical law."""

    species: VhsSpecies = ARGON
    law: str = "disabled"

    def rate(self, relative_speed):
        g = np.asarray(relative_speed, dtype=float)
        if np.any(~np.isfinite(g)) or np.any(g < 0):
            raise ValueError("relative speed must be finite and nonnegative")
        return np.zeros_like(g)


def elastic_pair(velocity, partner_velocity, direction):
    """Scatter identical elastic particles into the supplied unit direction.

    ``direction`` is an outgoing relative-velocity direction. The transport
    sampler must supply isotropic directions for the supported VHS law.
    Vectors may be broadcast over leading dimensions; final dimension is 3.
    The direction is checked, never silently normalized.
    """
    v, w, d = np.broadcast_arrays(np.asarray(velocity, dtype=float),
                                   np.asarray(partner_velocity, dtype=float),
                                   np.asarray(direction, dtype=float))
    if v.shape[-1:] != (3,) or np.any(~np.isfinite(v + w + d)):
        raise ValueError("finite Cartesian 3-vectors required")
    if np.any(np.abs(np.linalg.norm(d, axis=-1) - 1) > 1e-12):
        raise ValueError("outgoing direction must have unit norm within 1e-12")
    centre = (v + w) / 2
    relative = np.linalg.norm(v - w, axis=-1)[..., None] * d / 2
    return centre + relative, centre - relative
