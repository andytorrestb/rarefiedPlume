"""Gas property table.

Provided so a case can name a species instead of restating its constants. Note
that ``cases/3d-inflow`` currently overrides this from ``case.yaml`` with N2
values while its solver simulates argon -- see :data:`SPECIES` below and finding
SM-03.
"""

from __future__ import annotations

from dataclasses import dataclass

from plumetools.constants import AMU_TO_KG


@dataclass(frozen=True)
class Species:
    """Properties of one gas species.

    Attributes:
        name: dsmcFoam+ type id, matching ``typeIdList`` in ``dsmcProperties``.
        gamma: ratio of specific heats [-].
        molar_mass_g_per_mol: molar mass [g/mol].
        mass_kg: mass of one molecule [kg].
        diameter_m: VHS reference diameter [m].
        omega: VHS viscosity-temperature exponent [-].
        alpha: VSS scattering parameter [-].
        rotational_dof: rotational degrees of freedom [-].
    """

    name: str
    gamma: float
    molar_mass_g_per_mol: float
    mass_kg: float
    diameter_m: float
    omega: float
    alpha: float
    rotational_dof: int


#: Known species. Values for argon match ``constant/dsmcProperties`` in every case
#: in this repository.
#:
#: SM-03: the source-flow model in ``cases/3d-inflow`` is evaluated with the *N2*
#: entry (gamma 1.4, M_w 28.0134) while ``dsmcProperties`` declares ``typeIdList
#: (Ar)`` with mass 6.63e-26 and ``rotationalDegreesOfFreedom 0`` (so gamma = 5/3),
#: and the generated field is named ``boundaryNumberDensity_Ar``. At 300 K that is
#: a limiting velocity of 788.16 m/s instead of 558.89 m/s. The mismatch is
#: preserved for now via explicit ``gas:`` overrides in ``case.yaml``; correcting
#: it is a reviewed physics change, not part of the extraction.
SPECIES: dict[str, Species] = {
    "Ar": Species(
        name="Ar",
        gamma=5.0 / 3.0,
        molar_mass_g_per_mol=39.948,
        mass_kg=6.63e-26,
        diameter_m=4.17e-10,
        omega=0.74,
        alpha=1.36,
        rotational_dof=0,
    ),
    "N2": Species(
        name="N2",
        gamma=1.4,
        molar_mass_g_per_mol=28.0134,
        mass_kg=28.0134 * AMU_TO_KG,
        diameter_m=4.17e-10,
        omega=0.74,
        alpha=1.0,
        rotational_dof=2,
    ),
}


def get(name: str) -> Species:
    """Look up a species by name, or raise listing what is known."""
    try:
        return SPECIES[name]
    except KeyError:
        raise KeyError(f"unknown species {name!r}; known: {sorted(SPECIES)}") from None
