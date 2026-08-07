"""Physical constants and unit conversions for the AIAA 99-3455 case family.

These are **CODATA 2018 / SI-2019 exact** values, deliberately different from
:mod:`plumetools.constants`.

That module carries the 4-5 significant-figure values the legacy
``processInflowData.py`` used, and its docstring says why they cannot be updated:
the regression golden in ``tests/regression`` is bit-exact against them, so
substituting accurate constants there would break a freeze that exists to prove
the legacy extraction was behaviour-preserving.

This case family has no such freeze. It is new work with no golden to protect, so
it uses the accurate values. The two modules must therefore never be mixed in one
calculation -- the difference is about 5e-5 relative on ``k_B`` and 3.6e-4 on the
molecular mass, which is small but is not round-off, and tracing a discrepancy
back to a mis-imported constant is unpleasant.

Nothing here is a value printed by the paper. The paper's own inputs live in
:mod:`plumetools.markelov1999.geometry` and in each case's ``case.yaml``.
"""

from __future__ import annotations

#: Boltzmann constant [J/K]. Exact by the 2019 SI redefinition.
BOLTZMANN_J_PER_K = 1.380649e-23

#: Avogadro constant [1/mol]. Exact by the 2019 SI redefinition.
AVOGADRO_PER_MOL = 6.02214076e23

#: Molar gas constant [J/(mol K)] = k_B * N_A. Exact.
GAS_CONSTANT_J_PER_MOL_K = BOLTZMANN_J_PER_K * AVOGADRO_PER_MOL

#: Pounds per square inch to pascals [Pa/psi].
#:
#: Exact: 1 lbf = 4.4482216152605 N and 1 in = 0.0254 m, both defined, so
#: 1 psi = 4.4482216152605 / 0.0254**2 Pa. The legacy module rounds this to
#: 6894.75729, which is 4.2e-8 relative low -- irrelevant physically, but it means
#: a legacy 475 psi and a new 475 psi differ in the last few digits.
PSI_TO_PA = 4.4482216152605 / (0.0254 ** 2)

#: Inches to metres [m/in]. Exact by definition.
INCH_TO_M = 0.0254


def psi_to_pa(psi: float) -> float:
    """Convert a reservoir pressure from psi to Pa."""
    return float(psi) * PSI_TO_PA


def inch_to_m(inches: float) -> float:
    """Convert a length from inches to metres."""
    return float(inches) * INCH_TO_M


def molecular_mass_kg(molar_mass_g_per_mol: float) -> float:
    """Mass of one molecule [kg] from the molar mass [g/mol].

    ``m = M / (1000 * N_A)``. This is the *only* mass conversion this package
    uses; the legacy model instead multiplies the molar mass by an ``amu``
    constant, which is the same thing to within the precision of its constants
    but is a second, independent path to the same number.

    For N2 (28.0134 g/mol) this gives 4.651735e-26 kg.
    """
    if molar_mass_g_per_mol <= 0.0:
        raise ValueError(f"molar mass must be positive, got {molar_mass_g_per_mol}")
    return float(molar_mass_g_per_mol) / 1000.0 / AVOGADRO_PER_MOL
