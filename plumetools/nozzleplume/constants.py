"""Physical constants for the nozzle-plume family.

CODATA 2018 / SI-2019 **exact** values, restated rather than imported, for the
reason :mod:`plumetools.cai2012.constants` gives: :mod:`plumetools.constants`
carries legacy 4-5 significant-figure values a regression golden is bit-exact
against, and importing another case family's constants to get ``k_B`` would
couple two packages that share no code.
"""

from __future__ import annotations

#: Boltzmann constant [J/K]. Exact by the 2019 SI redefinition.
BOLTZMANN_J_PER_K = 1.380649e-23

#: Avogadro constant [1/mol]. Exact by the 2019 SI redefinition.
AVOGADRO_PER_MOL = 6.02214076e23

#: Standard gravity [m/s^2]. Exact by definition; used only to quote a specific
#: impulse next to the thrust, which is how an engine data sheet states it.
STANDARD_GRAVITY_M_PER_S2 = 9.80665
