"""Physical constants for the Cai & Wang 2012 case family.

CODATA 2018 / SI-2019 **exact** values.

Deliberately not imported from elsewhere in the repository:

* :mod:`plumetools.constants` carries the 4-5 significant-figure values the
  legacy ``processInflowData.py`` used and cannot be updated, because
  ``tests/regression`` is bit-exact against them.
* :mod:`plumetools.markelov1999.constants` has the right values, but it also
  carries that paper's psi and inch conversions. Importing another case family's
  constants module to get ``k_B`` would couple two families that share no
  physics, and the coupling would be invisible from either side.

Two exact defined constants are cheaper to restate than to entangle.
"""

from __future__ import annotations

#: Boltzmann constant [J/K]. Exact by the 2019 SI redefinition.
BOLTZMANN_J_PER_K = 1.380649e-23

#: Avogadro constant [1/mol]. Exact by the 2019 SI redefinition.
AVOGADRO_PER_MOL = 6.02214076e23
