"""Physical constants and unit conversions.

WARNING -- these are the values the legacy processInflowData.py used, kept to the
same precision. They are deliberately NOT the current CODATA values: substituting
more accurate constants would change every number the source-flow model produces
and break the regression golden. Any update is a physics change and belongs in a
dedicated, reviewed commit.
"""

from __future__ import annotations

#: Boltzmann constant [J/K]. Legacy 5 s.f. value; CODATA 2018 is 1.380649e-23.
BOLTZMANN_J_PER_K = 1.3806e-23

#: Avogadro constant [1/mol]. Legacy 4 s.f. value; CODATA 2018 is 6.02214076e23.
AVOGADRO_PER_MOL = 6.023e23

#: Atomic mass unit [kg]. Legacy value; CODATA 2018 is 1.66053907e-27.
AMU_TO_KG = 1.66605e-27

#: Pounds per square inch to pascals [Pa/psi]. Exact to the digits given.
PSI_TO_PA = 6894.75729
