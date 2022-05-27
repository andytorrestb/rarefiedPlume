#!/bin/bash

dsmcInitialise+
python3 processInflowData.py 743.397 802.2867 6.0201e20
decomposePar
mpirun -np 8 dsmcFoam+
postProcess -func sampleDict
python3 graphCaseValidation.py