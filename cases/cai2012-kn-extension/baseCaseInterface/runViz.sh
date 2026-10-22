#PBS -S /bin/bash
#PBS -q shortq
#PBS -l walltime=0:10:00
#PBS -l nodes=1:ppn=1
#PBS -N knx-iviz
#
# INTERFACE MODE: the family's imagery (../../AllpostCases --viz-only) reads
# ../../manifest.yaml, which lists only the direct-exit cases, and its viz
# specs assume the exit is in the mesh. Not adapted for interface cases in this
# increment; this job refuses rather than drawing a misleading picture.

printf 'FAILED: imagery is not adapted for interface cases (see the comment in runViz.sh)\n' >&2
exit 1
