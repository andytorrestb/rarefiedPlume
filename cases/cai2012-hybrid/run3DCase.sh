#!/usr/bin/env bash
#PBS -S /bin/bash
#PBS -q longq
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=192
#PBS -l mem=1200gb
#PBS -N hybrid-3d-ellipsoid

# Submit from the repository root or cases/cai2012-hybrid. PBS copies this
# script, so locate the submitted checkout through PBS_O_WORKDIR.
set -eo pipefail
SUBMIT_DIR=${PBS_O_WORKDIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}
SUBMIT_DIR=$(cd -- "$SUBMIT_DIR" && pwd)
if [ -f "$SUBMIT_DIR/cases/cai2012-hybrid/Allrun" ]; then
    SCRIPT_DIR="$SUBMIT_DIR/cases/cai2012-hybrid"
elif [ -f "$SUBMIT_DIR/Allrun" ] && [ -f "$SUBMIT_DIR/hybrid.yaml" ]; then
    SCRIPT_DIR="$SUBMIT_DIR"
else
    printf 'run3DCase.sh: submit from the repository root or cases/cai2012-hybrid\n' >&2
    exit 1
fi
REPO=$(cd -- "$SCRIPT_DIR/../.." && pwd)
PROFILE=${PROFILE:-hpc-kn1em3-ellipsoid}
CASE_DIR=${1:-${CASE_DIR:-"$SCRIPT_DIR/runs/$PROFILE"}}
[[ "$CASE_DIR" = /* ]] || CASE_DIR="$SUBMIT_DIR/$CASE_DIR"
STAGE=${STAGE:-all}
HPC_ENV=${HYBRID_HPC_ENV:-"$REPO/cases/cai2012-kn-extension/hpc_env.sh"}
cd "$REPO"

# This shared environment guards OpenFOAM's optional startup probes itself.
# Keep errexit for module/venv failures; never enable nounset around Lmod/FOAM.
source "$HPC_ENV"
: "${FOAM_USER_LIBBIN:?HPC environment did not define FOAM_USER_LIBBIN}"
: "${FOAM_APPBIN:?HPC environment did not define FOAM_APPBIN}"
# Always run the incremental build: an existing but stale library (e.g. a
# libhybridPlumeThermo.so without cappedTabulated) would otherwise only fail
# at thermo selection after the job has queued. Allwmake is a no-op when the
# sources are unchanged; a build failure aborts the job (errexit).
for app in thermophysicalModels dsmcBoundaryModels; do
    bash "$REPO/applications/$app/Allwmake"
done
export PATH="$FOAM_APPBIN:$PATH"
args=("$CASE_DIR" "$STAGE" --config "$SCRIPT_DIR/hybrid.yaml" --profile "$PROFILE" -np 192)
if [ "${ALLOW_UNSUITABLE:-0}" = 1 ]; then
    args+=(--allow-unsuitable)
fi
printf '>>> 3D profile %s, stage %s, 192 ranks, case %s\n' "$PROFILE" "$STAGE" "$CASE_DIR"
bash "$SCRIPT_DIR/Allrun" "${args[@]}"
