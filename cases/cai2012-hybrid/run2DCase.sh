#!/usr/bin/env bash
# Run the planar case from this checkout, regardless of the caller's directory.
set -eo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd -- "$SCRIPT_DIR/../.." && pwd)
cd "$REPO"

# OpenFOAM probes optional packages; do not apply errexit inside its setup.
if ! source /usr/lib/openfoam/openfoam2512/etc/bashrc; then
    printf 'run2DCase.sh: failed to initialize OpenFOAM v2512\n' >&2
    exit 1
fi
if [ -f "$REPO/.venv/bin/activate" ]; then
    source "$REPO/.venv/bin/activate"
fi

C="$SCRIPT_DIR/runs/planar-laptop-ellipse"
R="$SCRIPT_DIR/Allrun"

if [ ! -f "$C/hybrid.resolved.yaml" ]; then
    bash "$R" "$C" setup --profile planar-laptop-kn0p01
fi

bash "$R" "$C" cfd
bash "$R" "$C" post-cfd
bash "$R" "$C" prepare-dsmc
