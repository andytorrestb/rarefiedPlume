#!/bin/bash
# Queue the plate-unigas cases on PBS: each case's mesh job, then its run job
# once the mesh job has succeeded (afterok). Run from cases/plate-unigas.
#
#   bash submit_all.sh                    every case under Cases/
#   bash submit_all.sh Kn100_zD25 ...     only these
#
# A case that already holds results is skipped: resume it with
#   cd Cases/<case> && qsub runCase.sh
set -uo pipefail
cd "$(dirname "$0")" || exit 1
command -v qsub >/dev/null || { echo "submit_all: qsub is not on PATH" >&2; exit 1; }

if [ $# -gt 0 ]; then CASES=("$@"); else CASES=($(ls Cases)); fi
for case in "${CASES[@]}"; do
    dir=Cases/$case
    [ -f "$dir/case.yaml" ] || { echo "skip $case: no $dir/case.yaml" >&2; continue; }
    if ls -d "$dir"/processor0/0.[0-9]* "$dir"/processor0/[1-9]* "$dir"/0.[0-9]* "$dir"/[1-9]* >/dev/null 2>&1; then
        echo "skip $case: it holds results (resume with: cd $dir && qsub runCase.sh)"
        continue
    fi
    mesh=$(cd "$dir" && qsub -N "pm-$case" runMesh.sh) || { echo "FAILED to queue $case mesh" >&2; continue; }
    run=$(cd "$dir" && qsub -N "pr-$case" -W depend=afterok:"$mesh" runCase.sh) \
        || { echo "FAILED to queue $case run (mesh job $mesh is queued)" >&2; continue; }
    echo "$case: mesh $mesh, run $run"
done
