#PBS -S /bin/bash
#PBS -q shortq
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -m abe
#PBS -N viz-r-4d

set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO_ROOT/cases/cai2012-kn-extension/hpc_env.sh"

# printf '>>> Starting: run Allmesh\n'
# ./Allmesh
# printf '<<< Completed: run Allmesh\n'
# bash ./Allrun -np 192 > output.log 2> error.log
# ./Allpost > output_post.log 2> error_post.log
./Allviz > output_viz.log 2> error_viz.log
