#PBS -S /bin/bash
#PBS -q shortq
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N ugf-meshall
#
# Generate, validate and mesh all seven cases, in one job. Submit from here:
#
#     cd cases/cai2012-unigas && qsub runMeshAll.sh
#
#   1. python generate_cases.py      Cases/ and manifest-unigas.yaml
#   2. python validate_cases.py      the cases are the kn-extension's
#   3. ./Allmesh in every case       dictionaries, blockMesh, nozzle, checkMesh
#
# Any failure fails the job. Then one job per case: cd Cases/<case> && qsub
# runCase.sh. To re-mesh a single case: Cases/<case>/runMesh.sh.
set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with qsub from cases/cai2012-unigas\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
[ -f study.yaml ] && [ -f generate_cases.py ] \
    || { printf 'FAILED: %s is not the family directory\n' "$PWD" >&2; exit 1; }
source ./hpc_env.sh

printf '>>> Starting: generate the cases\n'
python generate_cases.py
printf '<<< Completed: generate the cases\n'

printf '>>> Starting: validate against the kn-extension\n'
python validate_cases.py
printf '<<< Completed: validate against the kn-extension\n'

mapfile -t CASES < <(python -c "import yaml
for c in yaml.safe_load(open('manifest-unigas.yaml'))['cases']: print(c['path'])")
[ "${#CASES[@]}" -gt 0 ] || { printf 'FAILED: the manifest lists no cases\n' >&2; exit 1; }
for case in "${CASES[@]}"; do
    printf '>>> Starting: mesh %s\n' "$case"
    (cd "$case" && bash ./Allmesh)
    printf '<<< Completed: mesh %s\n' "$case"
done
