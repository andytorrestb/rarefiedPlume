# cases/cai2012-kn-extension/hpc_env.sh -- the cluster environment, in one place.
#
# SOURCED by every job script of this family (baseCase/run*.sh, copied into each
# case, and runMeshAll.sh / runAnalyse.sh here); never executed. Each script
# finds it by looking in ./, ../ and ../../ from $PBS_O_WORKDIR.
#
# Never source it under `set -u`: Lmod's and OpenFOAM's start-up files read
# variables that are not set. The job scripts use `set -eo pipefail`.

if [[ $- == *u* ]]; then
    printf 'FAILED: hpc_env.sh was sourced under set -u; Lmod and OpenFOAM read unset variables\n' >&2
    return 1
fi

printf '>>> Starting: load Spack environment\n'
source /home2/appman/spack/share/spack/setup-env.sh
printf '<<< Completed: load Spack environment\n'

printf '>>> Starting: initialize Lmod\n'
source "$(spack location -i lmod)/lmod/lmod/init/bash"
printf '<<< Completed: initialize Lmod\n'

printf '>>> Starting: purge modules\n'
module purge
printf '<<< Completed: purge modules\n'
printf '>>> Starting: add Spack module path\n'
module use /apps/spack/share/spack/modules/linux-rocky8-zen4
printf '<<< Completed: add Spack module path\n'
printf '>>> Starting: load Python module\n'
module load python/3.11.7-gcc-8.5.0-wfpoppf
printf '<<< Completed: load Python module\n'
printf '>>> Starting: load OpenMPI module\n'
module load openmpi/5.0.2-gcc-13.2.0-4pbtkdf
printf '<<< Completed: load OpenMPI module\n'
printf '>>> Starting: load FFmpeg module\n'
module load ffmpeg/6.1.1-gcc-13.2.0-75kx2tp
printf '<<< Completed: load FFmpeg module\n'

# A file sourced without arguments of its own sees its caller's, and OpenFOAM's
# bashrc reads file-like ones as preference files. Hand it none, then put the
# caller's back.
printf '>>> Starting: load OpenFOAM environment\n'
hpc_env_args=("$@")
set --
if source /home2/torresfa/scratch/OpenFOAM/OpenFOAM-v2512/etc/bashrc; then
    set -- "${hpc_env_args[@]}"
    unset hpc_env_args
    printf '<<< Completed: load OpenFOAM environment\n'
else
    status=$?
    printf 'FAILED: loading OpenFOAM environment returned exit %d\n' "$status" >&2
    exit "$status"
fi

printf '>>> Starting: activate Python virtual environment\n'
source /home2/torresfa/scratch/OpenFOAM/rarefiedPlume/.venv/bin/activate
printf '<<< Completed: activate Python virtual environment\n'

printf '>>> Starting: configure ParaView environment\n'
export PVPYTHON_PATH="${PVPYTHON_PATH:-/home2/torresfa/src/paraview/ParaView-6.2.0-MPI-Linux-Python3.12-x86_64/bin/pvpython}"
export VTK_DEFAULT_OPENGL_WINDOW=vtkOSOpenGLRenderWindow
printf '<<< Completed: configure ParaView environment\n'
