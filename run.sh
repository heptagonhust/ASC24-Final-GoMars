#!/bin/bash
#SBATCH -p GPU_nodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=./output/gomars_test.%j.out
base_dir="$(pwd)"

export UCX_RC_PATH_MTU=2048
export I_MPI_HYDRA_RMK=slurm
export I_MPI_PIN=off
export OMP_NUM_THREADS=1
# export I_MPI_STACKSIZE=8388608
# export KMP_STACKSIZE=83886080

spack load intel-oneapi-mkl@2024.0.0
spack load intel-oneapi-mpi@2021.11.0
spack load intel-oneapi-compilers@2024.0.1
spack load hdf5/fxhrrhv

# ulimit -s unlimited

echo "====="
echo "cat $0"
cat $0
echo "====="


# mpirun ./cpu_bind.sh ./bind_test.sh
# mpirun build/gmcore_adv_driver.exe ./run/GMCORE-TESTBED/adv_mv.360x180/namelist 
# mpirun ./build/gmcore_driver.exe ./run/GMCORE-TESTBED/mz.1440x720/namelist
mpirun ./build/gmcore_driver.exe ./run/GMCORE-TESTBED/mz.360x180/namelist
# export LD_LIBRARY_PATH=/opt/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/hdf5-1.14.3-fxhrrhv46hbchtxp255okz3r2dotzmng/lib:$LD_LIBRARY_PATH
# mpirun ./build/gmcore_driver.exe ./run/GMCORE-TESTBED/mz.3600x1800/namelist
# mpirun ./build/gmcore_swm_driver.exe ./run/GMCORE-TESTBED/swm_mz.1440x720/namelist
