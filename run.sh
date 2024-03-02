#!/bin/bash
#SBATCH -p GPU_nodes
#SBATCH -N 3
#SBATCH -n 12
#SBATCH --output=./output/gomars_test.%j.out
base_dir="$(pwd)"

# you should run this under dir gmcore_main
export UCX_RC_PATH_MTU=2048
export I_MPI_HYDRA_RMK=slurm
export I_MPI_PIN=off
export OMP_NUM_THREADS=1
spack load intel-oneapi-mkl@2024.0.0
spack load intel-oneapi-mpi@2021.11.0
spack load intel-oneapi-compilers@2024.0.1/xbteted

echo "====="
echo "cat $0"
cat $0
echo "====="


# mpirun ./cpu_bind.sh ./bind_test.sh
mpirun build/gmcore_driver.exe ./run/GMCORE-TESTBED/bw.360x180/namelist 