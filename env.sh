spack load cmake@3.27.7
spack load ninja@1.11.1
spack load intel-oneapi-compilers@2024.0.2
spack load intel-oneapi-mpi@2021.11.0
spack load intel-oneapi-mkl@2024.0.0
spack load libxml2
spack load curl
spack load hdf5
# spack load hdf5 ~shared

export CC=mpiicx
export FC=mpiifx
export F77=mpiifx
