spack load cmake@3.24.4
spack load intel-oneapi-compilers@2024.0.1/xbteted
spack load intel-oneapi-mpi@2021.11.0+envmods~external-libfabric+generic-names~ilp64
spack load intel-oneapi-mkl@2024.0.0
spack load libxml2/q66mtbb
spack load curl
spack load hdf5/fxhrrhv
# spack load hdf5 ~shared

export CC=mpiicx
export FC=mpiifx
export F77=mpiifx
