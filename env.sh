
#================Intel compilers + OpenMPI(BAD)===================
# spack load cmake@3.27.7
# spack load ninja@1.11.1
# spack load libxml2
# spack load curl
# spack load hdf5 # introduces intel mpi
# spack load openmpi
# #spack load intel-oneapi-compilers@2024.0.2
# export CC=mpicc
# export FC=mpif90
# export F77=mpif77
#================intel compilers + OpenMPI===================

#================GCC + OpenMPI===================
# spack load gcc@12
# spack load cmake@3.27.7
# spack load ninja@1.11.1
# spack load libxml2
# spack load curl
# spack load hdf5 # introduces intel mpi
# spack load openmpi
# export CC=mpicc
# export FC=mpif90
# export F77=mpif77
#================GCC + OpenMPI===================

#================GCC + Intel MPI(BAD)===================
# spack load gcc@12 #which is bad and not support intel mpi, and doesn't have nvptx-none
spack load cmake@3.27.7
spack load ninja@1.11.1
spack load libxml2
spack load curl
spack load hdf5 # introduces intel mpi
spack load intel-oneapi-mpi@2021.11.0
export CC=mpicc
export FC=mpif90
export F77=mpif77
#================GCC + Intel MPI===================

#================Intel compiler + Intel MPI===================
# spack load cmake@3.27.7
# spack load ninja@1.11.1
# spack load libxml2
# spack load curl
# spack load hdf5 # introduces intel mpi
# spack load intel-oneapi-mpi@2021.11.0
# #spack load intel-oneapi-compilers@2023.2.1 #can't build gptl (xfortcom: Unknown command line argument)
# #spack load intel-oneapi-compilers@2024.0.2
# export CC=mpiicx
# export FC=mpiifx
# export F77=mpiifx
#================Intel compiler + Intel MPI===================


# $FC -v
echo '==========CC INFO=========='
which $CC
which gcc
$CC -v
$CC --version
echo '==========FC INFO=========='
which $FC
which gfortran
$FC -v
$FC --version
echo '===================='