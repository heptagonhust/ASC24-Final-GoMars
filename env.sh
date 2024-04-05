spack load cmake@3.27.7
spack load ninja@1.11.1
spack load libxml2
spack load curl
# spack load hdf5 # note intel-llvm-mpi inside
# spack load hdf5 ~shared

# spack load openmpi
spack load nvhpc

export OPAL_PREFIX=$(spack location -i nvhpc)/Linux_x86_64/23.9/comm_libs/mpi

# spack load gcc@12

# export CXXFLAGS="-fPIC -shared" FCFLAGS="-fPIC -shared"

# export CC=mpicc
# export FC=mpifort
# export F77=mpifort

export CC=mpicc
export FC=mpifort
export F77=mpifort

echo "============CC INFO============="
$CC --version
echo "==========END CC INFO==========="
