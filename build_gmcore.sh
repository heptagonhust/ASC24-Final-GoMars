#!/bin/bash

cd gmcore/

./pull_libs.py

spack load cmake@3.24.4
spack load intel-oneapi-compilers@2024.0.1
spack load intel-oneapi-mkl@2024.0.0

spack load intel-oneapi-mpi@2021.11.0

export CC=mpiicx
# export CXX=mpiicpx
export FC=mpiifx
# can also be mpiifx

# spack load netcdf-c
# spack load netcdf-fortran

# export NETCDF_ROOT=$(spack location -i netcdf-fortran)
current_dir=$(pwd)
export NETCDF_ROOT=$current_dir/netcdf
# export GPTL_ROOT=/opt/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/gptl-8.1.1-3hh72awwp2aabzqetd6yloj7xtjvuhu6
export GPTL_ROOT=$current_dir/gptl
# export NETCDF_ROOT=$(spack location -i netcdf-fortran)
# export NETCDF_ROOT=/data/asc24caeporo/xyw/finalasc/raw/netcdf/testf
# export NETCDF_ROOT=/data/asc24caeporo/xyw/finalasc/raw/netcdf/testc
export LD_LIBRARY_PATH=/opt/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/hdf5-1.14.3-fxhrrhv46hbchtxp255okz3r2dotzmng/lib:$LD_LIBRARY_PATH


# export FC=gfortran

# export CFLAGS=${CFLAGS:=}
# export FFLAGS=${FFLAGS:=}
# export CXXFLAGS=${CXXFLAGS:=}

# export CC=icx
# export CXX=icpx
# export FC=ifx

# add_compile_flag() {
# 	export CFLAGS="$1 $CFLAGS"
# 	export CXXFLAGS="$1 $CXXFLAGS"
# 	export FFLAGS="$1 $FFLAGS"
# }

# openmpi_base="/usr/mpi/gcc/openmpi-4.0.3rc4/"
# export PKG_CONFIG_PATH="${openmpi_base}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
# openmpi_FLAGS=$(pkg-config --cflags ompi)
# openmpi_LIBS="$(pkg-config --static --libs ompi)"
# add_compile_flag "${openmpi_FLAGS}"
# export LDFLAGS="${openmpi_LIBS} ${LDFLAGS:-}"

rm -rf build
cmake -B build
cd build
make -j8
