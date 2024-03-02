#!/bin/bash

cd gmcore_main/

./pull_libs.py

spack load cmake@3.24.4
spack load intel-oneapi-compilers@2024.0.1/xbteted
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
# export NETCDF_ROOT=$(spack location -i netcdf-fortran)
# export NETCDF_ROOT=/data/asc24caeporo/xyw/finalasc/raw/netcdf/testf
# export NETCDF_ROOT=/data/asc24caeporo/xyw/finalasc/raw/netcdf/testc
export LD_LIBRARY_PATH=/opt/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/hdf5-1.14.3-fxhrrhv46hbchtxp255okz3r2dotzmng/lib:$LD_LIBRARY_PATH


rm -rf build
cmake -B build
cd build
make -j8
