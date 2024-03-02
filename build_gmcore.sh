#!/bin/bash

cd gmcore/

./pull_libs.py

spack load cmake@3.24.4
spack load intel-oneapi-compilers@2024.0.1
spack load intel-oneapi-mkl@2024.0.0

spack load intel-oneapi-mpi@2021.11.0
# spack load hdf5 ~shared
spack load hdf5/fxhrrhv

export CC=mpiicx
export FC=mpiifx


current_dir=$(pwd)
export NETCDF_ROOT=$current_dir/netcdf
export GPTL_ROOT=$current_dir/gptl


rm -rf build
cmake -B build
cd build
make -j8
