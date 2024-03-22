#!/bin/bash

cd gmcore/

./pull_libs.py

spack load cmake@3.24.4
spack load intel-oneapi-compilers@2024.0.1
spack load intel-oneapi-mkl@2024.0.0

spack load intel-oneapi-mpi@2021.11.0

current_dir=$(pwd)
export NETCDF_ROOT=$current_dir/netcdf
export GPTL_ROOT=$current_dir/gptl
export LD_LIBRARY_PATH=/opt/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/hdf5-1.14.3-fxhrrhv46hbchtxp255okz3r2dotzmng/lib:$LD_LIBRARY_PATH

rm -rf build
cmake -B build
cd build
make -j8
