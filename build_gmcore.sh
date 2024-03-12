#!/bin/bash

set -e
cd "$(dirname $0)" || exit 1

source ./env.sh

cd gmcore

mkdir -p ./lib
if [ x"$(ls -A lib)" = x"" ]; then
  cp -r /data/gomars_libs/gmcore_libs/* ./lib
fi

# ./pull_libs.py

# spack load cmake@3.24.4
# spack load intel-oneapi-compilers@2024.0.1/xbteted 
# spack load intel-oneapi-mkl@2024.0.0
# spack load intel-oneapi-mpi@2021.11.0
# # spack load hdf5 ~shared
# spack load hdf5/fxhrrhv
# export CC=mpiicx
# export FC=mpiifx

export H5DIR=$(spack location -i hdf5)
export CURLDIR=$(spack location -i curl)
export XML2DIR=$(spack location -i libxml2)
export OPENMP_STATIC_LIB_PATH="$(spack location -i intel-oneapi-compilers@2024.0.2)/compiler/latest/lib/libiomp5.a"

echo "CC: $CC"
echo "FC: $FC"
echo "F77: $F77"

export NETCDF_ROOT="$(pwd)/netcdf"
export GPTL_ROOT="$(pwd)/gptl"

if [ x"$1" = xrebuild ]; then
  rm -rf build
fi

if [ ! -d build ]; then
  cmake -B build -G Ninja 
fi
cd build
# make -j8
ninja
