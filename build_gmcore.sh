#!/bin/bash

set -e
cd "$(dirname $0)" || exit 1

source ./env_xyw.sh

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
export OMPDIR="$(spack location -i intel-oneapi-compilers@2024.0.2)/compiler/latest/"
# source /data/spack/opt/spack/linux-ubuntu22.04-icelake/gcc-11.4.0/intel-oneapi-compilers-2024.0.2-lvfe6ufintzu3ibq3loire4oz62soeqe/setvars.sh

echo "CC: $CC"
echo "FC: $FC"
echo "F77: $F77"

export NETCDF_ROOT="$(pwd)/netcdf"
export GPTL_ROOT="$(pwd)/gptl"
export PIO_ROOT="$(pwd)/ParallelIO"
export ESMF_ROOT="$(pwd)/esmf"
export MCT_ROOT="$(pwd)/MCT"
export FOX_ROOT="$(pwd)/fox"

# get lapack from mkl
export LAPACK_ROOT=$(spack location -i intel-oneapi-mkl)/lib

echo "NETCDF_ROOT: $NETCDF_ROOT"
echo "GPTL_ROOT: $GPTL_ROOT"
echo "PIO_ROOT: $PIO_ROOT"
echo "LAPACK_ROOT: $LAPACK_ROOT"
echo "ESMF_ROOT: $ESMF_ROOT"
echo "MCT_ROOT: $MCT_ROOT"
echo "FOX_ROOT: $FOX_ROOT"


if [ x"$1" = xrebuild ]; then
  rm -rf build
fi

if [ ! -d build ]; then
  cmake -B build -DCMAKE_BUILD_TYPE=relwithdebinfo \
    -DCMAKE_RANLIB=/data/spack/opt/spack/linux-ubuntu22.04-icelake/gcc-11.4.0/intel-oneapi-compilers-2024.0.2-lvfe6ufintzu3ibq3loire4oz62soeqe/compiler/2024.0/bin/compiler/llvm-ranlib
fi
cd build
make -j64 VERBOSE=1
# ninja
