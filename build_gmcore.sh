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

export H5DIR="/usr/local/HDF_Group/HDF5/1.14.3"
export CURLDIR=$(spack location -i curl)
export XML2DIR=$(spack location -i libxml2)

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
ninja
