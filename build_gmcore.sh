#!/bin/bash
source ./env.sh

cd gmcore/

if [ ! -d "lib" ]; then
	mkdir -p ./lib
	cp -r /data/gomars_libs/gmcore_libs/* ./lib
fi
./pull_libs.py

# spack load cmake@3.24.4
# spack load intel-oneapi-compilers@2024.0.1/xbteted 
# spack load intel-oneapi-mkl@2024.0.0
# spack load intel-oneapi-mpi@2021.11.0
# # spack load hdf5 ~shared
# spack load hdf5/fxhrrhv
# export CC=mpiicx
# export FC=mpiifx
echo "CC: $CC"
echo "FC: $FC"
echo "F77: $F77"

current_dir=$(pwd)
export NETCDF_ROOT=$current_dir/netcdf
export GPTL_ROOT=$current_dir/gptl

target_dir="$current_dir/build"

if [ ! -d "$target_dir" ]; then
	rm -rf build
	cmake -B build -G Ninja 
	echo "hello"
fi
cd build
# make -j8
ninja
