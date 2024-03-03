#!/bin/bash
# spack load cmake@3.24.4
# spack load intel-oneapi-compilers@2024.0.1/xbteted
# spack load intel-oneapi-mpi@2021.11.0
# # spack load libxml2
# # spack load curl
# # spack load hdf5/fxhrrhv
# export CC=mpiicx
# export FC=mpiifort
# export F77=mpiifort
source ./env.sh

current_dir=$(pwd)
target_dir="$current_dir/gmcore/gptl"

if [ ! -d "$target_dir" ]; then
    mkdir -p "$target_dir"
fi

pushd gmcore
pushd gptl
CDIR=$(pwd)
# wget https://github.com/jmrosinski/GPTL/releases/download/v8.1.1/gptl-8.1.1.tar.gz
# tar -zxvf gptl-8.1.1.tar.gz

cp -r /data/gomars_data/gptl-8.1.1 .
pushd gptl-8.1.1

# wget https://gist.githubusercontent.com/bonfus/21dec6b966859f5f509b935f8b055a7f/raw/macros.make

./configure --enable-pmpi --prefix=${CDIR}

# make check

make install

popd
popd
popd