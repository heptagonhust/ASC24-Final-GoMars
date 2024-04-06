#!/bin/bash

set -e
cd "$(dirname $0)" || exit 1

source ./env_xyw.sh

cd gmcore

if [ ! -d MCT ]; then
  # git clone https://github.com/MCSclimate/MCT.git
  cp -r /data/gomars_data/MCT ./
fi

cd MCT

pwd

export MPILIBS="-L${I_MPI_ROOT}/lib"
export MPIHEADER="-I${I_MPI_ROOT}/include"
export MPIFC=mpiifx
export FC=mpiifx
export CC=mpiicx
export CXX=mpiicpx
export OPT=-O3

./configure --prefix=$(pwd)

make clean

make -j64

make install
