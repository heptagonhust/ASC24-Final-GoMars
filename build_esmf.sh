#!/bin/bash

set -e
cd "$(dirname $0)" || exit 1

source ./env_xyw.sh

cd gmcore

if [ ! -d esmf ]; then
  # git clone https://github.com/esmf-org/esmf.git
  cp -r /data/gomars_data/esmf ./
fi

cd esmf

pwd

export ESMF_DIR=$(pwd)
export ESMF_COMPILER=intel
export ESMF_COMM=intelmpi
export ESMF_CXXSTD=17
export I_MPI_CC=icx
export I_MPI_CXX=icpx
export I_MPI_F77=ifx

# Dependencies
export NETCDF_PATH=/home/asterich/ASC24-Final-GoMars/gmcore/netcdf
export ESMF_NETCDF=split
export ESMF_NETCDF_INCLUDE=${NETCDF_PATH}/include
export ESMF_NETCDF_LIBPATH=${NETCDF_PATH}/lib
export PIO_PATH=/home/asterich/ASC24-Final-GoMars/gmcore/ParallelIO
export ESMF_PIO=external
export ESMF_PIO_INCLUDE=${PIO_PATH}/include
export ESMF_PIO_LIBPATH=${PIO_PATH}/lib

echo "ESMF_DIR: $ESMF_DIR"

make VERBOSE=1 info
make VERBOSE=1 lib -j64

# install library

ESMF_INSTALL_PREFIX=${ESMF_DIR} \
ESMF_INSTALL_HEADERDIR=include/ \
ESMF_INSTALL_LIBDIR=lib/ \
ESMF_INSTALL_MODDIR=include/ \
ESMF_INSTALL_BINDIR=bin/ \
make install
