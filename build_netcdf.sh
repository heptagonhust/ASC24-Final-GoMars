#!/bin/bash

set -e
cd "$(dirname $0)" || exit 1

source ./env.sh

# export H5DIR=$(spack location -i hdf5 ~shared)
export H5DIR=$(spack location -i hdf5)
export CURLDIR=$(spack location -i curl)
export XML2DIR=$(spack location -i libxml2)

(
  cd gmcore
  netcdf_dir="netcdf"
  mkdir -p "$netcdf_dir"
  cd "$netcdf_dir"
  if [ x"$(ls -A .)" = x"" ]; then
    # wget https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
    # wget https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.1/netcdf-fortran-4.6.1.tar.gz
    # tar -zxvf netcdf-c-4.9.2.tar.gz
    # tar -zxvf netcdf-fortran-4.6.1.tar.gz
    cp -r /data/gomars_data/netcdf/netcdf-c-4.9.2/ .
    cp -r /data/gomars_data/netcdf/netcdf-fortran-4.6.1/ .
  fi
)

cd ./gmcore/netcdf/

CDIR=$(pwd)
NCDIR=$(pwd)
NFDIR=$(pwd)

(
  cd ./netcdf-c-4.9.2
  [ -f Makefile ] || [ -f makefile ] && make clean
  CPPFLAGS="-I${H5DIR}/include -I${CURLDIR}/include -I${XML2DIR}/include " \
    LDFLAGS="-L${H5DIR}/lib -L${CURLDIR}/lib -L${XML2DIR}/lib" \
    ./configure --prefix=${CDIR} --enable-shared=yes --enable-static=yes
  make install -j64
)

export LD_LIBRARY_PATH=${NCDIR}/lib:${LD_LIBRARY_PATH}
(
  cd ./netcdf-fortran-4.6.1
  [ -f Makefile ] || [ -f makefile ] && make clean
  CPPFLAGS="-I${NCDIR}/include -I${H5DIR}/include -I${CURLDIR}/include -I${XML2DIR}/include" \
    LDFLAGS="-L${NCDIR}/lib -L${H5DIR}/lib -L${CURLDIR}/lib -L${XML2DIR}/lib" \
    ./configure --prefix=${NFDIR} --enable-shared=no --enable-static=yes
  make install -j64
)

# wordaround: remove all shared libs
# rm ./lib/*so*
