#!/bin/bash

set -e

# spack load cmake@3.24.4
# spack load intel-oneapi-compilers@2024.0.1/xbteted
# spack load intel-oneapi-mkl@2024.0.0
# spack load intel-oneapi-mpi@2021.11.0
# spack load libxml2/q66mtbb
# spack load curl
# # spack load hdf5 ~shared
# spack load hdf5/fxhrrhv
# export CC=mpiicx
# export FC=mpiifort
# export F77=mpiifort

source ./env.sh

export CC=mpiicx
export FC=mpiifort
export F77=mpiifort


# export H5DIR=$(spack location -i hdf5 ~shared)
export H5DIR=$(spack location -i hdf5)
export CURLDIR=$(spack location -i curl)
export XML2DIR=$(spack location -i libxml2)

target_dir="netcdf"
pushd gmcore
if [ ! -d "$target_dir" ]; then
    mkdir "$target_dir"
    echo "created: $target_dir"
    pushd netcdf
    # wget https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
    # wget https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.1/netcdf-fortran-4.6.1.tar.gz
    # tar -zxvf netcdf-c-4.9.2.tar.gz
    # tar -zxvf netcdf-fortran-4.6.1.tar.gz
    cp -r /data/gomars_data/netcdf/netcdf-c-4.9.2/ .
    cp -r /data/gomars_data/netcdf/netcdf-fortran-4.6.1/ .
    popd
    popd
    # ./build_netcdf.sh
    # I copied the content of gomars_all_in_one.sh to this file and here it caused an bad recursion
else
    echo "existed: $target_dir"
    popd
fi

cd ./gmcore/netcdf/

CDIR=$(pwd)
NCDIR=$(pwd)
NFDIR=$(pwd)

cd netcdf-c-4.9.2/

# NMSLDIR=$(pwd)

CPPFLAGS="-I${H5DIR}/include -I${CURLDIR}/include -I${XML2DIR}/include " \
    LDFLAGS="-L${H5DIR}/lib -L${CURLDIR}/lib -L${XML2DIR}/lib" \
    ./configure --prefix=${CDIR}


make install -j64

export LD_LIBRARY_PATH=${NCDIR}/lib:${LD_LIBRARY_PATH}

cd  ../netcdf-fortran-4.6.1/

[ -f Makefile ] || [ -f makefile ] && make clean
CPPFLAGS="-I${NCDIR}/include -I${H5DIR}/include -I${CURLDIR}/include -I${XML2DIR}/include" \
    LDFLAGS="-L${NCDIR}/lib -L${H5DIR}/lib -L${CURLDIR}/lib -L${XML2DIR}/lib" \
    ./configure --prefix=${NFDIR}

make install -j64
