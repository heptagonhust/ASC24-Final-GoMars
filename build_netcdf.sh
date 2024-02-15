#!/bin/bash

spack load cmake@3.24.4
spack load intel-oneapi-compilers@2024.0.1
spack load intel-oneapi-mkl@2024.0.0
spack load intel-oneapi-mpi@2021.11.0
spack load libxml2
spack load curl
spack load hdf5/fxhrrhv

export CC=mpiicx
export FC=mpiifort
export F77=mpiifort

# export H5DIR=$(spack location -i hdf5)
export H5DIR=$(spack location -i hdf5/fxhrrhv)
export CURLDIR=$(spack location -i curl)
export XML2DIR=$(spack location -i libxml2)

cd ./gmcore/netcdf/

CDIR=$(pwd)
NCDIR=$(pwd)
NFDIR=$(pwd)

cd netcdf-c-4.9.2/

# export H5DIR=/opt/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/hdf5-1.14.3-fxhrrhv46hbchtxp255okz3r2dotzmng/
# export CURLDIR=/opt/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/curl-8.4.0-aayk5xsgzo4cyrucspyiktovue3mwnmn/

CPPFLAGS="-I${H5DIR}/include -I${CURLDIR}/include -I${XML2DIR}/include" \
    LDFLAGS="-L${H5DIR}/lib -L${CURLDIR}/lib -L${XML2DIR}/lib" \
    ./configure --prefix=${CDIR}

# CPPFLAGS="-I${H5DIR}/include" \
#     LDFLAGS="-L${H5DIR}/lib" \
#     ./configure --enable-parallel --disable-byterange
# make check
make install

# cd ./gmcore/netcdf/

# CDIR=$(pwd)


export LD_LIBRARY_PATH=${NCDIR}/lib:${LD_LIBRARY_PATH}

# echo $NFDIR

cd  ../netcdf-fortran-4.6.1/

# make distclean
make clean
CPPFLAGS="-I${NCDIR}/include -I${H5DIR}/include -I${CURLDIR}/include -I${XML2DIR}/include" \
    LDFLAGS="-L${NCDIR}/lib -L${H5DIR}/lib -L${CURLDIR}/lib -L${XML2DIR}/lib" \
    ./configure --prefix=${NFDIR}

# CPPFLAGS="-I${NCDIR}/include" \
#     LDFLAGS="-L${NCDIR}/lib" \ 
#     ./configure --prefix=${NFDIR}

# CPPFLAGS="-I${NCDIR}/include"
# LDFLAGS="-L${NCDIR}/lib" 
# ./configure --prefix=${NFDIR}

# CPPFLAGS="-I${NCDIR}/include" \
#     LDFLAGS="-L${NCDIR}/lib" \ 
#     ./configure --prefix=${NFDIR}

# make check
make install