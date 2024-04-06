#!/bin/bash

set -e
cd "$(dirname $0)" || exit 1

source ./env_xyw.sh

export HDF5_PATH=$(spack location -i hdf5)
export CURL_PATH=$(spack location -i curl)
export XML2_PATH=$(spack location -i libxml2)


export LD_LIBRARY_PATH="/home/asterich/ASC24-Final-GoMars/gmcore/netcdf/lib/libnetcdf.a:${LD_LIBRARY_PATH}"
export NETCDF_PATH=/home/asterich/ASC24-Final-GoMars/gmcore/netcdf
export CPPFLAGS="-I/home/asterich/ASC24-Final-GoMars/gmcore/netcdf/include"
export NetCDF_C_LIBRARY=/home/asterich/ASC24-Final-GoMars/gmcore/netcdf/lib/libnetcdf.a
export NetCDF_Fortran_LIBRARY=/home/asterich/ASC24-Final-GoMars/gmcore/netcdf/lib/libnetcdff.a
export NetCDF_INCLUDE_DIR=/home/asterich/ASC24-Final-GoMars/gmcore/netcdf/include

cd gmcore
if [ ! -d ParallelIO ]; then
  # git clone https://github.com/NCAR/ParallelIO.git
  cp -r /data/gomars_data/ParallelIO ./
fi

cd ParallelIO
pwd


### These commands are autotools, so it is necessary
### only when building the project for the first time.
# libtoolize --force
# aclocal
# autoheader
# autoconf
# automake --force-missing --add-missing


[ -f Makefile ] || [ -f makefile ] && make clean

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${HDF5_PATH}/lib/libhdf5.a:${HDF5_PATH}/lib/libhdf5_hl.a:${XML2_PATH}/lib/libxml2.so"

export LIBS="-lhdf5 -lhdf5_hl -lcurl -lxml2 -lc -lm"

CPPFLAGS="${CPPFLAGS} -I${HDF5_PATH}/include -I${CURL_PATH}/include -I${XML2_PATH}/include" \
 LDFLAGS="-L${NETCDF_PATH}/lib -L${HDF5_PATH}/lib -L${CURL_PATH}/lib -L${XML2_PATH}/lib" \
  ./configure --prefix="$(pwd)" --enable-fortran --disable-test-runs --disable-pnetcdf --enable-shared=no --enable-static=libnetcdf,netcdf
make check
make install -j64

# if [ -d build ]; then
#     rm -rf build
# fi
# mkdir build && cd build
# echo "current build path: $(pwd)"
# echo "current LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"

# export LIBXML2_PATH=$(spack location -i libxml2)

# CC=mpiicx FC=mpiifx cmake --log-level=VERBOSE \
#     -DNetCDF_PATH=/home/asterich/ASC24-Final-GoMars/gmcore/netcdf \
#     -DHDF5_PATH=$(spack location -i hdf5) \
#     -DCMAKE_SHARED_LINKER_FLAGS="-I${LIBXML2_PATH}/include -L${LIBXML2_PATH}/lib" \
#     -DWITH_PNETCDF=OFF ../

# make VERBOSE=1 >build.log 2>build.err