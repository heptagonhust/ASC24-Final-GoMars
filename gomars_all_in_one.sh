#!/bin/bash

current_dir=$(pwd)

if [ ! -d "$current_dir/gmcore" ]; then
    git clone https://gitee.com/dongli85/GMCORE gmcore
fi

# pull libs and data
pushd gmcore
if [ ! -d "$current_dir/gmcore/data" ]; then
    ./pull_libs.py
    ./pull_data.py -p earth
    ./pull_data.py -p mars
fi

# download netcdf
current_dir=$(pwd)
target_dir="$current_dir/netcdf"
if [ ! -d "$target_dir" ]; then
    mkdir -p "$target_dir"
    echo "created: $target_dir"
    pushd netcdf
    wget https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
    wget https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.1/netcdf-fortran-4.6.1.tar.gz
    tar -zxvf netcdf-c-4.9.2.tar.gz
    tar -zxvf netcdf-fortran-4.6.1.tar.gz
    popd
    popd
    # build netcdf
    # build netcdf-c && netcdf-fortran
    # cd ../../
    ./build_netcdf.sh
else
    echo "existed: $target_dir"
    popd
fi



# build gmcore
./build_gmcore.sh
