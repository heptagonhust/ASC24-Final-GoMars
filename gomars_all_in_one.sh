#!/bin/bash

current_dir=$(pwd)
current_user=$(whoami)
echo "The current user is: $current_user"

if [ ! -d "/data/gomars_output/$current_user" ]; then
    mkdir "/data/gomars_output/$current_user"
    echo "The current user is: me"
fi

# if [ ! -d "$current_dir/gmcore" ]; then
#     git clone https://gitee.com/dongli85/GMCORE gmcore
# fi

# pull libs and data
pushd gmcore
if [ ! -d "$current_dir/gmcore/data" ]; then
    mkdir "$current_dir/gmcore/data"
    # ./pull_data.py -p earth
    # ./pull_data.py -p mars
    cp -r /data/gomars_data/data/* ./data/
fi

if [ ! -d "$current_dir/gmcore/output" ]; then
    # ./pull_libs.py
    mkdir output/
fi

popd
if [ ! -d "$current_dir/gmcore/netcdf" ]; then
    ./build_netcdf.sh
fi

if [ ! -d "$current_dir/gmcore/gptl" ]; then
    ./build_gptl.sh
fi

if [ ! -d "$current_dir/gmcore/run" ]; then
    pushd gmcore
    mkdir run
    pushd run
    # git clone https://gitee.com/dongli85/GMCORE-TESTBED
    cp -r /data/gomars_data/namelist/GMCORE-TESTBED .
    popd
    popd
fi
# download netcdf
# current_dir=$(pwd)
# target_dir="$current_dir/netcdf"
# if [ ! -d "$target_dir" ]; then
#     mkdir -p "$target_dir"
#     echo "created: $target_dir"
#     pushd netcdf
#     wget https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
#     wget https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.1/netcdf-fortran-4.6.1.tar.gz
#     tar -zxvf netcdf-c-4.9.2.tar.gz
#     tar -zxvf netcdf-fortran-4.6.1.tar.gz
#     popd
#     popd
#     # build netcdf
#     # build netcdf-c && netcdf-fortran
#     # cd ../../
#     ./build_netcdf.sh
# else
#     echo "existed: $target_dir"
#     popd
# fi



# build gmcore
./build_gmcore.sh

# pushd gmcore

# run a case
# sbatch select_run.sh -c swm_mz.360x180
