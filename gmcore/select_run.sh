#!/bin/bash
#SBATCH -p GPU_nodes
#SBATCH -N 1
#SBATCH -n 15
#SBATCH --output=./output/gomars_test.%j.out



display_help() {
    echo "Options:"
    echo "-h: for help"
    echo "-c <case_name> : run selected case"
    echo "-l show cases"
}

example_list=('adv_dc4.360x180' 'adv_dcmip12.360x180' 'adv_mv.360x180' 'adv_sr.360x180' 'bw.150x75' 'bw.180x90' 'bw.360x180' 'dcmip31' 'hs.128x72' 'ksp15_01' 'mz.180x90' 'mz.360x180' 'mz.1440x720' 
'pgf.360x181' 'rh.180x90' 'rh.360x180' 'ss.180x90' 'ss.360x180' 'swm_cp.360x180' 'swm_jz.180x90' 'swm_jz.360x180' 'swm_jz.720x360' 'swm_mz.180x90' 'swm_mz.360x180' 'swm_mz.720x360' 'swm_rh.180x90'
'swm_rh.360x180' 'swm_rh.720x360' 'swm_rh.3600x1800' 'swm_sg.360x180' 'swm_sp.360x180' 'swm_vr.180x90' 'swm_vr.360x180' 'swm_vr.512x256' 'tc.360x180')
current_user=$(whoami)


print_example_list() {
    local i
    for ((i = 0; i < ${#example_list[@]}; i+=3)); do
        echo "${example_list[i]} ${example_list[i+1]} ${example_list[i+2]}"
    done
}

run_case() {
    case_name="$1"
    echo "run case"
    spack load intel-oneapi-mpi@2021.11.0+envmods~external-libfabric+generic-names~ilp64
    spack load intel-oneapi-compilers@2024.0.1/xbteted
    spack load hdf5/fxhrrhv

    # ulimit -s unlimit
    # edexample_list=('adv_dc4.360x180' 'adv_dcmip12.360x180' 'adv_mv.360x180' 'adv_sr.360x180' 'bw.150x75' 'bw.180x90' 'bw.360x180' 'dcmip31' 'hs.128x72' 'ksp15_01' 'mz.180x90' 'mz.360x180' 'mz.1440x720' 
    # 'pgf.360x181' 'rh.180x90' 'rh.360x180' 'ss.180x90' 'ss.360x180' 'swm_cp.360x180' 'swm_jz.180x90' 'swm_jz.360x180' 'swm_jz.720x360' 'swm_mz.180x90' 'swm_mz.360x180' 'swm_mz.720x360' 'swm_rh.180x90'
    # 'swm_rh.360x180' 'swm_rh.720x360' 'swm_rh.3600x1800' 'swm_sg.360x180' 'swm_sp.360x180' 'swm_vr.180x90' 'swm_vr.360x180' 'swm_vr.512x256' 'tc.360x180')

    

    adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
    swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
    normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)

    prefix="./run/GMCORE-TESTBED/"
    suffix="/namelist"
    namelist_relative_path="${prefix}${case_name}/${suffix}"
    namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
    data_path="/data/gomars_output/$current_user/ncdata/${case_name}"

    if [ ! -d ${data_path} ]; then    
        mkdir -p ${data_path}
    fi

    #set datapath to yours

    cd ${data_path}
    if [[ $namelist_absolute_path == *"adv"* ]]; then
        echo "adv case"
        exe_absolute_path=$adv_exe_absolute_path
    elif [[ $namelist_absolute_path == *"swm"* ]]; then
        echo "swm case"
        exe_absolute_path=$swm_exe_absolute_path
    else
        echo "normal case"
        exe_absolute_path=$normal_exe_absolute_path
    fi

    mpirun $exe_absolute_path $namelist_absolute_path
}



base_dir="$(pwd)"
node_name=$(hostname)

if [ $node_name != "hustcpu02" ]; then
    export UCX_RC_PATH_MTU=2048
    export I_MPI_HYDRA_RMK=slurm
    export I_MPI_PIN=off
    export OMP_NUM_THREADS=1
fi
# export I_MPI_STACKSIZE=8388608
# export KMP_STACKSIZE=83886080
input="$1"

echo "====="
echo "cat $0"
cat $0
echo "====="

if [ "$input" == "-h" ]; then
    display_help
elif [ "$input" == "-l" ]; then
    print_example_list
elif [ "$input" == "-c" ]; then
    if [ -z "$2" ]; then
        echo "Usage: $0 -c case_name"
        exit 1
    fi
    run_case $2 
else
    display_help
fi
