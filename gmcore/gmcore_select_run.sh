#!/bin/bash
#SBATCH -p GPU_nodes
#SBATCH -N 2
#SBATCH -n 60
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
 
print_example_list() {
    local i
    for ((i = 0; i < ${#example_list[@]}; i+=3)); do
        echo "${example_list[i]} ${example_list[i+1]} ${example_list[i+2]}"
    done
}



base_dir="$(pwd)"

export UCX_RC_PATH_MTU=2048
export I_MPI_HYDRA_RMK=slurm
export I_MPI_PIN=off
export OMP_NUM_THREADS=1
# export I_MPI_STACKSIZE=8388608
# export KMP_STACKSIZE=83886080
input="$1"

if [ "$input" == "-h" ]; then
    display_help
elif [ "$input" == "-l" ]; then
    print_example_list
elif [ "$input" == "-c" ]; then
    case_name="$2"
    echo "run case"
if [ -z "$2" ]; then
    echo "Usage: $0 case_name"
    exit 1
fi
spack load intel-oneapi-mkl@2024.0.0
spack load intel-oneapi-mpi@2021.11.0
spack load intel-oneapi-compilers@2024.0.1
spack load hdf5/fxhrrhv

# ulimit -s unlimited

echo "====="
echo "cat $0"
cat $0
echo "====="

adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)
# gmcore_adv_driver.exe

prefix="./run/GMCORE-TESTBED/"
suffix="/namelist"
namelist_relative_path="${prefix}${case_name}/${suffix}"

namelist_absolute_path=$(readlink -f ${namelist_relative_path} )

#set datapath to yours
data_path="/data/gomars_output/xyw/ncdata/${case_name}"


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
# cd /data/gomars_output/xyw/ncdata/swm_rh.360x180


mpirun $exe_absolute_path $namelist_absolute_path
# export LD_LIBRARY_PATH=/opt/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/hdf5-1.14.3-fxhrrhv46hbchtxp255okz3r2dotzmng/lib:$LD_LIBRARY_PATH
# mpirun ./build/gmcore_driver.exe ./run/GMCORE-TESTBED/mz.3600x1800/namelist
# mpirun ./build/gmcore_swm_driver.exe ./run/GMCORE-TESTBED/swm_mz.1440x720/namelist
else
    display_help
fi