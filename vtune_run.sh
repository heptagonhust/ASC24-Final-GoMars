#!/bin/bash
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -w hepnode4
#SBATCH --exclusive
#SBATCH --output=./output/slurm-%j.out
spack load intel-oneapi-vtune
echo "******batch_run.sh*******"
cat $0
echo "******batch_run.sh*******"
source ./env.sh
message=$2
# days=$3
run ( ) {
    if [ $(hostname) != "hepnode0" ]; then
        export UCX_RC_PATH_MTU=2048
        export I_MPI_HYDRA_RMK=slurm
        export I_MPI_PIN=off
        export OMP_NUM_THREADS=1
    fi
    case_name=$1
    node=$2
    proc=$3
    # days=$4
    adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
    swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
    normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)
    prefix="./run/GMCORE-TESTBED/"
    suffix="/namelist"
    namelist_relative_path="${prefix}${case_name}/${suffix}"
    namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
    data_path="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")"
    # now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
    days=$(grep 'run_days' ${namelist_absolute_path} | sed 's/.*= *\([0-9]*\).*/\1/')
    check_file="/data/gomars_output/public/N${2}n${3}/${case_name}_${days}days/baseline.nc" 
    check_dir="/data/gomars_output/public/N${2}n${3}/" 
    if [ -f "$check_file" ]; then
        echo "exist!"
    else
        echo "not exist, use default path at N1n16!"
        echo "!!! You should notice days!"
        check_file="/data/gomars_output/public/N1n16/${case_name}_${days}days/baseline.nc"
    fi
    cd ..
    current_dir=$(pwd)
    cd gmcore/
    if [ ! -d ${data_path} ]; then    
        mkdir -p ${data_path}
    fi
    pushd ${data_path}
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
    mpirun -n $3 -ppn $( expr $3 / $2 ) \
    vtune -collect hotspot -knob enable-stack-collection=true \
    -knob sampling-mode=hw -knob stack-size=0 -knob sampling-interval=100 \
    -r /home/xyw/ascgomars/ASC24-Final-GoMars/vtune/${case_name}_r001.hs -finalization-mode=deferred -- \
    ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
    rm -rf opt.nc
    mv *.nc opt.nc
    now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
    popd
    fd1="fd1=\"${check_file}\""
    fd2="fd2=\"${now_dir}\""
    echo $fd1 
    echo $fd2
    source activate ncl_stable
    if [[ $namelist_absolute_path == *"adv"* ]]; then
        echo "adv case"
        exe_absolute_path=$adv_exe_absolute_path
        ncl ../script/ncl/adv_verify_answer.ncl $fd1 $fd2
    elif [[ $namelist_absolute_path == *"swm"* ]]; then
        echo "swm case"
        exe_absolute_path=$swm_exe_absolute_path
        ncl ../script/ncl/swm_verify_answer.ncl $fd1 $fd2
    else
        echo "normal case"
        exe_absolute_path=$normal_exe_absolute_path
        ncl ../script/ncl/normal_verify_answer.ncl $fd1 $fd2
    fi
}
if [ -z $1 ]; then
    echo "Usage: $0 <case_list> <message>"
    exit 1
fi
cmd_array=()
while IFS= read -r cmd
do
    cmd_array+=("$cmd")
done < $1
for cmd in "${cmd_array[@]}"
do
    pushd ./gmcore
    run $cmd $2
    popd
done