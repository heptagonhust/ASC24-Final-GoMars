#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -w hepnode[1-3]
#SBATCH --exclusive
#SBATCH --output=./output/slurm-%j.out


echo "******batch_run.sh*******"
cat $0
echo "******batch_run.sh*******"


source ./env.sh

message=$2

run ( ) {
	
	if [ $(hostname) != "hepnode0" ]; then
		export UCX_RC_PATH_MTU=2048
		export I_MPI_HYDRA_RMK=slurm
		export I_MPI_PIN=off
		export OMP_NUM_THREADS=4
	fi
	case_name=$1
	node=$2
	proc=$3

	adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
	swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
	normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)

	prefix="./run/GMCORE-TESTBED/"
	suffix="/namelist"
	namelist_relative_path="${prefix}${case_name}/${suffix}"
	namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
	data_path="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")"
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
	# bash -c "mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path" #doesn't work
	# mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path
	mpirun -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
	popd
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
