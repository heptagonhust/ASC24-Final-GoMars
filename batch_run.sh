#!/bin/bash
if [  $(hostname) != "hustcpu02" ]; then
	export UCX_RC_PATH_MTU=2048
	export I_MPI_HYDRA_RMK=slurm
	export OMP_NUM_THREADS=1
fi

source ./env.sh

run ( ) {
	
	if [ $(hostname) != "hustcpu02" ]; then
		export UCX_RC_PATH_MTU=2048
		export I_MPI_HYDRA_RMK=slurm
		export I_MPI_PIN=off
		export OMP_NUM_THREADS=1
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
	data_path="/data/gomars_output/$(whoami)/ncdata/${case_name}/N${2}n${3}"

	if [ ! -d ${data_path} ]; then    
		mkdir -p ${data_path}
	fi

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
	# echo \$1:$1 \$2:$2 \$3:$3
	# echo exe_absolute_path:$exe_absolute_path namelist_absolute_path:$namelist_absolute_path
	echo "srun -N $2 -n $3 $exe_absolute_path $namelist_absolute_path"
	echo $(pwd)
}


if [ -z $1 ]; then
	echo "Usage: $0 <case_list>"
	exit 1
fi

while IFS= read -r cmd
do
	pushd ./gmcore
	case_name=$(echo $cmd | cut -d' ' -f1)
	node=$(echo $cmd | cut -d' ' -f2)
	proc=$(echo $cmd | cut -d' ' -f3)
	run $case_name $node $proc
	popd
done < $1

