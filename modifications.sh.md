diff --git a/1Nbatch_run.sh b/1Nbatch_run.sh
deleted file mode 100755
index a197089..0000000
--- a/1Nbatch_run.sh
+++ /dev/null
@@ -1,152 +0,0 @@
-#!/bin/bash
-#SBATCH -N 1
-#SBATCH -n 64
-#SBATCH -w hepnode[0-4]
-#SBATCH --exclusive
-#SBATCH --output=./output/slurm-%j.out
-
-
-echo "******batch_run.sh*******"
-cat $0
-echo "******batch_run.sh*******"
-
-
-source ./env.sh
-
-message=$2
-# days=$3
-
-run ( ) {
-	
-	if [ $(hostname) != "hepnode0" ]; then
-		# export I_MPI_FABRICS=shm:ofa
-		export I_MPI_SHM=clx_avx512
-		export UCX_RC_PATH_MTU=4096
-		export I_MPI_HYDRA_RMK=slurm
-		export I_MPI_PIN=off
-		export OMP_NUM_THREADS=1
-		# export I_MPI_ASYNC_PROGRESS=1
-		export I_MPI_DEBUG=10
-		export I_MPI_VAR_CHECK_SPELLING=1
-		# export I_MPI_INTRANODE_EAGER_THRESHLOD=1024
-  		# export FI_OFI_RXM_RX_SIZE=4096
-  		# export FI_OFI_RXM_TX_SIZE=4096
-
-
-		# export I_MPI_STATS=2
-		# export I_MPI_SHM_HEAP_CSIZE=-1
-		# export I_MPI_CACHE_BYPASS=1
-		# export I_MPI_WAIT_=1
-		# export I_MPI_SHM_NUM_BUFFERS=-1
-		# export I_MPI_SHM_BUFFER_SIZE=102400
-		export I_MPI_MALLOC=1
-	fi
-	export I_MPI_PIN=off
-
-
-	# export I_MPI_SHM_HEAP=1
-	# export I_MPI_PLATFORM=auto
-	# export I_MPI_TUNING_MODE=auto
-	# export I_MPI_TUNING_AUTO_POLICY=max
-	# export I_MPI_EAGER_THRESHOLD=1024
-	case_name=$1
-	node=$2
-	proc=$3
-	# days=$4
-	
-
-	adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
-	swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
-	normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)
-
-	prefix="./run/GMCORE-TESTBED/"
-	suffix="/namelist"
-	namelist_relative_path="${prefix}${case_name}/${suffix}"
-	namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
-	data_path="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")"
-	# now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-	days=$(grep 'run_days' ${namelist_absolute_path} | sed 's/.*= *\([0-9]*\).*/\1/')
-
-	check_file="/data/gomars_output/public/N${2}n${3}/${case_name}_${days}days/baseline.nc" 
-	check_dir="/data/gomars_output/public/N${2}n${3}/" 
-
-	if [ -f "$check_file" ]; then
-		echo "exist!"
-	else
-		echo "not exist, use default path at N1n16!"
-		echo "!!! You should notice days!"
-		check_file="/data/gomars_output/public/N1n16/${case_name}_${days}days/baseline.nc"
-	fi
-	cd ..
-	current_dir=$(pwd)
-	# mpitune_fast -f ./hostfile
-	cd gmcore/
-
-	if [ ! -d ${data_path} ]; then    
-		mkdir -p ${data_path}
-	fi
-
-	pushd ${data_path}
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-	fi
-	# bash -c "mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path" #doesn't work
-	# mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path
-	# mpitune_fast -hf hosts
-	# mpirun -genv I_MPI_PIN_PROCESSOR_LIST map=scatter -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-	mpirun -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-
-	rm -rf opt.nc
-	mv *.nc opt.nc
-	now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-
-	popd
-	fd1="fd1=\"${check_file}\""
-	fd2="fd2=\"${now_dir}\""
-
-	echo $fd1 
-	echo $fd2
-	source activate ncl_stable
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-		ncl ../script/ncl/adv_verify_answer.ncl $fd1 $fd2
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-		ncl ../script/ncl/swm_verify_answer.ncl $fd1 $fd2
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-		ncl ../script/ncl/normal_verify_answer.ncl $fd1 $fd2
-	fi
-
-	
-}
-
-
-if [ -z $1 ]; then
-	echo "Usage: $0 <case_list> <message>"
-	exit 1
-fi
-
-cmd_array=()
-
-while IFS= read -r cmd
-do
-	cmd_array+=("$cmd")
-done < $1
-
-for cmd in "${cmd_array[@]}"
-do
-	pushd ./gmcore
-	run $cmd $2
-	popd
-done
diff --git a/1Nfinal_run.sh b/1Nfinal_run.sh
deleted file mode 100755
index 561813d..0000000
--- a/1Nfinal_run.sh
+++ /dev/null
@@ -1,105 +0,0 @@
-#!/bin/bash
-#SBATCH -N 1
-#SBATCH -n 64
-#SBATCH -w hepnode0
-#SBATCH --exclusive
-#SBATCH --output=./output/slurm-%j.out
-
-
-echo "******final_run.sh*******"
-cat $0
-echo "******final_run.sh*******"
-
-
-source ./env.sh
-
-message=$2
-# days=$3
-
-run ( ) {
-	
-	if [ $(hostname) != "hepnode0" ]; then
-		export UCX_RC_PATH_MTU=2048
-		export I_MPI_HYDRA_RMK=slurm
-		export I_MPI_PIN=off
-		export OMP_NUM_THREADS=1
-	fi
-	export I_MPI_PIN=off
-	case_name=$1
-	node=$2
-	proc=$3
-	# days=$4
-	
-
-	adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
-	swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
-	normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)
-
-	prefix="./run/GMCORE-TESTBED/"
-	suffix="/namelist"
-	namelist_relative_path="${prefix}${case_name}/${suffix}"
-	namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
-	echo namelist_absolute_path: ${namelist_absolute_path}
-	data_path="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")"
-	# now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-	days=$(grep 'run_days' ${namelist_absolute_path} | sed 's/.*= *\([0-9]*\).*/\1/')
-
-	check_file="/data/gomars_output/public/N${2}n${3}/${case_name}_${days}days/baseline.nc" 
-	check_dir="/data/gomars_output/public/N${2}n${3}/" 
-
-	if [ -f "$check_file" ]; then
-		echo "exist!"
-	else
-		echo "not exist, use default path at N1n16!"
-		echo "!!! You should notice days!"
-		check_file="/data/gomars_output/public/N1n16/${case_name}_${days}days/baseline.nc"
-	fi
-	cd ..
-	current_dir=$(pwd)
-	cd gmcore/
-
-	if [ ! -d ${data_path} ]; then    
-		mkdir -p ${data_path}
-	fi
-
-	pushd ${data_path}
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-	fi
-	# bash -c "mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path" #doesn't work
-	mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path
-	# mpirun -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-
-
-	now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-
-	popd
-	
-}
-
-
-if [ -z $1 ]; then
-	echo "Usage: $0 <case_list> <message>"
-	exit 1
-fi
-
-cmd_array=()
-
-while IFS= read -r cmd
-do
-	cmd_array+=("$cmd")
-done < $1
-
-for cmd in "${cmd_array[@]}"
-do
-	pushd ./gmcore
-	run $cmd $2
-	popd
-done
diff --git a/2Nbatch_run.sh b/2Nbatch_run.sh
deleted file mode 100755
index aa27978..0000000
--- a/2Nbatch_run.sh
+++ /dev/null
@@ -1,125 +0,0 @@
-#!/bin/bash
-#SBATCH -N 2
-#SBATCH -n 120
-#SBATCH -w hepnode[0-4]
-#SBATCH --exclusive
-#SBATCH --output=./output/slurm-%j.out
-
-
-echo "******batch_run.sh*******"
-cat $0
-echo "******batch_run.sh*******"
-
-
-source ./env.sh
-
-message=$2
-# days=$3
-
-run ( ) {
-	
-	if [ $(hostname) != "hepnode0" ]; then
-		export UCX_RC_PATH_MTU=2048
-		export I_MPI_HYDRA_RMK=slurm
-		export I_MPI_PIN=off
-		export OMP_NUM_THREADS=1
-	fi
-	export I_MPI_PIN=off
-	case_name=$1
-	node=$2
-	proc=$3
-	# days=$4
-	
-
-	adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
-	swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
-	normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)
-
-	prefix="./run/GMCORE-TESTBED/"
-	suffix="/namelist"
-	namelist_relative_path="${prefix}${case_name}/${suffix}"
-	namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
-	data_path="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")"
-	# now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-	days=$(grep 'run_days' ${namelist_absolute_path} | sed 's/.*= *\([0-9]*\).*/\1/')
-
-	check_file="/data/gomars_output/public/N${2}n${3}/${case_name}_${days}days/baseline.nc" 
-	check_dir="/data/gomars_output/public/N${2}n${3}/" 
-
-	if [ -f "$check_file" ]; then
-		echo "exist!"
-	else
-		echo "not exist, use default path at N1n16!"
-		echo "!!! You should notice days!"
-		check_file="/data/gomars_output/public/N1n16/${case_name}_${days}days/baseline.nc"
-	fi
-	cd ..
-	current_dir=$(pwd)
-	cd gmcore/
-
-	if [ ! -d ${data_path} ]; then    
-		mkdir -p ${data_path}
-	fi
-
-	pushd ${data_path}
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-	fi
-	# bash -c "mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path" #doesn't work
-	# mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path
-	mpirun -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-
-	rm -rf opt.nc
-	mv *.nc opt.nc
-	now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-
-	popd
-	fd1="fd1=\"${check_file}\""
-	fd2="fd2=\"${now_dir}\""
-
-	echo $fd1 
-	echo $fd2
-	source activate ncl_stable
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-		ncl ../script/ncl/adv_verify_answer.ncl $fd1 $fd2
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-		ncl ../script/ncl/swm_verify_answer.ncl $fd1 $fd2
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-		ncl ../script/ncl/normal_verify_answer.ncl $fd1 $fd2
-	fi
-
-	
-}
-
-
-if [ -z $1 ]; then
-	echo "Usage: $0 <case_list> <message>"
-	exit 1
-fi
-
-cmd_array=()
-
-while IFS= read -r cmd
-do
-	cmd_array+=("$cmd")
-done < $1
-
-for cmd in "${cmd_array[@]}"
-do
-	pushd ./gmcore
-	run $cmd $2
-	popd
-done
diff --git a/2Nfinal_run.sh b/2Nfinal_run.sh
deleted file mode 100755
index 5d80061..0000000
--- a/2Nfinal_run.sh
+++ /dev/null
@@ -1,125 +0,0 @@
-#!/bin/bash
-#SBATCH -N 2
-#SBATCH -n 120
-#SBATCH -w hepnode[0,3]
-#SBATCH --exclusive
-#SBATCH --output=./output/slurm-%j.out
-
-
-echo "******2Nfinal_run.sh*******"
-cat $0
-echo "******2Nfinal_run.sh*******"
-
-
-source ./env.sh
-
-message=$2
-# days=$3
-
-run ( ) {
-	
-	if [ $(hostname) != "hepnode0" ]; then
-		export UCX_RC_PATH_MTU=2048
-		export I_MPI_HYDRA_RMK=slurm
-		export I_MPI_PIN=off
-		export OMP_NUM_THREADS=1
-	fi
-	export I_MPI_PIN=off
-	case_name=$1
-	node=$2
-	proc=$3
-	# days=$4
-	
-
-	adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
-	swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
-	normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)
-
-	prefix="./run/GMCORE-TESTBED/"
-	suffix="/namelist"
-	namelist_relative_path="${prefix}${case_name}/${suffix}"
-	namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
-	data_path="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")"
-	# now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-	days=$(grep 'run_days' ${namelist_absolute_path} | sed 's/.*= *\([0-9]*\).*/\1/')
-
-	check_file="/data/gomars_output/public/N${2}n${3}/${case_name}_${days}days/baseline.nc" 
-	check_dir="/data/gomars_output/public/N${2}n${3}/" 
-
-	if [ -f "$check_file" ]; then
-		echo "exist!"
-	else
-		echo "not exist, use default path at N1n16!"
-		echo "!!! You should notice days!"
-		check_file="/data/gomars_output/public/N1n16/${case_name}_${days}days/baseline.nc"
-	fi
-	cd ..
-	current_dir=$(pwd)
-	cd gmcore/
-
-	if [ ! -d ${data_path} ]; then    
-		mkdir -p ${data_path}
-	fi
-
-	pushd ${data_path}
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-	fi
-	# bash -c "mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path" #doesn't work
-	# mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path
-	mpirun -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-
-	# rm -rf opt.nc
-	# mv *.nc opt.nc
-	now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-
-	popd
-	# fd1="fd1=\"${check_file}\""
-	# fd2="fd2=\"${now_dir}\""
-
-	# echo $fd1 
-	# echo $fd2
-	# source activate ncl_stable
-	# if [[ $namelist_absolute_path == *"adv"* ]]; then
-	# 	echo "adv case"
-	# 	exe_absolute_path=$adv_exe_absolute_path
-	# 	ncl ../script/ncl/adv_verify_answer.ncl $fd1 $fd2
-	# elif [[ $namelist_absolute_path == *"swm"* ]]; then
-	# 	echo "swm case"
-	# 	exe_absolute_path=$swm_exe_absolute_path
-	# 	ncl ../script/ncl/swm_verify_answer.ncl $fd1 $fd2
-	# else
-	# 	echo "normal case"
-	# 	exe_absolute_path=$normal_exe_absolute_path
-	# 	ncl ../script/ncl/normal_verify_answer.ncl $fd1 $fd2
-	# fi
-
-	
-}
-
-
-if [ -z $1 ]; then
-	echo "Usage: $0 <case_list> <message>"
-	exit 1
-fi
-
-cmd_array=()
-
-while IFS= read -r cmd
-do
-	cmd_array+=("$cmd")
-done < $1
-
-for cmd in "${cmd_array[@]}"
-do
-	pushd ./gmcore
-	run $cmd $2
-	popd
-done
diff --git a/3Nfinal_run.sh b/3Nfinal_run.sh
deleted file mode 100755
index 8ae86c2..0000000
--- a/3Nfinal_run.sh
+++ /dev/null
@@ -1,125 +0,0 @@
-#!/bin/bash
-#SBATCH -N 4
-#SBATCH -n 240
-#SBATCH -w hepnode[1-4]
-#SBATCH --exclusive
-#SBATCH --output=./output/slurm-%j.out
-
-
-echo "******3Nfinal_run.sh*******"
-cat $0
-echo "******3Nfinal_run.sh*******"
-
-
-source ./env.sh
-
-message=$2
-# days=$3
-
-run ( ) {
-	
-	if [ $(hostname) != "hepnode0" ]; then
-		export UCX_RC_PATH_MTU=2048
-		export I_MPI_HYDRA_RMK=slurm
-		export I_MPI_PIN=off
-		export OMP_NUM_THREADS=1
-	fi
-	export I_MPI_PIN=off
-	case_name=$1
-	node=$2
-	proc=$3
-	# days=$4
-	
-
-	adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
-	swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
-	normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)
-
-	prefix="./run/GMCORE-TESTBED/"
-	suffix="/namelist"
-	namelist_relative_path="${prefix}${case_name}/${suffix}"
-	namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
-	data_path="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")"
-	# now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-	days=$(grep 'run_days' ${namelist_absolute_path} | sed 's/.*= *\([0-9]*\).*/\1/')
-
-	check_file="/data/gomars_output/public/N${2}n${3}/${case_name}_${days}days/baseline.nc" 
-	check_dir="/data/gomars_output/public/N${2}n${3}/" 
-
-	if [ -f "$check_file" ]; then
-		echo "exist!"
-	else
-		echo "not exist, use default path at N1n16!"
-		echo "!!! You should notice days!"
-		check_file="/data/gomars_output/public/N1n16/${case_name}_${days}days/baseline.nc"
-	fi
-	cd ..
-	current_dir=$(pwd)
-	cd gmcore/
-
-	if [ ! -d ${data_path} ]; then    
-		mkdir -p ${data_path}
-	fi
-
-	pushd ${data_path}
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-	fi
-	# bash -c "mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path" #doesn't work
-	# mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path
-	mpirun -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-
-	# rm -rf opt.nc
-	# mv *.nc opt.nc
-	now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-
-	popd
-	# fd1="fd1=\"${check_file}\""
-	# fd2="fd2=\"${now_dir}\""
-
-	# echo $fd1 
-	# echo $fd2
-	# source activate ncl_stable
-	# if [[ $namelist_absolute_path == *"adv"* ]]; then
-	# 	echo "adv case"
-	# 	exe_absolute_path=$adv_exe_absolute_path
-	# 	ncl ../script/ncl/adv_verify_answer.ncl $fd1 $fd2
-	# elif [[ $namelist_absolute_path == *"swm"* ]]; then
-	# 	echo "swm case"
-	# 	exe_absolute_path=$swm_exe_absolute_path
-	# 	ncl ../script/ncl/swm_verify_answer.ncl $fd1 $fd2
-	# else
-	# 	echo "normal case"
-	# 	exe_absolute_path=$normal_exe_absolute_path
-	# 	ncl ../script/ncl/normal_verify_answer.ncl $fd1 $fd2
-	# fi
-
-	
-}
-
-
-if [ -z $1 ]; then
-	echo "Usage: $0 <case_list> <message>"
-	exit 1
-fi
-
-cmd_array=()
-
-while IFS= read -r cmd
-do
-	cmd_array+=("$cmd")
-done < $1
-
-for cmd in "${cmd_array[@]}"
-do
-	pushd ./gmcore
-	run $cmd $2
-	popd
-done
diff --git a/4Nfinal_run.sh b/4Nfinal_run.sh
deleted file mode 100755
index a6f7014..0000000
--- a/4Nfinal_run.sh
+++ /dev/null
@@ -1,125 +0,0 @@
-#!/bin/bash
-#SBATCH -N 4
-#SBATCH -n 240
-#SBATCH -w hepnode[1-4]
-#SBATCH --exclusive
-#SBATCH --output=./output/slurm-%j.out
-
-
-echo "******4Nfinal_run.sh*******"
-cat $0
-echo "******4Nfinal_run.sh*******"
-
-
-source ./env.sh
-
-message=$2
-# days=$3
-
-run ( ) {
-	
-	if [ $(hostname) != "hepnode0" ]; then
-		export UCX_RC_PATH_MTU=2048
-		export I_MPI_HYDRA_RMK=slurm
-		export I_MPI_PIN=off
-		export OMP_NUM_THREADS=1
-	fi
-	export I_MPI_PIN=off
-	case_name=$1
-	node=$2
-	proc=$3
-	# days=$4
-	
-
-	adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
-	swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
-	normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)
-
-	prefix="./run/GMCORE-TESTBED/"
-	suffix="/namelist"
-	namelist_relative_path="${prefix}${case_name}/${suffix}"
-	namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
-	data_path="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")"
-	# now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-	days=$(grep 'run_days' ${namelist_absolute_path} | sed 's/.*= *\([0-9]*\).*/\1/')
-
-	check_file="/data/gomars_output/public/N${2}n${3}/${case_name}_${days}days/baseline.nc" 
-	check_dir="/data/gomars_output/public/N${2}n${3}/" 
-
-	if [ -f "$check_file" ]; then
-		echo "exist!"
-	else
-		echo "not exist, use default path at N1n16!"
-		echo "!!! You should notice days!"
-		check_file="/data/gomars_output/public/N1n16/${case_name}_${days}days/baseline.nc"
-	fi
-	cd ..
-	current_dir=$(pwd)
-	cd gmcore/
-
-	if [ ! -d ${data_path} ]; then    
-		mkdir -p ${data_path}
-	fi
-
-	pushd ${data_path}
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-	fi
-	# bash -c "mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path" #doesn't work
-	# mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path
-	mpirun -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-
-	# rm -rf opt.nc
-	# mv *.nc opt.nc
-	now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-
-	popd
-	# fd1="fd1=\"${check_file}\""
-	# fd2="fd2=\"${now_dir}\""
-
-	# echo $fd1 
-	# echo $fd2
-	# source activate ncl_stable
-	# if [[ $namelist_absolute_path == *"adv"* ]]; then
-	# 	echo "adv case"
-	# 	exe_absolute_path=$adv_exe_absolute_path
-	# 	ncl ../script/ncl/adv_verify_answer.ncl $fd1 $fd2
-	# elif [[ $namelist_absolute_path == *"swm"* ]]; then
-	# 	echo "swm case"
-	# 	exe_absolute_path=$swm_exe_absolute_path
-	# 	ncl ../script/ncl/swm_verify_answer.ncl $fd1 $fd2
-	# else
-	# 	echo "normal case"
-	# 	exe_absolute_path=$normal_exe_absolute_path
-	# 	ncl ../script/ncl/normal_verify_answer.ncl $fd1 $fd2
-	# fi
-
-	
-}
-
-
-if [ -z $1 ]; then
-	echo "Usage: $0 <case_list> <message>"
-	exit 1
-fi
-
-cmd_array=()
-
-while IFS= read -r cmd
-do
-	cmd_array+=("$cmd")
-done < $1
-
-for cmd in "${cmd_array[@]}"
-do
-	pushd ./gmcore
-	run $cmd $2
-	popd
-done
diff --git a/5Nbatch_run.sh b/5Nbatch_run.sh
deleted file mode 100755
index ec1f011..0000000
--- a/5Nbatch_run.sh
+++ /dev/null
@@ -1,125 +0,0 @@
-#!/bin/bash
-#SBATCH -N 5
-#SBATCH -n 320
-#SBATCH -w hepnode[0-4]
-#SBATCH --exclusive
-#SBATCH --output=./output/slurm-%j.out
-
-
-echo "******batch_run.sh*******"
-cat $0
-echo "******batch_run.sh*******"
-
-
-source ./env.sh
-
-message=$2
-# days=$3
-
-run ( ) {
-	
-	if [ $(hostname) != "hepnode0" ]; then
-		export UCX_RC_PATH_MTU=2048
-		export I_MPI_HYDRA_RMK=slurm
-		export I_MPI_PIN=off
-		export OMP_NUM_THREADS=1
-	fi
-	export I_MPI_PIN=off
-	case_name=$1
-	node=$2
-	proc=$3
-	# days=$4
-	
-
-	adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
-	swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
-	normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)
-
-	prefix="./run/GMCORE-TESTBED/"
-	suffix="/namelist"
-	namelist_relative_path="${prefix}${case_name}/${suffix}"
-	namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
-	data_path="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")"
-	# now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-	days=$(grep 'run_days' ${namelist_absolute_path} | sed 's/.*= *\([0-9]*\).*/\1/')
-
-	check_file="/data/gomars_output/public/N${2}n${3}/${case_name}_${days}days/baseline.nc" 
-	check_dir="/data/gomars_output/public/N${2}n${3}/" 
-
-	if [ -f "$check_file" ]; then
-		echo "exist!"
-	else
-		echo "not exist, use default path at N1n16!"
-		echo "!!! You should notice days!"
-		check_file="/data/gomars_output/public/N1n16/${case_name}_${days}days/baseline.nc"
-	fi
-	cd ..
-	current_dir=$(pwd)
-	cd gmcore/
-
-	if [ ! -d ${data_path} ]; then    
-		mkdir -p ${data_path}
-	fi
-
-	pushd ${data_path}
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-	fi
-	# bash -c "mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path" #doesn't work
-	# mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path
-	mpirun -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-
-	rm -rf opt.nc
-	mv *.nc opt.nc
-	now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-
-	popd
-	fd1="fd1=\"${check_file}\""
-	fd2="fd2=\"${now_dir}\""
-
-	echo $fd1 
-	echo $fd2
-	source activate ncl_stable
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-		ncl ../script/ncl/adv_verify_answer.ncl $fd1 $fd2
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-		ncl ../script/ncl/swm_verify_answer.ncl $fd1 $fd2
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-		ncl ../script/ncl/normal_verify_answer.ncl $fd1 $fd2
-	fi
-
-	
-}
-
-
-if [ -z $1 ]; then
-	echo "Usage: $0 <case_list> <message>"
-	exit 1
-fi
-
-cmd_array=()
-
-while IFS= read -r cmd
-do
-	cmd_array+=("$cmd")
-done < $1
-
-for cmd in "${cmd_array[@]}"
-do
-	pushd ./gmcore
-	run $cmd $2
-	popd
-done
diff --git a/5Nfinal_run.sh b/5Nfinal_run.sh
deleted file mode 100755
index d586359..0000000
--- a/5Nfinal_run.sh
+++ /dev/null
@@ -1,125 +0,0 @@
-#!/bin/bash
-#SBATCH -N 5
-#SBATCH -n 320
-#SBATCH -w hepnode[0-4]
-#SBATCH --exclusive
-#SBATCH --output=./output/slurm-%j.out
-
-
-echo "******5Nfinal_run.sh*******"
-cat $0
-echo "******5Nfinal_run.sh*******"
-
-
-source ./env.sh
-
-message=$2
-# days=$3
-
-run ( ) {
-	
-	if [ $(hostname) != "hepnode0" ]; then
-		export UCX_RC_PATH_MTU=2048
-		export I_MPI_HYDRA_RMK=slurm
-		export I_MPI_PIN=off
-		export OMP_NUM_THREADS=1
-	fi
-	export I_MPI_PIN=off
-	case_name=$1
-	node=$2
-	proc=$3
-	# days=$4
-	
-
-	adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
-	swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
-	normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)
-
-	prefix="./run/GMCORE-TESTBED/"
-	suffix="/namelist"
-	namelist_relative_path="${prefix}${case_name}/${suffix}"
-	namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
-	data_path="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")"
-	# now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-	days=$(grep 'run_days' ${namelist_absolute_path} | sed 's/.*= *\([0-9]*\).*/\1/')
-
-	check_file="/data/gomars_output/public/N${2}n${3}/${case_name}_${days}days/baseline.nc" 
-	check_dir="/data/gomars_output/public/N${2}n${3}/" 
-
-	if [ -f "$check_file" ]; then
-		echo "exist!"
-	else
-		echo "not exist, use default path at N1n16!"
-		echo "!!! You should notice days!"
-		check_file="/data/gomars_output/public/N1n16/${case_name}_${days}days/baseline.nc"
-	fi
-	cd ..
-	current_dir=$(pwd)
-	cd gmcore/
-
-	if [ ! -d ${data_path} ]; then    
-		mkdir -p ${data_path}
-	fi
-
-	pushd ${data_path}
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-	fi
-	# bash -c "mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path" #doesn't work
-	# mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path
-	mpirun -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-
-	# rm -rf opt.nc
-	# mv *.nc opt.nc
-	now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-
-	popd
-	# fd1="fd1=\"${check_file}\""
-	# fd2="fd2=\"${now_dir}\""
-
-	# echo $fd1 
-	# echo $fd2
-	# source activate ncl_stable
-	# if [[ $namelist_absolute_path == *"adv"* ]]; then
-	# 	echo "adv case"
-	# 	exe_absolute_path=$adv_exe_absolute_path
-	# 	ncl ../script/ncl/adv_verify_answer.ncl $fd1 $fd2
-	# elif [[ $namelist_absolute_path == *"swm"* ]]; then
-	# 	echo "swm case"
-	# 	exe_absolute_path=$swm_exe_absolute_path
-	# 	ncl ../script/ncl/swm_verify_answer.ncl $fd1 $fd2
-	# else
-	# 	echo "normal case"
-	# 	exe_absolute_path=$normal_exe_absolute_path
-	# 	ncl ../script/ncl/normal_verify_answer.ncl $fd1 $fd2
-	# fi
-
-	
-}
-
-
-if [ -z $1 ]; then
-	echo "Usage: $0 <case_list> <message>"
-	exit 1
-fi
-
-cmd_array=()
-
-while IFS= read -r cmd
-do
-	cmd_array+=("$cmd")
-done < $1
-
-for cmd in "${cmd_array[@]}"
-do
-	pushd ./gmcore
-	run $cmd $2
-	popd
-done
diff --git a/batch_nh.sh b/batch_nh.sh
deleted file mode 100755
index c218df3..0000000
--- a/batch_nh.sh
+++ /dev/null
@@ -1,152 +0,0 @@
-#!/bin/bash
-#SBATCH -N 4
-#SBATCH -n 200
-#SBATCH -w hepnode[0-4]
-#SBATCH --exclusive
-#SBATCH --output=./output/slurm-%j.out
-
-
-echo "******batch_run.sh*******"
-cat $0
-echo "******batch_run.sh*******"
-
-
-source ./env.sh
-
-message=$2
-# days=$3
-
-run ( ) {
-	
-	if [ $(hostname) != "hepnode0" ]; then
-		# export I_MPI_FABRICS=shm:ofa
-		export I_MPI_SHM=clx_avx512
-		export UCX_RC_PATH_MTU=4096
-		export I_MPI_HYDRA_RMK=slurm
-		export I_MPI_PIN=off
-		export OMP_NUM_THREADS=1
-		# export I_MPI_ASYNC_PROGRESS=1
-		export I_MPI_DEBUG=10
-		export I_MPI_VAR_CHECK_SPELLING=1
-		# export I_MPI_INTRANODE_EAGER_THRESHLOD=1024
-  		# export FI_OFI_RXM_RX_SIZE=4096
-  		# export FI_OFI_RXM_TX_SIZE=4096
-
-
-		# export I_MPI_STATS=2
-		# export I_MPI_SHM_HEAP_CSIZE=-1
-		# export I_MPI_CACHE_BYPASS=1
-		# export I_MPI_WAIT_=1
-		# export I_MPI_SHM_NUM_BUFFERS=-1
-		# export I_MPI_SHM_BUFFER_SIZE=102400
-		export I_MPI_MALLOC=1
-	fi
-	export I_MPI_PIN=off
-
-
-	# export I_MPI_SHM_HEAP=1
-	# export I_MPI_PLATFORM=auto
-	# export I_MPI_TUNING_MODE=auto
-	# export I_MPI_TUNING_AUTO_POLICY=max
-	# export I_MPI_EAGER_THRESHOLD=1024
-	case_name=$1
-	node=$2
-	proc=$3
-	# days=$4
-	
-
-	adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
-	swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
-	normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)
-
-	prefix="./run/GMCORE-TESTBED/"
-	suffix="/namelist"
-	namelist_relative_path="${prefix}${case_name}/${suffix}"
-	namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
-	data_path="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")"
-	# now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-	days=$(grep 'run_days' ${namelist_absolute_path} | sed 's/.*= *\([0-9]*\).*/\1/')
-
-	check_file="/data/gomars_output/public/N${2}n${3}/${case_name}_${days}days/baseline.nc" 
-	check_dir="/data/gomars_output/public/N${2}n${3}/" 
-
-	if [ -f "$check_file" ]; then
-		echo "exist!"
-	else
-		echo "not exist, use default path at N1n16!"
-		echo "!!! You should notice days!"
-		check_file="/data/gomars_output/public/N1n16/${case_name}_${days}days/baseline.nc"
-	fi
-	cd ..
-	current_dir=$(pwd)
-	# mpitune_fast -f ./hostfile
-	cd gmcore/
-
-	if [ ! -d ${data_path} ]; then    
-		mkdir -p ${data_path}
-	fi
-
-	pushd ${data_path}
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-	fi
-	# bash -c "mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path" #doesn't work
-	# mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path
-	# mpitune_fast -hf hosts
-	# mpirun -genv I_MPI_PIN_PROCESSOR_LIST map=scatter -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-	mpirun -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-
-	rm -rf opt.nc
-	mv *.nc opt.nc
-	now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-
-	popd
-	fd1="fd1=\"${check_file}\""
-	fd2="fd2=\"${now_dir}\""
-
-	echo $fd1 
-	echo $fd2
-	source activate ncl_stable
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-		ncl ../script/ncl/adv_verify_answer.ncl $fd1 $fd2
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-		ncl ../script/ncl/swm_verify_answer.ncl $fd1 $fd2
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-		ncl ../script/ncl/normal_verify_answer.ncl $fd1 $fd2
-	fi
-
-	
-}
-
-
-if [ -z $1 ]; then
-	echo "Usage: $0 <case_list> <message>"
-	exit 1
-fi
-
-cmd_array=()
-
-while IFS= read -r cmd
-do
-	cmd_array+=("$cmd")
-done < $1
-
-for cmd in "${cmd_array[@]}"
-do
-	pushd ./gmcore
-	run $cmd $2
-	popd
-done
diff --git a/batch_run.sh b/batch_run.sh
deleted file mode 100755
index 3f5bfa4..0000000
--- a/batch_run.sh
+++ /dev/null
@@ -1,125 +0,0 @@
-#!/bin/bash
-#SBATCH -N 4
-#SBATCH -n 240
-#SBATCH -w hepnode[0-4]
-#SBATCH --exclusive
-#SBATCH --output=./output/slurm-%j.out
-
-
-echo "******batch_run.sh*******"
-cat $0
-echo "******batch_run.sh*******"
-
-
-source ./env.sh
-
-message=$2
-# days=$3
-
-run ( ) {
-	
-	if [ $(hostname) != "hepnode0" ]; then
-		export UCX_RC_PATH_MTU=2048
-		export I_MPI_HYDRA_RMK=slurm
-		export I_MPI_PIN=off
-		export OMP_NUM_THREADS=1
-	fi
-	export I_MPI_PIN=off
-	case_name=$1
-	node=$2
-	proc=$3
-	# days=$4
-	
-
-	adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
-	swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
-	normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)
-
-	prefix="./run/GMCORE-TESTBED/"
-	suffix="/namelist"
-	namelist_relative_path="${prefix}${case_name}/${suffix}"
-	namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
-	data_path="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")"
-	# now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-	days=$(grep 'run_days' ${namelist_absolute_path} | sed 's/.*= *\([0-9]*\).*/\1/')
-
-	check_file="/data/gomars_output/public/N${2}n${3}/${case_name}_${days}days/baseline.nc" 
-	check_dir="/data/gomars_output/public/N${2}n${3}/" 
-
-	if [ -f "$check_file" ]; then
-		echo "exist!"
-	else
-		echo "not exist, use default path at N1n16!"
-		echo "!!! You should notice days!"
-		check_file="/data/gomars_output/public/N1n16/${case_name}_${days}days/baseline.nc"
-	fi
-	cd ..
-	current_dir=$(pwd)
-	cd gmcore/
-
-	if [ ! -d ${data_path} ]; then    
-		mkdir -p ${data_path}
-	fi
-
-	pushd ${data_path}
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-	fi
-	# bash -c "mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path" #doesn't work
-	# mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path
-	mpirun -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-
-	rm -rf opt.nc
-	mv *.nc opt.nc
-	now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-
-	popd
-	fd1="fd1=\"${check_file}\""
-	fd2="fd2=\"${now_dir}\""
-
-	echo $fd1 
-	echo $fd2
-	source activate ncl_stable
-	if [[ $namelist_absolute_path == *"adv"* ]]; then
-		echo "adv case"
-		exe_absolute_path=$adv_exe_absolute_path
-		ncl ../script/ncl/adv_verify_answer.ncl $fd1 $fd2
-	elif [[ $namelist_absolute_path == *"swm"* ]]; then
-		echo "swm case"
-		exe_absolute_path=$swm_exe_absolute_path
-		ncl ../script/ncl/swm_verify_answer.ncl $fd1 $fd2
-	else
-		echo "normal case"
-		exe_absolute_path=$normal_exe_absolute_path
-		ncl ../script/ncl/normal_verify_answer.ncl $fd1 $fd2
-	fi
-
-	
-}
-
-
-if [ -z $1 ]; then
-	echo "Usage: $0 <case_list> <message>"
-	exit 1
-fi
-
-cmd_array=()
-
-while IFS= read -r cmd
-do
-	cmd_array+=("$cmd")
-done < $1
-
-for cmd in "${cmd_array[@]}"
-do
-	pushd ./gmcore
-	run $cmd $2
-	popd
-done
diff --git a/bind_cpu.sh b/bind_cpu.sh
deleted file mode 100755
index cbb53e7..0000000
--- a/bind_cpu.sh
+++ /dev/null
@@ -1,25 +0,0 @@
-#!/bin/bash
-# LOCAL_RANK=$SLURM_LOCALID
-LOCAL_RANK=$MPI_LOCALRANKID # for Intel MPI
-# LOCAL_SIZE=$SLURM_TASKS_PER_NODE 
-# LOCAL_SIZE=${LOCAL_SIZE//(x3)/}
-LOCAL_SIZE=$MPI_LOCALNRANKS # for Intel MPI
-NCPUS=64 # Number of logical cores
-NUM_NUMA=8
-
-
-CORES_PER_PROCESS=$(($NCPUS / $LOCAL_SIZE))
-HARDCORES_PER_PROCESS=$(($CORES_PER_PROCESS))
-
-NUMA_ID=$(($NUM_NUMA * $LOCAL_RANK / $LOCAL_SIZE))
-
-CORE_START1=$(( $HARDCORES_PER_PROCESS * $LOCAL_RANK ))
-CORE_END1=$(( $HARDCORES_PER_PROCESS * $LOCAL_RANK + $HARDCORES_PER_PROCESS - 1 ))
-# CORE_START2=$(( $HARDCORES_PER_PROCESS * $LOCAL_RANK + 64 ))
-# CORE_END2=$(( $HARDCORES_PER_PROCESS * $LOCAL_RANK + $HARDCORES_PER_PROCESS + 63))
-
-CORES=$(seq -s, $CORE_START1 $CORE_END1)
-# CORES=$CORES",$(seq -s, $CORE_START2 $CORE_END2)"
-
-echo "Process $LOCAL_RANK on $(hostname) bound to cores $CORES"
-exec numactl -m "$NUMA_ID" -C "$CORES" $@
\ No newline at end of file
diff --git a/build_gmcore.sh b/build_gmcore.sh
index 02a073c..893eb59 100755
--- a/build_gmcore.sh
+++ b/build_gmcore.sh
@@ -1,51 +1,56 @@
 #!/bin/bash
 
-set -e
-cd "$(dirname $0)" || exit 1
+cd gmcore/
 
-source ./env.sh
+./pull_libs.py
 
-cd gmcore
+spack load cmake@3.24.4
+spack load intel-oneapi-compilers@2024.0.1
+spack load intel-oneapi-mkl@2024.0.0
 
-mkdir -p ./lib
-if [ x"$(ls -A lib)" = x"" ]; then
-  cp -r /data/gomars_libs/gmcore_libs/* ./lib
-fi
+spack load intel-oneapi-mpi@2021.11.0
 
-# ./pull_libs.py
+export CC=mpiicx
+# export CXX=mpiicpx
+export FC=mpiifx
+# can also be mpiifx
 
-# spack load cmake@3.24.4
-# spack load intel-oneapi-compilers@2024.0.1/xbteted 
-# spack load intel-oneapi-mkl@2024.0.0
-# spack load intel-oneapi-mpi@2021.11.0
-# # spack load hdf5 ~shared
-# spack load hdf5/fxhrrhv
-# export CC=mpiicx
-# export FC=mpiifx
+# spack load netcdf-c
+# spack load netcdf-fortran
 
-export H5DIR=$(spack location -i hdf5)
-export CURLDIR=$(spack location -i curl)
-export XML2DIR=$(spack location -i libxml2)
-export OMPDIR="$(spack location -i intel-oneapi-compilers@2024.0.2)/compiler/latest/"
+# export NETCDF_ROOT=$(spack location -i netcdf-fortran)
+current_dir=$(pwd)
+export NETCDF_ROOT=$current_dir/netcdf
+# export NETCDF_ROOT=$(spack location -i netcdf-fortran)
+# export NETCDF_ROOT=/data/asc24caeporo/xyw/finalasc/raw/netcdf/testf
+# export NETCDF_ROOT=/data/asc24caeporo/xyw/finalasc/raw/netcdf/testc
+export LD_LIBRARY_PATH=/opt/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/hdf5-1.14.3-fxhrrhv46hbchtxp255okz3r2dotzmng/lib:$LD_LIBRARY_PATH
 
-echo "CC: $CC"
-echo "FC: $FC"
-echo "F77: $F77"
 
-export NETCDF_ROOT="$(pwd)/netcdf"
-# export GPTL_ROOT="$(pwd)/gptl"
+# export FC=gfortran
 
-if [ x"$1" = xrebuild ]; then
-  rm -rf build
-fi
+# export CFLAGS=${CFLAGS:=}
+# export FFLAGS=${FFLAGS:=}
+# export CXXFLAGS=${CXXFLAGS:=}
 
-if [ ! -d build ]; then
-  cmake -B build -G Ninja 
-fi
+# export CC=icx
+# export CXX=icpx
+# export FC=ifx
 
-# cmake -B build -G Ninja \
-  # -DCMAKE_RANLIB=/data/spack/opt/spack/linux-ubuntu22.04-icelake/gcc-11.4.0/intel-oneapi-compilers-2024.0.2-lvfe6ufintzu3ibq3loire4oz62soeqe/compiler/2024.0/bin/compiler/llvm-ranlib
+# add_compile_flag() {
+# 	export CFLAGS="$1 $CFLAGS"
+# 	export CXXFLAGS="$1 $CXXFLAGS"
+# 	export FFLAGS="$1 $FFLAGS"
+# }
 
+# openmpi_base="/usr/mpi/gcc/openmpi-4.0.3rc4/"
+# export PKG_CONFIG_PATH="${openmpi_base}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
+# openmpi_FLAGS=$(pkg-config --cflags ompi)
+# openmpi_LIBS="$(pkg-config --static --libs ompi)"
+# add_compile_flag "${openmpi_FLAGS}"
+# export LDFLAGS="${openmpi_LIBS} ${LDFLAGS:-}"
+
+rm -rf build
+cmake -B build
 cd build
-#make -j64 VERBOSE=1
-ninja
+make -j8
diff --git a/build_gptl.sh b/build_gptl.sh
deleted file mode 100755
index 1c034d9..0000000
--- a/build_gptl.sh
+++ /dev/null
@@ -1,35 +0,0 @@
-#!/bin/bash
-
-# spack load cmake@3.24.4
-# spack load intel-oneapi-compilers@2024.0.1/xbteted
-# spack load intel-oneapi-mpi@2021.11.0
-# # spack load libxml2
-# # spack load curl
-# # spack load hdf5/fxhrrhv
-# export CC=mpiicx
-# export FC=mpiifort
-# export F77=mpiifort
-
-cd "$(dirname $0)" || exit 1
-set -e
-
-source ./env.sh
-
-gptl_dir="$(pwd)/gmcore/gptl"
-mkdir -p "$gptl_dir"
-
-cd "$gptl_dir"
-pwd
-if [ x"$(ls -A .)" = x"" ]; then
-  # wget https://github.com/jmrosinski/GPTL/releases/download/v8.1.1/gptl-8.1.1.tar.gz
-  # tar -zxvf gptl-8.1.1.tar.gz
-  cp -r /data/gomars_data/gptl-8.1.1 .
-fi
-
-cd gptl-8.1.1
-# wget https://gist.githubusercontent.com/bonfus/21dec6b966859f5f509b935f8b055a7f/raw/macros.make
-./configure --enable-pmpi --disable-openmp --prefix=$gptl_dir
-# make check
-make install
-# wordaround: remove all shared libs
-rm ${gptl_dir}/lib/*so*
diff --git a/build_netcdf.sh b/build_netcdf.sh
index b10c3b2..7a3d0ef 100755
--- a/build_netcdf.sh
+++ b/build_netcdf.sh
@@ -1,70 +1,71 @@
 #!/bin/bash
 
-set -e
-cd "$(dirname $0)" || exit 1
-
-# spack load cmake@3.24.4
-# spack load intel-oneapi-compilers@2024.0.1/xbteted
-# spack load intel-oneapi-mkl@2024.0.0
-# spack load intel-oneapi-mpi@2021.11.0
-# spack load libxml2/q66mtbb
-# spack load curl
-# # spack load hdf5 ~shared
-# spack load hdf5/fxhrrhv
-# export CC=mpiicx
-# export FC=mpiifort
-# export F77=mpiifort
-
-source ./env.sh
+spack load cmake@3.24.4
+spack load intel-oneapi-compilers@2024.0.1
+spack load intel-oneapi-mkl@2024.0.0
+spack load intel-oneapi-mpi@2021.11.0
+spack load libxml2
+spack load curl
+spack load hdf5/fxhrrhv
 
 export CC=mpiicx
 export FC=mpiifort
 export F77=mpiifort
 
-# export H5DIR=$(spack location -i hdf5 ~shared)
-export H5DIR=$(spack location -i hdf5)
+# export H5DIR=$(spack location -i hdf5)
+export H5DIR=$(spack location -i hdf5/fxhrrhv)
 export CURLDIR=$(spack location -i curl)
 export XML2DIR=$(spack location -i libxml2)
 
-(
-  cd gmcore
-  netcdf_dir="netcdf"
-  mkdir -p "$netcdf_dir"
-  cd "$netcdf_dir"
-  if [ x"$(ls -A .)" = x"" ]; then
-    # wget https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
-    # wget https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.1/netcdf-fortran-4.6.1.tar.gz
-    # tar -zxvf netcdf-c-4.9.2.tar.gz
-    # tar -zxvf netcdf-fortran-4.6.1.tar.gz
-    cp -r /data/gomars_data/netcdf/netcdf-c-4.9.2/ .
-    cp -r /data/gomars_data/netcdf/netcdf-fortran-4.6.1/ .
-  fi
-)
-
 cd ./gmcore/netcdf/
 
 CDIR=$(pwd)
 NCDIR=$(pwd)
 NFDIR=$(pwd)
 
-(
-  cd ./netcdf-c-4.9.2
-  [ -f Makefile ] || [ -f makefile ] && make clean
-  CPPFLAGS="-I${H5DIR}/include -I${CURLDIR}/include -I${XML2DIR}/include " \
+cd netcdf-c-4.9.2/
+
+# export H5DIR=/opt/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/hdf5-1.14.3-fxhrrhv46hbchtxp255okz3r2dotzmng/
+# export CURLDIR=/opt/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/curl-8.4.0-aayk5xsgzo4cyrucspyiktovue3mwnmn/
+
+CPPFLAGS="-I${H5DIR}/include -I${CURLDIR}/include -I${XML2DIR}/include" \
     LDFLAGS="-L${H5DIR}/lib -L${CURLDIR}/lib -L${XML2DIR}/lib" \
-    ./configure --prefix=${CDIR} --enable-parallel-tests
-  make install -j64
-)
+    ./configure --prefix=${CDIR}
+
+# CPPFLAGS="-I${H5DIR}/include" \
+#     LDFLAGS="-L${H5DIR}/lib" \
+#     ./configure --enable-parallel --disable-byterange
+# make check
+make install
+
+# cd ./gmcore/netcdf/
+
+# CDIR=$(pwd)
+
 
 export LD_LIBRARY_PATH=${NCDIR}/lib:${LD_LIBRARY_PATH}
-(
-  cd ./netcdf-fortran-4.6.1
-  [ -f Makefile ] || [ -f makefile ] && make clean
-  CPPFLAGS="-I${NCDIR}/include -I${H5DIR}/include -I${CURLDIR}/include -I${XML2DIR}/include" \
+
+# echo $NFDIR
+
+cd  ../netcdf-fortran-4.6.1/
+
+# make distclean
+make clean
+CPPFLAGS="-I${NCDIR}/include -I${H5DIR}/include -I${CURLDIR}/include -I${XML2DIR}/include" \
     LDFLAGS="-L${NCDIR}/lib -L${H5DIR}/lib -L${CURLDIR}/lib -L${XML2DIR}/lib" \
-    ./configure --prefix=${NFDIR} --enable-parallel-tests
-  make install -j64
-)
+    ./configure --prefix=${NFDIR}
+
+# CPPFLAGS="-I${NCDIR}/include" \
+#     LDFLAGS="-L${NCDIR}/lib" \ 
+#     ./configure --prefix=${NFDIR}
+
+# CPPFLAGS="-I${NCDIR}/include"
+# LDFLAGS="-L${NCDIR}/lib" 
+# ./configure --prefix=${NFDIR}
+
+# CPPFLAGS="-I${NCDIR}/include" \
+#     LDFLAGS="-L${NCDIR}/lib" \ 
+#     ./configure --prefix=${NFDIR}
 
-# wordaround: remove all shared libs
-rm ./lib/*so*
+# make check
+make install
\ No newline at end of file
diff --git a/env.sh b/env.sh
deleted file mode 100755
index de837ee..0000000
--- a/env.sh
+++ /dev/null
@@ -1,13 +0,0 @@
-spack load cmake
-spack load ninja@1.11.1
-spack load intel-oneapi-compilers@2024.0.2
-spack load intel-oneapi-mpi@2021.11.0
-spack load intel-oneapi-mkl@2024.0.0
-spack load libxml2
-spack load curl
-spack load hdf5
-# spack load hdf5 ~shared
-
-export CC=mpiicx
-export FC=mpiifx
-export F77=mpiifx
diff --git a/find_so_deps.sh b/find_so_deps.sh
deleted file mode 100755
index 8a20504..0000000
--- a/find_so_deps.sh
+++ /dev/null
@@ -1,25 +0,0 @@
-#!/bin/bash
-
-# Usage:
-#   ./find_so_deps.sh path/to/execuable
-# Example:
-#   ./find_so_deps.sh ./gmcore/build/gmcore_driver.exe
-
-find_dependenies() {
-    name="$1"
-    path="$2"
-    echo "$name:"
-    objdump -x "$path" 2> /dev/null | grep NEEDED | grep -oP '(  \S*\.so\S*)'
-}
-
-[ -n "$1" ] || exit 1
-[ -f "$1" ] || exit 1
-
-IFS=$'\n'
-
-for line in $(ldd "$1" | grep -oP '.*\.so.* => \/.* \(.*\)'); do
-    name="$(grep -oP '(lib.*so.*)(?= =>)' <<< "$line")"
-    path="$(grep -oP '(?<= => )(.*so.*)(?= \()' <<< "$line")"
-    find_dependenies "$name" "$path"
-done
-
diff --git a/gomars_all_in_one.sh b/gomars_all_in_one.sh
index 0bb7693..0271e92 100755
--- a/gomars_all_in_one.sh
+++ b/gomars_all_in_one.sh
@@ -1,78 +1,42 @@
 #!/bin/bash
 
 current_dir=$(pwd)
-current_user=$(whoami)
-echo "The current user is: $current_user"
 
-if [ ! -d "/data/gomars_output/$current_user" ]; then
-    mkdir "/data/gomars_output/$current_user"
-    echo "The current user is: me"
+if [ ! -d "$current_dir/gmcore" ]; then
+    git clone https://gitee.com/dongli85/GMCORE gmcore
 fi
 
-# if [ ! -d "$current_dir/gmcore" ]; then
-#     git clone https://gitee.com/dongli85/GMCORE gmcore
-# fi
-
 # pull libs and data
 pushd gmcore
 if [ ! -d "$current_dir/gmcore/data" ]; then
-    mkdir "$current_dir/gmcore/data"
-    # ./pull_data.py -p earth
-    # ./pull_data.py -p mars
-    cp -r /data/gomars_data/data/* ./data/
-fi
-
-if [ ! -d "$current_dir/gmcore/output" ]; then
-    # ./pull_libs.py
-    mkdir output/
-fi
-
-popd
-if [ ! -d "$current_dir/gmcore/netcdf" ]; then
-    ./build_netcdf.sh
-fi
-
-if [ ! -d "$current_dir/gmcore/gptl" ]; then
-    ./build_gptl.sh
+    ./pull_libs.py
+    ./pull_data.py -p earth
+    ./pull_data.py -p mars
 fi
 
-if [ ! -d "$current_dir/gmcore/run" ]; then
-    pushd gmcore
-    mkdir run
-    pushd run
-    # git clone https://gitee.com/dongli85/GMCORE-TESTBED
-    cp -r /data/gomars_data/namelist/GMCORE-TESTBED .
+# download netcdf
+current_dir=$(pwd)
+target_dir="$current_dir/netcdf"
+if [ ! -d "$target_dir" ]; then
+    mkdir -p "$target_dir"
+    echo "created: $target_dir"
+    pushd netcdf
+    wget https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
+    wget https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.1/netcdf-fortran-4.6.1.tar.gz
+    tar -zxvf netcdf-c-4.9.2.tar.gz
+    tar -zxvf netcdf-fortran-4.6.1.tar.gz
     popd
     popd
+    # build netcdf
+    # build netcdf-c && netcdf-fortran
+    # cd ../../
+    ./build_netcdf.sh
+else
+    echo "existed: $target_dir"
+    popd
 fi
-# download netcdf
-# current_dir=$(pwd)
-# target_dir="$current_dir/netcdf"
-# if [ ! -d "$target_dir" ]; then
-#     mkdir -p "$target_dir"
-#     echo "created: $target_dir"
-#     pushd netcdf
-#     wget https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
-#     wget https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.1/netcdf-fortran-4.6.1.tar.gz
-#     tar -zxvf netcdf-c-4.9.2.tar.gz
-#     tar -zxvf netcdf-fortran-4.6.1.tar.gz
-#     popd
-#     popd
-#     # build netcdf
-#     # build netcdf-c && netcdf-fortran
-#     # cd ../../
-#     ./build_netcdf.sh
-# else
-#     echo "existed: $target_dir"
-#     popd
-# fi
 
 
 
 # build gmcore
 ./build_gmcore.sh
-
-# pushd gmcore
-
-# run a case
-# sbatch select_run.sh -c swm_mz.360x180
diff --git a/install_conda.sh b/install_conda.sh
deleted file mode 100755
index 05acc8b..0000000
--- a/install_conda.sh
+++ /dev/null
@@ -1,4 +0,0 @@
-#!/bin/bash
-mkdir -p ~/miniconda3
-wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
-bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
\ No newline at end of file
diff --git a/ramfsrun.sh b/ramfsrun.sh
deleted file mode 100755
index 0517a16..0000000
--- a/ramfsrun.sh
+++ /dev/null
@@ -1,158 +0,0 @@
-#!/bin/bash
-#SBATCH -N 1
-#SBATCH -n 64
-#SBATCH -w hepnode1
-#SBATCH --exclusive
-#SBATCH --output=./output/slurm-%j.out
-
-
-echo "******ramfsrun.sh*******"
-cat $0
-echo "******ramfsrun.sh*******"
-
-
-source ./env.sh
-
-message=$2
-# days=$3
-
-run ( ) {
-	
-	if [ $(hostname) != "hepnode0" ]; then
-		# export I_MPI_FABRICS=ofa:ofa
-		export I_MPI_SHM=clx_avx512
-		export UCX_RC_PATH_MTU=4096
-		export I_MPI_HYDRA_RMK=slurm
-		export I_MPI_PIN=off
-		export OMP_NUM_THREADS=1
-		# export I_MPI_ASYNC_PROGRESS=1
-		export I_MPI_DEBUG=10
-		export I_MPI_VAR_CHECK_SPELLING=1
-		# export I_MPI_INTRANODE_EAGER_THRESHLOD=1024
-  		# export FI_OFI_RXM_RX_SIZE=4096
-  		# export FI_OFI_RXM_TX_SIZE=4096
-
-
-		# export I_MPI_STATS=2
-		# export I_MPI_SHM_HEAP_CSIZE=-1
-		# export I_MPI_CACHE_BYPASS=1
-		# export I_MPI_WAIT_=1
-		# export I_MPI_SHM_NUM_BUFFERS=-1
-		# export I_MPI_SHM_BUFFER_SIZE=102400
-		export I_MPI_MALLOC=1
-	fi
-	export I_MPI_PIN=off
- 	# export I_MPI_FILESYSTEM=1
-	# export I_MPI_FILESYSTEM_CB_CONFIG_LIST="hepnode0:1" 
-	# export I_MPI_FILESYSTEM_CB_CONFIG_LIST="*:1"
-
-	# export I_MPI_SHM_HEAP=1
-	# export I_MPI_PLATFORM=auto
-	# export I_MPI_TUNING_MODE=auto
-	# export I_MPI_TUNING_AUTO_POLICY=max
-	# export I_MPI_EAGER_THRESHOLD=1024
-	case_name=$1
-	node=$2
-	proc=$3
-	# days=$4
-	
-
-	adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
-	swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
-	normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)
-
-	prefix="./run/GMCORE-TESTBED/"
-	suffix="/namelist"
-	namelist_relative_path="${prefix}${case_name}/${suffix}"
-	namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
-	data_path="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")"
-	# now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-	days=$(grep 'run_days' ${namelist_absolute_path} | sed 's/.*= *\([0-9]*\).*/\1/')
-
-	check_file="/data/gomars_output/public/N${2}n${3}/${case_name}_${days}days/baseline.nc" 
-	check_dir="/data/gomars_output/public/N${2}n${3}/" 
-
-	if [ -f "$check_file" ]; then
-		echo "exist!"
-	else
-		echo "not exist, use default path at N1n16!"
-		echo "!!! You should notice days!"
-		check_file="/data/gomars_output/public/N1n16/${case_name}_${days}days/baseline.nc"
-	fi
-	cd ..
-	current_dir=$(pwd)
-	# mpitune_fast -f ./hostfile
-	cd gmcore/
-
-	if [ ! -d ${data_path} ]; then    
-		mkdir -p ${data_path}
-	fi
-
-	# pushd ${data_path}
-	# pushd /ramdisk
-	pushd /data/tmpfs
-	# pushd /tmp
-	# touch xywnmmsl
-	# if [[ $namelist_absolute_path == *"adv"* ]]; then
-	# 	echo "adv case"
-	# 	exe_absolute_path=$adv_exe_absolute_path
-	# elif [[ $namelist_absolute_path == *"swm"* ]]; then
-	# 	echo "swm case"
-	# 	exe_absolute_path=$swm_exe_absolute_path
-	# else
-	# 	echo "normal case"
-	# 	exe_absolute_path=$normal_exe_absolute_path
-	# fi
-	# bash -c "mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path" #doesn't work
-	# mpirun -n $3 -ppn $( expr $3 / $2 ) $exe_absolute_path $namelist_absolute_path
-	# mpitune_fast -hf hosts
-	# mpirun -genv I_MPI_PIN_PROCESSOR_LIST map=scatter -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-	mpirun -n $3 -ppn $( expr $3 / $2 ) ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-
-	# rm -rf opt.nc
-	# mv *.nc opt.nc
-	now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-
-	popd
-	fd1="fd1=\"${check_file}\""
-	fd2="fd2=\"${now_dir}\""
-
-	echo $fd1 
-	echo $fd2
-	# source activate ncl_stable
-	# if [[ $namelist_absolute_path == *"adv"* ]]; then
-	# 	echo "adv case"
-	# 	exe_absolute_path=$adv_exe_absolute_path
-	# 	ncl ../script/ncl/adv_verify_answer.ncl $fd1 $fd2
-	# elif [[ $namelist_absolute_path == *"swm"* ]]; then
-	# 	echo "swm case"
-	# 	exe_absolute_path=$swm_exe_absolute_path
-	# 	ncl ../script/ncl/swm_verify_answer.ncl $fd1 $fd2
-	# else
-	# 	echo "normal case"
-	# 	exe_absolute_path=$normal_exe_absolute_path
-	# 	ncl ../script/ncl/normal_verify_answer.ncl $fd1 $fd2
-	# fi
-
-	
-}
-
-
-if [ -z $1 ]; then
-	echo "Usage: $0 <case_list> <message>"
-	exit 1
-fi
-
-cmd_array=()
-
-while IFS= read -r cmd
-do
-	cmd_array+=("$cmd")
-done < $1
-
-for cmd in "${cmd_array[@]}"
-do
-	pushd ./gmcore
-	run $cmd $2
-	popd
-done
diff --git a/run_conda.sh b/run_conda.sh
deleted file mode 100755
index df739a6..0000000
--- a/run_conda.sh
+++ /dev/null
@@ -1,2 +0,0 @@
-#!/bin/bash
-~/miniconda3/bin/conda init bash
diff --git a/setup_ncl.sh b/setup_ncl.sh
deleted file mode 100755
index 270eca3..0000000
--- a/setup_ncl.sh
+++ /dev/null
@@ -1,7 +0,0 @@
-#!/bin/bash
-
-conda create -n ncl_stable -c conda-forge ncl
-source activate ncl_stable
-ncl -V
-
-# ncl xxx.ncl
diff --git a/vtune_run.sh b/vtune_run.sh
deleted file mode 100755
index fe9aca2..0000000
--- a/vtune_run.sh
+++ /dev/null
@@ -1,103 +0,0 @@
-#!/bin/bash
-#SBATCH -N 1
-#SBATCH -n 64
-#SBATCH -w hepnode4
-#SBATCH --exclusive
-#SBATCH --output=./output/slurm-%j.out
-spack load intel-oneapi-vtune
-echo "******batch_run.sh*******"
-cat $0
-echo "******batch_run.sh*******"
-source ./env.sh
-message=$2
-# days=$3
-run ( ) {
-    if [ $(hostname) != "hepnode0" ]; then
-        export UCX_RC_PATH_MTU=2048
-        export I_MPI_HYDRA_RMK=slurm
-        export I_MPI_PIN=off
-        export OMP_NUM_THREADS=1
-    fi
-    case_name=$1
-    node=$2
-    proc=$3
-    # days=$4
-    adv_exe_absolute_path=$(readlink -f ./build/gmcore_adv_driver.exe)
-    swm_exe_absolute_path=$(readlink -f ./build/gmcore_swm_driver.exe)
-    normal_exe_absolute_path=$(readlink -f ./build/gmcore_driver.exe)
-    prefix="./run/GMCORE-TESTBED/"
-    suffix="/namelist"
-    namelist_relative_path="${prefix}${case_name}/${suffix}"
-    namelist_absolute_path=$(readlink -f ${namelist_relative_path} )
-    data_path="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")"
-    # now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-    days=$(grep 'run_days' ${namelist_absolute_path} | sed 's/.*= *\([0-9]*\).*/\1/')
-    check_file="/data/gomars_output/public/N${2}n${3}/${case_name}_${days}days/baseline.nc" 
-    check_dir="/data/gomars_output/public/N${2}n${3}/" 
-    if [ -f "$check_file" ]; then
-        echo "exist!"
-    else
-        echo "not exist, use default path at N1n16!"
-        echo "!!! You should notice days!"
-        check_file="/data/gomars_output/public/N1n16/${case_name}_${days}days/baseline.nc"
-    fi
-    cd ..
-    current_dir=$(pwd)
-    cd gmcore/
-    if [ ! -d ${data_path} ]; then    
-        mkdir -p ${data_path}
-    fi
-    pushd ${data_path}
-    if [[ $namelist_absolute_path == *"adv"* ]]; then
-        echo "adv case"
-        exe_absolute_path=$adv_exe_absolute_path
-    elif [[ $namelist_absolute_path == *"swm"* ]]; then
-        echo "swm case"
-        exe_absolute_path=$swm_exe_absolute_path
-    else
-        echo "normal case"
-        exe_absolute_path=$normal_exe_absolute_path
-    fi
-    mpirun -n $3 -ppn $( expr $3 / $2 ) \
-    vtune -collect hotspot -knob enable-stack-collection=true \
-    -knob sampling-mode=hw -knob stack-size=0 -knob sampling-interval=100 \
-    -r /home/xyw/ascgomars/ASC24-Final-GoMars/vtune/${case_name}_r001.hs -finalization-mode=deferred -- \
-    ${current_dir}/bind_cpu.sh $exe_absolute_path $namelist_absolute_path
-    rm -rf opt.nc
-    mv *.nc opt.nc
-    now_dir="/data/gomars_output/$(whoami)/${case_name}/N${2}n${3}/"${message}"-$(date +"%y-%m-%d")/opt.nc"
-    popd
-    fd1="fd1=\"${check_file}\""
-    fd2="fd2=\"${now_dir}\""
-    echo $fd1 
-    echo $fd2
-    source activate ncl_stable
-    if [[ $namelist_absolute_path == *"adv"* ]]; then
-        echo "adv case"
-        exe_absolute_path=$adv_exe_absolute_path
-        ncl ../script/ncl/adv_verify_answer.ncl $fd1 $fd2
-    elif [[ $namelist_absolute_path == *"swm"* ]]; then
-        echo "swm case"
-        exe_absolute_path=$swm_exe_absolute_path
-        ncl ../script/ncl/swm_verify_answer.ncl $fd1 $fd2
-    else
-        echo "normal case"
-        exe_absolute_path=$normal_exe_absolute_path
-        ncl ../script/ncl/normal_verify_answer.ncl $fd1 $fd2
-    fi
-}
-if [ -z $1 ]; then
-    echo "Usage: $0 <case_list> <message>"
-    exit 1
-fi
-cmd_array=()
-while IFS= read -r cmd
-do
-    cmd_array+=("$cmd")
-done < $1
-for cmd in "${cmd_array[@]}"
-do
-    pushd ./gmcore
-    run $cmd $2
-    popd
-done
\ No newline at end of file
