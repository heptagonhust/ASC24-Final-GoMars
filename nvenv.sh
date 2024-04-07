./env.sh

NVARCH=`uname -s`_`uname -m`; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
export PATH=$NVCOMPILERS/$NVARCH/24.3/comm_libs/mpi/bin:$PATH
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/24.3/com
echo $MANPATH

which mpirun