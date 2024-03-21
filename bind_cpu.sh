#!/bin/bash
# LOCAL_RANK=$SLURM_LOCALID
LOCAL_RANK=$MPI_LOCALRANKID # for Intel MPI
# LOCAL_SIZE=$SLURM_TASKS_PER_NODE 
# LOCAL_SIZE=${LOCAL_SIZE//(x3)/}
LOCAL_SIZE=$MPI_LOCALNRANKS # for Intel MPI
NCPUS=$(nproc --all) # Number of logical cores
NUM_NUMA=4

echo NCPUS=$NCPUS
echo LOCAL_SIZE=$LOCAL_SIZE
echo CORES_PER_PROCESS=$CORES_PER_PROCESS
echo LOCAL_RANK=$LOCAL_RANK
CORES_PER_PROCESS=$(($NCPUS / $LOCAL_SIZE))
HARDCORES_PER_PROCESS=$(($CORES_PER_PROCESS / 2))

NUMA_ID=$(($NUM_NUMA * $LOCAL_RANK / $LOCAL_SIZE))

CORE_START1=$(( $HARDCORES_PER_PROCESS * $LOCAL_RANK ))
CORE_END1=$(( $HARDCORES_PER_PROCESS * $LOCAL_RANK + $HARDCORES_PER_PROCESS - 1 ))
# CORE_START2=$(( $HARDCORES_PER_PROCESS * $LOCAL_RANK + 64 ))
# CORE_END2=$(( $HARDCORES_PER_PROCESS * $LOCAL_RANK + $HARDCORES_PER_PROCESS + 63))

CORES=$(seq -s, $CORE_START1 $CORE_END1)
# CORES=$CORES",$(seq -s, $CORE_START2 $CORE_END2)"

echo "Process $LOCAL_RANK on $(hostname) bound to cores $CORES"
exec numactl -m "$NUMA_ID" -C "$CORES" $@