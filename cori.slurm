#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 64
#SBATCH -C haswell
#SBATCH -t 00:05:00
#SBATCH -J mel
#SBATCH --exclusive
####SBATCH --contiguous
#SBATCH --mail-type=ALL
#SBATCH -A <> 

# set -e

ulimit -c unlimited

echo "Running on:"
echo $SLURM_NODELIST

export TASKS_PER_SOCKET=16
export CORES_PER_SOCKET=16
export SOCKETS_PER_NODE=2
export CORES_PER_TASK=$(( $CORES_PER_SOCKET / $TASKS_PER_SOCKET ))

cd $SCRATCH
mkdir -p job.$SLURM_JOB_ID
cd job.$SLURM_JOB_ID

export SLURM_ARGS="--ntasks-per-node=$(( $TASKS_PER_SOCKET * $SOCKETS_PER_NODE )) \
    --cpus-per-task=$CORES_PER_TASK \
    --cpu_bind=cores \
    --cpu-freq=High \
    --mem-bind=local"

export GASNET_BACKTRACE=1
export UPCXX_SEGMENT_MB=512

# # To run with geometric graphs
# for NVERTICES in $((8 * 1024 * 1024)); do
#     for EXE in match match_rma; do
#         echo "===== $EXE ====="
#         srun $SLURM_ARGS $HOME/upcxx_matching/upcxx_matching/$EXE -n $NVERTICES -b
#         echo
#     done
# done


# To run with social network graphs
# for GRAPH in com-Orkut.bin com-friendster.ungraph.bin; do
for GRAPH in com-friendster.ungraph.bin; do
    # for EXE in match_rpc match_rma match_nrm match_p2p; do
    for EXE in match_rma; do
        echo "===== $EXE $GRAPH ====="
        sbcast --compress=lz4 $HOME/upcxx_matching/upcxx_matching/$EXE /tmp/$EXE
        srun $SLURM_ARGS /tmp/$EXE -f $SCRATCH/match-graphs/$GRAPH
        echo
    done
done
