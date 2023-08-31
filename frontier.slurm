#!/bin/bash
#SBATCH -J UPXML
#SBATCH -N 64
#SBATCH -p batch
#SBATCH -t 02:00:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# download MPI versions from: https://github.com/sg0/mel
# to compare with UPCXX RMA, build MPI version with -DUSE_MPI_RMA 
export BIN_PATH_RMA="$HOME/projs/mel-sg0"

export BIN_PATH="$HOME/projs/mel-upx"
export DATA_PATH="/lustre/orion/csc318/world-shared/graph_inputs/social"

export FI_CXI_RX_MATCH_MODE=hybrid

export UPCXX_SHARED_HEAP_SIZE=1GB

export GASNET_SPAWN_CONTROL=mpi
export GASNET_MAX_SEGSIZE=1GB/P
export GASNET_OFI_RECEIVE_BUFF_SIZE=recv

n=64
for ppn in 8 16 32; do 
pes=$(($n*$ppn))
for file in $DATA_PATH/*; do
echo "Half-approx Graph Matching using MPI RMA"
echo "------------------------------------------"
srun -N$n -n$pes -c1 --cpu-bind=threads --threads-per-core=1 $BIN_PATH_RMA/./match -b -f $file
echo "Half-approx Graph Matching using UPC++ RMA"
echo "------------------------------------------"
srun -N$n -n$pes -c1 --cpu-bind=threads --threads-per-core=1 $BIN_PATH/./match_rma -b -f $file
echo "Half-approx Graph Matching using UPC++ RPC"
echo "------------------------------------------"
srun -N$n -n$pes -c1 --cpu-bind=threads --threads-per-core=1 $BIN_PATH/./match_rpc -b -f $file
done
done
