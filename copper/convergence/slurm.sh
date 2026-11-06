#!/bin/bash

# 16 runs, 16 MPI processes each. One 48-core node runs 3 runs at once,
# so 6 array tasks cover everything (the last one only has 1 run).

#SBATCH --job-name=CuConv
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3850
#SBATCH --partition=compute
#SBATCH --time=48:00:00
#SBATCH --array=0-5
#SBATCH --account=su007-dq
#SBATCH --output=logs/out_%A_%a.txt
#SBATCH --error=logs/err_%A_%a.txt

module purge
module load GCC/13.2.0 OpenMPI/4.1.6 netCDF-Fortran/4.6.1 OpenBLAS/0.3.24

BASE_DIR="/home/p/phunsc/BraWl/copper/convergence"
RUNS_PER_NODE=3
PROCS_PER_RUN=16
mapfile -t dirs < "$BASE_DIR/run_dirs.txt"

for ((i = 0; i < RUNS_PER_NODE; i++)); do
    idx=$((SLURM_ARRAY_TASK_ID * RUNS_PER_NODE + i))
    [ "$idx" -ge "${#dirs[@]}" ] && break
    d="${dirs[$idx]}"

    echo "Running in $BASE_DIR/$d"
    cd "$BASE_DIR/$d" || { echo "Failed to enter $BASE_DIR/$d"; exit 1; }

    # No core binding, otherwise the concurrent runs all pin to the same cores
    mpirun -n $PROCS_PER_RUN --bind-to none "$BASE_DIR/../../brawl.run" > "$BASE_DIR/logs/run_${idx}.txt" 2>&1 &
done
wait
