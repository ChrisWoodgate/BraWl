#!/bin/bash

# Quick simulated annealing of every run (4 processes each, 12 at once on a
# 48-core node) to find the ground state energy for the WL energy windows.
# Usage: python3 anneal_windows.py setup; sbatch anneal.sh; python3 anneal_windows.py apply

#SBATCH --job-name=CuAnneal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3850
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --account=su007-dq
#SBATCH --output=logs/anneal_%j.txt

module purge
module load GCC/13.2.0 OpenMPI/4.1.6 netCDF-Fortran/4.6.1 OpenBLAS/0.3.24

BASE_DIR="/home/p/phunsc/BraWl/copper/convergence"

run_one() {
    cd "$BASE_DIR/$1/anneal" || exit 1
    mpirun -n 4 --bind-to none "$BASE_DIR/../../brawl.run" > run.log 2>&1
}
export -f run_one
export BASE_DIR

xargs -a "$BASE_DIR/run_dirs.txt" -P 12 -I{} bash -c 'run_one {}'
