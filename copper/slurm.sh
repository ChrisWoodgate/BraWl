#!/bin/bash

#SBATCH --job-name=brawl_array
#SBATCH --output=brawl_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3850
#SBATCH --partition=compute
#SBATCH --time=48:00:00
#SBATCH --array=0-19
#SBATCH --account=su007-dq
#SBATCH --job-name=AlNiCu
#SBATCH --output=logs/out_%A_%a.txt
#SBATCH --error=logs/err_%A_%a.txt

module purge
module load GCC/13.2.0 OpenMPI/4.1.6 netCDF-Fortran/4.6.1 OpenBLAS/0.3.24

BASE_DIR="/home/p/phunsc/BraWl/copper"
dirs=(Cu0.00 Cu0.05 Cu0.10 Cu0.15 Cu0.20 Cu0.25 Cu0.30 Cu0.35 Cu0.40 Cu0.45 Cu0.50 Cu0.55 Cu0.60 Cu0.65 Cu0.70 Cu0.75 Cu0.80 Cu0.85 Cu0.90 Cu0.95)

# Get directory for this array task
d="${dirs[$SLURM_ARRAY_TASK_ID]}"

echo "Running in $BASE_DIR/$d"
cd "$BASE_DIR/$d" || { echo "Failed to enter $BASE_DIR/$d"; exit 1; }

mpirun -n 32 ../../brawl.run
