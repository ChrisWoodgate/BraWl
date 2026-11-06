#!/bin/bash

# Run every Wang-Landau convergence run without slurm.
# Usage: BRAWL=/path/to/brawl.run ./run_wl.sh
# Optional: PROCS (MPI processes per run, default 16), PARALLEL (runs at once, default 3)

if [ -z "$BRAWL" ] || [ ! -x "$BRAWL" ]; then
    echo "Set BRAWL to the path of brawl.run, e.g. BRAWL=/path/to/brawl.run ./run_wl.sh"
    exit 1
fi

PROCS=${PROCS:-16}
PARALLEL=${PARALLEL:-3}
BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$BASE_DIR/logs"

run_one() {
    d="$1"
    cd "$BASE_DIR/$d" || { echo "Failed to enter $BASE_DIR/$d"; return 1; }
    mkdir -p data asro
    echo "Running $d"
    # No core binding, otherwise the concurrent runs all pin to the same cores
    mpirun -n "$PROCS" --bind-to none "$BRAWL" > "$BASE_DIR/logs/$(echo "$d" | tr '/' '_').txt" 2>&1
    echo "Finished $d"
}
export -f run_one
export BASE_DIR PROCS BRAWL

xargs -P "$PARALLEL" -I{} bash -c 'run_one {}' < "$BASE_DIR/run_dirs.txt"
