#!/bin/bash

# Rerun the copper production Wang-Landau runs (Cu0.00 ... Cu0.95) without slurm.
# Usage: BRAWL=/path/to/brawl.run ./run_wl.sh [Cu0.00 Cu0.50 ...]
# With no arguments every Cu* directory is run.
# Optional: PROCS (MPI processes per run, default 16), PARALLEL (runs at once, default 3),
#           ARCHIVE (1 = move existing data/asro/load_balance output aside first, default 1)

if [ -z "$BRAWL" ] || [ ! -x "$BRAWL" ]; then
    echo "Set BRAWL to the path of brawl.run, e.g. BRAWL=/path/to/brawl.run ./run_wl.sh"
    exit 1
fi

PROCS=${PROCS:-16}
PARALLEL=${PARALLEL:-3}
ARCHIVE=${ARCHIVE:-1}
BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BASE_DIR/logs"

if [ "$#" -gt 0 ]; then
    DIRS=("$@")
else
    DIRS=()
    for d in "$BASE_DIR"/Cu[0-9]*/; do
        DIRS+=("$(basename "$d")")
    done
fi

run_one() {
    d="$1"
    cd "$BASE_DIR/$d" || { echo "Failed to enter $BASE_DIR/$d"; return 1; }

    # Keep the previous results so they can be compared with the new ones
    if [ "$ARCHIVE" = "1" ]; then
        for sub in data asro load_balance; do
            if [ -d "$sub" ] && [ -n "$(ls -A "$sub" 2>/dev/null)" ]; then
                mkdir -p "previous_$STAMP"
                mv "$sub" "previous_$STAMP/$sub"
            fi
        done
    fi

    mkdir -p data asro load_balance
    echo "Running $d ($PROCS processes)"
    # No core binding, otherwise concurrent runs all pin to the same cores
    mpirun -n "$PROCS" --bind-to none "$BRAWL" > "$BASE_DIR/logs/${d}_$STAMP.txt" 2>&1
    echo "Finished $d"
}
export -f run_one
export BASE_DIR PROCS BRAWL ARCHIVE STAMP

printf '%s\n' "${DIRS[@]}" | xargs -P "$PARALLEL" -I{} bash -c 'run_one {}'
