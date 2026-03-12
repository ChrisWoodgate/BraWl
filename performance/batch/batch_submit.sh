#!/bin/bash

# Parameters
windows=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16)
folders=(1 2 3 4 5)
brawl_type=(0 1 2 3 4 5)
walkers=(01 02 03 04 05 06)
overlaps=(00 10 25 50 75)

# Optional overrides
walkers=(01)
#overlaps=(50)
#brawl_type=(0)

# Set job limits
MAX_JOBS=300          # Max jobs allowed in the queue
RESUME_THRESHOLD=200  # Wait until job count drops below this before submitting more
SLEEP_INTERVAL=60     # Seconds between queue checks
BATCH_SIZE=50        # Number of jobs to submit before sleeping

# Function to get current number of jobs in the queue
get_queued_jobs() {
  squeue -u "$USER" --noheader | wc -l
}

submit_counter=0

# Loop over combinations and control submission rate
for brawl in "${brawl_type[@]}"; do
  for window in "${windows[@]}"; do
    for folder in "${folders[@]}"; do
      for walker in "${walkers[@]}"; do
        for overlap in "${overlaps[@]}"; do

          # Wait if too many jobs are already queued
          while true; do
            current_jobs=$(get_queued_jobs)
            echo "Current jobs in queue: $current_jobs"
            if (( current_jobs < MAX_JOBS )); then
              break
            else
              # Wait until job count drops below RESUME_THRESHOLD
	      submit_counter=0
              while (( current_jobs >= RESUME_THRESHOLD )); do
                echo "Queue saturated ($current_jobs jobs), sleeping $SLEEP_INTERVAL seconds..."
                sleep $SLEEP_INTERVAL
                current_jobs=$(get_queued_jobs)
              done
              break
            fi
          done

          # Submit job
          batch_file="batch_${brawl}_${walker}_${window}_${overlap}_${folder}.sh"
          sbatch "$batch_file"
          ((submit_counter++))

          # Sleep after submitting a batch of jobs
          if (( submit_counter % BATCH_SIZE == 0 )); then
            echo "Submitted $submit_counter jobs, sleeping for $SLEEP_INTERVAL seconds..."
	    submit_counter=0
            sleep $SLEEP_INTERVAL
          fi

        done
      done
    done
  done
done
