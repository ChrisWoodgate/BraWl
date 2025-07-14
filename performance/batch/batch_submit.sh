#!/bin/bash

# List of core values (from your 'cores' array)
windows=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16)
folders=(1 2 3)
brawl_type=(0 1 2 3 4 5)
walkers=(01 02 03 04 05 06)
overlaps=(05 10 25 50 75)

walkers=(01)
overlaps=(50)
brawl_type=(0)

# Loop over each value in the cores array
for brawl in "${brawl_type[@]}"
do
  for window in "${windows[@]}"
  do
      for folder in "${folders[@]}"
      do
          for walker in "${walkers[@]}"
          do
            for overlap in "${overlaps[@]}"
            do
              # Create a new batch file name based on the current core value
              batch_file="batch_${brawl}_${walker}_${window}_${overlap}_${folder}.sh"

              # Submit the job with the modified batch file name
              sbatch "$batch_file"
            done
          done
      done
  done
done
