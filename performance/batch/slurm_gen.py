import numpy as np
import sys

cores_per_node = 128

windows = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
brawl_type = [0, 1, 2, 3, 4, 5]
walkers = [1, 2, 3, 4, 5, 6]
overlaps = [5, 10, 25, 50, 75]

# Overrides
#walkers = [1]
overlaps = [25]
brawl_type = [0, 2, 4]



base_script_content = '''#!/bin/bash

#SBATCH --output={brawl_type}_{walkers:02d}_{window:02d}_{overlap:02d}_{brawl_folder}.out
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={ntasks_compute}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3850
#SBATCH --partition=compute
#SBATCH --time=2:00:00
#SBATCH --account=su007-dq

module purge
module load GCC/13.2.0 OpenMPI/4.1.6 netCDF-Fortran/4.6.1 OpenBLAS/0.3.24

cd /home/p/phunsc/BraWl/performance/{brawl_type}_{walkers:02d}_{window:02d}_{overlap:02d}_{brawl_folder}

export SLURM_JOB_ID=${{SLURM_JOB_ID}}
srun -n {ntasks} ../../brawl.run
'''
for brawl in brawl_type:
  for walker in walkers:
    for window in windows:
      for overlap in overlaps:
        for brawl_folder in range(1,5+1):
          ntasks = walker*window
          nodes=int(np.ceil(max(ntasks/cores_per_node, 1)))
          script_content = base_script_content.format(nodes=nodes, ntasks=ntasks, ntasks_compute=int(np.ceil(ntasks/nodes)), brawl_type=brawl, walkers=walker, window=window, overlap=overlap, brawl_folder=brawl_folder)

          # Specify the file name
          file_name = 'batch_{brawl_type}_{walkers:02d}_{window:02d}_{overlap:02d}_{brawl_folder}.sh'.format(brawl_type=brawl, walkers=walker, window=window, overlap=overlap, brawl_folder=brawl_folder)

          # Write the content to the file
          with open(file_name, 'w') as file:
              file.write(script_content)

          print(f"The script has been written to {file_name}")

# cd /home/p/phunsc/BraWl/performance/{brawl_type}_{walkers:02d}_{window:02d}_{brawl_folder}
