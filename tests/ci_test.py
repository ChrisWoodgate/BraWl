import numpy as np
import sys
import os
import netCDF4 as nc

cwd = os.path.basename(os.getcwd())
test_dir = f"../99_ref/{cwd}"
ref_names = list(os.listdir(test_dir))
ref_files = [os.path.join(test_dir, f) for f in ref_names]

test_files = []
for root, dirs, files in os.walk("."):
  for file in files:
    if file in ref_names:
      test_files.append(os.path.join(root, os.path.basename(file)))

test_files = sorted(test_files, key=os.path.basename)
ref_files = sorted(ref_files, key=os.path.basename)

for i, _ in enumerate(ref_files):
  file_name = os.path.basename(ref_files[i])
  # --- Case 1: NetCDF (.nc) files ---
  if file_name.endswith(".nc"):
    try:
      with nc.Dataset(ref_files[i], "r") as ref_nc, nc.Dataset(test_files[i], "r") as test_nc:
        # Check dimensions
        if ref_nc.dimensions.keys() != test_nc.dimensions.keys():
          print(f"[DIM MISMATCH] {os.path.basename(test_files[i])}")
          sys.exit(1)

        # Check variable names
        if ref_nc.variables.keys() != test_nc.variables.keys():
          print(f"[VAR MISMATCH] {os.path.basename(test_files[i])}")
          sys.exit(1)

        # Check variable values
        for var_name in ref_nc.variables:
          ref_var = ref_nc.variables[var_name][:]
          test_var = test_nc.variables[var_name][:]

          # if cwd is 04_parallel_wang-landau
          if cwd == "04_parallel_wang-landau":
            ref_flat = ref_var.flatten()
            test_flat = test_var.flatten()

            nrmse = np.sqrt(np.mean((ref_flat - test_flat) ** 2))/np.mean(np.abs(ref_flat))

            if nrmse > 0.01:
              print(f"[DATA MISMATCH NRMSE >1%] {os.path.basename(test_files[i])}, variable: {var_name}")
              sys.exit(1)
          else:
            if not np.array_equal(ref_var, test_var):
              print(f"[DATA MISMATCH] {os.path.basename(test_files[i])}, variable: {var_name}")
              sys.exit(1)

    except Exception as e:
      print(f"[FILE ERROR] {file_name}: {e}")
      sys.exit(1)

  # --- Case 2: Text (.dat) files ---
  elif file_name.endswith(".dat"):
    try:
      with open(ref_files[i], 'r') as ref_f, open(test_files[i], 'r') as test_f:
        ref_lines = ref_f.readlines()
        test_lines = test_f.readlines()

        if ref_lines != test_lines:
          print(f"[DATA MISMATCH] {file_name} (differs)")
          sys.exit(1)

    except Exception as e:
      print(f"[FILE ERROR] {file_name}: {e}")
      sys.exit(1)

  # --- Case 3: Nested sampling (.energies) files ---
  # One header line (n_walkers, 1, 0, 'False', n_at) followed by one
  # line per iteration holding the iteration number and the energy
  # limit. Serial nested sampling with a fixed seed is deterministic,
  # so the header must match token for token and the energies exactly.
  elif file_name.endswith(".energies"):
    try:
      with open(ref_files[i], 'r') as ref_f, open(test_files[i], 'r') as test_f:
        ref_header = ref_f.readline().split()
        test_header = test_f.readline().split()

      if ref_header != test_header:
        print(f"[HEADER MISMATCH] {file_name}")
        sys.exit(1)

      ref_data = np.loadtxt(ref_files[i], skiprows=1, ndmin=2)
      test_data = np.loadtxt(test_files[i], skiprows=1, ndmin=2)

      if ref_data.shape != test_data.shape:
        print(f"[SHAPE MISMATCH] {file_name}: {ref_data.shape} vs {test_data.shape}")
        sys.exit(1)

      if not np.array_equal(ref_data, test_data):
        n_diff = int(np.sum(ref_data[:, 1] != test_data[:, 1]))
        max_diff = float(np.max(np.abs(ref_data[:, 1] - test_data[:, 1])))
        print(f"[DATA MISMATCH] {file_name}: {n_diff} iterations differ, max |dE| = {max_diff:.3e}")
        sys.exit(1)

    except Exception as e:
      print(f"[FILE ERROR] {file_name}: {e}")
      sys.exit(1)

  # --- Any other extension is not compared ---
  else:
    print(f"[SKIPPED] {file_name}: no comparison rule for this file type")

print("Pass")
sys.exit(0)
