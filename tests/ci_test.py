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

          ref_flat = ref_var.flatten()
          test_flat = test_var.flatten()

          nrmse = np.sqrt(np.mean((ref_flat - test_flat) ** 2))/np.mean(np.abs(ref_flat))

          if nrmse > 0.01:
            print(f"[DATA MISMATCH NRMSE >1%] {os.path.basename(test_files[i])}, variable: {var_name}")
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

print("Pass")
sys.exit(0)
