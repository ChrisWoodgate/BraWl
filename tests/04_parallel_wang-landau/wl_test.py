import netCDF4 as nc
import numpy as np
import sys

filename = "data/wl_dos.dat"
logdos_run = nc.Dataset(filename)
logdos_run = np.array(logdos_run["grid data"][:], dtype=np.float64)

filename = "../99_ref/wl_dos.nc"
logdos_ver = nc.Dataset(filename)
logdos_ver = np.array(logdos_ver["grid data"][:], dtype=np.float64)

rmse_per = np.sqrt(np.mean((logdos_ver-logdos_run)**2))/np.mean(logdos_ver)

if rmse_per < 0.01:
    print("Pass")
    sys.exit(0)
else:
    print("Fail")
    sys.exit(1)
