'''
------------------------------------------------------------------
Example analysis script for bontewarlo simulated annealing output.

C. D. Woodgate, Warwick                                      2023
------------------------------------------------------------------
'''
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append("/home/epp/phuqtp/PhD/code/monte_carlo_analysis/")
from  netCDF4 import Dataset
from pathlib import Path

# Names of elements
elements = ['Ni', 'Co', 'Cr']

# Concentrations of species
concentrations = 1.0/3.0**np.ones(3)

# Arrays into which to extract data
densities_array = []
temperatures_array = []
wcs_array = []
cs_array = []

# Number of calculations in the ensemble
n_calcs = 8

# Path to directory (can be changed)
path = './'

# Loop over calculations
for i in range(n_calcs):
    
    # Read radial density file
    file = path + 'radial_densities/proc_00' + str(i) + '_rho_of_T.nc'
    nc_data_file = Path(file)
    nc_data = Dataset(nc_data_file, "r", format="NETCDF4")

    # Get radial densities
    rho = nc_data.variables['rho data'][:]

    # Get temperatures
    T_data = nc_data.variables['T data'][:]
    
    # Append to the relevant arrays
    densities_array.append(rho)
    temperatures_array.append(T_data)
    
    # Close the file
    nc_data.close()
    
    # Read from the energy diagnostics file
    file = path + 'diagnostics/proc_000' + str(i) + 'energy_diagnostics.dat'
    data = np.genfromtxt(file, dtype=float, delimiter=',', names=True)

    # Get temperature and SHC
    temps = data['T']
    C = data['C']
    
    # Append to the relevant array
    cs_array.append(C)
    
    # Compute the Warren-Cowley ASRO parameters
    wcs = []
    for i in range(len(temps)):
        wcoft = np.zeros((3,3,2))
        for j in range(len(elements)):
            for k in range(len(elements)):
                wcoft[j,k,0] = 1.0 - 1.0/12.0*rho[i,1,j,k]/concentrations[j]
                wcoft[j,k,1] = 1.0 - 1.0/6.0*rho[i,2,j,k]/concentrations[j]
        wcs.append(wcoft)
        
    # Append these to the array
    wcs_array.append(wcs)
    
# Make them all numpy arrays
densities_array = np.array(densities_array)
cs_array = np.array(cs_array)
wcs_array = np.array(wcs_array)

# Average across the ensemble of simulations
average_wcs = np.sum(wcs_array, axis=0)/float(wcs_array.shape[0])
average_cs = np.sum(cs_array, axis=0)/float(cs_array.shape[0])

# Convert to k_b/atom
average_cs /= (8.617333262e-5/13.6)
cs_array /= (8.617333262e-5/13.6)

# Compute the variance of the SHC
cs_var = np.sum(np.square(cs_array), axis=0)/float(cs_array.shape[0]) - np.square(average_cs)

# Make the plot
fig, ax1 = plt.subplots(figsize=(5,3.75))
ax1.set_xlabel('$T$ (K)')
ax1.set_ylabel('$C$ ($k_B$/atom)')
ax1.plot(temps, average_cs, label='$C$', color='black')
ax1.fill_between(temps, average_cs-np.sqrt(cs_var), average_cs+np.sqrt(cs_var), color='grey', alpha=0.5, label='Error')
ax1.tick_params(axis='y')
ax1.set_ylim(0.0, 10.0)
ax1.set_xlim(0.0, 1200.0)
ax2 = ax1.twinx()
ax2.set_ylabel(r'$\alpha^{pq}_1$')
labels = ['Ni', 'Co', 'Cr']
pairs = [[0,1], [0,2], [1,2]]
for pair in pairs:
    ax2.plot(temps, average_wcs[:,pair[0],pair[1],0], label = labels[pair[1]] + '-' + labels[pair[0]]+' 1st')
ax1.legend(loc='lower right')
ax2.legend(loc='upper right')
ax2.set_ylim(-2.01,1.01)
chandles, clabels = ax1.get_legend_handles_labels()
rhohandles, rholabels = ax2.get_legend_handles_labels()
plt.title('CrCoNi ASRO', fontsize=16)
plt.tight_layout() 
plt.savefig('ASRO.pdf')
