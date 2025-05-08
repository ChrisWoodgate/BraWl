'''
Example script for plotting heat capacity and ASRO as function of temperature (Metropolis-Hastings simulated annealing)

C. D. Woodgate
'''
import numpy as np
import matplotlib.pyplot as plt
from  netCDF4 import Dataset

# CODATA values for conversion of eV to Ry, etc.
ry_in_ev = 13.605693122990
kb_in_ev = 8.617333262e-5

# Names of elements
elements = ['Fe', 'Ni']

# Concentrations of species
concentrations = [0.5, 0.5]

# Pull out data relating to temperature and energy
energy_data = np.loadtxt('energies/av_energy_diagnostics.dat')
temps = energy_data[:,0]
C = energy_data[:,2]/(kb_in_ev/ry_in_ev)

# Pull out data relating to atomic short-range order (ASRO)
nc_data = Dataset('asro/av_radial_density.nc', "r", format="NETCDF4")
rho = nc_data.variables['rho data'][:]

# Turn the radial densities into Warren-Cowley parameters
wcs = []
for i in range(len(temps)):
    wcoft = np.zeros((3,3,2))
    for j in range(len(elements)):
        for k in range(len(elements)):
            wcoft[j,k,0] = 1.0-1.0/12.0*rho[i,1,j,k]/concentrations[j]
            wcoft[j,k,1] = 1.0-1.0/6.0*rho[i,2,j,k]/concentrations[j]
    wcs.append(wcoft)
wcs = np.array(wcs)

# Close the NetCDF file
nc_data.close()

# Make the plot
fig, ax1 = plt.subplots(figsize=(6,4))

# Use the first axis for the heat capacity data
ax1.set_xlabel(r'Temperature, $T$ (K)')
ax1.set_ylabel(r'Heat capacity, $C_V$ ($k_B$/atom)')
ax1.plot(temps, C, label=r'$C_V$', color='black')
ax1.set_ylim(0.0, 2.5)
ax1.set_xlim(0.0, 1200.0)
ax1.legend(loc='upper right')

labels = ['Fe', 'Ni']

# And the second for  the ASRO
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set_ylabel(r'Warren-Cowley ASRO parameters, $\alpha^{pq}_n$')  # we already handled the x-label with ax1
ax2.plot(temps, wcs[:,0,1,0], label = labels[0] + '-' + labels[1] + r' $n=1$')
ax2.plot(temps, wcs[:,0,1,1], label = labels[0] + '-' + labels[1] + r' $n=2$')
ax2.legend(loc='lower right')
ax2.set_ylim(-1.01, 1.01) # Make the bounds slightly large so we can see the data more clearly

plt.title('FeNi Simulated Annealing')
fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.savefig('FeNi_simulated_annealing_results.pdf')
