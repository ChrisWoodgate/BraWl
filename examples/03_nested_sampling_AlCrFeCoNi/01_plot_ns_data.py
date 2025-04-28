'''
Example script for plotting results of an NS simulation

C. D. Woodgate, L. B. Partay
'''
import numpy as np
import matplotlib.pyplot as plt

# CODATA value for conversion of eV to Ry
ry_in_ev = 13.605693122990

# Load the NS data
ns_data = np.loadtxt('HC_fcc_al_1.00_crfeconi_K100.txt')

# Pull out the temperature, energy, and Cv data
temperature = ns_data[:,0]
energy = ns_data[:,3]
heat_capacity = ns_data[:,4]

# Rescale the energy so that the highest energy is set as the 'zero' of the energy scale
energy = energy - energy[-1]
energy = energy/108.0*1000.0*ry_in_ev

# TODO Need to check unit conversion on Cv. Would like k_B/atom
heat_capacity = heat_capacity * 1.0 # Replace 1.0 with correct conversion factor

# Make the plots
fig, ax = plt.subplots(figsize=(5,3.25))

# Title for the whole figure
fig.suptitle('Nested Sampling on fcc AlCrFeCoNi')

# Set the x limits
ax.set_xlim(0.0, temperature[-1])

# Plot the energy data on one axis
ax.plot(temperature, energy, label='$E$', color='black', linestyle='dashed')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Simulation energy, $E$ (meV/atom)')

# And the heat capacity data on another
ax2 = ax.twinx()
ax2.plot(temperature, heat_capacity, label='$C_V$')
ax2.set_ylabel('Heat Capacity, $C_V$ ($k_B$/atom)')
ax2.set_ylim(0.0, 1.1*np.max(heat_capacity))

fig.legend(loc=(0.65, 0.45))

plt.tight_layout()
plt.savefig('ns_alcrfeconi.pdf')
