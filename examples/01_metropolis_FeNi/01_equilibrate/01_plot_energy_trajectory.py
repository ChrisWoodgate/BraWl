'''
Example script for plotting simulation energy and ASRO as function of step number

C. D. Woodgate
'''
import numpy as np
import matplotlib.pyplot as plt

# CODATA value for conversion of eV to Ry
ry_in_ev = 13.605693122990

# Names of elements
elements = ['Fe', 'Ni']

# Concentrations of species
concentrations = [0.5, 0.5]

# File containing list of step numbers and associated energies
energy_data = np.loadtxt('trajectories/proc_0000_energy_trajectory_at_T_0300.0.dat')

# Pull out the relevant data
step_numbers = energy_data[:,0]
energies = energy_data[:,1] - energy_data[0,1]

# Convert total simulation energy (in Ry) to meV/atom
energies = energies/256.0*1000.0*ry_in_ev

# One 'sweep' is one trial move per atom
sweeps = step_numbers/256

# File containing list of step numbers and associated conditional pair probabilities
asro_data = np.loadtxt('trajectories/proc_0000_asro_trajectory_at_T_0300.0.dat')

# Indices 6 and 10 give me number of Fe-Ni pairs on first and second coordination shell, respectively
# Divide by the number of atoms on each of those coordination shells
asro_1 = asro_data[:,6]/12.0
asro_2 = asro_data[:,10]/6.0

# Plot the data
fig, ax = plt.subplots(2, 1, figsize=(4,5), sharex=True)

ax[0].plot(sweeps, energies, color='black')
ax[0].set_ylabel('Simulation energy (meV/atom)')
ax[0].ticklabel_format(axis='x', style='sci', scilimits=(10, 1e9))

ax[1].plot(sweeps, asro_1, label='Fe-Ni 1st nearest-neighbours')
ax[1].plot(sweeps, asro_2, label='Fe-Ni 2nd nearest-neighbours')
ax[1].set_xlabel('# of Monte Carlo sweeps')
ax[1].set_ylabel('Pair Probabilities')
ax[1].legend(fontsize='small')

fig.suptitle('Equilibration of FeNi at $T=300$ K')
plt.tight_layout()
plt.savefig('FeNi_energy_equilibration_300K.pdf')
