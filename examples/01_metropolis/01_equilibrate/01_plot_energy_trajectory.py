'''
Example script for plotting simulation energy and ASRO as function of step number

C. D. Woodgate
'''
import numpy as np
import matplotlib.pyplot as plt

# Names of elements
elements = ['Fe', 'Ni']

# Concentrations of species
concentrations = [0.5, 0.5]

energy_data = np.loadtxt('trajectories/proc_0000_energy_trajectory_at_T_0300.0.dat')

step_numbers = energy_data[:,0]

energies = energy_data[:,1]

plt.plot(step_numbers, energies)

plt.xlabel('# of trial Monte Carlo Moves')
plt.ylabel('Simulation energy (meV/atom)')
plt.title('Equilibration of FeNi')

plt.show()
