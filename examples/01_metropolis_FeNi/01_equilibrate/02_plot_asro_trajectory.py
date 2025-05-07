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

asro_data = np.loadtxt('trajectories/proc_0000_asro_trajectory_at_T_0300.0.dat')

step_numbers = asro_data[:,0]

# Indices 6 and 10 give me number of Fe-Ni pairs on first and second coordination shell, respectively
asro_1 = 1.0 - 2.0*asro_data[:,6]/12.0
asro_2 = 1.0 - 2.0*asro_data[:,10]/6.0

asro_1 = asro_data[:,6]/12.0
asro_2 = asro_data[:,10]/6.0

plt.plot(step_numbers, asro_1)
plt.plot(step_numbers, asro_2)

plt.xlabel('# of trial Monte Carlo Moves')
plt.ylabel('ASRO (meV/atom)')
plt.title('Equilibration of FeNi')

plt.show()
