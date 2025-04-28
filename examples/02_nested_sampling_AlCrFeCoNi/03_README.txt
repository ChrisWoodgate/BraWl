Example 3: Nested Sampling of the AlCrFeCoNi system, in the fcc structure.

This example demonstrates a nested sampling simulation of the AlCrFeCoNi system, in the fcc structure. As a quick demonstration, the input file is set up to use 108 atoms, and only a low-resolution sampling using 100 walkers. 
Atom-atom interactions are taken from:
C. D. Woodgate, G. A. Marchant, L. B. Pártay, J. B. Staunton, npj Comput. Mater. 10, 271 (2024)
DOI: https://doi.org/10.1038/s41524-024-01445-w 

A nested sampling run needs the following input files:
 - The file 'inputbrawl.txt' contains all relevant parameters for the simulation, while 'NiCoCr_V_ijs.txt' contains the atom-atom interaction parameters.
 - The file 'ns_input.txt' contains the relevant nested sampling parameters.

Run './brawl.run control=input_brawl.txt ns_input.txt' to launch the simulation. Once it has finished, analyse the results with 01_analyse.py

During the run the following information is printed:
The initial energies of the randomly generated configurations (these are the initial walkers)
Once the main nested sampling cycle starts, the following values are printed at regular intervals:
NS iteration number, current energy, current MC swap acceptance rate, total number of accepted MC swaps in the current iteration

** post-processing **

The '.energies' file can be used for post-processing calculation of the partition function. The format of the file such that it enables the 
usage of the ns_analyse python code from the pymatnest package.
(pymatnest is an open source python package for materials` nested sampling calculations, available on github:
https://github.com/libAtoms/pymatnest)

For the above test use ns_analyse as the following (-M minimum temperature in K, -D temperature increment in K, -n number of temperatures to evaluate thermodynamic properties at, -k Boltzmann-constant in Ry/K units)
./ns_analyse fcc_al_1.00_crfeconi_K100.energies -M 10 -D 5 -n 800 -k 6.3336374823584e-6 > HC_fcc_al_1.00_crfeconi_K100.energies.txt

The output of ns_analyse (HC_fcc_al_1.00_crfeconi_K100.txt) contains in columns, the temperature, log(partition function), total free energy, potential energy, heat capacity...etc. values as a function of temperature.
