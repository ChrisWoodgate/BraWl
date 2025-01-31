# BraWl

A code for performing lattice-based atomistic simulations of alloys with an internal energy given by a Bragg-Williams Hamiltonian, implementing both Monte Carlo methods (simulated annealing, transition matrix Monte Carlo, Wang-Landau sampling) and nested sampling. The code will be periodically updated so it is best to check the GitHub repository for the latest version.

Copyright (C) The Authors 2019-2025. Released under the GNU Lesser General Public License, version 3.

## Background

The Bragg-Williams Hamiltonian is an on-lattice Ising-like Hamiltonian describing the internal energy of a general substitutional alloy. The configuration of the alloy is specified by the *site occupation numbers*, $\{\xi_{i\alpha}\}$, where $\xi_{i\alpha}=1$ if site $i$ is occupied by an atom of species $\alpha$, and $\xi_{i\alpha}=0$ otherwise. Each lattice site must be constrained to have one (and only one) atom sitting on it, expressed as $\sum_\alpha \xi_{i\alpha}=1$ for all lattice sites $i$. The overall concentration of a chemical species, $c_\alpha$ is given by $c_\alpha = \frac{1}{N} \sum_i \xi_{i\alpha}$, where $N$ is the total number of lattice sites in the system. The energy associated with an atom of species $\alpha$ on site $i$ interacting with an atom of species $\alpha'$ on site $j$, referred to as an *effective pair interaction* is denoted $V_{i\alpha; j\alpha'}$. The Bragg-Williams Hamiltonian is then written
$$H(\{\xi_{i\alpha}\}) = \frac{1}{2}\sum_{i \alpha; j\alpha'} V_{i\alpha; j\alpha'} \xi_{i \alpha} \xi_{j \alpha'},$$
where the factor of $\frac{1}{2}$ accounts for double-counting.

## Citations
ANY publications/presentations/further work resulting from the use of this software should cite the original publication for which it was developed:
* C. D. Woodgate, J. B. Staunton, [Phys. Rev. B **105**, 115124 (2023)](https://doi.org/10.1103/PhysRevB.105.115124)

# Relevant Publications
A list of publications obtained using this code is:
1. G. A. Marchant, C. D. Woodgate, C. E. Patrick, J. B. Staunton, [Phys. Rev. B **103**, 094414 (2021)](https://doi.org/10.1103/PhysRevB.103.094414).
2. C. D. Woodgate, J. B. Staunton, [Phys. Rev. B **105**, 115124 (2022)](https://doi.org/10.1103/PhysRevB.105.115124).
3. C. D. Woodgate, J. B. Staunton, [Phys. Rev. Mater. **7**, 013801 (2023)](https://doi.org/10.1103/PhysRevMaterials.7.013801).
4. C. D. Woodgate, D. Hedlund, L. H. Lewis, J. B. Staunton, [Phys. Rev. Mater. **7**, 053801 (2023)](https://doi.org/10.1103/PhysRevMaterials.7.053801).
5. C. D. Woodgate, J. B. Staunton, [J. Appl. Phys. **135**, 135106 (2024)](https://doi.org/10.1063/5.0200862).
6. L. Shenoy, C. D. Woodgate, J. B. Staunton, A. P. Bartók, C. S. Becquart, C. Domain, J. R. Kermode, [Phys. Rev. Mater. **8**, 033804 (2024)](https://doi.org/10.1103/PhysRevMaterials.8.033804).
7. C. D. Woodgate, G. A. Marchant, L. B. Pártay, J. B. Staunton, [npj Comput. Mater. **10**, 271 (2024)](https://doi.org/10.1038/s41524-024-01445-w).
8. C. D. Woodgate, L. H. Lewis, J. B. Staunton, [npj Comput. Mater. **10**, 272 (2024)](https://doi.org/10.1038/s41524-024-01435-y).

## Compilation
At the moment the code is only tested with gfortran and OpenMPI. Put the code in a directory like `~/codes/BraWl`. It is my intention to test other compilers in future: watch this space!

On the Warwick SCRTP system, which uses environment modules, run
```
module purge
module load GCC/11.3.0 OpenMPI/4.1.4 netCDF-Fortran/4.6.0
```
then you should be able to build the code with
```
make compiler=gfortran
```
and run it in a directory with a suitable input file via
```
~/codes/BraWl/brawl.run
```
(Note that the code, by default, will look for a file called `input.txt` as its input.)

## Running the code
If you navigate to the `examples` subdirectory, you should find two examples demonstrating the code's usage which can be run inside those directories.

Most of the options specified in the input file are fairly self-explanatory. The least obvious is the `mode' option. Because it is my intention to include a 2D (and potentially 1D) option in future, the first digit indicates the number of spatial dimensions for the simulation. Then the last two digits the mode. At present, the implemented options are:
- 01: Simulated Annealing. Uses the Metropolis Monte Carlo algorithm with Kawasaki dynamics to perform simulated annealing on a system in an initially random configuration.
- 02: Draw Decorellated Samples. Optionally performs simulated annealing then draws samples of the grid N Monte Carlo steps apart. Good for generating supercell configurations for use other methods. *E.g.* this recent reference where the code was used to generate training/test data for a machine-learned interatomic potential: L. Shenoy, C. D. Woodgate, J. B. Staunton, A. P. Bartók, C. S. Becquart, C. Domain, J. R. Kermode, [Phys. Rev. Mater. **8**, 033804 (2024)](https://doi.org/10.1103/PhysRevMaterials.8.033804).
- 03: Nested sampling. Uses the nested sampling algorithm to sample the configuration space from random initial configurations, allowing to calculate the partition function at an arbitrary temperature during the post-processing step. This procedure is outlined in a recent publication: C. D. Woodgate, G. A. Marchant, L. B. Pártay, J. B. Staunton, [npj Comput. Mater. **10**, 271 (2024)](https://doi.org/10.1038/s41524-024-01445-w).

## Authors
- Hubert J. Naguszewski
- Livia B. Pártay
- Christopher D. Woodgate

## Contributing
Any/all contributions are welcome via pull requests. 

## License
This software is released under the LGPL-3.0 license. See the file LICENSE.txt for details.

## Questions, Brickbats, and Bouquets
[Drop me an email](mailto:christopher.woodgate@physics.org) if you have any questions, comments, or feedback!
