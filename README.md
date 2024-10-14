# Bonte Warlo

A code for performing lattice-based atomistic simulations of alloys with an internal energy given by a Bragg-Williams Hamiltonian. The code will be periodically updated so it is best to check the GitHub repository.

Copyright (C) Christopher Woodgate 2019-2024. Released under the GNU Lesser General Public License, version 3.

## Citations
ANY publications/presentations/further work resulting from the use of this software should cite the original publication for which it was developed:
* C. D. Woodgate, J. B. Staunton, Phys. Rev. B **105**, 115124 (2023)
DOI: [10.1103/PhysRevB.105.115124](https://doi.org/10.1103/PhysRevB.105.115124)

In addition, if you use the bcc implementation, you should cite the original paper for which that lattice type was implemented:
* C. D. Woodgate, J. B. Staunton, Phys. Rev. Mater. **7**, 013801 (2023)
DOI: [10.1103/PhysRevMaterials.7.013801](https://doi.org/10.1103/PhysRevMaterials.7.013801)

# Relevant Publications
A full list of publications obtained using this code is:
* C. D. Woodgate, J. B. Staunton, Phys. Rev. B **105**, 115124 (2022). DOI: [https://doi.org/10.1103/PhysRevB.105.115124](https://doi.org/10.1103/PhysRevB.105.115124)
* C. D. Woodgate, J. B. Staunton, Phys. Rev. Mater. **7**, 013801 (2023). DOI: [https://doi.org/10.1103/PhysRevMaterials.7.013801](https://doi.org/10.1103/PhysRevMaterials.7.013801)
* C. D. Woodgate, D. Hedlund, L. H. Lewis, J. B. Staunton, Phys. Rev. Mater. **7**, 053801 (2023). DOI: [https://doi.org/10.1103/PhysRevMaterials.7.053801](https://doi.org/10.1103/PhysRevMaterials.7.053801)
* C. D. Woodgate, J. B. Staunton, J. Appl. Phys. **135**, 135106 (2024). DOI: [https://doi.org/10.1063/5.0200862](https://doi.org/10.1063/5.0200862)
* L. Shenoy, C. D. Woodgate, J. B. Staunton, A. P. Bartók, C. S. Becquart, C. Domain, J. R. Kermode, [Phys. Rev. Mater. **8**, 033804 (2024)](https://doi.org/10.1103/PhysRevMaterials.8.033804).
* C. D. Woodgate, L. H. Lewis, J. B. Staunton, arXiv:2401:02809. DOI: [https://doi.org/10.48550/arXiv.2401.02809](https://doi.org/10.48550/arXiv.2401.02809)
* C. D. Woodgate, G. A. Marchant, L. B. Pártay, J. B. Staunton, [arXiv:2404.13173](https://doi.org/10.48550/arXiv.2404.13173).

## Compilation
At the moment the code is only tested with gfortran and OpenMPI. Put the code in a directory like `~/codes/bontewarlo`. It is my intention to test other compilers in future: watch this space!

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
~/codes/bontewarlo/bontewarlo.run
```

## Running the code
If you navigate to the `examples` subdirectory, you should find two examples demonstrating the code's usage which can be run inside those directories.

Most of the options specified in the input file are fairly self-explanatory. The least obvious is the `mode' option. Because it is my intention to include a 2D (and potentially 1D) option in future, the first digit indicates the number of spatial dimensions for the simulation. Then the last two digits the mode. At present, the implemented options are:
- 01: Simulated Annealing. Uses the Metropolis Monte Carlo algorithm with Kawasaki dynamics to perform simulated annealing on a system in an initially random configuration.
- 02: Draw Decorellated Samples. Optionally performs simulated annealing then draws samples of the grid N Monte Carlo steps apart. Good for generating supercell configurations for use other methods. *E.g.* this recent reference where the code was used to generate training/test data for a machine-learned interatomic potential: L. Shenoy, C. D. Woodgate, J. B. Staunton, A. P. Bartók, C. S. Becquart, C. Domain, J. R. Kermode, [Phys. Rev. Mater. **8**, 033804 (2024)](https://doi.org/10.1103/PhysRevMaterials.8.033804).
- 03: Nested sampling. (Implemented by Livia Bartók-Pártay, Department of Chemistry, University of Warwick, UK.) Uses the nested sampling algorithm to sample the configuration space from random initial configurations, facilitating calculation of the partition function at an arbitrary temperature during the post-processing step. This procedure is outlined in a recent preprint: C. D. Woodgate, G. A. Marchant, L. B. Pártay, J. B. Staunton, [arXiv:2404.13173](https://doi.org/10.48550/arXiv.2404.13173).
- 04: Transition Matrix Monte Carlo Method. (Implemented by Hubert Nagusewski, Department of Physics, University of Warwick.) 
- 05: Wang-Landau Sampling. (Implemented by Hubert Nagusewski, Department of Physics, University of Warwick.) 
- 06: Energy Spetrum algorithm. (Implemented by Hubert Nagusewski, Department of Physics, University of Warwick.) Needed for modes 04 and 05 to determine the range of possible simulation energies.

## Author
Christopher D. Woodgate

## Contributors
- Livia Bartók-Pártay (Nested Sampling implementation).
- Hubert Naguszewski (Wang-Landau sampling and Transition Matrix Monte Carlo (TMMC) implementations).

## Contributing
Any/all contributions are welcome via pull requests. 

## License
This software is released under the GPL license. See the file LICENSE.txt for details.

## Questions, Brickbats, and Bouquets
[Drop me an email](mailto:christopher.woodgate@physics.org) if you have any questions, comments, or feedback!
