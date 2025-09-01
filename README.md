# BraWl

[![Ubuntu (x86_64)](https://github.com/ChrisWoodgate/BraWl/actions/workflows/testingOnPush_Ubuntu.yaml/badge.svg)](https://github.com/ChrisWoodgate/BraWl/actions/workflows/testingOnPush_Ubuntu.yaml)
[![MacOS (ARM64)](https://github.com/ChrisWoodgate/BraWl/actions/workflows/testingOnPush_Apple.yaml/badge.svg)](https://github.com/ChrisWoodgate/BraWl/actions/workflows/testingOnPush_Apple.yaml)
[![License: MIT](https://img.shields.io/badge/License-LGPLv3-yellow.svg)](https://opensource.org/license/lgpl-3-0)

`BraWl` — a package for performing lattice-based atomistic simulations of alloys with an internal energy given by a Bragg-Williams Hamiltonian. The package implements a range of conventional and enhanced sampling techniques, including Metropolis-Hastings Monte Carlo, Nested Sampling, and Wang-Landau sampling. This `README` contains key information and installation instructions, while developer documentation can be found at: [https://chriswoodgate.github.io/BraWl/](https://chriswoodgate.github.io/BraWl/)

Copyright (C) The Authors 2019-2025. Released under the GNU Lesser General Public License, version 3.

## Background

The Bragg-Williams Hamiltonian is an on-lattice Ising-like Hamiltonian describing the internal energy of a general substitutional alloy. The configuration of the alloy is specified by the *site occupation numbers*, $\{\xi_{i\alpha}\}$, where $\xi_{i\alpha}=1$ if site $i$ is occupied by an atom of species $\alpha$, and $\xi_{i\alpha}=0$ otherwise. (These can be thought of as a little like the 'spins' of the Ising model.) Each lattice site must be constrained to have one (and only one) atom sitting on it, expressed as

$$\sum_\alpha \xi_{i\alpha}=1$$

for all lattice sites $i$. The overall concentration of a chemical species, $c_\alpha$ is given by

$$c_\alpha = \frac{1}{N} \sum_i \xi_{i\alpha},$$

where $N$ is the total number of lattice sites in the system. The energy associated with an atom of species $\alpha$ on site $i$ interacting with an atom of species $\alpha'$ on site $j$, referred to as an atom-atom *effective pair interaction* is denoted $V_{i\alpha; j\alpha'}$. The Bragg-Williams Hamiltonian is then written

$$H(\{\xi_{i\alpha}\}) = \frac{1}{2}\sum_{i \alpha; j\alpha'} V_{i\alpha; j\alpha'} \xi_{i \alpha} \xi_{j \alpha'},$$

where the factor of $\frac{1}{2}$ accounts for double-counting.

In practice, it is assumed that interactions are isotropic (direction-independent) and homogeneous (translationally invariant). Moreover, it is assumed that interactions have some finite range (a 'cutoff' radius), beyond which they are set to zero. The Hamiltonian is therefore re-written as

$$H(\{\xi_{i\alpha}\}) = \frac{1}{2}\sum_{i} \sum_{n} \left(\sum_{j \in n(i)} \sum_{\alpha \alpha'} V_{\alpha \alpha'}^{(n)} \xi_{i \alpha} \xi_{j \alpha'} \right),$$

where $n(i)$ denotes the set of lattice sites which are $n$ th nearest-neighbours to site $i$, and $V_{\alpha \alpha'}^{(n)}$ denotes the interaction between species $\alpha$ and $\alpha'$ on coordination shell $n$. This is the description which with BraWL works internally.

## Dependencies

### Main package
To build and run the main package with full functionality, you will require:
- A Fortran compiler.
- A Message Passing Interface (MPI) implementation.
- A LAPACK/BLAS implementation, such as `OpenBLAS`.
- The NetCDF-Fortran library. (This is _usually_ dependent on the standard NetCDF library, so you may need to explicitly install this, too.)
- GNU Make, for building the package.

Optionally, to build a serial version of the code with more limited functionality (Wang-Landau sampling is currently _only_ a parallel implementation due to the nature of the algorithm), it is sufficient simply to have GNU Make, a Fortran compiler and the NetCDF-Fortran library.

The code has been explicitly tested and verified to function with the following combination of compiler and library versions:
```
GCC/13.2.0
OpenMPI/4.1.6
netCDF-Fortran/4.6.1
```
From experience, compatibility with earlier/later versions of these libraries is strongly expected, but not explicitly guaranteed. Additionally, we have also had success compiling the code using the Intel Fortran compiler and Intel's MPI implementation.

### Nested Sampling analysis
In order to _analyse_ the output of a Nested Sampling simulation, you will require the `pymatnest` package, for which a list of dependencies and installation instructions can be found at the GitHub repository: [https://github.com/libAtoms/pymatnest](https://github.com/libAtoms/pymatnest).

### Example plotting scripts
For the example plotting scripts to function the following dependencies are required:
```
Python/3.11.5
Numpy/2.2.5
Matplotlib/3.10.1
Cycler/0.12.1
netCDF/1.7.2
```
The Python modules can be installed using `pip install -r requirements.txt`, executed in the main directory. It is recommended to install these within a virtual environment which can be created within `BraWl` using `python -m venv venv` and the activated using `source venv/bin/activate`. (Compatibility of the plotting scripts with later/earlier versions of the listed packages is anticipated, but not guaranteed.)

## Compilation
At the moment the code has been tested with GCC and OpenMPI, versions as specified above. Put the code in a directory like `~/codes/BraWl` and navigate to that directory. Assuming the dependencies mentioned above have been installed and configured correctly, building the code should be as simple as running the commands
```
make compiler=mpifort
```
if a parallel compilation is desired, or
```
make compiler=gfortran
```
if a serial compilation is desired.

### System-specific instructions
On the Warwick SCRTP system, which uses environment modules, the relevant modules (for use on the `taskfarm` and/or interactive SCRTP machines) can be loaded using the commands
```
module purge
module load GCC/13.2.0 OpenMPI/4.1.6 netCDF-Fortran/4.6.1
```

On Bristol's `bluepebble` cluster, the relevant modules (as of 2025/03) can be loaded using
```
module add gcc/12.3.0 openmpi/5.0.3 netcdf-c/4.9.2 netcdf-fortran
```

## Running the code
The code can be run from any directory containing the relevant input files, simply by running
```
mpirun -np <num_processors> /path/to/BraWl/brawl.run
```
if built in parallel, or
```
/path/to/BraWl/brawl.run
```
in serial.

By default, the code will look for input files named `brawl.inp` (defining the model), an interaction file specified in `brawl.inp`, and then an additional input file relevant to the sampling algorithm being used. For Metropolis-Hastings, Wang-Landau sampling, and Nested Sampling, the default behaviour is to look for files named `metropolis.inp`, `wang-landau.inp`, and `nested_sampling.inp` for these, respectively. If you wish to use files with names different to these, the code accepts command-line arguments which allow the user to specify the names of these files, _e.g._
```
mpirun -np <num_processors> /path/to/BraWl/brawl.run input=<brawl_input_name> metropolis=<metropolis_input_name>
```

## Examples
If you navigate to the `examples` subdirectory, you should find some example input files demonstrating the code's usage which can be run inside those directories. These input files are commented to explain what the various parameters mean and do.

Most options specified in the input files are fairly self-explanatory. Commented examples of input files can be found in the `examples` subdirectory. The least obvious is the `mode' option of `brawl.inp`. Because it is our intention to include 2D (and potentially 1D) options in future, the first digit indicates the number of spatial dimensions for the simulation. Then the last two digits the mode. At present, the implemented (and fully tested) options are:
- 01: Metropolis-Hastings Monte Carlo. Uses the Metropolis-Hastings algortithm to equilibrate a system then perform sampling. Can also be used to perform simulated annealing, _e.g._ as used in [npj Comput. Mater. **10** 272 (2024)](https://doi.org/10.1038/s41524-024-01435-y), or to draw decorrelated samples for use in other modelling approaches, such as training machine-learned interatomic potentials, _e.g._ as used in [Phys. Rev. Mater. **8**, 033804 (2024)](https://doi.org/10.1103/PhysRevMaterials.8.033804).
- 02: Wang-Landau sampling. Uses Wang-Landau sampling to compute the simulation density of states in energy, which can then be used to sample thermodynamic quantities as a post-processing step. This procedure is outlined in [arXiv:2503.13235](https://arxiv.org/abs/2503.13235)
- 03: Nested sampling. Uses the nested sampling algorithm to sample the configuration space from random initial configurations, allowing to calculate the partition function at an arbitrary temperature during the post-processing step. This procedure is outlined in [npj Comput. Mater. **10**, 271 (2024)](https://doi.org/10.1038/s41524-024-01445-w).

## Tests
The `tests` directory contains some test inputs/outputs to verify the code's core functionality once compiled. These tests are run automatically using GitHub Actions. You can find full details in `.github/workflows/`

## Documentation
The in addition to this README and the provided examples, the code also has (searchable) documentation which is auto-generated using [Doxygen](https://www.doxygen.nl). This documentation contains information about all modules, functions, subroutines, and derived types, so is particularly useful if you are looking to develop a new feature. There are two options available for viewing this documentation:
1. View the auto-generated documentation hosted online at [https://chriswoodgate.github.io/BraWl/](https://chriswoodgate.github.io/BraWl/). (Intended for end-users and developers looking to implement new features.)
2. Generate and 'host' your own version of the documentation on your local machine. (Intended for developers wanting to check that their comments/documentation is being rendered as expected.)

### Local generation and 'hosting' of documentation
To generate and host your own version of the documentation, you should:
1. Obtain [Doxygen](https://www.doxygen.nl), which is typically available through a package manager.
2. Run `doxygen docs/doxyfile` from the code's main directory.
3. Navigate to the freshly-generated `docs/html` directory and `open` index.html. (Alternatively, use your system's filemanager to open this in your web browser.)

## Citations
If you use `BraWl` in your research, please cite our preprint documenting the package and its capabilities:
* H. J. Naguszewski, L. B. Partay, D. Quigley, C. D. Woodgate, [arXiv:2505.05393](https://doi.org/10.48550/arXiv.2505.05393).

## List of publications
A (hopefully fairly complete) list of publications obtained using this code is:
1. G. A. Marchant, C. D. Woodgate, C. E. Patrick, J. B. Staunton, [Phys. Rev. B **103**, 094414 (2021)](https://doi.org/10.1103/PhysRevB.103.094414).
2. C. D. Woodgate, J. B. Staunton, [Phys. Rev. B **105**, 115124 (2022)](https://doi.org/10.1103/PhysRevB.105.115124).
3. C. D. Woodgate, J. B. Staunton, [Phys. Rev. Mater. **7**, 013801 (2023)](https://doi.org/10.1103/PhysRevMaterials.7.013801).
4. C. D. Woodgate, D. Hedlund, L. H. Lewis, J. B. Staunton, [Phys. Rev. Mater. **7**, 053801 (2023)](https://doi.org/10.1103/PhysRevMaterials.7.053801).
5. C. D. Woodgate, J. B. Staunton, [J. Appl. Phys. **135**, 135106 (2024)](https://doi.org/10.1063/5.0200862).
6. L. Shenoy, C. D. Woodgate, J. B. Staunton, A. P. Bartók, C. S. Becquart, C. Domain, J. R. Kermode, [Phys. Rev. Mater. **8**, 033804 (2024)](https://doi.org/10.1103/PhysRevMaterials.8.033804).
7. C. D. Woodgate, _"Modelling Atomic Arrangements in Multicomponent Alloys: A Perturbative, First-Principles-Based Approach"_ [Springer Series in Materials Science, Vol. 346 (Springer Nature Switzerland, Cham, 2024)](https://doi.org/10.1007/978-3-031-62021-8).
8. C. D. Woodgate, G. A. Marchant, L. B. Pártay, J. B. Staunton, [npj Comput. Mater. **10**, 271 (2024)](https://doi.org/10.1038/s41524-024-01445-w).
9. C. D. Woodgate, L. H. Lewis, J. B. Staunton, [npj Comput. Mater. **10**, 272 (2024)](https://doi.org/10.1038/s41524-024-01435-y).
10. C. D. Woodgate, H. J. Naguszewski, D. Redka, J. Minar, D. Quigley, J. B. Staunton, [J. Phys.: Mater. **8**, 045002 (2025)](https://doi.org/10.1088/2515-7639/adf468).

## Authors
- Hubert J. Naguszewski
- Livia B. Pártay
- Christopher D. Woodgate

## Contributors
- Heather Ratcliffe
- David Quigley
- Adam M. Krajewski

## Contributing
Contributions are welcome via pull requests. 

Features which are particularly welcome are:
- Implementations of new lattice types.
- Implementations of new/interesting sampling algorithms.
- Interfaces with other packages, _e.g._ writing output configurations to formats specific to other materials simulation codes.

Features which we do not currently plan on implementing are:
- Anything off-lattice. (The way configurations are represented would make this challenging!)

Additionally, as [developer documentation for the package](https://chriswoodgate.github.io/BraWl/) is auto-generated using [Doxygen](https://www.doxygen.nl), we ask that all contributions to the codebase follow the [Doxygen guide for formatting comment blocks](https://www.doxygen.nl/manual/docblocks.html#fortranblocks) and uses the relevant [Doxygen commands](https://www.doxygen.nl/manual/commands.html) to specify the details of new features/routines.

## License
This software is released under the LGPL-3.0 license. See the file LICENSE.txt for details.

## Questions, Brickbats, and Bouquets
[Drop me an email](mailto:christopher.woodgate@physics.org) if you have any questions, comments, or feedback!
