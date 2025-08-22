## Simulating the thermodynamics and phase stability of multicomponent alloys using conventional and enhanced sampling techniques

The Fortran package **BraWl** (named after the Bragg-Williams model)  facilitates simulation of the thermodynamics and phase stability of both binary and multicomponent alloys. It achieves this by providing implementation of both the Bragg-Williams Hamiltonian (a lattice based model expressing the internal energy of an alloy as a sum of atom-atom effective pair interactions) concurrently with a range of conventional and enhanced sampling techniques for exploration of the alloy configuration space. The result is a package which can determine phase equilibria as a function of both temperature and alloy composition, which leads to the construction of alloy phase diagrams. Additionally, the package can be used for extraction of representative equilibrated atomic configurations for visualisation, as well as for use in complementary modelling approaches. It provides a lightweight, fast and flexible foundation for a range of simulations relating to alloy thermodynamics and phase diagrams.

For an overview of the capabilities of the package, you can check out our recent preprint: H. J. Naguszewski, L. B. Partay, D. Quigley, C. D. Woodgate, [arXiv:2505.05393](https://doi.org/10.48550/arXiv.2505.05393).

The latest version of this documentation can be found at: [chriswoodgate.github.io/BraWl/](https://chriswoodgate.github.io/BraWl/)

The GitHub repo where the source code can be found is: [github.com/ChrisWoodgate/BraWl](https://github.com/ChrisWoodgate/BraWl)

---

## 🚀 Features

- 💡 Clean Fortran90 module structure
- 📦 Organized components for reusability
- ⚙️  MPI support for parallelism
- 🔬 Designed for scientific extensibility

---

## 📝 Citation

If you use **BraWl** in your research, please cite our preprint documenting the package and its capabilities:
* H. J. Naguszewski, L. B. Partay, D. Quigley, C. D. Woodgate, [arXiv:2505.05393](https://doi.org/10.48550/arXiv.2505.05393).

---

## 📚 Publications

A (hopefully fairly complete) list of publications obtained using the package is as follows:

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

---

## 💻 Authors

- Hubert J. Naguszewski
- Livia B. Pártay
- Christopher D. Woodgate

---

## 💻 Contributors

- Heather Ratcliffe
- David Quigley

---

## Contributing

Contributions are welcome via pull requests. 

Features which are particularly welcome are:
- Implementations of new lattice types.
- Implementations of new/interesting sampling algorithms.
- Interfaces with other packages, _e.g._ writing output configurations to formats specific to other materials simulation codes.

Features which we do not currently plan on implementing are:
- Anything off-lattice. (The way configurations are represented internally would make this challenging!)

Additionally, as [developer documentation for the package](https://chriswoodgate.github.io/BraWl/) is auto-generated using [Doxygen](https://www.doxygen.nl), we ask that all contributions to the codebase follow the [Doxygen guide for formatting comment blocks](https://www.doxygen.nl/manual/docblocks.html#fortranblocks) and uses the relevant [Doxygen commands](https://www.doxygen.nl/manual/commands.html) to specify the details of new features/routines.


---

## 🪪 License

This software is released under the LGPL-3.0 license. See the file LICENSE.txt for details.

---

## 🧭 Documentation Guide

- **Modules** — High-level structure of the project
    - **Modules List** – List of modules
    - **Module Members** – Complete list of members of *all* modules
- **Data Types** — Derived types/classes used in the project
- **Files** — Browse the project's file structure

---

## 📁 Project Structure

<pre>
BraWl/
├── src/               # Core source code
├── include/           # Any header files
├── examples/          # Examples of use
├── tests/             # Test cases
├── docs/              # Documentation (this!)
├── README.md          # High-level overview and README
├── Makefile           # Makefile (for building the code!)
├── LICENSE.txt        # Copyright statement and license (LGPL-3.0)
└── CODE_OF_CONDUCT.md # Contributor Covenant Code of Conduct
</pre>

