# 🔩 BraWl: Simulating the thermodynamics and phase stability of multicomponent alloys using conventional and enhanced sampling techniques

The Fortran package **BraWl** (named after the Bragg-Williams model)  facilitates simulation of the thermodynamics and phase stability of both binary and multicomponent alloys. It achieves this by providing implementation of both the Bragg-Williams Hamiltonian (a lattice based model expressing the internal energy of an alloy as a sum of atom-atom effective pair interactions) concurrently with a range of conventional and enhanced sampling techniques for exploration of the alloy configuration space. The result is a package which can determine phase equilibria as a function of both temperature and alloy composition, which leads to the construction of alloy phase diagrams. Additionally, the package can be used for extraction of representative equilibrated atomic configurations for visualisation, as well as for use in complementary modelling approaches.

It provides a lightweight, fast and flexible foundation for a range of simulations relating to alloy thermodynamics and phase diagrams.

For an overview of the capabilities of the package, you can check out our recent preprint: H. J. Naguszewski, L. B. Partay, D. Quigley, C. D. Woodgate, [arXiv:2505.05393](https://doi.org/10.48550/arXiv.2505.05393).

---

## 🚀 Features

- 💡 Clean Fortran90 module structure
- 📦 Organized components for reusability
- ⚙️  MPI support for parallelism
- 🔬 Designed for scientific extensibility

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
