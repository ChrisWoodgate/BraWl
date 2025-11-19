# BraWl: Simulating the thermodynamics and phase stability of multicomponent alloys using conventional and enhanced sampling techniques

The Fortran package **BraWl** (named after the Bragg-Williams model)  facilitates simulation of the thermodynamics and phase stability of both binary and multicomponent alloys. It achieves this by providing implementation of both the Bragg-Williams Hamiltonian (a lattice based model expressing the internal energy of an alloy as a sum of atom-atom effective pair interactions) concurrently with a range of conventional and enhanced sampling techniques for exploration of the alloy configuration space. The result is a package which can determine phase equilibria as a function of both temperature and alloy composition, which leads to the construction of alloy phase diagrams. Additionally, the package can be used for extraction of representative equilibrated atomic configurations for visualisation, as well as for use in complementary modelling approaches. It provides a lightweight, fast and flexible foundation for a range of simulations relating to alloy thermodynamics and phase diagrams.

For an overview of the theoretical backgroun and capabilities of the package, you can check out our recent preprint: H. J. Naguszewski, L. B. Partay, D. Quigley, C. D. Woodgate, [arXiv:2505.05393](https://doi.org/10.48550/arXiv.2505.05393). You can also checkout the at the following subpage(s):

* @ref theory


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

@page theory Theoretical Background

Here we briefly review some elementary details of alloy thermodynamics and statistical mechanics. You can also find much of this information on the arXiv: [arXiv:2505.05393](https://doi.org/10.48550/arXiv.2505.05393).

## Alloy Thermodynamics

In a substitutional alloy with fixed underlying lattice, a particular arrangement of atoms can be specified by a discrete set of site occupation numbers, \f${ \xi_{i\gamma} }\f$, where \f$\xi_{i \gamma} = 1\f$ if site \f$i\f$ is occupied by an atom of species \f$\gamma\f$ and \f$\xi_{i \gamma} = 0\f$ otherwise. The lattice index \f$i\f$ takes values in the range 1 to \f$N\f$, where \f$N\f$ is the total number of lattice sites in the system, while the species index \f$\gamma\f$ takes values in the range 1 to \f$s\f$, where \f$s\f$ is the number of chemical species (elements) present in the alloy composition. The physical constraint that each lattice site is occupied by one (and only one)  atom is expressed as \f$\sum_{\gamma} \xi_{i\gamma} = 1\f$, while the total concentration of a given chemical species is given by \f$c_\gamma = \frac{1}{N} \sum_i \xi_{i\gamma}\f$. (Naturally, vacancies can be treated in this framework by considering them as an additional chemical species present at a very low concentration.) If we consider an ensemble of alloy configurations, we can then define the site-wise concentrations, \f$\{ c_{i\gamma}\}\f$, as the average value of the site occupation numbers across the ensemble, \f$c_{i\gamma} = \langle \xi_{i\gamma} \rangle\f$, where \f$\langle \cdot \rangle\f$ denotes the average taken over the ensemble. (The precise meaning of `ensemble' will be defined presently.)

In the context of substitutional alloys, the word _phase_ refers to a physically homogeneous area of the material with uniform chemical composition and physical properties. At the atomic level, this means that, on average, the site-wise concentrations in a particular spatial region are homogeneous and describe some repeating motif of atoms arranged on a lattice. In general, a particular single-phase region of an alloy can be either a `solid solution', in which all lattice sites have equal partial lattice site occupancies and elements mix randomly and homogeneously, or an ordered intermetallic phase (sometimes referred to as an `intermetallic compound'), where a crystallographically ordered, repeating motif of elements can be identified. Naturally, an intermetallic phase can also admit substitutional disorder on one or more of the atomic sites in the repeating motif. In thermal equilibrium, it is possible for more than one phase to be present in an alloy, with the maximum permitted number of phases present determined, in general, by the Gibbs phase rule. An alloy which decomposes into multiple coexisting phases in thermal equilibrium is said to undergo `phase segregation' or `phase decomposition'. (Sometimes this process can also be referred to as one or more phases `precipitating' out of the solid solution.)

For a system with fixed underlying lattice, in the high temperature limit, entropy dictates that atoms should have no lattice site preference. However, as the temperature is lowered, this site symmetry will eventually be broken and either a partial or full site ordering (or decomposition) will become established. The relevant equilibrium phase(s) of an alloy are, in principle, fully determined once the pressure/volume, temperature, and alloy composition (_i.e._ set of alloying elements and their concentrations) have been specified. Freezing the underlying crystal lattice removes the degree of freedom provided by the pressure/volume, and it remains to inspect the probabilities of a given arrangement of atoms occurring at the specified temperature. 

In thermal equilibrium, the relevant probability distribution for determining the likelihood of a given configuration occurring is the Boltzmann distribution. The partition function, \f$Z\f$, is written as
\f[
    Z = \sum_{ \{ \xi_{i \gamma} \} } e^{-\beta E( \{ \xi_{i \gamma} \} )},
\f]
where \f$E( \{ \xi_{i \gamma} \} )\f$ is the energy associated with a configuration, \f$\{ \xi_{i \gamma} \}\f$, and the sum is taken over all possible atomic configurations. The symbol \f$\beta\f$ is defined via \f$\beta=1/k_B T\f$, where \f$T\f$ is temperature, and \f$k_B\f$ is the usual Boltzmann constant. The probability of a given arrangement of atoms occurring is then
\f[
    P(\{ \xi_{i\gamma} \}) = \frac{e^{-\beta E( \{ \xi_{i \gamma} \} )}}{Z}.
\f]
This probability distribution defines an ensemble of configurations distributed according to that probability. If the distribution is known, it is possible to recover various thermodynamic quantities of interest, such as the free energy and specific heat.

However, for any physically realistic form of the alloy internal energy, \f$E( \{ \xi_{i \gamma} \} )\f$, and for all but the smallest of simulation supercells, direct evaluation of the partition function is computationally intractable due to the huge number of configurations which must be considered. This is a particular problem in the context of alloys containing multiple elements, such as high-entropy alloys---those alloys containing four or more elements alloyed in near-equal ratios---as the size of the configuration space grows combinatorially both with the total number of atoms in the simulation cell, as well as with the number of elements present in a given composition. Consequently, it is necessary both to find means by which to evaluate the energy associated with a given arrangement of atoms which are accurate and computationally efficient, as well as to use sampling algorithms to reliably estimate thermodynamic quantities and determine equilibrium phases as a function of temperature without evaluating the partition function in a brute-force manner. Evaluation of the internal energy of the alloy will be discussed now, while specific sampling algorithms (those which are implemented in `BraWl`) will be discussed later in this article.

## Evaluation of the alloy internal energy

In the context of simulations performed on alloy supercells, the energy of a given configuration can be evaluated in a variety of ways. These include first-principles electronic structure calculations using density functional theory (DFT); interatomic potentials, including machine-learned interatomic potentials (MLIPs); and lattice-based models such as cluster expansions (CEs). All of these methods have their advantages and disadvantages. DFT calculations on alloy supercells, while highly accurate, are computationally expensive, rendering this option inviable when a large number of alloy energy evaluations are required for sampling. Calculations using interatomic potentials, MLIPs and CEs are significantly cheaper than direct DFT calculations, but they still frequently require a large DFT training dataset. Additionally, the number of fitting parameters required for these models grows significantly once a multicomponent alloy is considered, making them challenging to train and leading to concerns regarding their accuracy. Finally, it should be stressed that the computational cost of these models often remains prohibitive when a large number of energy evaluations is required.

An alternative, physically intuitive model for the internal energy of an alloy is the Bragg-Williams model, which assumes that the internal energy of an alloy takes a simple, pairwise form. The Bragg-Williams {\color{blue} model} Hamiltonian has the form
\f[
    H(\{\xi_{i\gamma}\}) = \frac{1}{2}\sum_{i \gamma; j\gamma'} V_{i\gamma; j\gamma'} \xi_{i \gamma} \xi_{j \gamma'},
\f]
where \f$V_{i\gamma; j\gamma'}\f$ denotes the effective pair interaction (EPI) between an atom of chemical species \f$\gamma\f$ on lattice site \f$i\f$ and an atom of chemical species \f$\gamma'\f$ on lattice site \f$j\f$. (The factor of \f$\frac{1}{2}\f$ eliminates double-counting in the summation.) For a system of finite size, it is assumed that periodic boundary conditions are applied in all coordinate directions. We note that the form of this Hamiltonian is very similar to that of the Lenz-Ising model, an elementary model in magnetism. 

Generally, the assumption is made that that the EPIs are spatially homogeneous and isotropic, and Eq.~\ref{eq:b-w1} is therefore rewritten as
\f[
    H(\{\xi_{i\gamma}\}) = \frac{1}{2}\sum_{i \gamma} \xi_{i \gamma} ( \sum_{n} \sum_{j \in n(i)} \sum_{\gamma'} V^{(n)}_{\gamma \gamma'} \xi_{j \gamma'} ),
\f]
where the sum over \f$i\f$ remains a sum over lattice sites, but the sum over \f$n\f$ denotes a sum over the coordination shells (nearest-neighbours, next-nearest-neighbours, _etc._) of the lattice. The notation \f$n(i)\f$ is then used to denote the set of lattice sites which are \f$n\f$th nearest-neighbours to site \f$i\f$. Then \f$V^{(n)}_{\gamma \gamma'}\f$ denotes the effective pair interaction between chemical species \f$\gamma\f$ and \f$\gamma'\f$ on coordination shell \f$n\f$. It is reasonable to assume that, for most alloys, the strength of EPIs will tail off quickly with decreasing distance, and the sum over \f$n\f$ can be taken over the first few coordination shells of the underlying lattice type being considered. (This is, of course, equivalent to imposing some radial `cutoff' on an interatomic potential.)

EPIs for the Bragg-Williams Hamiltonian can be obtained using a variety of methods, generally those based on density functional theory.  Similarly to the CE method, it is naturally possible to fit EPIs for a given alloy composition to a set of DFT total energy evaluations on alloy supercells with success. However, most frequently, such EPIs are obtained using the Korringa--Kohn--Rostoker (KKR) formulation of density functional theory, where the coherent potential approximation (CPA) can be used to describe the average electronic structure and consequent internal energy of the disordered alloy. There are then a variety of suitable techniques available for assessing the energetic cost of applied, inhomogeneous chemical perturbations to the CPA reference medium which naturally lead to extraction of EPIs. Such techniques include both the generalised perturbation method (GPM), as well as techniques using the language of concentration waves to describe the atomic-scale chemical fluctuations. Approaches based on concentration waves have been derived for alloys both in the binary and multicomponent settings. Once the EPIs for a given alloy composition are obtained, the phase stability of a particular alloy can be examined using sampling techniques applied to the Bragg-Williams model. This is the purpose of `BraWl` as presented in this work.

## Model limitations

It should be emphasised that there are two key limitations to the above discussion, and to the applicability of the Bragg--Williams model more generally.

The first limitation is the lack of consideration of entropic contributions beyond that made by the configurational entropy, such as vibrational, magnetic, and electronic entropies. Perhaps most important is the consideration of vibrational entropy, which is most relevant when a system transitions between two phases with different elastic properties. For example, if (upon cooling) a system undergoes a phase transition from a disordered phase to an ordered phase, the latter of which is elastically stiffer than the former, the difference in vibrational entropy between the two phases can thermodynamically stabilise the disordered phase and lower the disorder-order transition temperature. Such effects are not accounted for in the Bragg-Williams model directly. However, they can be included in some approximate way by modifying EPIs to account for these effects. More generally, the effect of other entropic contributions, such as magnetic and electronic contributions, should also be considered as appropriate to a given system of study.

The second limitation is the assumption of a fixed underlying crystal lattice. In an alloy where there are atomic size mismatches, there will often be local distortions to the underlying crystal lattice to accommodate the mismatches. (This is because configurations with such local lattice distortions represent true minima of the potential energy surface.) It is understood that these effects can affect predicted disorder-order transition temperatures in alloys.

Neither of the above effects (additional entropic contributions or local lattice distortions) are explicitly considered in the fixed-lattice Bragg-Williams model as implemented in `BraWl`, and users should therefore take due care when interpreting simulation results.
