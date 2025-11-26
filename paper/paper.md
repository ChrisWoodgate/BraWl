---
title: 'BraWl: Simulating the thermodynamics and phase stability of multicomponent alloys using conventional and enhanced sampling techniques'
tags:
  - Physics
  - Chemistry
  - Materials Science
  - Alloys
  - Statistical Mechanics
  - Thermodynamics
authors:
  - name: Hubert J. Naguszewski
    orcid: 0009-0006-1670-5437
    affiliation: 1
  - name: Livia B. Pártay
    orcid: 0000-0003-3249-3586
    affiliation: 2
  - name: David Quigley
    orcid: 0000-0003-4750-4372
    affiliation: 1
  - name: Christopher D. Woodgate
    corresponding: true
    orcid: 0000-0003-4876-5558
    affiliation: 3
affiliations:
 - name: Department of Physics, University of Warwick, United Kingdom
   index: 1
   ror: 01a77tt86
 - name: Department of Chemistry, University of Warwick, UK
   index: 2
   ror: 01a77tt86
 - name: School of Physics, University of Bristol, UK
   index: 3
   ror: 0524sp257
date: 13 May 2025
bibliography: brawl.bib
---

# Summary

Many technologically relevant materials, both structural and functional, are 'alloys' - systems in which two or more (typically) metallic elements are combined to produce a new material with desirable physical properties.
In a substitutional alloy, there is a fixed underlying crystal lattice, while the probability of a given constituent element of the alloy occupying a particular lattice site is determined by thermodynamic considerations.
In accordance with these considerations, atoms in an alloy can arrange themselves differently depending on the precise chemical composition and processing conditions.
Frequently, a mixture of elements will form a regular crystalline lattice with substitutional disorder (a 'solid solution') at high temperature, before atomic short- and long-range order emerges as the material is cooled.
The nature of atomic arrangements in a material determines many important physical properties.
For a given combination of elements, it is therefore crucial to understand the nature of atomic ordering in a material, as well as the temperature at which it emerges upon cooling, to guide materials processing strategies.
One physically intuitive model for the internal energy of an alloy is the Bragg-Williams model, which assumes that atoms in the alloy interact in a pairwise manner.
Crucially, the atom-atom effective pair interactions (EPIs) which appear in the Bragg-Williams Hamiltonian can be obtained _ab initio_ using density functional theory (DFT) calculations.
When appropriate sampling techniques are applied to the Bragg-Williams model, it is possible to explore the configuration space of a given alloy in detail and determine equilibrium phases as a function of temperature, leading to construction of phase diagrams.
Here, we present `BraWl`, a Fortran package implementing a range of conventional and enhanced sampling algorithms for exploration of the phase space of the Bragg-Williams model, facilitating study of diffusional solid-solid transformations in binary and multicomponent alloys.
These sampling algorithms include Metropolis-Hastings Monte Carlo, Wang-Landau sampling, and Nested Sampling.
We demonstrate the capabilities of the package by applying it to some prototypical binary and multicomponent alloys, including high-entropy alloys.

# Statement of need

The Fortran package `BraWl` facilitates simulation of the thermodynamics and phase stability of both binary and multicomponent alloys.
It achieves this by providing an implementation of the Bragg-Williams Hamiltonian (a lattice-based model expressing the internal energy of an alloy as a sum of atom-atom effective pair interactions) concurrently with implementations of a range of conventional and enhanced sampling techniques for exploration of the alloy configuration space.
The result is a package which can determine phase equilibria as a function of both temperature and alloy composition, which leads to the construction of alloy phase diagrams.
Additionally, the package can be used for extraction of representative equilibrated atomic configurations for visualisation, as well as for use in complementary modelling approaches.
An in-depth discussion of the background, underlying theory, and technical details of the sampling algorithms implemented in the package is given by @naguszewskibrawl2025.
Much of this discussion is also available in the package's documentation.

There are a range of existing packages capable of simulating alloy phase equilibria, both open- and closed-source.
Examples of widely-used such packages include ATAT [@vandewallealloy2002], ICET [@angqvisticet2019] and CELL [@rigamonticell2024], though all of these focus primarily on implementation of a generalised cluster expansion, rather than the simpler form of the Bragg-Williams Hamiltonian.
To our knowledge, there is no open-source package specifically focussing on the implementation of a range of sampling algorithms applied to the Bragg-Williams model.
We therefore believe that `BraWl` fills a gap in the capabilities of the current alloy software ecosystem.
Additionally, we hope that the modular way in which the package is constructed could enable implementation of more complex Hamiltonians, as well as further sampling algorithms in addition to those detailed below, in due course.

# Sampling algorithms

`BraWl` implements a range of conventional and enhanced sampling algorithms for exploration of the alloy configuration space.
Throughout, the package defaults to performing _swaps_ of pairs of atoms in the simulation cell, to conserve the overall concentration of each chemical species present in the simulation.
At present, the implemented sampling algorithms are:

1. **The Metropolis-Hastings Monte Carlo algorithm.** The Metropolis-Hastings algorithm allows for a system of interest to follow a chain of states which evolve to, and sample, an equilibrium ensemble [@metropolisequation1953; @landauguide2014]. It is a useful method for obtaining the equilibrium state of a system at a given temperature.
2. **Wang-Landau sampling.** Wang-Landau sampling is a flat histogram method which provides a means for high-throughput calculation of phase diagrams for atomistic/lattice model systems [@wangefficient2001]. The method allows for direct computation of an estimate of the density of states in energy, $g(E)$, and hence the partition function of a system of interest. From the partition function, thermodynamic quantities at any temperature of interest can then be obtained provided one has prior knowledge of the minimum and maximum energy relevant to those temperatures.
3. **Nested Sampling.** Nested sampling (NS) is powerful Bayesian inference technique [@skillingnested2004; @skillingnested2006] adapted to sample the potential energy surface of atomistic systems [@ashtonnested2022], giving direct access to the partition function at arbitrary temperatures for comprehensive thermodynamic analysis, without relying on advance knowledge of relevant structures or the range of energies accessible to them [@partayefficient2010; @partaynested2021].

# Physical quantities of interest

`BraWl` can extract a range of quantities of interest from a given alloy simulation.
The relevant quantities are: 

- **Internal energy.** For a given lattice type, system size, alloy composition, and set of atom-atom effective pair interactions, `BraWl` can evaluate the total energy associated with the alloy configuration. At the time of writing, for speed, common lattice types (fcc, bcc, simple cubic, _etc._) are hard-coded, with the intention that the range of implemented lattice types will be expanded over time as necessary.
- **Specific heat.** The isochoric (fixed volume) specific heat at a given temperature, $C_V(T)$, is a useful quantity for identifying phase transitions, as a plot of the specific heat as a function of temperature is expected to show a local peak at the temperature at which the transition occurs. Within `BraWl`, the specific heat is calculated via $$C_V(T) = \frac{\langle E^2\rangle - \langle E\rangle^2}{k_\textrm{B}T^2},$$ where $k_\textrm{B}$ is the usual Boltzmann constant, $E$ is the simulation energy, and $\langle \cdot \rangle$ denote thermodynamic averages obtained using the relevant sampling algorithm.
- **Atomic short-range order parameters.** To assess local atom-atom correlations in a simulation, `BraWl` can calculate the Warren-Cowley atomic short-range order parameters [@cowleyapproximate1950; @cowleyshortrange1965], adapted to the multicomponent setting, defined as $$\alpha^{\gamma \gamma'}_n=1-\frac{P^{\gamma \gamma'}_n}{c_{\gamma'}},$$ where $n$ refers to the $n$th coordination shell, $P^{\gamma \gamma'}_n$ is the conditional probability of an atom of type $q$ neighbouring an atom of type $p$ on shell $n$, and $c_q$ is the total concentration of atom type $q$. When $\alpha^{\gamma \gamma'}_n > 0$, $p$-$q$ pairs are disfavoured on shell $n$ and when $\alpha^{\gamma \gamma'}_n < 0$ they are favoured. The value $\alpha^{\gamma \gamma'}_n = 0$ corresponds to the ideal, maximally disordered solid solution.
- **Atomic long-range order parameters.** Over a simulation run (for example using the Metropolis algorithm), `BraWl` can calculate the average occupancy of a lattice site, _i.e._ the probability of an atom of a particular chemical species sitting on that site. This capability was first demonstrated on simulations of Fe-Ga alloys [@marchantab2021].

# Example Applications

BraWL has been used, with success, to study the phase behaviour of a range of binary and multicomponent alloys, for example the binary Fe-Ga system (Galfenol) [@marchantab2021], the Fe-Ni system [@woodgateint2024], the Cantor-Wu medium- and high-entropy alloys [@woodgatecompositional2022; @woodgateinterplay2023], the refractory high-entropy alloys [@woodgateshortrange2023; @woodgatecompetition2024], the Al$_x$CrFeCoNi system [@woodgatestructure2024], and the AlTiVNb and AlTiCrMo refractory high-entropy superalloys [@woodgateemergent2025].
Additionally, the package has been used to generate atomic configurations with physically motivated atomic short-range order (ASRO) and/or atomic long-range order (ALRO) for subsequent study using a range of other simulation techniques.
We highlight examples of its use in generating a training dataset for a machine-learned interatomic potential for the prototypical austenitic stainless steel, Fe$_7$Cr$_2$Ni [@shenoycollinearspin2024], as well as its use in generating configurations for use in a transition state study for ferromagnetic Fe-Ni alloys [@fisherlattice2025].
In this work, we explicitly consider several illustrative examples of the results which can be obtained using the sampling algorithms outlined above applied to the Bragg-Williams model as implemented in the package. Throughout these examples, the atom-atom EPIs for the Bragg-Williams model are obtained using the $S^{(2)}$ theory for multicomponent alloys [@khanstatistical2016; @woodgatemodelling2024].

As an example of the Metropolis-Hastings Monte Carlo algorithm, we consider its application to the binary FeNi alloy, first discussed by @woodgateint2024.
\autoref{fig:feniequilibration} shows the internal energy and conditional pair probabilities (quantifying ASRO) of a simulation cell containing 256 atoms as a function of the number of Metropolis-Hastings 'sweeps', where a sweep refers to performing a number of trial Metropolis-Hastings moves equal to the number of atoms in the simulation cell.
The simulation is performed at 300 K, below the alloy's $\textrm{L}1_0$ disorder-order transition temperature.
The $\textrm{L}1_0$ phase is a structure where $2/3$ of the nearest neighbours of Fe (Ni) atoms are Ni (Fe) atoms, and where none of the next-nearest neighbours of Fe (Ni) atoms are Ni (Fe) atoms.
It can be seen that this ordering is swiftly established as the number of Monte Carlo sweeps increases, albeit with some remaining thermal noise.

![Evolution of the simulation internal energy (top panel) and conditional pair probabilities (bottom panel) for an Fe$_{0.5}$Ni$_{0.5}$ alloy as a function of the number of Metropolis-Hastings sweeps at a simulation temperature of $T=300$ K. One 'sweep' is one trial move per atom in the system. Beyond approximately 100 sweeps, the system can be seen to have reached equilibrium, with $\textrm{L}1_0$ order established.\label{fig:feniequilibration}](feniequilibration.pdf){ width=70% }

As an example of Wang-Landau sampling, we consider its application to the AlTiCrMo refractory high-entropy superalloy, first discussed by @woodgateemergent2025, for which results are shown in \autoref{fig:wlAlTiCrMo}.
The top panel shows calculated energy probability distributions (histograms) at various temperatures, while the bottom panel shows the specific heat and evolution of the Warren-Cowley ASRO parameters as a function of temperature.
The high-temperature peak in the specific heat data is associated with the experimentally observed B2 crystallographic ordering.

Finally, as an example of application of the Nested Sampling algorithm, we consider its application to the AlCrFeCoNi high-entropy alloy, first discussed by @woodgatestructure2024.
\autoref{fig:nsalcrfeconi} plots the internal energy, $E$, and isochoric specific heat, $C_V$, obtained for the equiatomic, fcc, AlCrFeCoNi system.
The simulation cell contained 108 atoms.
The initial peak in the specific heat encountered upon cooling from high temperature is associated with an $\textrm{L}1_2$ ordering driven by Al, with subsequent peaks indicating eventual decomposition into multiple competing phases, which is understood to be consistent with experimental data for this system.

![Plots of energy probability distributions, Warren-Cowley ASRO parameters ($\alpha_n^{pq}$) and specific heat ($C_V$) as a function of temperature for AlTiCrMo obtained using lattice-based Monte Carlo simulations employing Wang-Landau sampling. Here, show $\alpha_n^{pq}$ only for $n = 1$. The zero of the energy scale for the energy histograms is set to be equal to the average internal energy of the alloy obtained at a simulation temperature of 3000 K.\label{fig:wlAlTiCrMo}](wlAlTiCrMo.pdf){ width=70% }

![Internal energy, $E$, and isochoric specific heat, $C_V$, obtained using the Nested Sampling algorithm applied to the equiatomic, fcc, AlCrFeCoNi high-entropy alloy. The simulation cell contained 108 atoms. Upon cooling, the initial peak in the specific heat is associated with an $\textrm{L}1_2$ ordering driven by Al, with subsequent peaks indicating eventual decomposition into multiple competing phases.\label{fig:nsalcrfeconi}](nsalcrfeconi.pdf){ width=70% }

# Performance

Fundamentally, all of the outlined sampling algorithms (Metropolis-Hastings Monte Carlo, Wang-Landau sampling, Nested Sampling) require a large number of evaluations of the alloy internal energy during the sampling process.
The time taken for evaluation of the model Hamiltonian is therefore the key factor in determining the computational cost of a sampling run on a given system.

In the context of the present work, two key observations should be made regarding how the computational cost of a simulation scales with the size of a simulation cell containing $N$ atoms.
The first observation relates to how evaluation of the Bragg-Williams model Hamiltonian scales with respect to system size.
As the atom-atom EPIs are assumed to be of finite range and each atom therefore only interacts with other atoms within some finite cutoff radius, evaluation of the model Hamiltonian scales approximately linearly with $N$ for systems which are sufficiently large.
The second relates to how the number of configurations accessible to a simulation cell scales with $N$.
Formally, for an $s$-species alloy with $N_1$ atoms of species 1, $N_2$ atoms of species 2, and so on, the number of possible atomic arrangements of these atoms on the $N$ available lattice sites, $n_\textrm{configs}$, is $$ n_\textrm{configs} = \frac{N!}{N_1!\times N_2! \times \dots \times N_s!}.$$ For an equiatomic binary alloy ($N_1 = N_2 = N/2$) this means that for a 32-atom cell there are approximately $6.01 \times 10^{8}$ possible configurations, while for a 256-atom cell this number becomes approximately $5.79 \times 10^{75}$.
For an equatomic quarternary system ($N_1 = N_2 = N_3 = N_4 = N/4$), these numbers become approximately $7.93 \times 10^{29}$ and $3.31 \times 10^{150}$ respectively.
Accordingly, for large systems with many elements in the alloy composition, longer sampling runs are anticipated to be needed to ensure the configuration space has been adequately explored and for any calculated quantities to be well-converged.
We emphasise that, naturally, any given sampling algorithm is not anticipated to visit all possible system configurations, only a representative sample that is converged.

To assess the performance of the `BraWl` package, we measure the time taken for a Metropolis-Hastings run at a single temperature, as this quantifies the rate of sampling and thus is crucial in determining the performance of all of the considered sampling algorithms.
We also focus exclusively on serial performance.
This is because the Metropolis-Hastings algorithm in not inherently parallelisable unless several independent simulations (_e.g._ at different temperatures, compositions, or using different random number seeds) are run independently, while parallel implementation/performance of the Wang-Landau sampling  [@naguszewskioptimal2025] and Nested Sampling algorithms [@unglertreplica2025] is discussed elsehwhere.
(Also, at the time of writing, only conventional, serial Nested Sampling is implemented in `BraWl`.)

As a benchmark, we consider a simulation on the equiatomic AlTiVNb alloy using the EPIs obtained by @woodgateemergent2025, which include interactions up to and including sixth-nearest neighbours, and a simulation cell consisting of $8\times8\times8$ bcc unit cells with $N=1024$ in total.
The sampling run consisted of a Metropolis-Hastings run with 1,024,000 trial Metropolis-Hastings swaps ($10^3$ trial moves per atom) during an initial burn-in phase, followed by a sampling run of 10,240,000 trial atomic swaps ($10^4$ trial moves per atom) during which statistics are gathered.
The energy of the simulation was stored every 1,024 trial moves for a total of 10,000 decorrelated energy samples, while the Warren-Cowley parameters up to and including the second coordination shell of the lattice were evaluated and stored every 10,240 trial moves for a total of 1,000 decorrelated ASRO samples.
(Typically, fewer samples are required to reliably converge ASRO results than are required to estimate quantities such as as the specific heat.) For a simulation performed locally in serial on an Apple M3 Pro CPU (ARM architecture) with `BraWl` built using the GNU compiler collection `v15.1.0` at `-O3` compiler optimisation, this simulation run took an average of 12.5 s across three attempts, indicating an average of approximately 900,000 trial atomic swaps per second.
On an 4.9 GHz 12th generation Intel Core i7-12700 CPU (x86 architecture) with `BraWl` built using the 2022 Intel Fortran compiler at `-O2` compiler optimisation, this simulation run took an average 31.8 s across three attempts, indicating an average of approximately 350,000 trial atomic swaps per second.

# Acknowledgements

H.J.N. and C.D.W. acknowledge invaluable training in scientific software development provided by Dr C. S. Brady and Dr H. Ratcliffe in their Introduction to Scientific Software Development course at the University of Warwick.
L.B.P. and C.D.W. acknowledge useful discussions with Dr Georgia A. Marchant (Department of Physics and Astronomy, Uppsala University, Sweden).
Finally, C.D.W. acknowledges guidance from Prof. Julie B. Staunton (Department of Physics, University of Warwick, UK) during the early stages of this work.

H.J.N. was supported by a studentship within the Centre for Doctoral Training in Modelling of Heterogeneous Systems, funded by the UK Engineering and Physical Sciences Research Council (EPSRC), Grant EP/S022848/1.
L.B.P. acknowledges support from the EPSRC through the individual Early Career Fellowship, Grant EP/T000163/1.
C.D.W. acknowledges support from an EPSRC Doctoral Prize Fellowship at the University of Bristol, Grant EP/W524414/1. 

Computing facilities for development and testing of the package were provided by the [Advanced Computing Research Centre (ACRC)](https://www.bristol.ac.uk/acrc/) of the University of Bristol, and by the [Scientific Computing Research Technology Platform (SCRTP)](https://warwick.ac.uk/research/rtp/sc/) of the University of Warwick.
We also acknowledge use of the Sulis Tier 2 HPC platform.
Sulis is funded by EPSRC Grant EP/T022108/1 and the HPC Midlands+ consortium.

# References

