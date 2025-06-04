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
 - name: Department of Physics, University of Warwick, UK
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
Crucially, the effective pair interactions (EPIs) which appear in the Bragg-Williams Hamiltonian can be obtained _ab initio_ using density functional theory (DFT) calculations.
When appropriate sampling techniques are applied to the Bragg-Williams model, it is possible to explore the configuration space of a given alloy in detail and determine equilibrium phases as a function of temperature, leading to construction of phase diagrams.
Here, we present `BraWl`, a Fortran package implementing a range of conventional and enhanced sampling algorithms for exploration of the phase space of the Bragg-Williams model, facilitating study of diffusional solid-solid transformations in binary and multicomponent alloys.
These sampling algorithms include Metropolis-Hastings Monte Carlo, Wang-Landau sampling, and Nested Sampling.
We demonstrate the capabilities of the package by applying it to some prototypical binary and multicomponent alloys, including high-entropy alloys.

# Statement of need

The Fortran package `BraWl` facilitates simulation of the thermodynamics and phase stability of both binary and multicomponent alloys.
It achieves this by providing implementation of both the Bragg-Williams Hamiltonian (a lattice based model expressing the internal energy of an alloy as a sum of atom-atom effective pair interactions) concurrently with a range of conventional and enhanced sampling techniques for exploration of the alloy configuration space.
The result is a package which can determine phase equilibria as a function of both temperature and alloy composition, which leads to the construction of alloy phase diagrams.
Additionally, the package can be used for extraction of representative equilibrated atomic configurations for visualisation, as well as for use in complementary modelling approaches.
An in-depth discussion of the background, underlying theory, and technical details of the sampling algorithms implemented in the package is given by @naguszewskibrawl2025.

There are a range of existing packages capable of simulating alloy phase equilibria, both open- and closed-source.
Examples of widely-used such packages include ATAT [@vandewallealloy2002], ICET [@angqvisticet2019] and CELL [@rigamonticell2024], though all of these focus primarily on implementation of a generalised cluster expansion, rather than the simpler form of the Bragg-Williams Hamiltonian.
To our knowledge, there is no open-source package specifically focussing on the implementation of a range of sampling algorithms applied to the Bragg-Williams model.
We therefore believe that BraWl fills a gap in the capabilities of the current alloy software ecosystem.
Additionally, we hope that the modular way in which the package is constructed could enable implementation of more complex Hamiltonians, as well as further sampling algorithms in addition to those detailed below, in due course.

## The Bragg-Williams Hamiltonian

A physically intuitive, lattice-based model for the internal energy of a substitutional alloy is the Bragg-Williams model [@braggeffect1934; @braggeffect1935], which assumes that the internal energy of an alloy takes a simple, pairwise form.
We specify a particular arrangement of atoms by a discrete set of site occupation numbers, $\left\{ \xi_{i\gamma} \right\}$, where $\xi_{i \gamma} = 1$ if lattice site $i$ is occupied by an atom of species $\gamma$ and $\xi_{i \gamma} = 0$ otherwise.
The Bragg-Williams Hamiltonian then has the form
\begin{equation}
    H(\{\xi_{i\gamma}\}) = \frac{1}{2}\sum_{i \gamma; j\gamma^{\prime}} V_{i\gamma; j\gamma^{\prime}} \xi_{i \gamma} \xi_{j \gamma^{\prime}},
    \label{eq:b-w1}
\end{equation}
where $V_{i\gamma; j\gamma^{\prime}}$ denotes the effective pair interaction (EPI) between an atom of chemical species $\gamma$ on lattice site $i$ and an atom of chemical species $\gamma^{\prime}$ on lattice site $j$.
(The factor of $\frac{1}{2}$ eliminates double-counting in the summation.)
For a system of finite size, it is assumed that periodic boundary conditions are applied in all coordinate directions.

Generally, the assumption is made that that the EPIs are spatially homogeneous and isotropic, and \autoref{eq:b-w1} is therefore rewritten as
\begin{equation}
    H(\{\xi_{i\gamma}\}) = \frac{1}{2}\sum_{i \gamma} \xi_{i \gamma} \left( \sum_{n} \sum_{j \in n(i)} \sum_{\gamma^{\prime}} V^{(n)}_{\gamma \gamma^{\prime}} \xi_{j \gamma^{\prime}} \right),
    \label{eq:b-w2}
\end{equation}
where the sum over $i$ remains a sum over lattice sites, but the sum over $n$ denotes a sum over the coordination shells (nearest-neighbours, next-nearest-neighbours, _etc._) of the lattice.
The notation $n(i)$ is then used to denote the set of lattice sites which are $n$th nearest-neighbours to site $i$.
Then $V^{(n)}_{\gamma \gamma^{\prime}}$ denotes the effective pair interaction between chemical species $\gamma$ and $\gamma^{\prime}$ on coordination shell $n$.
It is reasonable to assume that, for most alloys, the strength of EPIs will tail off quickly with decreasing distance, and the sum over $n$ can be taken over the first few coordination shells of the underlying lattice type being considered.
(This is, of course, equivalent to imposing some radial 'cutoff' on an interatomic potential.)

EPIs for the Bragg-Williams Hamiltonian can be obtained using a variety of methods, generally those based on density functional theory.
As examples of such methods, we highlight the recovery of such interactions obtained by fitting to a set of DFT total energy evaluations on alloy supercells [@zhangrobust2020; @liumonte2021], the generalised perturbation method (GPM) [@ducastellegeneralized1976; @rubanatomic2004], and techniques using the language of concentration waves to describe the atomic-scale chemical fluctuations [@singhatomic2015; @khanstatistical2016].
Once the EPIs for a given alloy composition are obtained, the phase stability of a particular alloy can be examined using sampling techniques applied to the Bragg-Williams model.
This is the purpose of `BraWl` as presented in this work.

# Sampling algorithms

`BraWl` implements a range of conventional and enhanced sampling algorithms for exploration of the alloy configuration space.
Throughout, the package defaults to performing _swaps_ of pairs of atoms in the simulation cell, to conserve the overall concentration of each chemical species present in the simulation.
At present, the implemented sampling algorithms are:

1. **The Metropolis-Hastings Monte Carlo algorithm.** The Metropolis-Hastings algorithm allows for a system of interest to follow a chain of states which evolve to, and sample, an equilibrium ensemble [@metropolisequation1953; landauguide2014]. It is a useful method for obtaining the equilibrium state of a system at a given temperature.
2. **Wang-Landau sampling.** Wang-Landau sampling is a flat histogram method which provides a means for high throughput calculation of phase diagrams for atomistic/lattice model systems [@wangefficient2001]. The method allows for direct computation of an estimate of the density of states in energy, $g(E)$, and hence the partition function of a system of interest. From the partition function, thermodynamic quantities at any temperature of interest can then be obtained provided one has prior knowledge of the minimum and maximum energy relevant to those temperatures.
3. **Nested Sampling.** Nested sampling (NS) is powerful Bayesian inference technique [@skillingnested2004; @skillingnested2006] adapted to sample the potential energy surface of atomistic systems [@ashtonnested2022], giving direct access to the partition function at arbitrary temperatures for comprehensive thermodynamic analysis, without relying on advance knowledge of relevant structures or the range of energies accessible to them [@partayefficient2010; @partaynested2021].

# Physical quantities of interest

`BraWl` can extract a range of quantities of interest from a given alloy simulation.
The relevant quantities are: 

- **Internal energy.** For a given lattice type, system size, alloy composition, and set of atom-atom effective pair interactions, `BraWl` can evaluate the total energy associated with the alloy configuration (\autoref{eq:b-w1}). At the time of writing, for speed, common lattice types (fcc, bcc, simple cubic, _etc._) are hard-coded, with the intention that the range of implemented lattice types will be expanded over time as necessary.
- **Heat capacity.** The isochoric (fixed volume) heat capacity at a given temperature, $C_V(T)$, is a useful quantity for identifying phase transitions, as a plot of the simulation heat capacity as a function of temperature is expected to show a local peak at the temperature at which the transition occurs. Within `BraWl`, the heat capacity is calculated via $$C_V(T) = \frac{\langle E^2\rangle - \langle E\rangle^2}{k_BT^2},$$ where $k_B$ is the usual Boltzmann constant, $E$ is the simulation energy, and $\langle \cdot \rangle$ denote thermodynamic averages obtained using the relevant sampling algorithm.
- **Atomic short-range order parameters.** To assess local atom-atom correlations in a simulation, `BraWl` can calculate the Warren-Cowley atomic short-range order parameters [@cowleyapproximate1950; @cowleyshortrange1965], adapted to the multicomponent setting, defined as $$\alpha^{\gamma \gamma'}_n=1-\frac{P^{\gamma \gamma'}_n}{c_{\gamma'}},$$ where $n$ refers to the $n$th coordination shell, $P^{\gamma \gamma'}_n$ is the conditional probability of an atom of type $q$ neighbouring an atom of type $p$ on shell $n$, and $c_q$ is the total concentration of atom type $q$. When $\alpha^{\gamma \gamma'}_n > 0$, $p$-$q$ pairs are disfavoured on shell $n$ and, when $\alpha^{\gamma \gamma'}_n < 0$ they are favoured. The value $\alpha^{\gamma \gamma'}_n = 0$ corresponds to the ideal, maximally disordered solid solution.
- **Atomic long-range order parameters.** Over a simulation run (for example using the Metropolis algorithm), `BraWl` can calculate the average occupancy of a lattice site, _i.e._ the probability of an atom of a particular chemical species sitting on that site. This capability was first demonstrated on simulations of Fe-Ga alloys [@marchantab2021].

# Example Applications

BraWL has been used, with success, to study the phase behaviour of a range of binary and multicomponent alloys, for example the binary Fe-Ga system (Galfenol) [@marchantab2021], the Fe-Ni system [@woodgateint2024], the Cantor-Wu medium- and high-entropy alloys [@woodgatecompositional2022; @woodgateinterplay2023], the refractory high-entropy alloys [@woodgateshortrange2023; @woodgatecompetition2024], the Al$_x$CrFeCoNi system [@woodgatestructure2024], and the AlTiVNb and AlTiCrMo refractory high-entropy superalloys [@woodgateemergent2025].
The package has also been used to generate atomic configurations for training datasets for machine-learned interatomic potentials, for example for the prototypical austenitic stainless steel, Fe$_7$Cr$_2$Ni [@shenoycollinearspin2024].

As an example of the Metropolis-Hastings Monte Carlo algorithm, we consider its application to the binary FeNi alloy, first discussed by @woodgateint2024.
\autoref{fig:feniequilibration} shows the internal energy and conditional pair probabilities (quantifying ASRO) of a simulation cell containing 256 atoms as a function of the number of Metropolis-Hastings 'sweeps', where a sweep refers to performing a number of trial Metropolis-Hastings moves equal to the number of atoms in the simulation cell.
The simulation is performed at 300 K, below the alloy's L1$_0$ disorder-order transition temperature.
The L1$_0$ phase is a structure where $2/3$ of the nearest neighbours of Fe (Ni) atoms are Ni (Fe) atoms, and where none of the next-nearest neighbours of Fe (Ni) atoms are Ni (Fe) atoms.
It can be seen that this ordering is swiftly established as the number of Monte Carlo sweeps increases, albeit with some remaining thermal noise.

![Evolution of the simulation internal energy (top panel) and conditional pair probabilities (bottom panel) for an Fe$_{0.5}$Ni$_{0.5}$ alloy as a function of the number of Metropolis-Hastings sweeps at a simulation temperature of $T=300$ K. One 'sweep' is one trial move per atom in the system. Beyond approximately 100 sweeps, the system can be seen to have reached equilibrium, with L1$_0$ order established.\label{fig:feniequilibration}](feniequilibration.pdf){ width=55% }

As an example of Wang-Landau sampling, we consider its application to the AlTiCrMo refractory high-entropy superalloy, first discussed by @woodgateemergent2025, for which results are shown in \autoref{fig:wlAlTiCrMo}.
The top panel shows calculated energy probability distributions (histograms) at various temperatures, while the bottom panel shows the simulation heat capacity and evolution of the Warren-Cowley ASRO parameters as a function of temperature.
The high-temperature peak in the heat capacity data is associated with the experimentally observed B2 crystallographic ordering.

![Plots of energy probability distributions, Warren-Cowley ASRO parameters ($\alpha_n^{pq}$) and simulation heat capacity ($C$) as a function of temperature for AlTiCrMo obtained using lattice-based Monte Carlo simulations employing Wang-Landau sampling. Here, show $\alpha_n^{pq}$ only for $n = 1$. The zero of the energy scale for the energy histograms is set to be equal to the average internal energy of the alloy obtained at a simulation temperature of 3000 K.\label{fig:wlAlTiCrMo}](wlAlTiCrMo.pdf){ width=55% }

Finally, as an example of application of the Nested Sampling algorithm, we consider its application to the AlCrFeCoNi high-entropy alloy, first discussed by @woodgatestructure2024.
\autoref{fig:nsalcrfeconi} plots the internal energy, $E$, and isochoric heat capacity, $C_V$, obtained for the equiatomic, fcc, AlCrFeCoNi system.
The simulation cell contained 108 atoms.
The initial peak in the heat capacity encountered upon cooling from high temperature is associated with an L1$_2$ ordering driven by Al, with subsequent peaks indicating eventual decomposition into multiple competing phases, which is understood to be consistent with experimental data for this system.

![Internal energy, $E$, and isochoric heat capacity, $C_V$, obtained using the Nested Sampling algorithm applied to the equiatomic, fcc, AlCrFeCoNi high-entropy alloy. The simulation cell contained 108 atoms. Upon cooling, the initial peak in the heat capacity is associated with an L1$_2$ ordering driven by Al, with subsequent peaks indicating eventual decomposition into multiple competing phases.\label{fig:nsalcrfeconi}](nsalcrfeconi.pdf){ width=55% }

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

