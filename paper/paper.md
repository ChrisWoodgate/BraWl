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
  - name: Livia B. Partay
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
 - name: Department of Physics, University of Warwick, Coventry, CV4 7AL, UK
   index: 1
   ror: 01a77tt86
 - name: Department of Chemistry, University of Warwick, Coventry, CV4 7AL, UK
   index: 2
   ror: 01a77tt86
 - name: H.H. Wills Physics Laboratory, University of Bristol, Royal Fort, Bristol, BS8 1TL, United Kingdom
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
Before outlining the specific need for this package, we first review some elementary details of alloy thermodynamics and statistical mechanics.

## Alloy thermodynamics

In a substitutional alloy with fixed underlying lattice, a particular arrangement of atoms can be specified by a discrete set of site occupation numbers, $`\left\{ \xi_{i\gamma} \right\}`$, where $`\xi_{i \gamma} = 1`$ if site $`i`$ is occupied by an atom of species $`\gamma`$ and $`\xi_{i \gamma} = 0`$ otherwise.
The lattice index $`i`$ takes values in the range 1 to $`N`$, where $`N`$ is the total number of lattice sites in the system, while the species index $`\gamma`$ takes values in the range 1 to $`s`$, where $`s`$ is the number of chemical species (elements) present in the alloy composition.
The physical constraint that each lattice site is occupied by one (and only one)  atom is expressed as $`\sum_{\gamma} \xi_{i\gamma} = 1`$, while the total concentration of a given chemical species is given by $`c_\gamma = \frac{1}{N} \sum_i \xi_{i\gamma}`$.
(Naturally, vacancies can be treated in this framework by considering them as an additional chemical species present at a very low concentration.)
If we consider an ensemble of alloy configurations, we can then define the site-wise concentrations, $`\left\{ c_{i\gamma}\right\}`$, as the average value of the site occupation numbers across the ensemble, $`c_{i\gamma} = \langle \xi_{i\gamma} \rangle`$, where $`\langle \cdot \rangle`$ denotes the average taken over the ensemble.

![Illustrations of potential states of a substitutional alloy in thermal equilibrium. A solid solution (a) is a state where lattice sites are occupied at random by elements of different chemical species. An ordered intermetallic compound (b) has an identifiable regular, repeating motif of atoms. A system may also undergo phase decomposition (c), where pairs of elements phase segregate from one another. In the multicomponent setting (d) there can be many possible competing phases.\label{fig:alloy_phase_cartoon}](brawlalloyillustrations.pdf)

In the context of substitutional alloys, the word _phase_ refers to a physically homogeneous area of the material with uniform chemical composition and physical properties.
At the atomic level, this means that, on the average, the site-wise concentrations in a particular spatial region are homogeneous and describe some repeating motif of atoms arranged on a lattice.
In general, a particular single-phase region of an alloy can be either a 'solid solution', in which all lattice sites have equal partial lattice site occupancies and elements mix randomly and homogeneously, or an ordered intermetallic phase (sometimes referred to as an 'intermetallic compound'), where a crystallographically ordered, repeating motif of elements can be identified.
Naturally, an intermetallic phase can also admit substitutional disorder on one or more of the atomic sites in the repeating motif.
In thermal equilibrium, it is possible for more than one phase to be present in an alloy, with the maximum permitted number of phases present determined, in general, by the Gibbs phase rule.
An alloy which decomposes into multiple coexisting phases in thermal equilibrium is said to undergo 'phase segregation' or 'phase decomposition'.
(Sometimes this process can also be referred to as one or more phases 'precipitating' out of the solid solution.)
Some illustrations of potential alloy phases for a toy, 2D alloy in both binary and multicomponent settings, are provided in Fig. \autoref{fig:alloyphasecartoon}


For a system with fixed underlying lattice, in the high temperature limit, entropy dictates that atoms should have no lattice site preference.
However, as the temperature is lowered, this site symmetry will eventually be broken and either a partial or full site ordering (or decomposition) will become established.
The relevant equilibrium phase(s) of an alloy are, in principle, fully determined once the pressure/volume, temperature, and alloy composition (_i.e._ set of alloying elements and their concentrations) have been specified.
Freezing the underlying crystal lattice removes the degree of freedom provided by the pressure/volume, and it remains to inspect the probabilities of a given arrangement of atoms occurring at the specified temperature.


In the canonical ensemble (fixed numbers of particles) in thermal equilibrium, the relevant probability distribution is the Boltzmann distribution.
The partition function, $`Z`$, is written as
\begin{equation}
    Z = \sum_{\left\{ \xi_{i \gamma} \right\}} e^{-\beta E\left( \left\{ \xi_{i \gamma} \right\} \right)},
\end{equation}
where $`E\left( \left\{ \xi_{i \gamma} \right\} \right)`$ is the energy associated with a configuration, $`\left\{ \xi_{i \gamma} \right\}`$, and the sum is taken over all possible atomic configurations.
The symbol $`\beta`$ is defined via $`\beta=1/k_B T`$, where $`T`$ is temperature, and $`k_B`$ is the usual Boltzmann constant.
The probability of a given arrangement of atoms occurring is then
\begin{equation}
    P\left(\left\{ \xi_{i\gamma} \right\}\right) = \frac{e^{-\beta E\left( \left\{ \xi_{i \gamma} \right\} \right)}}{Z}.
\end{equation}
If this probability distribution is known, it is then also possible to recover various thermodynamic quantities of interest, such as the free energy and heat capacity.


However, for any physically realistic form of the alloy internal energy, $`E\left( \left\{ \xi_{i \gamma} \right\} \right)`$, and for all but the smallest of simulation supercells, direct evaluation of the partition function is computationally intractable due to the huge number of configurations which must be considered.
This is a particular problem in the context of alloys containing multiple elements, such as high-entropy alloys - those alloys containing four or more elements alloyed in near-equal ratios `[@georgehighentropy2019]` - as the size of the configuration space grows combinatorially both with the total number of atoms in the simulation cell, as well as with the number of elements present in a given composition `[@zhangroadmap2025]`.
Consequently, it is necessary both to find means by which to evaluate the energy associated with a given arrangement of atoms which are accurate and computationally efficient, as well as to use sampling algorithms to reliably estimate thermodynamic quantities and determine equilibrium phases as a function of temperature `[@ferrarifrontiers2020; @ferrarisimulating2023]` without evaluating the partition function in a brute-force manner.
Evaluation of the internal energy of the alloy will be discussed now, while specific sampling algorithms (those which are implemented in `BraWl`) will be discussed later in this article.

## Evaluation of the alloy internal energy

In the context of simulations performed on alloy supercells, the energy of a given configuration can be evaluated in a variety of ways `[@widommodeling2018]`.
These include first-principles electronic structure calculations using density functional theory (DFT) `[@eisenbachfirstprinciples2019]`; interatomic potentials, including machine-learned interatomic potentials (MLIPs) `[@rosenbrockmachinelearned2021]`; and lattice-based models such as cluster expansions (CEs) `[@ekborgtannerconstruction2024]`.
All of these methods have their advantages and disadvantages.
DFT calculations on alloy supercells, while highly accurate, are computationally expensive, rendering this option inviable when a large number of alloy energy evaluations are required for sampling.
Calculations using interatomic potentials, MLIPs and CEs are significantly cheaper than direct DFT calculations, but they still frequently require a large DFT training dataset.
Additionally, the number of fitting parameters required for these models grows significantly once a multicomponent alloy is considered, making them challenging to train and leading to concerns regarding their accuracy.
Finally, it should be stressed that the computational cost of these models often remains non-negligible when a large number of energy evaluations is required.


An alternative, physically intuitive model for the internal energy of an alloy is the Bragg-Williams model `[@braggeffect1934, @braggeffect1935]`, which assumes that the internal energy of an alloy takes a simple, pairwise form.
The Bragg-Williams Hamiltonian has the form
\begin{equation}
    H(\{\xi_{i\gamma}\}) = \frac{1}{2}\sum_{i \gamma; j\gamma^{\prime}} V_{i\gamma; j\gamma^{\prime}} \xi_{i \gamma} \xi_{j \gamma^{\prime}},
    \label{eq:b-w1}
\end{equation}
where $`V_{i\gamma; j\gamma^{\prime}}`$ denotes the effective pair interaction (EPI) between an atom of chemical species $`\gamma`$ on lattice site $`i`$ and an atom of chemical species $\gamma^{\prime}$ on lattice site $`j`$.
(The factor of $`\frac{1}{2}`$ eliminates double-counting in the summation.)
For a system of finite size, it is assumed that periodic boundary conditions are applied in all coordinate directions.
We note that the form of this Hamiltonian is very similar to that of the Lenz-Ising model, an elementary model in magnetism `[@brushhistory1967]`. 

Generally, the assumption is made that that the EPIs are spatially homogeneous and isotropic, and Eq. \autoref{eq:b-w1} is therefore rewritten as
\begin{equation}
    H(\{\xi_{i\gamma}\}) = \frac{1}{2}\sum_{i \gamma} \xi_{i \gamma} \left( \sum_{n} \sum_{j \in n(i)} \sum_{\gamma^{\prime}} V^{(n)}_{\gamma \gamma^{\prime}} \xi_{j \gamma^{\prime}} \right),
    \label{eq:b-w3}
\end{equation}
where the sum over $`i`$ remains a sum over lattice sites, but the sum over $`n`$ denotes a sum over the coordination shells (nearest-neighbours, next-nearest-neighbours, _etc._) of the lattice.
The notation $`n(i)`$ is then used to denote the set of lattice sites which are $`n`$th nearest-neighbours to site $`i`$.
Then $`V^{(n)}_{\gamma \gamma^{\prime}}`$ denotes the effective pair interaction between chemical species $`\gamma`$ and $`\gamma^{\prime}`$ on coordination shell $`n`$.
It is reasonable to assume that, for most alloys, the strength of EPIs will tail off quickly with decreasing distance, and the sum over $`n`$ can be taken over the first few coordination shells of the underlying lattice type being considered.
(This is, of course, equivalent to imposing some radial 'cutoff' on an interatomic potential.)


EPIs for the Bragg-Williams Hamiltonian can be obtained using a variety of methods, generally those based on density functional theory.
Similarly to the CE method, it is naturally possible to fit EPIs for a given alloy composition to a set of DFT total energy evaluations on alloy supercells with success `[@liumachinenodate; @zhangrobust2020; @liumonte2021]`.
However, most frequently, such EPIs are obtained using the Korringa--Kohn--Rostoker (KKR) formulation of density functional theory `[@korringacalculation1947; @kohnsolution1954; @ebertcalculating2011]`, where the coherent potential approximation (CPA) can be used to describe the average electronic structure and consequent internal energy of the disordered alloy `[@sovencoherentpotential1967; @gyorffycoherentpotential1972; @stockscomplete1978]`.
There are then a variety of suitable techniques available for assessing the energetic cost of applied, inhomogeneous chemical perturbations to the CPA reference medium which naturally lead to extraction of EPIs.
Such techniques include both the generalised perturbation method (GPM) `[@ducastellegeneralized1976; @rubanatomic2004]`, as well as techniques using the language of concentration waves to describe the atomic-scale chemical fluctuations `[@khachaturyanordering1978; @gyorffyconcentration1983]`.
Approaches based on concentration waves have been derived for alloys both in the binary `[@stauntoncompositional1994; @johnsonfirst-principles1994]` and multicomponent `[@singhatomic2015; @khanstatistical2016; @woodgatemodelling2024]` settings.
Once the EPIs for a given alloy composition are obtained, the phase stability of a particular alloy can be examined using sampling techniques applied to the Bragg-Williams model.
This is the purpose of `BraWl` as presented in this work.

## The purpose of `BraWl`

There are a range of existing packages capable of simulating alloy phase equilibria, both open- and closed-source.
Examples of widely-used such packages include ATAT `[@vandewallealloy2002]`, ICET `[@angqvisticet2019]` and CELL `[@rigamonticell2024]`, though all of these focus primarily on implementation of a general cluster expansion, rather than the simpler form of the Bragg-Williams Hamiltonian.
To our knowledge, there is no open-source package specifically focussing on the implementation of a range of sampling algorithms applied to the Bragg-Williams model.
We therefore believe that `BraWl` fills a gap in the capabilities of the current alloy software ecosystem.
Additionally, we hope that the modular way in which the package is constructed could enable implementation of more complex Hamiltonians, as well as further sampling algorithms in addition to those detailed below, in due course.

# Sampling algorithms

`BraWl` implements a range of conventional and enhanced sampling algorithms for exploration of the alloy configuration space.
At present, these are the Metropolis-Hastings algorithm, Wang-Landau sampling, and Nested Sampling.
We briefly outline the details of each of these algorithms below.

## Metropolis-Hastings Monte Carlo

The Metropolis--Hastings algorithm is a useful method for obtaining the equilibrium state of a system.
It achieves this by allowing for a system of interest to follow a chain of states which evolve to, and sample, an equilibrium ensemble `[@metropolisequation1953; @landauguide2014]`.
For the Bragg-Williams model within the canonical ensemble, the algorithm functions by proposing a position swap between two randomly selected atoms in the simulation cell, and calculating the change in energy, $`\Delta E`$, induced by the swap.
The swap is accepted with a probability given by
\begin{equation}
    P_{n\rightarrow m} = 
    \begin{cases}
        \; \text{exp}\left(-\Delta E/k_BT\right), &\quad \Delta E > 0\\
        \; 1, &\quad \Delta E \leq 0,
    \end{cases}
\label{eq:metropolistransitionprobability}
\end{equation}
where $`n`$ labels the initial state, $`m`$ labels the proposed (swapped) state, $`T`$ denotes the simulation temperature, and $`k_B`$ is the usual Boltzmann constant.
These atom swaps can be performed according to Kawasaki dynamics `[@kawasakidiffusion1966]` (nearest neighbour swaps only) or performed between any two atoms in the system.
(The latter option, while less physical, typically allows the system to reach equilibrium in fewer trial Monte Carlo moves.)
Given enough trial moves, the system will reach equilibrium for a fixed simulation temperature.
Once at equilibrium, decorrelated samples can be drawn to obtain thermodynamic averages of various quantities.
Within `BraWl`, this functionality can be used to determine a range of quantities of interest, including atomic short- and long-range order parameters, as well as the simulation heat capacity.
(Definitions of these quantities are provided in Sec. \autoref{sec:physicalquantities}.)
These quantities can be plotted as a function of temperature by considering simulations performed across a range of temperatures.
It is also possible to perform 'simulated annealing' where, starting at high temperature the system is equilibrated at a given temperature and statistics drawn, before the simulation temperature is decreased and the cycle repeated until a desired target temperature is reached.
It also allows for determining the lowest available energy state of the alloy being studied to parametrise a Wang-Landau sampling, discussed below.
Finally, for a given simulation temperature, it is possible to draw decorrelated atomic configurations which can then be used for visualisation and as inputs to other modelling techniques.

## Wang-Landau sampling

Wang-Landau sampling is a flat histogram method which provides a means for high throughput calculation of phase diagrams for atomistic/lattice model systems `[@wangefficient2001]`.
The method allows for direct computation of an estimate of the density of states in energy $`g(E_i)`$, and hence the partition function 
\begin{equation*}
Z = \sum_i g(E_i) e^{-\beta E_i},
\end{equation*}
where the index $`i`$ runs over the appropriately discretised energy macrostates of a given Hamiltonian.
Thermodynamic quantities at any temperature of interest can then be obtained provided one has prior knowledge of the minimum and maximum energy relevant to those temperatures.
The method achieves this by starting from an initial 'guess' of the density of states, which is used to generate a Markov chain of configurations which, as the guess is iteratively refined, tends toward a uniform sampling over energy.
The uniformity is quantified by a  'flatness' criterion on a histogram of visited energies.
The algorithm starts with making initially large modifications to the estimated density of states.
Once the flatness criterion is achieved, the modification factor is reduced and sampling begins again, a process which is repeated over multiple iterations until a desired tolerance on the modification factor is achieved.

Within the Bragg-Williams model, Wang-Landau sampling performs atom swaps as in the Metropolis--Hastings method but with the following acceptance criterion
\begin{equation}
    P_{n\rightarrow m} = 
    \begin{cases}
        \; \frac{g(E_n)}{g(E_m)}, &\quad g(E_n) < g(E_m)\\
        \; 1, &\quad g(E_n) \geq g(E_m),
    \end{cases}
\label{eq:wl_transition_probability}
\end{equation}
where $`g`$ is the density of states, $`E_n`$ is the initial energy and $`E_m`$ is the energy associated with the configuration where the proposed swap has been made.
After each proposed swap, the density of states is updated according to
\begin{equation}
    g(E_i) \rightarrow g(E_i)f_k,
\end{equation}
where $`E_i`$ is the energy of the resultant state, $`f_k`$ is a modification factor initially ($`k=0`$) greater than 1, and $`k`$ is the current iteration index of the Wang-Landau sampling algorithm.
A histogram of the energies visited is maintained, $`H(E)`$, as is a measure of the 'flatness', $`F`$ of the histogram,
\begin{equation}
    F = \frac{\min(H(E))}{\frac{1}{N_b}\sum_i^N H(E_i)},
\end{equation}
where $`N_b`$ denotes the number of bins used in the histogram.
Once $`F`$ is above a given tolerance, sampling is interrupted and $`f`$ is reduced for the next sampling iteration, _e.g._ $`f_{k+1}=\sqrt{f_k}`$.
The visit histogram $`H(E)`$ is set to zero and a new Wang-Landau iteration begins. This is repeated until $`f`$ falls within some desired tolerance, _i.e._ sufficiently close to unity.

Within `BraWl`, the Wang-Landau sampling algorithm features a selection of parallelisation schemes using the message passing interface (MPI) as well as sampling enhancements.
These schemes include energy domain decomposition with dynamic domain sizing and multiple random walkers per domain as well as replica exchange.
This portion of `BraWl` is intended to be used to compute the simulation density of states (in energy) for a given alloy, from which a variety of data can be obtained such as energy distribution histograms at a given temperature as well as the specific heat capacity across a desired temperature range.

## Nested Sampling

Nested sampling (NS) is powerful Bayesian inference technique `[@skillingnested2004; @skillingnested2006; @ashtonnested2022]` adapted to sample the potential energy surface of atomistic systems, giving direct access to the partition function at arbitrary temperatures for comprehensive thermodynamic analysis, without relying on advance knowledge of relevant structures or the range of energies accessible to them `[@partayefficient2010; @partaynested2021]`.


Nested sampling is a top-down iterative approach, starting the sampling with random high-energy configurations and progressing towards the global minimum-energy structure through a series of consecutive nested energy levels. 
Within the Bragg-Williams model, these random high-energy configurations are configurations where the desired ratio of atomic species are assigned randomly to the available lattice sites.
When the sampling is initialised, an integer number, $`K`$, of random alloy configurations are produced - these are often referred to as _walkers_ or the _live set_ - with $`K`$ controlling the resolution (and also the computational cost) of the sampling.
During the iterative process, the configuration with the highest energy is recorded then removed from the live set, and substituted with a new configuration that has a lower energy, but uniformly randomly selected.
The uniform distribution of the $`K`$ walkers allows the estimation of the relative phase space volume, $`w_i= \{[K/(K+1)]^i-[K/(K+1)]^{i+1}\}`$, at iteration $`i`$, and thus the partition function can be evaluated simply during a post-processing step at any arbitrary $`\beta`$, as `[@partaynested2021]`
\begin{equation}
    Z = \sum_{i} w_i e^{-\beta E\left( \left\{ \xi_{i \gamma} \right\} \right)},
\end{equation}
where $`E( \{ \xi_{i \gamma}\})`$ is the energy of the configuration discarded in the $`i`$th NS iteration.
Thermodynamic quantities and weighted average of observables can be evaluated from this.
Since a simple rejection sampling quickly becomes unaffordable as lower energy regions need to be sampled, the new configurations are generated by a random walk in practice: one of the existing $`K`$ samples is selected randomly, cloned, and a series of species swap moves are proposed between randomly selected lattice sites, accepting every swap unless it would cause the energy to increase above the limit.
However, as swap moves cannot be adjusted as the sampling progresses (unlike, _e.g._, atomic displacement steps), the acceptance ratio of swaps decreases.
Thus, to ensure that new configurations are different from the starting structure, the number of proposed particle swap moves is doubled each time the acceptance probability (_i.e._ the ratio of the number of atoms moved during a set of swaps compared to the total number of atoms in the simulation cell) falls below 5\%.


In this work, since the lattice sites are fixed and the Hamiltonian is discretised, it is possible to create multiple different configurations with numerically the same energy value.
However, as the NS algorithm must be able to select the unique highest energy configuration during the iterative sampling process, we have to avoid such degeneracy.
Thus, the energy of each configuration is perturbed by a positive, uniform random number of value less than $`10^{-8}`$ Ry, making each energy value numerically distinct without effecting the uniform distribution of samples or any thermodynamic properties.

# Physical quantities of interest

`BraWl` can extract a range of quantities of interest from a given alloy simulation.
The relevant quantities are as follows.

## Internal energy

For a given lattice type, system size, alloy composition, and set of atom-atom effective pair interactions, `BraWl` can evaluate the total energy associated with the alloy configuration (Eq. \autoref{eq:b-w1}).
At the time of writing, for speed, common lattice types (fcc, bcc, simple cubic, _etc._) are hard-coded, with the intention that the range of implemented lattice types will be expanded over time as necessary.
Where only the relative _change_ of energy induced by swapping/substituting atoms is considered, `BraWl` takes advantage of the mathematical form of the Bragg-Williams Hamiltonian to only evaluate the relevant terms in the summation in Eq. \autoref{eq:b-w1} which are changed as a result of the swap/substitution.
This substantially reduces the cost of evaluation of change in simulation energy induced by a particular atomic swap.

## Heat capacity

The isochoric (fixed volume) heat capacity at a given temperature, $`C_V(T)`$, is a useful quantity for identifying phase transitions, as a plot of the simulation heat capacity as a function of temperature is expected to show a local peak at the temperature at which the transition occurs.
Within `BraWl`, the heat capacity is calculated via
\begin{equation}
    C_V(T) = \frac{\langle E^2\rangle - \langle E\rangle^2}{k_BT^2},
\end{equation}
where $`k_B`$ is the usual Boltzmann constant, $`E`$ is the simulation energy, and $`\langle \cdot \rangle`$ denote thermodynamic averages obtained using the relevant sampling algorithm.

## Atomic short-range order (ASRO) parameters

To assess local atom-atom correlations in a simulation, `BraWl` can calculate the Warren-Cowley atomic short-range order parameters `[@cowleyapproximate1950; @cowleyshort-range1965]`, adapted to the multicomponent setting, defined as
\begin{equation}
    \alpha^{\gamma \gamma^{\prime}}_n=1-\frac{P^{\gamma \gamma^{\prime}}_n}{c_{\gamma^{\prime}}}
\end{equation}
where $`n`$ refers to the $`n`$th coordination shell, $`P_n^{\gamma \gamma^{\prime}}`$ is the conditional probability of an atom of type $`q`$ neighbouring an atom of type $`p`$ on shell $`n`$, and $`c_q`$ is the total concentration of atom type $`q`$. When $`\alpha_n^{\gamma \gamma^{\prime}} > 0`$, $`p`$-$`q`$ pairs are disfavoured on shell $`n`$ and, when $`\alpha_n^{\gamma \gamma^{\prime}} < 0`$ they are favoured. The value $`\alpha_n^{\gamma \gamma^{\prime}} = 0`$ corresponds to the ideal, maximally disordered solid solution.

For maximal flexibility, `BraWl` outputs the conditional probabilities, $`P_n^{\gamma \gamma^{\prime}}`$, and the user can then choose whether or not to rescale and obtain the Warren-Cowley ASRO parameters as a post-processing step. The `BraWl` package can calculate these parameters averaged across a single configuration, or averaged at a particular temperature using one of the sampling algorithms implemented in the code.

## Atomic long-range order (ALRO) parameters

Over a simulation run (for example using the Metropolis algorithm), `BraWl` can calculate the average partial occupancies of a lattice site, $`\langle \xi_{i\gamma}\rangle = c_{i\gamma}`$, which are natural atomic long-range order parameters describing chemically ordered phases.
This capability was first demonstrated on simulations of Fe-Ga alloys `[@marchantab2021]`.

# Example Applications

`BraWL` has been used, with success, to study the phase behaviour of a range of binary and multicomponent alloys, for example the binary Fe-Ga system (Galfenol) `[@marchantab2021]`, the Fe-Ni system `[@woodgateintegrated2024]`, the Cantor-Wu medium- and high-entropy alloys `[@woodgatecompositional2022; woodgateinterplay2023]`, the refractory high-entropy alloys `[@woodgateshortrange2023; woodgatecompetition2024]`, the Al$`_x`$CrFeCoNi system `[@woodgatestructure2024]`, and the AlTiVNb and AlTiCrMo refractory high-entropy superalloys `[@woodgateemergentnodate]`.
The package has also been used to generate atomic configurations for training datasets for machine-learned interatomic potentials, for example for the prototypical austenitic stainless steel, Fe$`_7`$Cr$`_2`$Ni `[@shenoycollinearspin2024]`.

![Evolution of the simulation internal energy (top panel) and conditional pair probabilities (bottom panel) for an Fe$`_{0.5}`$Ni$`_{0.5}`$ alloy as a function of the number of Metropolis-Hastings sweeps at a simulation temperature of $`T=300`$ K. One 'sweep' is one trial move per atom in the system. Beyond approximately 100 sweeps, the system can be seen to have reached equilibrium, with L1$`_0`$ order established.\label{fig:feniequilibration}](feniequilibration.pdf){ width=40% }

As an example of the Metropolis-Hastings Monte Carlo algorithm, we consider its application to the binary FeNi alloy, first discussed by `@woodgateintegrated2024`.
Fig. \autoref{fig:feniequilibration} shows the internal energy and conditional pair probabilities (quantifying ASRO) of a simulation cell containing 256 atoms as a function of the number of Metropolis-Hastings 'sweeps', where a sweep refers to performing a number of trial Metropolis-Hastings moves equal to the number of atoms in the simulation cell.
The simulation is performed at 300 K, below the alloy's L1$_0$ disorder-order transition temperature.
The L1$_0$ phase is a structure where $`2/3`$ of the nearest neighbours of Fe (Ni) atoms are Ni (Fe) atoms, and where none of the next-nearest neighbours of Fe (Ni) atoms are Ni (Fe) atoms.
It can be seen that this ordering is swiftly established as the number of Monte Carlo sweeps increases, albeit with some remaining thermal noise.

![Plots of energy probability distributions, Warren-Cowley ASRO parameters ($`\alpha_n^{pq}`$) and simulation heat capacity ($C$) as a function of temperature for AlTiCrMo obtained using lattice-based Monte Carlo simulations employing Wang-Landau sampling. Here, show $`\alpha_n^{pq}`$ only for $`n = 1`$. The zero of the energy scale for the energy histograms is set to be equal to the average internal energy of the alloy obtained at a simulation temperature of 3000 K.\label{fig:wlAlTiCrMo}](wlAlTiCrMo.pdf){ width=40% }

As an example of Wang-Landau sampling, we consider its application to the AlTiCrMo refractory high-entropy superalloy, first discussed by `@woodgateemergentnodate`, for which results are shown in Fig. \autoref{fig:wlAlTiCrMo}.
The top panel shows calculated energy probability distributions (histograms) at various temperatures, while the bottom panel shows the simulation heat capacity and evolution of the Warren-Cowley ASRO parameters as a function of temperature.
The high-temperature peak in the heat capacity data is associated with the experimentally observed B2 crystallographic ordering.

![Internal energy, $`E`$, and isochoric heat capacity, $`C_V`$, obtained using the Nested Sampling algorithm applied to the equiatomic, fcc, AlCrFeCoNi high-entropy alloy. The simulation cell contained 108 atoms. Upon cooling, the initial peak in the heat capacity is associated with an L1$_2$ ordering driven by Al, with subsequent peaks indicating eventual decomposition into multiple competing phases.\label{fig:nsalcrfeconi}](nsalcrfeconi.pdf){ width=40% }

Finally, as an example of application of the Nested Sampling algorithm, we consider its application to the AlCrFeCoNi high-entropy alloy, first discussed in `@woodgatestructure2024`.
Fig. \autoref{fig:nsalcrfeconi} plots the internal energy, $`E`$, and isochoric heat capacity, $`C_V`$, obtained for the equiatomic, fcc, AlCrFeCoNi system.
The simulation cell contained 108 atoms.
The initial peak in the heat capacity encountered upon cooling from high temperature is associated with an L1$_2$ ordering driven by Al, with subsequent peaks indicating eventual decomposition into multiple competing phases, which is understood to be consistent with experimental data for this system.

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

