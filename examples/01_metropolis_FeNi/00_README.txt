# Example 1: Metropolis-Hastings Monte Carlo.

These examples demonstrate applications of Metropolis-Hastings Monte Carlo to the FeNi system. 

Atom-atom interactions are taken from:
C. D. Woodgate, L. H. Lewis, J. B. Staunton, [npj Comput. Mater. 10, 272 (2024)](https://doi.org/10.1038/s41524-024-01435-y).

These so-called 'effective pair interactions' are material-specific and must be chosen/derived suitably for the system being studied. It is common to recover them by fitting to an _ab initio_ DFT dataset, or from some form of perturbative concentration wave analysis applied to a given reference medium, typically that described by the coherent potential approximation (CPA). **It is therefore important to realise that `BraWl` is intended to be used as part of a multistage modelling workflow, rather than simply as a standalone package.

The first example, `01_equilibrate`, demonstrates how to use the Metropolis-Hasting algorithm to inspect equilibration of a simulation system at a given temperature. Navigate to that directory, then run `path/to/BraWl/brawl.run` (using `mpirun -np 1` if a parallel build). The simulation should run pretty quickly and then a range of new files/directories containing simulation data will appear. You can then analyse the simulation trajectory using the provided Python script, `01_visualise_equilibration.py`. You should see that the simulation reaches equiligrium fairly quickly, usually after around 100 Metropolis-Hastings 'sweeps' for a system of this size.

The second example, `02_simulated_annealing`, uses the Metropolis-Hastings algorithm and simulated annealing to inspect the equilibrium ASRO and simulation heat capacity as a function of temperature. Running BraWl in this directory takes a little longer (at each temperature we draw statistics over 256e4 trial 'moves'), but the simulation should finish in half an hour or so, usually less. Again, a Python script (`01_plot_results.py`) lets you see the simulation output. The data may be a little noisy (we have probably slightly under-sampled), but you should see the known $\textrm{L}1_0$ order establish itself as the simulation is cooled from high temperature, indicated by a peak in the heat capacity data and associated change in the Warren-Cowley ASRO parameters.

N.B. The 'acceptance' rate counts a trial swap of two atoms of the same chemical species as an 'accepted' move, so at low temperatures for a 50:50 binary alloy, this rate tends to 0.5. (This choice ensures that the acceptance rate tends to 1 as temperature tends to infinity.)
