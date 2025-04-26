Example 1: Metropolis-Hastings Monte Carlo.

These examples demonstrate applications of Metropolis-Hastings Monte Carlo to the  FeNi system. 

Atom-atom interactions are taken from:
C. D. Woodgate, L. H. Lewis, J. B. Staunton, npj Comput. Mater. 10, 272 (2024).
DOI: https://doi.org/10.1038/s41524-024-01435-y

The first example, 01_equilibrate, demonstrates how to use the Metropolis-Hasting algorithm to inspect equilibration of a simulation system at a given temperature.

Navigate to that directory, then run `path/to/BraWl/brawl.run` (using mpirun if a parallel build). The simulation should run pretty quickly and then a range of new files/directories will appear.

You can analyse the simulation trajectory using the provided Python script, `01_visualise_equilibration.py`

Currently this directory is 'work in progress'. Please check back later for a full example!
