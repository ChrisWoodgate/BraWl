Example 3: Wang-Landau sampling of the AlTiCrMo system, in the fcc structure.

This example demonstrates a Wang-Landau simulation of the AlTiCrMo system, in the fcc structure. As a quick demonstration, the input file is set up to use 128 atoms with 4 windows each with 2 walkers. 
Atom-atom interactions are taken from:
`???`

A Wang-Landau sampling run needs the following input files:
 - The file `brawl.inp` contains all relevant parameters for the simulation, while `AlTiCrMo.vij` contains the atom-atom interaction parameters.
 - The file `wl_input.inp` contains the relevant Wang-Landau parameters.

Run `mpirun -n 8 ../../brawl.run` to launch the simulation. Once it has finished, analyse the results with `wl_vis.py`.

The simulation should take around 30 minutes.
