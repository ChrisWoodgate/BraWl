# Example 3: Wang-Landau sampling of the AlTiCrMo system, in the bcc structure.

This example demonstrates a Wang-Landau simulation of the AlTiCrMo system, in the bcc structure. As a quick demonstration, the input file is set up to use 128 atoms with 4 windows each with 2 walkers. 
Atom-atom interactions are taken from:
C. D. Woodgate, H. J. Naguszewski, D. Redka, J. Minar, D. Quigley, J. B. Staunton, [arXiv:2503.13235](https://arxiv.org/abs/2503.13235).

A Wang-Landau sampling run needs the following input files:
 - The file 'input.inp' contains all relevant parameters for the simulation, while 'AlTiCrMo.vij' contains the atom-atom interaction parameters.
 - The file 'wl\_input.inp' contains the relevant Wang-Landau parameters. An optional `dos_tolerance=<percent>` line stops the run early once the RMSE of the change in ln g(E) between WL iterations (uniform offset removed), as a percentage of the ln g(E) range, falls below that value (e.g. `dos_tolerance=0.5`). If omitted (or <= 0) only the `tolerance` on the modification factor ln f is used.

Run 'mpirun -n 8 /path/to/brawl/brawl.run' to launch the simulation. Once it has finished, analyse the results with wl\_vis.py. (Note that this example requires a parallel build.)

The simulation should take at most 5 minutes. You may want to run it in the background. A good tool for this (and other terminal multiplexing activities) is GNU Screen.
