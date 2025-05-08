# Example 3: Wang-Landau sampling of the AlTiCrMo system, in the bcc structure.

This example demonstrates a Wang-Landau simulation of the AlTiCrMo system, in the bcc structure. As a quick demonstration, the input file is set up to use 128 atoms with 4 windows each with 2 walkers. 
Atom-atom interactions are taken from:
C. D. Woodgate, H. J. Naguszewski, D. Redka, J. Minar, D. Quigley, J. B. Staunton, [arXiv:2503.13235](https://arxiv.org/abs/2503.13235).

A Wang-Landau sampling run needs the following input files:
 - The file 'input.inp' contains all relevant parameters for the simulation, while 'AlTiCrMo.vij' contains the atom-atom interaction parameters.
 - The file 'wl\_input.inp' contains the relevant Wang-Landau parameters.

Run 'mpirun -n 8 /path/to/brawl/brawl.run' to launch the simulation. Once it has finished, analyse the results with wl\_vis.py. (Note that this example requires a parallel build.)

The simulation should take around 45 minutes. You may want to run it in the background. A good tool for this (and other terminal multiplexing activities) is GNU Screen.
