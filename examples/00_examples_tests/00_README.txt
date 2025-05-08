# Example 0: Examples and tests

*N.B.* _This is intended for developers of the code to understand core functionality. If you are simply an end-user, you can skip this example._

This example demonstrates how to run the `examples` main program, which calls routines kept in `src/howto_examples.f90`.

To build the `examples` main, use the option `make compiler=<compiler> example`

The file `brawl.inp` contains all relevant parameters for setting up the simulation, while 'epis.vij' contains the atom-atom interaction parameters. Then `metropolis.inp` contains some details of parameters relating to Metropolis-Hastings Monte Carlo which will be demonstrated.

Then run `mpirun -n 1 /path/to/BraWL/example.run` (in parallel) or `/path/to/brawl/example.run` (in serial) from this directory to launch the simulation. There will be graphical output to explain what happens. 🙂
