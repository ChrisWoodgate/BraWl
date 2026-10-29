!> @file    example.f90
!>
!> @brief   Example program showing the code's functionalities for
!>          developers
!>
!> @details This file contains an example program showcasing the code's
!>          functionalities. It is built with
!>          `make compiler=<compiler> example` and run from the
!>          directory `examples/00_examples_tests`. The sequence of
!>          calls mirrors that of the main program in main.f90: the
!>          global (cross-rank) Metropolis arrays only exist in the
!>          parallel build, so the calls which set them up and tear
!>          them down are guarded by the USE_MPI preprocessor macro.
!>
!> @author  C. D. Woodgate
!> @author  Claude Fable 5.1
!> @date    2019-2026
program example
  
  use initialise
  use comms
  use shared_data
  use metropolis
  use io
  use kinds
  use c_functions
  use netcdf_io
  use write_xyz
  use metropolis_output
  use display
  use howto_examples

  implicit none

  ! Runtime parameters type
  type(run_params) :: setup

  ! Metropolis parameters type
  type(metropolis_params) :: metropolis

  ! Start MPI
  call comms_initialise()

  ! Print software info to the screen
  if(my_rank == 0) call write_info('s')

  ! Parse inputs
  call parse_inputs(setup, my_rank)

  ! Allocate space for atom-atom interaction parameters
  call initialise_interaction(setup)

  ! Read in atom-atom interaction
  call read_exchange(setup, my_rank)

  ! Read the Metropolis control file
  call parse_metropolis_inputs(metropolis, my_rank)

  ! Initialise PNRG
  ! static_seed is true if we would like to use a fixed seed and false
  ! otherwise
  call initialise_prng(setup%static_seed)

  ! Initialise some function pointers
  call initialise_function_pointers(setup)

  ! Make directories for data
  call make_data_directories(my_rank)

  ! Make sure the required directories have been made before any rank
  ! tries to write to them
  call comms_wait()

#ifdef USE_MPI

  ! Initialise some global arrays. These are only used (and only
  ! defined) in the parallel build.
  call initialise_global_metropolis_arrays(setup, metropolis)

#endif

  ! Initialise some local arrays
  call initialise_local_arrays(setup)
  call initialise_local_metropolis_arrays(setup, metropolis)

  !-----------------------------------------!
  ! Main examples routine in howto_examples !
  !-----------------------------------------!
  call examples(setup, metropolis, my_rank)

  ! Clean up
  call local_metropolis_clean_up()

#ifdef USE_MPI

  ! Clean up
  call global_metropolis_clean_up()

#endif

  ! Clean up
  call clean_up_interaction()

  ! Clean up
  call local_clean_up(setup)

  ! Print software info to the screen
  if(my_rank == 0) call write_info('f')

  ! Finalise MPI
  call comms_finalise()

end program example
