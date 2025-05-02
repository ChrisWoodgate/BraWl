!> @file    tests.F90
!>
!> @brief   Test program for verifying core code functionality against
!>          established test cases with known results.
!>
!> @TODO    Verify:
!>           - Input read-in
!>           - Hamiltonian evaluation
!>           - Energy trajectory for Metropolis run
!>           - Generated configuration at start/end of run
!>
!> @author  C. D. Woodgate
!>
!> @date    2025
program test
  
  use initialise
  use comms
  use shared_data
  use metropolis
  use nested_sampling
  use io
  use kinds
  use c_functions
  use write_netcdf
  use write_xyz
  use metropolis_output
  use display

#ifdef USE_MPI

  use tmmc
  use wang_landau

#endif

  implicit none

  ! Runtime parameters type
  type(run_params) :: setup

  ! Metropolis parameters type
  type(metropolis_params) :: metropolis

  ! Nested Sampling parameters type
  type(ns_params) :: ns_setup

#ifdef USE_MPI

  ! Tmmc parameters type
  type(tmmc_params) :: tmmc_setup

  ! Wang Landau parameters type
  type(wl_params) :: wl_setup

#endif

  ! Start MPI
  call comms_initialise()

  ! Print software info to the screen
  if(my_rank == 0) call write_info('s')

  ! Print that this is a test run
  call print_centered_message('Test cases', '-')

  ! We will test both fcc and bcc implementations
  call print_centered_message('Testing FCC example', '-')

  ! Parse inputs
  ! call parse_inputs(setup, my_rank)

  ! Allocate space for atom-atom interaction parameters
  ! call initialise_interaction(setup)

  ! Read in atom-atom interaction
  ! call read_exchange(setup, my_rank)

  ! Initialise PNRG
  ! static_seed is true if we would like to use a fixed seed and false
  ! otherwise
  ! call initialise_prng(setup%static_seed)

  ! Initialise some function pointers
  ! call initialise_function_pointers(setup)

  ! Initialise some local arrays
  ! call initialise_local_arrays(setup)

  ! Clean up
  ! call clean_up_interaction()

  ! Clean up
  ! call local_clean_up(setup)

  ! Print software info to the screen
  if(my_rank == 0) call write_info('f')

  ! Finalise MPI
  call comms_finalise()

end program test
