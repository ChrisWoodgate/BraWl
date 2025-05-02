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
  use io
  use kinds
  use c_functions
  use write_netcdf
  use write_xyz
  use metropolis_output
  use display
  use tests

  implicit none

  ! Runtime parameters type
  type(run_params) :: setup

  ! Start MPI
  call comms_initialise()

  ! Print software info to the screen
  if(my_rank == 0) call write_info('s')

  ! Print that this is a test run
  call print_centered_message('Test cases', '-')

  ! We will test both fcc and bcc implementations
  call print_centered_message('Testing FCC example', '-')

  setup%n_1 = 4
  setup%n_2 = 4
  setup%n_3 = 4
  setup%n_basis = 1
  setup%n_species=3
  setup%lattice='fcc'
  allocate(setup%species_concentrations(3))
  setup%species_concentrations=(/0.333, 0.333, 0.334/)
  setup%interaction_file = 'fcc.vij'
  setup%interaction_range = 4
  setup%static_seed = .True.
  setup%wc_range = 3

  call test_suite(setup, my_rank, 256, 'generate')

  ! Print software info to the screen
  if(my_rank == 0) call write_info('f')

  ! Finalise MPI
  call comms_finalise()

end program test
