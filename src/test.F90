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
  use netcdf_io
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
  call print_centered_message('Test cases', '-', .True.)

  ! We will test both fcc and bcc implementations
  write(6,'(72("-"))')
  call print_centered_message('Testing FCC example', '-')
  write(6,'(72("-"),/)')

  !-----------------------!
  ! Tests for fcc lattice !
  !-----------------------!
  setup%n_1 = 4
  setup%n_2 = 4
  setup%n_3 = 4
  setup%n_basis = 1
  setup%n_species=5
  setup%lattice='fcc'
  allocate(setup%species_concentrations(0:5))
  setup%species_concentrations = 0.0_real64
  allocate(setup%species_numbers(5))
  setup%species_numbers = 0
  allocate(setup%species_names(5))
  setup%species_names = (/ 'Al', 'Cr', 'Fe', 'Co', 'Ni' /)
  ! The '0th' value of species_concentrations should be zero---used in some array initialisations
  setup%species_concentrations=(/0.0_real64, 0.2_real64, 0.2_real64, 0.2_real64, 0.2_real64, 0.2_real64/)
  setup%interaction_file = 'fcc_epi.vij'
  setup%interaction_range = 4
  setup%static_seed = .True.
  setup%wc_range = 3

  call test_suite(setup, my_rank, 256, 'fcc', 'test')

  deallocate(setup%species_concentrations,setup%species_numbers,setup%species_names)

  ! We will test both fcc and bcc implementations
  write(6,'(72("-"))')
  call print_centered_message('Testing BCC example', '-')
  write(6,'(72("-"),/)')

  !-----------------------!
  ! Tests for bcc lattice !
  !-----------------------!
  setup%n_1 = 4
  setup%n_2 = 4
  setup%n_3 = 4
  setup%n_basis = 1
  setup%n_species=4
  setup%lattice='bcc'
  allocate(setup%species_concentrations(0:4))
  setup%species_concentrations = 0.0_real64
  allocate(setup%species_numbers(4))
  setup%species_numbers = 0
  allocate(setup%species_names(4))
  setup%species_names = (/ 'Al', 'Ti', ' V', 'Nb'/)
  ! The '0th' value of species_concentrations should be zero---used in some array initialisations
  setup%species_concentrations=(/0.0_real64, 0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64/)
  setup%interaction_file = 'bcc_epi.vij'
  setup%interaction_range = 6
  setup%static_seed = .True.
  setup%wc_range = 3

  call test_suite(setup, my_rank, 128, 'bcc', 'test')

  ! Print software info to the screen
  if(my_rank == 0) call write_info('f')

  ! Finalise MPI
  call comms_finalise()

end program test
