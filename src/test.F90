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

  ! Flag for passing tests for each lattice type
  logical :: fcc_pass, bcc_pass

  ! Start MPI
  call comms_initialise()

  ! Print software info to the screen
  if(my_rank == 0) call write_info('s')

  ! If we try to run the test suite in parallel, this is not currently handled
  if (p .gt. 1) then
    call comms_finalise()
    stop 'Please run the test suite on one processor: mpirun -np 1 /path/to/brawl/tests.run'
  end if

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

  fcc_pass = test_suite(setup, my_rank, 256, 'fcc', 'test')

  ! If we pass, print to the screen
  if (fcc_pass) then
    print*, 'Passed fcc lattice tests', new_line('a')
  end if

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

  bcc_pass = test_suite(setup, my_rank, 128, 'bcc', 'test')

  ! If we pass, print to the screen
  if (bcc_pass) then
    print*, 'Passed bcc lattice tests', new_line('a')
  end if

  deallocate(setup%species_concentrations,setup%species_numbers,setup%species_names)

  ! If all tests are passed, print to the screen
  if ((fcc_pass) .and. (bcc_pass)) then
    print*, 'Passed all tests successfully!'
  ! Otherwise, tell the user that there is a problem
  else
    print*, 'Some (or all) tests failed :-('
    print*, 'Check the output above to see what the WARNINGs were'
  end if

  ! Print software info to the screen
  if(my_rank == 0) call write_info('f')

  ! Finalise MPI
  call comms_finalise()

end program test
