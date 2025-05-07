!> @file    tests.F90
!>
!> @brief   Routines for testing core functionality
!>
!> @details This module contains routines which run tests of core
!>          functionality
!>
!> @author  C. D. Woodgate
!>
!> @date    2025
module tests

  use initialise
  use kinds
  use shared_data
  use io
  use display, only : print_centered_message
  use c_functions
  use bw_hamiltonian
  use random_site
  use analytics
  use netcdf_io
  use write_xyz
  use metropolis_output
  use metropolis, only : monte_carlo_step_lattice
  
  implicit none

  private

  public :: test_suite

  contains

  !> @brief   Subroutine for testing core functionality of code
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  my_rank Rank of this process
  !> @param  n_steps Number of trial Metropolis-Hastings moves to run
  !> @param  lattice Lattice type to test
  !> @param  mode If 'generate', will generate data. If 'test', will test.
  !>
  !> @return True if all tests pass, False if not
  function test_suite(setup, my_rank, n_steps, lattice, mode) result(passed)

    ! Derived type describing simulation setup
    type(run_params) :: setup

    ! Rank of this processor
    integer, intent(in) :: my_rank

    ! How many trial MC steps to run
    integer, intent(in)  :: n_steps

    ! Lattice type to test
    character(len=*)  :: lattice

    ! Mode, whether to generate data or test data
    character(len=*)  :: mode

    ! Whether or not we pass the test
    logical :: passed

    ! Integers used in calculations
    integer :: i, accept, ierr

    ! Integers to count any disagreements
    integer :: energy_disagreements, asro_disagreements

    ! Temperature and temperature steps
    real(real64) :: beta, temp, sim_temp

    ! Array for loading configuration arrays and cross-checking
    integer(array_int), allocatable, dimension(:,:,:,:) :: test_config

    ! Arrays for checking energies of trajectories
    real(real64), dimension(:), allocatable :: trajectory_energy,      &
                                               trajectory_energy_test, &
                                               indices

    ! Arrays for checking energies of trajectories
    real(real64), dimension(:,:,:,:), allocatable :: trajectory_asro, trajectory_asro_test

    ! Start by assuming we have passed
    ! Will switch to .False. if a test fails
    passed=.True.

    ! Allocate the config array
    if (allocated(config)) deallocate(config)
    allocate(config(setup%n_basis, 2*setup%n_1, 2*setup%n_2, 2*setup%n_3))

    ! Allocate the array for lattice shells
    if (allocated(shells)) deallocate(shells)
    if (.not. allocated(shells)) allocate(shells(setup%wc_range))

    ! Initialise some function pointers
    call initialise_function_pointers(setup)

    ! Allocate space for atom-atom interaction parameters
    call initialise_interaction(setup)

    ! Read in atom-atom interaction
    call read_exchange(setup, my_rank)

    ! Initialise PNRG
    ! static_seed is true if we would like to use a fixed seed and false
    ! otherwise
    call initialise_prng(setup%static_seed)

    ! Print to the screen that we have successfully initialised the lattice
    print*, ' '
    call print_centered_message('Initialising lattice', '-', .True.)

    ! Do an initial setup of the lattice
    call initial_setup(setup, config)

    ! Get the coordination shells of the lattice
    call lattice_shells(setup, shells, config)

    !------------------------------!
    ! Lattice initialisation check !
    !------------------------------!

    ! Now either read in the reference data (if in 'test' mode) or generate
    ! this data, if in 'generate' mode.
    if (trim(mode) .eq. 'test') then
      call ncdf_config_reader('99_ref/'//lattice//'_start_config.nc', test_config, setup)
    else if (trim(mode) .eq. 'generate') then
      call ncdf_grid_state_writer('99_ref/'//lattice//'_start_config.nc', ierr, config, setup)
      test_config=config
    end if

    ! Print to the screen if we do
    call print_centered_message('Checking initial random '//lattice//' configurations', '-', .True.)

    ! Print to the screen whether or not we agree
    if (.not.(configs_equal(config, test_config))) then
      print*, 'WARNING: Initial random '//lattice//' configuration generated different from reference'
      passed=.False.
    else
      print*, 'Initial random '//lattice//' configurations (fixed seed) are identical!', new_line('a')
    end if

    !---------------------------!
    ! Metropolis-Hastings sweep !
    !---------------------------!

    ! Now start a Metropolis-Hastings sweep
    call print_centered_message('Starting trial Metropolis-Hastings sweep', '-', .True.)

    ! And allocate some memory for storing output data
    allocate(trajectory_energy(n_steps))
    allocate(indices(n_steps))
    allocate(trajectory_asro(setup%n_species, setup%n_species, setup%wc_range, n_steps))

    ! Set the simulation temperature
    temp=300.0
    sim_temp = temp*k_b_in_Ry
    beta = 1.0_real64/sim_temp

    ! Do some Monte Carlo moves, storing energies
    print*, 'Will run ', n_steps, ' trial Metropolis-Hastings moves', new_line('a')
    do i=1, n_steps
      accept = monte_carlo_step_lattice(setup, config, beta)
      trajectory_energy(i) = setup%full_energy(config)
      trajectory_asro(:,:,:,i) = radial_densities(setup, config, setup%wc_range, shells)
      indices(i) = real(i, kind=real64)
    end do
    call print_centered_message('Completed trial Metropolis-Hastings sweep', '-', .True.)

    !------------------------------!
    ! Energy/ASRO trajectory check !
    !------------------------------!

    ! Now either read in the reference data (if in 'test' mode) or generate
    ! this data, if in 'generate' mode.
    if (trim(mode) .eq. 'test') then
      call ncdf_radial_density_reader('99_ref/'//lattice//'_data.nc',  &
                                      trajectory_asro_test,            &
                                      trajectory_energy_test,          &
                                      setup,                           &
                                      n_steps)
    else if (trim(mode) .eq. 'generate') then
      call ncdf_radial_density_writer('99_ref/'//lattice//'_data.nc',  &
                                      trajectory_asro,                 &
                                      shells,                          &
                                      indices,                         &
                                      trajectory_energy,               &
                                      setup)

      trajectory_energy_test = trajectory_energy
      trajectory_asro_test = trajectory_asro
    end if

    call print_centered_message('Checking '//lattice//' energy trajectory', '-', .True.)
    energy_disagreements = array_equal_1D(trajectory_energy,           &
                                          trajectory_energy_test)
    if (energy_disagreements .gt. 0) then
      print*, 'WARNING: There are ', energy_disagreements,             &
            ' outside tolerance in calculated '//lattice//' energies', &
            new_line('a')
      passed = .False.
    else
      print*, 'Energy trajectory for '//lattice//' lattice agrees with reference', &
               new_line('a')
    end if

    call print_centered_message('Checking '//lattice//' ASRO trajectory', '-', .True.)
    asro_disagreements = array_equal_4D(trajectory_asro, trajectory_asro_test)
    if (asro_disagreements .gt. 0) then
      print*, 'WARNING: there are ', asro_disagreements,               &
              ' outside tolerance in calculated '//lattice//' asro',   &
              new_line('a')
      passed=.False.
    else
      print*, 'ASRO trajectory for '//lattice//' lattice agrees with reference', new_line('a')
    end if

    !---------------------------!
    ! Final configuration check !
    !---------------------------!

    ! Now either read in the reference data (if in 'test' mode) or generate
    ! this data, if in 'generate' mode.
    if (trim(mode) .eq. 'test') then
      call ncdf_config_reader('99_ref/'//lattice//'_end_config.nc', test_config, setup)
    else if (trim(mode) .eq. 'generate') then
      call ncdf_grid_state_writer('99_ref/'//lattice//'_end_config.nc', ierr, config, setup)
      test_config=config
    end if

    call print_centered_message('Checking final '//lattice//' configurations', '-', .True.)

    if (.not.(configs_equal(config, test_config))) then
      print*, 'WARNING: Final '//lattice//' configurations are different'
      passed = .False.
    else
      print*, 'Final '//lattice//' configurations are identical!', new_line('a')
    end if

    !----------------------------------!
    ! Clean up at the end of the tests !
    !----------------------------------!

    call clean_up_interaction()
    if(allocated(config)) deallocate(config)
    if(allocated(test_config)) deallocate(test_config)
    if(allocated(shells)) deallocate(shells)
    if(allocated(indices)) deallocate(indices)
    if(allocated(trajectory_energy)) deallocate(trajectory_energy)
    if(allocated(trajectory_asro)) deallocate(trajectory_asro)
    if(allocated(trajectory_energy_test)) deallocate(trajectory_energy_test)
    if(allocated(trajectory_asro_test)) deallocate(trajectory_asro_test)

  end function test_suite

  !> @brief   Function for comparing two 4D arrays
  !>
  !> @author  C. D. Woodgate
  !>
  !> @date    2025
  !>
  !> @param  array1 First array
  !> @param  array2 Second array
  !> @param  tolerance Tolerance. Defaults to epsilon(real32), i.e. single-bit precision
  !>
  !> @return The number of indices where the array elements are not
  !>         within the specified tolerance.
  function array_equal_4D(array1, array2, tolerance) result(disagreements)

    real(real64), dimension(:,:,:,:), allocatable, intent(in) :: array1, array2
    real(real64), optional :: tolerance
    real(real32) :: single
    real(real64) :: tol
    integer, dimension(4) :: shape1, shape2
    integer :: disagreements
    integer :: i, j, k, l

    if (.not.present(tolerance)) then
      tol = real(epsilon(single), kind=real64)
    else
      tol = tolerance
    end if

    shape1 = shape(array1)
    shape2 = shape(array2)

    ! Check that two configs are of the same shape
    if (.not. all(shape1 .eq. shape2)) then
      disagreements = -1
      print*, 'Warning: input arrays not of same shape in array_equal_4D'
      return
    end if

    disagreements = 0

    ! Compare the two arrays element by element
    do l=1, shape1(4)
      do k=1, shape1(3)
        do j=1, shape1(2)
          do i=1, shape1(1)
            if (abs(array1(i,j,k,l)-array2(i,j,k,l)) .gt. tol) then
              disagreements = disagreements + 1
            end if
          end do
        end do
      end do
    end do

  end function array_equal_4D

  !> @brief   Function for comparing two 1D arrays
  !>
  !> @author  C. D. Woodgate
  !>
  !> @date    2025
  !>
  !> @param  array1 First array
  !> @param  array2 Second array
  !> @param  tolerance Tolerance. Defaults to epsilon(real32), i.e. single-bit precision
  !>
  !> @return The number of indices where the array elements are not
  !>         within the specified tolerance.
  function array_equal_1D(array1, array2, tolerance) result(disagreements)

    real(real64), dimension(:), allocatable, intent(in) :: array1, array2
    real(real64), optional :: tolerance
    real(real32) :: single
    real(real64) :: tol
    integer :: len1, len2
    integer :: disagreements
    integer :: i

    if (.not.present(tolerance)) then
      tol = real(epsilon(single), kind=real64)
    else
      tol = tolerance
    end if

    len1 = size(array1)
    len2 = size(array2)

    ! Check that two configs are of the same shape
    if (.not. (len1 .eq. len2)) then
      disagreements = -1
      print*, 'Warning: input arrays not of same length in array_equal_1D'
      return
    end if

    disagreements = 0

    ! Compare the two arrays element by element
    do i=1, len1
      if (abs(array1(i)-array2(i)) .gt. tol) then
        disagreements = disagreements + 1
      end if
    end do

  end function array_equal_1D

  !> @brief   Function for comparing two configurations
  !>
  !> @author  C. D. Woodgate
  !>
  !> @date    2025
  !>
  !> @param  config1 First configuration
  !> @param  config2 Second configuration
  !>
  !> @return True if configurations are the same, False if not
  function configs_equal(config1, config2) result(equal)

    integer(array_int), dimension(:,:,:,:), allocatable, intent(in) :: config1, config2
    integer, dimension(4) :: shape1, shape2
    logical :: equal
    integer :: i, j, k, l

    equal=.True.

    shape1 = shape(config1)
    shape2 = shape(config2)

    ! Check that two configs are of the same shape
    if (.not. all(shape1 .eq. shape2)) then
      equal = .False.
      return
    end if

    ! Compare the two configurations element by element
    do l=1, shape1(4)
      do k=1, shape1(3)
        do j=1, shape1(2)
          do i=1, shape1(1)
            if (config1(i,j,k,l) .ne. config2(i,j,k,l)) then
              equal = .False.
              return
            end if
          end do
        end do
      end do
    end do

  end function configs_equal

end module tests
