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
  subroutine test_suite(setup, my_rank, n_steps, lattice, mode)

    ! Mode, whether to generate data or test data
    character(len=*)  :: lattice

    ! Mode, whether to generate data or test data
    character(len=*)  :: mode

    ! Rank of this processor
    integer, intent(in) :: my_rank

    ! Derived type describing simulation setup
    type(run_params) :: setup

    ! Integers used in calculations
    integer :: i, accept, n_steps, ierr

    ! Array for loading configuration arrays and cross-checking
    integer(array_int), allocatable, dimension(:,:,:,:) :: test_config

    ! Temperature and temperature steps
    real(real64) :: beta, temp, sim_temp

    ! Arrays for checking energies of trajectories
    real(real64), dimension(:), allocatable :: trajectory_energy,      &
                                               trajectory_energy_test, &
                                               indices

    ! Arrays for checking energies of trajectories
    real(real64), dimension(:,:,:,:), allocatable :: trajectory_asro, trajectory_asro_test

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

    call print_centered_message('Initialising lattice', '-')

    ! Do an initial setup of the lattice
    call initial_setup(setup, config)

    ! Get the coordination shells of the lattice
    call lattice_shells(setup, shells, config)

    if (trim(mode) .eq. 'test') then
      call ncdf_config_reader('99_ref/'//lattice//'_start_config.nc', test_config, setup)
    else if (trim(mode) .eq. 'generate') then
      call ncdf_grid_state_writer('99_ref/'//lattice//'_start_config.nc', ierr, config, setup)
      test_config=config
    end if

    allocate(trajectory_energy(n_steps))
    allocate(indices(n_steps))
    allocate(trajectory_asro(setup%n_species, setup%n_species, setup%wc_range, n_steps))

    if (.not.(configs_equal(config, test_config))) then
      stop 'Initial random configuration generated different from reference'
    else
      print*, 'Initial random configurations (fixed seed) are identical!'
    end if

    call print_centered_message('Starting trial Metropolis-Hastings sweep', '-')

    ! Choose a simulation temperature
    temp=300.0
    sim_temp = temp*k_b_in_Ry
    beta = 1.0_real64/sim_temp

    do i=1, n_steps
      ! Do some Monte Carlo moves, storing energies

      accept = monte_carlo_step_lattice(setup, config, beta)

      trajectory_energy(i) = setup%full_energy(config)
      trajectory_asro(:,:,:,i) = radial_densities(setup, config, setup%wc_range, shells)
      indices(i) = real(i, kind=real64)
    end do

    call print_centered_message('Completed trial Metropolis-Hastings sweep', '-')

    ! Check the energy trajectories
    if (trim(mode) .eq. 'test') then
      call ncdf_radial_density_reader('99_ref/'//lattice//'_data.nc', trajectory_asro_test, trajectory_energy_test, setup, 256)
    else if (trim(mode) .eq. 'generate') then
      call ncdf_radial_density_writer('99_ref/'//lattice//'_data.nc', trajectory_asro, shells, indices, trajectory_energy, setup)
    end if

    if (trim(mode) .eq. 'test') then
      call ncdf_config_reader('99_ref/'//lattice//'_end_config.nc', test_config, setup)
    else if (trim(mode) .eq. 'generate') then
      call ncdf_grid_state_writer('99_ref/'//lattice//'_end_config.nc', ierr, config, setup)
      test_config=config
    end if

    if (.not.(configs_equal(config, test_config))) then
      stop 'Final configurations are different'
    else
      print*, 'Final random configurations are identical!'
    end if

    call print_centered_message('Tests complete!', '-')

    ! Clean up
    call clean_up_interaction()

    if(allocated(config)) deallocate(config)
    if(allocated(test_config)) deallocate(test_config)
    if(allocated(shells)) deallocate(shells)
    if(allocated(indices)) deallocate(indices)
    if(allocated(trajectory_energy)) deallocate(trajectory_energy)
    if(allocated(trajectory_asro)) deallocate(trajectory_asro)
    if(allocated(trajectory_energy_test)) deallocate(trajectory_energy_test)
    if(allocated(trajectory_asro_test)) deallocate(trajectory_asro_test)

  end subroutine test_suite

  !> @brief   Function for comparing two 1D arrays
  !>
  !> @author  C. D. Woodgate
  !>
  !> @date    2025
  !>
  !> @param  config1 First configuration
  !> @param  config2 Second configuration
  !> @param  tol Tolerance. Defaults to tiny(real32), i.e. single-bit precision
  !>
  !> @return The number of indices where the array elements are not
  !>         within the specified tolerance.
  function array_equal_1D(array1, array2, tol) result(disagreements)

    real(real64), dimension(:), allocatable, intent(in) :: array1, array2
    real(real64), optional :: tol
    real(real32) :: single
    integer :: len1, len2
    integer :: disagreements
    integer :: i

    if (.not.present(tol)) then
      tol = real(epsilon(single), kind=real64)
    end if

    len1 = size(array1)
    len2 = size(array2)

    ! Check that two configs are of the same shape
    if (.not. (len1 .eq. len2)) then
      disagreements = len1
      print*, 'Warning: input arrays not of same length in array_equal_1D'
      return
    end if

    disagreements = 0

    ! Compare the two arrays element by element
    do i=1, len1
      if (abs(array1(i)-array2(i)) .gt. tol) then
        disagreements = disagreements + 1
        return
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
