!> @file    initialise.f90
!>
!> @brief   Assorted routines and tools used when initialising a
!>          simulation
!>
!> @details This module contains routines used for memory allocations,
!>          for initalising the state of a configuration, for setting
!>          and for setting the seed of the random number generator.
!>
!> @author  C. D. Woodgate
!> @date    2019-2025
module initialise

  use kinds
  use constants
  use derived_types
  use shared_data
  use io
  use c_functions
  use bw_hamiltonian
  use random_site
  use neighbour_tables
  use comms, only : comms_finalise, my_rank, comms_wait
  use iso_c_binding
  
  implicit none

  private

  public :: initialise_prng, initialise_function_pointers,             &
            initialise_interaction, initialise_local_arrays,           &
            clean_up_interaction, initialise_local_metropolis_arrays,  &
            local_clean_up, local_metropolis_clean_up, initial_setup

#ifdef USE_MPI

  public :: initialise_global_metropolis_arrays,                       &
            global_metropolis_clean_up

#endif

  contains

  !> @brief   Subroutine to initalise the PRNG
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  seedtime If not 0, time-based random seed(s) are used. If
  !>                  0, static seed is used (for testing).
  !>
  !> @return None
  subroutine initialise_prng(static_seed)

    logical :: static_seed
    integer :: seedtime = 0, status, job_id
    character(len=20) :: slurm_job_id_str

    if(my_rank == 0) then
      write(6,'(18("-"),x,"Seeding random number generator(s)",x,18("-"),/)')
      if (static_seed) then
        write(6,'(2x,"Using fixed seed (for testing)",/)')
      else
        write(6,'(2x,"Using time-based random seed",/)')
      end if
    end if

    if (static_seed) then
      seedtime = 0
    else
      seedtime = 1
    end if

    ! Read slurm job id if it exists
    call GET_ENVIRONMENT_VARIABLE('SLURM_JOB_ID', slurm_job_id_str, status=status)
    if (status == 0) then
      read(slurm_job_id_str, *) job_id
    else
      job_id = 0
    endif

    ! Initialise the prng
    seed = f90_init_genrand(seedtime, int(my_rank, kind=C_INT), job_id)

    call comms_wait()

#ifdef USE_MPI

    print*, 'Thread ', my_rank, ' has seed ', seed

#else

    print*, 'Using seed '

#endif

    call flush(6)

    call comms_wait()

  end subroutine initialise_prng

#ifdef USE_MPI

  !> @brief   Subroutine to initalise the global Metropolis arrays
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  metropolis Derived type containing Metroplis parameters
  !>
  !> @return None
  subroutine initialise_global_metropolis_arrays(setup, metropolis)

    type(run_params), intent(inout) :: setup
    type(metropolis_params), intent(inout) :: metropolis

    if (metropolis%calculate_energies) then
      ! Array for storing energy as a function of temperature
      allocate(av_energies_of_T(metropolis%T_steps))
      ! Array for storing heat capacity as a function of temperature
      allocate(av_C_of_T(metropolis%T_steps))
    end if
    ! Array for storing energy as a function of temperature
    allocate(av_acceptance_of_T(metropolis%T_steps))

    if (metropolis%calculate_asro) then
      ! Radial densities as a function of temperature
      allocate(av_rho_of_T(setup%n_species, setup%n_species, &
                       setup%wc_range, metropolis%T_steps))
      av_rho_of_T = 0.0_real64
    end if

    if (metropolis%calculate_alro) then
      ! Radial densities as a function of temperature
      allocate(av_order_of_T(setup%n_species, setup%n_basis,           &
                             setup%n_grid_1, setup%n_grid_2,           &
                             setup%n_grid_3, metropolis%T_steps))
      av_order_of_T = 0.0_real64
    end if

  end subroutine initialise_global_metropolis_arrays

#endif

  !> @brief   Subroutine to initalise function pointers and the
  !>          lattice-dependent sizes of the simulation cell
  !>
  !> @details Sets the number of atoms in the basis, the total number of
  !>          atoms, and the extent of the configuration array in each
  !>          direction for the requested lattice type, then points the
  !>          energy, random site, and random neighbour function
  !>          pointers at the relevant hard-coded implementations.
  !>
  !> @author  C. D. Woodgate
  !> @author  Claude Fable 5.1
  !> @date    2019-2026
  !>
  !> @param  setup Derived type containing simulation parameters
  !>
  !> @return None
  subroutine initialise_function_pointers(setup)

    type(run_params), intent(inout) :: setup

    setup%full_energy => total_energy

    ! Setup functions used in calculations
    if(trim(setup%lattice) == 'simple_cubic') then
      setup%n_basis = 1
      setup%n_atoms = 8*setup%n_1*setup%n_2*setup%n_3
      setup%n_grid_1 = 2*setup%n_1
      setup%n_grid_2 = 2*setup%n_2
      setup%n_grid_3 = 2*setup%n_3
      ! Every grid point is a lattice site, one lattice parameter apart
      setup%lattice_vectors = reshape( (/ 1.0, 0.0, 0.0,     &
                                          0.0, 1.0, 0.0,     &
                                          0.0, 0.0, 1.0  /), &
                                       (/ 3, 3 /))
      allocate(setup%basis_vectors(3,1))
      setup%basis_vectors(:,1)   = (/ 0.0, 0.0, 0.0 /)
      if (setup%interaction_range .eq. 1) then
        setup%nbr_energy => simple_cubic_energy_1shells
      else if (setup%interaction_range .eq. 2) then
        setup%nbr_energy => simple_cubic_energy_2shells
      else
        call comms_finalise()
        stop 'Unsupported number of shells'
      end if
      setup%rdm_site => simple_cubic_random_site
      setup%rdm_nbr => simple_cubic_random_nbr
    else if(trim(setup%lattice) == 'bcc') then
      setup%n_basis = 1
      setup%n_atoms = 2*setup%n_1*setup%n_2*setup%n_3
      setup%n_grid_1 = 2*setup%n_1
      setup%n_grid_2 = 2*setup%n_2
      setup%n_grid_3 = 2*setup%n_3
      setup%lattice_vectors = reshape( (/ 0.5, 0.0, 0.0,     &
                                          0.0, 0.5, 0.0,     &
                                          0.0, 0.0, 0.5  /), &
                                       (/ 3, 3 /))
!      setup%lattice_vectors = reshape( (/-0.5, 0.5, 0.5,     &
!                                          0.5,-0.5, 0.5,     &
!                                          0.5, 0.5,-0.5  /), &
!                                       (/ 3, 3 /))
      allocate(setup%basis_vectors(3,1))
      setup%basis_vectors(:,1)   = (/ 0.0, 0.0, 0.0 /)
      if (setup%interaction_range .eq. 1) then
        setup%nbr_energy => bcc_energy_1shells
      else if (setup%interaction_range .eq. 2) then
        setup%nbr_energy => bcc_energy_2shells
      else if (setup%interaction_range .eq. 3) then
        setup%nbr_energy => bcc_energy_3shells
      else if (setup%interaction_range .eq. 4) then
        setup%nbr_energy => bcc_energy_4shells
      else if (setup%interaction_range .eq. 5) then
        setup%nbr_energy => bcc_energy_5shells
      else if (setup%interaction_range .eq. 6) then
        setup%nbr_energy => bcc_energy_6shells
      else if (setup%interaction_range .eq. 7) then
        setup%nbr_energy => bcc_energy_7shells
      else if (setup%interaction_range .eq. 8) then
        setup%nbr_energy => bcc_energy_8shells
      else if (setup%interaction_range .eq. 9) then
        setup%nbr_energy => bcc_energy_9shells
      else if (setup%interaction_range .eq. 10) then
        setup%nbr_energy => bcc_energy_10shells
      else
        call comms_finalise()
        stop 'Unsupported number of shells'
      end if
      setup%rdm_site => bcc_random_site
      setup%rdm_nbr => bcc_random_nbr
    else if(trim(setup%lattice) == 'fcc') then

      setup%n_basis = 1
      setup%n_atoms = 4*setup%n_1*setup%n_2*setup%n_3
      setup%n_grid_1 = 2*setup%n_1
      setup%n_grid_2 = 2*setup%n_2
      setup%n_grid_3 = 2*setup%n_3
      setup%lattice_vectors = reshape( (/ 0.5, 0.0, 0.0,     &
                                          0.0, 0.5, 0.0,     &
                                          0.0, 0.0, 0.5  /), &
                                       (/ 3, 3 /))
!      setup%lattice_vectors = reshape( (/ 0.0, 0.5, 0.5,     &
!                                          0.5, 0.0, 0.5,     &
!                                          0.5, 0.5, 0.0  /), &
!                                       (/ 3, 3 /))
      allocate(setup%basis_vectors(3,1))
      setup%basis_vectors(:,1)   = (/ 0.0, 0.0, 0.0 /)
      if (setup%interaction_range .eq. 1) then
        setup%nbr_energy => fcc_energy_1shells
      else if (setup%interaction_range .eq. 2) then
        setup%nbr_energy => fcc_energy_2shells
      else if (setup%interaction_range .eq. 3) then
        setup%nbr_energy => fcc_energy_3shells
      else if (setup%interaction_range .eq. 4) then
        setup%nbr_energy => fcc_energy_4shells
      else if (setup%interaction_range .eq. 5) then
        setup%nbr_energy => fcc_energy_5shells
      else if (setup%interaction_range .eq. 6) then
        setup%nbr_energy => fcc_energy_6shells
      else
        call comms_finalise()
        stop 'Unsupported number of shells'
      end if
      setup%rdm_site => fcc_random_site
      setup%rdm_nbr => fcc_random_nbr
    else if(trim(setup%lattice) == 'ideal_hcp') then

      ! The axial ratio is fixed at the ideal value; the pair tables
      ! are built at that ratio with distance-shell indexing and serve
      ! the ASRO analysis, while the energy uses the hard-coded kernels
      if (setup%c_over_a .gt. 0.0_real64) then
        call comms_finalise()
        stop "c_over_a cannot be set for 'ideal_hcp'; use lattice 'hcp'"
      end if
      if (trim(setup%epi_indexing) .ne. 'shells') then
        call comms_finalise()
        stop "epi_indexing must be 'shells' for 'ideal_hcp'"
      end if
      setup%c_over_a = ideal_hcp_axial_ratio
      call set_hcp_geometry(setup)
      call build_hcp_tables(setup, max(setup%interaction_range,        &
                                       setup%wc_range - 1))

      if (setup%interaction_range .eq. 1) then
        setup%nbr_energy => ideal_hcp_energy_1shells
      else if (setup%interaction_range .eq. 2) then
        setup%nbr_energy => ideal_hcp_energy_2shells
      else if (setup%interaction_range .eq. 3) then
        setup%nbr_energy => ideal_hcp_energy_3shells
      else if (setup%interaction_range .eq. 4) then
        setup%nbr_energy => ideal_hcp_energy_4shells
      else if (setup%interaction_range .eq. 5) then
        setup%nbr_energy => ideal_hcp_energy_5shells
      else if (setup%interaction_range .eq. 6) then
        setup%nbr_energy => ideal_hcp_energy_6shells
      else if (setup%interaction_range .eq. 7) then
        setup%nbr_energy => ideal_hcp_energy_7shells
      else if (setup%interaction_range .eq. 8) then
        setup%nbr_energy => ideal_hcp_energy_8shells
      else if (setup%interaction_range .eq. 9) then
        setup%nbr_energy => ideal_hcp_energy_9shells
      else if (setup%interaction_range .eq. 10) then
        setup%nbr_energy => ideal_hcp_energy_10shells
      else
        call comms_finalise()
        stop 'Unsupported number of shells'
      end if
      setup%rdm_site => ideal_hcp_random_site
      setup%rdm_nbr => ideal_hcp_random_nbr
      call print_pair_tables(setup)
    else if(trim(setup%lattice) == 'hcp') then

      ! General axial ratio: the energy, the ASRO analysis and the
      ! indexing of the interaction file all follow the runtime pair
      ! tables built for this ratio
      if (setup%c_over_a .le. 0.0_real64) then
        call comms_finalise()
        stop "Lattice 'hcp' requires c_over_a in the input file"
      end if
      call set_hcp_geometry(setup)
      call build_hcp_tables(setup, max(setup%interaction_range,        &
                                       setup%wc_range - 1))
      setup%nbr_energy => hcp_energy_tables
      ! Swaps of 'nearest neighbours' use the twelve topological
      ! nearest neighbours of the hcp lattice whatever the axial ratio
      setup%rdm_site => ideal_hcp_random_site
      setup%rdm_nbr => ideal_hcp_random_nbr
      call print_pair_tables(setup)
    else
      call comms_finalise()
      stop 'Lattice type not yet implemented!'
    end if

    ! Check on concentrations and numbers of atoms
    if ((abs(sum(setup%species_concentrations)-1.0_real64) &
                                        .gt. 0.001_real64) &
    .and. (sum(setup%species_numbers) .ne. setup%n_atoms)) &
    then
      call comms_finalise()
      stop 'Invalid numbers of atoms or concentrations specified'
    end if

  end subroutine initialise_function_pointers

  !> @brief   Subroutine to set the sizes, lattice vectors and basis
  !>          positions of the hcp lattice for the axial ratio in
  !>          setup%c_over_a
  !>
  !> @details The hcp lattice is stored on a hexagonal primitive cell
  !>          with a two-atom basis: one grid point per cell, with the
  !>          basis index distinguishing the two atoms. Lattice vectors
  !>          are a_1 = (1, 0, 0), a_2 = (-1/2, sqrt(3)/2, 0) and
  !>          a_3 = (0, 0, c/a) in units of the in-plane lattice
  !>          parameter a, and the basis atoms sit at fractional
  !>          coordinates (1/3, 2/3, 1/4) and (2/3, 1/3, 3/4).
  !>
  !> @author  C. D. Woodgate
  !> @author  Claude Fable 5.1
  !> @date    2026
  !>
  !> @param  setup Derived type containing simulation parameters
  !>
  !> @return None
  subroutine set_hcp_geometry(setup)

    type(run_params), intent(inout) :: setup
    real(real64) :: c_over_a

    c_over_a = setup%c_over_a

    setup%n_basis = 2
    setup%n_atoms = 2*setup%n_1*setup%n_2*setup%n_3
    setup%n_grid_1 = setup%n_1
    setup%n_grid_2 = setup%n_2
    setup%n_grid_3 = setup%n_3

    ! Primitive lattice vectors in units of the in-plane lattice
    ! parameter a
    setup%lattice_vectors(:,1) = (/  1.0_real64,                       &
                                     0.0_real64,                       &
                                     0.0_real64 /)
    setup%lattice_vectors(:,2) = (/ -0.5_real64,                       &
                                     0.5_real64*sqrt(3.0_real64),      &
                                     0.0_real64 /)
    setup%lattice_vectors(:,3) = (/  0.0_real64,                       &
                                     0.0_real64,                       &
                                     c_over_a /)

    ! Cartesian positions of the two basis atoms in units of a,
    ! corresponding to fractional coordinates (1/3, 2/3, 1/4) and
    ! (2/3, 1/3, 3/4)
    if (allocated(setup%basis_vectors)) deallocate(setup%basis_vectors)
    allocate(setup%basis_vectors(3,2))
    setup%basis_vectors(:,1) = (/ 0.0_real64,                          &
                                  1.0_real64/sqrt(3.0_real64),         &
                                  0.25_real64*c_over_a /)
    setup%basis_vectors(:,2) = (/ 0.5_real64,                          &
                                  0.5_real64/sqrt(3.0_real64),         &
                                  0.75_real64*c_over_a /)

  end subroutine set_hcp_geometry

  !> @brief   Subroutine to allocate memory for storing effective pair
  !>          interactions
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !>
  !> @return None
  subroutine initialise_interaction(setup)

    type(run_params), intent(inout) :: setup

    ! Allocate space for atom-atom interchange parameters
    allocate(V_ex(setup%n_species, setup%n_species, setup%interaction_range))

  end subroutine initialise_interaction

  !> @brief   Subroutine to allocate memory used by all processors,
  !>          for all types of sim, even if sim is parallel
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !>
  !> @return None
  subroutine initialise_local_arrays(setup)

    type(run_params), intent(inout) :: setup

    ! On-shell distances
    allocate(shells(setup%wc_range))
    shells = 0.0_real64

    ! Allocate array for storing configuration
    allocate(config(setup%n_basis, setup%n_grid_1, setup%n_grid_2,     &
                    setup%n_grid_3))
    config = 0_array_int

  end subroutine initialise_local_arrays

  !> @brief   Subroutine to allocate memory used by all processors for
  !>          a parallel Metropolis run
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  metropolis Derived type containing Metroplis parameters
  !>
  !> @return None
  subroutine initialise_local_metropolis_arrays(setup, metropolis)

    type(run_params), intent(inout) :: setup
    type(metropolis_params), intent(inout) :: metropolis

    if (metropolis%calculate_energies) then
      ! Array for storing energy as a function of temperature
      allocate(energies_of_T(metropolis%T_steps))
      energies_of_T = 0.0_real64

      ! Array for storing energy as a function of temperature
      allocate(C_of_T(metropolis%T_steps))
      C_of_T = 0.0_real64
    end if

    ! Array for storing energy as a function of temperature
    allocate(acceptance_of_T(metropolis%T_steps))
    acceptance_of_T = 0.0_real64

    if (metropolis%calculate_asro) then
      ! Radial densities as a function of temperature
      allocate(rho_of_T(setup%n_species, setup%n_species, &
                       setup%wc_range, metropolis%T_steps))
      rho_of_T = 0.0_real64
    end if

    if (metropolis%calculate_alro) then
      ! Radial densities as a function of temperature
      allocate(order_of_T(setup%n_species, setup%n_basis,              &
                          setup%n_grid_1, setup%n_grid_2,              &
                          setup%n_grid_3, metropolis%T_steps))
      order_of_T = 0.0_real64
    end if

    ! Array for storing temperatures
    allocate(temperature(metropolis%T_steps))
    temperature = 0.0_real64

  end subroutine initialise_local_metropolis_arrays

  !> @brief   Subroutine to clean up memory used for storing effective 
  !>          pair interactions
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @return None
  subroutine clean_up_interaction()

    deallocate(V_ex)

  end subroutine clean_up_interaction

  !> @brief   Subroutine to clean up memory used by all modes of
  !>          simulation, regardless of parallelism used
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !>
  !> @return None
  subroutine local_clean_up(setup)

    type(run_params), intent(inout) :: setup

    call clean_up_tables(setup)
    if (allocated(setup%basis_vectors)) deallocate(setup%basis_vectors)
    if (allocated(shells)) deallocate(shells)
    if (allocated(config)) deallocate(config)
    if (allocated(setup%species_names)) deallocate(setup%species_names)
    if (allocated(setup%species_concentrations)) deallocate(setup%species_concentrations)

  end subroutine local_clean_up

  !> @brief   Subroutine to clean up memory used by all processors for
  !>          a Metropolis simulation
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !>
  !> @return None
  subroutine local_metropolis_clean_up()

    if (allocated(rho_of_T)) deallocate(rho_of_T)
    if (allocated(order_of_T)) deallocate(order_of_T)
    if (allocated(temperature)) deallocate(temperature)
    if (allocated(energies_of_T)) deallocate(energies_of_T)
    if (allocated(C_of_T)) deallocate(C_of_T)
    if (allocated(acceptance_of_T)) deallocate(acceptance_of_T)

  end subroutine local_metropolis_clean_up

#ifdef USE_MPI

  !> @brief   Subroutine to clean up memory used by process 1 for a
  !>          parallel Metropolis simulation
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @return None
  subroutine global_metropolis_clean_up()

    if (allocated(av_rho_of_T)) deallocate(av_rho_of_T)
    if (allocated(av_order_of_T)) deallocate(av_order_of_T)
    if (allocated(av_energies_of_T)) deallocate(av_energies_of_T)
    if (allocated(av_C_of_T)) deallocate(av_C_of_T)
    if (allocated(av_acceptance_of_T)) deallocate(av_acceptance_of_T)

  end subroutine global_metropolis_clean_up

#endif

  !> @brief   Subroutine to initialise the simulation in a random
  !>          configuration with the correct overall concentration of
  !>          (or numbers of particles of) each species
  !>
  !> @author  C. D. Woodgate
  !> @author  Claude Fable 5.1
  !> @date    2019-2026
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  config Array for simulation configuration
  !>
  !> @return None
  subroutine initial_setup(setup, config)

    type(run_params) :: setup
    integer(array_int), allocatable, dimension(:,:,:,:) :: config
    integer(int32), dimension(4) :: grid_dims
    integer :: i, j, k, b, n_sites, idx
    integer(int64) :: round_err, round_check, inc
    integer(array_int) :: l, n_species
    real(real64) :: rand
    integer(int64), dimension(setup%n_species) :: species_count, check

    check = 0
    if (setup%lattice == 'simple_cubic') then
      n_sites = setup%n_1*setup%n_2*setup%n_3*8
    else if (setup%lattice == 'bcc') then
      n_sites = setup%n_1*setup%n_2*setup%n_3*2
    else if (setup%lattice == 'fcc') then
      n_sites = setup%n_1*setup%n_2*setup%n_3*4
    else if ((setup%lattice == 'ideal_hcp') .or.                      &
             (setup%lattice == 'hcp')) then
      n_sites = setup%n_1*setup%n_2*setup%n_3*2
    else
      call comms_finalise()
      stop 'Lattice type not supported!'
    end if

    ! Make sure my arrays have been allocated and bomb if not
    if (.not. allocated(config)) then
      call comms_finalise()
      stop 'config not allocated in function initial_setup'
    end if

    ! Get the dimensions of the grid I am using
    grid_dims = shape(config)
    n_species = int(setup%n_species, kind=array_int)

    ! Case user specified number of each species
    if (sum(setup%species_numbers) .eq. setup%n_atoms) then
      species_count = setup%species_numbers
      do i=1, setup%n_species
        setup%species_concentrations(i) = real(species_count(i))/real(setup%n_atoms)
      end do
    ! Case user specified concentrations of each species
    !else
    !  do i=1, n_species-1
    !    species_count(i) = floor(real(n_sites)*setup%species_concentrations(i))
    !  end do
    !  species_count(n_species) = n_sites - sum(species_count(1:(n_species-1)))
    !end if
    else
      do i=1, n_species
        species_count(i) = nint(real(n_sites)*setup%species_concentrations(i))
      end do

      ! Work out the rounding error from the concentrations
      round_err = sum(species_count)-n_sites

      ! +1 if round_err is positive
      ! -1 if round_err is negative
      inc = sign(1_int64, -round_err)

      ! Save a copy for incrementing and checking
      round_check = round_err

      ! If we have fewer atoms than we should, gradually add one to each
      ! species, starting from index 1, until correct, avoiding any zero
      ! species. Do the opposite if we have more atoms than we should.
      ! The variable 'inc' determines whether we are adding or removing.
      do while (round_check .ne. 0)
        do j=1, n_species
          if (species_count(j) .eq. 0) cycle
          if (round_check .eq. 0) cycle
          species_count(j) = species_count(j) + inc
          round_check = round_check + inc
        end do
      end do

    end if

    ! Make sure species counts make sense
    if (.not. (sum(species_count).eq.n_sites)) then
      call comms_finalise()
      print*, 'Species counts do not sum to the number of lattice sites'
      print*, 'Species counts: ', species_count
      stop 'Error in initial_setup()'
    end if

    ! Set configuration to be zero
    config = 0_array_int

    ! Set up the lattice
    ! Deal with simple cubic case
    if(trim(setup%lattice) == 'simple_cubic') then
      ! Loop over lattice sites
      do k=1, grid_dims(4)
        do j=1, grid_dims(3)
          do i=1, grid_dims(2)
            do while (config(1,i,j,k) .eq. 0_int32)
              ! Get a random number
              rand = genrand()
              ! Loop over species
              do l=1, n_species
                ! Decide which species to sit on that site
                if ((rand .ge. sum(setup%species_concentrations(0:(l-1)))) .and. &
                    (rand .le. sum(setup%species_concentrations(0:l)))) then
                  if (check(l) .lt. species_count(l)) then
                    config(1,i,j,k) = l
                    check(l) = check(l) + 1
                  end if
                end if
              end do ! Species
            end do ! While
          end do ! i
        end do ! j
      end do ! k
    ! Deal with bcc case
    else if (trim(setup%lattice) == 'bcc') then
      ! Loop over lattice sites
      do k=1, grid_dims(4)
        do j=1, grid_dims(3)/2
          do i=1, grid_dims(2)/2
            do while (config(1,2*i-modulo(k,2),2*j-modulo(k,2),k) .eq. 0_int32)
              ! Get a random number
              rand = genrand()
              ! Loop over species
              do l=1, n_species
                ! Decide which species to sit on that site
                if ((rand .ge. sum(setup%species_concentrations(0:(l-1)))) .and. &
                    (rand .le. sum(setup%species_concentrations(0:l)))) then
                  if (check(l) .lt. species_count(l)) then
                    if( modulo(k, 2) == 1) then
                      config(1,2*i-1, 2*j-1,k) = l
                    else
                      config(1,2*i, 2*j,k) = l
                    end if
                    check(l) = check(l) + 1
                  end if
                end if
              end do ! Species
            end do ! While
          end do ! i
        end do ! j
      end do ! k
    ! Deal with fcc case
    else if (trim(setup%lattice) == 'fcc') then
      ! Loop over lattice sites
      do k=1, grid_dims(4)
        do j=1, grid_dims(3)
          do i=1, grid_dims(2)/2
            idx = 2*i - modulo(k,2)*modulo(j,2) - modulo(k+1,2)*modulo(j+1,2)
            do while (config(1, idx,j,k) .eq. 0_int32)
              ! Get a random number
              rand = genrand()
              ! Loop over species
              do l=1, n_species
                ! Decide which species to sit on that site
                if ((rand .ge. sum(setup%species_concentrations(0:(l-1)))) .and. &
                    (rand .le. sum(setup%species_concentrations(0:l)))) then
                  if (check(l) .lt. species_count(l)) then
                    if( modulo(k, 2) == 1) then
                      if( modulo(j,2) == 1) then
                        config(1, 2*i-1, j, k) = l
                      else
                        config(1, 2*i, j, k) = l
                      end if
                    else
                      if( modulo(j,2) == 1) then
                        config(1, 2*i, j, k) = l
                      else
                        config(1, 2*i-1, j, k) = l
                      end if
                    end if
                    check(l) = check(l) + 1
                  end if
                end if
              end do ! Species
            end do ! While
          end do ! i
        end do ! j
      end do ! k
    ! Deal with the hcp cases: every grid point is a lattice site
    else if((trim(setup%lattice) == 'ideal_hcp') .or.                  &
            (trim(setup%lattice) == 'hcp')) then
      ! Loop over lattice sites
      do k=1, grid_dims(4)
        do j=1, grid_dims(3)
          do i=1, grid_dims(2)
            do b=1, grid_dims(1)
              do while (config(b,i,j,k) .eq. 0_int32)
                ! Get a random number
                rand = genrand()
                ! Loop over species
                do l=1, n_species
                  ! Decide which species to sit on that site
                  if ((rand .ge. sum(setup%species_concentrations(0:(l-1)))) .and. &
                      (rand .le. sum(setup%species_concentrations(0:l)))) then
                    if (check(l) .lt. species_count(l)) then
                      config(b,i,j,k) = l
                      check(l) = check(l) + 1
                    end if
                  end if
                end do ! Species, l
              end do ! While
            end do ! b
          end do ! i
        end do ! j
      end do ! k
    else
      call comms_finalise()
      print*, 'Lattice type: ', setup%lattice, ' not yet implemented'
      stop
    end if

  end subroutine initial_setup


end module initialise
