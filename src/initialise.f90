!----------------------------------------------------------------------!
! initialise.f90                                                       !
!                                                                      !
! Module with routines associated with initialising a simulation.      !
!                                                                      !
! C. D. Woodgate,  Warwick                                        2023 !
!----------------------------------------------------------------------!
module initialise

  use kinds
  use shared_data
  use io
  use c_functions
  use energetics
  use random_site
  use comms
  
  implicit none

  contains

  !--------------------------------------------------------------------!
  ! Subroutine to initialise the PNRG                                  !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine initialise_pnrg(seedtime)
    integer :: seedtime

    if(my_rank == 0) then
      write(6,'(17("-"),x,"Initialising random number generators",x,16("-"),/)')
    end if

    ! Initialise the prng
    seed = f90_init_genrand(seedtime, int(my_rank, kind=C_INT))

    call comms_wait()
    print*, 'Thread ', my_rank, ' has seed ', seed
    call flush(6)

    call comms_wait()

  end subroutine initialise_pnrg

  !--------------------------------------------------------------------!
  ! Subroutine to initialise arrays for storing averages (rank 1 only) !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine initialise_global_arrays(setup)
    type(run_params), intent(inout) :: setup
    ! Array for storing energy as a function of temperature
    allocate(av_energies_of_T(setup%T_steps))
    ! Array for storing energy as a function of temperature
    allocate(av_C_of_T(setup%T_steps))
    ! Array for storing energy as a function of temperature
    allocate(av_acceptance_of_T(setup%T_steps))
    ! Radial densities as a function of temperature
    allocate(av_rho_of_T(setup%n_species, setup%n_species, &
                     setup%wc_range, setup%T_steps))
    av_rho_of_T = 0.0_real64
  end subroutine initialise_global_arrays

  !--------------------------------------------------------------------!
  ! Subroutine to initialise function pointers                         !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine initialise_function_pointers(setup)
    type(run_params), intent(inout) :: setup

    setup%full_energy => total_energy

    ! Setup functions used in calculations
    if(trim(setup%lattice) == 'simple_cubic') then
      setup%n_basis = 1
      setup%n_atoms = 8*setup%n_1*setup%n_2*setup%n_3
      if (setup%interaction_range .eq. 1) then
        setup%nbr_energy => simple_cubic_energy_1shells
      else
        print*, 'Unsupported number of shells'
        stop
      end if
      setup%rdm_site => simple_cubic_random_site
      setup%rdm_nbr => simple_cubic_random_nbr
    else if(trim(setup%lattice) == 'bcc') then
      setup%n_basis = 1
      setup%n_atoms = 2*setup%n_1*setup%n_2*setup%n_3
      setup%lattice_vectors = reshape( (/-0.5, 0.5, 0.5,     &
                                          0.5,-0.5, 0.5,     &
                                          0.5, 0.5,-0.5  /), &
                                       (/ 3, 3 /))
      setup%basis_vectors   = (/ 0.0, 0.0, 0.0 /)
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
      else
        print*, 'Unsupported number of shells'
        stop
      end if
      setup%rdm_site => bcc_random_site
      setup%rdm_nbr => bcc_random_nbr
    else if(trim(setup%lattice) == 'fcc') then
      setup%n_basis = 1
      setup%n_atoms = 4*setup%n_1*setup%n_2*setup%n_3
      setup%lattice_vectors = reshape( (/ 0.0, 0.5, 0.5,     &
                                          0.5, 0.0, 0.5,     &
                                          0.5, 0.5, 0.0  /), &
                                       (/ 3, 3 /))
      setup%basis_vectors   = (/ 0.0, 0.0, 0.0 /)
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
        print*, 'Unsupported number of shells'
        stop
      end if
      setup%rdm_site => fcc_random_site
      setup%rdm_nbr => fcc_random_nbr
    ! Bomb if we ask for a lattice type we don't have!
    else
      print*, 'Lattice type not yet implemented!'
      stop
    end if
  end subroutine initialise_function_pointers

  !--------------------------------------------------------------------!
  ! Subroutine to allocate memory for storing interaction parameters   !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine initialise_interaction(setup)
    type(run_params), intent(inout) :: setup

    ! Allocate space for atom-atom interchange parameters
    allocate(V_ex(setup%n_species, setup%n_species, setup%interaction_range))

  end subroutine initialise_interaction

  !--------------------------------------------------------------------!
  ! Subroutine to allocate memory for all processors                   !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine initialise_local_arrays(setup)
    type(run_params), intent(inout) :: setup
    ! Array for storing energy as a function of temperature
    allocate(energies_of_T(setup%T_steps))
    ! Array for storing energy as a function of temperature
    allocate(C_of_T(setup%T_steps))
    ! Array for storing energy as a function of temperature
    allocate(acceptance_of_T(setup%T_steps))
    ! Radial densities as a function of temperature
    allocate(rho_of_T(setup%n_species, setup%n_species, &
                     setup%wc_range, setup%T_steps))
    rho_of_T = 0.0_real64
    ! On-shell distances
    allocate(shells(setup%wc_range))
    shells = 0.0_real64
    ! Array for storing temperatures
    allocate(temperature(setup%T_steps))
    ! Allocate array for storing configuration
    allocate(config(setup%n_basis, 2*setup%n_1, 2*setup%n_2, 2*setup%n_3))
  end subroutine initialise_local_arrays

  !--------------------------------------------------------------------!
  ! Subroutine to clean up interaction array                           !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine clean_up_interaction(setup)
    type(run_params), intent(inout) :: setup

    deallocate(V_ex)

  end subroutine clean_up_interaction

  !--------------------------------------------------------------------!
  ! Subroutine to clean up local arrays                                !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine local_clean_up(setup)
    type(run_params), intent(inout) :: setup
    deallocate(rho_of_T)
    deallocate(shells)
    deallocate(temperature)
    deallocate(config)
    deallocate(energies_of_T, C_of_T, acceptance_of_T)
    deallocate(setup%species_names, setup%species_concentrations)
  end subroutine local_clean_up

  !--------------------------------------------------------------------!
  ! Subroutine to clean up average arrays                              !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine global_clean_up()
    deallocate(av_rho_of_T)
    deallocate(av_energies_of_T, av_C_of_T, av_acceptance_of_T)
  end subroutine global_clean_up

  !--------------------------------------------------------------------!
  ! Subroutine to initialise the grid with a random configuration      !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine initial_setup(setup, config)
    type(run_params) :: setup
    integer(int16), allocatable, dimension(:,:,:,:) :: config
    integer(int32), dimension(4) :: grid_dims
    integer :: i, j, k, n_sites, idx
    integer(int16) :: l, n_species
    real(real64) :: rand
    integer(int32), dimension(setup%n_species) :: species_count, check

    check = 0
    if (setup%lattice == 'simple_cubic') then
      n_sites = setup%n_1*setup%n_2*setup%n_3*8
    else if (setup%lattice == 'bcc') then
      n_sites = setup%n_1*setup%n_2*setup%n_3*2
    else if (setup%lattice == 'fcc') then
      n_sites = setup%n_1*setup%n_2*setup%n_3*4
    else
      stop 'Lattice type not supported!'
    end if

    ! Make sure my arrays have been allocated and bomb if not
    if (.not. allocated(config)) then
      stop 'config not allocated in function initial_setup'
    end if

    ! Get the dimensions of the grid I am using
    grid_dims = shape(config)
    n_species = int(setup%n_species, kind=int16)

    do i=1, n_species-1
      species_count(i) = floor(real(n_sites)*setup%species_concentrations(i))
    end do
    species_count(n_species) = n_sites - sum(species_count(1:(n_species-1)))

    ! Set configuration to be zero
    config = 0_int16

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
    else
      print*, 'Lattice type: ', setup%lattice, ' not yet implemented'
      stop
    end if

  end subroutine initial_setup

end module initialise
