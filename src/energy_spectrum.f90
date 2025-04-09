!----------------------------------------------------------------------!
! Energy Spectrum Module                                               !
!                                                                      !
! This program runs MC for the given system and attempts to find       !
! lowest energy state.                                                 !
!                                                                      !
! H. Naguszewski, Warwick                                         2024 !
!----------------------------------------------------------------------!

module energy_spectrum
  use initialise
  use kinds
  use shared_data
  use c_functions
  use random_site
  use metropolis
  use mpi

  implicit none

  contains

  subroutine es_main(setup, es_setup, my_rank)
    ! Rank of this processor
    integer, intent(in) :: my_rank

    ! Arrays for storing data
    type(run_params) :: setup
    type(es_params) :: es_setup

    ! Integers used in calculations
    integer :: ierr, unique_energy, i, energy_min_loc, energy_max_loc, ierror

    ! Temperature and temperature steps
    real(real64) :: acceptance, step, energy_to_ry, min_energy, min_change, max_change
    real(real64) :: min_energy_buf, min_change_buf, max_change_buf

    energy_to_ry = setup%n_atoms/(Ry_to_eV*1000)

    ! Set up the lattice
    call initial_setup(setup, config)

    call lattice_shells(setup, shells, config)

    min_energy = setup%full_energy(config)
    min_change = HUGE(min_change)
    max_change = TINY(max_change)

    if (my_rank == 0) then
      write (6, '(/,72("-"),/)')
      write (6, '(24("-"),x,"Commencing Simulation!",x,24("-"),/)')
      print *, "Number of atoms", setup%n_atoms
    end if

    acceptance = run_es_sweeps(setup, es_setup, config, min_energy, min_change, max_change)
    
    call MPI_REDUCE(min_energy, min_energy_buf, 1, MPI_DOUBLE_PRECISION, &
    MPI_MIN, 0, MPI_COMM_WORLD, ierror)
    call MPI_REDUCE(min_change, min_change_buf, 1, MPI_DOUBLE_PRECISION, &
    MPI_MIN, 0, MPI_COMM_WORLD, ierror)
    call MPI_REDUCE(max_change, max_change_buf, 1, MPI_DOUBLE_PRECISION, &
    MPI_MAX, 0, MPI_COMM_WORLD, ierror)
    if (my_rank == 0) then
      print *, "Energy Reached", min_energy/energy_to_ry, "meV Min change", min_change/energy_to_ry, &
      "meV Max change", max_change/energy_to_ry, "meV"
    end if

    call comms_wait()
    if (my_rank == 0) then
      write (*, *)
      write (6, '(25("-"),x,"Simulation Complete!",x,25("-"))')
    end if

  end subroutine es_main

  function run_es_sweeps(setup, es_setup, config, min_energy, min_change, max_change) result(acceptance)
    integer(int16), dimension(:, :, :, :) :: config
    class(run_params), intent(in) :: setup
    class(es_params), intent(in) :: es_setup
    real(real64), intent(inout) :: min_energy, min_change, max_change

    integer, dimension(4) :: rdm1, rdm2
    real(real64) :: e_swapped, e_unswapped, pair_swapped, pair_unswapped, eps, delta_e, beta
    integer :: acceptance, i, cycle_loop, reject
    integer(int16) :: site1, site2

    ! Establish total energy before any moves
    e_unswapped = setup%full_energy(config)
    e_swapped = e_unswapped

    acceptance = 0.0_real64
    beta = 1.0_real64/(k_b_in_Ry*10.0_real64)

    do i = 1, es_setup%mc_sweeps*setup%n_atoms
      ! Make one MC trial
      ! Generate random numbers
      rdm1 = setup%rdm_site()
      rdm2 = setup%rdm_site()

      ! Get what is on those sites
      site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
      site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

      ! Calculate energy if different species
      if (site1 /= site2) then
        pair_unswapped = pair_energy(setup, config, rdm1, rdm2)

        call pair_swap(config, rdm1, rdm2)

        pair_swapped = pair_energy(setup, config, rdm1, rdm2)
        e_swapped = e_unswapped - pair_unswapped + pair_swapped

        delta_e = e_swapped - e_unswapped
        if (delta_e < min_change) then
          min_change = ABS(delta_e)
        end if
        if (delta_e > max_change) then
          max_change = ABS(delta_e)
        end if

        ! Accept or reject move
        if (genrand() .lt. exp(-delta_e*beta)) then
          e_unswapped = e_swapped
        else
          call pair_swap(config, rdm1, rdm2)
        end if

        if (e_unswapped < min_energy) then
          min_energy = e_unswapped
        end if
      end if
    end do

    !print*, "Energy Drift: ", ABS(setup%full_energy(config) - e_unswapped)/ABS(setup%full_energy(config))
  end function run_es_sweeps
end module energy_spectrum

