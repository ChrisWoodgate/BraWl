!> @file    metropolis.f90
!>
!> @brief   Routines for running Metropolis Monte Carlo simulations
!>
!> @details This module contains routines which run Metropolis Monte
!>          Carlo simulations and extract various quantities of
!>          interest.
!>
!> @author  C. D. Woodgate
!> @date    2019-2025
module metropolis

  use initialise
  use kinds
  use shared_data
  use io
  use comms
  use c_functions
  use bw_hamiltonian
  use random_site
  use analytics
  use netcdf_io
  use write_xyz
  use metropolis_output
  
  implicit none

  private

  public :: metropolis_simulated_annealing,                            &
            metropolis_decorrelated_samples,                           &
            monte_carlo_step_nbr,                                      &
            monte_carlo_step_lattice

  contains

  !> @brief   Subroutine for performing simulated annealing using
  !>          the Metropolis Monte Carlo algorithm
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  metropolis Derived type containing Metropolis MC parameters
  !> @param  my_rank Rank of this process
  !>
  !> @return None
  subroutine metropolis_simulated_annealing(setup, metropolis, my_rank)

    ! Rank of this processor
    integer, intent(in) :: my_rank

    ! Derived type describing simulation setup
    type(run_params) :: setup
    type(metropolis_params) :: metropolis

    ! Integers used in calculations
    integer :: i,j,k, accept, ierr, n_sweeps, n_sweep_steps
    integer :: n_save_energy, n_save_asro, n_save_alro, n_save_trajectory

    ! Temperature and temperature steps
    real(real64) :: beta, temp, sim_temp, current_energy, step_E,     &
                    step_Esq, C, acceptance

    ! Name of file for grid state and xyz file at this temperature
    character(len=72) :: grid_file, xyz_file, xyz_trajectory_file

    ! Name of file for writing diagnostics at the end
    character(len=72) :: energy_trajectory_file, asro_trajectory_file

    ! Name of file for writing diagnostics at the end
    character(len=43) :: diagnostics_file
  
    ! Name of file for writing radial densities at the end
    character(len=37) :: radial_file, alro_file

    ! Radial densities at each temperature step
    real(real64), allocatable, dimension(:,:,:) :: r_densities, asro

    ! Long-range order parameters at each temperature step
    real(real64), allocatable, dimension(:,:,:,:,:) :: order

    ! Read the Metropolis control file
    call parse_metropolis_inputs(metropolis, my_rank)

    ! Make the relevant directories
    if ((metropolis%write_final_config_xyz).or.(metropolis%write_final_config_nc)) then
      if(my_rank == 0) call execute_command_line('mkdir -p configs')
    end if
    if (metropolis%calculate_energies) then
      if(my_rank == 0) call execute_command_line('mkdir -p energies')
    end if
    if (metropolis%calculate_asro) then
      if(my_rank == 0) call execute_command_line('mkdir -p asro')
    end if
    if (metropolis%calculate_alro) then
      if(my_rank == 0) call execute_command_line('mkdir -p alro')
    end if
    if ((metropolis%write_trajectory_energy).or.(metropolis%write_trajectory_asro).or.(metropolis%write_trajectory_xyz)) then
      if(my_rank == 0) call execute_command_line('mkdir -p trajectories')
    end if

#ifdef USE_MPI

    ! Initialise some global arrays
    call initialise_global_metropolis_arrays(setup, metropolis)

#endif

    ! Initialise some local arrays
    call initialise_local_metropolis_arrays(setup, metropolis)

    if (metropolis%read_start_config_nc) then
      call ncdf_config_reader(metropolis%start_config_file, config, setup)
    else
      call initial_setup(setup, config)
    end if

    call lattice_shells(setup, shells, config)

    n_save_energy = floor(real(metropolis%n_mc_steps)                  &
                          /real(metropolis%n_sample_steps))

    n_save_asro = floor(real(metropolis%n_mc_steps)                    &
                        /real(metropolis%n_sample_steps_asro))

    n_save_alro = floor(real(metropolis%n_mc_steps)                    &
                        /real(metropolis%n_sample_steps_alro))

    n_save_trajectory = floor(real(metropolis%n_mc_steps)              &
                            /real(metropolis%n_sample_steps_trajectory))

    ! Are we swapping neighbours or on the whole lattice?
    if (metropolis%nbr_swap) then
      setup%mc_step => monte_carlo_step_nbr
    else
      setup%mc_step => monte_carlo_step_lattice
    end if

    !if (metropolis%calculate_asro) then
      ! Allocate memory for radial densities
      allocate(r_densities(setup%n_species, setup%n_species, setup%wc_range))
      allocate(asro(setup%n_species, setup%n_species, setup%wc_range))
    !end if

    !if (metropolis%calculate_alro) then
      ! Allocate memory for order_parameters
      allocate(order(setup%n_species, setup%n_basis, 2*setup%n_1, 2*setup%n_2, 2*setup%n_3))
    !end if

    if(my_rank == 0) then
      write(6,'(24("-"),x,"Commencing Simulation!",x,24("-"),/)')
    end if

    ! Loop over temperature steps
    do j=1, metropolis%T_steps
  
      ! Zero relevant values at each temperature step as needed
      if (allocated(order)) order = 0.0_real64
      if (allocated(r_densities)) r_densities = 0.0_real64
      if (allocated(asro)) r_densities = 0.0_real64

      step_E = 0.0_real64; step_Esq=0.0_real64
      acceptance = 0.0_real64
    
      ! Work out the temperature and corresponding beta
      temp = metropolis%T + real(j-1, real64)*metropolis%delta_T
      sim_temp = temp*k_b_in_Ry
      beta = 1.0_real64/sim_temp
    
      ! Store this in an array
      temperature(j) = temp
    
      !---------!
      ! Burn in !
      !---------!
      if (j==1) then
        if (metropolis%burn_in_start) then
          do i=1, metropolis%n_burn_in_steps
            ! Make one MC move
            accept = setup%mc_step(config, beta)
            acceptance = acceptance + accept
          end do
        end if
      else if (metropolis%T_steps .eq. 1) then
        if (metropolis%burn_in) then
          do i=1, metropolis%n_burn_in_steps
            ! Make one MC move
            accept = setup%mc_step(config, beta)
            acceptance = acceptance + accept
          end do
        end if
      else
        if (metropolis%burn_in) then
          do i=1, metropolis%n_burn_in_steps
            ! Make one MC move
            accept = setup%mc_step(config, beta)
            acceptance = acceptance + accept
          end do
        end if
      end if

      if ((metropolis%burn_in) .or. ((metropolis%burn_in_start).and.(j==1))) then
          if (my_rank ==0) then
            write(6,'(a,f7.2,a)',advance='yes') &
            " Burn-in complete at temperature ", temp, " on process 0."
            write(6,'(a,i7,a)',advance='yes') &
            " Attempted", int(metropolis%n_burn_in_steps), " trial Monte Carlo moves,"
            write(6,'(a,i7,a)',advance='yes') &
            " of which ", int(acceptance), " were accepted,"
            write(6,'(a,f7.2,a,/)',advance='yes') &
            " corresponding to an acceptance rate of ", &
            100.0*acceptance/float(metropolis%n_burn_in_steps), " %"
          end if
      end if

      if (metropolis%write_trajectory_xyz) then

        ! Name of xyz file to which we will write trajectory
        write(xyz_trajectory_file, '(A,I4.4,A,I4.4,F2.1,A)') 'trajectories/proc_', &
        my_rank, '_trajectory_at_T_', int(temp), temp-int(temp),'.xyz'

        ! Delete this file if it already exists
        open(unit=100, iostat=ierr, file=xyz_trajectory_file, status='old')
        if (ierr == 0) close(100, status='delete')

      end if

      if (metropolis%write_trajectory_energy) then

        ! Name of file to which we will write trajectory energies
        write(energy_trajectory_file, '(A,I4.4,A,I4.4,F2.1,A)') 'trajectories/proc_', &
        my_rank, '_energy_trajectory_at_T_', int(temp), temp-int(temp),'.dat'

        ! Delete this file if it already exists
        open(unit=101, iostat=ierr, file=energy_trajectory_file, status='old')
        if (ierr == 0) close(101, status='delete')

      end if

      if (metropolis%write_trajectory_asro) then

        ! Name of file to which we will write trajectory ASRO parameters
        write(asro_trajectory_file, '(A,I4.4,A,I4.4,F2.1,A)') 'trajectories/proc_', &
        my_rank, '_asro_trajectory_at_T_', int(temp), temp-int(temp),'.dat'

        ! Delete this file if it already exists
        open(unit=102, iostat=ierr, file=asro_trajectory_file, status='old')
        if (ierr == 0) close(102, status='delete')

      end if

      !-----------------------------------------------------------------!
      ! Gather information about initial state of trajectory, if needed !
      !-----------------------------------------------------------------!
      ! Storing of data to do with energies
      if (metropolis%calculate_energies) then
        current_energy = setup%full_energy(config)
        if (metropolis%write_trajectory_energy) then
          call energy_trajectory_writer(energy_trajectory_file, 0, current_energy)
        end if
      end if

      ! Storing of data to do with ASRO
      if (metropolis%calculate_asro) then
        asro = radial_densities(setup, config, setup%wc_range, shells)
        if (metropolis%write_trajectory_asro) then
          call asro_trajectory_writer(asro_trajectory_file, 0, asro)
        end if
      end if

      ! Write (or append) trajectory configuration to .xyz file
      if (metropolis%write_trajectory_xyz) then
        call xyz_writer(trim(xyz_trajectory_file), config, setup, .True.)
      end if

      ! Chop the run up into sweeps during which we need no draw no data
      ! (Avoids unneccessary branching)
      n_sweeps = metropolis%n_mc_steps/metropolis%n_sample_steps
      n_sweep_steps = metropolis%n_mc_steps/n_sweeps

      ! Set acceptance rate back to zero for main MC loop
      acceptance = 0.0_real64

      !-----------------------!
      ! Main Monte Carlo loop !
      !-----------------------!
      do i=1, n_sweeps

        do k=1, n_sweep_steps
          ! Make one MC move
          accept = setup%mc_step(config, beta)
          acceptance = acceptance + accept
        end do

        ! Storing of data to do with energies
        if (metropolis%calculate_energies) then
          ! Current energy
          current_energy = setup%full_energy(config)

          ! Add this to total for averaging
          step_E   = step_E + current_energy

          ! Add square to total for averaging
          step_Esq = step_Esq + current_energy**2

          ! Write (or append) trajectory energy to file
          if (metropolis%write_trajectory_energy) then
            if (mod(i, metropolis%n_sample_steps_trajectory) .eq. 0) then
              call energy_trajectory_writer(energy_trajectory_file, i, current_energy)
            end if
          end if
        end if

        ! Storing of data to do with ASRO
        if (metropolis%calculate_asro) then

          if (mod(i, metropolis%n_sample_steps_asro) .eq. 0) then
            asro = radial_densities(setup, config, setup%wc_range, shells)
            ! Add radial densities for averaging
            r_densities = r_densities + asro
          end if

          if (metropolis%write_trajectory_asro) then
            if (mod(i, metropolis%n_sample_steps_trajectory) .eq. 0) then
              call asro_trajectory_writer(asro_trajectory_file, i, asro)
            end if
          end if
        end if

        ! Storing data to do with ALRO
        if (metropolis%calculate_alro) then
          if (mod(i, metropolis%n_sample_steps_alro) .eq. 0) then
            ! Add radial densities for averaging
            call store_state(order, config, setup)
          end if
        end if

        ! Write (or append) trajectory configuration to .xyz file
        if (metropolis%write_trajectory_xyz) then
          if (mod(i, metropolis%n_sample_steps_trajectory) .eq. 0) then
            ! Write xyz trajectory file
            call xyz_writer(trim(xyz_trajectory_file), config, setup, .True.)
          end if
        end if

      end do

      ! Acceptance rate at this temperature
      acceptance_of_T(j) = acceptance/real(metropolis%n_mc_steps)

      if (metropolis%calculate_energies) then
        ! Store the average energy per atom at this temperature
        energies_of_T(j) = step_E/n_save_energy/setup%n_atoms

        ! Heat capacity (per atom) at this temperature
        C = (step_Esq/n_save_energy - (step_E/n_save_energy)**2) &
            /(sim_temp*temp)/setup%n_atoms

        ! Store the specific heat capacity at this temperature
        C_of_T(j) = C
      end if

      if (metropolis%calculate_asro) then
        ! Store the average radial densities at this temperature
        rho_of_T(:,:,:,j) = r_densities/n_save_asro
      end if

      if (metropolis%calculate_alro) then
        ! Store the average radial densities at this temperature
        order_of_T(:,:,:,:,:,j) = order/real(n_save_alro)
      end if

      ! Dump grids as xyz files if needed
      if (metropolis%write_final_config_xyz) then
          write(xyz_file, '(A11 I4.4 A12 I4.4 F2.1 A4)') 'configs/proc_', &
          my_rank, '_config_at_T_', int(temp), temp-int(temp),'.xyz'

          ! Write xyz file
          call xyz_writer(trim(xyz_file), config, setup)
      end if
  
      ! Dump grids as NetCDF files if needed
      if (metropolis%write_final_config_nc) then
        write(grid_file, '(A11 I4.4 A11 I4.4 F2.1 A3)') 'configs/proc_', my_rank, '_grid_at_T_', &
                                             int(temp), temp-int(temp), '.nc'
        ! Write grid to file
        call ncdf_grid_state_writer(trim(grid_file), ierr, config, setup)
      end if

      if (my_rank ==0) then
        ! Write that we have completed a particular temperature
        write(6,'(a,f7.2,a)',advance='yes') &
        " Sampling at temperature ", temp, " complete on process 0."
        write(6,'(a,I10,a)',advance='yes') &
        " Attempted", int(metropolis%n_mc_steps), " trial Monte Carlo moves,"
        write(6,'(a,i10,a)',advance='yes') &
        " of which ", int(acceptance), " were accepted,"
        write(6,'(a,f7.2,a)',advance='yes') &
        " corresponding to an acceptance rate of ", &
        100.0*acceptance/float(metropolis%n_mc_steps), " %"
        write(6,'(a,f7.2,a)',advance='yes') &
        " Average internal energy was ", 13.606_real64*1000*energies_of_T(j), " meV/atom"
        write(6,'(a,f9.4,a,/)',advance='yes') &
        " Estimate of heat capacity is ", C_of_T(j)/k_B_in_Ry, " kB/atom"
      end if
    
    end do ! Loop over temperature

  
    
    ! Write energy diagnostics
    if (metropolis%calculate_energies) then
      write(diagnostics_file, '(A,I4.4,A)') 'energies/proc_', my_rank, &
                                           'energy_diagnostics.dat'
      call diagnostics_writer(trim(diagnostics_file), temperature, &
                              energies_of_T, C_of_T, acceptance_of_T)
    end if

    ! Write radial densities
    if (metropolis%calculate_asro) then
      write(radial_file, '(A,I3.3,A)') 'asro/proc_', my_rank, '_rho_of_T.nc'
      ! Write the radial densities to file
      call ncdf_radial_density_writer(trim(radial_file), rho_of_T, &
                                      shells, temperature, energies_of_T, setup)
    end if

    ! Write radial densities
    if (metropolis%calculate_alro) then
      write(alro_file, '(A22 I3.3 A12)') 'alro/proc_', my_rank, '_rho_of_T.nc'
      ! Write the radial densities to file
      call ncdf_order_writer(trim(alro_file), ierr, order_of_T, temperature, setup)
    end if


#ifdef USE_MPI

    ! Average results across the simulation
    call comms_reduce_metropolis_results(setup, metropolis)

    ! Write averaged results to file if needed
    if (my_rank .eq. 0) then
      ! Write energy diagnostics
      if (metropolis%calculate_energies) then
        call diagnostics_writer('energies/av_energy_diagnostics.dat', temperature, &
                                av_energies_of_T, av_C_of_T, av_acceptance_of_T)
      end if

      ! Write radial densities
      if (metropolis%calculate_asro) then
        !Write the radial densities to file
        call ncdf_radial_density_writer('asro/av_radial_density.nc', av_rho_of_T, &
                                      shells, temperature, av_energies_of_T, setup)
      end if
    end if

#endif

    if (allocated(r_densities)) deallocate(r_densities)
    if (allocated(asro)) deallocate(asro)
    if (allocated(order)) deallocate(order)

    ! Clean up
    call local_metropolis_clean_up()

#ifdef USE_MPI

    ! Clean up
    call global_metropolis_clean_up()

#endif

    if(my_rank == 0) then
      write(6,'(25("-"),x,"Simulation Complete!",x,25("-"))')
    end if
  
  end subroutine metropolis_simulated_annealing

  !> @brief   Subroutine for performing simulated annealing using
  !>          the Metropolis Monte Carlo algorithm, then drawing a set
  !>          of equilibrated samples each N steps apart
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  metropolis Derived type containing Metropolis MC parameters
  !> @param  my_rank Rank of this process
  !>
  !> @return None
  subroutine metropolis_decorrelated_samples(setup, metropolis, my_rank)

    ! Rank of this processor
    integer, intent(in) :: my_rank

    ! Derived type describing simulation setup
    type(run_params) :: setup
    type(metropolis_params) :: metropolis
  
    ! Integers used in calculations
    integer :: i,j, accept
    integer :: n_save_energy, n_save_asro, n_save_alro, n_save_trajectory

    ! Temperature and temperature steps
    real(real64) :: beta, temp, sim_temp, current_energy, acceptance
  
    ! Name of xyz file
    character(len=42) :: xyz_file

    ! Read the Metropolis control file
    call parse_metropolis_inputs(metropolis, my_rank)

    ! Set up the lattice
    call initial_setup(setup, config)

#ifdef USE_MPI

    ! Initialise some global arrays
    call initialise_global_metropolis_arrays(setup, metropolis)

#endif

    ! Initialise some local arrays
    call initialise_local_metropolis_arrays(setup, metropolis)

    call lattice_shells(setup, shells, config)

    n_save_energy = floor(real(metropolis%n_mc_steps)                  &
                          /real(metropolis%n_sample_steps))

    n_save_asro = floor(real(metropolis%n_mc_steps)                    &
                        /real(metropolis%n_sample_steps_asro))

    n_save_alro = floor(real(metropolis%n_mc_steps)                    &
                        /real(metropolis%n_sample_steps_alro))

    n_save_trajectory = floor(real(metropolis%n_mc_steps)              &
                            /real(metropolis%n_sample_steps_trajectory))

    ! Are we swapping neighbours or on the whole lattice?
    if (metropolis%nbr_swap) then
      setup%mc_step => monte_carlo_step_nbr
    else
      setup%mc_step => monte_carlo_step_lattice
    end if

    if(my_rank == 0) then
      write(6,'(/,72("-"),/)')
      write(6,'(24("-"),x,"Commencing Simulation!",x,24("-"),/)')
    end if

    !---------------------------------------------------!
    ! Burn-in at each temperature (simulated annealing) !
    !---------------------------------------------------!
    do j=1, metropolis%T_steps
  
      ! Work out the temperature and corresponding beta
      temp = metropolis%T + real(j-1, real64)*metropolis%delta_T
      sim_temp = temp*k_b_in_Ry
      beta = 1.0_real64/sim_temp
    
      ! Burn in
      if (metropolis%burn_in) then

        acceptance = 0.0_real64

        do i=1, metropolis%n_burn_in_steps
          ! Make one MC move
          accept = setup%mc_step(config, beta)
          acceptance = acceptance + accept
        end do

        if(my_rank == 0) then
          write(6,'(a,f7.2,a)',advance='yes') &
          " Burn-in complete at temperature ", temp, " on process 0."
          write(6,'(a,i7,a)',advance='yes') &
          " Accepted ", int(acceptance), " Monte Carlo moves at this temperature,"
          write(6,'(a,f7.2,a,/)',advance='yes') &
          " Corresponding to an acceptance rate of ", &
          100.0*acceptance/float(metropolis%n_burn_in_steps), " %"
        end if

      end if

    end do ! Loop over temperature

    !--------------------!
    ! Target Temperature !
    !--------------------!
 
    acceptance=0

    do i=1, metropolis%n_mc_steps
    
        ! Make one MC move
        accept = setup%mc_step(config, beta)
  
        acceptance = acceptance + accept

        ! Draw samples
        if (mod(i, metropolis%n_sample_steps_asro) .eq. 0) then

          ! Get the energy of this configuration
          current_energy = setup%full_energy(config)

          write(xyz_file, '(A11 I3.3 A8 I4.4 A6 I4.4 F2.1 A4)') &
          'configs/proc_', my_rank, '_config_',                   &
          int(i/metropolis%n_sample_steps_asro), '_at_T_', int(temp),&
          temp-int(temp),'.xyz'

          ! Write xyz file
          call xyz_writer(xyz_file, config, setup)
 
          if (my_rank == 0) then
            write(6,'(a,i7,a,/)',advance='yes') &
            " Accepted an additional ", int(acceptance), " Monte Carlo moves before sample."
          end if

          acceptance=0
          
        end if
    
      end do

    ! Clean up
    call local_metropolis_clean_up()

#ifdef USE_MPI

    ! Clean up
    call global_metropolis_clean_up()

#endif

    if(my_rank == 0) then
      write(6,'(25("-"),x,"Simulation Complete!",x,25("-"))')
    end if
  
  end subroutine metropolis_decorrelated_samples

  !> @brief   Subroutine for performing a trial Metropolis-Kawasaki swap
  !>          assuming pair can be swapped across the entire lattice
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  config Current grid state
  !> @param  beta Inverse temperature to use in acceptance rate
  !>
  !> @return None
  function monte_carlo_step_lattice(setup, config, beta) result(accept)

    integer(array_int), dimension(:,:,:,:) :: config
    class(run_params), intent(in) :: setup
    integer :: accept
    integer, dimension(4) :: rdm1, rdm2
    real(real64) , intent(in) :: beta
    real(real64) :: e_unswapped, e_swapped, delta_e
    integer(array_int) :: site1, site2

    ! Pick two random sites on the lattice
    rdm1 = setup%rdm_site()
    rdm2 = setup%rdm_site()

    ! Find out which atoms are on those sites
    site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
    site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

    ! If they are the same chemical species, we don't need to proceed
    ! further
    ! Note: previously this was counted as an 'unaccepted' move. Now we
    !       accept the move as this recovers a 100% acceptance rate in
    !       the limit T->\infty.
    if (site1 == site2) then
      accept = 1
      return
    end if

    ! Get the energy associated with the current configuration
    ! Note: This is a 'local' evaluation - as the interaction is
    !       short-ranged in real space, we do not need to evaluate the
    !       energy of the entire cell for large cells.
    e_unswapped = pair_energy(setup, config, rdm1, rdm2)

    ! Trial a swap of the selected pair of atoms
    call pair_swap(config, rdm1, rdm2)

    ! Evaluate the energy if the pair of atoms is swapped
    e_swapped = pair_energy(setup, config, rdm1, rdm2)   

    ! Assess the change in energy as a result of the swap
    delta_e = e_swapped - e_unswapped

    ! Metropolis condition
    ! If the change in energy is negative, always accept it
    if(delta_e .lt. 0.0_real64) then
      accept = 1
      return

    ! Else, if the change in energy is positive, accept it 
    ! if exp(-beta*deltaE) > chosen random number in [0,1]
    else if (genrand() .lt. exp(-beta*delta_E)) then
      accept = 1
      return

    ! Otherwise, swap is rejected
    else
      accept = 0
      ! As swap has been rejected, swap the pair of atoms back
      call pair_swap(config, rdm1, rdm2)
    end if

  end function monte_carlo_step_lattice

  !> @brief   Subroutine for performing a trial Metropolis-Kawasaki swap
  !>          assuming only nearest-neighbour pair can be swapped.
  !>
  !> @details More physical than allowing swaps across the entire
  !>          lattice, but much slower to reach equilibrium.
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  config Current grid state
  !> @param  beta Inverse temperature to use in acceptance rate
  !>
  !> @return None
  function monte_carlo_step_nbr(setup, config, beta) result(accept)

    integer(array_int), dimension(:,:,:,:) :: config
    class(run_params), intent(in) :: setup
    integer :: accept
    integer, dimension(4) :: rdm1, rdm2
    real(real64) , intent(in) :: beta
    real(real64) :: e_unswapped, e_swapped, delta_e
    integer(array_int) :: site1, site2

    ! Pick a random site on the lattice and its neighbour
    rdm1 = setup%rdm_site()
    rdm2 = setup%rdm_nbr(rdm1)

    ! Find out which atoms are on those sites
    site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
    site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

    ! If they are the same chemical species, we don't need to proceed
    ! further
    ! Note: previously this was counted as an 'unaccepted' move. Now we
    !       accept the move as this recovers a 100% acceptance rate in
    !       the limit T->\infty.
    if (site1 == site2) then
      accept = 1
      return
    end if

    ! Get the energy associated with the current configuration
    ! Note: This is a 'local' evaluation - as the interaction is
    !       short-ranged in real space, we do not need to evaluate the
    !       energy of the entire cell for large cells.
    e_unswapped = pair_energy(setup, config, rdm1, rdm2)

    ! Trial a swap of the selected pair of atoms
    call pair_swap(config, rdm1, rdm2)

    ! Evaluate the energy if the pair of atoms is swapped
    e_swapped = pair_energy(setup, config, rdm1, rdm2)   

    ! Assess the change in energy as a result of the swap
    delta_e = e_swapped - e_unswapped

    ! Metropolis condition
    ! If the change in energy is negative, always accept it
    if(delta_e .lt. 0.0_real64) then
      accept = 1
      return

    ! Else, if the change in energy is positive, accept it 
    ! if exp(-beta*deltaE) > chosen random number in [0,1]
    else if (genrand() .lt. exp(-beta*delta_E)) then
      accept = 1
      return

    ! Otherwise, swap is rejected
    else
      accept = 0
      ! As swap has been rejected, swap the pair of atoms back
      call pair_swap(config, rdm1, rdm2)
    end if

  end function monte_carlo_step_nbr

end module metropolis
