!> @file    howto_examples.f90
!>
!> @brief   Example subroutines showing the code's functionalities for
!>          developers
!>
!> @details This file contains some example subroutines showcasing the
!>          code's functionalities. They are called by the test program
!>          in example.f90
!>
!> @author  C. D. Woodgate
!> @date    2021-2024
module howto_examples 

  use initialise
  use kinds
  use shared_data
  use io
  use comms
  use c_functions
  use bw_hamiltonian
  use random_site
  use analytics
  use write_netcdf
  use write_xyz
  use metropolis_output
  use metropolis
  
  implicit none

  contains

  !> @brief   Subroutine showing examples of code's functionalities
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  my_rank Rank of current MPI process
  !>
  !> @return None
  subroutine examples(setup, my_rank)

    ! Rank of this processor
    integer, intent(in) :: my_rank

    ! Arrays for storing data
    type(run_params) :: setup
  
    ! Integers used in calculations
    integer :: i,j, div_steps, accept, n_save
    
    ! Temperature and temperature steps
    real(real64) :: beta, temp, sim_temp, current_energy, step_E,     &
                    step_Esq, C, acceptance
  
    ! Name of file for grid state and xyz file at this temperature
    character(len=36) :: xyz_file

    ! Name of file for writing diagnostics at the end
    character(len=43) :: diagnostics_file
  
    ! Name of file for writing radial densities at the end
    character(len=37) :: radial_file

    ! Radial densities array (to populate)
    real(real64), allocatable, dimension(:,:,:) :: r_densities

    ! This is for calculating radial density functions later
    ! (It populates the array "shells".)
    call lattice_shells(setup, shells, config)

    ! Are we swapping just nearest-neighbours or 
    ! allowing any pair of lattice sites?
    if (setup%nbr_swap) then
      setup%mc_step => monte_carlo_step_nbr
    else
      setup%mc_step => monte_carlo_step_lattice
    end if

    !----------!
    ! EXAMPLES !
    !----------!
  
    !--------------------------------------------!
    ! 1. Initialise a random grid and display it !
    !--------------------------------------------!

    ! Set up the initial state of the grid (random occupancies)
    ! on this rank. (Each rank works on different grid currently.)
    call initial_setup(setup, config)

    if(my_rank == 0) then
      print*, ' '
      print*, ' This is what the grid initially looks like: '
      print*, ' '
    end if
  
    ! You can see what it looks like (sort of) using this routine
    call display_grid(config)

    !---------------------------------!
    ! 2. Perform a Metropolis MC step !
    !---------------------------------!

    ! Would normally use j to loop over temperatures
    j = 1

    ! Work out the temperature and corresponding beta
    temp = setup%T + real(j-1, real64)*setup%delta_T
    sim_temp = temp*k_b_in_Ry
    beta = 1.0_real64/sim_temp

    accept=0.0_real64
    
    do while (accept .lt. 0.5_real64)
          accept = setup%mc_step(config, beta)
    end do

    if(my_rank == 0) then
      print*, ' '
      print*, ' This is what the grid looks like after one successful step: '
      print*, ' '
    end if

    call display_grid(config)
  
    !---------------------------------------------------!
    ! 3. Get the energy associated with a configuration !
    !---------------------------------------------------!

    current_energy = setup%full_energy(config)

    if(my_rank == 0) then
      print*, ' '
      print*, ' The energy of that configuration is: ', current_energy, ' Ry'
      print*, ' or: ', current_energy*13.606*1000.0, ' meV'
      print*, ' or: ', current_energy*13.606*1000.0/setup%n_atoms, ' meV/atom'
      print*, ' '
    end if

    !------------------------------------!
    ! 4. Run some kind of equillibration !
    !------------------------------------!
    n_save=floor(real(setup%mc_steps)/real(setup%sample_steps))

    acceptance = 0.0_real64
    step_E = 0.0_real64
    step_Esq = 0.0_real64

    do i=1, setup%mc_steps
    
      ! Make one MC move
      accept = setup%mc_step(config, beta)
  
      acceptance = acceptance + accept

      ! Write percentage progress to screen
      if (mod(i, setup%sample_steps) .eq. 0) then
        current_energy = setup%full_energy(config)
        step_E   = step_E + current_energy
        step_Esq = step_Esq + current_energy**2
      end if
    
    end do

    ! Store the average energy per atom at this temperature
    energies_of_T(j) = step_E/n_save/setup%n_atoms

    ! Heat capacity (per atom) at this temperature
    C = (step_Esq/n_save - (step_E/n_save)**2)/(sim_temp*temp)/setup%n_atoms

    ! Acceptance rate at this temperature
    acceptance_of_T(j) = acceptance/real(setup%mc_steps)
  
    ! Store the specific heat capacity at this temperature
    C_of_T(j) = C
    
    if (my_rank ==0) then
      ! Write that we have completed a particular temperature
      write(6,'(a,f7.2,a)',advance='yes') &
      " Temperature ", temp, " complete on process 0."
      write(6,'(a,f7.2,a)',advance='yes') &
      " Internal energy was ", 13.606_real64*1000*energies_of_T(j), " meV/atom"
      write(6,'(a,f9.4,a)',advance='yes') &
      " Heat capacity was ", C_of_T(j)/k_B_in_Ry, " kB/atom"
      write(6,'(a,f7.2,a,/)',advance='yes') &
      " Swap acceptance rate was ", 100.0*acceptance_of_T(j), "%"
    end if
    
    !------------------------------------!
    ! 5. Write the grid to a .xyz file   !
    !------------------------------------!

    write(xyz_file, '(A11 I3.3 A12 I4.4 F2.1 A4)') 'grids/proc_', &
    my_rank, 'config_at_T_', int(temp), temp-int(temp),'.xyz'
  
    ! Write xyz file
    call xyz_writer(xyz_file, config, setup)
  
    !----------------------------------------------------------------!
    ! 6. Compute the Warren-Cowley ASRO parameters and write to file !
    !----------------------------------------------------------------!

    ! Work out how to find each shell
    call lattice_shells(setup, shells, config)

    ! Compute the radial densities
    r_densities = radial_densities(setup, config, setup%wc_range, shells)
  
    write(radial_file, '(A22 I3.3 A12)') 'radial_densities/proc_', my_rank, '_rho_of_T.nc'

    ! Write the radial densities to file
    call ncdf_radial_density_writer_once(radial_file, r_densities, shells, setup)

    if(my_rank == 0) then
      write(6,'(25("-"),x,"Simulation Complete!",x,25("-"))')
    end if
  
  end subroutine examples

end module howto_examples 
