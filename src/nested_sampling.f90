!> @file    nested_sampling.f90                                  
!>
!> @brief   Assorted routines and tools to perform Nested Sampling                                                                      !
!>
!> @details This module contains routines necessary for the Nested Sampling
!>          calculations.                                                           
!>
!> @author  L. B. Partay
!> @date    2024  
module nested_sampling

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

  !> @brief   Main nested sampling routine.  
  !>
  !> @details This routine performs the nested sampling calculation. Input parameters (such as
  !>          the number of walkers, random walk parameters...etc.) are read from the "ns_input.txt" 
  !>          file. The first section of the routine generate the initial random configurations, which is
  !>          followed by the main nested sampling cycle: recording and discarding the highest energy
  !>          configuration, replacing it with a new one, generated from an existing walker via a series
  !>          of random swaps.
  !>          The routine makes use of the pair_swap subroutine to calculate energies.
  !> 
  !> @param  setup Derived type containing simulation parameters
  !> @param  nested_sampling Derived type containing nested sampling parameters
  !> 
  !> @return None
  !>
  !> @author  L. B. Partay
  !> @date    2024  
  subroutine nested_sampling_main(setup, nested_sampling)

    ! Arrays for storing instances of the system, as nested sampling walkers
    type(run_params) :: setup
    type(ns_params) :: nested_sampling
  
    ! Energy before and after swap and change
    real(real64) :: e_unswapped, e_swapped, delta_e, ener_limit, rnd
    real(real64) :: acc_ratio, rnde
    real(real64), dimension(:), allocatable :: walker_energies
    ! Array for all the walkers configurations, with indices of n_base, n1, n2, n3, i_walker
    integer(int16), dimension(:,:,:,:,:), allocatable :: ns_walkers
    
    ! Occupancy of each site
    integer(int16) :: site1, site2
    integer :: i_walker, i_step, i_iter, irnd, n_at
    integer :: n_acc, i_max(1), extra_steps
  
    ! Random lattice sites
    integer, dimension(4) :: rdm1, rdm2

    ! open and parse input parameters needed for nested sampling
    call read_ns_file("ns_input.inp", nested_sampling)
    open(35,file=nested_sampling%outfile_ener)

    ! creat arrays for storing all the NS walker configurations and their energies
    allocate(ns_walkers(setup%n_basis, 2*setup%n_1, 2*setup%n_2, 2*setup%n_3,nested_sampling%n_walkers))
    allocate(walker_energies(nested_sampling%n_walkers))
    walker_energies=0.0d0

    ! initialise all walkers 
    do i_walker = 1,nested_sampling%n_walkers
    
       ! Set up the energies file for ns_analyse post-procesing (script available from pymatnest package) 
       if (i_walker == 1) then
          n_at = setup%n_1*setup%n_2*setup%n_3*setup%n_basis*setup%n_species
          ! header for energies file
          write(35,*) nested_sampling%n_walkers, 1, 0, 'False', n_at

       endif
      
       ! set up and generate initial random configurations for each walker
       call initial_setup(setup, config)
       ns_walkers(:,:,:,:,i_walker) = config
       
       ! Check the average concentrations
       if (i_walker == 1) call print_particle_count(setup, ns_walkers(:,:,:,:,i_walker))
    
       ! Calculate all the inital energies and print to screen
       ! Add small random number to the energy to avoid configurations to be degenerate in energy
       call random_number(rnde)
       walker_energies(i_walker)=setup%full_energy(ns_walkers(:,:,:,:,i_walker))+rnde*1e-8
       print*, 'Energy', i_walker, 'is: ', walker_energies(i_walker)
    
    enddo
    
    extra_steps=0
    
    !---------------------------------!
    ! Main NS iteration cycle         !
    !---------------------------------!
    print*, '********* MAIN NS CYCLE STARTS ************'
    print*, 'Will do',nested_sampling%n_iter,'iterations, starting with', nested_sampling%n_steps,'MC steps initially.'
    
    do i_iter = 1, nested_sampling%n_iter
    
       ! find highest energy, which is saved and discarded from the list of walkers
       i_max=maxloc(walker_energies)
       ener_limit=walker_energies(i_max(1))
    
       ! save energy to NS ener output: nested_sampling%outfile_ener 
       write(35,*) i_iter, ener_limit
       ! save config to NS xyz output: nested_sampling%outfile_traj
       if (mod(i_iter,nested_sampling%traj_freq)==0) then
          call xyz_writer(nested_sampling%outfile_traj, ns_walkers(:,:,:,:,i_max(1)), setup, trajectory=.true.)
       endif
    
       ! print sampling stage info to screen in every 0.5*number of walkers step
       ! also adjust the number of attempted steps if the acceptance ratio is too low. 
       if (mod(i_iter,int(nested_sampling%n_walkers/2.0))==0) then

          ! acceptance ratio of the previous random walk      
          acc_ratio=real(n_acc)/real(nested_sampling%n_steps+extra_steps)
          
          ! number of accepted steps should be such that 5% of atoms are moved
          if ((n_acc < n_at*0.05).and.(extra_steps<nested_sampling%n_steps*100)) then
             ! increase the extra steps to increase attempted steps
             extra_steps = extra_steps+int(nested_sampling%n_steps)
             print*, i_iter, ener_limit, acc_ratio, n_acc,'n_steps increased to', &
                    & nested_sampling%n_steps+extra_steps
          else 
             print*, i_iter, ener_limit, acc_ratio, n_acc
          endif

       endif
    
       ! Generate a new sample, starting from an existing walker and doing a set of random steps

       ! pick configuration for cloning
       call random_number(rnd)
       irnd=ceiling(rnd*nested_sampling%n_walkers)
       ns_walkers(:,:,:,:,i_max(1)) = ns_walkers(:,:,:,:,irnd)
       walker_energies(i_max(1)) = walker_energies(irnd)
    
       n_acc=0
       i_step=0
       ! do random walk
       do i_step = 1, nested_sampling%n_steps+extra_steps
       
          ! Pick random site 1 and get what atom type it accopies
          rdm1 = setup%rdm_site()
          site1 = ns_walkers(rdm1(1), rdm1(2), rdm1(3), rdm1(4),i_max(1))
          do 
             ! Generate random site 2 that has different atom type
             rdm2 = setup%rdm_site()
             !if (acc_ratio < 0.01) then
             !   rdm2 = ns_setup(i_max(1))%rdm_nbr(rdm1)
             !else
             !   rdm2 = ns_setup(i_max(1))%rdm_site()
             !endif
            
             site2 = ns_walkers(rdm2(1), rdm2(2), rdm2(3), rdm2(4),i_max(1))
             if (site1 /= site2) exit ! if same type, try again
          enddo
    
          ! Compute the energies before and after and work out the change
          e_unswapped = pair_energy(setup, ns_walkers(:,:,:,:,i_max(1)), rdm1, rdm2)
          call pair_swap(ns_walkers(:,:,:,:,i_max(1)), rdm1, rdm2)
          e_swapped = pair_energy(setup, ns_walkers(:,:,:,:,i_max(1)), rdm1, rdm2)
    
          ! since we are always swapping different atom types, the delta_e cannot be zero
          delta_e = e_swapped - e_unswapped
    
          if (walker_energies(i_max(1))+delta_e < ener_limit) then
             !accept
             walker_energies(i_max(1))=walker_energies(i_max(1))+delta_e
             n_acc = n_acc+1
          else
             !reject (swap pairs back)
             call pair_swap(ns_walkers(:,:,:,:,i_max(1)), rdm1, rdm2)
          endif
    
       enddo
       if (n_acc==0) then
         write(*,*) "WARNING! None of the attempted swaps were accepted, hence the cloned configuration", &
                    & "remained identical to the parent."
         write(*,*) "Consider: increasing the number of steps, use a different strategy for swapping pair or", &
                    & "decreasing the number of NS iterations (might have reached low enough energies?)"         
       endif        
    
    enddo
    print*, 'All NS iterations are done, sampling finished.'
    close(35)
    
  
  end subroutine nested_sampling_main

end module nested_sampling
