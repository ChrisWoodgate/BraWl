!> @file    analytics.f90
!>
!> @brief   Assorted routines and tools for analysing a simulation
!>
!> @details This module contains routines for analysing various aspects
!>          of a simulation, including counting the number of particles,
!>          calculating the neighbours of a given site on the lattice,
!>          and calculating the atomic short-range order parameters.
!>
!> @author  C. D. Woodgate
!> @date    2019-2025
module analytics

  use kinds
  use derived_types
  use constants
  use shared_data
  use io
  use display
  
  implicit none

  private

  public :: store_state, average_state, total_particle_count,          &
            print_particle_count, lattice_shells, radial_densities

  contains

  !> @brief   Subroutine to add this configuration to the average
  !>
  !> @details This is needed to compute the average occupancy of each
  !>          lattice site at the end of a simulation.
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  averages Array of floats where averages are being stored
  !> @param  config Current atomic configuration
  !> @param  setup Derived type containing simulation parameters
  !>
  !> @return None
  subroutine store_state(averages, config, setup)
    !integer(array_int), allocatable, dimension(:,:,:), intent(in) :: config
    integer(array_int), dimension(:,:,:), intent(in) :: config
    type(run_params), intent(in) :: setup
    real(real64), dimension(:,:,:,:), intent(inout), allocatable :: averages
    integer :: i,j,k,l

    do i=1, setup%n_species
      do l=1, 2*setup%n_3
        do k=1, 2*setup%n_2
          do j=1, 2*setup%n_1
            if (config(j,k,l) == i) then
              averages(i,j,k,l) = averages(i,j,k,l) + 1.0_real64
            end if
          end do
        end do
      end do
    end do
  end subroutine store_state

  !> @brief   Subroutine to compute average occupancies
  !>
  !> @details At the end of the simulation, divide the stored
  !>          occupancies by the number of times we stored them, to get
  !>          the average occupancy of each lattice site.
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  averages Array of floats where averages are being stored
  !> @param  setup Derived type containing simulation parameters
  !> @param  n_steps Number of steps at which the state was stored
  !>
  !> @return None
  subroutine average_state(averages, setup, n_steps)
    type(run_params), intent(in) :: setup
    integer, intent(in) :: n_steps
    real(real64), dimension(:,:,:,:), intent(inout), allocatable :: averages
    integer :: i

    do i=1, setup%n_species
      averages(i,:,:,:) = (1.0/real(n_steps, real64))*averages(i,:,:,:)
    end do
  end subroutine average_state
  
  !> @brief   Function to count the total number of particles in the box
  !>
  !> @details This function was mainly used for testing during code
  !>          development, to make sure no particles were dissapearing.
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  config Current atomic configuration
  !>
  !> @return The total number of particles in the simulation cell
  function total_particle_count(setup, config) result(total_count)
    !integer(array_int), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(array_int), dimension(:,:,:,:), intent(in) :: config
    type(run_params) :: setup
    integer :: total_count
    integer :: i,j,k,b

    total_count = 0

    do b=1, setup%n_basis
      do k=1, 2*setup%n_3
        do j=1, 2*setup%n_2
          do i=1, 2*setup%n_1
            if (config(b,i,j,k) .ne. 0_array_int) then
              total_count = total_count + 1
            end if
          end do
        end do
      end do
    end do

  end function total_particle_count

  !> @brief   Subroutine to print the number of atoms of each species
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2023
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  config Current atomic configuration
  !>
  !> @return None
  subroutine print_particle_count(setup, config, my_rank)
    !integer(array_int), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(array_int), dimension(:,:,:,:), intent(in) :: config
    type(run_params) :: setup
    integer, dimension(4) :: sizes
    integer, dimension(:), allocatable :: species_count
    integer :: i,j,k, n, my_rank

    if(my_rank == 0) then
      write(6,'(16("-"),x,"Checking contents of simulation cell(s)",x, 17("-"),/)')
    end if

    sizes = shape(config)

    allocate(species_count(setup%n_species))

    species_count = 0
    n=0

    do k=1, sizes(4)
      do j=1, sizes(3)
        do i=1, sizes(2)
          if (config(1,i,j,k) .ne. 0_array_int) then
            n = n+1
            species_count(config(1,i,j,k)) = &
              species_count(config(1,i,j,k)) + 1
          end if
        end do
      end do
    end do

    if(my_rank == 0) then
      write(6,'(x,"Simulation cell is",x,A,/)') setup%lattice
      write(6,'(x,"There are:",/)')
      write(6,'(x,I3,x,"cells in direction 1,")') setup%n_1
      write(6,'(x,I3,x,"cells in direction 2,")') setup%n_2
      write(6,'(x,I3,x,"cells in direction 3,",/)') setup%n_3
      write(6,'(x,"for a total of ",I5,x,"atoms in the cell.",/)') n

      write(6,'(x,"The breakdown is:",/)')

      write(6,'(x,"Index | Element | Number of Atoms")')
      write(6,'(x,33("-"))')
      do i=1, setup%n_species
        write(6,'(x,I5," |      ", A, " | ", I9)')                         &
              i, setup%species_names(i), species_count(i)
      end do
    end if

    deallocate(species_count)

    if(my_rank == 0) then
      write(6,'(/,17("-"),x,"End of info about simulation cell(s)",x, 17("-"),/)')
    end if

  end subroutine print_particle_count

  !> @brief   Subroutine to compute lattice shell distances
  !>
  !> @details Could just calculate these once, but this routine is good
  !>          for a variety of lattice types.
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  shells Array where shell distances will be stored
  !> @param  config Current atomic configuration
  !>
  !> @return None
  subroutine lattice_shells(setup, shells, configuration)

    integer(array_int), dimension(:,:,:,:) :: configuration
    type(run_params), intent(in) :: setup
    integer :: i,j,k,b,l
    real(real64) :: dist
    !real(real64), dimension(3) :: r_vec
    real(real64), dimension(:), allocatable :: all_shells, shells

    ! Factor of eight to account for the fact that simulation
    ! doubles number of cells in each direction to build lattice
    allocate(all_shells(8*setup%n_1*setup%n_2*setup%n_3*setup%n_basis+1))

    all_shells = 0.0_real64
    shells = 0.0_real64

    l = 1

!    ! Loop over all lattice sites
!    do k=1, 2*setup%n_3
!      do j=1, 2*setup%n_2
!        do i=1, 2*setup%n_1
!          do b=1, setup%n_basis
!            r_vec = real(i)*setup%lattice_vectors(:,1) + &
!                    real(j)*setup%lattice_vectors(:,2) + &
!                    real(k)*setup%lattice_vectors(:,3) + &
!                    real(b)*setup%basis_vectors
!            dist = norm2(r_vec)
!            all_shells(l) = dist
!            l=l+1
!         end do
!        end do
!      end do
!    end do

    ! Loop over all lattice sites
    do k=1, 2*setup%n_3
      do j=1, 2*setup%n_2
        do i=1, 2*setup%n_1
          do b=1, setup%n_basis
            ! Cycle if this lattice site is empty
            if (configuration(b,i,j,k) .eq. 0_array_int) cycle
            dist     = sqrt(real((k-1)**2) + &
                            real((j-1)**2) + &
                            real((i-1)**2))
            all_shells(l) = dist
            l=l+1
         end do
        end do
      end do
    end do

    ! Order the list
    call quicksort(all_shells)

    ! Counter for how many non-repeated distances
    ! we have counted
    l=1

    ! Count the non-repeated distances
    do i=1, size(all_shells)-1
      if (abs(all_shells(i)-all_shells(i+1)) .lt. 1e-3_real64) cycle
      shells(l) = all_shells(i)
      l=l+1
      if (l .gt. setup%wc_range) exit
    end do

    ! Deallocate the array of all distances
    deallocate(all_shells)

  end subroutine lattice_shells

  !> @brief   Subroutine to compute radial densities (ASRO parameters)
  !>
  !> @details This routine computes the conditional probabilities of
  !>          one type of atom neighbouring another type of atom. These
  !>          are *not* the Warren-Cowley ASRO parameters, but you can
  !>          convert to them using a simple rescaling.
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  configuration Current atomic configuration
  !> @param  n_shells Number of shells on which to compute probabilities
  !> @param  shell_distances Array of lattice shell distances
  !>
  !> @return The computed radial densities (conditional probabilities)
  function radial_densities(setup, configuration, n_shells,            &
                            shell_distances) result(r_densities)
    type(run_params), intent(in) :: setup
    integer(array_int), dimension(:,:,:,:) :: configuration
    real(real64), dimension(:), allocatable :: shell_distances
    real(real64), dimension(setup%n_species,setup%n_species,           &
                            n_shells) :: r_densities
    integer, intent(in) :: n_shells
    integer :: i_1,i_2,i_3,j_1,j_2,j_3,j_b,jj_1,jj_2,jj_3, &
               l, species_i, species_j, i,j, i_b
    integer :: loop_1, loop_2, loop_3
    integer, dimension(setup%n_species) :: particle_counts
    real(real64) :: distance, d_x, d_y, d_z

    ! Array for counting the number of each species
    ! Initialise it to zero
    particle_counts = 0

    ! Output array of radial densities
    ! Initialise it to zero
    r_densities = 0.0_real64

    ! Count how many of each species there are
    do i_3=1, setup%n_3*2
      do i_2=1, setup%n_2*2
        do i_1=1, setup%n_1*2
          do i_b=1, setup%n_basis
            do l=1, setup%n_species
              if (configuration(i_b, i_1, i_2, i_3) .eq.              &
                                int(l, kind=array_int)) then
                particle_counts(l) = particle_counts(l) + 1
              end if
            end do
          end do
        end do
      end do
    end do

    ! Check that we won't divide by zero later
    do l=1, setup%n_species
      if (particle_counts(l) .eq. 0) then
        print*, 'Warning, one or more particle counts are zero in radial_densities()'
      end if
    end do

    ! For small systems, limit loop size to prevent double counting
    loop_1 = min(setup%n_1, 5)
    loop_2 = min(setup%n_2, 5)
    loop_3 = min(setup%n_3, 5)

    ! Loop over all lattice sites
    do i_3=1, 2*setup%n_3
      do i_2=1, 2*setup%n_2
        do i_1=1, 2*setup%n_1
          do i_b=1, setup%n_basis
          ! Cycle if this site is empty
          if (configuration(i_b, i_1, i_2, i_3) .eq. 0_array_int) cycle
            ! Loop over neighbouring sites, accounting for
            ! P.B.C.s
            do jj_3=i_3-loop_3, i_3+loop_3, 1
              j_3 = modulo(jj_3-1, 2*setup%n_3) + 1
              do jj_2=i_2-loop_2, i_2+loop_2, 1
                j_2 = modulo(jj_2-1, 2*setup%n_2) + 1
                do jj_1=i_1-loop_1, i_1+loop_1, 1
                  j_1 = modulo(jj_1-1, 2*setup%n_1) + 1
                  do j_b=1, setup%n_basis
                    if (configuration(j_b, j_1, j_2, j_3) .eq. 0_array_int) cycle
                    ! Compute the distance to this site, accounting
                    ! for PBCs
                    d_x = real(i_1-j_1)
                    d_y = real(i_2-j_2)
                    d_z = real(i_3-j_3)

                    d_x = d_x - float(2*setup%n_1)* &
                                nint(d_x/float(2*setup%n_1))
                    d_y = d_y - float(2*setup%n_2)* &
                                nint(d_y/float(2*setup%n_2))
                    d_z = d_z - float(2*setup%n_3)* &
                                nint(d_z/float(2*setup%n_3))

                    distance = sqrt(d_x**2 + d_y**2 + d_z**2)

                    ! Loop over and find which shell this sits in
                    do l=1, n_shells
                      if (abs(distance - shell_distances(l)) &
                          .lt. 1e-3_real64) then

                        ! Add it to the relevant entry for this shell
                        species_i = configuration(i_b,i_1, i_2, i_3)
                        species_j = configuration(j_b,j_1, j_2, j_3)
                        r_densities(species_i, species_j, l) = &
                          r_densities(species_i, species_j, l) + 1.0_real64 
                      end if
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    ! Nice nested do loop...

    ! Average them
    do i=1, n_shells
      do j=1, setup%n_species
        r_densities(j,:,i) = r_densities(j,:,i)/particle_counts(j)
      end do
    end do
    
  end function radial_densities


  !> @brief   Implementation of the quicksort algorithm for arrays
  !>
  !> @details Puts an array of real numbers into size order from
  !>          smallest to largest. Note that this operates *on* the
  !>          input array---if you would like to keep the unsorted
  !>          array, make a copy of it.
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2023
  !>
  !> @param  array Array of real (real64) numbers to sort
  !>
  !> @return None
  recursive subroutine quicksort(array)

    real(real64), intent(inout)::array(:)
    real(real64) :: temp,pivot
    integer :: i,j,last,left,right

    last=size(array)

    if (last.lt.50) then ! use insertion sort on small arrays
       do i=2,last
          temp=array(i)
          do j=i-1,1,-1
             if (array(j).le.temp) exit
             array(j+1)=array(j)
          enddo
          array(j+1)=temp
       enddo
       return
    endif
    ! find median of three pivot
    ! and place sentinels at first and last elements
    temp=array(last/2)
    array(last/2)=array(2)
    if (temp.gt.array(last)) then
       array(2)=array(last)
       array(last)=temp
    else
       array(2)=temp
    endif
    if (array(1).gt.array(last)) then
       temp=array(1)
       array(1)=array(last)
       array(last)=temp
    endif
    if (array(1).gt.array(2)) then
       temp=array(1)
       array(1)=array(2)
       array(2)=temp
    endif
    pivot=array(2)

    left=3
    right=last-1
    do
       do while(array(left).lt.pivot)
          left=left+1
       enddo
       do while(array(right).gt.pivot)
          right=right-1
       enddo
       if (left.ge.right) exit
       temp=array(left)
       array(left)=array(right)
       array(right)=temp
       left=left+1
       right=right-1
    enddo
    if (left.eq.right) left=left+1
    call quicksort(array(1:left-1))
    call quicksort(array(left:))

  end subroutine quicksort

end module analytics
