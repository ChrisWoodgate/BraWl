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
!> @author  Claude Fable 5.1
!> @date    2019-2026
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

  ! Half-width, in cells, of the window of neighbouring cells scanned
  ! by lattice_shells() and radial_densities() for the ideal hcp
  ! lattice. Offsets of up to five cells cover every coordination
  ! shell closer than about 4a, well beyond the ten shells for which
  ! energy kernels exist.
  integer, parameter :: ideal_hcp_asro_window = 5

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

    integer(array_int), dimension(:,:,:,:), intent(in) :: config
    type(run_params), intent(in) :: setup
    real(real64), dimension(:,:,:,:,:), intent(inout), allocatable :: averages
    integer :: i,j,k,l,m

    do m=1, setup%n_grid_3
      do l=1, setup%n_grid_2
        do k=1, setup%n_grid_1
          do j=1, setup%n_basis
            do i=1, setup%n_species
              if (config(j,k,l,m) == i) then
                averages(i,j,k,l,m) = averages(i,j,k,l,m) + 1.0_real64
              end if
            end do
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
      do k=1, setup%n_grid_3
        do j=1, setup%n_grid_2
          do i=1, setup%n_grid_1
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
  !> @author  Claude Fable 5.1
  !> @date    2019-2026
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
    integer :: i,j,k,b, n, my_rank

    if(my_rank == 0) then
      write(6,'(16("-"),x,"Checking contents of simulation cell(s)",x, 15("-"),/)')
    end if

    sizes = shape(config)

    allocate(species_count(setup%n_species))

    species_count = 0
    n=0

    do k=1, sizes(4)
      do j=1, sizes(3)
        do i=1, sizes(2)
          do b=1, sizes(1)
            if (config(b,i,j,k) .ne. 0_array_int) then
              n = n+1
              species_count(config(b,i,j,k)) = &
                species_count(config(b,i,j,k)) + 1
            end if
          end do
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
  !> @details For the cubic lattices, distances are measured in units
  !>          of the simulation grid spacing (half the conventional
  !>          lattice parameter) from the first grid point to every
  !>          occupied grid point. For ideal_hcp, real-space distances
  !>          in units of the lattice parameter are computed from the
  !>          lattice and basis vectors, from a basis-1 site to every
  !>          site within the same window of cell offsets scanned by
  !>          radial_densities(), so that the two routines agree by
  !>          construction. In both cases the first shell is the
  !>          on-site distance of zero.
  !>
  !> @author  C. D. Woodgate
  !> @author  Claude Fable 5.1
  !> @date    2019-2026
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  shells Array where shell distances will be stored
  !> @param  configuration Current atomic configuration
  !>
  !> @return None
  subroutine lattice_shells(setup, shells, configuration)

    integer(array_int), dimension(:,:,:,:) :: configuration
    type(run_params), intent(in) :: setup
    integer :: i,j,k,b,l
    integer :: loop_1, loop_2, loop_3
    real(real64) :: dist
    real(real64), dimension(:), allocatable :: all_shells, shells

    shells = 0.0_real64

    l = 1

    if (trim(setup%lattice) == 'ideal_hcp') then

      ! Window of cell offsets to scan (same as radial_densities)
      loop_1 = min(setup%n_grid_1, ideal_hcp_asro_window)
      loop_2 = min(setup%n_grid_2, ideal_hcp_asro_window)
      loop_3 = min(setup%n_grid_3, ideal_hcp_asro_window)

      allocate(all_shells(setup%n_basis*(2*loop_1+1)*(2*loop_2+1)      &
                          *(2*loop_3+1)))
      all_shells = 0.0_real64

      ! Real-space distances from a basis-1 site to every site in the
      ! window
      do k=-loop_3, loop_3
        do j=-loop_2, loop_2
          do i=-loop_1, loop_1
            do b=1, setup%n_basis
              all_shells(l) = basis_site_separation(setup, 1, b,       &
                                                    i, j, k)
              l=l+1
            end do
          end do
        end do
      end do

    else

      ! One entry for every point of the simulation grid (plus one)
      allocate(all_shells(setup%n_grid_1*setup%n_grid_2                &
                          *setup%n_grid_3*setup%n_basis+1))
      all_shells = 0.0_real64

      ! Loop over all lattice sites
      do k=1, setup%n_grid_3
        do j=1, setup%n_grid_2
          do i=1, setup%n_grid_1
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

    end if

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
  !>          convert to them using a simple rescaling. Pairs are
  !>          assigned to shells by distance: in units of the grid
  !>          spacing for the cubic lattices, and in real space from the
  !>          lattice and basis vectors for ideal_hcp (see
  !>          lattice_shells()).
  !>
  !> @author  C. D. Woodgate
  !> @author  Claude Fable 5.1
  !> @date    2019-2026
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
    logical :: hcp

    ! Does this lattice use real-space distances?
    hcp = (trim(setup%lattice) == 'ideal_hcp')

    ! Array for counting the number of each species
    ! Initialise it to zero
    particle_counts = 0

    ! Output array of radial densities
    ! Initialise it to zero
    r_densities = 0.0_real64

    ! Count how many of each species there are
    do i_3=1, setup%n_grid_3
      do i_2=1, setup%n_grid_2
        do i_1=1, setup%n_grid_1
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
    if (hcp) then
      loop_1 = min(setup%n_grid_1, ideal_hcp_asro_window)
      loop_2 = min(setup%n_grid_2, ideal_hcp_asro_window)
      loop_3 = min(setup%n_grid_3, ideal_hcp_asro_window)
    else
      loop_1 = min(setup%n_1, 5)
      loop_2 = min(setup%n_2, 5)
      loop_3 = min(setup%n_3, 5)
    end if

    ! Loop over all lattice sites
    do i_3=1, setup%n_grid_3
      do i_2=1, setup%n_grid_2
        do i_1=1, setup%n_grid_1
          do i_b=1, setup%n_basis
          ! Cycle if this site is empty
          if (configuration(i_b, i_1, i_2, i_3) .eq. 0_array_int) cycle
            ! Loop over neighbouring sites, accounting for
            ! P.B.C.s
            do jj_3=i_3-loop_3, i_3+loop_3, 1
              j_3 = modulo(jj_3-1, setup%n_grid_3) + 1
              do jj_2=i_2-loop_2, i_2+loop_2, 1
                j_2 = modulo(jj_2-1, setup%n_grid_2) + 1
                do jj_1=i_1-loop_1, i_1+loop_1, 1
                  j_1 = modulo(jj_1-1, setup%n_grid_1) + 1
                  do j_b=1, setup%n_basis
                    if (configuration(j_b, j_1, j_2, j_3) .eq. 0_array_int) cycle
                    if (hcp) then
                      ! Real-space distance from the (unwrapped) cell
                      ! offset and the basis positions
                      distance = basis_site_separation(setup, i_b,     &
                                   j_b, jj_1-i_1, jj_2-i_2, jj_3-i_3)
                    else
                      ! Compute the distance to this site, accounting
                      ! for PBCs
                      d_x = real(i_1-j_1)
                      d_y = real(i_2-j_2)
                      d_z = real(i_3-j_3)

                      d_x = d_x - float(setup%n_grid_1)* &
                                  nint(d_x/float(setup%n_grid_1))
                      d_y = d_y - float(setup%n_grid_2)* &
                                  nint(d_y/float(setup%n_grid_2))
                      d_z = d_z - float(setup%n_grid_3)* &
                                  nint(d_z/float(setup%n_grid_3))

                      distance = sqrt(d_x**2 + d_y**2 + d_z**2)
                    end if

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


  !> @brief   Function to compute the real-space separation of two
  !>          sites of a lattice with a multi-atom basis
  !>
  !> @details The second site is specified relative to the first by the
  !>          offset between their cells, which is used as given (no
  !>          periodic wrapping), together with the basis index of each
  !>          site. The separation is formed from the lattice and
  !>          basis vectors stored in setup.
  !>
  !> @author  C. D. Woodgate
  !> @author  Claude Fable 5.1
  !> @date    2026
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  b_i Basis index of the first site
  !> @param  b_j Basis index of the second site
  !> @param  d_1 Cell offset of the second site along a_1
  !> @param  d_2 Cell offset of the second site along a_2
  !> @param  d_3 Cell offset of the second site along a_3
  !>
  !> @return The separation in units of the lattice parameter
  pure function basis_site_separation(setup, b_i, b_j, d_1, d_2,      &
                                      d_3) result(dist)

    type(run_params), intent(in) :: setup
    integer, intent(in) :: b_i, b_j, d_1, d_2, d_3
    real(real64) :: dist
    real(real64), dimension(3) :: r_vec

    r_vec = real(d_1, real64)*setup%lattice_vectors(:,1)               &
          + real(d_2, real64)*setup%lattice_vectors(:,2)               &
          + real(d_3, real64)*setup%lattice_vectors(:,3)               &
          + setup%basis_vectors(:,b_j) - setup%basis_vectors(:,b_i)

    dist = norm2(r_vec)

  end function basis_site_separation

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
