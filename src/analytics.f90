!----------------------------------------------------------------------!
! analytics.f90                                                        !
!                                                                      !
! Various routines and tools for analysing the simulation              !
!                                                                      !
! C. D. Woodgate,  Warwick                                        2023 !
!----------------------------------------------------------------------!
module analytics

  use kinds
  use shared_data
  use io
  use display
  
  implicit none

  contains

  !--------------------------------------------------------------------!
  ! Routine to add this state to the average (for LRO)                 !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine store_state(averages, config, setup)
    integer(int16), allocatable, dimension(:,:,:), intent(in) :: config
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

  !--------------------------------------------------------------------!
  ! Routine to compute the average (for LRO)                           !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine average_state(averages, setup, n_steps)
    type(run_params), intent(in) :: setup
    integer, intent(in) :: n_steps
    real(real64), dimension(:,:,:,:), intent(inout), allocatable :: averages
    integer :: i

    do i=1, setup%n_species
      averages(i,:,:,:) = (1.0/real(n_steps, real64))*averages(i,:,:,:)
    end do
  end subroutine average_state
  
  !--------------------------------------------------------------------!
  ! Routine to count the number of particles (used for testing)        !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  function total_particle_count(setup, config) result(total_count)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    type(run_params) :: setup
    integer :: total_count
    integer :: i,j,k,b

    total_count = 0

    do b=1, setup%n_basis
      do k=1, setup%n_3
        do j=1, setup%n_2
          do i=1, setup%n_1
            if (config(b,i,j,k) .ne. 0_int16) then
              total_count = total_count + 1
            end if
          end do
        end do
      end do
    end do

  end function total_particle_count

  !--------------------------------------------------------------------!
  ! Routine to print the number of particles of each species (used for !
  ! testing)                                                           !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine print_particle_count(setup, config)
    integer(int16), allocatable, dimension(:,:,:), intent(in) :: config
    type(run_params) :: setup
    integer, dimension(3) :: sizes
    integer, dimension(:), allocatable :: species_count
    integer :: i,j,k, n

    sizes = shape(config)

    allocate(species_count(setup%n_species))

    species_count = 0
    n=0

    do k=1, sizes(3)
      do j=1, sizes(2)
        do i=1, sizes(1)
          if (config(i,j,k) .ne. 0_int16) then
            n = n+1
            species_count(config(i,j,k)) = &
              species_count(config(i,j,k)) + 1
          end if
        end do
      end do
    end do

    print*, 'Particle counts are: '

    do i=1, setup%n_species-1
      print*, 'Species ', i ,species_count(i)
    end do

    print*, 'Species ', setup%n_species,                   &
             species_count(setup%n_species), new_line('a')
    
    deallocate(species_count)

  end subroutine print_particle_count

  !--------------------------------------------------------------------!
  ! Subroutine to compute on-shell distances. Could just calculate     !
  ! these once, but this is good for a variety of lattices.            !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine lattice_shells(setup, shells, configuration)
    integer(int16), dimension(:,:,:,:), allocatable :: configuration
    type(run_params), intent(in) :: setup
    integer :: i,j,k,b,l
    real(real64) :: dist
    real(real64), dimension(3) :: r_vec
    real(real64), dimension(:), allocatable :: all_shells, shells

    ! Factor of eight to account for the fact that simulation
    ! doubles number of cells in each direction to build lattice
    allocate(all_shells(8*setup%n_1*setup%n_2*setup%n_3*setup%n_basis))

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
            if (configuration(b,i,j,k) .eq. 0_int16) cycle
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
    do i=1, size(all_shells)
      if (abs(all_shells(i)-all_shells(i+1)) .lt. 1e-3_real64) cycle
      shells(l) = all_shells(i)
      l=l+1
      if (l .gt. setup%wc_range) exit
    end do

    ! Deallocate the array of all distances
    deallocate(all_shells)

  end subroutine lattice_shells

  !--------------------------------------------------------------------!
  ! Subroutine to compute radial densities, i.e. atomic short-range    !
  ! order parameters.                                                  !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine radial_densities(setup, configuration, n_shells,   &
                              shell_distances, r_densities, T)
    type(run_params), intent(in) :: setup
    integer(int16), dimension(:,:,:,:), allocatable :: configuration
    real(real64), dimension(:), allocatable :: shell_distances
    real(real64), dimension(:,:,:,:), allocatable :: r_densities
    integer, intent(in) :: n_shells, T
    integer :: i_1,i_2,i_3,j_1,j_2,j_3,j_b,jj_1,jj_2,jj_3, &
               l, species_i, species_j, i,j, i_b
    integer, dimension(setup%n_species) :: particle_counts
    real(real64) :: distance, d_x, d_y, d_z

    ! Array for counting the number of each species
    particle_counts = 0

    ! Count how many of each species there are
    do i_3=1, setup%n_3*2
      do i_2=1, setup%n_2*2
        do i_1=1, setup%n_1*2
          do i_b=1, setup%n_basis
            do l=1, setup%n_species
              if (configuration(i_b, i_1, i_2, i_3) .eq. int(l, kind=int16)) then
                particle_counts(l) = particle_counts(l) + 1
              end if
            end do
          end do
        end do
      end do
    end do

    ! Loop over all lattice sites
    do i_3=1, 2*setup%n_3
      do i_2=1, 2*setup%n_2
        do i_1=1, 2*setup%n_1
          do i_b=1, setup%n_basis
          ! Cycle if this site is empty
          if (configuration(i_b, i_1, i_2, i_3) .eq. 0_int16) cycle
            ! Loop over neighbouring sites, accounting for
            ! P.B.C.s
            do jj_3=i_3-4, i_3+4, 1
              j_3 = modulo(jj_3-1, 2*setup%n_3) + 1
              do jj_2=i_2-4, i_2+4, 1
                j_2 = modulo(jj_2-1, 2*setup%n_2) + 1
                do jj_1=i_1-4, i_1+4, 1
                  j_1 = modulo(jj_1-1, 2*setup%n_1) + 1
                  do j_b=1, setup%n_basis
                    if (configuration(j_b, j_1, j_2, j_3) .eq. 0_int16) cycle
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
                        r_densities(species_i, species_j, l, T) = &
                          r_densities(species_i, species_j, l, T) + 1.0_real64 
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
        r_densities(j,:,i, T) = r_densities(j,:,i, T)/particle_counts(j)
      end do
    end do
    
  end subroutine radial_densities


  !--------------------------------------------------------------------!
  ! Quicksort routine.                                                 !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
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
