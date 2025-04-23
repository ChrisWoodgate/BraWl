!> @file    write_xyz.f90
!>
!> @brief   Module for writing xyz file                             
!>
!> @details This module contains a routine to produce a configuration file in xyz format 
!>
!> @author  C. D. Woodgate
!> @date    2019-2025
! C. D. Woodgate,  Warwick                                        2023 !
module write_xyz

  use kinds
  use constants
  use shared_data
  use analytics

  implicit none

  contains

  !> @brief   Subroutine to write xyz file                         
  !>
  !> @details This is routine writes an xyz file of the current configuration, in the extended xyz format, 
  !>          (incl. lattice, element type and coordinates) with option to append configurations in a 
  !>          single trajectory file.
  !>
  !> @author  C. D. Woodgate & L. B. Partay
  !> @date    2024
  !>
  !> @param  filename Name of file to be created                                          
  !> @param  configuration Current atomic configuration
  !> @param  setup Derived type containing simulation parameters
  !> @param  trajectory Logical type, if true, configuration will be appended to the end of the file
  !>
  !> @return None
  subroutine xyz_writer(filename, configuration, setup, trajectory)
    ! Simulation setup information
    type(run_params), intent(in) :: setup

    ! Data to write to file
    !integer(array_int), dimension(:,:,:,:), allocatable, intent(in) :: configuration
    integer(array_int), dimension(:,:,:,:), intent(in) :: configuration

    ! Option to open existing xyz file and append new configuration to the end
    logical, optional :: trajectory

    ! Sizes of grid
    integer, dimension(4) :: grid_sizes

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    integer :: n_particles, i, j, k, l
    logical :: exist
    real(real64), dimension(3) :: pos

    n_particles = total_particle_count(setup, configuration)

    grid_sizes = shape(configuration)

    ! Check if we want to open new xyz, overwrite existing xyz or append existing xyz file
    if (present(trajectory)) then

      if (trajectory) then
        inquire(file=filename, exist=exist)
        if (exist) then
          open(7, file=filename, status="old", position="append", action="write")
        else
          open(7, file=filename, status="new", action="write")
        end if      
      end if

    else ! open xyz (if exists, will ovrewrite)

      open(unit=7, file=filename)

    endif  

    write(7, *) n_particles

    write(7,*) 'Lattice="',setup%lattice_parameter*setup%n_1, 0.0, 0.0 , &
                   &   0.0, setup%lattice_parameter*setup%n_2, 0.0 , &
                   &   0.0, 0.0, setup%lattice_parameter*setup%n_3,     '"'

!    if (present(trajectory)) then
!      if (trajectory) then
!        write(7,*) 'Lattice="',setup%lattice_parameter*setup%n_1, 0.0, 0.0 , &
!                       &   0.0, setup%lattice_parameter*setup%n_2, 0.0 , &
!                       &   0.0, 0.0, setup%lattice_parameter*setup%n_3,     '"'
!      end if          
!    else           
!      write(7,*) ''
!    end if

    do i=1, grid_sizes(1)
      do j=1, grid_sizes(2)
        do k=1, grid_sizes(3)
          do l=1, grid_sizes(4)
            if (configuration(i,j,k,l) .eq. 0) cycle
            pos = real(j-1)*setup%lattice_vectors(1,:) + &
                  real(k-1)*setup%lattice_vectors(2,:) + &
                  real(l-1)*setup%lattice_vectors(3,:) + &
                  real(i)*setup%basis_vectors
            pos = pos*setup%lattice_parameter
            write(7,*) setup%species_names(configuration(i,j,k,l)), pos(1),pos(2),pos(3)
          end do
        end do
      end do
    end do

    close(unit=7)

  end subroutine xyz_writer
   
end module write_xyz
