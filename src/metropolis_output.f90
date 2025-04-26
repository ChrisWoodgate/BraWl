!> @file    metropolis_output.f90
!>
!> @brief   Routines for writing plain text outputs from a Metropolis
!>          simulation
!>
!> @author  C. D. Woodgate
!> @date    2019-2025
module metropolis_output

  use kinds
  use shared_data
  use analytics

  implicit none

  private

  public :: diagnostics_writer, energy_trajectory_writer,              &
            asro_trajectory_writer

  contains

  !> @brief   Subroutine to write energy of a Metropolis simulation
  !>          trajectory to a plain text file
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  filename Name of file to which to write
  !> @param  step_number Index of current MC trial step
  !> @param  energy Simulation energy to write
  !>
  !> @return None
  subroutine energy_trajectory_writer(filename, step_number, energy)

    character(len=*), intent(in) :: filename
    integer :: step_number
    real(real64) :: energy
    logical :: exist

    inquire(file=filename, exist=exist)
    if (exist) then
      open(7, file=filename, status="old", position="append", action="write")
    else
      open(7, file=filename, status="new", action="write")
      write(7, *) '# step_number E'
    end if

    write(7,"(I13,x,f20.10,x)") step_number, energy

    close(7)

  end subroutine energy_trajectory_writer

  !> @brief   Subroutine to write ASRO of a Metropolis simulation
  !>          trajectory to a plain text file
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  filename Name of file to which to write
  !> @param  step_number Index of current MC trial step
  !> @param  energy Simulation energy to write
  !>
  !> @return None
  subroutine asro_trajectory_writer(filename, step_number, asro)

    character(len=*), intent(in) :: filename
    integer :: step_number, i, j, k
    integer, dimension(3) :: asro_size
    real(real64), dimension(:,:,:) :: asro
    logical :: exist

    asro_size = shape(asro)

    inquire(file=filename, exist=exist)
    if (exist) then
      open(7, file=filename, status="old", position="append", action="write")
    else
      open(7, file=filename, status="new", action="write")
      write(7, *) '# step_number ASRO'
    end if

    write(7,"(I13,x,f20.10,x)", advance='no') step_number

    do k=1, asro_size(3)
      do j=1, asro_size(2)
        do i=1, asro_size(1)
          write(7,"(f8.5,x)", advance='no') asro(i,j,k)
        end do
      end do
    end do

    write(7,"(x)", advance='yes')

    close(7)

  end subroutine asro_trajectory_writer

  !> @brief   Subroutine to write outputs of a Metropolis simulation to
  !>          plain text file
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2025
  !>
  !> @param  filename Name of file to which to write
  !> @param  temps Array of simulation temperatures
  !> @param  energies Array of simulation energies
  !> @param  C Array of simulation heat capacities
  !> @param  acceptance Array of simulation heat capacities
  !>
  !> @return None
  subroutine diagnostics_writer(filename, temps, energies, C, acceptance)

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    real(real64), allocatable, dimension(:), intent(in) :: &
                            temps, energies, C, acceptance

    integer :: i, n_steps

    n_steps = size(energies)

    open(unit=7, file=filename)

    write(7, *) 'T, E, C, Acceptance Rate'

    do i=1, n_steps
      write(7,*) temps(i), ', ', energies(i), ', ', C(i), &
                 ', ', acceptance(i)
    end do

    close(7)

  end subroutine diagnostics_writer
   
end module metropolis_output
