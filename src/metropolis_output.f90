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

  public :: diagnostics_writer

  contains

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
