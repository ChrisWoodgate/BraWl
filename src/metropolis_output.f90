!----------------------------------------------------------------------!
! write_diagnostics.f90                                                !
!                                                                      !
! Module for writing energy data                                       !
!                                                                      !
! C. D. Woodgate,  Warwick                                        2023 !
!----------------------------------------------------------------------!
module write_diagnostics

  use kinds
  use shared_data
  use analytics

  implicit none

  contains

  !--------------------------------------------------------------------!
  ! Routine to write energy and diagnostic data                        !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
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
   
end module write_diagnostics
