!> @file    c_functions.f90
!>
!> @brief   Interfaces for routines written in C
!>
!> @details This module contains interfaces for routines written in C.
!>          Currently this is just the Mersenne Twister PNRG, but other
!>          interfaces can be added in future as needed.
!>
!> @author  C. D. Woodgate
!> @date    2019-2025
module c_functions

  use iso_c_binding

  implicit none

  private

  public :: genrand, f90_init_genrand

  interface
    ! Mersenne Twister PNRG
    ! Needs to have seed set with f90_init_genrand()
    pure function genrand() bind(C, name='genrand')
      use, intrinsic :: iso_c_binding, only : C_double
      real(C_double) :: genrand
    end function genrand

    ! Routine to set the seed for genrand
    function f90_init_genrand(seedtime, my_rank, job_id) bind(C, name='f90_init_genrand')
      use, intrinsic :: iso_c_binding, only : C_int
      integer(C_int) :: f90_init_genrand
      integer(C_int), value :: my_rank
      integer(C_int), value :: seedtime
      integer(C_int), value :: job_id
    end function f90_init_genrand
  end interface

end module c_functions
