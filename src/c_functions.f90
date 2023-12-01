!----------------------------------------------------------------------!
! c_functions.f90                                                      !
!                                                                      !
! Module to interface with Mersenne Twiseter PNRG. (Used to also       !
! interface with some other routines that have now been ported to      !
! Fortran                                                              !
!                                                                      !
! C. D. Woodgate,  Warwick                                        2023 !
!----------------------------------------------------------------------!
module c_functions
  use iso_c_binding

  implicit none

  interface
    ! Mersenne Twister PNRG
    ! Needs to have seed set with f90_init_genrand()
    pure function genrand() bind(C, name='genrand')
      use, intrinsic :: iso_c_binding, only : C_double
      real(C_double) :: genrand
    end function genrand

    ! Routine to set the seed for genrand
    function f90_init_genrand(seedtime, my_rank) bind(C, name='f90_init_genrand')
      use, intrinsic :: iso_c_binding, only : C_int
      integer(C_int) :: f90_init_genrand
      integer(C_int), value :: my_rank
      integer(C_int), value :: seedtime
    end function f90_init_genrand
  end interface

end module c_functions
