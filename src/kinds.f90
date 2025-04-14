!> @file    kinds.f90
!>
!> @brief   Standard KINDS for numeric data types
!>
!> @details Code originally written by H. Ratcliffe and C. S. Brady
!>          for the Research Software Engineering course delivered as
!>          part of the taught component of the EPSRC-supported Centre
!>          for Doctoral Training (CDT) in Modelling of Heterogeneous
!>          Systems (HetSys) at the University of Warwick. Re-used here
!>          with permission.
!>
!>          Original description (verbatim):
!>          This module provides the numeric kinds we think you might
!>          need, using names that match those in F2008. This makes it
!>          smooth to change to std F2008 by simply removing this
!>          module.
!>
!>          Do remember that not every length of type has to exist on a
!>          particular system - for those which don't you will get AT
!>          LEAST the requested length if this does not exceed the
!>          longest the system offers.
!>
!>
!> @author  H. Ratcliffe
!> @author  C. S. Brady
!> @author  C. D. Woodgate
!>
!> @date    2019-2025
module kinds

  implicit none

  !> Very short integer
  !> (8 bit, -128 to 127)
  integer, parameter :: int8 = selected_int_kind(2)

  !> Short integer
  !> (16 bit, -32 768 to 32 767)
  integer, parameter :: int16 = selected_int_kind(4)

  !> "Normal" integer
  !> (32 bit == 4 byte, -2 147 483 648 to 2 147 483 647)
  integer, parameter :: int32 = selected_int_kind(9)

  !> Long integer
  !> (64 bit, -9 223 372 036 854 775 808 to 9 223 372 036 854 775 807)
  integer, parameter :: int64 = selected_int_kind(15)

  !> Normal "float"
  !> (32 bit = 4 bytes, approx -3.4e38 to 3.4e38 and covering values
  !>  down to about 1e-38 magnitude)
  integer, parameter :: real32 = selected_real_kind(6, 37)

  !> Longer "double"
  !> (64 bit, approx -1.8e308 to 1.8e308 and covering values down to
  !>  about 2e-308 magnitude)
  integer, parameter :: real64 = selected_real_kind(15, 307)

  !> Extra long "real"
  !> (128 bit = 16 bytes). Few systems have numbers this long at all!
  integer, parameter :: real128 = selected_real_kind(33, 4931)

end module
