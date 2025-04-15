!> @file    constants.f90
!>
!> @brief   Module containing relevant constants.
!>
!> @details 
!>
!> @author  C. D. Woodgate
!> @date    2025
module constants

  use kinds

  implicit none

  private

  public :: k_b_in_eV, k_b_in_Ry, Ry_to_eV

  !--------------------------------------------------------------------!
  ! Hard-coded physical constants. Mainly for conversion of internal   !
  ! units. (Internally the code works with Rydbergs.)                  !
  !                                                                    !
  ! Currently taken from the CODATA recommended values 2022            !
  ! https://physics.nist.gov/cuu/Constants/index.html                  !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!

  ! 2022 CODATA Value for k_b in eV/K
  real(real64), parameter :: k_b_in_eV &
                             =8.167333262e-5_real64

  ! 2022 CODATA Value for k_b in eV/K 
  real(real64), parameter :: k_b_in_Ry &
                             =8.167333262e-5_real64/13.605693122990_real64

  ! 2022 CODATA value for Rydberg in eV
  real(real64), parameter :: Ry_to_eV = 13.605693122_real64

end module constants
