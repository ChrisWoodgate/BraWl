!> @file    shared_data.f90
!>
!> @brief   Module containing important shared data.
!>
!> @details 
!>
!> @author  C. D. Woodgate
!> @date    2019-2025
module shared_data

  use kinds
  use constants
  
  implicit none

  save

  ! Random number seed
  integer(8) :: seed

  ! Arrays that are only used by rank 1
  real(real64), allocatable, dimension(:) :: av_energies_of_T,        &
                                             av_C_of_T,               &
                                             av_acceptance_of_T
  real(real64), allocatable, dimension(:,:,:,:) :: av_rho_of_T
  real(real64), allocatable, dimension(:,:,:,:,:,:) :: av_order_of_T

  ! Arrays used on all processors
  ! Indices run (basis, x, y, z)
  integer(array_int), allocatable, dimension(:,:,:,:) :: config
  real(real64), allocatable, dimension(:) :: energies_of_T, C_of_T,   &
                                             acceptance_of_T,         &
                                             temperature
  real(real64), allocatable, dimension(:,:,:,:,:,:) :: order_of_T
  real(real64), allocatable, dimension(:,:,:,:) :: rho_of_T
  real(real64), dimension(:,:,:,:), allocatable :: order
  real(real64), dimension(:,:,:), allocatable :: V_ex
  real(real64), allocatable, dimension(:) :: shells
  
end module shared_data
