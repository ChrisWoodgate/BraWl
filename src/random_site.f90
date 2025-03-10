!----------------------------------------------------------------------!
! random_site.f90                                                      !
!                                                                      !
! Module for getting random lattice site and neighbour.                !
!                                                                      !
! C. D. Woodgate,  Bristol                                        2025 !
!----------------------------------------------------------------------!
module random_site

  use kinds
  use shared_data
  use c_functions
  
  implicit none

  ! Array for neighbours on the sc lattice
  integer, parameter, dimension(3,6) :: &
  sc_nbrs = reshape((/  0,  0,  1, &
                        0,  1,  0, &
                        1,  0,  0, &
                        0,  0, -1, &
                        0, -1,  0, &
                       -1,  0,  0 /), (/3, 6/))

  ! Array for neighbours on the bcc lattice
  integer, parameter, dimension(3,8) :: &
  bcc_nbrs = reshape((/  1,  1,  1, &
                         1,  1, -1, &
                         1, -1,  1, &
                         1, -1, -1, &
                        -1,  1,  1, &
                        -1,  1, -1, &
                        -1, -1,  1, &
                        -1, -1, -1 /), (/3, 8/))

  ! Array for neighbours on the fcc lattice
  integer, parameter, dimension(3,12) :: &
  fcc_nbrs = reshape((/  0,  1,  1, &
                         0,  1, -1, &
                         0, -1,  1, &
                         0, -1, -1, &
                         1,  1,  0, &
                         1, -1,  0, &
                         1,  0,  1, &
                         1,  0, -1, &
                        -1,  1,  0, &
                        -1, -1,  0, &
                        -1,  0,  1, &
                        -1,  0, -1 /), (/3, 12/))
  contains

  !-------------------------------------------------!
  ! Function to get a random site on the sc lattice !
  !-------------------------------------------------!
  pure function simple_cubic_random_site(setup) result(site)
    class(run_params), intent(in) :: setup
    integer, dimension(4) :: site
    site(1) = 1
    site(2) = floor(genrand()*2.0_real64*real(setup%n_1)) +1
    site(3) = floor(genrand()*2.0_real64*real(setup%n_2)) +1
    site(4) = floor(genrand()*2.0_real64*real(setup%n_3)) +1
  end function simple_cubic_random_site
  
  !------------------------------------------------------!
  ! Function to get a random neighbour on the sc lattice !
  !------------------------------------------------------!
  pure function simple_cubic_random_nbr(setup, site) result(nbr)
    class(run_params), intent(in) :: setup
    integer, dimension(4), intent(in) :: site
    integer, dimension(4) :: nbr 
    integer :: n

    ! Generate random integer corresponding to neighbour
    n = floor(6.0_real64*genrand())+1

    ! Make this the neighbour
    nbr(2:) = site(2:) + sc_nbrs(:,n)

    ! Wrap coordinates into box
    nbr(1) = 1
    nbr(2) = modulo(nbr(2)-1, 2*setup%n_1) + 1
    nbr(3) = modulo(nbr(3)-1, 2*setup%n_2) + 1
    nbr(4) = modulo(nbr(4)-1, 2*setup%n_3) + 1
  end function simple_cubic_random_nbr
  
  !--------------------------------------------------!
  ! Function to get a random site on the bcc lattice !
  !--------------------------------------------------!
  pure function bcc_random_site(setup) result(site)
    class(run_params), intent(in) :: setup
    integer, dimension(4) :: site
    site(1) = 1
    site(4) = floor(2.0_real64*genrand()*real(setup%n_3, real64)) +1
    site(2) = 2 * floor(genrand()*real(setup%n_1, real64)) &
                + 2 - modulo(site(4), 2)
    site(3) = 2 * floor(genrand()*real(setup%n_2, real64)) &
                + 2 - modulo(site(4), 2)
  end function bcc_random_site
  
  !-------------------------------------------------------!
  ! Function to get a random neighbour on the bcc lattice !
  !-------------------------------------------------------!
  pure function bcc_random_nbr(setup, site) result(nbr)
    class(run_params), intent(in) :: setup
    integer, dimension(4), intent(in) :: site
    integer, dimension(4) :: nbr 
    integer :: n

    ! Generate random integer corresponding to neighbour
    n = floor(8.0_real64*genrand())+1

    ! Make this the neighbour
    nbr(2:) = site(2:) + bcc_nbrs(:,n)

    ! Wrap coordinates into box
    nbr(1) = 1
    nbr(2) = modulo(nbr(2)-1, 2*setup%n_1) + 1
    nbr(3) = modulo(nbr(3)-1, 2*setup%n_2) + 1
    nbr(4) = modulo(nbr(4)-1, 2*setup%n_3) + 1
  end function bcc_random_nbr
  
  !--------------------------------------------------!
  ! Function to get a random site on the fcc lattice !
  !--------------------------------------------------!
  pure function fcc_random_site(setup) result(site)
    class(run_params), intent(in) :: setup
    integer, dimension(4) :: site
    site(1) = 1
    site(4) = floor(2.0_real64*genrand()*real(setup%n_3, real64)) + 1
    site(2) = floor(2.0_real64*genrand()*real(setup%n_1, real64)) + 1
    site(3) = 2 * floor(genrand()*real(setup%n_2, real64)) + 1 &
                    + modulo((site(2)-modulo(site(4),2)), 2)
  end function fcc_random_site

  !-------------------------------------------------------!
  ! Function to get a random neighbour on the fcc lattice !
  !-------------------------------------------------------!
  pure function fcc_random_nbr(setup, site) result(nbr)
    class(run_params), intent(in) :: setup
    integer, dimension(4), intent(in) :: site
    integer, dimension(4) :: nbr 
    integer :: n

    ! Generate random integer corresponding to neighbour
    n = floor(12.0_real64*genrand())+1

    ! Make this the neighbour
    nbr(2:) = site(2:) + fcc_nbrs(:,n)

    ! Wrap coordinates into box
    nbr(1) = 1
    nbr(2) = modulo(nbr(2)-1, 2*setup%n_1) + 1
    nbr(3) = modulo(nbr(3)-1, 2*setup%n_2) + 1
    nbr(4) = modulo(nbr(4)-1, 2*setup%n_3) + 1
  end function fcc_random_nbr

  !--------------------------------------------------------------------!
  ! Function to swap a pair of atoms on the lattice.                   !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  subroutine pair_swap(config, idx1, idx2)
    !integer(int16), allocatable, dimension(:,:,:,:) :: config
    integer(int16), dimension(:,:,:,:) :: config
    integer, dimension(4), intent(in) :: idx1, idx2
    integer(int16) :: species1, species2

    species1 = config(idx1(1), idx1(2), idx1(3), idx1(4))
    species2 = config(idx2(1), idx2(2), idx2(3), idx2(4))

    config(idx1(1), idx1(2), idx1(3), idx1(4)) = species2
    config(idx2(1), idx2(2), idx2(3), idx2(4)) = species1

  end subroutine pair_swap
  
end module random_site
