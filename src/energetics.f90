!----------------------------------------------------------------------!
! energetics.f90                                                       !
!                                                                      !
! Module implementing the Bragg-Williams Hamiltonian.                  !
!                                                                      !
! C. D. Woodgate,  Bristol                                        2025 !
!----------------------------------------------------------------------!
module energetics

  use kinds
  use shared_data
  use io
  use c_functions
  
  implicit none

  contains

  !--------------------------------------------------------------------!
  ! Function to compute the total energy of the simulation             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function total_energy(setup,config) result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer :: i, j, k, l

    energy=0.0_real64

    do l=1, 2*setup%n_3
      do k=1, 2*setup%n_2
        do j=1, 2*setup%n_1
          do i=1, setup%n_basis
            if (config(i,j,k,l) .eq. 0_int16) cycle
            energy = energy + setup%nbr_energy(config, i, j, k, l)
          end do
        end do
      end do
    end do

    ! Divide by two to correct double-counting
    energy = 0.5_real64*energy

  end function total_energy

  !--------------------------------------------------------------------!
  ! Function to compute the energy associated with just one pair of    !
  ! atoms. (This means we can avoid re-evaluating the full Hamiltonian !
  ! for things like trial Metropolis-Kawasaki swaps                    !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function pair_energy(setup, config, idx1, idx2)&
       result(energy)
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    type(run_params), intent(in) :: setup
    integer, dimension(4), intent(in) :: idx1, idx2
    real(real64) :: energy

    ! The function nbr_energy will only evaluate the portion of the
    ! energy coming from the atoms around a site interacting with the
    ! atom on that site. There is NO factor of 1/2 because there is no
    ! double-counting in this routine.
    energy = setup%nbr_energy(config, idx1(1), idx1(2), idx1(3), idx1(4)) &
           + setup%nbr_energy(config, idx2(1), idx2(2), idx2(3), idx2(4))
  end function pair_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 1st coordination     !
  ! shell to the energy for the BCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_shell1_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, ip1, im1, jp1, jm1, kp1, km1, ib

    energy=0.0_real64
    
    ! Compute where my neighbours are
    ip1 = modulo(site_i, 2*setup%n_1) + 1
    im1 = modulo(site_i-2, 2*setup%n_1) + 1
    jp1 = modulo(site_j, 2*setup%n_2) + 1
    jm1 = modulo(site_j-2, 2*setup%n_2) + 1
    kp1 = modulo(site_k, 2*setup%n_3) + 1
    km1 = modulo(site_k-2, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b
      
    allocate(nbrs(8))
    nbrs(1) = config(ib, ip1, jp1, kp1)
    nbrs(2) = config(ib, ip1, jp1, km1)
    nbrs(3) = config(ib, ip1, jm1, kp1)
    nbrs(4) = config(ib, ip1, jm1, km1)
    nbrs(5) = config(ib, im1, jp1, kp1)
    nbrs(6) = config(ib, im1, jp1, km1)
    nbrs(7) = config(ib, im1, jm1, kp1)
    nbrs(8) = config(ib, im1, jm1, km1)
    do i=1, 8
      energy = energy + V_ex(species, nbrs(i), 1)
    end do
    deallocate(nbrs)
  end function bcc_shell1_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 2nd coordination     !
  ! shell to the energy for the BCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_shell2_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, ip2, im2, jp2, jm2, kp2, km2, ib

    energy=0.0_real64
    
    ! Compute where my neighbours are
    ip2 = modulo(site_i+1, 2*setup%n_1) + 1
    im2 = modulo(site_i-3, 2*setup%n_1) + 1
    jp2 = modulo(site_j+1, 2*setup%n_2) + 1
    jm2 = modulo(site_j-3, 2*setup%n_2) + 1
    kp2 = modulo(site_k+1, 2*setup%n_3) + 1
    km2 = modulo(site_k-3, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b
      
    allocate(nbrs(6))
    nbrs(1) = config(ib, ip2, site_j  , site_k  )
    nbrs(2) = config(ib, im2, site_j  , site_k  )
    nbrs(3) = config(ib, site_i  , jm2, site_k  )
    nbrs(4) = config(ib, site_i  , jp2, site_k  )
    nbrs(5) = config(ib, site_i  , site_j  , kp2)
    nbrs(6) = config(ib, site_i  , site_j  , km2)
    do i=1, 6
      energy = energy + V_ex(species, nbrs(i), 2)
    end do
    deallocate(nbrs)
  end function bcc_shell2_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 3rd coordination     !
  ! shell to the energy for the BCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_shell3_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, ip2, im2, jp2, jm2, kp2, km2, ib

    energy=0.0_real64
    
    ! Compute where my neighbours are
    ip2 = modulo(site_i+1, 2*setup%n_1) + 1
    im2 = modulo(site_i-3, 2*setup%n_1) + 1
    jp2 = modulo(site_j+1, 2*setup%n_2) + 1
    jm2 = modulo(site_j-3, 2*setup%n_2) + 1
    kp2 = modulo(site_k+1, 2*setup%n_3) + 1
    km2 = modulo(site_k-3, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b
      
    allocate(nbrs(12))
    nbrs(1)  = config(ib,site_i,  jm2,  km2)
    nbrs(2)  = config(ib, ip2, site_j,  km2)
    nbrs(3)  = config(ib, im2, site_j,  km2)
    nbrs(4)  = config(ib,site_i,  jp2,  km2)
    nbrs(5)  = config(ib, ip2,  jm2, site_k)
    nbrs(6)  = config(ib, im2,  jm2, site_k)
    nbrs(7)  = config(ib, ip2,  jp2, site_k)
    nbrs(8)  = config(ib, im2,  jp2, site_k)
    nbrs(9)  = config(ib,site_i,  jm2,  kp2)
    nbrs(10) = config(ib, ip2, site_j,  kp2)
    nbrs(11) = config(ib, im2, site_j,  kp2)
    nbrs(12) = config(ib,site_i,  jp2,  kp2)
    do i=1, 12
      energy = energy + V_ex(species, nbrs(i), 3)
    end do
    deallocate(nbrs)
  end function bcc_shell3_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 4th coordination     !
  ! shell to the energy for the BCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_shell4_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, up, dn, fw, bw, lt, rt, upupup, dndndn, &
               fwfwfw, bwbwbw, ltltlt, rtrtrt, ib

    energy=0.0_real64
    
    up = modulo(site_i, 2*setup%n_1) + 1
    dn = modulo(site_i-2, 2*setup%n_1) + 1
    lt = modulo(site_j, 2*setup%n_2) + 1
    rt = modulo(site_j-2, 2*setup%n_2) + 1
    fw = modulo(site_k, 2*setup%n_3) + 1
    bw = modulo(site_k-2, 2*setup%n_3) + 1
    up = modulo(site_i, 2*setup%n_1) + 1
    dn = modulo(site_i-2, 2*setup%n_1) + 1
    upupup = modulo(site_i+2, 2*setup%n_1) + 1
    dndndn = modulo(site_i-4, 2*setup%n_1) + 1
    ltltlt = modulo(site_j+2, 2*setup%n_2) + 1
    rtrtrt = modulo(site_j-4, 2*setup%n_2) + 1
    fwfwfw = modulo(site_k+2, 2*setup%n_3) + 1
    bwbwbw = modulo(site_k-4, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b

    allocate(nbrs(24))
    nbrs(1)   = config(ib, up, lt, fwfwfw)
    nbrs(2)   = config(ib, dn, lt, fwfwfw)
    nbrs(3)   = config(ib, up, rt, fwfwfw)
    nbrs(4)   = config(ib, dn, rt, fwfwfw)
    nbrs(5)   = config(ib, up, lt, bwbwbw)
    nbrs(6)   = config(ib, dn, lt, bwbwbw)
    nbrs(7)   = config(ib, up, rt, bwbwbw)
    nbrs(8)   = config(ib, dn, rt, bwbwbw)
    nbrs(9)   = config(ib, up, ltltlt, fw)
    nbrs(10)  = config(ib, dn, ltltlt, fw)
    nbrs(11)  = config(ib, up, ltltlt, bw)
    nbrs(12)  = config(ib, dn, ltltlt, bw)
    nbrs(13)  = config(ib, up, rtrtrt, fw)
    nbrs(14)  = config(ib, dn, rtrtrt, fw)
    nbrs(15)  = config(ib, up, rtrtrt, bw)
    nbrs(16)  = config(ib, dn, rtrtrt, bw)
    nbrs(17)  = config(ib, upupup, lt, fw)
    nbrs(18)  = config(ib, upupup, rt, fw)
    nbrs(19)  = config(ib, upupup, lt, bw)
    nbrs(20)  = config(ib, upupup, rt, bw)
    nbrs(21)  = config(ib, dndndn, lt, fw)
    nbrs(22)  = config(ib, dndndn, rt, fw)
    nbrs(23)  = config(ib, dndndn, lt, bw)
    nbrs(24)  = config(ib, dndndn, rt, bw)
     
    ! Sum them
    do i=1, 24
      energy = energy + V_ex(species, nbrs(i),4)
    end do
    deallocate(nbrs)
  end function bcc_shell4_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 5th coordination     !
  ! shell to the energy for the BCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_shell5_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, upup, dndn, fwfw, bwbw, ltlt, rtrt, ib

    energy=0.0_real64
    
    upup = modulo(site_i+1, 2*setup%n_1) + 1
    dndn = modulo(site_i-3, 2*setup%n_1) + 1
    ltlt = modulo(site_j+1, 2*setup%n_2) + 1
    rtrt = modulo(site_j-3, 2*setup%n_2) + 1
    fwfw = modulo(site_k+1, 2*setup%n_3) + 1
    bwbw = modulo(site_k-3, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b

    allocate(nbrs(8))
    nbrs(1)  = config(ib, upup, ltlt, fwfw)
    nbrs(2)  = config(ib, dndn, ltlt, fwfw)
    nbrs(3)  = config(ib, upup, rtrt, fwfw)
    nbrs(4)  = config(ib, dndn, rtrt, fwfw)
    nbrs(5)  = config(ib, upup, ltlt, bwbw)
    nbrs(6)  = config(ib, dndn, ltlt, bwbw)
    nbrs(7)  = config(ib, upup, rtrt, bwbw)
    nbrs(8)  = config(ib, dndn, rtrt, bwbw)
     
    ! Sum them
    do i=1, 8
      energy = energy + V_ex(species, nbrs(i),5)
    end do
    deallocate(nbrs)
  end function bcc_shell5_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 6th coordination     !
  ! shell to the energy for the BCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_shell6_energy(setup, site_b, site_i, site_j, site_k, &
                             config,  species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, upupupup, dndndndn, fwfwfwfw, bwbwbwbw, &
               ltltltlt, rtrtrtrt, ib

    energy=0.0_real64
    
    upupupup = modulo(site_i+3, 2*setup%n_1) + 1
    dndndndn = modulo(site_i-5, 2*setup%n_1) + 1
    ltltltlt = modulo(site_j+3, 2*setup%n_2) + 1
    rtrtrtrt = modulo(site_j-5, 2*setup%n_2) + 1
    fwfwfwfw = modulo(site_k+3, 2*setup%n_3) + 1
    bwbwbwbw = modulo(site_k-5, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b

    allocate(nbrs(6))
    nbrs(1)  = config(ib, upupupup, site_j, site_k)
    nbrs(2)  = config(ib, dndndndn, site_j, site_k)
    nbrs(3)  = config(ib, site_i, ltltltlt, site_k)
    nbrs(4)  = config(ib, site_i, rtrtrtrt, site_k)
    nbrs(5)  = config(ib, site_i, site_j, fwfwfwfw)
    nbrs(6)  = config(ib, site_i, site_j, bwbwbwbw)
     
    ! Sum them
    do i=1, 6
      energy = energy + V_ex(species, nbrs(i),6)
    end do
    deallocate(nbrs)
  end function bcc_shell6_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 7th coordination     !
  ! shell to the energy for the BCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_shell7_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, up, dn, fw, bw, lt, rt, upupup, dndndn, &
               fwfwfw, bwbwbw, ltltlt, rtrtrt, ib

    energy=0.0_real64

    up = modulo(site_i, 2*setup%n_1) + 1
    dn = modulo(site_i-2, 2*setup%n_1) + 1
    lt = modulo(site_j, 2*setup%n_2) + 1
    rt = modulo(site_j-2, 2*setup%n_2) + 1
    fw = modulo(site_k, 2*setup%n_3) + 1
    bw = modulo(site_k-2, 2*setup%n_3) + 1
    upupup = modulo(site_i+2, 2*setup%n_1) + 1
    dndndn = modulo(site_i-4, 2*setup%n_1) + 1
    ltltlt = modulo(site_j+2, 2*setup%n_2) + 1
    rtrtrt = modulo(site_j-4, 2*setup%n_2) + 1
    fwfwfw = modulo(site_k+2, 2*setup%n_3) + 1
    bwbwbw = modulo(site_k-4, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b

    allocate(nbrs(24))

    nbrs(1)   = config(ib,     up, ltltlt, fwfwfw)
    nbrs(2)   = config(ib,     dn, ltltlt, fwfwfw)
    nbrs(3)   = config(ib, upupup,     lt, fwfwfw)
    nbrs(4)   = config(ib, dndndn,     lt, fwfwfw)
    nbrs(5)   = config(ib, upupup,     rt, fwfwfw)
    nbrs(6)   = config(ib, dndndn,     rt, fwfwfw)
    nbrs(7)   = config(ib,     up, rtrtrt, fwfwfw)
    nbrs(8)   = config(ib,     dn, rtrtrt, fwfwfw)
    nbrs(9)   = config(ib, upupup, ltltlt,     fw)
    nbrs(10)  = config(ib, dndndn, ltltlt,     fw)
    nbrs(11)  = config(ib, upupup, rtrtrt,     fw)
    nbrs(12)  = config(ib, dndndn, rtrtrt,     fw)
    nbrs(13)  = config(ib, upupup, ltltlt,     bw)
    nbrs(14)  = config(ib, dndndn, ltltlt,     bw)
    nbrs(15)  = config(ib, upupup, rtrtrt,     bw)
    nbrs(16)  = config(ib, dndndn, rtrtrt,     bw)
    nbrs(17)  = config(ib,     up, ltltlt, bwbwbw)
    nbrs(18)  = config(ib,     dn, ltltlt, bwbwbw)
    nbrs(19)  = config(ib, upupup,     lt, bwbwbw)
    nbrs(20)  = config(ib, dndndn,     lt, bwbwbw)
    nbrs(21)  = config(ib, upupup,     rt, bwbwbw)
    nbrs(22)  = config(ib, dndndn,     rt, bwbwbw)
    nbrs(23)  = config(ib,     up, rtrtrt, bwbwbw)
    nbrs(24)  = config(ib,     dn, rtrtrt, bwbwbw)

    ! Sum them
    do i=1, 24
      energy = energy + V_ex(species, nbrs(i),7)
    end do

    deallocate(nbrs)

  end function bcc_shell7_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 8th coordination     !
  ! shell to the energy for the BCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_shell8_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, upup, dndn, fwfw, bwbw, ltlt, rtrt, upupupup, &
               dndndndn, fwfwfwfw, bwbwbwbw, ltltltlt, rtrtrtrt, ib

    energy=0.0_real64

    upup = modulo(site_i+1, 2*setup%n_1) + 1
    dndn = modulo(site_i-3, 2*setup%n_1) + 1
    ltlt = modulo(site_j+1, 2*setup%n_2) + 1
    rtrt = modulo(site_j-3, 2*setup%n_2) + 1
    fwfw = modulo(site_k+1, 2*setup%n_3) + 1
    bwbw = modulo(site_k-3, 2*setup%n_3) + 1
    upupupup = modulo(site_i+3, 2*setup%n_1) + 1
    dndndndn = modulo(site_i-5, 2*setup%n_1) + 1
    ltltltlt = modulo(site_j+3, 2*setup%n_2) + 1
    rtrtrtrt = modulo(site_j-5, 2*setup%n_2) + 1
    fwfwfwfw = modulo(site_k+3, 2*setup%n_3) + 1
    bwbwbwbw = modulo(site_k-5, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b

    allocate(nbrs(24))

    nbrs(1)   = config(ib,   site_i,     ltlt, fwfwfwfw)
    nbrs(2)   = config(ib,     upup,   site_j, fwfwfwfw)
    nbrs(3)   = config(ib,     dndn,   site_j, fwfwfwfw)
    nbrs(4)   = config(ib,   site_i,     rtrt, fwfwfwfw)
    nbrs(5)   = config(ib,   site_i, ltltltlt,     fwfw)
    nbrs(6)   = config(ib, dndndndn,   site_j,     fwfw)
    nbrs(7)   = config(ib, upupupup,   site_j,     fwfw)
    nbrs(8)   = config(ib,   site_i, rtrtrtrt,     fwfw)
    nbrs(9)   = config(ib,     upup, ltltltlt,   site_k)
    nbrs(10)  = config(ib,     dndn, ltltltlt,   site_k)
    nbrs(11)  = config(ib, upupupup,     ltlt,   site_k)
    nbrs(12)  = config(ib, dndndndn,     ltlt,   site_k)
    nbrs(13)  = config(ib, upupupup,     rtrt,   site_k)
    nbrs(14)  = config(ib, dndndndn,     rtrt,   site_k)
    nbrs(15)  = config(ib,   site_i, ltltltlt,   site_k)
    nbrs(16)  = config(ib, dndndndn,   site_j,   site_k)
    nbrs(17)  = config(ib,   site_i, ltltltlt,     bwbw)
    nbrs(18)  = config(ib, dndndndn,   site_j,     bwbw)
    nbrs(19)  = config(ib, upupupup,   site_j,     bwbw)
    nbrs(20)  = config(ib,   site_i, rtrtrtrt,     bwbw)
    nbrs(21)  = config(ib,   site_i,     ltlt, bwbwbwbw)
    nbrs(22)  = config(ib,     upup,   site_j, bwbwbwbw)
    nbrs(23)  = config(ib,     dndn,   site_j, bwbwbwbw)
    nbrs(24)  = config(ib,   site_i,     rtrt, bwbwbwbw)

    ! Sum them
    do i=1, 24
      energy = energy + V_ex(species, nbrs(i),8)
    end do

    deallocate(nbrs)

  end function bcc_shell8_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 9th coordination     !
  ! shell to the energy for the BCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_shell9_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, upup, dndn, fwfw, bwbw, ltlt, rtrt, upupupup, &
               dndndndn, fwfwfwfw, bwbwbwbw, ltltltlt, rtrtrtrt, ib

    energy=0.0_real64

    upup = modulo(site_i+1, 2*setup%n_1) + 1
    dndn = modulo(site_i-3, 2*setup%n_1) + 1
    ltlt = modulo(site_j+1, 2*setup%n_2) + 1
    rtrt = modulo(site_j-3, 2*setup%n_2) + 1
    fwfw = modulo(site_k+1, 2*setup%n_3) + 1
    bwbw = modulo(site_k-3, 2*setup%n_3) + 1
    upupupup = modulo(site_i+3, 2*setup%n_1) + 1
    dndndndn = modulo(site_i-5, 2*setup%n_1) + 1
    ltltltlt = modulo(site_j+3, 2*setup%n_2) + 1
    rtrtrtrt = modulo(site_j-5, 2*setup%n_2) + 1
    fwfwfwfw = modulo(site_k+3, 2*setup%n_3) + 1
    bwbwbwbw = modulo(site_k-5, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b

    allocate(nbrs(24))

    nbrs(1)   = config(ib,     upup,     ltlt, fwfwfwfw)
    nbrs(2)   = config(ib,     dndn,     ltlt, fwfwfwfw)
    nbrs(3)   = config(ib,     upup,     rtrt, fwfwfwfw)
    nbrs(4)   = config(ib,     dndn,     rtrt, fwfwfwfw)
    nbrs(5)   = config(ib,     upup, ltltltlt,     fwfw)
    nbrs(6)   = config(ib,     dndn, ltltltlt,     fwfw)
    nbrs(7)   = config(ib, upupupup,     ltlt,     fwfw)
    nbrs(8)   = config(ib, dndndndn,     ltlt,     fwfw)
    nbrs(9)   = config(ib, upupupup,     rtrt,     fwfw)
    nbrs(10)  = config(ib, dndndndn,     rtrt,     fwfw)
    nbrs(11)  = config(ib,     upup, rtrtrtrt,     fwfw)
    nbrs(12)  = config(ib,     dndn, rtrtrtrt,     fwfw)
    nbrs(13)  = config(ib,     upup, ltltltlt,     bwbw)
    nbrs(14)  = config(ib,     dndn, ltltltlt,     bwbw)
    nbrs(15)  = config(ib, upupupup,     ltlt,     bwbw)
    nbrs(16)  = config(ib, dndndndn,     ltlt,     bwbw)
    nbrs(17)  = config(ib, upupupup,     rtrt,     bwbw)
    nbrs(18)  = config(ib, dndndndn,     rtrt,     bwbw)
    nbrs(19)  = config(ib,     upup, rtrtrtrt,     bwbw)
    nbrs(20)  = config(ib,     dndn, rtrtrtrt,     bwbw)
    nbrs(21)  = config(ib,     upup,     ltlt, bwbwbwbw)
    nbrs(22)  = config(ib,     dndn,     ltlt, bwbwbwbw)
    nbrs(23)  = config(ib,     upup,     rtrt, bwbwbwbw)
    nbrs(24)  = config(ib,     dndn,     rtrt, bwbwbwbw)

    ! Sum them
    do i=1, 24
      energy = energy + V_ex(species, nbrs(i),9)
    end do

    deallocate(nbrs)

  end function bcc_shell9_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 10th coordination    !
  ! shell to the energy for the BCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_shell10_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, up, dn, fw, bw, lt, rt, upupup, dndndn, fwfwfw, &
               bwbwbw, ltltlt, rtrtrt, upupupupup, dndndndndn,    &
               fwfwfwfwfw, bwbwbwbwbw, ltltltltlt, rtrtrtrtrt, ib

    energy=0.0_real64

    up = modulo(  site_i, 2*setup%n_1) + 1
    dn = modulo(site_i-2, 2*setup%n_1) + 1
    lt = modulo(  site_j, 2*setup%n_2) + 1
    rt = modulo(site_j-2, 2*setup%n_2) + 1
    fw = modulo(site_k, 2*setup%n_3) + 1
    bw = modulo(site_k-2, 2*setup%n_3) + 1
    upupup = modulo(site_i+2, 2*setup%n_1) + 1
    dndndn = modulo(site_i-4, 2*setup%n_1) + 1
    ltltlt = modulo(site_j+2, 2*setup%n_2) + 1
    rtrtrt = modulo(site_j-4, 2*setup%n_2) + 1
    fwfwfw = modulo(site_k+2, 2*setup%n_3) + 1
    bwbwbw = modulo(site_k-4, 2*setup%n_3) + 1
    upupupupup = modulo(site_i+4, 2*setup%n_1) + 1
    dndndndndn = modulo(site_i-6, 2*setup%n_1) + 1
    ltltltltlt = modulo(site_j+4, 2*setup%n_2) + 1
    rtrtrtrtrt = modulo(site_j-6, 2*setup%n_2) + 1
    fwfwfwfwfw = modulo(site_k+4, 2*setup%n_3) + 1
    bwbwbwbwbw = modulo(site_k-6, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b

    allocate(nbrs(32))

    nbrs(1)   = config(ib,         up,         lt, fwfwfwfwfw)
    nbrs(2)   = config(ib,         dn,         lt, fwfwfwfwfw)
    nbrs(3)   = config(ib,         up,         rt, fwfwfwfwfw)
    nbrs(4)   = config(ib,         dn,         rt, fwfwfwfwfw)
    nbrs(5)   = config(ib,     upupup,     ltltlt,     fwfwfw)
    nbrs(6)   = config(ib,     dndndn,     ltltlt,     fwfwfw)
    nbrs(7)   = config(ib,     upupup,     rtrtrt,     fwfwfw)
    nbrs(8)   = config(ib,     dndndn,     rtrtrt,     fwfwfw)
    nbrs(9)   = config(ib,         up, rtrtrtrtrt,         fw)
    nbrs(10)  = config(ib,         dn, rtrtrtrtrt,         fw)
    nbrs(11)  = config(ib, upupupupup,         rt,         fw)
    nbrs(12)  = config(ib, dndndndndn,         rt,         fw)
    nbrs(13)  = config(ib, upupupupup,         lt,         fw)
    nbrs(14)  = config(ib, dndndndndn,         lt,         fw)
    nbrs(15)  = config(ib,         up, rtrtrtrtrt,         fw)
    nbrs(16)  = config(ib,         dn, rtrtrtrtrt,         fw)
    nbrs(17)  = config(ib,         up, rtrtrtrtrt,         bw)
    nbrs(18)  = config(ib,         dn, rtrtrtrtrt,         bw)
    nbrs(19)  = config(ib, upupupupup,         rt,         bw)
    nbrs(20)  = config(ib, dndndndndn,         rt,         bw)
    nbrs(21)  = config(ib, upupupupup,         lt,         bw)
    nbrs(22)  = config(ib, dndndndndn,         lt,         bw)
    nbrs(23)  = config(ib,         up, rtrtrtrtrt,         bw)
    nbrs(24)  = config(ib,         dn, rtrtrtrtrt,         bw)
    nbrs(25)  = config(ib,     upupup,     ltltlt,     bwbwbw)
    nbrs(26)  = config(ib,     dndndn,     ltltlt,     bwbwbw)
    nbrs(27)  = config(ib,     upupup,     rtrtrt,     bwbwbw)
    nbrs(28)  = config(ib,     dndndn,     rtrtrt,     bwbwbw)
    nbrs(29)  = config(ib,         up,         lt, bwbwbwbwbw)
    nbrs(30)  = config(ib,         dn,         lt, bwbwbwbwbw)
    nbrs(31)  = config(ib,         up,         rt, bwbwbwbwbw)
    nbrs(32)  = config(ib,         dn,         rt, bwbwbwbwbw)

    ! Sum them
    do i=1, 32
      energy = energy + V_ex(species, nbrs(i),10)
    end do

    deallocate(nbrs)

  end function bcc_shell10_energy

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 1st    !
  ! coordination shell on the BCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_energy_1shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy= bcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species)
    
  end function bcc_energy_1shells

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 2nd    !
  ! coordination shell on the BCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_energy_2shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = bcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell2_energy(setup, site_b, site_i, site_j, site_k, config, species)
    
  end function bcc_energy_2shells

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 3rd    !
  ! coordination shell on the BCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_energy_3shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = bcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell2_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell3_energy(setup, site_b, site_i, site_j, site_k, config, species)
    
  end function bcc_energy_3shells

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 4th    !
  ! coordination shell on the BCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_energy_4shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = bcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell2_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell3_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell4_energy(setup, site_b, site_i, site_j, site_k, config, species)
    
  end function bcc_energy_4shells

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 5th    !
  ! coordination shell on the BCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_energy_5shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = bcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell2_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell3_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell4_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell5_energy(setup, site_b, site_i, site_j, site_k, config, species)
    
  end function bcc_energy_5shells

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 6th    !
  ! coordination shell on the BCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_energy_6shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = bcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell2_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell3_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell4_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell5_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell6_energy(setup, site_b, site_i, site_j, site_k, config, species)

  end function bcc_energy_6shells

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 7th    !
  ! coordination shell on the BCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_energy_7shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = bcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell2_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell3_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell4_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell5_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell6_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell7_energy(setup, site_b, site_i, site_j, site_k, config, species)

  end function bcc_energy_7shells

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 8th    !
  ! coordination shell on the BCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_energy_8shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = bcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell2_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell3_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell4_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell5_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell6_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell7_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell8_energy(setup, site_b, site_i, site_j, site_k, config, species)

  end function bcc_energy_8shells

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 9th    !
  ! coordination shell on the BCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_energy_9shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = bcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell2_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell3_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell4_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell5_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell6_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell7_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell8_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell9_energy(setup, site_b, site_i, site_j, site_k, config, species)

  end function bcc_energy_9shells

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 10th   !
  ! coordination shell on the BCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function bcc_energy_10shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = bcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell2_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell3_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell4_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell5_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell6_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell7_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell8_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell9_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + bcc_shell10_energy(setup, site_b, site_i, site_j, site_k, config, species)

  end function bcc_energy_10shells

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 1st coordination     !
  ! shell to the energy for the FCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function fcc_shell1_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, up, dn, fw, bw, lt, rt, ib

    energy=0.0_real64
    
    ! Compute where my neighbours are
    up = modulo(site_i, 2*setup%n_1) + 1
    dn = modulo(site_i-2, 2*setup%n_1) + 1
    lt = modulo(site_j, 2*setup%n_2) + 1
    rt = modulo(site_j-2, 2*setup%n_2) + 1
    fw = modulo(site_k, 2*setup%n_3) + 1
    bw = modulo(site_k-2, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b
      
    allocate(nbrs(12))
    nbrs(1)  = config(ib, site_i, rt, fw)
    nbrs(2)  = config(ib, site_i, rt, bw)
    nbrs(3)  = config(ib, site_i, lt, fw)
    nbrs(4)  = config(ib, site_i, lt, bw)
    nbrs(5)  = config(ib, up, rt, site_k)
    nbrs(6)  = config(ib, up, lt, site_k)
    nbrs(7)  = config(ib, up, site_j, fw)
    nbrs(8)  = config(ib, up, site_j, bw)
    nbrs(9)  = config(ib, dn, rt, site_k)
    nbrs(10) = config(ib, dn, lt, site_k)
    nbrs(11) = config(ib, dn, site_j, fw)
    nbrs(12) = config(ib, dn, site_j, bw)
    do i=1, 12
      energy = energy + V_ex(species, nbrs(i), 1)
    end do
    deallocate(nbrs)
  end function fcc_shell1_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 2nd coordination     !
  ! shell to the energy for the FCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function fcc_shell2_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i
    integer :: upup, dndn, fwfw, bwbw, ltlt, rtrt, ib

    energy=0.0_real64
    
    upup = modulo(site_i+1, 2*setup%n_1) + 1
    dndn = modulo(site_i-3, 2*setup%n_1) + 1
    ltlt = modulo(site_j+1, 2*setup%n_2) + 1
    rtrt = modulo(site_j-3, 2*setup%n_2) + 1
    fwfw = modulo(site_k+1, 2*setup%n_3) + 1
    bwbw = modulo(site_k-3, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b

    allocate(nbrs(6))
    nbrs(1)  = config(ib, upup, site_j, site_k)
    nbrs(2)  = config(ib, dndn, site_j, site_k)
    nbrs(3)  = config(ib, site_i, ltlt, site_k)
    nbrs(4)  = config(ib, site_i, rtrt, site_k)
    nbrs(5)  = config(ib, site_i, site_j, fwfw)
    nbrs(6)  = config(ib, site_i, site_j, bwbw)
    do i=1, 6
      energy = energy + V_ex(species, nbrs(i), 2)
    end do
    deallocate(nbrs)
  end function fcc_shell2_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 3rd coordination     !
  ! shell to the energy for the FCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function fcc_shell3_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, up, dn, fw, bw, lt, rt
    integer :: upup, dndn, fwfw, bwbw, ltlt, rtrt, ib

    energy=0.0_real64
    
    up = modulo(site_i, 2*setup%n_1) + 1
    dn = modulo(site_i-2, 2*setup%n_1) + 1
    lt = modulo(site_j, 2*setup%n_2) + 1
    rt = modulo(site_j-2, 2*setup%n_2) + 1
    fw = modulo(site_k, 2*setup%n_3) + 1
    bw = modulo(site_k-2, 2*setup%n_3) + 1
      
    upup = modulo(site_i+1, 2*setup%n_1) + 1
    dndn = modulo(site_i-3, 2*setup%n_1) + 1
    ltlt = modulo(site_j+1, 2*setup%n_2) + 1
    rtrt = modulo(site_j-3, 2*setup%n_2) + 1
    fwfw = modulo(site_k+1, 2*setup%n_3) + 1
    bwbw = modulo(site_k-3, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b

    allocate(nbrs(24))
    nbrs(1)   = config(ib, dndn, lt, fw)
    nbrs(2)   = config(ib, dndn, lt, bw)
    nbrs(3)   = config(ib, dndn, rt, fw)
    nbrs(4)   = config(ib, dndn, rt, bw)
    nbrs(5)   = config(ib, upup, lt, fw)
    nbrs(6)   = config(ib, upup, lt, bw)
    nbrs(7)   = config(ib, upup, rt, fw)
    nbrs(8)   = config(ib, upup, rt, bw)
    nbrs(9)   = config(ib, up, ltlt, fw)
    nbrs(10)  = config(ib, dn, ltlt, fw)
    nbrs(11)  = config(ib, up, ltlt, bw)
    nbrs(12)  = config(ib, dn, ltlt, bw)
    nbrs(13)  = config(ib, up, rtrt, fw)
    nbrs(14)  = config(ib, dn, rtrt, fw)
    nbrs(15)  = config(ib, up, rtrt, bw)
    nbrs(16)  = config(ib, dn, rtrt, bw)
    nbrs(17)  = config(ib, up, lt, fwfw)
    nbrs(18)  = config(ib, dn, lt, fwfw)
    nbrs(19)  = config(ib, up, rt, fwfw)
    nbrs(20)  = config(ib, dn, rt, fwfw)
    nbrs(21)  = config(ib, up, lt, bwbw)
    nbrs(22)  = config(ib, dn, lt, bwbw)
    nbrs(23)  = config(ib, up, rt, bwbw)
    nbrs(24)  = config(ib, dn, rt, bwbw)
    do i=1, 24
      energy = energy + V_ex(species, nbrs(i), 3)
    end do
    deallocate(nbrs)
  end function fcc_shell3_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 4th coordination     !
  ! shell to the energy for the FCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function fcc_shell4_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i
    integer :: upup, dndn, fwfw, bwbw, ltlt, rtrt, ib

    energy=0.0_real64
    
    upup = modulo(site_i+1, 2*setup%n_1) + 1
    dndn = modulo(site_i-3, 2*setup%n_1) + 1
    ltlt = modulo(site_j+1, 2*setup%n_2) + 1
    rtrt = modulo(site_j-3, 2*setup%n_2) + 1
    fwfw = modulo(site_k+1, 2*setup%n_3) + 1
    bwbw = modulo(site_k-3, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b

    allocate(nbrs(12))
    nbrs(1)  = config(ib, upup, site_j, bwbw)
    nbrs(2)  = config(ib, dndn, site_j, bwbw)
    nbrs(3)  = config(ib, site_i, ltlt, bwbw)
    nbrs(4)  = config(ib, site_i, rtrt, bwbw)
    nbrs(5)  = config(ib, upup, ltlt, site_k)
    nbrs(6)  = config(ib, dndn, ltlt, site_k)
    nbrs(7)  = config(ib, upup, rtrt, site_k)
    nbrs(8)  = config(ib, dndn, rtrt, site_k)
    nbrs(9)  = config(ib, upup, site_j, fwfw)
    nbrs(10) = config(ib, dndn, site_j, fwfw)
    nbrs(11) = config(ib, site_i, ltlt, fwfw)
    nbrs(12) = config(ib, site_i, rtrt, fwfw)
    do i=1, 12
      energy = energy + V_ex(species, nbrs(i), 4)
    end do
    deallocate(nbrs)
  end function fcc_shell4_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 5th coordination     !
  ! shell to the energy for the FCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function fcc_shell5_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, up, dn, fw, bw, lt, rt
    integer :: upupup, dndndn, fwfwfw, bwbwbw, ltltlt, rtrtrt, ib

    energy=0.0_real64
    
    up = modulo(site_i, 2*setup%n_1) + 1
    dn = modulo(site_i-2, 2*setup%n_1) + 1
    lt = modulo(site_j, 2*setup%n_2) + 1
    rt = modulo(site_j-2, 2*setup%n_2) + 1
    fw = modulo(site_k, 2*setup%n_3) + 1
    bw = modulo(site_k-2, 2*setup%n_3) + 1
      
    upupup = modulo(site_i+2, 2*setup%n_1) + 1
    dndndn = modulo(site_i-4, 2*setup%n_1) + 1
    ltltlt = modulo(site_j+2, 2*setup%n_2) + 1
    rtrtrt = modulo(site_j-4, 2*setup%n_2) + 1
    fwfwfw = modulo(site_k+2, 2*setup%n_3) + 1
    bwbwbw = modulo(site_k-4, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b

    allocate(nbrs(24))
    nbrs(1)   = config(ib, up, site_j, bwbwbw)
    nbrs(2)   = config(ib, dn, site_j, bwbwbw)
    nbrs(3)   = config(ib, site_i, rt, bwbwbw)
    nbrs(4)   = config(ib, site_i, lt, bwbwbw)
    nbrs(5)   = config(ib, site_i, ltltlt, bw)
    nbrs(6)   = config(ib, upupup, site_j, bw)
    nbrs(7)   = config(ib, dndndn, site_j, bw)
    nbrs(8)   = config(ib, site_i, rtrtrt, bw)
    nbrs(9)   = config(ib, up, ltltlt, site_k)
    nbrs(10)  = config(ib, dn, ltltlt, site_k)
    nbrs(11)  = config(ib, upupup, lt, site_k)
    nbrs(12)  = config(ib, dndndn, lt, site_k)
    nbrs(13)  = config(ib, upupup, rt, site_k)
    nbrs(14)  = config(ib, dndndn, rt, site_k)
    nbrs(15)  = config(ib, up, rtrtrt, site_k)
    nbrs(16)  = config(ib, dn, rtrtrt, site_k)
    nbrs(17)  = config(ib, site_i, ltltlt, fw)
    nbrs(18)  = config(ib, upupup, site_j, fw)
    nbrs(19)  = config(ib, dndndn, site_j, fw)
    nbrs(20)  = config(ib, site_k, rtrtrt, fw)
    nbrs(21)  = config(ib, up, site_j, fwfwfw)
    nbrs(22)  = config(ib, dn, site_j, fwfwfw)
    nbrs(23)  = config(ib, site_i, rt, fwfwfw)
    nbrs(24)  = config(ib, site_i, lt, fwfwfw)
    do i=1, 24
      energy = energy + V_ex(species, nbrs(i), 5)
    end do
    deallocate(nbrs)
  end function fcc_shell5_energy

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 6th coordination     !
  ! shell to the energy for the FCC lattice                            !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function fcc_shell6_energy(setup, site_b, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, upup, dndn, fwfw, bwbw, ltlt, rtrt, ib

    energy=0.0_real64

    upup = modulo(site_i+1, 2*setup%n_1) + 1
    dndn = modulo(site_i-3, 2*setup%n_1) + 1
    ltlt = modulo(site_j+1, 2*setup%n_2) + 1
    rtrt = modulo(site_j-3, 2*setup%n_2) + 1
    fwfw = modulo(site_k+1, 2*setup%n_3) + 1
    bwbw = modulo(site_k-3, 2*setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b
    
    allocate(nbrs(8))
    nbrs(1)   = config(ib, upup, ltlt, bwbw)
    nbrs(2)   = config(ib, dndn, ltlt, bwbw)
    nbrs(3)   = config(ib, upup, rtrt, bwbw)
    nbrs(4)   = config(ib, dndn, rtrt, bwbw)
    nbrs(5)   = config(ib, upup, ltlt, fwfw)
    nbrs(6)   = config(ib, dndn, ltlt, fwfw)
    nbrs(7)   = config(ib, upup, rtrt, fwfw)
    nbrs(8)   = config(ib, dndn, rtrt, fwfw)
    do i=1, 8
      energy = energy + V_ex(species, nbrs(i), 6)
    end do
    deallocate(nbrs)
  end function fcc_shell6_energy

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 1st    !
  ! coordination shell on the FCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function fcc_energy_1shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy= fcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species)
    
  end function fcc_energy_1shells

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 2nd    !
  ! coordination shell on the FCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function fcc_energy_2shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = fcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell2_energy(setup, site_b, site_i, site_j, site_k, config, species)
    
  end function fcc_energy_2shells

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 3rd    !
  ! coordination shell on the FCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function fcc_energy_3shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = fcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell2_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell3_energy(setup, site_b, site_i, site_j, site_k, config, species)
    
  end function fcc_energy_3shells

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 4th    !
  ! coordination shell on the FCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function fcc_energy_4shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = fcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell2_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell3_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell4_energy(setup, site_b, site_i, site_j, site_k, config, species)
    
  end function fcc_energy_4shells

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 5th    !
  ! coordination shell on the FCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function fcc_energy_5shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = fcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell2_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell3_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell4_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell5_energy(setup, site_b, site_i, site_j, site_k, config, species)
    
  end function fcc_energy_5shells

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 6th    !
  ! coordination shell on the FCC lattice.                             !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function fcc_energy_6shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = fcc_shell1_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell2_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell3_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell4_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell5_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           + fcc_shell6_energy(setup, site_b, site_i, site_j, site_k, config, species)
    
  end function fcc_energy_6shells

  !--------------------------------------------------------------------!
  ! Function to compute the contribution from the 1st coordination     !
  ! shell to the energy for the simple cubic lattice                   !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function simple_cubic_1shell_energy(setup, site_b, site_i, site_j, site_k, config, species) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, up, dn, fw, bw, lt, rt, ib

    energy=0.0_real64
    
    ! Compute where my neighbours are
    up = modulo(  site_i, setup%n_1) + 1
    dn = modulo(site_i-2, setup%n_1) + 1
    lt = modulo(  site_j, setup%n_2) + 1
    rt = modulo(site_j-2, setup%n_2) + 1
    fw = modulo(  site_k, setup%n_3) + 1
    bw = modulo(site_k-2, setup%n_3) + 1

    ! Basis index (always =1 for this lattice implementation,
    ! but keep here for generality)
    ib = site_b
      
    allocate(nbrs(6))

    ! Compute the energies of neighbours
    nbrs(1) = config(ib,    up, site_j, site_k)
    nbrs(2) = config(ib,    dn, site_j, site_k)
    nbrs(3) = config(ib,site_i,     lt, site_k)
    nbrs(4) = config(ib,site_i,     rt, site_k)
    nbrs(5) = config(ib,site_i, site_j,     fw)
    nbrs(6) = config(ib,site_i, site_j,     bw)
    
    ! Sum them
    do i=1, 6
      energy = energy + V_ex(species, nbrs(i),1)
    end do
 
    deallocate(nbrs)
  end function simple_cubic_1shell_energy

  !--------------------------------------------------------------------!
  ! Function to compute the energy for an interaction up to the 1st    !
  ! coordination shell on the simple cubic lattice.                    !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  function simple_cubic_energy_1shells(setup, config, site_b, site_i, site_j, site_k) &
           result(energy)
    !integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    integer(int16), dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_b, site_i, site_j, site_k
    integer(int16) :: species

    species = config(site_b, site_i, site_j, site_k)

    energy = simple_cubic_1shell_energy(setup, site_b, site_i, site_j, site_k, config, species)
    
  end function simple_cubic_energy_1shells

end module energetics
