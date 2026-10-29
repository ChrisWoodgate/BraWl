!> @file    neighbour_tables.f90
!>
!> @brief   Runtime neighbour tables for lattices with a multi-atom
!>          basis, built from the lattice geometry at start-up
!>
!> @details The energy of the Bragg-Williams model depends only on
!>          which pairs of sites share an effective pair interaction
!>          (EPI), never on the distance between them. For a lattice
!>          such as hcp the set of neighbour offsets is fixed by the
!>          topology of the lattice, while the axial ratio c/a only
!>          decides how those offsets are ordered by distance and which
!>          of them happen to coincide. This module therefore works
!>          with *pair classes*: sets of neighbour offsets that are
!>          equivalent under the crystal symmetry. For hcp a class is
!>          identified by the in-plane separation of the two atoms and
!>          the number of layers (in units of c/2) between them, both
!>          of which are independent of c/a.
!>
!>          Classes are enumerated once at start-up, sorted by distance
!>          for the requested axial ratio, and assigned to *groups*,
!>          which are the blocks of the interaction file (and the
!>          shells of the ASRO analysis). With epi_indexing = 'shells'
!>          classes whose distances agree within shell_tolerance are
!>          merged into one group, which recovers the conventional
!>          coordination shells (and, at the ideal axial ratio, exactly
!>          the shells of the hard-coded ideal_hcp kernels). With
!>          epi_indexing = 'classes' every class is its own group, which
!>          allows symmetry-inequivalent neighbours at the same distance
!>          (for example the in-plane and out-of-plane nearest
!>          neighbours of hcp) to carry different interactions.
!>
!>          The resulting per-sublattice neighbour lists, ordered by
!>          group, drive both the energy kernel (hcp_energy_tables in
!>          bw_hamiltonian.f90) and the ASRO analysis (analytics.f90),
!>          so there is a single definition of which pairs belong to
!>          which group.
!>
!>          The class key (in-plane separation, layers apart) coincides
!>          with the symmetry orbits of pair vectors for all in-plane
!>          separations below 7a, which is the first distance at which
!>          two inequivalent in-plane vectors share a length. The
!>          builder refuses to tabulate beyond that distance.
!>
!> @author  C. D. Woodgate
!> @author  Claude Fable 5.1
!>
!> @date    2026
module neighbour_tables

  use kinds
  use constants
  use derived_types
  use comms, only : comms_finalise, my_rank

  implicit none

  private

  public :: build_hcp_tables, print_pair_tables, clean_up_tables,     &
            ideal_hcp_axial_ratio

  !> Axial ratio c/a of the ideal hcp lattice, sqrt(8/3)
  real(real64), parameter :: ideal_hcp_axial_ratio =                   &
                             sqrt(8.0_real64/3.0_real64)

  ! Fractional coordinates of the two hcp basis atoms: in-plane in
  ! thirds of a_1 and a_2 (first index), and along c in quarters of
  ! a_3, for basis atoms 1 and 2 (second index / element)
  integer, parameter, dimension(2,2) :: hcp_basis_thirds =             &
                                        reshape((/ 1, 2,              &
                                                   2, 1 /), (/2, 2/))
  integer, parameter, dimension(2) :: hcp_basis_quarters = (/ 1, 3 /)

  ! Largest in-plane separation (units of a) up to which the class key
  ! is guaranteed to coincide with the symmetry orbits of pair vectors
  real(real64), parameter :: max_safe_in_plane = 7.0_real64

  ! Largest enumeration box (cells) the builder will try
  integer, parameter :: max_reach = 15

  contains

  !> @brief   Subroutine to build the pair classes, groups and
  !>          neighbour lists of the hcp lattice for a given axial ratio
  !>
  !> @details Enumerates every pair vector from each basis atom within
  !>          a box of cells that is enlarged until at least
  !>          n_groups_needed groups lie safely inside it, identifies
  !>          the pair classes by their exact integer key, sorts them by
  !>          distance (ties broken by layers apart, then in-plane
  !>          separation, so that in-plane classes precede out-of-plane
  !>          ones at equal distance), merges them into groups according
  !>          to setup%epi_indexing and setup%shell_tolerance, and
  !>          fills setup%tables. The axial ratio is taken from
  !>          setup%c_over_a and the number of neighbours used by the
  !>          energy kernel and by the ASRO analysis are set from
  !>          setup%interaction_range and setup%wc_range.
  !>
  !> @author  C. D. Woodgate
  !> @author  Claude Fable 5.1
  !> @date    2026
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  n_groups_needed Number of groups to tabulate
  !>
  !> @return None
  subroutine build_hcp_tables(setup, n_groups_needed)

    type(run_params), intent(inout) :: setup
    integer, intent(in) :: n_groups_needed

    ! Half-widths of the enumeration box (cells) in-plane and along c
    integer :: reach, reach_z

    ! Enumerated pair vectors: (central basis atom, neighbour basis
    ! atom, offset along a_1, a_2, a_3, class) per entry
    integer, allocatable, dimension(:,:) :: ent
    integer :: n_ent, n_max_ent

    ! Classes: exact key (9 rho^2, layers), squared distance, rank
    integer, allocatable, dimension(:) :: cls_k1, cls_z, rank, grp
    real(real64), allocatable, dimension(:) :: cls_d2
    integer :: n_cls, n_grp, n_used

    ! Per basis atom: number of neighbours in the used classes
    integer, dimension(2) :: n_list
    integer, allocatable, dimension(:) :: key

    real(real64) :: gamma, tol, d2, d2_max
    integer :: b0, b1, di, dj, dk, x, y, z4, k1, zeta, c, i, j, n, m
    logical :: by_shells

    gamma = setup%c_over_a
    tol = setup%shell_tolerance

    if (gamma .le. 0.0_real64) then
      call comms_finalise()
      stop 'build_hcp_tables: c_over_a must be positive'
    end if

    if (trim(setup%epi_indexing) .eq. 'shells') then
      by_shells = .true.
    else if (trim(setup%epi_indexing) .eq. 'classes') then
      by_shells = .false.
    else
      call comms_finalise()
      stop "epi_indexing must be 'shells' or 'classes'"
    end if

    if (n_groups_needed .lt. 1) then
      call comms_finalise()
      stop 'build_hcp_tables: need at least one group'
    end if

    call clean_up_tables(setup)

    reach = 3
    reach_z = 2

    !------------------------------------------------------------------!
    ! Enumerate pair vectors in a box, enlarging it until enough groups !
    ! lie safely inside                                                !
    !------------------------------------------------------------------!
    do

      ! Every pair closer than sqrt(d2_max) lies inside the box: the
      ! in-plane metric satisfies rho^2 >= 3 x^2 / 4 for either lattice
      ! coordinate x, and the basis offsets add at most one cell, so
      ! in-plane separations up to (sqrt(3)/2)(reach-1) need
      ! |di|,|dj| <= reach; along c, (reach_z-1) c suffices.
      d2_max = min(0.75_real64*real(reach-1, real64)**2,              &
                   (real(reach_z-1, real64)*gamma)**2)

      n_max_ent = 4*(2*reach+1)**2*(2*reach_z+1)
      allocate(ent(6, n_max_ent))
      allocate(cls_k1(n_max_ent), cls_z(n_max_ent), cls_d2(n_max_ent))
      n_ent = 0
      n_cls = 0

      do b0=1, 2
        do b1=1, 2
          do dk=-reach_z, reach_z
            do dj=-reach, reach
              do di=-reach, reach

                ! Pair vector in thirds (in-plane) and quarters (c)
                x  = 3*di + hcp_basis_thirds(1,b1)                     &
                          - hcp_basis_thirds(1,b0)
                y  = 3*dj + hcp_basis_thirds(2,b1)                     &
                          - hcp_basis_thirds(2,b0)
                z4 = 4*dk + hcp_basis_quarters(b1)                     &
                          - hcp_basis_quarters(b0)

                ! Exact class key: 9 rho^2 and layers apart (c/2)
                k1 = x*x + y*y - x*y
                zeta = abs(z4)/2

                ! Skip the site itself
                if ((k1 .eq. 0) .and. (zeta .eq. 0)) cycle

                d2 = real(k1, real64)/9.0_real64                       &
                   + real(zeta*zeta, real64)*gamma*gamma/4.0_real64
                if (d2 .gt. d2_max) cycle

                ! Find or create the class of this pair vector
                c = 0
                do i=1, n_cls
                  if ((cls_k1(i) .eq. k1) .and. (cls_z(i) .eq. zeta)) then
                    c = i
                    exit
                  end if
                end do
                if (c .eq. 0) then
                  n_cls = n_cls + 1
                  c = n_cls
                  cls_k1(c) = k1
                  cls_z(c) = zeta
                  cls_d2(c) = d2
                end if

                n_ent = n_ent + 1
                ent(:, n_ent) = (/ b0, b1, di, dj, dk, c /)

              end do
            end do
          end do
        end do
      end do

      ! Rank the classes by distance, with ties broken by layers apart
      ! and then by in-plane separation
      allocate(rank(n_cls), grp(n_cls))
      call rank_classes(n_cls, cls_d2, cls_z, cls_k1, tol, rank)

      ! Assign groups in rank order
      n_grp = 0
      do i=1, n_cls
        c = rank(i)
        if (i .eq. 1) then
          n_grp = 1
        else if (.not. by_shells) then
          n_grp = n_grp + 1
        else if (sqrt(cls_d2(c)) .gt.                                  &
                 sqrt(cls_d2(rank(i-1)))*(1.0_real64 + tol)) then
          ! Compare with the previous class rather than the first of
          ! the group so that the merge is the same for all members
          n_grp = n_grp + 1
        end if
        grp(c) = n_grp
      end do

      ! Enough groups, with the last needed one safely inside the box?
      if (n_grp .ge. n_groups_needed) then
        do i=1, n_cls
          if (grp(rank(i)) .eq. n_groups_needed) then
            if (sqrt(cls_d2(rank(i)))*(1.0_real64 + 2.0_real64*tol)   &
                .lt. sqrt(d2_max)) exit
          end if
        end do
        if (i .le. n_cls) exit
      end if

      deallocate(ent, cls_k1, cls_z, cls_d2, rank, grp)
      reach = reach + 2
      reach_z = reach_z + 2

      if (reach .gt. max_reach) then
        call comms_finalise()
        stop 'build_hcp_tables: too many groups requested'
      end if

    end do

    !------------------------------------------------------------------!
    ! Fill the tables with the classes in the needed groups            !
    !------------------------------------------------------------------!

    ! Classes used, in rank order
    n_used = 0
    do i=1, n_cls
      if (grp(rank(i)) .le. n_groups_needed) n_used = i
    end do

    ! The class key is only a symmetry-orbit key below a limiting
    ! in-plane separation
    do i=1, n_used
      if (sqrt(real(cls_k1(rank(i)), real64))/3.0_real64 .ge.          &
          max_safe_in_plane) then
        call comms_finalise()
        stop 'build_hcp_tables: requested range beyond tabulated limit'
      end if
    end do

    setup%tables%n_classes = n_used
    setup%tables%n_groups = n_groups_needed
    allocate(setup%tables%class_key(2, n_used))
    allocate(setup%tables%class_distance(n_used))
    allocate(setup%tables%class_group(n_used))
    allocate(setup%tables%class_count(n_used, 2))
    allocate(setup%tables%group_distance(n_groups_needed))
    allocate(setup%tables%group_count(n_groups_needed, 2))
    setup%tables%class_count = 0
    setup%tables%group_count = 0
    setup%tables%group_distance = -1.0_real64

    do i=1, n_used
      c = rank(i)
      setup%tables%class_key(:, i) = (/ cls_k1(c), cls_z(c) /)
      setup%tables%class_distance(i) = sqrt(cls_d2(c))
      setup%tables%class_group(i) = grp(c)
      if (setup%tables%group_distance(grp(c)) .lt. 0.0_real64) then
        setup%tables%group_distance(grp(c)) = sqrt(cls_d2(c))
      end if
    end do

    ! Count the neighbours of each basis atom in each used class, and
    ! sort the entries by rank so that the lists are ordered by group
    allocate(key(n_ent))
    do n=1, n_ent
      key(n) = 0
      do i=1, n_used
        if (rank(i) .eq. ent(6, n)) key(n) = i
      end do
    end do
    n_list = 0
    do n=1, n_ent
      if (key(n) .eq. 0) cycle
      b0 = ent(1, n)
      n_list(b0) = n_list(b0) + 1
      setup%tables%class_count(key(n), b0) =                           &
        setup%tables%class_count(key(n), b0) + 1
    end do
    do i=1, n_used
      do b0=1, 2
        setup%tables%group_count(setup%tables%class_group(i), b0) =    &
          setup%tables%group_count(setup%tables%class_group(i), b0)    &
          + setup%tables%class_count(i, b0)
      end do
    end do

    m = maxval(n_list)
    allocate(setup%tables%nbrs(4, m, 2))
    allocate(setup%tables%nbr_group(m, 2))
    allocate(setup%tables%n_nbrs_energy(2))
    allocate(setup%tables%n_nbrs_asro(2))
    setup%tables%nbrs = 0
    setup%tables%nbr_group = 0
    setup%tables%max_offset = 0

    ! Fill the lists in order of class rank, and within a class in
    ! order of enumeration, so that the tables are reproducible
    do b0=1, 2
      j = 0
      do i=1, n_used
        do n=1, n_ent
          if ((key(n) .ne. i) .or. (ent(1, n) .ne. b0)) cycle
          j = j + 1
          setup%tables%nbrs(:, j, b0) = ent(2:5, n)
          setup%tables%nbr_group(j, b0) = setup%tables%class_group(i)
          setup%tables%max_offset = max(setup%tables%max_offset,       &
                                        maxval(abs(ent(3:5, n))))
        end do
      end do
      setup%tables%n_nbrs_energy(b0) = 0
      setup%tables%n_nbrs_asro(b0) = 0
      do n=1, j
        if (setup%tables%nbr_group(n, b0) .le. setup%interaction_range) &
          setup%tables%n_nbrs_energy(b0) = n
        if (setup%tables%nbr_group(n, b0) .le. setup%wc_range - 1)     &
          setup%tables%n_nbrs_asro(b0) = n
      end do
    end do

    setup%tables%built = .true.

    deallocate(ent, cls_k1, cls_z, cls_d2, rank, grp, key)

  end subroutine build_hcp_tables

  !> @brief   Subroutine to rank pair classes by distance
  !>
  !> @details Classes whose distances agree within the relative
  !>          tolerance are ordered by layers apart and then by
  !>          in-plane separation, so that the order is reproducible
  !>          at degenerate distances such as those of the ideal axial
  !>          ratio. A simple insertion sort suffices for the few dozen
  !>          classes involved.
  !>
  !> @author  C. D. Woodgate
  !> @author  Claude Fable 5.1
  !> @date    2026
  !>
  !> @param  n Number of classes
  !> @param  d2 Squared distance of each class
  !> @param  z Layers apart of each class
  !> @param  k1 Nine times the squared in-plane separation of each class
  !> @param  tol Relative tolerance on distance
  !> @param  rank Class indices in increasing order (output)
  !>
  !> @return None
  subroutine rank_classes(n, d2, z, k1, tol, rank)

    integer, intent(in) :: n
    real(real64), dimension(:), intent(in) :: d2
    integer, dimension(:), intent(in) :: z, k1
    real(real64), intent(in) :: tol
    integer, dimension(:), intent(out) :: rank
    integer :: i, j, c

    do i=1, n
      c = i
      j = i - 1
      do while (j .ge. 1)
        if (.not. class_before(c, rank(j))) exit
        rank(j+1) = rank(j)
        j = j - 1
      end do
      rank(j+1) = c
    end do

    contains

    !> True if class a should precede class b
    logical function class_before(a, b)
      integer, intent(in) :: a, b
      if (abs(sqrt(d2(a)) - sqrt(d2(b))) .le.                          &
          tol*max(sqrt(d2(a)), sqrt(d2(b)))) then
        class_before = (z(a) .lt. z(b)) .or.                           &
                       ((z(a) .eq. z(b)) .and. (k1(a) .lt. k1(b)))
      else
        class_before = d2(a) .lt. d2(b)
      end if
    end function class_before

  end subroutine rank_classes

  !> @brief   Subroutine to print the pair classes and groups to the
  !>          screen
  !>
  !> @details Prints, for each group, its distance, multiplicity, and
  !>          the classes it contains, followed by the ranges used for
  !>          the interactions and for the ASRO analysis, so that the
  !>          user can match the blocks of the interaction file to the
  !>          groups. Groups closer than a thousand times the merging
  !>          tolerance are flagged, as they are easy to confuse.
  !>
  !> @author  C. D. Woodgate
  !> @author  Claude Fable 5.1
  !> @date    2026
  !>
  !> @param  setup Derived type containing simulation parameters
  !>
  !> @return None
  subroutine print_pair_tables(setup)

    type(run_params), intent(in) :: setup
    integer :: g, i, n_close
    real(real64) :: gap

    if (my_rank .ne. 0) return
    if (.not. setup%tables%built) return

    write(6,'(16("-"),x,"Pair classes and groups for this lattice",x,  &
             &15("-"),/)')
    write(6,'(x,"Lattice: ",A,"   c/a = ",F8.5,"  (ideal ",F8.5,")")') &
          trim(setup%lattice), setup%c_over_a, ideal_hcp_axial_ratio
    if (trim(setup%epi_indexing) .eq. 'shells') then
      write(6,'(x,"EPI blocks indexed by distance shells; classes ",   &
               &"merged within",/,x,"a relative tolerance of ",ES9.2,  &
               &" on the distance",/)') setup%shell_tolerance
    else
      write(6,'(x,"EPI blocks indexed by pair classes (no merging ",   &
               &"of equal distances)",/)')
    end if

    write(6,'(x,"Group    d/a   Mult.   Classes (in-plane d/a, ",      &
             &"layers of c/2, mult.)")')
    write(6,'(x,66("-"))')
    do g=1, setup%tables%n_groups
      write(6,'(x,I5,x,F7.4,x,I5,3x)', advance='no')                   &
            g, setup%tables%group_distance(g),                         &
            setup%tables%group_count(g, 1)
      do i=1, setup%tables%n_classes
        if (setup%tables%class_group(i) .ne. g) cycle
        write(6,'("(",F6.4,", ",I1,", ",I2,") ")', advance='no')       &
              sqrt(real(setup%tables%class_key(1,i), real64))          &
              /3.0_real64,                                             &
              setup%tables%class_key(2,i),                             &
              setup%tables%class_count(i, 1)
      end do
      write(6,'(a)') ''
    end do
    write(6,'(x,66("-"))')

    write(6,'(x,"Interactions use groups 1 to ",I3,                    &
             &" (distances up to ",F7.4," a)")')                       &
          setup%interaction_range,                                     &
          setup%tables%group_distance(setup%interaction_range)
    if (setup%wc_range .gt. 1) then
      write(6,'(x,"ASRO analysis uses groups 1 to ",I3,                &
               &" after the on-site term")') setup%wc_range - 1
    end if

    ! Flag groups that are easy to confuse
    n_close = 0
    do g=2, setup%tables%n_groups
      gap = (setup%tables%group_distance(g)                            &
             - setup%tables%group_distance(g-1))                       &
            /setup%tables%group_distance(g)
      if (gap .lt. 1000.0_real64*setup%shell_tolerance) then
        if (n_close .eq. 0) write(6,'(a)') ''
        n_close = n_close + 1
        write(6,'(x,"Note: groups ",I3," and ",I3," differ in distance", &
                 &" by only ",ES9.2," (relative)")') g-1, g, gap
      end if
    end do

    write(6,'(/,13("-"),x,"End of pair classes and groups for this",   &
             &x,"lattice",x,12("-"),/)')

  end subroutine print_pair_tables

  !> @brief   Subroutine to free the memory held by the pair tables
  !>
  !> @author  C. D. Woodgate
  !> @author  Claude Fable 5.1
  !> @date    2026
  !>
  !> @param  setup Derived type containing simulation parameters
  !>
  !> @return None
  subroutine clean_up_tables(setup)

    type(run_params), intent(inout) :: setup

    if (allocated(setup%tables%class_key))                             &
      deallocate(setup%tables%class_key)
    if (allocated(setup%tables%class_distance))                        &
      deallocate(setup%tables%class_distance)
    if (allocated(setup%tables%class_group))                           &
      deallocate(setup%tables%class_group)
    if (allocated(setup%tables%class_count))                           &
      deallocate(setup%tables%class_count)
    if (allocated(setup%tables%group_distance))                        &
      deallocate(setup%tables%group_distance)
    if (allocated(setup%tables%group_count))                           &
      deallocate(setup%tables%group_count)
    if (allocated(setup%tables%nbrs))                                  &
      deallocate(setup%tables%nbrs)
    if (allocated(setup%tables%nbr_group))                             &
      deallocate(setup%tables%nbr_group)
    if (allocated(setup%tables%n_nbrs_energy))                         &
      deallocate(setup%tables%n_nbrs_energy)
    if (allocated(setup%tables%n_nbrs_asro))                           &
      deallocate(setup%tables%n_nbrs_asro)
    setup%tables%built = .false.
    setup%tables%n_classes = 0
    setup%tables%n_groups = 0
    setup%tables%max_offset = 0

  end subroutine clean_up_tables

end module neighbour_tables
