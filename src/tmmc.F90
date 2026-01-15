!> @file    tmmc.f90                              
!>
!> @brief   Assorted routines and tools to perform Transition Matric Monte Carlo (TMMC)
!>
!> @details WIP. DO NOT USE. This module does not function correctly.                                                            
!>
!> @author  H. J. Naguszewski
!> @date    2024 
module tmmc

#ifdef USE_MPI

  use initialise
  use kinds
  use shared_data
  use derived_types
  use netcdf_io
  use comms
  use analytics
  use c_functions
  use random_site
  use metropolis
  use mpi

  implicit none

  contains

  !------------------------------------------------------------------!
  ! Main TMMC routine.                                               !
  !                                                                  !
  ! H. J. Naguszewski,  Warwick                                 2024 !
  !------------------------------------------------------------------!
  subroutine tmmc_main(setup, tmmc_setup, my_rank)

    ! Rank of this processor
    integer, intent(in) :: my_rank
    integer :: ierror, num_proc

    ! Input file data types
    type(run_params) :: setup
    type(tmmc_params) :: tmmc_setup

    ! Loop integers and error handling variable
    integer :: i, j, ierr

    ! Temperature and temperature steps
    real(real64) :: temp, beta

    ! tmmc variables and arrays
    integer :: bins
    real(real64) :: bin_width, energy_to_ry, target_energy, acceptance
    real(real64), allocatable :: bin_edges(:), probability_dist(:), bin_probability(:), bin_probability_buffer(:), &
                                 energy_bias(:)
    real(real64), allocatable :: trans_matrix(:, :), norm_trans_matrix(:, :), trans_matrix_buffer(:, :), energy_bias_all(:, :)

    ! MPI variables
    integer :: mpi_bins, mpi_start_idx, mpi_end_idx, mpi_index
    integer, allocatable :: intervals(:,:)
    real(real64) :: start, end, reduce_time, bias_time, tmmc_time, bin_overlap
    real(real64), allocatable :: mpi_bin_edges(:), energy_bias_mpi(:)

    ! window variables
    integer, allocatable :: window_indices(:, :)
    integer :: num_windows, num_walkers

    if (my_rank == 0) then
      write (6, '(72("~"))')
      write (6, '(5("~"),x,"Error: WIP module. Does not function.",x,6("~"))')
      write (6, '(72("~"))')
    end if
    call MPI_FINALIZE(ierror)
    call EXIT(0)

    ! Get number of MPI processes
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_proc, ierror)
    tmmc_setup%mc_sweeps = tmmc_setup%mc_sweeps/tmmc_setup%weight_update

    ! Check if number of MPI processes is divisible by number of windows
    bins = tmmc_setup%bins
    num_windows = tmmc_setup%num_windows
    num_walkers = num_proc/num_windows
    if (MOD(num_proc, num_windows) /= 0) then
      if (my_rank == 0) then
        write (6, '(72("~"))')
        write (6, '(5("~"),x,"Error: Number of MPI processes not divisible by num_windows",x,6("~"))')
        write (6, '(72("~"))')
      end if
      call MPI_FINALIZE(ierror)
      call EXIT(0)
    end if

    ! Get start and end indices for energy windows
    allocate (window_indices(num_windows, 2))
    allocate(intervals(num_windows, 2))
    
    call divide_range(1, bins, num_windows, intervals)
    bin_overlap = tmmc_setup%bin_overlap

    window_indices(1, 1) = intervals(1,1)
    window_indices(1,2) = INT(intervals(1,2) + ABS(intervals(2,1)-intervals(2,2))*bin_overlap)
    do i = 2, num_windows-1
      window_indices(i, 1) = INT(intervals(i,1) - ABS(intervals(i-1,1)-intervals(i-1,2))*bin_overlap)
      window_indices(i, 2) = INT(intervals(i,2) + ABS(intervals(i+1,1)-intervals(i+1,2))*bin_overlap)
    end do
    window_indices(num_windows, 1) = INT(intervals(num_windows,1) - ABS(intervals(num_windows-1,1) &
                                    -intervals(num_windows-1,2))*bin_overlap)
    window_indices(num_windows,2) = intervals(num_windows,2)

    mpi_index = my_rank/num_walkers + 1
    mpi_start_idx = window_indices(mpi_index, 1)
    mpi_end_idx = window_indices(mpi_index, 2)
    mpi_bins = mpi_end_idx - mpi_start_idx + 1

    if (my_rank == 0) then
      print*, "Window Indices"
      print*, window_indices(:, 1)
      print*, window_indices(:, 2)
      print*, "Window Bins"
      print*, window_indices(:, 2) - window_indices(:, 1) + 1
    end if

    ! Allocate arrays
    allocate (bin_edges(bins + 1))
    allocate (probability_dist(bins))
    allocate (bin_probability(bins))
    allocate (trans_matrix(bins, bins))
    allocate (norm_trans_matrix(bins, bins))
    allocate (energy_bias(bins))

    ! MPI arrays
    allocate (mpi_bin_edges(mpi_bins + 1))
    allocate (energy_bias_mpi(mpi_bins))
    if (my_rank == 0) then
      allocate (trans_matrix_buffer(bins, bins))
      allocate (bin_probability_buffer(bins))
      bin_probability_buffer = 0.0_real64
      allocate (energy_bias_all(bins, tmmc_setup%weight_update))
    end if

    ! Set temperature in appropriate units for simulation
    temp = tmmc_setup%T*k_b_in_Ry
    beta = 1.0_real64/temp

    ! Conversion meV/atom to Rydberg
    energy_to_ry = setup%n_atoms/(Ry_to_eV*1000)

    ! Create energy bins and set mpi bins
    j = 1
    bin_width = (tmmc_setup%energy_max - tmmc_setup%energy_min)/real(tmmc_setup%bins)*energy_to_ry
    do i = 1, bins + 1
      bin_edges(i) = tmmc_setup%energy_min*energy_to_ry + (i - 1)*bin_width
    end do
    do i = mpi_start_idx, mpi_end_idx + 1
      mpi_bin_edges(j) = bin_edges(i)
      j = j + 1
    end do

    ! Target energy for burn in
    target_energy = (mpi_bin_edges(1) + mpi_bin_edges(SIZE(mpi_bin_edges)))/2

    ! Initialise
    energy_bias = 0.0_real64; probability_dist = 0.0_real64
    trans_matrix = 0.0_real64; norm_trans_matrix = 0.0_real64
    bin_probability = 0.0_real64

    ! Set up the lattice
    call initial_setup(setup, config)
    call lattice_shells(setup, shells, config)

    if (my_rank == 0) then
      write (6, '(/,72("-"),/)')
      write (6, '(24("-"),x,"Commencing Simulation!",x,24("-"),/)')
      print *, "Number of atoms", setup%n_atoms
    end if
    call comms_wait()

    !---------!
    ! Burn in !
    !---------!
    call tmmc_burn_in(setup, config, target_energy, MINVAL(mpi_bin_edges), MAXVAL(mpi_bin_edges))
    print*, "Rank: ", my_rank, "Burn-in complete"
    call comms_wait()

    if (my_rank == 0) then
      write (*, *)
      write (6, '(27("-"),x,"Burn-in complete",x,27("-"),/)')
      write (*, *)
    end if

    !--------------------!
    ! Main Loop          !
    !--------------------!
    do i = 1, tmmc_setup%weight_update
      start = MPI_Wtime()
      bin_probability = 0.0_real64
      acceptance = run_tmmc_sweeps(setup, tmmc_setup, config, temp, bins, &
                                   bin_edges, mpi_start_idx, mpi_end_idx, energy_bias, trans_matrix, bin_probability)

      call comms_wait()
      tmmc_time = MPI_Wtime()

      call MPI_REDUCE(trans_matrix, trans_matrix_buffer, bins*bins, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      call MPI_REDUCE(bin_probability, bin_probability_buffer, bins, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      reduce_time = MPI_Wtime()

      if (my_rank == 0) then
        call bias_from_tm(energy_bias, probability_dist, norm_trans_matrix, trans_matrix_buffer, &
                          bins, bin_edges, bin_width, beta, 1)
        ! Write output files
        call ncdf_writer_1d("dos_bins.dat", ierr, bin_edges)
        call ncdf_writer_1d("dos_probability.dat", ierr, probability_dist)
        call ncdf_writer_1d("bin_probability.dat", ierr, bin_probability_buffer)
        call ncdf_writer_2d("energy_bias_all.dat", ierr, energy_bias_all)
        call ncdf_writer_2d("transition_matrix.dat", ierr, trans_matrix_buffer)
      end if

      call comms_wait()
      bias_time = MPI_Wtime()
      call MPI_Bcast(energy_bias, bins, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call comms_wait()

      end = MPI_Wtime()
      if (my_rank == 0) then
        energy_bias_all(:, i) = energy_bias
        write (6, '(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a)', advance='yes') "TMMC: ", tmmc_time - start, &
          "s | R:", reduce_time - tmmc_time, "s | Bias:", bias_time - reduce_time, "s | B:", end - bias_time, "s"
        write (6, '(a,i0,a,f8.2,a,f8.2,a)', advance='yes') "Weight Update ", i, ": Accepted ", &
          (acceptance/(tmmc_setup%mc_sweeps*setup%n_atoms)*100.0), "% of Monte Carlo moves. Time taken: ", end - start, "s"
      end if
    end do

    !----------------------------!
    ! Write final output files   !
    !----------------------------!
    call MPI_REDUCE(trans_matrix, trans_matrix_buffer, bins*bins, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    call MPI_REDUCE(bin_probability*mpi_bins, bin_probability_buffer, bins, &
                    MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierror)
    if (my_rank == 0) then
      call bias_from_tm(energy_bias, probability_dist, norm_trans_matrix, trans_matrix_buffer, &
                        bins, bin_edges, bin_width, temp, 0)
      ! Normalize bins visited array
      bin_probability = bin_probability_buffer
      bin_probability = bin_probability/sum(bin_probability)

      ! Write output files
      call ncdf_writer_1d("dos_bins.dat", ierr, bin_edges)
      call ncdf_writer_1d("dos_probability.dat", ierr, probability_dist)
      call ncdf_writer_1d("bin_probability.dat", ierr, bin_probability)
      call ncdf_writer_2d("energy_bias_all.dat", ierr, energy_bias_all)
      call ncdf_writer_2d("transition_matrix.dat", ierr, trans_matrix_buffer)
    end if

    if (my_rank == 0) then
      write (*, *)
      write (6, '(25("-"),x,"Simulation Complete!",x,25("-"))')
    end if

  end subroutine tmmc_main

  !------------------------------------------------------------------!
  ! Routine for obtaining index of particular bin.                   !
  !                                                                  !
  ! H. J. Naguszewski,  Warwick                                 2024 !
  !------------------------------------------------------------------!
  integer function bin_index(energy, bin_edges, bins) result(index)
    integer, intent(in) :: bins
    real(real64), intent(in) :: energy
    real(real64), dimension(:), intent(in) :: bin_edges
    real(real64) :: bin_range

    bin_range = bin_edges(bins + 1) - bin_edges(1)
    index = int(((energy - bin_edges(1))/(bin_range))*real(bins)) &
            + 1
  end function bin_index

  !------------------------------------------------------------------!
  ! Routine for constructing the bias from the transition matrix.    !
  !                                                                  !
  ! H. J. Naguszewski,  Warwick                                 2024 !
  !------------------------------------------------------------------!
  subroutine bias_from_tm(energy_bias, probability_dist, norm_trans_matrix, &
                          trans_matrix, bins, bin_edges, bin_width, beta, bias_min)
    real(real64), dimension(:), intent(inout) :: energy_bias, probability_dist
    real(real64), dimension(:), intent(in) :: bin_edges
    real(real64), dimension(:, :), intent(inout) :: norm_trans_matrix
    real(real64), dimension(:, :), intent(in) :: trans_matrix
    real(real64), intent(in) :: bin_width, beta
    integer, intent(in) :: bins, bias_min

    integer, parameter :: lwmax = 8000
    integer :: info, lwork
    real(real64), dimension(bins) :: wr, wi
    real(real64), dimension(lwmax) :: work
    real(real64), dimension(bins, bins) :: work_tm, vl, vr
    logical, dimension(bins) :: max_mask

    integer :: i, ii
    real(real64) :: Pnorm, bin_energy, min_bias, mincount
    external :: dgeev

    max_mask = .True.

    ! zero normalized transition matrix
    norm_trans_matrix = 0.0_real64

    ! Compute as appropriately normalised collection matrix
    do i = 1, bins
      Pnorm = sum(trans_matrix(:, i))
      do ii = 1, bins
        if (Pnorm > 0.0_real64) then
          norm_trans_matrix(ii, i) = trans_matrix(ii, i)/Pnorm
        end if
      end do
    end do

    !-------------------------------------------------------------!
    ! Find the dominant eigenvector and store as probability_dist !
    !-------------------------------------------------------------!
    work_tm = norm_trans_matrix

    ! query work space
    lwork = -1
    call dgeev('V', 'V', bins, work_tm, bins, wr, wi, vl, bins, &
               vr, bins, work, lwork, info)

    lwork = min(lwmax, int(work(1)))

    ! solve eigenproblem
    call dgeev('V', 'V', bins, work_tm, bins, wr, wi, vl, bins, &
               vr, bins, work, lwork, info)

    ! check convergence
    IF (info .gt. 0) THEN
      WRITE (*, *) 'The algorithm failed to compute eigenvalues.'
      STOP
    END IF

    probability_dist = abs(vr(:, maxloc(wr, 1, mask=max_mask)))

    ! In case there are bins for which we have no data (yet)
    ! replace zeros with minimum non-zero probability
    if (bias_min .eq. 1) then
      mincount = minval(probability_dist, &
                        MASK=(probability_dist > 0.0_real64))

      do i = 1, bins
        probability_dist(i) = max(probability_dist(i), mincount)
      end do
    end if

    probability_dist = probability_dist/sum(probability_dist)

    ! Construct energy_bias function needed for uniform energy
    ! sampling
    do i = 1, bins
      bin_energy = bin_edges(i) + 0.5*bin_width
      energy_bias(i) = log(probability_dist(i))/beta
    end do
    ! Shift energy_bias so that minimum value is zero
    min_bias = minval(energy_bias)
    energy_bias = energy_bias - min_bias
  end subroutine bias_from_tm

  !------------------------------------------------------------------!
  ! Routine to run TMMC sweeps.                                      !
  !                                                                  !
  ! H. J. Naguszewski,  Warwick                                 2024 !
  !------------------------------------------------------------------!
  function run_tmmc_sweeps(setup, tmmc_setup, config, temp, bins, bin_edges, mpi_start_idx, &
                           mpi_end_idx, energy_bias, trans_matrix, bin_probability) result(acceptance)
    integer(array_int), dimension(:, :, :, :) :: config
    class(run_params), intent(in) :: setup
    class(tmmc_params), intent(in) :: tmmc_setup
    integer, intent(in) :: bins, mpi_start_idx, mpi_end_idx
    real(real64), dimension(:), intent(in) :: energy_bias, bin_edges
    real(real64), dimension(:), intent(inout) :: bin_probability
    real(real64), intent(in) :: temp
    real(real64), dimension(:, :), intent(inout) :: trans_matrix

    integer, dimension(4) :: rdm1, rdm2
    real(real64) :: e_swapped, e_unswapped, delta_e, beta, unswapped_bias, swapped_bias
    integer :: acceptance, i, ibin, jbin
    integer(array_int) :: site1, site2

    ! inverse temp
    beta = 1.0_real64/temp

    ! Establish total energy before any moves
    e_unswapped = setup%full_energy(config)
    e_swapped = e_unswapped

    acceptance = 0.0_real64

    do i = 1, tmmc_setup%mc_sweeps*setup%n_atoms
      ! Make one MC trial
      ! Generate random numbers
      rdm1 = setup%rdm_site()
      rdm2 = setup%rdm_site()

      ! Get what is on those sites
      site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
      site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

      call pair_swap(config, rdm1, rdm2)

      ! Compute energy if sites have different elements
      if (site1 /= site2) then
        e_swapped = setup%full_energy(config)
      end if

      ibin = bin_index(e_unswapped, bin_edges, bins)
      jbin = bin_index(e_swapped, bin_edges, bins)

      !print*, ibin, jbin, mpi_start_idx - 1, mpi_end_idx + 1

      ! Only compute energy change if within limits where V is defined
      if (jbin > 0 .and. jbin < bins + 1) then
        ! Probability of staying in ibin, ignoring energy_bias
        trans_matrix(ibin, ibin) = trans_matrix(ibin, ibin) &
                                   + 1.0_real64 - min(1.0_real64, exp(-beta*(e_swapped - e_unswapped)))

        ! Probability of moving to jbin, ignoring energy_bias
        trans_matrix(jbin, ibin) = trans_matrix(jbin, ibin) &
                                   + min(1.0_real64, exp(-beta*(e_swapped - e_unswapped)))

        ! Only compute move if within mpi bin energy range
        if (jbin > mpi_start_idx - 1 .and. jbin < mpi_end_idx + 1) then
          ! Add change in V into diff_energy
          unswapped_bias = energy_bias(ibin)
          swapped_bias = energy_bias(jbin)
          delta_e = (e_swapped + swapped_bias) - (e_unswapped + unswapped_bias)

          ! Accept or reject move
          if (genrand() .lt. exp(-beta*delta_e)) then
            bin_probability(jbin) = bin_probability(jbin) + 1.0_real64
            acceptance = acceptance + 1
            e_unswapped = e_swapped
          else
            bin_probability(ibin) = bin_probability(ibin) + 1.0_real64
            call pair_swap(config, rdm1, rdm2)
          end if
        else
          call pair_swap(config, rdm1, rdm2)
        end if
      else
        ! reject and reset
        call pair_swap(config, rdm1, rdm2)
      end if
    end do
  end function run_tmmc_sweeps

  subroutine tmmc_burn_in(setup, config, target_energy, min_e, max_e)
    integer(array_int), dimension(:, :, :, :) :: config
    class(run_params), intent(in) :: setup
    real(real64), intent(in) :: target_energy, min_e, max_e

    integer, dimension(4) :: rdm1, rdm2
    real(real64) :: e_swapped, e_unswapped, delta_e
    integer(array_int) :: site1, site2

    ! Establish total energy before any moves
    e_unswapped = setup%full_energy(config)

    do while (.True.)
      if (e_unswapped > min_e .and. e_unswapped < max_e) then
        exit
      end if
      ! Make one MC trial
      ! Generate random numbers
      rdm1 = setup%rdm_site()
      rdm2 = setup%rdm_site()

      ! Get what is on those sites
      site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
      site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

      ! If sites different proceed
      if (site1 /= site2) then
        call pair_swap(config, rdm1, rdm2)
        e_swapped = setup%full_energy(config)

        delta_e = e_swapped - e_unswapped

        ! Accept or reject move
        if (e_swapped > target_energy .and. delta_e < 0) then
          e_unswapped = e_swapped
        else if (e_swapped < target_energy .and. delta_e > 0) then
          e_unswapped = e_swapped
        else if (genrand() .lt. 0.0001_real64) then ! to prevent getting stuck in local minimum
          e_unswapped = e_swapped
        else
          call pair_swap(config, rdm1, rdm2)
        end if
      end if
    end do
  end subroutine tmmc_burn_in

  subroutine divide_range(start, finish, num_intervals, intervals)
    integer, intent(in) :: num_intervals
    integer, intent(in) :: start, finish
    integer, intent(out) :: intervals(num_intervals, 2)
    real(real64) :: factor, index, power, b, n, g, idx, z
    integer :: i

    intervals(1,1) = start
    intervals(num_intervals,2) = finish

    power = 2.5
    b = finish
    n = num_intervals + 1
    g = 0.5
    idx = 3
    z = 2

    !factor = (b-1.0_real64)/((n+1.0_real64)**power-1.0_real64)
    !factor = (1.0_real64 - b)/(EXP(g)-EXP(g*n))
    factor = (idx - b)/(EXP(g*z)-EXP(g*n))

    do i = 2, num_intervals
      !index = FLOOR(factor*(i**power-1)+1)
      !index = FLOOR(factor*EXP(g*i)+1-factor*EXP(g))
      index = FLOOR(factor*EXP(g*i)+idx-factor*EXP(g*z))
      intervals(i-1,2) = INT(index)
      intervals(i,1) = INT(index + 1)
    end do
  end subroutine divide_range

#endif

end module tmmc
