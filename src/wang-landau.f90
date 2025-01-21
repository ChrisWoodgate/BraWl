!----------------------------------------------------------------------!
! Wang-Landau module                                                   !
!                                                                      !
! H. Naguszewski, Warwick                                         2024 !
!----------------------------------------------------------------------!

module wang_landau
  use initialise
  use kinds
  use shared_data
  use c_functions
  use random_site
  use metropolis
  use mpi

  implicit none

contains

  subroutine wl_main(setup, wl_setup, my_rank)
    ! Rank of this processor
    integer, intent(in) :: my_rank
    integer :: ierror, mpi_processes

    ! Input file data types
    type(run_params) :: setup
    type(wl_params) :: wl_setup

    ! Loop integers and error handling variable
    integer :: i, j, ierr

    ! Temperature and temperature steps
    real(real64) :: temp, beta

    ! wl variables and arrays
    integer :: bins, resets
    real(real64) :: bin_width, energy_to_ry, wl_f, wl_f_prev, tolerance, flatness_tolerance
    real(real64) :: target_energy, min_val, flatness, bins_buffer, bins_min, radial_min, radial_min_buffer
    real(real64), allocatable :: bin_edges(:), wl_hist(:), wl_logdos(:), bin_energy(:)
    logical :: pre_sampled, rho_saved

    ! Name of file for writing radial densities at the end
    character(len=37) :: radial_file

    ! window variables
    integer, allocatable :: window_indices(:, :)
    integer :: num_windows, num_walkers

    ! mpi variables
    real(real64) :: start, end, time_max, time_min
    integer :: mpi_bins, mpi_start_idx, mpi_end_idx, mpi_index, mpi_start_idx_buffer, mpi_end_idx_buffer, beta_index
    integer, allocatable :: window_intervals(:,:), window_rank_index(:,:)
    real(real64) :: scale_factor, scale_count, wl_logdos_min, bin_overlap, beta_diff, beta_original, beta_merge
    real(real64), allocatable :: mpi_bin_edges(:), mpi_wl_hist(:), wl_logdos_buffer(:), wl_logdos_write(:)
    real(real64), allocatable :: rank_time(:), rank_time_buffer(:,:)

    ! radial density across energy
    real(real64), allocatable :: rho_of_E(:,:,:,:), rho_of_E_buffer(:,:,:,:)
    integer, allocatable :: radial_record(:), radial_record_buffer(:)
    integer :: radial_mc_steps

    ! Minimum value in array to be considered
    min_val = wl_setup%tolerance*1e-1_real64

    ! Radial densities as a function of energy
    allocate(rho_of_E(setup%n_species, setup%n_species, setup%wc_range, wl_setup%bins))
    if (my_rank == 0) then
      allocate(rho_of_E_buffer(setup%n_species, setup%n_species, setup%wc_range, wl_setup%bins))
    end if
    rho_of_E = 0.0_real64

    ! Path to radial file and name of radial file
    radial_file = "radial_densities/rho_of_E.dat"

    ! Get number of MPI processes
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_processes, ierror)

    ! Check if number of MPI processes is divisible by number of windows
    bins = wl_setup%bins
    num_windows = wl_setup%num_windows
    num_walkers = mpi_processes/num_windows
    wl_setup%radial_samples = INT(wl_setup%radial_samples/num_walkers)
    wl_setup%radial_samples = MAX(INT(wl_setup%radial_samples*num_walkers), 1)
    if (MOD(mpi_processes, num_windows) /= 0) then
      if (my_rank == 0) then
        write (6, '(72("~"))')
        write (6, '(5("~"),x,"Error: Number of MPI processes not divisible by num_windows",x,6("~"))')
        write (6, '(72("~"))')
      end if
      call MPI_FINALIZE(ierror)
      call EXIT(0)
    end if

    ! Get start and end indices for energy windows
    allocate(window_indices(num_windows, 2))
    allocate(window_intervals(wl_setup%num_windows, 2))
    allocate(window_rank_index(wl_setup%num_windows, 2))
    
    call divide_range(my_rank, wl_setup, mpi_processes, window_intervals, window_rank_index, num_walkers)

    call create_window_intervals(wl_setup, num_walkers, window_intervals, window_indices, &
    mpi_index, mpi_start_idx, mpi_end_idx, mpi_bins)

    ! allocate arrays
    allocate(radial_record(wl_setup%bins))
    allocate(radial_record_buffer(wl_setup%bins))
    allocate(bin_edges(bins + 1))
    allocate(wl_hist(bins))
    allocate(wl_logdos(bins))
    allocate(wl_logdos_buffer(bins))
    allocate(bin_energy(bins))

    ! MPI arrays
    allocate(mpi_bin_edges(mpi_bins + 1))
    allocate(mpi_wl_hist(mpi_bins))
    allocate(rank_time(num_windows))
    allocate(rank_time_buffer(num_windows,3))

    if (my_rank == 0) then
      allocate(wl_logdos_write(bins))
    end if

    call create_energy_bins(setup, wl_setup, mpi_start_idx, mpi_end_idx, bin_edges, bin_energy, mpi_bin_edges)

    ! Target energy for burn in
    target_energy = (mpi_bin_edges(1) + mpi_bin_edges(SIZE(mpi_bin_edges)))/2

    ! Initialize
    wl_hist = 0.0_real64; wl_logdos = 0.0_real64; wl_f = wl_setup%wl_f; wl_f_prev = wl_f
    flatness = 0.0_real64; pre_sampled = .False.
    mpi_wl_hist = 0.0_real64; wl_logdos_buffer = 0.0_real64; rank_time = 0.0_real64
    radial_record = 0
    radial_record_buffer = 0
    radial_mc_steps = 0
    rho_saved = .False.

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
    call burn_in(setup, config, MINVAL(mpi_bin_edges), MAXVAL(mpi_bin_edges), &
                mpi_processes, num_walkers, my_rank, mpi_index, window_rank_index, .True.)
    print*, "Rank: ", my_rank, "Burn-in complete"
    call comms_wait()
    if (my_rank == 0) then
      write (*, *)
      write (6, '(27("-"),x,"Burn-in complete",x,27("-"),/)')
      write (*, *)
    end if

    !--------------------!
    ! Main Wang-Landau   !
    !--------------------!
    do while (wl_f > wl_setup%tolerance)
      if (pre_sampled .neqv. .True.) then ! First reset after system had time to explore
        if (minval(mpi_wl_hist) > 10.0_real64) then
          start = mpi_wtime()
          pre_sampled = .True.
          mpi_wl_hist = 0.0_real64
        end if
      end if

      call sweeps(setup, wl_setup, config, num_walkers, bin_edges, mpi_start_idx, mpi_end_idx, &
                  mpi_wl_hist, wl_logdos, wl_f, mpi_index, window_intervals, radial_record, rho_of_E, rho_saved, radial_mc_steps)

      flatness = minval(mpi_wl_hist)/(sum(mpi_wl_hist)/mpi_bins)
      bins_min = count(mpi_wl_hist > min_val)/REAL(mpi_bins)

      if (pre_sampled .neqv. .True.) then
        write (6, '(a,i0,a,f6.2,a)') "Rank: ", my_rank, " | bins visited: ", bins_min*100_real64, "%"
      end if

      if (pre_sampled .eqv. .True.) then
        !Check if we're wl_setup%flatness% flat
        if (flatness > wl_setup%flatness .and. minval(mpi_wl_hist) > 10) then
          ! End timer
          end = mpi_wtime()
          call comms_wait()

          !Reset the histogram
          mpi_wl_hist = 0.0_real64
          !Reduce f
          wl_f = wl_f*1/2
          wl_logdos = wl_logdos - minval(wl_logdos, MASK=(wl_logdos > 0.0_real64))!wl_logdos_min
          wl_logdos = ABS(wl_logdos * merge(0, 1, wl_logdos < 0.0_real64))
          !Average DoS across walkers
          call dos_average(wl_setup, num_walkers, mpi_index, wl_logdos, wl_logdos_buffer)

          rank_time = 0.0_real64
          rank_time(mpi_index) = end - start
          call MPI_ALLREDUCE(end - start, time_max, 1, MPI_DOUBLE_PRECISION, &
          MPI_MAX, MPI_COMM_WORLD, ierror)

          call MPI_REDUCE(end - start, time_min, 1, MPI_DOUBLE_PRECISION, &
          MPI_MIN, 0, MPI_COMM_WORLD, ierror)

          call MPI_REDUCE(rank_time/num_walkers, rank_time_buffer(:,1), wl_setup%num_windows, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, 0, MPI_COMM_WORLD, ierror)

          call MPI_REDUCE(rank_time, rank_time_buffer(:,3), wl_setup%num_windows, MPI_DOUBLE_PRECISION, &
          MPI_MAX, 0, MPI_COMM_WORLD, ierror)

          rank_time = time_max
          rank_time(mpi_index) = end - start

          call MPI_REDUCE(rank_time, rank_time_buffer(:,2), wl_setup%num_windows, MPI_DOUBLE_PRECISION, &
          MPI_MIN, 0, MPI_COMM_WORLD, ierror)

          call MPI_BCAST(rank_time_buffer, wl_setup%num_windows*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

          call MPI_REDUCE(radial_record, radial_record_buffer, wl_setup%bins, MPI_INT, &
          MPI_SUM, 0, MPI_COMM_WORLD, ierror)

          radial_min = 0
          do i=1,wl_setup%bins
            radial_min = radial_min + REAL(MIN(radial_record_buffer(i),wl_setup%radial_samples))
          end do
          radial_min = radial_min/REAL(wl_setup%radial_samples*wl_setup%bins)
          call MPI_BCAST(radial_min, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

          call comms_wait()
          if (my_rank == 0) then
            wl_logdos_write = wl_logdos
            write (6, '(a,f20.18,a,f8.2,a)', advance='no') "Flatness reached for W-L F of: ", wl_f_prev, &
                    " | Radial samples: ", radial_min*100_real64, "%"
              write (*, *)
            do i=1, wl_setup%num_windows
              write (6, '(a,i0,a,f8.2,a,f8.2,a,f8.2,a)') "MPI Window: ", i, " | Avg. time: ", rank_time_buffer(i,1), &
              "s | Time min: ", rank_time_buffer(i,2), "s Time max: " , rank_time_buffer(i,3), "s"
            end do
            write (*, *)
            wl_f_prev = wl_f
          end if

          ! MPI send and recieve calls for combining window DoS
          call dos_combine(wl_setup, mpi_index, my_rank, num_walkers, &
          window_indices, window_rank_index, &
          wl_logdos, wl_logdos_buffer, wl_logdos_write)

          !rho_saved
          if (rho_saved .eqv. .False.) then
            if (radial_min >= 1) then
              rho_saved = .True.
              call MPI_REDUCE(rho_of_E, rho_of_E_buffer, setup%n_species*setup%n_species*setup%wc_range*wl_setup%bins, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
              if (my_rank == 0) then 
                print*, "Radial densities saved"
                do i=1, wl_setup%bins
                  rho_of_E_buffer(:,:,:,i) = rho_of_E_buffer(:,:,:,i)/radial_record_buffer(i)
                end do
                call ncdf_radial_density_writer_across_energy(radial_file, rho_of_E_buffer, shells, bin_energy, setup)
              end if
            end if
          end if

          if (my_rank == 0) then
            ! Write output files
            call ncdf_writer_1d("wl_dos_bins.dat", ierr, bin_edges)
            call ncdf_writer_1d("wl_dos.dat", ierr, wl_logdos_write)
            call ncdf_writer_1d("wl_hist.dat", ierr, wl_hist)
            wl_logdos = wl_logdos_write
          end if

          call mpi_window_optimise(wl_setup, my_rank, bin_edges, window_intervals, window_indices, &
                                  mpi_bins, mpi_index, mpi_bin_edges, mpi_wl_hist, rank_time_buffer, num_walkers)
          mpi_start_idx = window_indices(mpi_index, 1)
          mpi_end_idx = window_indices(mpi_index, 2)
          mpi_bins = mpi_end_idx - mpi_start_idx + 1
          call burn_in(setup, config, MINVAL(mpi_bin_edges), MAXVAL(mpi_bin_edges), &
                        mpi_processes, num_walkers, my_rank, mpi_index, window_rank_index, .False.)

          call MPI_BCAST(wl_logdos, wl_setup%bins, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
          ! Zero elements not worked on
          wl_logdos(1:window_indices(mpi_index,1)-1) = 0.0_real64
          wl_logdos(window_indices(mpi_index,2)+1:wl_setup%bins) = 0.0_real64
          ! Subtract minimum value
          wl_logdos = wl_logdos - minval(wl_logdos, MASK=(wl_logdos > 0.0_real64))
          wl_logdos = ABS(wl_logdos * merge(0, 1, wl_logdos < 0.0_real64))
          call comms_wait()
          start = mpi_wtime()
        end if
      end if
    end do

    if (my_rank == 0) then
      write (*, *)
      write (6, '(25("-"),x,"Simulation Complete!",x,25("-"))')
    end if

  end subroutine wl_main

  integer function bin_index(energy, bin_edges, bins) result(index)
    integer, intent(in) :: bins
    real(real64), intent(in) :: energy
    real(real64), dimension(:), intent(in) :: bin_edges
    real(real64) :: bin_range

    bin_range = bin_edges(bins + 1) - bin_edges(1)
    index = int(((energy - bin_edges(1))/(bin_range))*real(bins)) + 1
  end function bin_index

  subroutine sweeps(setup, wl_setup, config, num_walkers, bin_edges, &
                         mpi_start_idx, mpi_end_idx, mpi_wl_hist, wl_logdos, wl_f, &
                         mpi_index, window_intervals, radial_record, rho_of_E, rho_saved, radial_mc_steps)
    integer(int16), dimension(:, :, :, :) :: config
    class(run_params), intent(in) :: setup
    class(wl_params), intent(in) :: wl_setup
    integer, intent(in) :: num_walkers, mpi_start_idx, mpi_end_idx, mpi_index
    integer, dimension(:,:), intent(inout) :: window_intervals
    integer, dimension(:), intent(inout) :: radial_record
    real(real64), dimension(:,:,:,:), intent(inout) :: rho_of_E
    real(real64), dimension(:), intent(in) :: bin_edges
    real(real64), dimension(:), intent(inout) :: mpi_wl_hist, wl_logdos
    real(real64), intent(in) :: wl_f
    logical, intent(in) :: rho_saved
    integer, intent(inout) :: radial_mc_steps

    integer, dimension(4) :: rdm1, rdm2
    real(real64) :: e_swapped, e_unswapped, delta_e
    integer :: i, ibin, jbin, iradial
    integer(int16) :: site1, site2

    ! Establish total energy before any moves
    e_unswapped = setup%full_energy(config)
    e_swapped = e_unswapped

    do i = 1, wl_setup%mc_sweeps*setup%n_atoms

      ! Make one MC trial
      ! Generate random numbers
      rdm1 = setup%rdm_site()
      rdm2 = setup%rdm_site()

      ! Get what is on those sites
      site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
      site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

      call pair_swap(config, rdm1, rdm2)

      ! Calculate energy if different species
      if (site1 /= site2) then
        e_swapped = setup%full_energy(config)
      end if

      ibin = bin_index(e_unswapped, bin_edges, wl_setup%bins)
      jbin = bin_index(e_swapped, bin_edges, wl_setup%bins)

      ! Calculate radial density and add to appropriate location in array
      if (jbin < wl_setup%bins + 1 .and. jbin > 0) then
        if (rho_saved .eqv. .False.) then
          if (radial_record(jbin) < MAX(wl_setup%radial_samples/num_walkers,1)) then
            radial_mc_steps = radial_mc_steps + 1
            if(radial_mc_steps >= setup%n_atoms) then
              radial_mc_steps = 0
              radial_record(jbin) = radial_record(jbin) + 1
              rho_of_E(:,:,:,jbin) = rho_of_E(:,:,:,jbin) + radial_densities(setup, config, setup%wc_range, shells)
            end if
          end if
        end if
      end if

      ! Only compute energy change if within limits where V is defined and within MPI region
      if (jbin > mpi_start_idx - 1 .and. jbin < mpi_end_idx + 1) then

        ! Add change in V into diff_energy
        delta_e = e_swapped - e_unswapped

        ! Accept or reject move
        if (genrand() .lt. exp((wl_logdos(ibin) - wl_logdos(jbin)))) then
          e_unswapped = e_swapped
        else
          call pair_swap(config, rdm1, rdm2)
          jbin = ibin
        end if
        mpi_wl_hist(jbin - mpi_start_idx + 1) = mpi_wl_hist(jbin - mpi_start_idx + 1) + 1.0_real64
        wl_logdos(jbin) = wl_logdos(jbin) + wl_f
      else
        ! reject and reset
        call pair_swap(config, rdm1, rdm2)
      end if
    end do

  end subroutine sweeps

  subroutine burn_in(setup, config, min_e, max_e, &
                        mpi_processes, num_walkers, my_rank, mpi_index, window_rank_index, comms)
    integer(int16), dimension(:, :, :, :) :: config
    class(run_params), intent(in) :: setup
    real(real64), intent(in) :: min_e, max_e
    integer, intent(in) :: mpi_processes, num_walkers, my_rank, mpi_index
    integer, dimension(:,:), intent(in) :: window_rank_index
    logical, intent(in) :: comms

    integer, dimension(4) :: rdm1, rdm2
    real(real64) :: e_swapped, e_unswapped, delta_e, target_energy
    integer(int16) :: site1, site2
    logical :: stop_burn_in, flag
    integer :: rank, rank_index, request, ierror, status

    stop_burn_in = .False.
    ! Target energy
    target_energy = (min_e + max_e)/2

    ! Establish total energy before any moves
    e_unswapped = setup%full_energy(config)

    ! Non-blocking MPI receive
    if (comms .eqv. .True.) then
      call MPI_IRECV(stop_burn_in, 1, MPI_LOGICAL, MPI_ANY_SOURCE, 10000, MPI_COMM_WORLD, request, ierror)
    end if

    do while(.True.)
      ! Check if MPI message received
      if (comms .eqv. .True.) then
        call MPI_TEST(request, flag, MPI_STATUS_IGNORE, ierror)
      end if

      ! Stop burn if other rank in window is burnt in
      ! or if burnt in send configuration to rest of window
      if (stop_burn_in .eqv. .True.) then
        call MPI_RECV(config, SIZE(config), MPI_SHORT, MPI_ANY_SOURCE, 10001, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror) ! Check if you can put an array of accepted values in the source variable
        exit
      else if (e_unswapped > min_e .and. e_unswapped < max_e) then
        if (comms .eqv. .True.) then
          if (flag .eqv. .False.) then
            stop_burn_in = .True.
            call MPI_CANCEL(request, ierror)
            call MPI_REQUEST_FREE(request, ierror)
            do rank=window_rank_index(mpi_index, 1), window_rank_index(mpi_index, 2)
              if (rank /= my_rank) then
                call MPI_SEND(stop_burn_in, 1, MPI_LOGICAL, rank, 10000, MPI_COMM_WORLD, ierror)
                call MPI_SEND(config, SIZE(config), MPI_SHORT, rank, 10001, MPI_COMM_WORLD, ierror)
              end if
            end do
          end if
        end if
        exit
      end if

      ! Make one MC trial
      ! Generate random numbers
      rdm1 = setup%rdm_site()
      rdm2 = setup%rdm_site()

      ! Get what is on those sites
      site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
      site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

      if (site1 /= site2) then
        call pair_swap(config, rdm1, rdm2)
        e_swapped = setup%full_energy(config)

        delta_e = e_swapped - e_unswapped

        ! Accept or reject move
        if (e_swapped > target_energy .and. delta_e < 0) then
          e_unswapped = e_swapped
        else if (e_swapped < target_energy .and. delta_e > 0) then
          e_unswapped = e_swapped
        else if (genrand() .lt. 0.001_real64) then ! to prevent getting stuck in local minimum (should adjust this later to something more scientific instead of an arbitrary number)
          e_unswapped = e_swapped
        else
          call pair_swap(config, rdm1, rdm2)
        end if
      end if
    end do
  end subroutine burn_in

  subroutine divide_range(my_rank, wl_setup, mpi_process, window_intervals, window_rank_index, num_walkers)
    ! Subroutine input
    type(wl_params) :: wl_setup
    integer, intent(in) :: my_rank, mpi_process

    ! Subroutine output
    integer, dimension(:, :), intent(out) :: window_intervals, window_rank_index
    integer, intent(out) :: num_walkers

    ! Loop indices
    integer :: i

    ! Calculation variables
    real(real64) :: b, n, power, factor
    ! Set first values outside of loop
    window_intervals(1,1) = 1
    window_intervals(wl_setup%num_windows,2) = wl_setup%bins

    ! values for equation to generate bin distribution
    ! factor derived from equation of form:
    ! y = Ax^B + C
    power = 2
    b = wl_setup%bins
    n = wl_setup%num_windows + 1 ! +1 since "edges" are being obtained for window_intervals array

    factor = (b-1.0_real64)/((n+1.0_real64)**power-1.0_real64)

    do i = 2, wl_setup%num_windows
      window_intervals(i-1,2) = INT(FLOOR(factor*(i**power-1)+1))
      window_intervals(i,1) = INT(FLOOR(factor*(i**power-1)+1) + 1)
    end do

    ! Distribute MPI processes
    window_rank_index(1,1) = 0
    window_rank_index(wl_setup%num_windows, 2) = mpi_process - 1
    do i = 2, wl_setup%num_windows
      window_rank_index(i,1) = INT(FLOOR(REAL(mpi_process/wl_setup%num_windows)))*(i-1)
      window_rank_index(i-1,2) = window_rank_index(i,1) - 1
    end do
  end subroutine divide_range

  subroutine create_window_intervals(wl_setup, num_walkers, window_intervals, window_indices, &
    mpi_index, mpi_start_idx, mpi_end_idx, mpi_bins)
    class(wl_params), intent(in) :: wl_setup
    integer, dimension(:,:), intent(in) :: window_intervals
    integer, intent(in) :: num_walkers

    integer, intent(out) :: mpi_start_idx, mpi_end_idx, mpi_bins, mpi_index
    integer, dimension(:,:), intent(inout) :: window_indices

    integer :: i

    window_indices(1, 1) = window_intervals(1,1)
    window_indices(1,2) = INT(window_intervals(1,2) + wl_setup%bin_overlap)
    do i = 2, wl_setup%num_windows-1
      window_indices(i, 1) = MAX(INT(window_intervals(i,1) - wl_setup%bin_overlap), 1)
      window_indices(i, 2) = MIN(INT(window_intervals(i,2) + wl_setup%bin_overlap), wl_setup%bins)
    end do
    window_indices(wl_setup%num_windows, 1) = INT(window_intervals(wl_setup%num_windows,1) &
                                    - ABS(window_intervals(wl_setup%num_windows-1,1) &
                                    - wl_setup%bin_overlap))
    window_indices(wl_setup%num_windows,2) = window_intervals(wl_setup%num_windows,2)

    mpi_index = my_rank/num_walkers + 1
    mpi_start_idx = window_indices(mpi_index, 1)
    mpi_end_idx = window_indices(mpi_index, 2)
    mpi_bins = mpi_end_idx - mpi_start_idx + 1
  end subroutine create_window_intervals

  subroutine create_energy_bins(setup, wl_setup, mpi_start_idx, mpi_end_idx, bin_edges, bin_energy, mpi_bin_edges)
    class(run_params), intent(in) :: setup
    class(wl_params), intent(in) :: wl_setup
    integer, intent(in) :: mpi_start_idx, mpi_end_idx

    real(real64), dimension(:), intent(inout) :: bin_edges, bin_energy, mpi_bin_edges
    
    ! Internal
    integer :: i, j
    real(real64) :: energy_to_ry, bin_width

    ! Conversion meV/atom to Rydberg
    energy_to_ry = setup%n_atoms/(eV_to_Ry*1000)
    ! Create energy wl_setup%bins and set mpi wl_setup%bins
    j = 1
    bin_width = (wl_setup%energy_max - wl_setup%energy_min)/real(wl_setup%bins)*energy_to_ry
    do i = 1, wl_setup%bins + 1
      bin_edges(i) = wl_setup%energy_min*energy_to_ry + (i - 1)*bin_width
    end do
    do i = 1, wl_setup%bins
      bin_energy(i) = wl_setup%energy_min*energy_to_ry + (i - 0.5)*bin_width
    end do
    do i = mpi_start_idx, mpi_end_idx + 1
      mpi_bin_edges(j) = bin_edges(i)
      j = j + 1
    end do

  end subroutine create_energy_bins

  subroutine dos_average(wl_setup, num_walkers, mpi_index, wl_logdos, wl_logdos_buffer)
    class(wl_params), intent(in) :: wl_setup
    integer, intent(in) :: num_walkers, mpi_index

    real(real64), dimension(:), intent(inout) :: wl_logdos, wl_logdos_buffer

    integer :: i, j, ierror

    do i = 1, wl_setup%num_windows
      if (mpi_index == i) then
        if (my_rank /= (i - 1)*num_walkers) then
          call MPI_Send(wl_logdos, wl_setup%bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers, i, MPI_COMM_WORLD, ierror)
          !print*, my_rank, "send", (i - 1)*num_walkers
        end if
        if (my_rank == (i - 1)*num_walkers) then
          do j = 1, num_walkers - 1
            call MPI_Recv(wl_logdos_buffer, wl_setup%bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers + j, i, MPI_COMM_WORLD, &
                          MPI_STATUS_IGNORE, ierror)
            !print*, my_rank, "recv", (i - 1)*num_walkers + j
            wl_logdos = wl_logdos + wl_logdos_buffer
          end do
        end if
        if (my_rank == (i - 1)*num_walkers) then
          wl_logdos = wl_logdos/num_walkers
          do j = 1, num_walkers - 1
            call MPI_Send(wl_logdos, wl_setup%bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers + j, i, MPI_COMM_WORLD, &
                          ierror)
            !print*, my_rank, "send", (i - 1)*num_walkers + j
          end do
        else
          call MPI_Recv(wl_logdos_buffer, wl_setup%bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers, i, MPI_COMM_WORLD, &
                        MPI_STATUS_IGNORE, ierror)
          !print*, my_rank, "recv", (i - 1)*num_walkers
        end if
      end if
    end do
  end subroutine dos_average

  subroutine dos_combine(wl_setup, mpi_index, my_rank, num_walkers, &
    window_indices, window_rank_index, &
    wl_logdos, wl_logdos_buffer, wl_logdos_write)

    ! Input
    type(wl_params), intent(in) :: wl_setup
    integer, intent(in) :: mpi_index, my_rank, num_walkers
    integer, dimension(:,:), intent(in) :: window_indices, window_rank_index

    ! Input-Output
    real(real64), dimension(:), intent(inout) :: wl_logdos, wl_logdos_buffer, wl_logdos_write

    ! Internal
    integer :: i, j, ierror, beta_index
    real(real64) :: beta_original, beta_merge, beta_diff, scale_factor
    integer :: mpi_start_idx, mpi_end_idx

    if (my_rank == 0) then
    wl_logdos_write = wl_logdos
    end if

    do i = 2, wl_setup%num_windows
      if (my_rank == window_rank_index(i,1)) then
        call MPI_Send(wl_logdos, wl_setup%bins, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierror)
        call MPI_Send(window_indices(mpi_index,1), 1, MPI_INT, 0, 1, MPI_COMM_WORLD, ierror)
        call MPI_Send(window_indices(mpi_index,2), 1, MPI_INT, 0, 2, MPI_COMM_WORLD, ierror)
      end if
      if (my_rank == 0) then
        call MPI_Recv(wl_logdos_buffer, wl_setup%bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers, 0, MPI_COMM_WORLD, &
          MPI_STATUS_IGNORE, ierror)
        call MPI_Recv(mpi_start_idx, 1, MPI_INT, (i - 1)*num_walkers, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        call MPI_Recv(mpi_end_idx, 1, MPI_INT, (i - 1)*num_walkers, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        scale_factor = 0.0_real64
        beta_diff = HUGE(beta_diff)
        do j = 0,  window_indices(i-1, 2) - window_indices(i, 1) - 1
          beta_original = wl_logdos_write(mpi_start_idx + j + 1) - wl_logdos_write(mpi_start_idx + j)
          beta_merge = wl_logdos_buffer(mpi_start_idx + j + 1) - wl_logdos_buffer(mpi_start_idx + j)
          if (ABS(beta_original - beta_merge) < beta_diff) then
            beta_diff = ABS(beta_original - beta_merge)
            beta_index = mpi_start_idx + j
          end if
        end do

        do j = beta_index, mpi_end_idx
          wl_logdos_write(j) = wl_logdos_buffer(j) + wl_logdos_write(beta_index) - wl_logdos_buffer(beta_index)
        end do
      end if
    end do
end subroutine dos_combine

  subroutine mpi_window_optimise(wl_setup, my_rank, bin_edges, window_intervals, window_indices, &
                                mpi_bins, mpi_index, mpi_bin_edges, mpi_wl_hist, rank_all_time, num_walkers)
    ! Input
    type(wl_params) :: wl_setup
    integer, intent(in) :: my_rank, num_walkers
    real(real64), dimension(:,:), intent(inout) :: rank_all_time
    real(real64), dimension(:), intent(in) :: bin_edges

    ! Output
    integer, dimension(:,:), intent(inout) :: window_indices

    ! Input-Output
    integer, dimension(:,:), intent(inout) :: window_intervals
    integer, intent(inout) :: mpi_bins, mpi_index
    real(real64), dimension(:), allocatable, intent(inout) :: mpi_bin_edges, mpi_wl_hist

    ! Internal
    integer :: i, ierror, bins_to_subtract, sum_below_2
    real(real64) :: diffusion(wl_setup%num_windows)
    integer :: bins(wl_setup%num_windows)

    ! Perform window size adjustment then broadcast
    if (my_rank == 0) then
      diffusion = (window_intervals(:,2) - window_intervals(:,1) + 1)**2/rank_all_time(:,1)
      bins = NINT(REAL(wl_setup%bins)*diffusion/SUM(diffusion))
      
      ! Compute bins to subtract
      sum_below_2 = 0
      do i = 1, wl_setup%num_windows
          if (bins(i) < 2) then
              sum_below_2 = sum_below_2 + bins(i)
          end if
      end do
      bins_to_subtract = count(bins < 2) * 2 - sum_below_2
      
      ! Set all bins less than 2 to 2
      do i = 1, wl_setup%num_windows
          if (bins(i) < 2) then
              bins(i) = 2
          end if
      end do
      
      ! Adjust the max bin
      bins(maxloc(bins)) = bins(maxloc(bins)) - bins_to_subtract

      print*, "Old"
      print*, window_intervals(:,1)
      print*, window_intervals(:,2)

      window_intervals(1, 2) = bins(1)
      do i=2, wl_setup%num_windows
        window_intervals(i, 1) = window_intervals(i-1, 2) + 1
        window_intervals(i, 2) = window_intervals(i, 1) + bins(i) - 1
      end do
      window_intervals(wl_setup%num_windows, 1) = window_intervals(wl_setup%num_windows-1, 2) + 1
      window_intervals(wl_setup%num_windows, 2) = wl_setup%bins
      print*, "New"
      print*, window_intervals(:,1)
      print*, window_intervals(:,2)
      !print*, 1+time_mult_array
    end if

    call MPI_BCAST(window_intervals, wl_setup%num_windows*2, MPI_INT, 0, MPI_COMM_WORLD, ierror)

    ! Populate MPI arrays and indlude MPI window overlap
    call mpi_arrays(wl_setup, my_rank, bin_edges, window_intervals, window_indices, &
                    mpi_bin_edges, mpi_wl_hist, mpi_bins, mpi_index, num_walkers)
  end subroutine mpi_window_optimise

  subroutine mpi_arrays(wl_setup, my_rank, bin_edges, window_intervals, window_indices, &
                        mpi_bin_edges, mpi_wl_hist, mpi_bins, mpi_index, num_walkers)
    ! Subroutine input
    type(wl_params) :: wl_setup
    integer, intent(in) :: my_rank, num_walkers
    integer, dimension(:,:), intent(in) :: window_intervals
    real(real64), dimension(:), intent(in) :: bin_edges

    ! Subroutine output
    integer, dimension(:,:), intent(inout) :: window_indices
    real(real64), dimension(:), allocatable, intent(inout) :: mpi_bin_edges, mpi_wl_hist
    integer, intent(out) :: mpi_bins, mpi_index

    ! Loop indices
    integer :: i, j

    deallocate(mpi_bin_edges)
    deallocate(mpi_wl_hist)

    window_indices(1, 1) = window_intervals(1,1)
    window_indices(1,2) = INT(window_intervals(1,2) + wl_setup%bin_overlap)
    do i = 2, wl_setup%num_windows-1
      window_indices(i, 1) = MAX(INT(window_intervals(i,1) - wl_setup%bin_overlap), 1)
      window_indices(i, 2) = MIN(INT(window_intervals(i,2) + wl_setup%bin_overlap), wl_setup%bins)
    end do
    window_indices(wl_setup%num_windows, 1) = INT(window_intervals(wl_setup%num_windows,1) &
                                    - ABS(window_intervals(wl_setup%num_windows-1,1) &
                                    - wl_setup%bin_overlap))
    window_indices(wl_setup%num_windows,2) = window_intervals(wl_setup%num_windows,2)

    mpi_index = my_rank/num_walkers + 1
    mpi_bins = window_indices(mpi_index,2) - window_indices(mpi_index,1) + 1

    allocate(mpi_bin_edges(mpi_bins + 1))
    allocate(mpi_wl_hist(mpi_bins))
    mpi_bin_edges = 0.0_real64
    mpi_wl_hist = 0.0_real64

    j = 1
    do i = window_indices(mpi_index, 1), window_indices(mpi_index, 2) + 1
      mpi_bin_edges(j) = bin_edges(i)
      j = j + 1
    end do
  end subroutine mpi_arrays

end module wang_landau
