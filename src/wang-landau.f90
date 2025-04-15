!> @file    wang-landau.f90                                  
!>
!> @brief   Assorted routines and tools to perform Wang Landau Sampling                                                                      !
!>
!> @details This module contains routines necessary for the Wang Landau Sampling
!>          calculations.                                                           
!>
!> @author  H. J. Naguszewski
!> @date    2024

module wang_landau
  use initialise
  use kinds
  use shared_data
  use c_functions
  use random_site
  use metropolis
  use mpi
  use display

  implicit none
 
  ! MPI variables
  integer :: mpi_processes
  real(real64) :: start, end, time_max, time_min, test_time_1, test_time_2
  integer :: mpi_bins, mpi_start_idx, mpi_end_idx, mpi_index, mpi_start_idx_buffer
  integer :: mpi_end_idx_buffer, beta_index, i_sweeps
  real(real64) :: scale_factor, scale_count, wl_logdos_min, bin_overlap, beta_diff, beta_original, beta_merge
  real(real64), allocatable :: mpi_bin_edges(:), mpi_wl_hist(:), wl_logdos_buffer(:), wl_logdos_combine(:)
  real(real64), allocatable :: rank_time(:), rank_time_buffer(:,:)

  logical :: flag
  
  ! Window variables
  integer, allocatable :: window_intervals(:,:), window_rank_index(:,:)
  integer, allocatable :: window_indices(:, :)
  integer :: num_walkers
  
  ! Setup types
  type(run_params) :: setup_internal
  type(wl_params) :: wl_setup_internal
  
  ! Temperature and temperature steps
  real(real64) :: temp, beta
  
  ! WL variables and arrays
  integer :: bins, resets, pre_sampled_state, ibin, itemp
  real(real64) :: bin_width, energy_to_ry, wl_f, wl_f_prev, tolerance, flatness_tolerance
  real(real64) :: target_energy, flatness, bins_buffer, bins_min, radial_min, radial_min_buffer
  real(real64), allocatable :: bin_edges(:), wl_hist(:), wl_logdos(:), bin_energy(:), mean_energy(:,:), prob(:)
  integer, allocatable :: pre_sampled(:), pre_sampled_buffer(:)
  logical :: rho_saved
  
  ! Load balancing metrics
  real(real64), allocatable :: lb_bins(:,:), lb_avg_time(:,:), lb_max_time(:,:), window_time(:), diffusion_prev(:)
  integer :: converged, converged_sum
  
  ! Radial density across energy
  real(real64), allocatable :: rho_of_E(:,:,:,:), rho_of_E_buffer(:,:,:,:)
  integer, allocatable :: radial_record(:), radial_record_buffer(:)
  logical, allocatable :: radial_record_bool(:)
  integer :: radial_mc_steps
  real(real64) :: radial_time
  
  ! Rank of this processor
  integer :: ierror, request
  
  ! Loop integers and error handling variable
  integer :: i, j
  integer :: iter, num_iter
  
  ! Name of file for writing radial densities at the end
  character(len=37) :: radial_file
  
  contains

  !> @brief   Main Wang Landau sampling routine.  
  !>
  !> @details This routine performs the Wang Landau sampling calculation. Input parameters (such as
  !>          the energy range, number of windows etc.) are read from the "wl_input.txt" 
  !>          file. The first section of the routine generates the initial random configurations, moves them
  !>          into approriate energy windows, performs pre-sampling to get an estimate for the density of states
  !>          before performing the main Wang Landau sampling cycle. The main cycle consists of converging
  !>          
  !>          The routine makes use of the pair_swap subroutine to calculate energies.
  !> 
  !> @param   setup Derived type containing simulation parameters
  !> @param   wl_setup Derived type containing Wang Landau sampling parameters
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine wl_main(setup, wl_setup)
    type(run_params) :: setup
    type(wl_params) :: wl_setup

    ! Makes internal copies that use module wide variables
    setup_internal = setup
    wl_setup_internal = wl_setup

    ! Check if number of MPI processes is divisible by number of windows
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_processes, ierror)
    bins = wl_setup_internal%bins
    if (MOD(mpi_processes, wl_setup_internal%num_windows) /= 0) then
      if (my_rank == 0) then
        write (6, '(72("~"))')
        write (6, '(5("~"),x,"Error: Number of MPI processes not divisible by num_windows",x,6("~"))')
        write (6, '(72("~"))')
      end if
      call MPI_FINALIZE(ierror)
      call EXIT(0)
    end if

    ! Initial setup
    call primary_array_allocation()

    call divide_range(window_intervals, window_rank_index)

    call create_window_intervals(window_intervals, window_indices, mpi_bins)

    call secondary_array_allocation(mpi_bins)

    call create_energy_bins(bin_edges)

    call initialise_variables()

    if (my_rank == 0) then
      write (6, '(/,72("-"),/)')
      call print_centered_message("Commencing Simulation!", "=")
      write(*,*)
      print *, "Number of atoms", setup_internal%n_atoms
      print *, "Number of iterations", num_iter
    end if
    call comms_wait()

    !---------!
    ! Burn in !
    !---------!
    
    call enter_energy_window()
    print*, "Rank: ", my_rank, "within energy window"
  
    flush(6)
    call comms_wait()
    if (my_rank == 0) then
      write (*, *)
      call print_centered_message("Ranks Within Energy Windows!", "=")
      write (*, *)
    end if

    !--------------------!
    ! Pre-Sampling       !
    !--------------------!
    call pre_sampling(wl_logdos, mpi_wl_hist)

    flush(6)
    call comms_wait()
    if (my_rank == 0) then
      call print_centered_message("Pre-sampling Complete!", "=")
      write (*, *)
    end if

    !--------------------!
    ! Main Wang-Landau   !
    !--------------------!
    iter = 0
    converged = 0
    converged_sum = 0
    radial_time = 0.0_real64
    mpi_wl_hist = 0.0_real64
    i = 0
    i_sweeps = 0
    start = mpi_wtime()
    do while (wl_f > wl_setup_internal%tolerance)
      i_sweeps = i_sweeps + 1
      if (MOD(i_sweeps, 1000) == 0) then
        i_sweeps = 0
        if (my_rank == 0) then
          print*, ""
          print*, "1000 Iterations performed. Converged ranks: ", REAL(converged_sum)/REAL(mpi_processes)
        end if
        call comms_wait()
        if (converged == 0) then
          print*, "Unconverged rank: ", my_rank, "Flatness: ", flatness, "Explored: ", &
          REAL(COUNT(mpi_wl_hist/=0))/REAL(SIZE(mpi_wl_hist))
        end if
      end if 
      call sweeps(wl_logdos, mpi_wl_hist)
      if (MOD(i_sweeps, 10) == 0) then        
        call replica_exchange(config)
      end if

      flatness = minval(mpi_wl_hist)/(sum(mpi_wl_hist)/mpi_bins)

      if (converged == 0 .and. flatness > wl_setup_internal%flatness .and. minval(mpi_wl_hist) > 10.0_real64) then
        ! End timer
        end = mpi_wtime()
        converged = 1
      end if

      call MPI_ALLREDUCE(converged, converged_sum, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)

      if (converged_sum == mpi_processes) then
        i = 0
        i_sweeps = 0
        converged = 0
        converged_sum = 0
        !Reset the histogram
        mpi_wl_hist = 0.0_real64
        ! Increase iteration
        iter = iter + 1
        !Reduce f
        wl_f = wl_f*1/2
        wl_logdos = wl_logdos - minval(wl_logdos, MASK=(wl_logdos > 0.0_real64))!wl_logdos_min
        wl_logdos = ABS(wl_logdos * merge(0, 1, wl_logdos < 0.0_real64))
        !Average DoS across walkers
        call dos_average(wl_logdos)
        call dos_combine(wl_logdos)

        call reduce_time(start, end, radial_time)

        call comms_wait()
        ! MPI send and recieve calls for combining window DoS
        
        if (my_rank == 0) then
          write (6, '(20("-"),x,a,i3,a,i3,x,20("-"))', advance='no') "Wang-Landau Iteration: ", iter, &
          "/", num_iter
          write (*, *)
          write (6, '(a,f20.18,a,f8.2,a)', advance='no') "Flatness reached f of: ", wl_f_prev, &
                  " | Radial samples: ", radial_min*100_real64, "%"
            write (*, *)
          do i=1, wl_setup_internal%num_windows
            write (6, '(a,i3,a,f12.2,a,f12.2,a,f12.2,a)') "MPI Window: ", i, " | Avg. time: ", rank_time_buffer(i,1), &
            "s | Time min: ", rank_time_buffer(i,2), "s Time max: " , rank_time_buffer(i,3), "s"
          end do
          wl_f_prev = wl_f
        end if

        call save_rho_E(rho_saved, radial_record)

        call save_wl_data(bin_edges, wl_logdos, wl_hist)
        call save_load_balance_data(window_indices, rank_time_buffer)
        
        call mpi_window_optimise(iter)

        call compute_mean_energy(wl_logdos)
        call enter_energy_window()

        call zero_subtract_logdos(wl_logdos)
        call comms_wait()
        if (my_rank == 0) then
          call print_centered_message("Load Balancing Complete!", "-")
          write (*, *)
        end if
        start = mpi_wtime()
      end if
    end do

    if (my_rank == 0) then
      write (*, *)
      call print_centered_message("Simulation Complete!", "=")
    end if

  end subroutine wl_main

  !> @brief   Reducing MPI timing across processes
  !>          
  !> @param   setup Derived type containing simulation parameters
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine reduce_time(start, end, radial_time)
    real(real64), intent(in) :: start, end
    real(real64), intent(inout) :: radial_time
    rank_time = 0.0_real64
    rank_time(mpi_index) = end - start - radial_time
    call MPI_ALLREDUCE(end - start, time_max, 1, MPI_DOUBLE_PRECISION, &
    MPI_MAX, MPI_COMM_WORLD, ierror)
    call MPI_REDUCE(end - start, time_min, 1, MPI_DOUBLE_PRECISION, &
    MPI_MIN, 0, MPI_COMM_WORLD, ierror)
    call MPI_REDUCE(rank_time/num_walkers, rank_time_buffer(:,1), wl_setup_internal%num_windows, MPI_DOUBLE_PRECISION, &
    MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    call MPI_REDUCE(rank_time, rank_time_buffer(:,3), wl_setup_internal%num_windows, MPI_DOUBLE_PRECISION, &
    MPI_MAX, 0, MPI_COMM_WORLD, ierror)
    rank_time = time_max
    rank_time(mpi_index) = end - start - radial_time
    radial_time = 0.0_real64
    call MPI_REDUCE(rank_time, rank_time_buffer(:,2), wl_setup_internal%num_windows, MPI_DOUBLE_PRECISION, &
    MPI_MIN, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(rank_time_buffer, wl_setup_internal%num_windows*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine reduce_time

  !> @brief   Saving radial density as a function of energy
  !>
  !> @details Routine that checks if the requisite number of radial density samples
  !>          has been recorded and then saves to output file. Once data has been saved
  !>          the routine is not executed.
  !>          
  !> @param   rho_saved Boolean that indicates whether radial density has been written to file
  !> @param   radial_record Integer 1D array Indivial MPI process storage array of how many radial density samples
  !>          have been taken in each energy bin
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine save_rho_E(rho_saved, radial_record)
    logical, intent(inout) :: rho_saved
    integer, intent(in) :: radial_record(:)
    
    if (.not. rho_saved) then
      call MPI_REDUCE(radial_record, radial_record_buffer, wl_setup_internal%bins, MPI_INT, &
      MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      if (my_rank == 0) then
        radial_min = 0
        do i=1,wl_setup_internal%bins
          radial_min = radial_min + REAL(MIN(radial_record_buffer(i),wl_setup_internal%radial_samples))
        end do
        do i=1,wl_setup_internal%bins
          if (radial_record_buffer(i) >= wl_setup_internal%radial_samples) then
            radial_record_bool(i) = .True.
          end if
        end do
        radial_min = radial_min/REAL(wl_setup_internal%radial_samples*wl_setup_internal%bins)
      end if

      call MPI_BCAST(radial_min, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_BCAST(radial_record_bool, wl_setup_internal%bins, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

      if (radial_min >= 1) then
        rho_saved = .True.
        call MPI_REDUCE(rho_of_E, rho_of_E_buffer, &
        setup_internal%n_species*setup_internal%n_species*setup_internal%wc_range*wl_setup_internal%bins, &
        MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        if (my_rank == 0) then 
            call print_centered_message("Radial Densities Saved", "=")
            write (*, *)
          do i=1, wl_setup_internal%bins
            rho_of_E_buffer(:,:,:,i) = rho_of_E_buffer(:,:,:,i)/REAL(radial_record_buffer(i))
          end do
          call ncdf_radial_density_writer_across_energy(radial_file, rho_of_E_buffer, shells, bin_energy, setup_internal)
        end if
      end if
    end if
  end subroutine save_rho_E

  !> @brief   Save Wang Landau data
  !>
  !> @details Routine that saves the energy binning, density of states and hit histogram
  !>          
  !> @param   bin_edges Double 1D array containing energy bin edges
  !> @param   wl_logdos Double 1D array containing density of states
  !> @param   wl_hist Double 1D array contining histogram of energy hits
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine save_wl_data(bin_edges, wl_logdos, wl_hist)
    real(real64), allocatable, intent(in) :: bin_edges(:), wl_logdos(:), wl_hist(:)
    if (my_rank == 0) then
      ! Write output files
      call ncdf_writer_1d("wl_dos_bins.dat", ierror, bin_edges)
      call ncdf_writer_1d("wl_dos.dat", ierror, wl_logdos)
      call ncdf_writer_1d("wl_hist.dat", ierror, wl_hist)
    end if
  end subroutine save_wl_data

  !> @brief   Save load balancing data
  !>
  !> @details Routine that saves data relevant to performance analysis.
  !>          
  !> @param   window_indices Integer 2D array containing energy bin edges
  !> @param   rank_time_buffer Double 2D array containing density of states
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine save_load_balance_data(window_indices, rank_time_buffer)
    integer, intent(in) :: window_indices(:, :)
    real(real64), intent(in) :: rank_time_buffer(:, :)
    if (my_rank == 0) then
      lb_bins(iter, :) = REAL(window_indices(:, 2) - window_indices(:, 1) + 1)
      lb_avg_time(iter, :) = rank_time_buffer(:,1)
      lb_max_time(iter, :) = rank_time_buffer(:,3)
      window_time(iter) = MAXVAL(rank_time_buffer(:,3))
      call ncdf_writer_2d("wl_lb_bins.dat", ierror, lb_bins)
      call ncdf_writer_2d("wl_lb_avg_time.dat", ierror, lb_avg_time)
      call ncdf_writer_2d("wl_lb_max_time.dat", ierror, lb_max_time)
      call ncdf_writer_1d("wl_window_time.dat", ierror, window_time)
    end if
  end subroutine save_load_balance_data

  !> @brief   Computes the mean energy of the system
  !>
  !> @details Routine that computes the mean energy of the sytem.
  !>          
  !> @param   wl_logdos Double 1D array containing density of states
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine compute_mean_energy(wl_logdos)
    real(real64), allocatable, intent(in) :: wl_logdos(:)
    wl_logdos_buffer = wl_logdos - maxval(wl_logdos)
    do itemp = 1, 300
      beta = 1.0_real64 / (k_b_in_Ry * itemp * 10.0_real64)

      ! Compute reweighted histogram
      do ibin = 1, wl_setup_internal%bins
          prob(ibin) = wl_logdos_buffer(ibin) - beta * (bin_edges(ibin) + 0.5_real64 * bin_width)
      end do

      ! Normalize probabilities
      prob = exp(prob - maxval(prob))
      prob = prob / sum(prob)

      ! Compute mean energy
      mean_energy(itemp, 1) = sum((bin_edges(:wl_setup_internal%bins) + 0.5_real64 * bin_width) * prob)
      mean_energy(itemp, 2) = beta
    end do
    wl_logdos_buffer = wl_logdos
  end subroutine compute_mean_energy

  !> @brief   Scales density of states for each process
  !>
  !> @details Routine that takes the density of states array and subtracts the minimum value
  !>          within the array from the array and then zeroes elements beyond energy window
  !>          of process
  !>          
  !> @param   wl_logdos Double 1D array containing density of states
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine zero_subtract_logdos(wl_logdos)
    real(real64), allocatable, intent(inout) :: wl_logdos(:)
    ! Zero elements not worked on
    wl_logdos(1:window_indices(mpi_index,1)-1) = 0.0_real64
    wl_logdos(window_indices(mpi_index,2)+1:wl_setup_internal%bins) = 0.0_real64
    ! Subtract minimum value
    wl_logdos = wl_logdos - minval(wl_logdos, MASK=(wl_logdos > 0.0_real64))
    wl_logdos = ABS(wl_logdos * merge(0, 1, wl_logdos < 0.0_real64))
  end subroutine zero_subtract_logdos

  !> @brief   Finds index of energy bin given energy
  !>
  !> @details Routine that outputs the index of the energy bin of given energy
  !>      
  !> @param   energy Double that contains current energy of system  
  !> @param   bin_edges Double 1D array containing energy bin edges
  !> @param   bins Integer that indicates number of energy bins
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  integer function bin_index(energy, bin_edges, bins) result(index)
    integer, intent(in) :: bins
    real(real64), intent(in) :: energy
    real(real64), dimension(:), intent(in) :: bin_edges
    real(real64) :: bin_range

    bin_range = bin_edges(bins + 1) - bin_edges(1)
    index = int(((energy - bin_edges(1))/(bin_range))*real(bins)) + 1
  end function bin_index

  !> @brief   Main routine that performs Monte Carlo sweeps
  !>
  !> @details The routine performs a user set number of Monte Carlo sweeps and
  !>          records and updates data relevant to the Wang Landau sampling process
  !>          (density of states and hit histogram).
  !>          
  !> @param   wl_logdos Double 1D array containing density of states
  !> @param   mpi_wl_hist Double 1D array containing energy bin visit information
  !>          for working MPI process
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine sweeps(wl_logdos, mpi_wl_hist)
    real(real64), allocatable, intent(inout) :: wl_logdos(:), mpi_wl_hist(:)
    integer, dimension(4) :: rdm1, rdm2
    real(real64) :: e_swapped, e_unswapped, pair_unswapped, pair_swapped, delta_e, radial_start, radial_end
    integer :: i, ibin, jbin, iradial
    integer(int16) :: site1, site2

    ! Establish total energy before any moves
    e_unswapped = setup_internal%full_energy(config)
    e_swapped = e_unswapped

    do i = 1, wl_setup_internal%mc_sweeps*setup_internal%n_atoms
      ! Make one MC trial
      ! Generate random numbers
      rdm1 = setup_internal%rdm_site()
      rdm2 = setup_internal%rdm_site()
      ! Get what is on those sites
      site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
      site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

      e_swapped = e_unswapped

      pair_unswapped = pair_energy(setup_internal, config, rdm1, rdm2)
      pair_swapped = pair_unswapped

      call pair_swap(config, rdm1, rdm2)
      ! Calculate energy if different species
      if (site1 /= site2) then
        pair_swapped = pair_energy(setup_internal, config, rdm1, rdm2)
        e_swapped = e_unswapped - pair_unswapped + pair_swapped
      end if
      ibin = bin_index(e_unswapped, bin_edges, wl_setup_internal%bins)
      jbin = bin_index(e_swapped, bin_edges, wl_setup_internal%bins)
      ! Calculate radial density and add to appropriate location in array
      radial_start = mpi_wtime()
      if (jbin < wl_setup_internal%bins + 1 .and. jbin > 0) then
        if (rho_saved .eqv. .False.) then
          if (radial_record(jbin) < MAX(wl_setup_internal%radial_samples/num_walkers,1)) then
            if (radial_record_bool(jbin) .eqv. .False.) then
              radial_mc_steps = radial_mc_steps + 1
              if(radial_mc_steps >= setup_internal%n_atoms) then
                radial_mc_steps = 0
                radial_record(jbin) = radial_record(jbin) + 1
                rho_of_E(:,:,:,jbin) = rho_of_E(:,:,:,jbin) + &
                radial_densities(setup_internal, config, setup_internal%wc_range, shells)
                radial_end = mpi_wtime()
                if (converged == 0) then
                  radial_time = radial_time + radial_end - radial_start
                end if
              end if
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
        if (MOD(i,INT(0.02_real64*REAL(setup_internal%n_atoms))) == 0) then
          mpi_wl_hist(jbin - mpi_start_idx + 1) = mpi_wl_hist(jbin - mpi_start_idx + 1) + 1.0_real64
        end if
        wl_logdos(jbin) = wl_logdos(jbin) + wl_f
      else
        ! reject and reset
        call pair_swap(config, rdm1, rdm2)
      end if
    end do

  end subroutine sweeps

  !> @brief   Routine for moving processes into assigned energy window
  !>
  !> @details Routine that performs Metropolis Monte Carlo moves in order
  !>          to move process atom configuration into assigned energy window
  !>          
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine enter_energy_window()
    integer, dimension(4) :: rdm1, rdm2
    real(real64) :: e_swapped, e_unswapped, pair_swapped, pair_unswapped, delta_e, target_energy, condition
    real(real64) :: beta, beta_end, beta_start, weight, min_e, max_e
    integer(int16) :: site1, site2
    logical :: stop_enter_energy_window, flag
    integer :: rank, rank_index, request, ierror, status, i_steps, i_sweeps, sweeps

    stop_enter_energy_window = .False.
    flag = .False.
    ! Target energy
    min_e = MINVAL(mpi_bin_edges)
    max_e = MAXVAL(mpi_bin_edges)
    target_energy = (min_e + max_e)/2.0_real64
    condition = ABS(max_e - min_e)*0.1_real64

    !beta = mean_energy(minloc(abs(mean_energy(:,1) - min_e), DIM=1),2)

    ! Establish total energy before any moves
    e_unswapped = setup_internal%full_energy(config)

    if (e_unswapped < target_energy) then
      beta = 1.0_real64/(k_B_in_Ry*30000.0_real64)
      beta_start = beta
      beta_end = beta!mean_energy(minloc(abs(mean_energy(:,1) - min_e), DIM=1),2)
    else 
      beta = 1.0_real64/(k_B_in_Ry*30000.0_real64)
      beta_start = beta
      beta_end = 1.0_real64/(k_B_in_Ry*1.0_real64)
    end if

    ! Non-blocking MPI receive
    call MPI_IRECV(stop_enter_energy_window, 1, MPI_LOGICAL, MPI_ANY_SOURCE, 10000, MPI_COMM_WORLD, request, ierror)

    i_steps = 0
    i_sweeps = 0
    sweeps = 100
    do while(.True.)
      i_steps = i_steps + 1
      if (MOD(i_steps, setup_internal%n_atoms) == 0) then
        i_steps = 0
        i_sweeps = i_sweeps + 1
        weight = MAX(REAL(i_sweeps) / REAL(sweeps), 1.0_real64)
        beta = (1.0_real64 - weight) * beta_start + weight * beta_end
      end if

      ! Check if MPI message received
      call MPI_TEST(request, flag, MPI_STATUS_IGNORE, ierror)

      ! Stop burn if other rank in window is burnt in
      ! or if burnt in send configuration to rest of window
      if (flag .eqv. .True.) then
        call MPI_RECV(config, SIZE(config), MPI_SHORT, MPI_ANY_SOURCE, 10001, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror) ! Check if you can put an array of accepted values in the source variable
        exit
      else if (e_unswapped < max_e-condition .and. e_unswapped > min_e+condition) then
        e_unswapped = setup_internal%full_energy(config)
        if (e_unswapped < max_e-condition .and. e_unswapped > min_e+condition) then
        else
          cycle
        end if
        if (flag .eqv. .False.) then
          stop_enter_energy_window = .True.
          call MPI_CANCEL(request, ierror)
          call MPI_REQUEST_FREE(request, ierror)
          do rank=window_rank_index(mpi_index, 1), window_rank_index(mpi_index, 2)
            if (rank /= my_rank) then
              call MPI_ISEND(stop_enter_energy_window, 1, MPI_LOGICAL, rank, 10000, MPI_COMM_WORLD, request, ierror)
              call MPI_ISEND(config, SIZE(config), MPI_SHORT, rank, 10001, MPI_COMM_WORLD, request, ierror)
            end if
          end do
        end if
        exit
      end if

      ! Make one MC trial
      ! Generate random numbers
      rdm1 = setup_internal%rdm_site()
      rdm2 = setup_internal%rdm_site()

      ! Get what is on those sites
      site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
      site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

      ! Calculate energy if different species
      if (site1 /= site2) then
        pair_unswapped = pair_energy(setup_internal, config, rdm1, rdm2)

        call pair_swap(config, rdm1, rdm2)

        pair_swapped = pair_energy(setup_internal, config, rdm1, rdm2)
        e_swapped = e_unswapped - pair_unswapped + pair_swapped

        delta_e = e_swapped - e_unswapped

        ! Accept or reject move
        if (genrand() .lt. exp(-delta_e*beta)) then ! to prevent getting stuck in local minimum (should adjust this later to something more scientific instead of an arbitrary number)
          e_unswapped = e_swapped
        else
          call pair_swap(config, rdm1, rdm2)
        end if
      end if
    end do

    call comms_wait()
    call comms_purge()
  end subroutine enter_energy_window

  !> @brief   Pre-sampling routine that performs initial sweeps
  !>
  !> @details Routine that performs pre-sampling such that each energy
  !>          bin has been visited and an initial density of states is obtained.
  !>          After sampling performs a window optimisation based on pre-sampling run time.
  !>          
  !> @param   wl_logdos Double 1D array containing density of states
  !> @param   mpi_wl_hist Double 1D array containing energy bin visit information
  !>          for working MPI process
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine pre_sampling(wl_logdos, mpi_wl_hist)
    real(real64), allocatable, intent(inout) :: wl_logdos(:), mpi_wl_hist(:)
    pre_sampled = 0
    pre_sampled_buffer = 0
    pre_sampled_state = 0
    radial_time = 0.0_real64
    start = mpi_wtime()
    end = start
    do while (SUM(pre_sampled_buffer) < wl_setup_internal%num_windows)
      call sweeps(wl_logdos, mpi_wl_hist)

      if (minval(mpi_wl_hist) > REAL(setup_internal%n_atoms) .and. pre_sampled_state == 0) then
        pre_sampled(mpi_index) = 1
      end if
      call MPI_ALLREDUCE(pre_sampled, pre_sampled_buffer, wl_setup_internal%num_windows, &
      MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
      if (pre_sampled_buffer(mpi_index) == 1) then
        if (pre_sampled_state == 0) then
          end = mpi_wtime() - radial_time
          pre_sampled_state = 1
          mpi_wl_hist = REAL(setup_internal%n_atoms) + 1.0_real64
          write (6, '(a,i3,a,f6.2,a)') "Rank: ", INT(my_rank), " | bins visited: ", 100.0_real64, "%"
        end if
      else
        bins_min = REAL(count(mpi_wl_hist > 1.0_real64))/REAL(mpi_bins)
        write (6, '(a,i3,a,f6.2,a)') "Rank: ", INT(my_rank), " | bins visited: ", REAL(bins_min*100.0_real64), "%"
      end if
    end do
    call save_rho_E(rho_saved, radial_record)
    call zero_subtract_logdos(wl_logdos)
    call dos_average(wl_logdos)
    call dos_combine(wl_logdos)
    
    call compute_mean_energy(wl_logdos)

    radial_time = 0.0_real64 ! radial time already accounted for
    call reduce_time(start, end, radial_time)
    
    call mpi_window_optimise(0)

    if (my_rank == 0) then
      print*, window_indices(:,1)
      print*, window_indices(:,2)
    end if
    if (my_rank == 0 ) then
      call print_centered_message("Pre-sampling Timings", "-")
      do i=1, wl_setup_internal%num_windows
        write (6, '(a,i3,a,f12.2,a,f12.2,a,f12.2,a)') "MPI Window: ", i, " | Avg. time: ", rank_time_buffer(i,1), &
        "s | Time min: ", rank_time_buffer(i,2), "s Time max: " , rank_time_buffer(i,3), "s"
      end do
      write (*, *)
    end if

    call enter_energy_window()

    call zero_subtract_logdos(wl_logdos)
  end subroutine pre_sampling

  !> @brief   Creates initial energy windows
  !>
  !> @details Routine that divides energy bins into windows using
  !>          y = Ax^B + C
  !>          
  !>          
  !> @param   window_intervals Integer 2D array containing energy bin indices for
  !>          each window
  !> @param   window_rank_index Integer 2D array containing the window index for
  !>          each MPI process (which window the process is in)
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine divide_range(window_intervals, window_rank_index)
    integer, intent(inout) :: window_intervals(:,:), window_rank_index(:,:)
    ! Loop indices
    integer :: i

    ! Calculation variables
    real(real64) :: b, n, power, factor
    ! Set first values outside of loop
    window_intervals(1,1) = 1
    window_intervals(wl_setup_internal%num_windows,2) = wl_setup_internal%bins

    ! values for equation to generate bin distribution
    ! factor derived from equation of form:
    ! y = Ax^B + C
    power = 1
    b = wl_setup_internal%bins
    n = wl_setup_internal%num_windows

    factor = (b-1.0_real64)/((n+1.0_real64)**power-1.0_real64)

    do i = 2, wl_setup_internal%num_windows
      window_intervals(i-1,2) = INT(FLOOR(factor*((i-1)**power)+1))
      window_intervals(i,1) = window_intervals(i-1,2) + 1
    end do

    if (my_rank == 0) then
      print*, window_intervals(:,2) - window_intervals(:,1) + 1
    end if

    ! Distribute MPI processes
    window_rank_index(1,1) = 0
    window_rank_index(wl_setup_internal%num_windows, 2) = mpi_processes - 1
    do i = 2, wl_setup_internal%num_windows
      window_rank_index(i,1) = INT(FLOOR(REAL(mpi_processes/wl_setup_internal%num_windows)))*(i-1)
      window_rank_index(i-1,2) = window_rank_index(i,1) - 1
    end do
  end subroutine divide_range

  !> @brief   Routine that determines window indices for use by MPI processes
  !>
  !> @details Routine that takes the window intervals and adds the user defined overlap
  !>          region. After doing so it defined the relevant MPI indices.
  !>          
  !> @param   window_intervals Integer 2D array containing energy bin indices for
  !>          each window
  !> @param   window_indices Integer 2D array containing the energy bin indices for
  !>          each MPI process including ovelap
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine create_window_intervals(window_intervals, window_indices, mpi_bins)
    integer, intent(inout) :: window_intervals(:,:), window_indices(:,:), mpi_bins
    integer :: i

    ! Halve overlap due to implementation
    wl_setup_internal%bin_overlap = MAX(wl_setup_internal%bin_overlap/2, 1)

    window_indices(1, 1) = window_intervals(1,1)
    window_indices(1,2) = INT(window_intervals(1,2) + wl_setup_internal%bin_overlap)
    do i = 2, wl_setup_internal%num_windows-1
      window_indices(i, 1) = MAX(INT(window_intervals(i,1) - wl_setup_internal%bin_overlap), 1)
      window_indices(i, 2) = MIN(INT(window_intervals(i,2) + wl_setup_internal%bin_overlap), wl_setup_internal%bins)
    end do
    window_indices(wl_setup_internal%num_windows, 1) = MAX(INT(window_intervals(wl_setup_internal%num_windows,1) &
                                    - wl_setup_internal%bin_overlap), 1)
    window_indices(wl_setup_internal%num_windows,2) = window_intervals(wl_setup_internal%num_windows,2)

    num_walkers = mpi_processes/wl_setup_internal%num_windows

    mpi_index = my_rank/num_walkers + 1
    mpi_start_idx = window_indices(mpi_index, 1)
    mpi_end_idx = window_indices(mpi_index, 2)
    mpi_bins = mpi_end_idx - mpi_start_idx + 1
  end subroutine create_window_intervals

  !> @brief   Routine that creates energy bins
  !>
  !> @details Routine that creates energy bins based on user input.
  !>          Takes minimum and maximum energy and creates a uniformly sized number of bins
  !>          based on user defined number.
  !>          
  !> @param   bin_edges Doulbe 1D array that store energy bin edges
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine create_energy_bins(bin_edges)  
    real(real64), intent(inout) :: bin_edges(:)
    ! Internal
    integer :: i, j
    real(real64) :: energy_to_ry, bin_width

    ! Conversion meV/atom to Rydberg
    energy_to_ry = setup_internal%n_atoms/(Ry_to_eV*1000)

    ! Create energy wl_setup_internal%bins and set mpi wl_setup_internal%bins
    j = 1
    bin_width = (wl_setup_internal%energy_max - wl_setup_internal%energy_min)/real(wl_setup_internal%bins)*energy_to_ry
    do i = 1, wl_setup_internal%bins + 1
      bin_edges(i) = wl_setup_internal%energy_min*energy_to_ry + (i - 1)*bin_width
    end do
    do i = 1, wl_setup_internal%bins
      bin_energy(i) = wl_setup_internal%energy_min*energy_to_ry + (i - 0.5)*bin_width
    end do
    do i = mpi_start_idx, mpi_end_idx + 1
      mpi_bin_edges(j) = bin_edges(i)
      j = j + 1
    end do

  end subroutine create_energy_bins

  !> @brief   First set of array allocations
  !>
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine primary_array_allocation()
    ! Number of WL iteration to be performed
    num_iter = 0
    wl_f = wl_setup_internal%wl_f
    do while (wl_f > wl_setup_internal%tolerance)
      num_iter = num_iter + 1
      wl_f = wl_f*0.5
    end do
    wl_f = wl_setup_internal%wl_f
    allocate(lb_bins(num_iter, wl_setup_internal%num_windows))
    allocate(lb_avg_time(num_iter, wl_setup_internal%num_windows))
    allocate(lb_max_time(num_iter, wl_setup_internal%num_windows))
    allocate(window_time(num_iter))
    allocate(diffusion_prev(wl_setup_internal%num_windows))
    allocate(pre_sampled(wl_setup_internal%num_windows))
    allocate(pre_sampled_buffer(wl_setup_internal%num_windows))

    ! Radial densities as a function of energy
    allocate(rho_of_E(setup_internal%n_species, setup_internal%n_species, setup_internal%wc_range, wl_setup_internal%bins))
    allocate(rho_of_E_buffer(setup_internal%n_species, setup_internal%n_species, setup_internal%wc_range, wl_setup_internal%bins))

    ! Get start and end indices for energy windows
    allocate(window_indices(wl_setup_internal%num_windows, 2))
    allocate(window_intervals(wl_setup_internal%num_windows, 2))
    allocate(window_rank_index(wl_setup_internal%num_windows, 2))
    
    ! allocate arrays
    allocate(radial_record(wl_setup_internal%bins))
    allocate(radial_record_buffer(wl_setup_internal%bins))
    allocate(radial_record_bool(wl_setup_internal%bins))
    allocate(bin_edges(bins + 1))
    allocate(wl_hist(bins))
    allocate(wl_logdos(bins))
    allocate(wl_logdos_buffer(bins))
    allocate(bin_energy(bins))
    allocate(mean_energy(300,2))
    allocate(prob(bins))
  end subroutine primary_array_allocation

  !> @brief   Second set of array allocations
  !>
  !> @details Occurs after primary allocation due to dependacy on mpi_bins variable
  !>          
  !> @param   mpi_bins Integer for number of bins assigned to each MPI process
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine secondary_array_allocation(mpi_bins)
    integer, intent(in) :: mpi_bins
    ! MPI arrays
    allocate(mpi_bin_edges(mpi_bins + 1))
    allocate(mpi_wl_hist(mpi_bins))
    allocate(rank_time(wl_setup_internal%num_windows))
    allocate(rank_time_buffer(wl_setup_internal%num_windows,3))

    if (my_rank == 0) then
      allocate(wl_logdos_combine(bins))
    end if
  end subroutine

  !> @brief   Initialise variables and setup lattice
  !>
  !> @details Intialises variables and arrays. Performs initial setup for lattice configuration
  !>          
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine initialise_variables()
    ! Target energy for burn in
    target_energy = (mpi_bin_edges(1) + mpi_bin_edges(SIZE(mpi_bin_edges)))/2

    ! Load balance diffusion
    diffusion_prev = (window_intervals(:,2) - window_intervals(:,1) + 1)**2

    ! Load balance arrays
    lb_bins = -1.0_real64
    lb_avg_time = -1.0_real64
    lb_max_time = -1.0_real64
    window_time = -1.0_real64

    rho_of_E_buffer = 0.0_real64
    rho_of_E = 0.0_real64

    ! Path to radial file and name of radial file
    radial_file = "radial_densities/rho_of_E.dat"

    num_walkers = mpi_processes/wl_setup_internal%num_windows
    wl_setup_internal%radial_samples = INT(wl_setup_internal%radial_samples/num_walkers)
    wl_setup_internal%radial_samples = MAX(INT(wl_setup_internal%radial_samples*num_walkers), 1)

    ! Initialize
    wl_hist = 0.0_real64; wl_logdos = 0.0_real64; wl_f = wl_setup_internal%wl_f; wl_f_prev = wl_f
    flatness = 0.0_real64;
    mpi_wl_hist = 0.0_real64; wl_logdos_buffer = 0.0_real64; rank_time = 0.0_real64
    radial_record = 0
    radial_record_buffer = 0
    radial_record_bool = .False.
    radial_mc_steps = 0
    rho_saved = .False.

    mean_energy = 1.0_real64 / (k_b_in_Ry * 10.0_real64)

    ! Set up the lattice
    call initial_setup(setup_internal, config)
    call lattice_shells(setup_internal, shells, config)
  end subroutine initialise_variables

  !> @brief   Routine that averages DoS within window
  !>        
  !> @param   wl_logdos Double 1D array containing density of states
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine dos_average(wl_logdos)
    real(real64), intent(inout) :: wl_logdos(:)
    integer :: i, j

    do i = 1, wl_setup_internal%num_windows
      if (mpi_index == i) then
        if (my_rank /= (i - 1)*num_walkers) then
          call MPI_Send(wl_logdos, wl_setup_internal%bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers, i, MPI_COMM_WORLD, ierror)
          !print*, my_rank, "send", (i - 1)*num_walkers
        end if
        if (my_rank == (i - 1)*num_walkers) then
          do j = 1, num_walkers - 1
            call MPI_Recv(wl_logdos_buffer, wl_setup_internal%bins, MPI_DOUBLE_PRECISION, &
            (i - 1)*num_walkers + j, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            !print*, my_rank, "recv", (i - 1)*num_walkers + j
            wl_logdos = wl_logdos + wl_logdos_buffer
          end do
        end if
        if (my_rank == (i - 1)*num_walkers) then
          wl_logdos = wl_logdos/num_walkers
          do j = 1, num_walkers - 1
            call MPI_Send(wl_logdos, wl_setup_internal%bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers + j, i, MPI_COMM_WORLD, &
                          ierror)
            !print*, my_rank, "send", (i - 1)*num_walkers + j
          end do
        else
          call MPI_Recv(wl_logdos_buffer, wl_setup_internal%bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers, i, MPI_COMM_WORLD, &
                        MPI_STATUS_IGNORE, ierror)
          wl_logdos = wl_logdos_buffer
          !print*, my_rank, "recv", (i - 1)*num_walkers
        end if
      end if
    end do
  end subroutine dos_average
  
  !> @brief   Routine that combines DoS across windows
  !>
  !> @details Routine that "stitches" together DoS across windows.
  !>          Finds point within overlap region where thermodynamic beta is 
  !>          most similar, scales DoS and combines.
  !>          
  !> @param   wl_logdos Double 1D array containing density of states
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine dos_combine(wl_logdos)
    real(real64), intent(inout) :: wl_logdos(:)
    ! Internal
    integer :: i, j, beta_index
    real(real64) :: beta_original, beta_merge, beta_diff, scale_factor
    integer :: mpi_start_idx, mpi_end_idx

    if (my_rank == 0) then
      wl_logdos_combine = wl_logdos
    end if

    do i = 2, wl_setup_internal%num_windows
      if (my_rank == window_rank_index(i,1)) then
        call MPI_Send(wl_logdos, wl_setup_internal%bins, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierror)
        call MPI_Send(window_indices(mpi_index,1), 1, MPI_INT, 0, 1, MPI_COMM_WORLD, ierror)
        call MPI_Send(window_indices(mpi_index,2), 1, MPI_INT, 0, 2, MPI_COMM_WORLD, ierror)
      end if
      if (my_rank == 0) then
        call MPI_Recv(wl_logdos_buffer, wl_setup_internal%bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers, 0, MPI_COMM_WORLD, &
          MPI_STATUS_IGNORE, ierror)
        call MPI_Recv(mpi_start_idx, 1, MPI_INT, (i - 1)*num_walkers, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        call MPI_Recv(mpi_end_idx, 1, MPI_INT, (i - 1)*num_walkers, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        scale_factor = 0.0_real64
        beta_diff = HUGE(beta_diff)
        do j = 0,  window_indices(i-1, 2) - window_indices(i, 1) - 1
          beta_original = wl_logdos_combine(mpi_start_idx + j + 1) - wl_logdos_combine(mpi_start_idx + j)
          beta_merge = wl_logdos_buffer(mpi_start_idx + j + 1) - wl_logdos_buffer(mpi_start_idx + j)
          if (ABS(beta_original - beta_merge) < beta_diff) then
            beta_diff = ABS(beta_original - beta_merge)
            beta_index = mpi_start_idx + j
          end if
        end do

        do j = beta_index, mpi_end_idx
          wl_logdos_combine(j) = wl_logdos_buffer(j) + wl_logdos_combine(beta_index) - wl_logdos_buffer(beta_index)
        end do
      end if
    end do

    if (my_rank == 0) then
      wl_logdos = wl_logdos_combine
    end if
    call MPI_BCAST(wl_logdos, wl_setup_internal%bins, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine dos_combine

  !> @brief   Performs window size optimisation
  !>
  !> @details Routine that adjusts the size of each energy window based on time taken
  !>          for window to first converge. Uses a diffusion coefficient approach with
  !>          a dampening factor that reduces the window size adjustment with each iteration.
  !>          Window sizes are adjusted based on the diffusion coefficients normalised to sum to
  !>          1. Each window is assigned a percentage of the total energy bins based on this
  !>          normalised diffusion coefficient.
  !>          
  !> @param   iter Integer for the current Wang Landau interation
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine mpi_window_optimise(iter)
    integer, intent(in) :: iter

    ! Internal
    integer :: i, j, ierror, min_bins
    real(real64) :: diffusion(wl_setup_internal%num_windows), diffusion_merge(wl_setup_internal%num_windows), factor, scaling
    integer :: bins(wl_setup_internal%num_windows)
    integer :: bins_max_indices(wl_setup_internal%num_windows)
    logical :: mk(wl_setup_internal%num_windows)

    ! Perform window size adjustment then broadcast
    if (wl_setup_internal%num_windows > 1) then
      if (my_rank == 0) then
        factor = 1.0_real64
        scaling = 0.8_real64
        factor = factor*(scaling**(iter))
        diffusion = (REAL((window_intervals(:,2) - window_intervals(:,1) + 1))) &
        /(rank_time_buffer(:,1))
        diffusion_merge = diffusion_prev/sum(diffusion_prev)*(1.0_real64-factor) + factor*diffusion/sum(diffusion)
        diffusion_prev = diffusion_merge
        
        bins = NINT(REAL(wl_setup_internal%bins)*diffusion_merge/SUM(diffusion_merge))
        
        ! Set all bins less than min_bins to min_bins
        min_bins = 1
        do i = 1, wl_setup_internal%num_windows
          if (bins(i) < min_bins) then
              bins(i) = min_bins
          end if
        end do

        mk = .True.
        do i = 1, wl_setup_internal%num_windows
          bins_max_indices(i) = MAXLOC(bins,dim=1,mask=mk)
          mk(MAXLOC(bins,mk)) = .FALSE.
        end do

        i = 0
        j = 0
        do while(SUM(bins) /= wl_setup_internal%bins)
          i = i + 1
          do j = 1, MIN(i, wl_setup_internal%num_windows)
            if (SUM(bins) > wl_setup_internal%bins .and. bins(bins_max_indices(j)) > 2) then
              bins(bins_max_indices(j)) = bins(bins_max_indices(j)) - 1
            end if
            if (SUM(bins) < wl_setup_internal%bins .and. bins(bins_max_indices(j)) > 2) then
              bins(bins_max_indices(j)) = bins(bins_max_indices(j)) + 1
            end if
          end do
        end do

        window_intervals(1, 2) = bins(1)
        do i=2, wl_setup_internal%num_windows
          window_intervals(i, 1) = window_intervals(i-1, 2) + 1
          window_intervals(i, 2) = window_intervals(i, 1) + bins(i) - 1
        end do
        window_intervals(wl_setup_internal%num_windows, 1) = window_intervals(wl_setup_internal%num_windows-1, 2) + 1
        window_intervals(wl_setup_internal%num_windows, 2) = wl_setup_internal%bins
      end if

      call MPI_BCAST(window_intervals, wl_setup_internal%num_windows*2, MPI_INT, 0, MPI_COMM_WORLD, ierror)

      ! Populate MPI arrays and indlude MPI window overlap
      call mpi_arrays(window_intervals, window_indices, mpi_bin_edges, mpi_wl_hist, mpi_bins)

    end if
    mpi_start_idx = window_indices(mpi_index, 1)
    mpi_end_idx = window_indices(mpi_index, 2)
    mpi_bins = mpi_end_idx - mpi_start_idx + 1
  end subroutine mpi_window_optimise

  !> @brief   Routine that creates MPI arrays
  !>
  !> @details This routine deallocates mpi_bin_edges and mpi_wl_hist and then
  !>          re-calculates the window energy bin indices before reallocating
  !>          the MPI arrays and filling them with appropriate values.
  !>          To be used in conjuction with mpi_window_optimise.
  !>          
  !> @param   window_intervals Integer 2D array containing energy bin indices for
  !>          each window
  !> @param   window_indices Integer 2D array containing the energy bin indices for
  !>          each MPI process including ovelap
  !> @param   mpi_bin_edges Double 1D array containing energy bin edges
  !>          for working MPI process
  !> @param   mpi_wl_hist Double 1D array containing energy bin visit information
  !>          for working MPI process
  !> @param   mpi_bins Integer for number of bins assigned to each MPI process
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine mpi_arrays(window_intervals, window_indices, mpi_bin_edges, mpi_wl_hist, mpi_bins)
    integer, dimension(:,:), intent(in) :: window_intervals
    integer, dimension(:,:), intent(inout) :: window_indices
    real(real64), dimension(:), allocatable, intent(inout) :: mpi_bin_edges, mpi_wl_hist
    integer, intent(out) :: mpi_bins

    ! Loop indices
    integer :: i, j

    window_indices(1, 1) = window_intervals(1,1)
    window_indices(1,2) = INT(window_intervals(1,2) + wl_setup_internal%bin_overlap)
    do i = 2, wl_setup_internal%num_windows-1
      window_indices(i, 1) = MAX(INT(window_intervals(i,1) - wl_setup_internal%bin_overlap), 1)
      window_indices(i, 2) = MIN(INT(window_intervals(i,2) + wl_setup_internal%bin_overlap), wl_setup_internal%bins)
    end do
    window_indices(wl_setup_internal%num_windows, 1) = MAX(INT(window_intervals(wl_setup_internal%num_windows,1) &
                                    - wl_setup_internal%bin_overlap), 1)
    window_indices(wl_setup_internal%num_windows,2) = window_intervals(wl_setup_internal%num_windows,2)

    mpi_index = my_rank/num_walkers + 1
    mpi_bins = window_indices(mpi_index,2) - window_indices(mpi_index,1) + 1

    if (allocated(mpi_bin_edges)) deallocate(mpi_bin_edges)
    if (allocated(mpi_wl_hist)) deallocate(mpi_wl_hist)
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

  !> @brief   Replica exchange
  !>
  !> @details Routine that performs a replica exchange if there are two process
  !>          from different windows within the same overlap region.
  !>          
  !> @param   config Short 4D array that stores lattice configuration 
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine replica_exchange(config)
  integer(int16), dimension(:,:,:,:) :: config
  ! Declare local variables
  integer, dimension(num_walkers, 2) :: overlap_lower, overlap_upper
  integer, dimension(num_walkers, 2) :: overlap_exchange
  integer :: i, j, k, exchange_index, ierror
  integer :: exchange_count, ibin, jbin
  logical :: exchange_match, accept
  integer :: overlap_loc, request
  integer :: overlap_mpi(mpi_processes, 2), overlap_mpi_buffer(mpi_processes, 2)
  real(real64) :: delta_e, e_swapped, e_unswapped
   
  ! Perform binning and initialize overlap_mpi
  e_unswapped = setup_internal%full_energy(config)
  ibin = bin_index(e_unswapped, bin_edges, wl_setup_internal%bins)
  jbin = ibin
  accept = .False.
  overlap_mpi = 0
  overlap_exchange = -1

  if (ibin > mpi_start_idx - 1 .and. &
  ibin < mpi_start_idx + wl_setup_internal%bin_overlap .and. mpi_index > 1) then
    overlap_loc = mpi_index - 1
  elseif (ibin > mpi_end_idx - wl_setup_internal%bin_overlap .and. &
   ibin < mpi_end_idx + 1 .and. mpi_index < wl_setup_internal%num_windows) then
    overlap_loc = mpi_index
  else
    overlap_loc = 0
  end if
 
  overlap_mpi(my_rank+1, 1) = my_rank
  overlap_mpi(my_rank+1, 2) = overlap_loc
 
  ! Reduce to find max overlap
  call MPI_REDUCE(overlap_mpi, overlap_mpi_buffer, mpi_processes*2, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ierror)
  overlap_mpi = overlap_mpi_buffer

  ! Exchange loop
  do i=1, wl_setup_internal%num_windows-1
    overlap_exchange = -1
    if (my_rank == 0) then
      overlap_lower = overlap_mpi((i-1)*num_walkers+1:i*num_walkers, :)
      overlap_upper = overlap_mpi((i)*num_walkers+1:(i+1)*num_walkers, :)

      exchange_index = 0
      exchange_count = 1
 
      call shuffle_rows(overlap_lower, num_walkers)
      call shuffle_rows(overlap_upper, num_walkers)
 
      do j = 1, num_walkers
        if (overlap_lower(j, 2) == 0) cycle  ! Skip rows with 0
        do k = 1, num_walkers
          if (overlap_upper(k, 2) == 0) cycle  ! Skip rows with 0
          if (overlap_lower(j, 2) == overlap_upper(k, 2)) then
            ! Matching rows found
            exchange_count = exchange_count + 1
            ! Store the matching row indices (first columns only)
            exchange_index = exchange_index + 1
            overlap_exchange(exchange_index, 1) = overlap_lower(j, 1)
            overlap_exchange(exchange_index, 2) = overlap_upper(k, 1)

            ! Set matching rows to 0
            overlap_lower(j, :) = 0
            overlap_upper(k, :) = 0
          end if
        end do
      end do
    end if

    ! Broadcast exchange data
    call MPI_BCAST(overlap_exchange, num_walkers*2, MPI_INT, 0, MPI_COMM_WORLD, ierror)

    ! MPI SEND RECV calls for replica exchange
    do j=1, COUNT(overlap_exchange(:,1) > -1)
      call comms_wait()
      if (my_rank == overlap_exchange(j,1)) then
        call MPI_RECV(e_swapped, 1, MPI_DOUBLE_PRECISION, overlap_exchange(j,2), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        jbin = bin_index(e_swapped, bin_edges, wl_setup_internal%bins)
        if (genrand() .lt. exp((wl_logdos(ibin) - wl_logdos(jbin)))) then
          accept = .True.
          call MPI_SEND(accept, 1, MPI_INT, overlap_exchange(j,2), 1, MPI_COMM_WORLD, ierror)
          call MPI_ISEND(config, SIZE(config), MPI_SHORT, overlap_exchange(j,2), 2, MPI_COMM_WORLD, request, ierror)
          call MPI_RECV(config, SIZE(config), MPI_SHORT, overlap_exchange(j,2), 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        else
          call MPI_SEND(accept, 1, MPI_INT, overlap_exchange(j,2), 1, MPI_COMM_WORLD, ierror)
        end if
      elseif (my_rank == overlap_exchange(j,2)) then
        call MPI_SEND(e_unswapped, 1, MPI_DOUBLE_PRECISION, overlap_exchange(j,1), 0, MPI_COMM_WORLD, ierror)
        call MPI_RECV(accept, 1, MPI_INT, overlap_exchange(j,1), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        if (accept) then
          call MPI_ISEND(config, SIZE(config), MPI_SHORT, overlap_exchange(j,1), 2, MPI_COMM_WORLD, request, ierror)
          call MPI_RECV(config, SIZE(config), MPI_SHORT, overlap_exchange(j,1), 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        end if
      end if
    end do
   end do
 end subroutine replica_exchange
  
  !> @brief   Routine that shuffles rows in integer array
  !>     
  !> @param   array Integer input array to have rows shuffles
  !> @param   num_walkers Integer for number of walkers within window
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  ! Subroutine to shuffle the rows of the array
  subroutine shuffle_rows(array, num_walkers)
    integer, dimension(num_walkers, 2), intent(inout) :: array
    integer, intent(in) :: num_walkers
    integer :: i, rand_index, temp(2)
    ! Shuffle the rows in the array randomly
    do i = num_walkers, 2, -1
        rand_index = MIN(1+int(genrand()*i), num_walkers)  ! Generate random index between 1 and i
        ! Swap rows i and rand_index
        temp = array(i, :)
        array(i, :) = array(rand_index, :)
        array(rand_index, :) = temp
    end do
  end subroutine

  !> @brief   Purge pending communication
  !>
  !> @details Routine that recieves any pending communications into a buffer array
  !>          that is allocated, filled then deallocated. Essentially deletes pending
  !>          communications.
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024 
  subroutine comms_purge()
    integer :: mpi_counter, status(MPI_STATUS_SIZE), ierror
    logical :: flag
    character, allocatable :: discard(:)
    do while(.true.)
      CALL MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, flag, status, ierror)
      if (.not. flag) exit  ! No more messages
      ! get the message size in bytes
      CALL MPI_GET_COUNT(status, MPI_BYTE, mpi_counter, ierror)
      ! allocate a dummy buffer
      allocate(discard(mpi_counter))
      ! receive and discard the message
      CALL MPI_RECV(discard, mpi_counter, MPI_BYTE, status(MPI_SOURCE), status(MPI_TAG), MPI_COMM_WORLD, status, ierror)
      deallocate(discard)  ! free memory
    end do
  end subroutine comms_purge

end module wang_landau
