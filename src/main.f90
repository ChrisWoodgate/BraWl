!> @file    main.f90
!>
!> @brief   Main BraWl program
!>
!> @author  H. J. Naguszewski
!> @author  L. B. Partay
!> @author  C. D. Woodgate
!>
!> @date    2019-2025
program main
  
  use initialise
  use comms
  use shared_data
  use metropolis
  use nested_sampling
  use tmmc
  use wang_landau
  use io
  use kinds
  use c_functions
  use write_netcdf
  use write_xyz
  use metropolis_output
  use display

  implicit none

  ! Runtime parameters type
  type(run_params) :: setup

  ! Metropolis parameters type
  type(metropolis_params) :: metropolis_setup

  ! Nested Sampling parameters type
  type(ns_params) :: ns_setup

  ! Tmmc parameters type
  type(tmmc_params) :: tmmc_setup

  ! Wang Landau parameters type
  type(wl_params) :: wl_setup

  !--------------------------------------------------------------------!
  !                          Initial Setup                             !
  !--------------------------------------------------------------------!

  ! Start MPI
  call comms_initialise()

  ! Print software info to the screen
  if(my_rank == 0) call write_info('s')

  ! Parse inputs
  call parse_inputs(setup, my_rank)

  ! Allocate space for atom-atom interaction parameters
  call initialise_interaction(setup)

  ! Read in atom-atom interaction
  call read_exchange(setup, my_rank)

  ! Initialise PNRG
  ! static_seed is true if we would like to use a fixed seed and false
  ! otherwise
  call initialise_prng(setup%static_seed)

  ! Initialise some function pointers
  call initialise_function_pointers(setup)

  ! Initialise some local arrays
  call initialise_local_arrays(setup)

  !--------------------------------------------------------------------!
  !                           Main Routines                            !
  !--------------------------------------------------------------------!

  !---------------------------------!
  ! Metropolis-Hastings Monte Carlo !
  !---------------------------------!
  if (setup%mode == 301) then

    ! Make the relevant directories
    if(my_rank == 0) call execute_command_line('mkdir -p grids')
    if(my_rank == 0) call execute_command_line('mkdir -p diagnostics')
    if(my_rank == 0) call execute_command_line('mkdir -p radial_densities')
    if(my_rank == 0) call execute_command_line('mkdir -p trajectories')
    if(my_rank == 0) call execute_command_line('mkdir -p alro')
    if(my_rank == 0) call execute_command_line('mkdir -p energies')

    ! Run Metropolis with Kawasaki dynamics
    call metropolis_simulated_annealing(setup, metropolis_setup, my_rank)

  else if (setup%mode == 302) then

    ! Make the relevant directories
    if(my_rank == 0) call execute_command_line('mkdir -p grids')
    if(my_rank == 0) call execute_command_line('mkdir -p diagnostics')
    if(my_rank == 0) call execute_command_line('mkdir -p radial_densities')

    ! Draw decorrelated samples
    call metropolis_decorrelated_samples(setup, metropolis_setup, my_rank)

  !---------------------------!
  ! Nested sampling algorithm !
  !---------------------------!
  else if (setup%mode == 303) then

    ! The variable 'p' is defined in comms.f90 and is the number of MPI
    ! processes. The nested sampling algorithm is only implemented in
    ! serial at present, so we terminate the programme if the user has
    ! tried to execute in parallel
    if (p .gt. 1) then
      call comms_finalise()
      stop ' Nested sampling is only implemented in serial'
    end if

    ! Nested Sampling algorithm
    call nested_sampling_main(setup, ns_setup)

  !-------------------------------!
  ! Transition matrix Monte Carlo !
  !-------------------------------!
  else if (setup%mode == 304) then

    ! Read TMMC input file
    call read_tmmc_file("tmmc_input.inp", tmmc_setup, my_rank)

    ! Run main TMMC routine
    ! Note: this is currently under development. The program will exit
    !       if you try to use this routine.
    call tmmc_main(setup, tmmc_setup, my_rank)

  !----------------------!
  ! Wang-Landau sampling !
  !----------------------!
  else if (setup%mode == 305) then

    ! Make the relevant directories
    if(my_rank == 0) call execute_command_line('mkdir -p radial_densities')

    ! Read the Wang-Landau input file
    call read_wl_file("wl_input.inp", wl_setup, my_rank)

    ! Run Wang Landau sampling
    call wl_main(setup, wl_setup)

  !-------------------------------------!
  ! Case of unrecognised mode requested !
  !-------------------------------------!
  else

   if (my_rank .eq. 1) then
     print*, ' Unrecognised mode', setup%mode
   end if
   call comms_finalise()
   stop 'Exiting as unrecognised mode requested'

  end if

  !--------------------------------------------------------------------!
  !                     Cleanup at end of run                          !
  !--------------------------------------------------------------------!

  ! Clean up
  call clean_up_interaction()

  ! Clean up
  call local_clean_up(setup)

  ! Print software info to the screen
  if(my_rank == 0) call write_info('f')

  ! Finalise MPI
  call comms_finalise()

end program main
