!----------------------------------------------------------------------!
! main.f90                                                             !
!                                                                      !
! Main bontewarlo program.                                             !
!                                                                      !
! C. D. Woodgate,  Warwick                                        2024 !
!----------------------------------------------------------------------!
program main
  
  use initialise
  use comms
  use shared_data
  use metropolis
  use nested_sampling
  use tmmc
  use wang_landau
  use energy_spectrum
  use config_output
  use io
  use kinds
  use c_functions
  use write_netcdf
  use write_xyz
  use write_diagnostics
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

  ! Energy Spectrum parameters type
  type(es_params) :: es_setup

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
  call initialise_pnrg(setup%seedtime)

  ! Initialise some function pointers
  call initialise_function_pointers(setup)

  ! Make directories for data
  call make_data_directories(my_rank)

  ! Initialise some local arrays
  call initialise_local_arrays(setup)

  !---------------!
  ! Main Routines !
  !---------------!
  if (setup%mode == 301) then

    ! Metropolis with Kawasaki dynamics
    call metropolis_simulated_annealing(setup, metropolis_setup, my_rank)

  else if (setup%mode == 302) then

    ! Draw decorrelated samples
    call metropolis_decorrelated_samples(setup, metropolis_setup, my_rank)

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
    call nested_sampling_main(setup, ns_setup, my_rank)

  else if (setup%mode == 304) then

    ! Tmmc algorithm
    call read_tmmc_file("tmmc_input.txt", tmmc_setup, my_rank)
    call tmmc_main(setup, tmmc_setup, my_rank)

  else if (setup%mode == 305) then

    ! Wang Landau algorithm
    call read_wl_file("wl_input.txt", wl_setup, my_rank)
    call wl_main(setup, wl_setup, my_rank)

  else if (setup%mode == 306) then

    ! Energy spectrum algorithm
    call read_es_file("es_input.txt", es_setup, my_rank)
    call es_main(setup, es_setup, my_rank)

  else if (setup%mode == 307) then

    ! Configuration at temperature output
    call co_main(setup, my_rank)

  else

   print*, ' Unrecognised mode', setup%mode

  end if

  ! Clean up
  call clean_up_interaction()

  ! Clean up
  call local_clean_up(setup)

  ! Print software info to the screen
  if(my_rank == 0) call write_info('f')

  ! Finalise MPI
  call comms_finalise()

end program main
