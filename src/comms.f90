!> @file    comms.f90
!>
!> @brief   Routines associated with calls to the MPI library
!>
!> @author  C. D. Woodgate
!> @date    2021-2025
module comms

  use mpi
  use kinds
  use shared_data
  use derived_types

  implicit none

  private

  public :: my_rank, p

  public :: comms_initialise, comms_wait,                              &
            comms_reduce_metropolis_results, comms_finalise, comms_purge

  save

  ! total number of processes
  integer :: p

  ! rank of my processor
  integer :: my_rank

  ! start and end times for mpi
  real(real64) :: t1, t2

  ! mpi status_info
  integer, dimension(mpi_status_size) :: status_info

  ! error variables
  integer :: ierr

  ! has the message been recieved
  logical :: flag

  ! size of message
  integer :: mpi_counter

  contains

  !> @brief   Subroutine to initialise MPI
  !>
  !> @author  C. D. Woodgate
  !> @date    2021-2023
  !>
  !> @return None
  subroutine comms_initialise()

    ! initialise mpi
    call mpi_init(ierr)

    ! set up the size and rank of the communicator
    call mpi_comm_rank(mpi_comm_world,my_rank,ierr)
    call mpi_comm_size(mpi_comm_world,p,ierr)

    ! start the clock
    t1 = mpi_wtime()

  end subroutine comms_initialise

  !> @brief   Subroutine to call MPI_BARRIER---all processes wait
  !>
  !> @author  C. D. Woodgate
  !> @date    2019-2023
  !>
  !> @return None
  subroutine comms_wait()

    ! Call MPI_BARRIER
    call mpi_barrier(mpi_comm_world, ierr)

  end subroutine comms_wait

  !> @brief   Reduces (averages) results of a parallel ensemble of 
  !>          Metropolis Monte Carlo simulations.
  !>
  !> @author  C. D. Woodgate
  !> @date    2021-2023
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  setup Derived type containing Metropolis parameters
  !>
  !> @return None
  subroutine comms_reduce_metropolis_results(setup, metropolis)

    ! Input contains information about simulation
    type(run_params), intent(in) :: setup
    ! Input contains information about Metropolis
    type(metropolis_params), intent(in) :: metropolis

    ! Bring all simulation energy arrays, <E>(T), to rank 0
    ! and sum them.
    call mpi_reduce(energies_of_T, av_energies_of_T,                   &
                    metropolis%T_steps, MPI_DOUBLE, MPI_SUM, 0,        &
                    mpi_comm_world, ierr)

    ! Divide by the number of simulations to get the average
    av_energies_of_T = av_energies_of_T/real(p)
  
    ! Do the same for the heat capacity data
    call mpi_reduce(C_of_T, av_C_of_T,metropolis%T_steps,              &
                 MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_comm_world, ierr)
    av_C_of_T = av_C_of_T/real(p)
  
    ! Do the same with the acceptance rates
    call mpi_reduce(acceptance_of_T, av_acceptance_of_T,               &
                    metropolis%T_steps, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, mpi_comm_world, ierr)
    av_acceptance_of_T = av_acceptance_of_T/real(p)
  
    ! Do the same with the radial densities
    call mpi_reduce(rho_of_T, av_rho_of_T,                             &
               metropolis%T_steps*(setup%n_species**2)*setup%wc_range, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_comm_world,  &
                    ierr)
    av_rho_of_T = av_rho_of_T/real(p)

  end subroutine comms_reduce_metropolis_results

  !> @brief   Subroutine to finalise MPI
  !>
  !> @details The block which is currently commented out can be used to
  !>          display the time taken.
  !>
  !> @author  C. D. Woodgate
  !> @date    2021-2023
  !>
  !> @return None
  subroutine comms_finalise()

    ! stop the clock
    t2 = mpi_wtime()

!    ! print time taken
!    if(my_rank .eq. 0) then
!      print '(a,f0.2,a)', "Total time elapsed for mpi run is ", t2-t1, &
!                          " seconds."
!    end if

    ! clean up mpi
    call mpi_finalize(ierr)

  end subroutine comms_finalise

  !> @brief   Purge pending communication
  !>
  !> @details Routine that recieves any pending communications into a buffer array
  !>          that is allocated, filled then deallocated. Essentially deletes pending
  !>          communications.
  !> 
  !> @return  None
  !>
  !> @author  H. J. Naguszewski
  !> @date    2025 
  subroutine comms_purge()
    character, allocatable :: discard(:)
    do while(.true.)
      CALL MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, flag, status_info, ierr)
      if (.not. flag) exit  ! No more messages
      ! get the message size in bytes
      CALL MPI_GET_COUNT(status_info, MPI_BYTE, mpi_counter, ierr)
      ! allocate a dummy buffer
      allocate(discard(mpi_counter))
      ! receive and discard the message
      CALL MPI_RECV(discard, mpi_counter, MPI_BYTE, status_info(MPI_SOURCE), status_info(MPI_TAG), &
      MPI_COMM_WORLD, status_info, ierr)
      deallocate(discard)  ! free memory
    end do
  end subroutine comms_purge

end module comms
