module comms

  use mpi
  use kinds
  use mpi_shared_data

  implicit none
  save

  !> total number of processes
  integer :: p
  !> rank of my processor
  integer :: my_rank
  !> start and end times for mpi
  real(real64) :: t1, t2
  !> mpi status_info
  integer, dimension(mpi_status_size) :: status_info
  !> error variables
  integer :: ierr
  !> our communicator
  integer :: cart_comm
  !> this processor coordinates
  integer, dimension(3) :: my_rank_coords
  !> neighbouring ranks
  integer, dimension(6) :: my_rank_neighbours
  integer :: east, west, north, south, up, down

  contains

  !--------------------------------------------------------
  !> @author
  !> c. woodgate
  !
  ! description:
  !> comms_initialise sets up the mpi communicator
  !
  ! revision history:
  ! todo_dd_mm_yyyy -
  !
  !> @param[in] error
  !--------------------------------------------------------
  subroutine comms_initialise()

    !> initialise mpi
    call mpi_init(ierr)

    !> set up the size and rank of the communicator
    call mpi_comm_rank(mpi_comm_world,my_rank,ierr)
    call mpi_comm_size(mpi_comm_world,p,ierr)

    !> start the clock
    t1 = mpi_wtime()

  end subroutine comms_initialise

  subroutine comms_wait()

  call mpi_barrier(mpi_comm_world, ierr)

  end subroutine comms_wait

  subroutine comms_reduce_results(setup)

    type(run_params), intent(in) :: setup

    call mpi_reduce(energies_of_T, av_energies_of_T, setup%T_steps,         &
                    MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, ierr)

    av_energies_of_T = av_energies_of_T/real(p)
  
    call mpi_reduce(C_of_T, av_C_of_T, setup%T_steps,                       &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_comm_world, ierr)
  
    av_C_of_T = av_C_of_T/real(p)
  
    call mpi_reduce(acceptance_of_T, av_acceptance_of_T, setup%T_steps,     &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_comm_world, ierr)
  
    av_acceptance_of_T = av_acceptance_of_T/real(p)
  
    call mpi_reduce(rho_of_T, av_rho_of_T,                                  &
                    setup%T_steps*(setup%n_species**2)*setup%radial_d,      &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_comm_world, ierr)
  
    av_rho_of_T = av_rho_of_T/real(p)
  

  end subroutine comms_reduce_results

  !--------------------------------------------------------
  !> @author
  !> c. woodgate
  !
  ! description:
  !> comms_finalize cleans up the mpi communicator
  !
  ! revision history:
  ! todo_dd_mm_yyyy -
  !
  !> @param[in] error
  !--------------------------------------------------------
  subroutine comms_finalise()

    !> stop the clock
    t2 = mpi_wtime()

    !> print time taken
    if(my_rank .eq. 0) then
      print '(a,f0.2,a)', "total time elapsed for mpi run is ", t2-t1, " seconds."
    end if

    !> clean up mpi
    call mpi_finalize(ierr)

  end subroutine comms_finalise

end module comms
