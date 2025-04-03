!----------------------------------------------------------------------!
! write_netcdf.f90                                                     !
!                                                                      !
! Module containing all calls to NetCDF library routines               !
!                                                                      !
! C. D. Woodgate,  Warwick                                        2023 !
!----------------------------------------------------------------------!
module write_netcdf

  use kinds
  use netcdf
  use shared_data

  implicit none

  contains

  !--------------------------------------------------------------------!
  ! Routine to write radial density array calculated once to file      !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2024 !
  !--------------------------------------------------------------------!
  subroutine ncdf_radial_density_writer_once(filename, rho, r, setup)

    integer, parameter :: rho_ndims = 3

    type(run_params), intent(in) :: setup

    ! Data to write to file
    real(real64), dimension(:,:,:), allocatable, intent(in) :: rho
    real(real64), dimension(:), allocatable, intent(in) :: r

    ! Number of dimensions of my grid data
    integer, dimension(rho_ndims) :: rho_sizes, rho_dim_ids
    integer :: r_size, r_dim_id, T_size, T_dim_id, U_size, U_dim_id

    ! Names of my dimensions
    character(len=1), dimension(rho_ndims) :: rho_dims=(/"i", "j", "r"/)
    character(len=3) :: r_dims = "r_i"

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    ! Variables used in writing process
    integer :: file_id, i

    ! Ids for variables
    integer :: rho_id, r_id

    ! Get the sizes of my incoming arrays
    rho_sizes  = shape(rho)
    r_size = size(r)

    ! Create the file
    call check(nf90_create(filename, nf90_clobber, file_id))

    ! Add information about global runtime data
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_1', setup%n_1))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_2', setup%n_2))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_3', setup%n_3))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Number of Species', setup%n_species))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Lattice Type', setup%lattice))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Interaction file', setup%interaction_file))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Concentrations', setup%species_concentrations))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Warren-Cowley Range', &
                            setup%wc_range))

    ! Define the 3D variables and dimensions
    do i = 1, rho_ndims
      call check(nf90_def_dim(file_id, rho_dims(i), &
                              rho_sizes(i), rho_dim_ids(i)))
    end do

    call check(nf90_def_var(file_id, "rho data", NF90_DOUBLE, &
                            rho_dim_ids, rho_id))

    call check(nf90_def_dim(file_id, r_dims, r_size, r_dim_id))

    call check(nf90_def_var(file_id, "r data", NF90_DOUBLE, &
                            r_dim_id, r_id))

    ! Finish defining metadata
    call check(nf90_enddef(file_id))

    ! Dump the variables to file
    call check(nf90_put_var(file_id, rho_id, rho))
    call check(nf90_put_var(file_id, r_id, r))

    ! Close the file
    call check(nf90_close(file_id))

  end subroutine ncdf_radial_density_writer_once

  !--------------------------------------------------------------------!
  ! Routine to write radial densities to file                          !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine ncdf_radial_density_writer(filename, rho, r, T, U_of_T, setup)

    integer, parameter :: rho_ndims = 4

    type(run_params), intent(in) :: setup

    ! Data to write to file
    real(real64), dimension(:,:,:,:), allocatable, intent(in) :: rho
    real(real64), dimension(:), allocatable, intent(in) :: r, T
    real(real64), dimension(:), allocatable, intent(in) :: U_of_T

    ! Number of dimensions of my grid data
    integer, dimension(rho_ndims) :: rho_sizes, rho_dim_ids
    integer :: r_size, r_dim_id, T_size, T_dim_id, U_size, U_dim_id

    ! Names of my dimensions
    character(len=1), dimension(rho_ndims) :: rho_dims=(/"i", "j", "r", "T"/)
    character(len=3) :: r_dims = "r_i"
    character(len=3) :: T_dims = "T_i"
    character(len=3) :: U_dims = "U_i"

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    ! Variables used in writing process
    integer :: file_id, i

    ! Ids for variables
    integer :: rho_id, r_id, T_id, U_id

    ! Get the sizes of my incoming arrays
    rho_sizes  = shape(rho)
    r_size = size(r)
    T_size = size(T)
    U_size = size(U_of_T)

    ! Create the file
    call check(nf90_create(filename, nf90_clobber, file_id))

    ! Add information about global runtime data
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_1', setup%n_1))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_2', setup%n_2))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_3', setup%n_3))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Number of Species', setup%n_species))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Lattice Type', setup%lattice))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Interaction file', setup%interaction_file))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Concentrations', setup%species_concentrations))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Warren-Cowley Range', &
                            setup%wc_range))

    ! Define the 3D variables and dimensions
    do i = 1, rho_ndims
      call check(nf90_def_dim(file_id, rho_dims(i), &
                              rho_sizes(i), rho_dim_ids(i)))
    end do 

    call check(nf90_def_var(file_id, "rho data", NF90_DOUBLE, &
                            rho_dim_ids, rho_id))

    call check(nf90_def_dim(file_id, r_dims, r_size, r_dim_id))

    call check(nf90_def_var(file_id, "r data", NF90_DOUBLE, &
                            r_dim_id, r_id))

    call check(nf90_def_dim(file_id, T_dims, T_size, T_dim_id))

    call check(nf90_def_var(file_id, "T data", NF90_DOUBLE, &
                            T_dim_id, T_id))

    call check(nf90_def_dim(file_id, U_dims, U_size, U_dim_id))

    call check(nf90_def_var(file_id, "U data", NF90_DOUBLE, &
                            U_dim_id, U_id))

    ! Finish defining metadata
    call check(nf90_enddef(file_id))

    ! Dump the variables to file
    call check(nf90_put_var(file_id, rho_id, rho))
    call check(nf90_put_var(file_id, r_id, r))
    call check(nf90_put_var(file_id, t_id, t))
    call check(nf90_put_var(file_id, U_id, U_of_T))

    ! Close the file
    call check(nf90_close(file_id))

  end subroutine ncdf_radial_density_writer

  !--------------------------------------------------------------------!
  ! Routine to write radial densities across energy to file            !
  !                                                                    !
  ! H. J. Naguszewski,  Warwick                                   2023 !
  !--------------------------------------------------------------------!
  subroutine ncdf_radial_density_writer_across_energy(filename, rho, r, U, setup)

    integer, parameter :: rho_ndims = 4

    type(run_params), intent(in) :: setup

    ! Data to write to file
    real(real64), dimension(:,:,:,:), allocatable, intent(in) :: rho
    real(real64), dimension(:), allocatable, intent(in) :: r
    real(real64), dimension(:), allocatable, intent(in) :: U

    ! Number of dimensions of my grid data
    integer, dimension(rho_ndims) :: rho_sizes, rho_dim_ids
    integer :: r_size, r_dim_id, U_size, U_dim_id

    ! Names of my dimensions
    character(len=1), dimension(rho_ndims) :: rho_dims=(/"i", "j", "r", "U"/)
    character(len=3) :: r_dims = "r_i"
    character(len=3) :: U_dims = "U_i"

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    ! Variables used in writing process
    integer :: file_id, i

    ! Ids for variables
    integer :: rho_id, r_id, U_id

    ! Get the sizes of my incoming arrays
    rho_sizes  = shape(rho)
    r_size = size(r)
    U_size = size(U)

    ! Create the file
    call check(nf90_create(filename, nf90_clobber, file_id))

    ! Add information about global runtime data
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_1', setup%n_1))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_2', setup%n_2))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_3', setup%n_3))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Number of Species', setup%n_species))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Lattice Type', setup%lattice))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Interaction file', setup%interaction_file))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Concentrations', setup%species_concentrations))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Warren-Cowley Range', &
                            setup%wc_range))

    ! Define the 3D variables and dimensions
    do i = 1, rho_ndims
      call check(nf90_def_dim(file_id, rho_dims(i), &
                              rho_sizes(i), rho_dim_ids(i)))
    end do 

    call check(nf90_def_var(file_id, "rho data", NF90_DOUBLE, &
                            rho_dim_ids, rho_id))

    call check(nf90_def_dim(file_id, r_dims, r_size, r_dim_id))

    call check(nf90_def_var(file_id, "r data", NF90_DOUBLE, &
                            r_dim_id, r_id))

    call check(nf90_def_dim(file_id, U_dims, U_size, U_dim_id))

    call check(nf90_def_var(file_id, "U data", NF90_DOUBLE, &
                            U_dim_id, U_id))

    ! Finish defining metadata
    call check(nf90_enddef(file_id))

    ! Dump the variables to file
    call check(nf90_put_var(file_id, rho_id, rho))
    call check(nf90_put_var(file_id, r_id, r))
    call check(nf90_put_var(file_id, U_id, U))

    ! Close the file
    call check(nf90_close(file_id))

  end subroutine ncdf_radial_density_writer_across_energy

  !--------------------------------------------------------------------!
  ! Routine to write order parameters to file                          !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine ncdf_order_writer(filename, ierr, order, temperature, setup)

    integer, parameter :: grid_ndims = 5

    type(run_params), intent(in) :: setup

    ! Data to write to file
    real(real64), dimension(:,:,:,:,:), allocatable, intent(in) :: order
    real(real64), dimension(:), allocatable, intent(in) :: temperature

    ! Number of dimensions of my grid data
    integer, dimension(grid_ndims) :: grid_sizes, grid_dim_ids
    integer :: temp_size, temp_dim_id

    ! Names of my dimensions
    character(len=1), dimension(grid_ndims) :: grid_dims=(/"x", "y" , "z", "s", "t"/)
    character(len=4) :: temp_dims = "temp"

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    ! Variables used in writing process
    integer :: ierr, file_id, i

    ! Ids for variables
    integer :: grid_id, temp_id

    ! Get the sizes of my incoming arrays
    grid_sizes  = shape(order)
    temp_size = size(temperature)

    ! Create the file, overwriting if it exists !
    ierr = nf90_create(filename, nf90_clobber, file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ! Add information about global runtime data
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_1', setup%n_1))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_2', setup%n_2))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_3', setup%n_3))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Number of Species', setup%n_species))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Lattice Type', setup%lattice))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Interaction file', setup%interaction_file))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Concentrations', setup%species_concentrations))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Warren-Cowley Range', &
                            setup%wc_range))

    ! Define the 4D variables and dimensions !
    do i = 1, grid_ndims
      ierr = nf90_def_dim(file_id, grid_dims(i), grid_sizes(i), grid_dim_ids(i))
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
    end do 

    ierr = nf90_def_var(file_id, "grid data", NF90_DOUBLE, grid_dim_ids, grid_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    if (allocated(temperature)) then
      ! Define the temperature variable and dimensions !
      ierr = nf90_def_dim(file_id, temp_dims, temp_size, temp_dim_id)
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
  
      ierr = nf90_def_var(file_id, "temperature data", NF90_DOUBLE, temp_dim_id, temp_id)
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
    end if

    ! Finish defining metadata !
    ierr = nf90_enddef(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ! Dump the variables to file !
    ierr = nf90_put_var(file_id, grid_id, order) 
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    if (allocated(temperature)) then
      ierr = nf90_put_var(file_id, temp_id, temperature) 
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
    end if

    ! Close the file !
    ierr = nf90_close(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

  end subroutine ncdf_order_writer

  !--------------------------------------------------------------------!
  ! Routine to write grid state to file                                !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine ncdf_grid_state_writer(filename, ierr, &
                                    state, temperature, setup)

    integer, parameter :: grid_ndims = 4

    type(run_params), intent(in) :: setup

    ! Data to write to file
    integer(int16), dimension(:,:,:,:), allocatable, intent(in) :: state
    real(real64), intent(in) :: temperature

    ! Number of dimensions of my grid data
    integer, dimension(grid_ndims) :: grid_sizes, grid_dim_ids

    ! Names of my dimensions
    character(len=1), dimension(grid_ndims) :: &
                      grid_dims=(/"b", "x", "y" ,"z"/)

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    ! Variables used in writing process
    integer :: ierr, file_id, i

    ! Ids for variables
    integer :: grid_id

    ! Get the sizes of my incoming arrays
    grid_sizes  = shape(state)

    ! Create the file, overwriting if it exists !
    ierr = nf90_create(filename, nf90_clobber, file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !> Add information about global runtime data
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_basis', setup%n_basis))
    ! Add information about global runtime data
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_1', setup%n_1))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_2', setup%n_2))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_3', setup%n_3))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Number of Species', setup%n_species))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'temperature', temperature))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Lattice Type', setup%lattice))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Interaction file', setup%interaction_file))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Concentrations', setup%species_concentrations))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Warren-Cowley Range', &
                            setup%wc_range))

    ! Define the 4D variables and dimensions !
    do i = 1, grid_ndims
      ierr = nf90_def_dim(file_id, grid_dims(i), grid_sizes(i), grid_dim_ids(i))
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
    end do 

    ierr = nf90_def_var(file_id, "grid data", NF90_SHORT, grid_dim_ids, grid_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ! Finish defining metadata !
    ierr = nf90_enddef(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ! Dump the variables to file !
    ierr = nf90_put_var(file_id, grid_id, state) 
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ! Close the file !
    ierr = nf90_close(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

  end subroutine ncdf_grid_state_writer

  !---------------------------------------------!
  ! Routine to write order parameters to a file !
  !---------------------------------------------!
  subroutine ncdf_grid_states_writer(filename, ierr, states, temperature, setup)

    integer, parameter :: grid_ndims = 4

    type(run_params), intent(in) :: setup

    ! Data to write to file
    integer(int16), dimension(:,:,:,:), allocatable, intent(in) :: states
    real(real64), dimension(:), allocatable, intent(in) :: temperature

    ! Number of dimensions of my grid data
    integer, dimension(grid_ndims) :: grid_sizes, grid_dim_ids
    integer :: temp_size, temp_dim_id

    ! Names of my dimensions
    character(len=1), dimension(grid_ndims) :: grid_dims=(/"x", "y" , "z", "t"/)
    character(len=4) :: temp_dims = "temp"

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    ! Variables used in writing process
    integer :: ierr, file_id, i

    ! Ids for variables
    integer :: grid_id, temp_id

    ! Get the sizes of my incoming arrays
    grid_sizes  = shape(states)
    temp_size = size(temperature)

    ! Create the file, overwriting if it exists !
    ierr = nf90_create(filename, nf90_clobber, file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ! Add information about global runtime data
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_1', setup%n_1))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_2', setup%n_2))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'N_3', setup%n_3))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Number of Species', setup%n_species))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Lattice Type', setup%lattice))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Interaction file', setup%interaction_file))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Concentrations', setup%species_concentrations))
    call check(nf90_put_att(file_id, NF90_GLOBAL, &
                            'Warren-Cowley Range', &
                            setup%wc_range))

    ! Define the 4D variables and dimensions !
    do i = 1, grid_ndims
      ierr = nf90_def_dim(file_id, grid_dims(i), grid_sizes(i), grid_dim_ids(i))
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
    end do 

    ierr = nf90_def_var(file_id, "grid data", NF90_SHORT, grid_dim_ids, grid_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    if (allocated(temperature)) then
      ! Define the temperature variable and dimensions !
      ierr = nf90_def_dim(file_id, temp_dims, temp_size, temp_dim_id)
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
  
      ierr = nf90_def_var(file_id, "temperature data", NF90_DOUBLE, temp_dim_id, temp_id)
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
    end if

    ! Finish defining metadata !
    ierr = nf90_enddef(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ! Dump the variables to file !
    ierr = nf90_put_var(file_id, grid_id, states) 
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    if (allocated(temperature)) then
      ierr = nf90_put_var(file_id, temp_id, temperature) 
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
    end if

    ! Close the file !
    ierr = nf90_close(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

  end subroutine ncdf_grid_states_writer

  subroutine ncdf_writer_1d(filename, ierr, grid_data)

    integer, parameter :: grid_ndims = 1

    ! Data to write to file
    real(real64), dimension(:), allocatable, intent(in) :: grid_data
    ! Number of dimensions of my rho data
    integer, dimension(grid_ndims) :: grid_sizes, grid_dim_ids

    ! Names of my dimensions
    character(len=1), dimension(grid_ndims) :: grid_dims=(/"x"/)

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    ! Variables used in writing process
    integer :: ierr, file_id, i

    ! Ids for variables
    integer :: grid_id

    ! Get the sizes of my incoming arrays
    grid_sizes  = shape(grid_data)

    !-------------------------------------------!
    ! Create the file, overwriting if it exists !
    !-------------------------------------------!
    ierr = nf90_create(filename, nf90_clobber, file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !----------------------------------------!
    ! Define the 2D variables and dimensions !
    !----------------------------------------!
    do i = 1, grid_ndims
      ierr = nf90_def_dim(file_id, grid_dims(i), grid_sizes(i), grid_dim_ids(i))
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
    end do 

    ! grid_data
    ierr = nf90_def_var(file_id, "grid data", NF90_DOUBLE, grid_dim_ids, grid_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !--------------------------!
    ! Finish defining metadata !
    !--------------------------!
    ierr = nf90_enddef(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ierr = nf90_put_var(file_id, grid_id, grid_data) 
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !----------------!
    ! Close the file !
    !----------------!
    ierr = nf90_close(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

  end subroutine ncdf_writer_1d

  subroutine ncdf_writer_2d(filename, ierr, grid_data)

    integer, parameter :: grid_ndims = 2

    ! Data to write to file
    real(real64), dimension(:,:), allocatable, intent(in) :: grid_data
    ! Number of dimensions of my rho data
    integer, dimension(grid_ndims) :: grid_sizes, grid_dim_ids

    ! Names of my dimensions
    character(len=1), dimension(grid_ndims) :: grid_dims=(/"x", "y" /)

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    ! Variables used in writing process
    integer :: ierr, file_id, i

    ! Ids for variables
    integer :: grid_id

    ! Get the sizes of my incoming arrays
    grid_sizes  = shape(grid_data)

    !-------------------------------------------!
    ! Create the file, overwriting if it exists !
    !-------------------------------------------!
    ierr = nf90_create(filename, nf90_clobber, file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !----------------------------------------!
    ! Define the 2D variables and dimensions !
    !----------------------------------------!
    do i = 1, grid_ndims
      ierr = nf90_def_dim(file_id, grid_dims(i), grid_sizes(i), grid_dim_ids(i))
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
    end do 

    ! grid_data
    ierr = nf90_def_var(file_id, "grid data", NF90_DOUBLE, grid_dim_ids, grid_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !--------------------------!
    ! Finish defining metadata !
    !--------------------------!
    ierr = nf90_enddef(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ierr = nf90_put_var(file_id, grid_id, grid_data) 
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !----------------!
    ! Close the file !
    !----------------!
    ierr = nf90_close(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

  end subroutine ncdf_writer_2d

  subroutine ncdf_writer_5d(filename, ierr, grid_data)

    integer, parameter :: grid_ndims = 5

    ! Data to write to file
    real(real64), dimension(:,:,:,:,:), allocatable, intent(in) :: &
      grid_data
    ! Number of dimensions of my rho data
    integer, dimension(grid_ndims) :: grid_sizes, grid_dim_ids

    ! Names of my dimensions
    character(len=2), dimension(grid_ndims) :: &
      grid_dims=(/"s1", "s2", "xx", "yy" , "zz"/)

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    ! Variables used in writing process
    integer :: ierr, file_id, i

    ! Ids for variables
    integer :: grid_id

    ! Get the sizes of my incoming arrays
    grid_sizes  = shape(grid_data)

    !-------------------------------------------!
    ! Create the file, overwriting if it exists !
    !-------------------------------------------!
    ierr = nf90_create(filename, nf90_clobber, file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !----------------------------------------!
    ! Define the 3D variables and dimensions !
    !----------------------------------------!
    do i = 1, grid_ndims
      ierr = nf90_def_dim(file_id, grid_dims(i), grid_sizes(i), grid_dim_ids(i))
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
    end do 

    ! grid_data
    ierr = nf90_def_var(file_id, "grid data", NF90_DOUBLE, grid_dim_ids, grid_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !--------------------------!
    ! Finish defining metadata !
    !--------------------------!
    ierr = nf90_enddef(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ierr = nf90_put_var(file_id, grid_id, grid_data) 
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !----------------!
    ! Close the file !
    !----------------!
    ierr = nf90_close(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

  end subroutine ncdf_writer_5d

  subroutine ncdf_writer_3d(filename, ierr, grid_data)

    integer, parameter :: grid_ndims = 3

    ! Data to write to file
    real(real64), dimension(:,:,:), allocatable, intent(in) :: grid_data
    ! Number of dimensions of my rho data
    integer, dimension(grid_ndims) :: grid_sizes, grid_dim_ids

    ! Names of my dimensions
    character(len=1), dimension(grid_ndims) :: grid_dims=(/"x", "y" , "z"/)

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    ! Variables used in writing process
    integer :: ierr, file_id, i

    ! Ids for variables
    integer :: grid_id

    ! Get the sizes of my incoming arrays
    grid_sizes  = shape(grid_data)

    !-------------------------------------------!
    ! Create the file, overwriting if it exists !
    !-------------------------------------------!
    ierr = nf90_create(filename, nf90_clobber, file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !----------------------------------------!
    ! Define the 3D variables and dimensions !
    !----------------------------------------!
    do i = 1, grid_ndims
      ierr = nf90_def_dim(file_id, grid_dims(i), grid_sizes(i), grid_dim_ids(i))
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
    end do 

    ! grid_data
    ierr = nf90_def_var(file_id, "grid data", NF90_DOUBLE, grid_dim_ids, grid_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !--------------------------!
    ! Finish defining metadata !
    !--------------------------!
    ierr = nf90_enddef(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ierr = nf90_put_var(file_id, grid_id, grid_data) 
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !----------------!
    ! Close the file !
    !----------------!
    ierr = nf90_close(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

  end subroutine ncdf_writer_3d

  subroutine ncdf_writer_4d(filename, ierr, grid_data)

    integer, parameter :: grid_ndims = 4

    ! Data to write to file
    real(real64), dimension(:,:,:,:), allocatable, intent(in) :: grid_data
    ! Number of dimensions of my rho data
    integer, dimension(grid_ndims) :: grid_sizes, grid_dim_ids

    ! Names of my dimensions
    character(len=1), dimension(grid_ndims) :: grid_dims=(/"x", "y" , "z", "t"/)

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    ! Variables used in writing process
    integer :: ierr, file_id, i

    ! Ids for variables
    integer :: grid_id

    ! Get the sizes of my incoming arrays
    grid_sizes  = shape(grid_data)

    !-------------------------------------------!
    ! Create the file, overwriting if it exists !
    !-------------------------------------------!
    ierr = nf90_create(filename, nf90_clobber, file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !----------------------------------------!
    ! Define the 3D variables and dimensions !
    !----------------------------------------!
    do i = 1, grid_ndims
      ierr = nf90_def_dim(file_id, grid_dims(i), grid_sizes(i), grid_dim_ids(i))
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
    end do 

    ! grid_data
    ierr = nf90_def_var(file_id, "grid data", NF90_DOUBLE, grid_dim_ids, grid_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !--------------------------!
    ! Finish defining metadata !
    !--------------------------!
    ierr = nf90_enddef(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ierr = nf90_put_var(file_id, grid_id, grid_data) 
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !----------------!
    ! Close the file !
    !----------------!
    ierr = nf90_close(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

  end subroutine ncdf_writer_4d

  subroutine ncdf_writer_3d_short(filename, ierr, grid_data, energies)

    real(real64), dimension(:), allocatable, intent(in) :: energies

    ! Number of dimensions of my grid data
    integer :: energy_size, energy_dim_id

    ! Names of my dimensions
    character(len=6) :: energy_dims = "energy"

    ! Ids for variables
    integer :: energy_id

    integer, parameter :: grid_ndims = 3

    ! Data to write to file
    integer(int16), dimension(:,:,:), allocatable, intent(in) :: grid_data
    ! Number of dimensions of my rho data
    integer, dimension(grid_ndims) :: grid_sizes, grid_dim_ids

    ! Names of my dimensions
    character(len=1), dimension(grid_ndims) :: grid_dims=(/"x", "y" , "z"/)

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    ! Variables used in writing process
    integer :: ierr, file_id, i

    ! Ids for variables
    integer :: grid_id

    ! Get the sizes of my incoming arrays
    grid_sizes  = shape(grid_data)
    energy_size = size(energies)

    !-------------------------------------------!
    ! Create the file, overwriting if it exists !
    !-------------------------------------------!
    ierr = nf90_create(filename, nf90_clobber, file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !----------------------------------------!
    ! Define the 3D variables and dimensions !
    !----------------------------------------!
    do i = 1, grid_ndims
      ierr = nf90_def_dim(file_id, grid_dims(i), grid_sizes(i), grid_dim_ids(i))
      if (ierr /= nf90_noerr) then
        print*, trim(nf90_strerror(ierr))
        return
      end if
    end do 

    ! grid_data
    ierr = nf90_def_var(file_id, "grid data", NF90_SHORT, grid_dim_ids, grid_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ! Define the energy variable and dimensions !
    ierr = nf90_def_dim(file_id, energy_dims, energy_size, energy_dim_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if
  
    ierr = nf90_def_var(file_id, "energy", NF90_DOUBLE, energy_dim_id, energy_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !--------------------------!
    ! Finish defining metadata !
    !--------------------------!
    ierr = nf90_enddef(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ierr = nf90_put_var(file_id, grid_id, grid_data) 
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    ierr = nf90_put_var(file_id, energy_id, energies) 
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

    !----------------!
    ! Close the file !
    !----------------!
    ierr = nf90_close(file_id)
    if (ierr /= nf90_noerr) then
      print*, trim(nf90_strerror(ierr))
      return
    end if

  end subroutine ncdf_writer_3d_short

  subroutine check(stat)
    !> integer error code
    integer, intent ( in) :: stat

    !> check for error
    if(stat /= nf90_noerr) then
      print *, trim(nf90_strerror(stat))
      stop "stopped"
    end if
  end subroutine check

  !--------------------------------------------------------------------!
  ! Subroutine to read and parse bin edge netcdf file                  !
  !                                                                    !
  ! H. J. Naguszewski, Warwick                                    2024 !
  !--------------------------------------------------------------------!
  subroutine read_1D_array(filename, varname, array)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: varname
    real(real64), intent(out) :: array(:)
    integer :: ncid, varid, ierr

    ! Open the NetCDF file
    ierr = nf90_open(filename, nf90_nowrite, ncid)
    if (ierr /= nf90_noerr) then
      print *, 'Error opening file:', nf90_strerror(ierr)
      return
    end if

    ! Get the variable ID
    ierr = nf90_inq_varid(ncid, varname, varid)
    if (ierr /= nf90_noerr) then
      print *, 'Error getting variable ID:', nf90_strerror(ierr)
    end if

    ! Read the data into the array
    ierr = nf90_get_var(ncid, varid, array)
    if (ierr /= nf90_noerr) then
      print *, 'Error reading variable:', nf90_strerror(ierr)
    end if

    ! Close the NetCDF file
    ierr = nf90_close(ncid)
    if (ierr /= nf90_noerr) then
      print *, 'Error closing file:', nf90_strerror(ierr)
    end if

  end subroutine read_1D_array
   
end module write_netcdf
