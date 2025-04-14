!----------------------------------------------------------------------!
! Configuration Output Module                                          !
!                                                                      !
! This program runs MC for the loaded temperatures and outputs         !
! the configuration once equilibrium is achieved                       !
!                                                                      !
! H. J. Naguszewski, Warwick                                      2024 !
!----------------------------------------------------------------------!

module config_output
  use initialise
  use kinds
  use shared_data
  use c_functions
  use random_site
  use metropolis
  use mpi

  implicit none

  contains

  subroutine co_main(setup, my_rank)
    ! Rank of this processor
    integer, intent(in) :: my_rank

    ! Arrays for storing data
    type(run_params) :: setup

    ! Integers used in calculations
    integer :: i, j, acceptance
    real(real64) :: beta

    ! Temperature and temperature steps
    real(real64) :: energy_to_ry

    ! Temperature file and arrays
    character(len=100) :: filename, xyz_file
    real(real64), allocatable :: temperatures(:)
    integer :: num_temp

    filename = "temperatures.txt"
    print*, '##################################'
    print*, '# Parsing temperature input file #'
    print*, '##################################'

    call read_temp_file(filename, temperatures, num_temp)

    print*, "# Temperatures Found"
    do i=1, num_temp
      write(6, '(f8.2,a)') temperatures(i), " K" 
    end do

    energy_to_ry = setup%n_atoms/(Ry_to_eV*1000)

    ! Set up the lattice
    call initial_setup(setup, config)

    call lattice_shells(setup, shells, config)

    setup%mc_step => monte_carlo_step_lattice

    if (my_rank == 0) then
      write (6, '(/,72("-"),/)')
      write (6, '(24("-"),x,"Commencing Simulation!",x,24("-"),/)')
      print *, "Number of atoms", setup%n_atoms
    end if

    do i=1, num_temp
      write(6, '(a,f8.2,a)') "Running for temperature: ", temperatures(i), " K"
      beta = 1.0_real64/(temperatures(i)*k_b_in_Ry)
      do j=1, setup%n_atoms*50000
        acceptance = setup%mc_step(config, beta)
      end do
      write(xyz_file, '(A12 I4.4 F3.2 A4)') 'grids/config_at_T_', int(temperatures(i)), &
            temperatures(i)-int(temperatures(i)),'.xyz'
      call xyz_writer(xyz_file, config, setup)
    end do

    if (my_rank == 0) then
      write (*, *)
      write (6, '(25("-"),x,"Simulation Complete!",x,25("-"))')
    end if

  end subroutine co_main

  subroutine read_temp_file(filename, array, n)
    real(real64), dimension(:), allocatable, intent(out) :: array
    integer, intent(out) :: n
    character(len=100), intent(in) :: filename
    integer :: i
    integer :: unit
    integer :: line_count
    character(len=100) :: line
  
    ! Initialize the line counter
    line_count = 0
    unit = 0
  
    ! First pass: Open the file and count the number of lines (elements)
    open(unit=unit, file=filename, status='old')
  
    do while (.true.)
       read(unit,'(A)', iostat=i) line
       if (i /= 0) exit  ! End of file reached
       line_count = line_count + 1
    end do
  
    ! Now allocate the array based on the line count
    n = line_count
    allocate(array(n))
  
    ! Rewind and read the file again to populate the array
    rewind(unit)
  
    ! Read the data into the array
    do i = 1, n
       read(unit, *) array(i)
    end do
  
    ! Close the file
    close(unit)
  
  end subroutine read_temp_file

end module config_output

