!> @file    io.f90
!>
!> @brief   Assorted routines and tools for file and data i/o
!>
!> @details This module contains routines for reading/writing
!>          information about a simulation (either to the screen or to
!>          file in plain text format. Data to be stored in binary
!>          format is written using the NetCDF library, for which the
!>          relevant routines can be found in write_netcdf.f90
!>
!> @author  C. D. Woodgate
!> @date    2020-2025
module io

  use kinds
  use derived_types
  use shared_data
  use command_line
  use display
  use comms
  
  implicit none

  private

  public :: write_info, make_data_directories, read_control_file,      &
            echo_control_file, read_exchange, parse_inputs,            &
            parse_metropolis_inputs, read_metropolis_file,             &
            echo_metropolis_file, read_ns_file, read_tmmc_file,        &
            read_wl_file

  ! Variables for keeping track of time
  real(real64) :: t_start, t_stop

  contains

  !> @brief   Subroutine to print software version, date, and time
  !>
  !> @details Some ASCII-art, names of developers/contributors, and
  !>          date/time when execution of the program started and
  !>          finished.
  !>
  !> @author  C. D. Woodgate
  !> @date    2024-2025
  !>
  !> @param  point Character string (either 's' or 'f') telling us
  !>               whether we are at the start or finish of simulation
  !>
  !> @return None
  subroutine write_info(point)

    character (len=1) point
    character (len=10) date,time

    ! Get date and time
    call date_and_time(date=date,time=time)

    write(6,'(/,72("="))')
    write(6,'(22x,"BraWl Version 0.4.1, 15.04.25")')
    write(6,'(72("-"))')

    if (point .eq. 's') then

      ! Start the clock
      call cpu_time(t_start)

      ! Print relevant info to screen
      write(6, '(20x,"    ____            _       ____")')
      write(6, '(20x,"   / __ )_________ | |     / / /")')
      write(6, '(20x,"  / __  / ___/ __ `/ | /| / / / ")')
      write(6, '(20x," / /_/ / /  / /_/ /| |/ |/ / /  ")')
      write(6, '(20x,"/_____/_/   \__,_/ |__/|__/_/   ")')
      write(6, '(20x,"                                ")')
      write(6,'(72("-"))')
      write(6, '("      Authors: Hubert J. Naguszewski,  ")')
      write(6, '("               Livia B. Partay,        ")')
      write(6, '("               Christopher D. Woodgate")')
      write(6,'(72("-"))')
      write(6, '(" Contributors: Heather Ratcliffe")')
      write(6, '("               David Quigley    ")')
      write(6,'(72("-"))')
      write(6,'(15x,"This run started at",1x,a," on",1x,a)')           &
               time(1:2)//":"//time(3:4)//":"//time(5:6),              &
               date(7:8)//"."//date(5:6)//"."//date(1:4)
    else if(point .eq. 'f') then
      write(6,'(15x,"This run finished at",1x,a," on",1x,a)')          &
               time(1:2)//":"//time(3:4)//":"//time(5:6),              &
               date(7:8)//"."//date(5:6)//"."//date(1:4)
      write(6,'(72("-"))')
      ! Start the clock
      call cpu_time(t_stop)
      write(6,'(20x,"Execution took",1x,f9.1 " seconds")') t_stop-t_start
    endif

    write(6,'(72("="),/)' )

  end subroutine write_info

  !> @brief   Subroutine to make directories for storing data
  !>
  !> @author  C. D. Woodgate
  !> @date    2021-2025
  !>
  !> @param  point Character string (either 's' or 'f') telling us
  !>               whether we are at the start or finish of simulation
  !>
  !> @return None
  subroutine make_data_directories(my_rank)

    integer :: my_rank

    ! make a directory for the grid states, diagnostics, 
    ! and radial_densities for each thread
    if(my_rank == 0) call execute_command_line('mkdir -p configs')
    if(my_rank == 0) call execute_command_line('mkdir -p energies')
    if(my_rank == 0) call execute_command_line('mkdir -p asro')

  end subroutine make_data_directories

  !> @brief   Subroutine to read the input file defining the simulation
  !>
  !> @details As of v0.4.0, need to read this file AND a separate input
  !>          file relevant to Metropolis, Nested Sampling, Wang-Landau
  !>          algorithms.
  !>
  !> @author  C. D. Woodgate
  !> @date    2020-2025
  !>
  !> @param  filename Name of file to read
  !> @param  parameters Derived type containing simulation parameters
  !> @param  my_rank Rank of current MPI process
  !>
  !> @return None
  subroutine read_control_file(filename, parameters, my_rank)

    character(len=*), intent(in) :: filename
    logical, dimension(11) :: check
    logical, dimension(2) :: cs_or_counts
    type(run_params) :: parameters
    character(len=144) :: buffer, label, first_char
    integer :: line, pos, ios, my_rank
    logical :: exists

    ! Set my checking_array to be false initially
    check = .false.
    cs_or_counts = .false.

    ios=0; line=0

    ! Initialise to default values---this means Valgrid will be happy
    parameters%mode = 301
    parameters%lattice = 'fcc'
    parameters%lattice_parameter = 3.57
    parameters%n_1 = 4
    parameters%n_2 = 4
    parameters%n_3 = 4
    parameters%n_basis = 1
    parameters%n_species = 4
    parameters%interaction_file = 'V_ijs.txt'
    parameters%wc_range = 2
    parameters%static_seed = .true.

    ! See if the relevant file exists
    inquire(file=trim(filename), exist=exists)

    ! Exit cleanly if we can't find it
    if (.not. exists) then
      call comms_finalise()
      stop 'Could not find input file ' // trim(filename)
    end if

    ! If we can find it, output that we are reading it
    if(my_rank == 0) then
      write(6,'(26("-"),x,"Parsing input file",x,26("-"),/)')
    end if

    ! Open it for reading
    open(15, file=filename, iostat=ios)

    ! Read line by line from the buffer until we hit the end of the file
    do while (ios==0)

      read(15, "(A)", iostat=ios) buffer

      if(ios==0) then
        line=line+1

        ! Check if the first non-whitespace character is a hash.
        ! If so, this is a comment line---ignore it.
        first_char = trim(buffer)
        if (first_char(1:1) .eq. '#') then
          continue
        end if

        pos = scan(buffer, '=')
        label=buffer(1:pos-1)
        buffer = buffer(pos+1:)

        select case (label)
        case ('mode')
          read(buffer, *, iostat=ios) parameters%mode
          check(1) = .true.
        case ('lattice')
          read(buffer, *, iostat=ios) parameters%lattice
          check(2) = .true.
        case ('lattice_parameter')
          read(buffer, *, iostat=ios) parameters%lattice_parameter
          check(3) = .true.
        case ('n_1')
          read(buffer, *, iostat=ios) parameters%n_1
          check(4) = .true.
        case ('n_2')
          read(buffer, *, iostat=ios) parameters%n_2
          check(5) = .true.
        case ('n_3')
          read(buffer, *, iostat=ios) parameters%n_3
          check(6) = .true.
        case ('n_species')
          read(buffer, *, iostat=ios) parameters%n_species
          check(7) = .true.
        case ('interaction_file')
          read(buffer, *, iostat=ios) parameters%interaction_file
          check(8) = .true.
        case ('interaction_range')
          read(buffer, *, iostat=ios) parameters%interaction_range
          check(9) = .true.
        case ('wc_range')
          read(buffer, *, iostat=ios) parameters%wc_range
        case ('static_seed')
          read(buffer, *, iostat=ios) parameters%static_seed
        case default
        end select
      end if
    end do

    close(15)

    ! Now we know how many species there are, we can go back and read
    ! the concentrations/n_atoms of each species, and their labels
    allocate(parameters%species_names(parameters%n_species))
    allocate(parameters%species_concentrations(0:parameters%n_species))
    allocate(parameters%species_numbers(parameters%n_species))

    ! Set these arrays to be zero initially
    parameters%species_concentrations = 0.0_real64
    parameters%species_numbers = 0

    line=0

    open(15, file=filename, iostat=ios)

    ! Read line by line from the buffer until we hit the end of the file
    do while (ios==0)

      read(15, "(A)", iostat=ios) buffer

      if(ios==0) then
        line=line+1

        ! Check if the first non-whitespace character is a hash.
        ! If so, this is a comment line---ignore it.
        first_char = trim(buffer)
        if (first_char(1:1) .eq. '#') then
          continue
        end if

        pos = scan(buffer, '=')
        label=buffer(1:pos-1)
        buffer = buffer(pos+1:)

        select case (label)
        case ('species_names')
          read(buffer, *, iostat=ios) parameters%species_names
          check(10) = .true.
        case ('species_concentrations')
          read(buffer, *, iostat=ios) parameters%species_concentrations(1:)
          cs_or_counts(1) = .true.
          check(11) = .true.
        case ('species_numbers')
          read(buffer, *, iostat=ios) parameters%species_numbers(:)
          cs_or_counts(2) = .true.
          check(11) = .true.
        case default
        end select
      end if
    end do

    close(15)

    ! Check that the user has specified either the concentrations of each species
    ! or the number of atoms of each species, but not both.
    if (cs_or_counts(1) .and. cs_or_counts(2)) then
      call comms_finalise()
      stop 'You cannot specify both chemical concentrations and numbers of atoms!'
    else if ((.not. cs_or_counts(1)) .and. (.not. cs_or_counts(2))) then
      call comms_finalise()
      stop 'You must specify either chemical concentrations or numbers of atoms.'
    end if

    ! Check that the user has provided all the necessary inputs
    ! Exit and tell them which is missing if needed
    if (.not. all(check)) then
      call comms_finalise()
      if (.not. check(1)) then
        stop "Missing 'mode' in system file"
      else if (.not. check(2)) then
        stop "Missing 'lattice' in system file"
      else if (.not. check(3)) then
        stop "Missing 'lattice_parameter' in system file"
      else if (.not. check(4)) then
        stop "Missing 'n_1' in system file"
      else if (.not. check(5)) then
        stop "Missing 'n_2' in system file"
      else if (.not. check(6)) then
        stop "Missing 'n_3' in system file"
      else if (.not. check(7)) then
        stop "Missing 'n_species' in system file"
      else if (.not. check(8)) then
        stop "Missing 'interaction_file' in system file"
      else if (.not. check(9)) then
        stop "Missing 'interaction_range' in system file"
      else if (.not. check(10)) then
        stop "Missing 'species_names' in system file"
      else
      stop 'Missing parameter in system file'
      end if
    end if

  end subroutine read_control_file

  !> @brief   Subroutine to echo the contents of the input file to
  !>          the screen
  !>
  !> @author  C. D. Woodgate
  !> @date    2020-2025
  !>
  !> @param  parameters Derived type containing simulation parameters
  !>
  !> @return None
  subroutine echo_control_file(parameters)

    type(run_params) :: parameters
    integer :: i

    print*, ' Read mode = ', parameters%mode
    print*, ' Read n_1 = ', parameters%n_1
    print*, ' Read n_2 = ', parameters%n_2
    print*, ' Read n_3 = ', parameters%n_3
    print*, ' Read n_basis = ', parameters%n_basis
    print*, ' Read n_species = ', parameters%n_species
    print*, ' Read lattice parameter = ', parameters%lattice_parameter
    print*, ' Read lattice = ', parameters%lattice
    print*, ' Read interaction_file = ', parameters%interaction_file
    print*, ' Read wc_range = ', parameters%wc_range

    ! Print specified concentrations/numbers of atoms
    if (abs(sum(parameters%species_concentrations)-1.0_real64) &
        .lt. 0.001) then
      do i=1, parameters%n_species
        print*, ' Read species ', i, ' = ', parameters%species_names(i), &
                ' at concentration ', parameters%species_concentrations(i)
      enddo
    else
      do i=1, parameters%n_species
        print*, ' Read species ', i, ' = ', parameters%species_names(i), &
                ' at ', parameters%species_numbers(i), ' atoms'
      enddo
    end if

  end subroutine echo_control_file

  !> @brief   Subroutine to read the atom-atom effective pair
  !>          interactions (EPIs) from file
  !>
  !> @details The file should represent the EPIs on each coordination
  !>          shell as an sxs matrix, with a blank line between the set
  !>          of EPIs for one coordination shell and those for the
  !>          next. They should be given in order from nearest to
  !>          furthest coordination shell.
  !>
  !> @author  C. D. Woodgate
  !> @date    2020-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  my_rank Rank of current MPI process
  !>
  !> @return None
  subroutine read_exchange(setup, my_rank)

    type(run_params) , intent(in) :: setup
    integer :: my_rank

    if(my_rank == 0) then
      !write(6,'(72("-"),/)', advance='no')
      write(6,'(8("-"),x,"Reading atom-atom effective pair interaction parameters",x, 7("-"),/)')
    end if

    V_ex = 0.0_real64
    open(16, file=setup%interaction_file)
    read(16,*) V_ex   
    close(16)

    if (my_rank .eq. 0) then
      ! Print it out as a sanity check for user
      call pretty_print_exchange(setup)
    end if

    if(my_rank == 0) then
      write(6,'(3("-"),x,"Read atom-atom effective pair interaction parameters successfully",x, &
                &2("-"),/)')!, advance='no')
      !write(6,'(72("-"),/)')
    end if

  end subroutine read_exchange

  !> @brief   Subroutine to parse command line arguments
  !>
  !> @details For routines which actually *parse* the command line
  !>          arguments, see command_line.f90
  !>
  !> @author  C. D. Woodgate
  !> @date    2020-2025
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  my_rank Rank of current MPI process
  !>
  !> @return None
  subroutine parse_inputs(setup, my_rank)

    type(run_params) :: setup
    integer :: my_rank
    character(len=30) :: control = ' '

    ! Parse the name of the input file
    if(my_rank == 0) then
      write(6,'(22("-"),x,"Parsing name of input file",x,22("-"),/)')
    end if
      if(.not. get_arg('control', control)) then
        if (my_rank == 0) then
          print*, 'Input file not specified with "control=<name>"'
          print*, 'Defaulting to searching for "input.txt"'
          print*, ' '
        end if
        control = 'input.txt'
      else
        if (my_rank == 0) then
          print*, 'Input file name is: ', control
          print*, ' '
        end if
      end if

    ! Read the input file
    call read_control_file(control, setup, my_rank)

    if(my_rank == 0) then
      call echo_control_file(setup)
      write(6,'(/,20("-"),x,"Parsed input file successfully",x,20("-"),/)')
    end if

  end subroutine parse_inputs

  !> @brief   Subroutine to parse the Metropolis MC input file
  !>
  !> @author  C. D. Woodgate
  !> @date    2020-2025
  !>
  !> @param  metropolis Derived type containing Metropolis MC parameters
  !> @param  my_rank Rank of current MPI process
  !>
  !> @return None
  subroutine parse_metropolis_inputs(metropolis, my_rank)

    type(metropolis_params) :: metropolis
    integer :: my_rank
    character(len=30) :: control = ' '

    ! Parse the name of the input file
    if(my_rank == 0) then
      write(6,'(/,17("-"),x,"Parsing name of Metropolis input file",x,16("-"),/)')
    end if
      if(.not. get_arg('metropolis', control)) then
        if (my_rank == 0) then
          print*, ' Metropolis input file not specified with "metropolis=<name>"'
          print*, ' '
          print*, ' Defaulting to searching for "metropolis.inp"'
        end if
        control = 'metropolis.inp'
      else
        if (my_rank == 0) then
          print*, 'Input file name is: ', control
          print*, ' '
        end if
      end if

    ! Read the input file
    call read_metropolis_file(control, metropolis, my_rank)

    write(6,'(x,"Parameters to be used are as follows",/)')

    if(my_rank == 0) then
      call echo_metropolis_file(metropolis)
      write(6,'(/,15("-"),x,"Parsed Metropolis input file successfully",x,14("-"),/)')
    end if

  end subroutine parse_metropolis_inputs

  !> @brief   Subroutine to parse the Metropolis MC input file
  !>
  !> @author  C. D. Woodgate
  !> @date    2025
  !>
  !> @param  metropolis Derived type containing Metropolis MC parameters
  !> @param  my_rank Rank of current MPI process
  !>
  !> @return None
  subroutine read_metropolis_file(filename, metropolis, my_rank)

    character(len=*), intent(in) :: filename
    logical, dimension(4) :: check
    type(metropolis_params) :: metropolis
    character(len=144) :: buffer, label, first_char
    integer :: line, pos, ios, my_rank
    logical :: exists

    ! Set my checking_array to be false initially
    check = .false.

    ios=0; line=0

    ! Defaults if these are not specified.
    metropolis%burn_in_start = .False.
    metropolis%burn_in = .False.
    metropolis%n_burn_in_steps = 0
    metropolis%calculate_energies = .true.
    metropolis%write_trajectory_energy = .false.
    metropolis%calculate_asro = .true.
    metropolis%calculate_alro = .false.
    metropolis%n_sample_steps_asro = 0
    metropolis%n_sample_steps_alro = 0
    metropolis%write_trajectory_xyz = .false.
    metropolis%write_trajectory_energy = .false.
    metropolis%write_trajectory_asro = .false.
    metropolis%n_sample_steps_trajectory = 0
    metropolis%write_final_config_xyz = .false.
    metropolis%read_start_config_nc = .false.
    metropolis%T_steps = 1
    metropolis%delta_T = 1
    metropolis%nbr_swap = .false.

    ! See if the relevant file exists
    inquire(file=trim(filename), exist=exists)

    ! Exit cleanly if we can't find it
    if (.not. exists) then
      call comms_finalise()
      stop 'Could not find Metropolis control file: ' // trim(filename)
    end if

    ! If we can find it, output that we are reading it
    if(my_rank == 0) then
      write(6,'(/,19("-"),x,"Parsing Metropolis control file",x,20("-"),/)')
    end if

    ! Open it for reading
    open(15, file=filename, iostat=ios)

    ! Read line by line from the buffer until we hit the end of the file
    do while (ios==0)

      read(15, "(A)", iostat=ios) buffer

      if(ios==0) then
        line=line+1

        ! Check if the first non-whitespace character is a hash.
        ! If so, this is a comment line---ignore it.
        first_char = trim(buffer)
        if (first_char(1:1) .eq. '#') then
          continue
        end if

        pos = scan(buffer, '=')
        label=buffer(1:pos-1)
        buffer = buffer(pos+1:)

        select case (label)
        case ('mode')
          read(buffer, *, iostat=ios) metropolis%mode
          check(1) = .true.
        case ('n_mc_steps')
          read(buffer, *, iostat=ios) metropolis%n_mc_steps
          check(2) = .true.
        case ('burn_in_start')
          read(buffer, *, iostat=ios) metropolis%burn_in_start
        case ('burn_in')
          read(buffer, *, iostat=ios) metropolis%burn_in
        case ('n_burn_in_steps')
          read(buffer, *, iostat=ios) metropolis%n_burn_in_steps
        case ('calculate_energies')
          read(buffer, *, iostat=ios) metropolis%calculate_energies
        case ('n_sample_steps')
          read(buffer, *, iostat=ios) metropolis%n_sample_steps
          check(3) = .true.
        case ('calculate_asro')
          read(buffer, *, iostat=ios) metropolis%calculate_asro
        case ('n_sample_steps_asro')
          read(buffer, *, iostat=ios) metropolis%n_sample_steps_asro
        case ('calculate_alro')
          read(buffer, *, iostat=ios) metropolis%calculate_alro
        case ('n_sample_steps_alro')
          read(buffer, *, iostat=ios) metropolis%n_sample_steps_alro
        case ('n_sample_steps_trajectory')
          read(buffer, *, iostat=ios) metropolis%n_sample_steps_trajectory
        case ('write_trajectory_xyz')
          read(buffer, *, iostat=ios) metropolis%write_trajectory_xyz
        case ('write_trajectory_energy')
          read(buffer, *, iostat=ios) metropolis%write_trajectory_energy
        case ('write_trajectory_asro')
          read(buffer, *, iostat=ios) metropolis%write_trajectory_asro
        case ('write_final_config_xyz')
          read(buffer, *, iostat=ios) metropolis%write_final_config_xyz
        case ('write_final_config_nc')
          read(buffer, *, iostat=ios) metropolis%write_final_config_nc
        case ('read_start_config_nc')
          read(buffer, *, iostat=ios) metropolis%read_start_config_nc
        case ('start_config_file')
          read(buffer, *, iostat=ios) metropolis%start_config_file
        case ('T')
          read(buffer, *, iostat=ios) metropolis%T
          check(4) = .true.
        case ('T_steps')
          read(buffer, *, iostat=ios) metropolis%T_steps
        case ('delta_T')
          read(buffer, *, iostat=ios) metropolis%delta_T
        case ('nbr_swap')
          read(buffer, *, iostat=ios) metropolis%nbr_swap
        case default
        end select
      end if
    end do

    close(15)

    if (metropolis%n_sample_steps_alro .eq. 0) then
      metropolis%n_sample_steps_asro = metropolis%n_sample_steps
    endif
    if (metropolis%n_sample_steps_alro .eq. 0) then
      metropolis%n_sample_steps_alro = metropolis%n_sample_steps
    endif
    if (metropolis%n_sample_steps_trajectory .eq. 0) then
      metropolis%n_sample_steps_trajectory = metropolis%n_sample_steps
    endif

    ! Check that the user has provided all the necessary inputs
    ! Exit and tell them which is missing if needed
    if (.not. all(check)) then
      call comms_finalise()
      if (.not. check(1)) then
        stop "Missing 'mode' in Metropolis input file"
      else if (.not. check(2)) then
        stop "Missing 'n_mc_steps' in Metropolis input file"
      else if (.not. check(3)) then
        stop "Missing 'n_sample_steps' in Metropolis input file"
      else if (.not. check(4)) then
        stop "Missing 'T' in Metropolis input file"
      else
        stop "Missing something in Metropolis input file"
      end if
    end if

  end subroutine read_metropolis_file

  !> @brief   Subroutine to echo the contents of the Metropolis input
  !>          file to the screen
  !>
  !> @author  C. D. Woodgate
  !> @date    2025
  !>
  !> @param  metropolis Derived type containing Metroplis parameters
  !>
  !> @return None
  subroutine echo_metropolis_file(metropolis)

    type(metropolis_params) :: metropolis

    print*, ' mode =                       ', metropolis%mode
    print*, ' n_mc_steps =                 ', metropolis%n_mc_steps
    print*, ' burn_in_start =              ', metropolis%burn_in_start
    print*, ' burn_in =                    ', metropolis%burn_in
    print*, ' n_burn_in_steps =            ', metropolis%n_burn_in_steps
    print*, ' n_sample_steps =             ', metropolis%n_sample_steps
    print*, ' calculate_energies =         ', metropolis%calculate_energies
    print*, ' write_trajectory_energy =    ', metropolis%write_trajectory_energy
    print*, ' write_trajectory_asro =      ', metropolis%write_trajectory_asro
    print*, ' calculate_asro =             ', metropolis%calculate_asro
    print*, ' n_sample_steps_asro =        ', metropolis%n_sample_steps_asro
    print*, ' calculate_alro =             ', metropolis%calculate_alro
    print*, ' n_sample_steps_alro =        ', metropolis%n_sample_steps_alro
    print*, ' write_trajectory_xyz =       ', metropolis%write_trajectory_xyz
    print*, ' n_sample_steps_trajectory =  ', metropolis%n_sample_steps_trajectory
    print*, ' write_final_config_xyz =     ', metropolis%write_final_config_xyz
    print*, ' write_final_config_nc =      ', metropolis%write_final_config_nc
    print*, ' read_start_config_nc =       ', metropolis%read_start_config_nc
    if (metropolis%read_start_config_nc) then
      print*, ' starting configuration file  ', metropolis%start_config_file
    end if
    print*, ' T =                          ', metropolis%T
    print*, ' T_steps =                    ', metropolis%T_steps
    print*, ' delta_T =                    ', metropolis%delta_T
    print*, ' nbr_swap =                   ', metropolis%nbr_swap

  end subroutine echo_metropolis_file

  !> @brief   Subroutine to read and parse nested sampling control file
  !>
  !> @author  L. B. Partay
  !> @date    2024
  !>
  !> @param  filename Name of the nested sampling input file (expected
  !>                  to be "ns_input.txt")
  !> @param  parameters Derived type of ns_params, containing nested
  !>                    sampling parameters
  !>
  !> @return None
  subroutine read_ns_file(filename, parameters)

    character(len=*), intent(in) :: filename
    logical, dimension(8) :: check
    type(ns_params) :: parameters
    character(len=144) :: buffer, label, first_char
    integer :: line, pos, ios

    check = .false.

    ios=0; line=0

    ! Print to screen that we are looking for the NS input file
    if(my_rank == 0) then
      write(6,'(/,18("-"),x,"Parsing Nested Sampling input file",x,18("-"),/)')
    end if

    print*, ' Looking for Nested Sampling input file named: ', filename, new_line('a')

    open(25, file=filename, iostat=ios)

    ! Exit cleanly if we cannot find it
    if (ios .ne. 0) then
      call comms_finalise()
      stop 'Could not parse input file. Aborting...'
    end if

    ! Otherwise, read it
    do while (ios==0)

      read(25, "(A)", iostat=ios) buffer

      if(ios==0) then
        line=line+1

        ! Check if the first non-whitespace character is a hash.
        ! If so, this is a comment line---ignore it.
        first_char = trim(buffer)
        if (first_char(1:1) .eq. '#') then
          continue
        end if

        pos = scan(buffer, '=')
        label=buffer(1:pos-1)
        buffer = buffer(pos+1:)

        select case (label)
        case ('n_walkers')
          read(buffer, *, iostat=ios) parameters%n_walkers
          print*, ' Read n_walkers = ', parameters%n_walkers
        case ('n_steps')
          read(buffer, *, iostat=ios) parameters%n_steps
          print*, ' Read n_steps = ', parameters%n_steps
        case ('n_iter')
          read(buffer, *, iostat=ios) parameters%n_iter
          print*, ' Read n_iter = ', parameters%n_iter
        case ('outfile_ener')
          read(buffer, *, iostat=ios) parameters%outfile_ener
          print*, ' Read outfile_ener = ', parameters%outfile_ener
        case ('outfile_traj')
          read(buffer, *, iostat=ios) parameters%outfile_traj
          print*, ' Read outfile_traj = ', parameters%outfile_traj
        case ('traj_freq')
          read(buffer, *, iostat=ios) parameters%traj_freq
          print*, ' Write configuration every n-th NS iteration = ', parameters%traj_freq
        case default
        end select
      end if
    end do

    ! If all has gone to plan, print that we have read this file
    if(my_rank == 0) then
      write(6,'(/,12("-"),x,"Successfully parsed Nested Sampling input file",x,12("-"),/)')
    end if

    close(25)

  end subroutine read_ns_file

  !> @brief   Subroutine to read and parse transition matrix Monte Carlo
  !>          (TMMC) input file
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024
  !>
  !> @param  filename Name of the TMMC input file
  !> @param  parameters Derived type of tmmc_params, containing TMMC
  !>                    parameters
  !>
  !> @return None
  subroutine read_tmmc_file(filename, parameters, my_rank)

    integer :: my_rank
    character(len=*), intent(in) :: filename
    logical, dimension(8) :: check
    type(tmmc_params) :: parameters
    character(len=100) :: buffer, label
    integer :: line, pos, ios

    check = .false.

    ios=0; line=0

    open(25, file=filename, iostat=ios)

    if (ios .ne. 0) then
      call comms_finalise()
      stop 'Could not parse tmmc input file. Aborting...'
    end if

    if (my_rank == 0) then
      write(*,'(a)', advance='no') new_line('a')
      print*, '###############################'
      print*, '#  Parsing tmmc input file    #'
      print*, '###############################'

      print*, '# tmmc input file name: ', filename
    end if

    do while (ios==0)

      read(25, "(A)", iostat=ios) buffer

      if(ios==0) then
        line=line+1

        pos = scan(buffer, '=')
        label=buffer(1:pos-1)
        buffer = buffer(pos+1:)

        select case (label)
        case ('mc_sweeps')
          read(buffer, *, iostat=ios) parameters%mc_sweeps
          if (my_rank == 0) then
            print*, '# Read mc_sweeps = ', parameters%mc_sweeps
          end if
          check(1) = .true.
        case ('bins')
          read(buffer, *, iostat=ios) parameters%bins
          if (my_rank == 0) then
            print*, '# Read bins = ', parameters%bins
          end if
          check(2) = .true.
        case ('num_windows')
          read(buffer, *, iostat=ios) parameters%num_windows
          if (my_rank == 0) then
            print*, '# Read num_windows = ', parameters%num_windows
          end if
          check(3) = .true.
        case ('bin_overlap')
          read(buffer, *, iostat=ios) parameters%bin_overlap
          if (my_rank == 0) then
            print*, '# Read bin_overlap = ', parameters%bin_overlap
          end if
          check(4) = .true.
        case ('weight_update')
          read(buffer, *, iostat=ios) parameters%weight_update
          if (my_rank == 0) then
            print*, '# Read weight_update = ', parameters%weight_update
          end if
          check(5) = .true.
        case ('energy_min')
          read(buffer, *, iostat=ios) parameters%energy_min
          if (my_rank == 0) then
            print*, '# Read energy_min = ', parameters%energy_min
          end if
          check(6) = .true.
        case ('energy_max')
          read(buffer, *, iostat=ios) parameters%energy_max
          if (my_rank == 0) then
            print*, '# Read energy_max = ', parameters%energy_max
          end if
          check(7) = .true.
        case ('T')
          read(buffer, *, iostat=ios) parameters%T
          if (my_rank == 0) then
            print*, '# Read T = ', parameters%T
          end if
          check(8) = .true.
        case default
          if (my_rank == 0) then
            print*, '# Skipping invalid label'
          end if
        end select
      end if
    end do

    if (my_rank == 0) then
      print*, '# Finished parsing tmmc input file #'
      print*, '####################################', new_line('a')
    end if
    close(25)

    if (.not. all(check)) then
      call comms_finalise()
      stop 'Missing parameter in tmmc input file'
    end if

  end subroutine read_tmmc_file

  !> @brief   Subroutine to read and parse Wang-Landau sampling input
  !>          file
  !>
  !> @author  H. J. Naguszewski
  !> @date    2024
  !>
  !> @param  filename Name of the Wang-Landau input file
  !> @param  parameters Derived type of wl_params, containing
  !>                    Wang-Landau parameters
  !>
  !> @return None
  subroutine read_wl_file(filename, parameters, my_rank)

    integer :: my_rank
    character(len=*), intent(in) :: filename
    logical, dimension(10) :: check
    type(wl_params) :: parameters
    character(len=100) :: buffer, label
    integer :: line, pos, ios

    check = .false.

    ios=0; line=0

    open(25, file=filename, iostat=ios)

    if (ios .ne. 0) then
      call comms_finalise()
      stop 'Could not parse wang landau input file. Aborting...'
    end if

    if (my_rank == 0) then
      write(*,'(a)', advance='no') new_line('a')
      print*, '##################################'
      print*, '# Parsing wang landau input file #'
      print*, '##################################'

      print*, '# wang landau input file name: ', filename
    end if

    do while (ios==0)

      read(25, "(A)", iostat=ios) buffer

      if(ios==0) then
        line=line+1

        pos = scan(buffer, '=')
        label=buffer(1:pos-1)
        buffer = buffer(pos+1:)

        select case (label)
        case ('mc_sweeps')
          read(buffer, *, iostat=ios) parameters%mc_sweeps
          if (my_rank == 0) then
            print*, '# Read mc_sweeps = ', parameters%mc_sweeps
          end if
          check(1) = .true.
        case ('bins')
          read(buffer, *, iostat=ios) parameters%bins
          if (my_rank == 0) then
            print*, '# Read bins = ', parameters%bins
          end if
          check(2) = .true.
        case ('num_windows')
          read(buffer, *, iostat=ios) parameters%num_windows
          if (my_rank == 0) then
            print*, '# Read num_windows = ', parameters%num_windows
          end if
          check(3) = .true.
        case ('bin_overlap')
          read(buffer, *, iostat=ios) parameters%bin_overlap
          if (my_rank == 0) then
            print*, '# Read bin_overlap = ', parameters%bin_overlap
          end if
          check(4) = .true.
        case ('tolerance')
          read(buffer, *, iostat=ios) parameters%tolerance
          if (my_rank == 0) then
            print*, '# Read tolerance = ', parameters%tolerance
          end if
          check(5) = .true.
        case ('flatness')
          read(buffer, *, iostat=ios) parameters%flatness
          if (my_rank == 0) then
            print*, '# Read flatness = ', parameters%flatness
          end if
          check(6) = .true.
        case ('wl_f')
          read(buffer, *, iostat=ios) parameters%wl_f
          if (my_rank == 0) then
            print*, '# Read wl_f = ', parameters%wl_f
          end if
          check(7) = .true.
        case ('energy_min')
          read(buffer, *, iostat=ios) parameters%energy_min
          if (my_rank == 0) then
            print*, '# Read energy_min = ', parameters%energy_min
          end if
          check(8) = .true.
        case ('energy_max')
          read(buffer, *, iostat=ios) parameters%energy_max
          if (my_rank == 0) then
            print*, '# Read energy_max = ', parameters%energy_max
          end if
          check(9) = .true.
        case ('radial_samples')
          read(buffer, *, iostat=ios) parameters%radial_samples
          if (my_rank == 0) then
            print*, '# Read radial_samples = ', parameters%radial_samples
          end if
          check(10) = .true.

        case default
          if (my_rank == 0) then
            print*, '# Skipping invalid label'
          end if
        end select
      end if
    end do

    if (my_rank == 0) then
      print*, '# Finished parsing wang landau input file #'
      print*, '###########################################', new_line('a')
    end if
    close(25)

    if (.not. all(check)) then
      call comms_finalise()
      stop 'Missing parameter in wang landau input file'
    end if

  end subroutine read_wl_file

end module io
