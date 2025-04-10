!> @file    command_line.f90
!>
!> @brief   Routines for parsing command line arguments.
!>
!> @details Code originally written by H. Ratcliffe and C. S. Brady
!>          for the Research Software Engineering course delivered as
!>          part of the taught component of the EPSRC-supported Centre
!>          for Doctoral Training (CDT) in Modelling of Heterogeneous
!>          Systems (HetSys) at the University of Warwick. Re-used here
!>          with permission.
!>
!> @author  H. Ratcliffe
!> @author  C. S. Brady
!> @author  C. D. Woodgate
!>
!> @date    2019-2025
module command_line

  use kinds

  implicit none

  ! length of character value string
  integer, parameter :: len = 30

  ! type containing a key-value pair
  ! init. to default values
  type cmd_arg
    character(len=20) :: name = "NULL"
    character(len=len) :: value = ""
  end type

  !> @brief Read arguments by name or number
  !>
  !> @details All of the members of this interface take either a 'name'
  !>          or a 'num' parameter, detailling which argument to look
  !>          up and the type of the 'val' argument dictates the type to
  !>          attempt parsing the argument as.
  !>
  !>          Interface dispatches to the correct actual function based
  !>          on arguments in call. This means we can convert to a
  !>          requested type. We call like e.g.: get_arg(10, var) or
  !>          get_arg("name", var) and we will try and interpret the
  !>          value as whatever type var is. We can supply optional arg
  !>          to those to capture whether or not the name exists.
  !>
  !>          Note that if a name was present, but the value could not
  !>          be parsed, exists will be TRUE but correct will be FALSE.
  !>
  !> @author  H. Ratcliffe
  !> @author  C. S. Brady
  !>
  !> @date    2019-2021
  !>
  !> @param name Name to look up (if used)
  !> @param num Arg. number to look up (if used)
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !>
  !> @return True if the name is found and parsed, False otherwise
  interface get_arg
    module procedure get_arg_num_int, get_arg_name_int
    module procedure get_arg_num_dbl, get_arg_name_dbl
    module procedure get_arg_num_str, get_arg_name_str
  end interface

  ! Module level variables - these are accessible
  ! to all functions in the module
  ! but we make them private to within the module only

  ! the argument list
  type(cmd_arg), dimension(:), allocatable :: all_args

  ! the number of arguments
  integer :: num_args = 0
  private :: all_args, num_args

  contains

  !> @brief Parse command line args
  !>
  !> @details This function can be called multiple times and will
  !>          freshly parse ALL arguments each time. We assume
  !>          these are entered as 'name=value' and if there is no '='
  !>          sign, then value is set to a sentinel. For most
  !>          flexibility, we assume all values are reals, and thus we
  !>          use 'HUGE(1.0)' as the no-value sentinel. This number is
  !>          greater than any other real value or a 'num' parameter,
  !>          detailling which argument to look up and the type of the
  !>          'val' argument dictates the type to attempt parsing the
  !>          argument as.
  !>
  !> @author  H. Ratcliffe
  !> @author  C. S. Brady
  !>
  !> @date    2019-2021
  !>
  !> @return None
  subroutine parse_args()

    ! Strictly we can't be sure 50 chars is enough
    ! but for command-line args it's enough if we're sensible
    ! We wont overflow, but our strings may get truncated
    CHARACTER(LEN=51) :: arg
    INTEGER :: i, indx

    num_args = command_argument_count()
    if(num_args > 0) then

      ! If this is not the first call, all_args may be already allocated
      ! Deallocate if needed, and allocate to correct size
      if(allocated(all_args)) deallocate(all_args)
      allocate(all_args(num_args))

      ! loop over all arguments
      do i = 1, num_args
        call get_command_argument(i, arg)
        ! Location of the '=' sign
        ! If not found, return value is 0
        indx = index(arg, '=')
        if(indx > 1) then
          ! All characters up to '=', not including it
          ! but with any leading spaces removed
          all_args(i)%name = adjustl(arg(1:indx-1))
          ! All characters after '='
          all_args(i)%value= adjustl(arg(indx+1:))
       else
          all_args(i)%name = trim(adjustl(arg))
          ! Value already has a default value, so leave it alone
        end if
      end do
    endif

  end subroutine parse_args

  !--------------------------------------------------------------------!
  ! The next functions let you access the parsed args. You need to     !
  ! have called parse_args first or there wont be any args! You don't  !
  ! need to call these explicitly, you should just use get_arg and the !
  ! INTERFACE declared above will handle the rest.                     !
  !--------------------------------------------------------------------!

  !> @brief Read by number for double precision values
  !>
  !> @author  H. Ratcliffe
  !> @author  C. S. Brady
  !>
  !> @date    2019-2021
  !>
  !> @param num Argument number to read
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !>
  !> @return True if the argument is found and parsed, False otherwise
  function get_arg_num_dbl(num, val, exists)

    logical :: get_arg_num_dbl
    integer, intent(in) :: num
    real(kind=real64), intent(out) :: val
    logical, intent(out), optional :: exists
    logical :: found
    integer :: ierr

    found = .false.
    ! Check requested number is in range
    if(num <= num_args .and. num > 0) then
      ! Read it from string into value
      ! We don't need to specify the format in general
      read(all_args(num)%value, *, iostat=ierr) val
      found = .true.
    end if

    if(present(exists)) then
      exists = found
    end if

    ! Return value is whether value is found and correctly parsed
    get_arg_num_dbl = (found .and. (ierr == 0))

  end function get_arg_num_dbl

  !> @brief Read by name for double precision values
  !>
  !> @author  H. Ratcliffe
  !> @author  C. S. Brady
  !>
  !> @date    2019-2021
  !>
  !> @param name Argument name to look up
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !>
  !> @return True if the name is found and parsed, False otherwise
  function get_arg_name_dbl(name, val, exists)

    logical :: get_arg_name_dbl
    character(len=*), intent(in) :: name
    real(kind=real64), intent(out) :: val
    integer :: i
    logical, intent(out), optional :: exists
    logical :: found
    integer :: ierr

    found = .false.
    ! Our cmd_arg type is already initialised to the sentinel
    do i = 1, num_args
      if(all_args(i)%name == trim(adjustl(name))) then
        found = .true.
        read(all_args(i)%value, *, iostat=ierr) val
        exit
      end if
    end do

    if(present(exists)) then
      exists = found
    end if

    ! Return value is whether value is found and correctly parsed
    get_arg_name_dbl = (found .and. (ierr == 0))

  end function get_arg_name_dbl

  !> @brief Read by number for long integer values
  !>
  !> @author  H. Ratcliffe
  !> @author  C. S. Brady
  !>
  !> @date    2019-2021
  !>
  !> @param num Argument number to read
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !>
  !> @return True if the argument is found and parsed, False otherwise
  function get_arg_num_int(num, val, exists)

    logical :: get_arg_num_int
    integer, intent(in) :: num
    integer(kind=int32), intent(out) :: val
    logical, intent(out), optional :: exists
    logical :: found
    integer :: ierr

    found = .false.
    ! Check requested number is in range
    if(num <= num_args .and. num > 0) then
      ! READ it from string into value
      ! We don't need to specify the format in general
      read(all_args(num)%value, *, iostat=ierr) val
      found = .true.
    end if

    if(present(exists)) then
      exists = found
    end if

    ! Return value is whether value is found and correctly parsed
    get_arg_num_int = (found .and. (ierr == 0))

  end function get_arg_num_int

  !> @brief Read by name for long integer values
  !>
  !> @author  H. Ratcliffe
  !> @author  C. S. Brady
  !>
  !> @date    2019-2021
  !>
  !> @param name Argument name to look up
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !>
  !> @return True if the name is found and parsed, False otherwise
  function get_arg_name_int(name, val, exists)

    logical :: get_arg_name_int
    character(len=*), intent(in) :: name
    integer(kind=int32), intent(out) :: val
    integer :: i
    logical, intent(out), optional :: exists
    logical :: found
    integer :: ierr

    found = .false.
    ! Our cmd_arg type is already initialised to the sentinel
    do i = 1, num_args
      if(all_args(i)%name == trim(adjustl(name))) then
        found = .true.
        read(all_args(i)%value, *, iostat=ierr) val
        exit
      end if
    end do

    if(present(exists)) then
      exists = found
    end if

    ! Return value is whether value is found and correctly parsed
    get_arg_name_int = (found .and. (ierr == 0))

  end function get_arg_name_int

  !> @brief Read by number for string/character values
  !>
  !> @author  H. Ratcliffe
  !> @author  C. S. Brady
  !>
  !> @date    2019-2021
  !>
  !> @param num Argument number to read
  !> @param val Value to read into
  !> @param exists Whether the name was found - this is already
  !>               contained in the return value, but is given
  !>               for consistency with the other members
  !>
  !> @return True if the argument is found and parsed, False otherwise
  function get_arg_num_str(num, val, exists)

    logical :: get_arg_num_str
    integer, intent(in) :: num
    character(len=*), intent(out) :: val
    logical, intent(out), optional :: exists
    logical :: found

    found = .false.
    ! Check requested number is in range
    if(num <= num_args .and. num > 0) then
      ! Read it from string into value
      ! We don't need to specify the format in general
      val = all_args(num)%value
      found = .true.
    end if

    if(present(exists)) then
      exists = found
    end if

    ! Return value is whether value is found and correctly parsed
    get_arg_num_str = found

  end function get_arg_num_str

  !> @brief Read by name for string values
  !>
  !> @author  H. Ratcliffe
  !> @author  C. S. Brady
  !>
  !> @date    2019-2021
  !>
  !> @param name Argument name to look up
  !> @param val Value to read into
  !> @param exists Whether the name was found - this is already
  !>               contained in the return value, but is given
  !>               for consistency with the other members
  !>
  !> @return True if the name is found and parsed, False otherwise
  function get_arg_name_str(name, val, exists)

    logical :: get_arg_name_str
    character(len=*), intent(in) :: name
    character(len=*), intent(out) :: val
    integer :: i
    logical, intent(out), optional :: exists
    logical :: found

    found = .false.
    ! Our cmd_arg type is already initialised to the sentinel
    do i = 1, num_args
      if(all_args(i)%name == trim(adjustl(name))) then
        found = .true.
        val = all_args(i)%value
        exit
      end if
    end do

    if(present(exists)) then
      exists = found
    end if

    ! Return value is whether value is found and correctly parsed
    get_arg_name_str = found

  end function get_arg_name_str

  !--------------------------------------------------------------------!
  ! This lets you get just the string value from an argument by name.  !
  ! You still need to have called parse_args first!                    !
  !--------------------------------------------------------------------!

  !> @brief Lookup an argument by name and return the value as a string
  !>
  !> @author  H. Ratcliffe
  !> @author  C. S. Brady
  !>
  !> @date    2019-2021
  !>
  !> @param name Argument name to look up
  !> @param exists Whether the name was found
  !>
  !> @return The string value associated with the given name
  function get_arg_value(name, exists)

    character(len=len) :: get_arg_value
    character(len=*), intent(in) :: name
    logical, intent(out), optional :: exists
    type(cmd_arg) :: tmp
    integer :: i
    logical :: found

    ! Initialise to the default value
    ! Use a temporary to get the default values
    ! This makes sure we match the expected sentinel
    get_arg_value = tmp%value
    found = .false.

    do i = 1, num_args
      if(all_args(i)%name .eq. trim(adjustl(name))) then
        get_arg_value = all_args(i)%value
        found = .true.
        exit
      end if
    end do

    if(present(exists)) then
      exists = found
    end if

  end function get_arg_value

end module command_line
