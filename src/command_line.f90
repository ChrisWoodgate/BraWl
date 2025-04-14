  !> @brief functions to parse command line arguments
  !>
  !> module to read command line arguments to a program
  !> we assume they are of the form name=value (spaces around '=' are ignored) or are a flag
  !> note: 'val=""' differs from 'val' - the latter is a flag, the former an empty string
  !> value can be extracted as a string, integer, a long-integer
  !> or a single or double-precision real, according to the
  !> type passed in.
  !> argument names are limited to 20 chars, and values
  !> to 30 chars as read.
  !>
  !> note that the only functions you should call from outside are
  !> get_arg and get_arg_value for 'key=value' arguments
  !> arg_present to check presence and state (flag or valued)
  !> arg_count and dump_names for general inquiries
  !> note: flags can also be read by name as a logical: true if present, false if not
  !> note: total count may not match command_argument_count due to
  !> parsing spaces out of key( )=()value syntax!

  !> @author h ratcliffe, senior research software engineer, university of warwick
  !> @date 1/11/24

module command_line

  use kinds

  implicit none

  private
  public :: get_arg, get_arg_value, arg_present ! main accessors
  public :: arg_count, dump_names, str_wrapper  ! helpers

  logical :: initial_parse_done = .false.

  !> type containing a key-value pair
  ! init. to default values
  type cmd_arg
    character(len=:), allocatable :: name
    character(len=:), allocatable :: value
    logical :: has_value = .true. !whether a value was supplied (c.f. empty)
  end type

  !> wrapper to allow array of allocatable strings
  type str_wrapper
    character(len=:), allocatable :: str
  end type

  !> @brief read arguments by name or number
  !>
  !> get arguments. if the first parameter is a string,
  !> this is interpreted as the name of the parameter, if
  !> an integer, it is the position of the argument in the input list.
  !> @param name name to look up (supply this or num)
  !> @param num arg. number to look up (supply this or name)
  !> @param val value to read into, with type matching that to parse as
  !> @param exists whether the name was found
  !> @return true if the name is found and parsed, false otherwise
  interface get_arg
    module procedure get_arg_num_logical, get_arg_name_logical
    module procedure get_arg_num_int, get_arg_name_int
    module procedure get_arg_num_long, get_arg_name_long
    module procedure get_arg_num_float, get_arg_name_float
    module procedure get_arg_num_dbl, get_arg_name_dbl
    module procedure get_arg_num_str, get_arg_name_str
  end interface

  !> the argument list
  type(cmd_arg), dimension(:), allocatable :: all_args
  !> the number of arguments
  integer :: num_args = 0
  integer, parameter :: max_string_len = 200
  private :: all_args, num_args

  contains

  !> @brief parse out command line args
  !>
  !> this function can be called multiple times
  !> and will freshly parse all arguments each time
  !> we assume these are entered as 'name=value'. if there is
  !> no '=' sign, we set an empty value
  subroutine parse_args()

    ! strictly we can't be sure any max_string_len is enough
    ! but for command-line args it's enough if we're sensible
    ! we wont overflow, but our strings may get truncated
    ! we'll print an warning, since this is probably unintended input

    ! note: some codes have reason to disable implicit re-allocation so we
    ! take the extra effort to allocate all our strings manually

    integer :: i_arg, i_tok, indx
    type(cmd_arg), dimension(:), allocatable :: all_args_tmp
    character(len=max_string_len) :: arg, tmp, tmp_name, tmp_val
    integer :: arg_in_length, tmp_len
    logical :: truncated

    truncated = .false.

    num_args = command_argument_count()
    if(num_args > 0) then

      ! if this is not the first call, all_args may be already allocated
      ! deallocate if needed, and allocate to sufficient size
      ! will be trimmed to actual size after parsing
      if(allocated(all_args)) deallocate(all_args)
      allocate(all_args(num_args))

      i_arg = 1 !index of current arg
      i_tok = 1 ! index of current input token
      ! loop over all arguments and extract
      do while (i_tok <= num_args)
        ! first extract name and value parts in all cases
        ! this consumes 1, 2 or 3 tokens depending on spaces

        call get_command_argument(i_tok, arg, length=arg_in_length)
        i_tok = i_tok + 1

        if(arg_in_length > max_string_len) truncated = .true.

        ! location of the '=' sign
        ! if not found, return value of this is 0
        indx = index(arg, '=')

        !look at next chars - remove all whitespace
        tmp = adjustl(arg(indx+1:))
        tmp_len = len_trim(tmp)
        if(indx > 1 .and. tmp_len > 0) then
          ! all characters after '='
          tmp_val = tmp
          ! all characters up to '=', not including it
          ! but with any leading spaces removed
          tmp_name = adjustl(arg(1:indx-1))
        else if(indx > 1) then
          ! have an '=' but no following value
          ! consume next token
          call get_command_argument(i_tok, tmp, length=arg_in_length)
          i_tok = i_tok + 1
          if(arg_in_length > max_string_len) truncated = .true.

          tmp_val = adjustl(tmp)
          tmp_name = adjustl(arg(1:indx-1))
        else   ! have not yet found the equals!
          ! set name, then hunt value...
          tmp_name = adjustl(arg)

          !peek next token - will need either 0, 1 or 2 more
          call get_command_argument(i_tok, tmp, length=arg_in_length)
          if(arg_in_length > max_string_len) truncated = .true.

          indx = index(adjustl(tmp), '=')
          if(indx /= 1) then
            ! next token does not lead with '=', assume this is a flag and
            ! do not consume next. set value for clarity, and mark
            tmp_val = ""
            all_args(i_arg)%has_value = .false.
          else
            ! consume this one and possibly one more
            i_tok = i_tok + 1
            if(len_trim(adjustl(tmp)) > 1) then
              !this token has content following the '='
              tmp_val = adjustl(tmp(2:))
            else
              ! consume another
              call get_command_argument(i_tok, tmp, length=arg_in_length)
              i_tok = i_tok + 1
              if(arg_in_length > max_string_len) truncated = .true.

              tmp_val = adjustl(tmp)
            end if
          end if
        end if

        ! explicitly allocate and set the values
        allocate(character(len=len_trim(tmp_name)) :: all_args(i_arg)%name)
        all_args(i_arg)%name = trim(tmp_name)
        allocate(character(len=len_trim(tmp_val)) :: all_args(i_arg)%value)
        all_args(i_arg)%value = trim(tmp_val)

        i_arg = i_arg + 1
      end do
      !i_arg is now the actual parsed count
      !shrink array to get rid of excess unfilled space
      num_args = i_arg-1
      call move_alloc(all_args, all_args_tmp)
      allocate(all_args(num_args))
      all_args = all_args_tmp(1:num_args)
      deallocate(all_args_tmp)
    endif

    if(truncated) print'(a,i0, a)', "warning: very long argument truncated. to support arguments&
   & longer than ", max_string_len, " increase the max_string_len parameter"


  end subroutine parse_args

  !> helper function - do a parse if it hasn't been done yet
  subroutine initial_parse

    if(.not. initial_parse_done) then
      call parse_args
      initial_parse_done = .true.
    end if

  end subroutine initial_parse

  !> get the number of arguments
  !> note: total count may not match command_argument_count due to
  !> parsing key=value syntax!
  function arg_count()
    integer :: arg_count

    call initial_parse
    arg_count = num_args
  end function

!------------------------------------------------------------------

  !> @brief read by number for logical values
  !> @param num argument number to read
  !> @param val value to read into
  !> @param exists whether the name was found
  !> @return true if the name is found and parsed, false otherwise
  function get_arg_num_logical(num, val, exists)

    logical :: get_arg_num_logical
    integer, intent(in) :: num
    logical, intent(inout) :: val
    logical, intent(out), optional :: exists
    logical :: found
    integer :: ierr

    call initial_parse

    found = .false.
    ! check requested number is in range
    if(num <= num_args .and. num > 0) then
      ! read it from string into value
      read(all_args(num)%value, *, iostat=ierr) val
      found = .true.
    end if

    if(present(exists)) then
      exists = found
    end if

    ! return value is whether value is found and correctly parsed
    get_arg_num_logical = (found .and. (ierr == 0))

  end function get_arg_num_logical

  !> @brief read by name for logical values
  !> note : a flag (name with no '=value' part) will parse as
  !> .true. if present, and .false. if not
  !> @param name argument name to look up
  !> @param val value to read into
  !> @param exists whether the name was found
  !> @return true if the name is found and parsed, false otherwise
  function get_arg_name_logical(name, val, exists)

    logical :: get_arg_name_logical
    character(len=*), intent(in) :: name
    logical, intent(inout) :: val
    integer :: i
    logical, intent(out), optional :: exists
    logical :: found
    integer :: ierr

    call initial_parse

    found = .false.
    val = .false.
    ierr = 0
    ! our cmd_arg type is already initialised to the sentinel
    do i = 1, num_args
      if(all_args(i)%name == trim(adjustl(name))) then
        found = .true.
        if(all_args(i)%has_value) then
          read(all_args(i)%value, *, iostat=ierr) val
        else
          val = .true.
        end if
        exit
      end if
    end do

    if(present(exists)) then
      exists = found
    end if

    ! return value is whether value is found and correctly parsed
    get_arg_name_logical = (found .and. (ierr == 0))

  end function get_arg_name_logical


  !> @brief read by number for double precision values
  !> @param num argument number to read
  !> @param val value to read into
  !> @param exists whether the name was found
  !> @return true if the name is found and parsed, false otherwise
  function get_arg_num_dbl(num, val, exists)

    logical :: get_arg_num_dbl
    integer, intent(in) :: num
    real(kind=real64), intent(inout) :: val
    logical, intent(out), optional :: exists
    logical :: found
    integer :: ierr

    call initial_parse

    found = .false.
    ! check requested number is in range
    if(num <= num_args .and. num > 0) then
      ! read it from string into value
      read(all_args(num)%value, *, iostat=ierr) val
      found = .true.
    end if

    if(present(exists)) then
      exists = found
    end if

    ! return value is whether value is found and correctly parsed
    get_arg_num_dbl = (found .and. (ierr == 0))

  end function get_arg_num_dbl

  !> @brief read by name for double precision values
  !> @param name argument name to look up
  !> @param val value to read into
  !> @param exists whether the name was found
  !> @return true if the name is found and parsed, false otherwise
  function get_arg_name_dbl(name, val, exists)

    logical :: get_arg_name_dbl
    character(len=*), intent(in) :: name
    real(kind=real64), intent(inout) :: val
    integer :: i
    logical, intent(out), optional :: exists
    logical :: found
    integer :: ierr

    call initial_parse

    found = .false.
    ! our cmd_arg type is already initialised to the sentinel
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

    ! return value is whether value is found and correctly parsed
    get_arg_name_dbl = (found .and. (ierr == 0))

  end function get_arg_name_dbl

  ! command line parsing should be avoided in performance critical code
  ! so extra overhead from double call and downcast is not a problem


  !> @brief read by number for single precision (float) values
  !> @param num argument number to read
  !> @param val value to read into
  !> @param exists whether the name was found
  !> @return true if the name is found and parsed, false otherwise
  function get_arg_num_float(num, val, exists)
    logical :: get_arg_num_float
    integer, intent(in) :: num
    real(kind=real32), intent(inout) :: val
    real(kind=real64) :: tmp
    logical, intent(out), optional :: exists

    get_arg_num_float = get_arg_num_dbl(num, tmp, exists)
    if( abs(tmp) < huge(val)) then
      !value in range. convert. note: there may be precision loss
      val = real(tmp, kind=real32)
    else
      !value out of range, can't be parsed
      get_arg_num_float  = .false.
    end if

  end function

  !> @brief read by name for single precision (float) values
  !> @param name argument name to look up
  !> @param val value to read into
  !> @param exists whether the name was found
  !> @return true if the name is found and parsed, false otherwise
  function get_arg_name_float(name, val, exists)
    logical :: get_arg_name_float
    character(len=*), intent(in) :: name
    real(kind=real32), intent(inout) :: val
    real(kind=real64) :: tmp
    logical, intent(out), optional :: exists

    get_arg_name_float = get_arg_name_dbl(name, tmp, exists)
    if( abs(tmp) < huge(val)) then
      !value in range. convert. note: there may be precision loss
      val = real(tmp, kind=real32)
    else
      !value out of range, can't be parsed
      get_arg_name_float  = .false.
    end if

  end function


  !> @brief read by number for integer values
  !> @param num argument number to read
  !> @param val value to read into
  !> @param exists whether the name was found
  !> @return true if the name is found and parsed, false otherwise
  function get_arg_num_int(num, val, exists)

    logical :: get_arg_num_int
    integer, intent(in) :: num
    integer(kind=int32), intent(inout) :: val
    logical, intent(out), optional :: exists
    logical :: found
    integer :: ierr

    call initial_parse

    found = .false.
    ! check requested number is in range
    if(num <= num_args .and. num > 0) then
      ! read it from string into value
      ! we don't need to specify the format in general
      read(all_args(num)%value, *, iostat=ierr) val
      found = .true.
    end if

    if(present(exists)) then
      exists = found
    end if

    ! return value is whether value is found and correctly parsed
    get_arg_num_int = (found .and. (ierr == 0))

  end function get_arg_num_int

  !> @brief read by name for integer values
  !> @param name argument name to look up
  !> @param val value to read into
  !> @param exists whether the name was found
  !> @return true if the name is found and parsed, false otherwise
  function get_arg_name_int(name, val, exists)

    logical :: get_arg_name_int
    character(len=*), intent(in) :: name
    integer(kind=int32), intent(inout) :: val
    integer :: i
    logical, intent(out), optional :: exists
    logical :: found
    integer :: ierr

    call initial_parse

    found = .false.
    ! our cmd_arg type is already initialised to the sentinel
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

    ! return value is whether value is found and correctly parsed
    get_arg_name_int = (found .and. (ierr == 0))

  end function get_arg_name_int

  !> @brief read by number for long integer values
  !> @param num argument number to read
  !> @param val value to read into
  !> @param exists whether the name was found
  !> @return true if the name is found and parsed, false otherwise
  function get_arg_num_long(num, val, exists)

    logical :: get_arg_num_long
    integer, intent(in) :: num
    integer(kind=int64), intent(inout) :: val
    logical, intent(out), optional :: exists
    logical :: found
    integer :: ierr

    call initial_parse

    found = .false.
    ! check requested number is in range
    if(num <= num_args .and. num > 0) then
      ! read it from string into value
      ! we don't need to specify the format in general
      read(all_args(num)%value, *, iostat=ierr) val
      found = .true.
    end if

    if(present(exists)) then
      exists = found
    end if

    ! return value is whether value is found and correctly parsed
    get_arg_num_long = (found .and. (ierr == 0))

  end function get_arg_num_long

  !> @brief read by name for long integer values
  !> @param name argument name to look up
  !> @param val value to read into
  !> @param exists whether the name was found
  !> @return true if the name is found and parsed, false otherwise
  function get_arg_name_long(name, val, exists)

    logical :: get_arg_name_long
    character(len=*), intent(in) :: name
    integer(kind=int64), intent(inout) :: val
    integer :: i
    logical, intent(out), optional :: exists
    logical :: found
    integer :: ierr

   call initial_parse

    found = .false.
    ! our cmd_arg type is already initialised to the sentinel
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

    ! return value is whether value is found and correctly parsed
    get_arg_name_long = (found .and. (ierr == 0))

  end function get_arg_name_long

  !> @brief read by number for string/character values
  !> @param num argument number to read
  !> @param val value to read into
  !> @param exists whether the name was found - this is already contained in
  !> the return value, but is given for consistency with the other members
  !> @return true if the name is found and parsed, false otherwise
  function get_arg_num_str(num, val, exists)

    logical :: get_arg_num_str
    integer, intent(in) :: num
    character(len=*), intent(inout) :: val
    logical, intent(out), optional :: exists
    logical :: found

    call initial_parse

    found = .false.
    ! check requested number is in range
    if(num <= num_args .and. num > 0) then
      ! read it from string into value
      ! we don't need to specify the format in general
      val = all_args(num)%value
      found = .true.
    end if

    if(present(exists)) then
      exists = found
    end if

    ! return value is whether value is found and correctly parsed
    get_arg_num_str = found

  end function get_arg_num_str

  !> @brief read by name for string values
  !> note: if the string passed is shorter than the value, it will be truncated
  !> if the length is not known use get_arg_value to get an allocatable string
  !> @param name argument name to look up
  !> @param val value to read into
  !> @param exists whether the name was found - this is already contained in
  !> the return value, but is given for consistency with the other members
  !> @return true if the name is found and parsed, false otherwise
  function get_arg_name_str(name, val, exists)

    logical :: get_arg_name_str
    character(len=*), intent(in) :: name
    character(len=*), intent(inout) :: val
    integer :: i
    logical, intent(out), optional :: exists
    logical :: found

    call initial_parse

    found = .false.
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

    ! return value is whether value is found and correctly parsed
    get_arg_name_str = found

  end function get_arg_name_str


!--------------------------------------------------------------------
  !> @brief check presence of an argument by name
  !> @param name argument name to look up
  !> @param has_value whether this argument has a defined value (also .false. if not present)
  !> @return true if present, false if not
  function arg_present(name, has_value) result(found)

   logical :: found
   character(len=*), intent(in) :: name
   logical, intent(out), optional :: has_value
   integer :: i

    call initial_parse

    found = .false.
    if(present(has_value)) has_value = .false.
    do i = 1, num_args
      if(all_args(i)%name == trim(adjustl(name))) then
        found = .true.
        if(present(has_value)) has_value = all_args(i)%has_value
        exit
      endif
    end do

  end function arg_present

  !> @brief get all the argument names (by copy)
  !> order will _probably_ match input order, but this is not guaranteed
  !> @return array of str_wrapper types containing all argument names
  function dump_names()
    type(str_wrapper), dimension(:), allocatable :: dump_names
    integer :: i

    ! do this the long way again to support non implicit allocations
    allocate(dump_names(num_args))
    do i = 1, num_args
      allocate(dump_names(i)%str, source =all_args(i)%name)
    end do

  end function dump_names

  !> @brief lookup an argument by name and return the value as an (allocatable) string
  !> if the name is not present, an empty string is returned
  !> @param name argument name to look up
  !> @param exists whether the name was found
  !> @return the string value associated with the given name
  function get_arg_value(name, exists)

    character(len=:), allocatable :: get_arg_value
    character(len=*), intent(in) :: name
    logical, intent(out), optional :: exists
    integer :: i
    logical :: found

    call initial_parse

    found = .false.

    do i = 1, num_args
      if(all_args(i)%name .eq. trim(adjustl(name))) then
        allocate(get_arg_value, source=all_args(i)%value)
        found = .true.
        exit
      end if
    end do

    ! return empty string, not unallocated one.
    if(.not. allocated(get_arg_value)) allocate(character(len=0) :: get_arg_value)
    if(present(exists)) then
      exists = found
    end if

  end function get_arg_value


end module command_line
