  !> @brief Functions to parse command line arguments
  !>
  !> Module to read command line arguments to a program
  !> We assume they are of the form name=value (spaces around '=' are ignored) or are a flag
  !> NOTE: 'val=""' differs from 'val' - the latter is a flag, the former an empty string
  !> Value can be extracted as a string, integer, a long-integer
  !> or a single or double-precision real, according to the
  !> type passed in.
  !> Argument names are limited to 20 chars, and values
  !> to 30 chars as read.
  !>
  !> Note that the only functions you should call from outside are
  !> get_arg and get_arg_value for 'key=value' arguments
  !> arg_present to check presence and state (flag or valued)
  !> arg_count and dump_names for general inquiries
  !> NOTE: flags can also be read by name as a logical: TRUE if present, FALSE if not
  !> NOTE: total count may not match COMMAND_ARGUMENT_COUNT due to
  !> parsing spaces out of key( )=()value syntax!

  !> @author H Ratcliffe, Senior Research Software Engineer, University of Warwick
  !> @date 1/11/24

MODULE command_line

  USE kinds

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_arg, get_arg_value, arg_present ! Main accessors
  PUBLIC :: arg_count, dump_names, str_wrapper  ! Helpers

  LOGICAL :: initial_parse_done = .FALSE.

  !> Type containing a key-value pair
  ! Init. to default values
  TYPE cmd_arg
    CHARACTER(LEN=:), ALLOCATABLE :: name
    CHARACTER(LEN=:), ALLOCATABLE :: value
    LOGICAL :: has_value = .TRUE. !Whether a value was supplied (c.f. empty)
  END TYPE

  !> Wrapper to allow array of allocatable strings
  TYPE str_wrapper
    CHARACTER(LEN=:), ALLOCATABLE :: str
  END TYPE

  !> @brief Read arguments by name or number
  !>
  !> Get arguments. If the first parameter is a string,
  !> this is interpreted as the name of the parameter, if
  !> an integer, it is the position of the argument in the input list.
  !> @param name Name to look up (supply this OR num)
  !> @param num Arg. number to look up (supply this OR name)
  !> @param val Value to read into, with type matching that to parse as
  !> @param exists Whether the name was found
  !> @return True if the name is found and parsed, False otherwise
  INTERFACE get_arg
    MODULE PROCEDURE get_arg_num_logical, get_arg_name_logical
    MODULE PROCEDURE get_arg_num_int, get_arg_name_int
    MODULE PROCEDURE get_arg_num_long, get_arg_name_long
    MODULE PROCEDURE get_arg_num_float, get_arg_name_float
    MODULE PROCEDURE get_arg_num_dbl, get_arg_name_dbl
    MODULE PROCEDURE get_arg_num_str, get_arg_name_str
  END INTERFACE

  !> The argument list
  TYPE(cmd_arg), DIMENSION(:), ALLOCATABLE :: all_args
  !> The number of arguments
  INTEGER :: num_args = 0
  INTEGER, PARAMETER :: max_string_len = 200
  PRIVATE :: all_args, num_args

  CONTAINS

  !> @brief Parse out command line args
  !>
  !> This function can be called multiple times
  !> and will freshly parse ALL arguments each time
  !> We assume these are entered as 'name=value'. If there is
  !> no '=' sign, we set an empty value
  SUBROUTINE parse_args()

    ! Strictly we can't be sure any max_string_len is enough
    ! but for command-line args it's enough if we're sensible
    ! We wont overflow, but our strings may get truncated
    ! We'll print an warning, since this is probably unintended input

    ! NOTE: some codes have reason to disable implicit re-allocation so we
    ! take the extra effort to allocate all our strings manually

    INTEGER :: i_arg, i_tok, indx
    TYPE(cmd_arg), DIMENSION(:), ALLOCATABLE :: all_args_tmp
    CHARACTER(LEN=max_string_len) :: arg, tmp, tmp_name, tmp_val
    INTEGER :: arg_in_length, tmp_len
    LOGICAL :: truncated

    truncated = .FALSE.

    num_args = COMMAND_ARGUMENT_COUNT()
    IF(num_args > 0) THEN

      ! If this is not the first call, all_args may be already allocated
      ! Deallocate if needed, and allocate to sufficient size
      ! Will be trimmed to actual size after parsing
      IF(ALLOCATED(all_args)) DEALLOCATE(all_args)
      ALLOCATE(all_args(num_args))

      i_arg = 1 !Index of current arg
      i_tok = 1 ! Index of current input token
      ! Loop over all arguments and extract
      DO WHILE (i_tok <= num_args)
        ! First extract name and value parts in all cases
        ! This consumes 1, 2 or 3 tokens depending on spaces

        CALL GET_COMMAND_ARGUMENT(i_tok, arg, length=arg_in_length)
        i_tok = i_tok + 1

        IF(arg_in_length > max_string_len) truncated = .TRUE.

        ! Location of the '=' sign
        ! If not found, return value of this is 0
        indx = INDEX(arg, '=')

        !Look at next chars - remove all whitespace
        tmp = ADJUSTL(arg(indx+1:))
        tmp_len = LEN_TRIM(tmp)
        IF(indx > 1 .AND. tmp_len > 0) THEN
          ! All characters after '='
          tmp_val = tmp
          ! All characters up to '=', not including it
          ! but with any leading spaces removed
          tmp_name = ADJUSTL(arg(1:indx-1))
        ELSE IF(indx > 1) THEN
          ! Have an '=' but no following value
          ! Consume next token
          CALL GET_COMMAND_ARGUMENT(i_tok, tmp, length=arg_in_length)
          i_tok = i_tok + 1
          IF(arg_in_length > max_string_len) truncated = .TRUE.

          tmp_val = ADJUSTL(tmp)
          tmp_name = ADJUSTL(arg(1:indx-1))
        ELSE   ! Have not yet found the equals!
          ! Set name, then hunt value...
          tmp_name = ADJUSTL(arg)

          !Peek next token - will need either 0, 1 or 2 more
          CALL GET_COMMAND_ARGUMENT(i_tok, tmp, length=arg_in_length)
          IF(arg_in_length > max_string_len) truncated = .TRUE.

          indx = INDEX(ADJUSTL(tmp), '=')
          IF(indx /= 1) THEN
            ! Next token does not lead with '=', assume this is a flag and
            ! DO NOT consume next. Set value for clarity, and mark
            tmp_val = ""
            all_args(i_arg)%has_value = .FALSE.
          ELSE
            ! Consume this one and possibly one more
            i_tok = i_tok + 1
            IF(LEN_TRIM(ADJUSTL(tmp)) > 1) THEN
              !This token has content following the '='
              tmp_val = ADJUSTL(tmp(2:))
            ELSE
              ! Consume another
              CALL GET_COMMAND_ARGUMENT(i_tok, tmp, length=arg_in_length)
              i_tok = i_tok + 1
              IF(arg_in_length > max_string_len) truncated = .TRUE.

              tmp_val = ADJUSTL(tmp)
            END IF
          END IF
        END IF

        ! Explicitly allocate and set the values
        ALLOCATE(CHARACTER(len=LEN_TRIM(tmp_name)) :: all_args(i_arg)%name)
        all_args(i_arg)%name = TRIM(tmp_name)
        ALLOCATE(CHARACTER(len=LEN_TRIM(tmp_val)) :: all_args(i_arg)%value)
        all_args(i_arg)%value = TRIM(tmp_val)

        i_arg = i_arg + 1
      END DO
      !i_arg is now the actual parsed count
      !Shrink array to get rid of excess unfilled space
      num_args = i_arg-1
      CALL MOVE_ALLOC(all_args, all_args_tmp)
      ALLOCATE(all_args(num_args))
      all_args = all_args_tmp(1:num_args)
      DEALLOCATE(all_args_tmp)
    ENDIF

    IF(truncated) PRINT'(A,I0, A)', "WARNING: Very long argument truncated. To support arguments&
   & longer than ", max_string_len, " increase the max_string_len parameter"


  END SUBROUTINE parse_args

  !> Helper function - do a parse if it hasn't been done yet
  SUBROUTINE initial_parse

    IF(.NOT. initial_parse_done) THEN
      CALL parse_args
      initial_parse_done = .TRUE.
    END IF

  END SUBROUTINE initial_parse

  !> Get the number of arguments
  !> NOTE: total count may not match COMMAND_ARGUMENT_COUNT due to
  !> parsing key=value syntax!
  FUNCTION arg_count()
    INTEGER :: arg_count

    CALL initial_parse
    arg_count = num_args
  END FUNCTION

!------------------------------------------------------------------

  !> @brief Read by number for logical values
  !> @param num Argument number to read
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !> @return True if the name is found and parsed, False otherwise
  FUNCTION get_arg_num_logical(num, val, exists)

    LOGICAL :: get_arg_num_logical
    INTEGER, INTENT(IN) :: num
    LOGICAL, INTENT(INOUT) :: val
    LOGICAL, INTENT(OUT), OPTIONAL :: exists
    LOGICAL :: found
    INTEGER :: ierr

    CALL initial_parse

    found = .FALSE.
    ! Check requested number is in range
    IF(num <= num_args .AND. num > 0) THEN
      ! READ it from string into value
      READ(all_args(num)%value, *, IOSTAT=ierr) val
      found = .TRUE.
    END IF

    IF(PRESENT(exists)) THEN
      exists = found
    END IF

    ! Return value is whether value is found and correctly parsed
    get_arg_num_logical = (found .AND. (ierr == 0))

  END FUNCTION get_arg_num_logical

  !> @brief Read by name for logical values
  !> NOTE : a flag (name with no '=value' part) will parse as
  !> .TRUE. if present, and .FALSE. if not
  !> @param name Argument name to look up
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !> @return True if the name is found and parsed, False otherwise
  FUNCTION get_arg_name_logical(name, val, exists)

    LOGICAL :: get_arg_name_logical
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL, INTENT(INOUT) :: val
    INTEGER :: i
    LOGICAL, INTENT(OUT), OPTIONAL :: exists
    LOGICAL :: found
    INTEGER :: ierr

    CALL initial_parse

    found = .FALSE.
    val = .FALSE.
    ierr = 0
    ! Our cmd_arg type is already initialised to the sentinel
    DO i = 1, num_args
      IF(all_args(i)%name == TRIM(ADJUSTL(name))) THEN
        found = .TRUE.
        IF(all_args(i)%has_value) THEN
          READ(all_args(i)%value, *, IOSTAT=ierr) val
        ELSE
          val = .TRUE.
        END IF
        EXIT
      END IF
    END DO

    IF(PRESENT(exists)) THEN
      exists = found
    END IF

    ! Return value is whether value is found and correctly parsed
    get_arg_name_logical = (found .AND. (ierr == 0))

  END FUNCTION get_arg_name_logical


  !> @brief Read by number for double precision values
  !> @param num Argument number to read
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !> @return True if the name is found and parsed, False otherwise
  FUNCTION get_arg_num_dbl(num, val, exists)

    LOGICAL :: get_arg_num_dbl
    INTEGER, INTENT(IN) :: num
    REAL(KIND=REAL64), INTENT(INOUT) :: val
    LOGICAL, INTENT(OUT), OPTIONAL :: exists
    LOGICAL :: found
    INTEGER :: ierr

    CALL initial_parse

    found = .FALSE.
    ! Check requested number is in range
    IF(num <= num_args .AND. num > 0) THEN
      ! READ it from string into value
      READ(all_args(num)%value, *, IOSTAT=ierr) val
      found = .TRUE.
    END IF

    IF(PRESENT(exists)) THEN
      exists = found
    END IF

    ! Return value is whether value is found and correctly parsed
    get_arg_num_dbl = (found .AND. (ierr == 0))

  END FUNCTION get_arg_num_dbl

  !> @brief Read by name for double precision values
  !> @param name Argument name to look up
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !> @return True if the name is found and parsed, False otherwise
  FUNCTION get_arg_name_dbl(name, val, exists)

    LOGICAL :: get_arg_name_dbl
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(KIND=REAL64), INTENT(INOUT) :: val
    INTEGER :: i
    LOGICAL, INTENT(OUT), OPTIONAL :: exists
    LOGICAL :: found
    INTEGER :: ierr

    CALL initial_parse

    found = .FALSE.
    ! Our cmd_arg type is already initialised to the sentinel
    DO i = 1, num_args
      IF(all_args(i)%name == TRIM(ADJUSTL(name))) THEN
        found = .TRUE.
        READ(all_args(i)%value, *, IOSTAT=ierr) val
        EXIT
      END IF
    END DO

    IF(PRESENT(exists)) THEN
      exists = found
    END IF

    ! Return value is whether value is found and correctly parsed
    get_arg_name_dbl = (found .AND. (ierr == 0))

  END FUNCTION get_arg_name_dbl

  ! Command line parsing should be avoided in performance critical code
  ! so extra overhead from double call and downcast is not a problem

! \TODO Use this approach, or just dupe. all the code?

  !> @brief Read by number for single precision (float) values
  !> @param num Argument number to read
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !> @return True if the name is found and parsed, False otherwise
  FUNCTION get_arg_num_float(num, val, exists)
    LOGICAL :: get_arg_num_float
    INTEGER, INTENT(IN) :: num
    REAL(KIND=REAL32), INTENT(INOUT) :: val
    REAL(KIND=REAL64) :: tmp
    LOGICAL, INTENT(OUT), OPTIONAL :: exists

    get_arg_num_float = get_arg_num_dbl(num, tmp, exists)
    IF( ABS(tmp) < HUGE(val)) THEN
      !Value in range. Convert. Note: there may be precision loss
      val = REAL(tmp, KIND=REAL32)
    ELSE
      !Value out of range, can't be parsed
      get_arg_num_float  = .FALSE.
    END IF

  END FUNCTION

  !> @brief Read by name for single precision (float) values
  !> @param name Argument name to look up
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !> @return True if the name is found and parsed, False otherwise
  FUNCTION get_arg_name_float(name, val, exists)
    LOGICAL :: get_arg_name_float
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(KIND=REAL32), INTENT(INOUT) :: val
    REAL(KIND=REAL64) :: tmp
    LOGICAL, INTENT(OUT), OPTIONAL :: exists

    get_arg_name_float = get_arg_name_dbl(name, tmp, exists)
    IF( ABS(tmp) < HUGE(val)) THEN
      !Value in range. Convert. Note: there may be precision loss
      val = REAL(tmp, KIND=REAL32)
    ELSE
      !Value out of range, can't be parsed
      get_arg_name_float  = .FALSE.
    END IF

  END FUNCTION


  !> @brief Read by number for integer values
  !> @param num Argument number to read
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !> @return True if the name is found and parsed, False otherwise
  FUNCTION get_arg_num_int(num, val, exists)

    LOGICAL :: get_arg_num_int
    INTEGER, INTENT(IN) :: num
    INTEGER(KIND=INT32), INTENT(INOUT) :: val
    LOGICAL, INTENT(OUT), OPTIONAL :: exists
    LOGICAL :: found
    INTEGER :: ierr

    CALL initial_parse

    found = .FALSE.
    ! Check requested number is in range
    IF(num <= num_args .AND. num > 0) THEN
      ! READ it from string into value
      ! We don't need to specify the format in general
      READ(all_args(num)%value, *, IOSTAT=ierr) val
      found = .TRUE.
    END IF

    IF(PRESENT(exists)) THEN
      exists = found
    END IF

    ! Return value is whether value is found and correctly parsed
    get_arg_num_int = (found .AND. (ierr == 0))

  END FUNCTION get_arg_num_int

  !> @brief Read by name for integer values
  !> @param name Argument name to look up
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !> @return True if the name is found and parsed, False otherwise
  FUNCTION get_arg_name_int(name, val, exists)

    LOGICAL :: get_arg_name_int
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER(KIND=INT32), INTENT(INOUT) :: val
    INTEGER :: i
    LOGICAL, INTENT(OUT), OPTIONAL :: exists
    LOGICAL :: found
    INTEGER :: ierr

    CALL initial_parse

    found = .FALSE.
    ! Our cmd_arg type is already initialised to the sentinel
    DO i = 1, num_args
      IF(all_args(i)%name == TRIM(ADJUSTL(name))) THEN
        found = .TRUE.
        READ(all_args(i)%value, *, IOSTAT=ierr) val
        EXIT
      END IF
    END DO

    IF(PRESENT(exists)) THEN
      exists = found
    END IF

    ! Return value is whether value is found and correctly parsed
    get_arg_name_int = (found .AND. (ierr == 0))

  END FUNCTION get_arg_name_int

  !> @brief Read by number for long integer values
  !> @param num Argument number to read
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !> @return True if the name is found and parsed, False otherwise
  FUNCTION get_arg_num_long(num, val, exists)

    LOGICAL :: get_arg_num_long
    INTEGER, INTENT(IN) :: num
    INTEGER(KIND=INT64), INTENT(INOUT) :: val
    LOGICAL, INTENT(OUT), OPTIONAL :: exists
    LOGICAL :: found
    INTEGER :: ierr

    CALL initial_parse

    found = .FALSE.
    ! Check requested number is in range
    IF(num <= num_args .AND. num > 0) THEN
      ! READ it from string into value
      ! We don't need to specify the format in general
      READ(all_args(num)%value, *, IOSTAT=ierr) val
      found = .TRUE.
    END IF

    IF(PRESENT(exists)) THEN
      exists = found
    END IF

    ! Return value is whether value is found and correctly parsed
    get_arg_num_long = (found .AND. (ierr == 0))

  END FUNCTION get_arg_num_long

  !> @brief Read by name for long integer values
  !> @param name Argument name to look up
  !> @param val Value to read into
  !> @param exists Whether the name was found
  !> @return True if the name is found and parsed, False otherwise
  FUNCTION get_arg_name_long(name, val, exists)

    LOGICAL :: get_arg_name_long
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER(KIND=INT64), INTENT(INOUT) :: val
    INTEGER :: i
    LOGICAL, INTENT(OUT), OPTIONAL :: exists
    LOGICAL :: found
    INTEGER :: ierr

   CALL initial_parse

    found = .FALSE.
    ! Our cmd_arg type is already initialised to the sentinel
    DO i = 1, num_args
      IF(all_args(i)%name == TRIM(ADJUSTL(name))) THEN
        found = .TRUE.
        READ(all_args(i)%value, *, IOSTAT=ierr) val
        EXIT
      END IF
    END DO

    IF(PRESENT(exists)) THEN
      exists = found
    END IF

    ! Return value is whether value is found and correctly parsed
    get_arg_name_long = (found .AND. (ierr == 0))

  END FUNCTION get_arg_name_long

  !> @brief Read by number for string/character values
  !> @param num Argument number to read
  !> @param val Value to read into
  !> @param exists Whether the name was found - this is already contained in
  !> the return value, but is given for consistency with the other members
  !> @return True if the name is found and parsed, False otherwise
  FUNCTION get_arg_num_str(num, val, exists)

    LOGICAL :: get_arg_num_str
    INTEGER, INTENT(IN) :: num
    CHARACTER(LEN=*), INTENT(INOUT) :: val
    LOGICAL, INTENT(OUT), OPTIONAL :: exists
    LOGICAL :: found

    CALL initial_parse

    found = .FALSE.
    ! Check requested number is in range
    IF(num <= num_args .AND. num > 0) THEN
      ! READ it from string into value
      ! We don't need to specify the format in general
      val = all_args(num)%value
      found = .TRUE.
    END IF

    IF(PRESENT(exists)) THEN
      exists = found
    END IF

    ! Return value is whether value is found and correctly parsed
    get_arg_num_str = found

  END FUNCTION get_arg_num_str

  !> @brief Read by name for string values
  !> NOTE: if the string passed is shorter than the value, it will be truncated
  !> If the length is not known use get_arg_value to get an allocatable string
  !> @param name Argument name to look up
  !> @param val Value to read into
  !> @param exists Whether the name was found - this is already contained in
  !> the return value, but is given for consistency with the other members
  !> @return True if the name is found and parsed, False otherwise
  FUNCTION get_arg_name_str(name, val, exists)

    LOGICAL :: get_arg_name_str
    CHARACTER(LEN=*), INTENT(IN) :: name
    CHARACTER(LEN=*), INTENT(INOUT) :: val
    INTEGER :: i
    LOGICAL, INTENT(OUT), OPTIONAL :: exists
    LOGICAL :: found

    CALL initial_parse

    found = .FALSE.
    DO i = 1, num_args
      IF(all_args(i)%name == TRIM(ADJUSTL(name))) THEN
        found = .TRUE.
        val = all_args(i)%value
        EXIT
      END IF
    END DO

    IF(PRESENT(exists)) THEN
      exists = found
    END IF

    ! Return value is whether value is found and correctly parsed
    get_arg_name_str = found

  END FUNCTION get_arg_name_str


!--------------------------------------------------------------------
  !> @brief Check presence of an argument by name
  !> @param name Argument name to look up
  !> @param has_value Whether this argument has a defined value (also .FALSE. if not present)
  !> @return True if present, False if not
  FUNCTION arg_present(name, has_value) RESULT(found)

   LOGICAL :: found
   CHARACTER(LEN=*), INTENT(IN) :: name
   LOGICAL, INTENT(OUT), OPTIONAL :: has_value
   INTEGER :: i

    CALL initial_parse

    found = .FALSE.
    IF(PRESENT(has_value)) has_value = .FALSE.
    DO i = 1, num_args
      IF(all_args(i)%name == TRIM(ADJUSTL(name))) THEN
        found = .TRUE.
        IF(PRESENT(has_value)) has_value = all_args(i)%has_value
        EXIT
      ENDIF
    END DO

  END FUNCTION arg_present

  !> @brief Get all the argument names (by copy)
  !> Order will _probably_ match input order, but this is not guaranteed
  !> @return Array of str_wrapper types containing all argument names
  FUNCTION dump_names()
    TYPE(str_wrapper), DIMENSION(:), ALLOCATABLE :: dump_names
    INTEGER :: i

    ! Do this the long way again to support non implicit allocations
    ALLOCATE(dump_names(num_args))
    DO i = 1, num_args
      ALLOCATE(dump_names(i)%str, SOURCE =all_args(i)%name)
    END DO

  END FUNCTION dump_names

  !> @brief Lookup an argument by name and return the value as an (allocatable) string
  !> If the name is NOT PRESENT, an empty string is returned
  !> @param name Argument name to look up
  !> @param exists Whether the name was found
  !> @return The string value associated with the given name
  FUNCTION get_arg_value(name, exists)

    CHARACTER(LEN=:), ALLOCATABLE :: get_arg_value
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL, INTENT(OUT), OPTIONAL :: exists
    INTEGER :: i
    LOGICAL :: found

    CALL initial_parse

    found = .FALSE.

    DO i = 1, num_args
      IF(all_args(i)%name .EQ. TRIM(ADJUSTL(name))) THEN
        ALLOCATE(get_arg_value, SOURCE=all_args(i)%value)
        found = .TRUE.
        EXIT
      END IF
    END DO

    ! Return empty string, not unallocated one.
    IF(.NOT. ALLOCATED(get_arg_value)) ALLOCATE(CHARACTER(LEN=0) :: get_arg_value)
    IF(PRESENT(exists)) THEN
      exists = found
    END IF

  END FUNCTION get_arg_value


END MODULE command_line
