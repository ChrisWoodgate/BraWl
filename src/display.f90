!> @file    display.f90
!>
!> @brief   Assorted routines and tools for printing simulation
!>          information to the screen.
!>
!> @author  C. D. Woodgate
!> @date    2020-2023
module display

  use kinds
  use shared_data
  use constants
  use derived_types
  
  implicit none

  private

  public :: pretty_print_exchange, display_grid, print_centered_message

  contains

  !> @brief   Subroutine to print atom-atom effective pair interactions
  !>          to the screen in a human-readable format.
  !>
  !> @details User has choice of units to use. Default is meV, but also
  !>          support mRy.
  !>
  !> @author  C. D. Woodgate
  !> @date    2020-2023
  !>
  !> @param  setup Derived type containing simulation parameters
  !> @param  units Optional argument specifying the units to use
  !>
  !> @return None
  subroutine pretty_print_exchange(setup, units)
    integer :: i,j,k
    character(len=*), optional :: units
    character(len=10) :: aunits
    type(run_params), intent(in) :: setup

    ! Can do mRy or meV, but default is meV
    if (present(units)) then
      aunits = units
    else
      aunits = 'meV'
    end if

    ! Case mRy
    if (trim(aunits) .eq. 'mRy') then
      write(*,"(x,'V_ij to be used in calculation (mRy):',/)") 
      do i=1, setup%interaction_range
        write(*,'(a, I2, a, a)') ' Interchange interaction on shell: ', i, new_line('a')
          write(*, '(A)', advance='no') ' '
        do j=1, setup%n_species
          write(*, '(A, A)', advance='no') '          ', setup%species_names(j)
        end do
          write(*, '(A)', advance='yes') ''
        do j=1, setup%n_species
          write(*, '(A, A)', advance='no') ' ', setup%species_names(j)
          do k=1, setup%n_species
            write(*, '(A, F8.3, A)', advance='no') '  ', 1000*V_ex(k,j,i), '  '
            if (k .eq. setup%n_species) write(*,'(a)') ''
          end do
          if (j .eq. setup%n_species) write(*,'(a)') ''
        end do
      end do
      write(*,'(a)') ''
    ! Case meV
    else if (trim(aunits) .eq. 'meV') then
      write(*,"(x,'V_ij to be used in calculation (meV):',/)") 
      do i=1, setup%interaction_range
        write(*,'(a, I2, a, a)') ' Interchange interaction on shell: ', i, new_line('a')
          write(*, '(A)', advance='no') ' '
        do j=1, setup%n_species
          write(*, '(A, A)', advance='no') '          ', setup%species_names(j)
        end do
          write(*, '(A)', advance='yes') ''
        do j=1, setup%n_species
          write(*, '(A, A)', advance='no') ' ', setup%species_names(j)
          do k=1, setup%n_species
            write(*, '(A, F8.3, A)', advance='no') '  ', 1000*Ry_to_eV*V_ex(k,j,i), '  '
            if (k .eq. setup%n_species) write(*,'(a)') ''
          end do
          if (j .eq. setup%n_species) write(*,'(a)') ''
        end do
      end do
      write(*,'(a)') ''
    end if

  end subroutine pretty_print_exchange

  !> @brief   Subroutine to print the current state of the grid to the
  !>          screen, layer by layer.
  !>
  !> @details Currently supports the cubic representation of the fcc and
  !           bcc lattice types, as well as simple cubic.
  !>
  !> @author  C. D. Woodgate
  !> @date    2020-2023
  !>
  !> @param  grid Current simulation configuration
  !> @param  show_borders Optional argument controlling display of 
  !>                      borders
  !> @param  title Optional argument specifying the title to display
  !>
  !> @return None
  subroutine display_grid(grid, show_borders, title)

    integer(kind=array_int), intent(in), dimension(:,:,:,:) :: grid
    logical, intent(in), optional :: show_borders
    character(len=*), optional :: title
    logical :: borders
    integer, dimension(4) :: sizes
    integer :: ix, iy, iz
    character(len=1) :: c
    character(len=4), parameter :: clrstr = char(27)//'[2j'

    borders = .true.

    if (present(show_borders)) borders = show_borders

    write(*,'(a)') clrstr

    sizes = shape(grid)
    if (present(title)) then
      write(*, '(a)') title
    end if

    do iz = sizes(4), 1, -1

      write(*, '(a, I2)') 'Layer for z = ', iz

      if (borders) write(*, '(a)') repeat('=', sizes(2)+2)

      do iy = sizes(3), 1, -1
        if (borders) write(*, '(a)', advance='no') '|'
        do ix = 1, sizes(2)
          if (grid(1,ix,iy,iz) .ne. 0) then
            write(c, '(I1)') grid(1,ix, iy, iz)
          else
            c = ' '
          end if
          write(*, '(a)', advance='no') c
        end do
        if (borders) write(*, '(a)', advance='no') '|'
        write(*, '(a)') ''
      end do
      if (borders) write(*, '(a)') repeat('=', sizes(2)+2)
    end do

  end subroutine display_grid

  !> @brief   Centered printing routine  
  !>
  !> @details This routine prints centered text with user specified filler either side of printed message.
  !>          e.g. "----- Output -----". Optional argument controls whether a newline is inserted.
  !>          
  !> @param  message String to be printed
  !> @param  fill_char String to be used to fill message (can be " ")
  !> @param  newline If True, a newline is printed after the message. If False, not.
  !> 
  !> @return None
  !>
  !> @author  H. J. Naguszewski
  !> @author  C. D. Woodgate
  !> @date    2025 
  subroutine print_centered_message(message, fill_char, newline)
    implicit none
    character(len=*), intent(in) :: message, fill_char
    logical, optional :: newline
    character(len=72) :: output_str
    logical :: anewline
    integer :: num_fill_chars

    if (present(newline)) then
      anewline = newline
    else
      anewline = .False.
    end if

    ! Calculate the number of fill characters needed on each side
    num_fill_chars = (72 - len_trim(message) - 2) / 2

    ! Construct the formatted string with fill characters on both sides
    output_str = repeat(fill_char, num_fill_chars) // " " // message // " " // repeat(fill_char, num_fill_chars)

    ! If there's an odd number of fill characters, we need to add one more to the right side
    if (mod(72 - len_trim(message) - 2, 2) /= 0) then
        output_str = output_str // fill_char
    end if

    ! Write formatted output
    if (anewline) then
      write(6, '(A,/)') output_str
    else
      write(6, '(A)') output_str
    end if

end subroutine print_centered_message

end module display
