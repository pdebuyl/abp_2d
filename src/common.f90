module common
  implicit none

  private

  public :: rk
  public :: pi
  public :: normal_distribution
  public :: normal_value
  public :: get_character_argument

  integer, parameter :: rk = selected_real_kind(15)

  real(kind=rk), parameter :: pi = 4*atan(1._rk)

  ! prng state
  logical :: has_gauss = .false.
  real(kind=rk) :: gauss

  interface normal_distribution
     module procedure normal_distribution_1
     module procedure normal_distribution_2
  end interface normal_distribution

contains

  subroutine normal_distribution_1(data, scale)
    real(kind=rk), intent(inout) :: data(:)
    real(kind=rk), intent(in) :: scale

    integer :: ni, i

    ni = size(data, dim=1)

    do i = 1, ni
       data(i) = normal_value()*scale
    end do

  end subroutine normal_distribution_1

  subroutine normal_distribution_2(data, scale)
    real(kind=rk), intent(inout) :: data(:,:)
    real(kind=rk), intent(in) :: scale

    integer :: ni, nj, i, j

    ni = size(data, dim=1)
    nj = size(data, dim=2)

    do j = 1, nj
       do i = 1, ni
          data(i, j) = normal_value()*scale
       end do
    end do

  end subroutine normal_distribution_2

  function normal_value() result(r)
    real(kind=rk) :: r

    logical :: found
    real(kind=rk) :: u1, u2, radius

    if ( has_gauss ) then
       r = gauss
       has_gauss = .false.
    else
       found = .false.
       do while (.not. found)
          call random_number(u1)
          u1 = 2*u1 - 1
          call random_number(u2)
          u2 = 2*u2 - 1
          radius = (u1**2+u2**2)
          if ( ( radius < 1 ) .and. (radius > 0) ) found = .true.
       end do
       gauss = u1 * sqrt( -2 * log(radius)/radius )
       has_gauss = .true.
       r = u2 * sqrt( -2 * log(radius)/radius )
    end if

  end function normal_value

  !> Return the position-th command-line argument as an integer
  function get_character_argument(position) result(arg)
    integer, intent(in) :: position
    character(len=:), allocatable :: arg
    character(len=128) :: r

    if (command_argument_count() < position) then
       stop 'missing argument for parameter file'
    end if

    call get_command_argument(position, r)

    arg = trim(r)

  end function get_character_argument

end module common
