module threefry_m
  use, intrinsic :: iso_c_binding
  use common
  implicit none

  type, bind(c) :: threefry_ctr_t
     integer(c_long) :: c0 = 0
     integer(c_long) :: c1 = 0
  end type threefry_ctr_t

  type, bind(c) :: threefry_key_t
     integer(c_long) :: c0 = 0
     integer(c_long) :: c1 = 0
  end type threefry_key_t

  type threefry_t
     type(threefry_ctr_t) :: c
     type(threefry_key_t) :: k
     real(kind=rk) :: normal_value
     logical :: has_normal = .false.
   contains
     procedure :: random_normal
     procedure :: urandom_seed
  end type threefry_t

  interface
     integer(kind=c_int) function threefry_int32(c, k) bind(c, name='threefry_int32')
       use, intrinsic :: iso_c_binding
       import threefry_ctr_t, threefry_key_t
       type(threefry_ctr_t) :: c
       type(threefry_key_t) :: k
     end function threefry_int32
  end interface

  interface
     real(c_double) function threefry_double(c, k) bind(c, name='threefry_c_double')
       use, intrinsic :: iso_c_binding
       import threefry_ctr_t, threefry_key_t
       type(threefry_ctr_t) :: c
       type(threefry_key_t) :: k
     end function threefry_double
  end interface

contains

  function random_normal(this) result(r)
    class(threefry_t), intent(inout) :: this
    real(kind=rk) :: r

    logical :: found
    real(kind=rk) :: u1, u2, radius

    if ( this%has_normal ) then
       r = this%normal_value
       this%has_normal = .false.
    else
       found = .false.
       do while (.not. found)
          u1 = 2*threefry_double(this%c, this%k) - 1
          u2 = 2*threefry_double(this%c, this%k) - 1
          radius = (u1**2+u2**2)
          if ( ( radius < 1 ) .and. (radius > 0) ) found = .true.
       end do
       this%normal_value = u1 * sqrt( -2 * log(radius)/radius )
       this%has_normal = .true.
       r = u2 * sqrt( -2 * log(radius)/radius )
    end if

  end function random_normal

  subroutine urandom_seed(this)
    class(threefry_t), intent(inout) :: this
    integer :: f_unit, error

    open(newunit=f_unit, file='/dev/urandom', status='old', access='stream', iostat=error)

    if (error /= 0) &
         stop 'error opening /dev/urandom in seed_rng'
    
    read(f_unit) this%k%c0

    close(f_unit)

  end subroutine urandom_seed

end module threefry_m
