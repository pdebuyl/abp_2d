module threefry_m
  use, intrinsic :: iso_c_binding
  use common
  implicit none

  type, bind(c) :: threefry_ctr_t
     integer(c_long) :: c0
     integer(c_long) :: c1
  end type threefry_ctr_t

  type, bind(c) :: threefry_key_t
     integer(c_long) :: c0
     integer(c_long) :: c1
  end type threefry_key_t

  type threefry_t
     type(threefry_ctr_t) :: c
     type(threefry_key_t) :: k
     real(kind=rk) :: normal
     logical :: has_normal = .false.
  end type threefry_t

  interface
     integer function threefry_int32(c, k) bind(c, name='threefry_int32')
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

end module threefry_m
