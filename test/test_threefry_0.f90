program test_threefry_0
  use threefry_m
  implicit none

  type(threefry_t) :: rng

  integer :: i
  real(kind=rk) :: x

  rng%c%c0 = 0
  rng%c%c1 = 0
  rng%k%c0 = z'deadbeef'
  rng%k%c1 = z'badcafe'

  do i = 1, 10
     write(*,*) threefry_int32(rng%c, rng%k)
  end do

end program test_threefry_0
