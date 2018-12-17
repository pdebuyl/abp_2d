program test_threefry_1
  use threefry_m
  implicit none

  type(threefry_t) :: rng

  integer :: i
  real(kind=rk) :: x

  rng%c%c0 = 0
  rng%c%c1 = 0
  rng%k%c0 = 123456789
  rng%k%c1 = 987654321

  do i = 1, 10000
     write(*,*) threefry_double(rng%c, rng%k)
     rng%c%c0 = rng%c%c0 + 1
  end do

end program test_threefry_1
