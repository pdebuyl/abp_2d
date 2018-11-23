program test_min_dist
  use abp_model_m
  implicit none

  integer :: i
  real(kind=rk) :: noise(2), dt, x(2)

  dt = 3._rk

  x = 0

  do i = 1, 1000
     call normal_distribution(noise, scale=sqrt(2*dt))
     x = x + noise
     write(*,*) noise, x
  end do

end program test_min_dist

