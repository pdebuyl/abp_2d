program test_min_dist
  use common
  use abp_model_m
  implicit none

  type(abp_t) :: abp

  integer :: i, j
  real(kind=rk) :: x1(2), x2(2)

  abp%box_l = [5._rk, 6._rk]
  abp%box_l_i = 1/abp%box_l

  do i = 1, 1

     call random_number(x1)
     x1 = -2*abp%box_l + 5*x1*abp%box_l
     call random_number(x2)
     x2 = -2*abp%box_l + 5*x2*abp%box_l

     write(*,*) '---------------------------------'
     write(*,*) 'coordinates'
     write(*,*) x1, x2
     write(*,*) 'delta'
     write(*,*) x1-x2
     write(*,*) 'min_dist'
     write(*,*) abp%min_dist(x1, x2)

     x1 = modulo(x1 - x2, abp%box_l)
     do j = 1, 2
        if (x1(j) > abp%box_l(j)/2) then
           x1(j) = x1(j) - abp%box_l(j)
        else if (x1(j) < -abp%box_l(j)/2) then
           x1(j) = x1(j) + abp%box_l(j)
        end if
     end do

     write(*,*) 'pedestrian version'
     write(*,*) x1

  end do

end program test_min_dist

