program run_abp_probe
  use common
  use abp_model_m
  implicit none

  integer, parameter :: N = 20

  type(abp_t) :: abp

  integer :: i
  real(kind=rk) :: dt

  integer :: pos_unit, theta_unit, vel_unit

  open(newunit=pos_unit, file='position.txt')
  open(newunit=vel_unit, file='velocity.txt')
  open(newunit=theta_unit, file='theta.txt')

  call abp%init(N, L=[16.6_rk, 16.6_rk])

  abp%v0 = 1._rk
  abp%Dr = 0.01_rk
  abp%D = 1._rk
  abp%sigma = 1._rk
  abp%epsilon = 1

  abp%D(1) = 0.1_rk
  abp%v0(1) = 0
  abp%sigma(1) = 5

  call abp%random_placement()

  dt = 0.00001_rk

  do i = 1, 100

     call abp%srk_step(100, dt)

  end do

  dt = 0.001_rk

  do i = 1, 5000

     call abp%srk_step(500, dt)

     write(pos_unit,'(40e15.5)') abp%x
     write(vel_unit,'(40e15.5)') abp%v
     write(theta_unit,'(20e15.5)') abp%theta

  end do


end program run_abp_probe
