program run_abp_probe
  use common
  use abp_model_m
  implicit none

  integer :: N

  type(abp_t) :: abp

  integer :: i
  real(kind=rk) :: dt

  integer :: pos_unit, theta_unit, vel_unit
  character(len=100) :: tmp
  character(len=:), allocatable :: coord_format, theta_format
  real(kind=rk) :: phi, sigma_abp, sigma_probe, l


  sigma_abp = 0.5_rk
  sigma_probe = 4.5_rk
  l = 30
  phi = 0.3_rk

  N = 1 + floor((phi*l**2 - pi*sigma_probe**2) / (pi*sigma_abp**2))

  write(*,*) 'ABP simulation'
  write(*,*) 'l', l
  write(*,*) 'phi', phi
  write(*,*) 'N', N

  write(*,*) pi*(sigma_probe**2 + (N-1)*sigma_abp**2) / l**2

  allocate(character(len=12) :: coord_format)
  write(coord_format, '(A,i5.5,A)') '(', 2*N, 'e15.5)'
  allocate(character(len=12) :: theta_format)
  write(theta_format, '(A,i5.5,A)') '(', N, 'e15.5)'

  open(newunit=pos_unit, file='position.txt')
  open(newunit=vel_unit, file='velocity.txt')
  open(newunit=theta_unit, file='theta.txt')

  call seed_rng

  call abp%init(N, L=[l, l])

  abp%v0 = 1._rk
  abp%Dr = 0.1_rk
  abp%D = 0.1_rk
  abp%mu = 1/abp%D
  abp%sigma = sigma_abp
  abp%epsilon = 1

  ! no noise and no self-prop for probe
  abp%D(1) = 0
  abp%v0(1) = 0
  abp%mu(1) = abp%mu(1)*sigma_abp/sigma_probe
  abp%sigma(1) = sigma_probe

  call abp%random_placement()

  dt = 0.00001_rk

  do i = 1, 100

     if (modulo(i, 100)==0) &
          write(*,*) 'thermalization', i

     call abp%srk_step(100, dt)

  end do

  dt = 0.001_rk

  do i = 1, 5000

     if (modulo(i, 100)==0) &
          write(*,*) 'integration', i

     call abp%srk_step(500, dt)

     write(pos_unit, coord_format) abp%x
     write(vel_unit, coord_format) abp%v
     write(theta_unit, theta_format) abp%theta

  end do


end program run_abp_probe
