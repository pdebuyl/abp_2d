program run_abp_probe
  use common
  use abp_model_m
  implicit none

  integer, parameter :: N = 20

  type(abp_t) :: abp
  real(kind=rk), allocatable :: x1(:,:), force1(:,:)
  real(kind=rk), allocatable :: noise(:,:), theta_noise(:)

  integer :: i, j
  real(kind=rk) :: dt, dist
  logical :: too_close

  allocate(abp%x(2, N))
  allocate(abp%theta(N))
  allocate(abp%force(2, N))
  allocate(abp%D(N))
  allocate(abp%Dr(N))
  allocate(abp%v0(N))
  allocate(abp%sigma(N))
  allocate(abp%epsilon(N))

  abp%N = N
  abp%box_l = 20
  abp%box_l_i = 1/abp%box_l
  abp%v0 = 0.2_rk
  abp%Dr = 0.1_rk
  abp%D = 1._rk
  abp%sigma = 1._rk
  abp%epsilon = 1

  abp%D(1) = 0.1_rk
  abp%v0(1) = 0
  abp%sigma(1) = 5

  do i = 1, N
     too_close = .true.
     do while (too_close)
        call random_number(abp%x(:,i))
        abp%x(:,i) = abp%x(:,i)*abp%box_l
        too_close = .false.
        do j = 1, i-1
           dist = norm2(abp%min_dist(abp%x(:,i), abp%x(:,j)))
           if (dist <= 0.9_rk*(abp%sigma(i)+abp%sigma(j))) then
              too_close = .true.
              exit
           end if
        end do
     end do
  end do

  call random_number(abp%theta)
  abp%theta = abp%theta * 2 * pi

  allocate(x1(2, N))
  allocate(noise(2, N))
  allocate(theta_noise(N))
  allocate(force1(2, N))

  dt = 0.00001_rk

  do i = 1, 100

     call abp%srk_step(100, dt)

  end do

  dt = 0.001_rk

  do i = 1, 5000

     call abp%srk_step(100, dt)

     write(22,'(40e15.5)') abp%x
     write(23,'(20e15.5)') abp%theta
     write(24,'(40e15.5)') evec(abp%theta)


  end do


end program run_abp_probe
