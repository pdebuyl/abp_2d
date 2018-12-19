program run_abp_probe
  use common
  use config_m
  use abp_model_m
  implicit none

  integer :: N

  type(abp_t) :: abp
  type(config_t) :: conf

  integer :: i, n_therm, n_loop, n_steps
  real(kind=rk) :: dt, probe_T

  integer :: pos_unit, theta_unit, vel_unit
  character(len=:), allocatable :: coord_format, theta_format
  real(kind=rk) :: phi, sigma_abp, sigma_probe, l
  integer, allocatable :: idx_data(:)
  real(kind=rk), allocatable :: io_x(:,:), io_v(:,:), io_th(:)
  integer, allocatable :: io_id(:)

  call conf%init(filename=get_character_argument(1))

  sigma_abp = conf%get_d('sigma')
  sigma_probe = conf%get_d('sigma_probe')
  l = conf%get_d('box_l')
  phi = conf%get_d('phi')

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

  allocate(idx_data(N))
  allocate(io_x(2,N))
  allocate(io_v(2,N))
  allocate(io_th(N))
  allocate(io_id(N))

  call abp%init(N, L=[l, l])
  do i = 1, size(abp%rng)
     call abp%rng(i)%urandom_seed()
  end do

  abp%v0 = conf%get_d('v0')
  abp%Dr = conf%get_d('Dr')
  abp%D = conf%get_d('D')
  abp%mu = 1/abp%D
  abp%sigma = sigma_abp
  abp%epsilon = 1

  dt = conf%get_d('dt')
  n_therm = conf%get_i('n_therm')
  n_loop = conf%get_i('n_loop')
  n_steps = conf%get_i('n_steps')

  ! no noise and no self-prop for probe
  probe_T = conf%get_d('probe_T')
  abp%D(1) = probe_T*abp%D(1)*sigma_abp/sigma_probe
  abp%v0(1) = 0
  abp%mu(1) = abp%mu(1)*sigma_abp/sigma_probe
  abp%sigma(1) = sigma_probe

  call abp%random_placement()

  do i = 1, n_therm

     if (modulo(i, 100)==0) &
          write(*,*) 'thermalization', i

     call abp%srk_step(n_steps, dt/10)

  end do

  do i = 1, n_loop

     if (modulo(i, 100)==0) &
          write(*,*) 'integration', i

     call abp%srk_step(n_steps, dt)

     call write_dumps

  end do

contains

  subroutine write_dumps

    integer :: i, idx

    do i = 1, N
       idx = abp%id(i)
       io_x(:,idx) = abp%x(:,i)
       io_v(:,idx) = abp%v(:,i)
       io_th(idx) = abp%theta(i)
       io_id(idx) = abp%id(i)
    end do

    write(pos_unit, coord_format) io_x
    write(vel_unit, coord_format) io_v
    write(theta_unit, theta_format) io_th

  end subroutine write_dumps

end program run_abp_probe
