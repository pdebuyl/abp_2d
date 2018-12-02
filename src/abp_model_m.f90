module abp_model_m
  use common
  implicit none

  private

  public :: abp_t
  public :: evec
  public :: cut_factor

  type abp_t
     real(kind=rk), allocatable :: x(:,:)
     real(kind=rk), allocatable :: v(:,:)
     real(kind=rk), allocatable :: theta(:)
     real(kind=rk), allocatable :: force(:,:)
     real(kind=rk), allocatable :: D(:)
     real(kind=rk), allocatable :: Dr(:)
     real(kind=rk), allocatable :: v0(:)
     real(kind=rk), allocatable :: sigma(:)
     real(kind=rk), allocatable :: epsilon(:)
     ! maybe a rng state
     integer :: N
     real(kind=rk) :: box_l(2)
     real(kind=rk) :: box_l_i(2)
     ! wall parameters
     ! trap parameters
     ! constant external force parameter
   contains
     procedure :: init
     procedure :: min_dist
     procedure :: compute_force
     procedure :: srk_step
  end type abp_t

  interface evec
     module procedure evec_s
     module procedure evec_1
  end interface evec

  real(kind=rk), parameter :: cut_factor = 2._rk**(1._rk/6._rk)

contains

  subroutine init(this, N, L)
    class(abp_t), intent(out) :: this
    integer, intent(in) :: N
    real(kind=rk), intent(in) :: L(2)

    allocate(this%x(2, N))
    allocate(this%v(2, N))
    allocate(this%theta(N))
    allocate(this%force(2, N))
    allocate(this%D(N))
    allocate(this%Dr(N))
    allocate(this%v0(N))
    allocate(this%sigma(N))
    allocate(this%epsilon(N))

    this%N = N
    this%box_l = L
    this%box_l_i = 1/L

  end subroutine init

  function min_dist(this, x1, x2) result(r)
    class(abp_t), intent(in) :: this
    real(kind=rk), intent(in) :: x1(2), x2(2)
    real(kind=rk) :: r(2)

    r = x1 - x2
    r = r - floor(r*this%box_l_i+0.5_rk)*this%box_l

  end function min_dist

  subroutine compute_force(this)
    class(abp_t), intent(inout) :: this

    integer :: i, j
    real(kind=rk) :: dist(2), r, f(2)
    real(kind=rk) :: sigma, epsilon


    this%force = 0

    do i = 1, this%N
       do j = i+1, this%N

          dist = this%min_dist(this%x(:,i), this%x(:,j))
          sigma = this%sigma(i)+this%sigma(j)
          epsilon = 1

          r = norm2(dist)
          if (r < cut_factor*sigma) then
             f = lj_force(dist, r, sigma, epsilon)
             this%force(:,i) = this%force(:,i) + f
             this%force(:,j) = this%force(:,j) - f
          end if

       end do
    end do

  end subroutine compute_force

  function lj_force(x, r, sigma, epsilon) result(f)
    real(kind=rk), intent(in) :: x(2)
    real(kind=rk), intent(in) :: r
    real(kind=rk), intent(in) :: sigma
    real(kind=rk), intent(in) :: epsilon
    real(kind=rk) :: f(2)

    real(kind=rk) :: r_sq, sig6_o_r6

    r_sq = r**2
    sig6_o_r6 = sigma**6/r_sq**3

    f = 24*epsilon* sig6_o_r6/r_sq * (2*sig6_o_r6 - 1) * x

  end function lj_force


  function evec_s(theta) result(r)
    real(kind=rk), intent(in) :: theta
    real(kind=rk) :: r(2)

    r(1) = cos(theta)
    r(2) = sin(theta)

  end function evec_s

  function evec_1(theta) result(r)
    real(kind=rk), intent(in) :: theta(:)
    real(kind=rk) :: r(2, size(theta))

    integer :: i

    do i = 1, size(theta)
       r(:,i) = evec_s(theta(i))
    end do

  end function evec_1

  subroutine srk_step(this, nsteps, dt)
    class(abp_t), intent(inout) :: this
    integer, intent(in) :: nsteps
    real(kind=rk), intent(in) :: dt

    logical, save :: first = .true.
    real(kind=rk), allocatable, save :: noise(:,:), theta_noise(:)
    real(kind=rk), allocatable, save :: x1(:,:), force1(:,:)

    integer :: ii, j

    if (first) then
       allocate(noise(2, this%N))
       allocate(theta_noise(this%N))
       allocate(x1(2, this%N))
       allocate(force1(2, this%N))
       first = .false.
    end if

     do ii = 1, nsteps
        ! Keep original positions for final update
        x1 = this%x

        ! Generate the noise
        call normal_distribution(noise, scale=sqrt(2*dt))
        call normal_distribution(theta_noise, scale=sqrt(2*dt))

        ! compute force at begin of timestep
        call this%compute_force
        force1 = this%force

        ! First update of the coordinates
        do j = 1, this%N
           this%x(:,j) = this%x(:,j) &
                + this%D(j)*force1(:,j)*dt &
                + this%v0(j) * evec(this%theta(j))*dt &
                + noise(:,j)*sqrt(this%D(j))
           this%theta(j) = this%theta(j) + theta_noise(j)*sqrt(this%Dr(j))
        end do

        ! compute force at end of timestep
        call this%compute_force

        do j = 1, this%N
           this%v(:,j) = this%D(j)*(force1(:,j)+this%force(:,j))/2 &
                + this%v0(j) * evec(this%theta(j))

           this%x(:,j) = x1(:,j) + this%v(:,j) * dt &
                + noise(:,j)*sqrt(this%D(j))

        end do

     end do

   end subroutine srk_step

end module abp_model_m
