module abp_model_m
  use common
  implicit none

  private

  public :: abp_t
  public :: evec
  public :: cut_factor

  type pair_list_t
     real(kind=rk) :: cell_l(2)
     integer, allocatable :: cell_list(:,:,:)
     integer, allocatable :: cell_count(:,:)
  end type pair_list_t

  type p_list_t
     integer :: count
     integer, allocatable :: idx(:)
  end type p_list_t

  type abp_t
     real(kind=rk), allocatable :: x(:,:)
     real(kind=rk), allocatable :: x_old(:,:)
     real(kind=rk), allocatable :: v(:,:)
     real(kind=rk), allocatable :: theta(:)
     real(kind=rk), allocatable :: force(:,:)
     real(kind=rk), allocatable :: D(:)
     real(kind=rk), allocatable :: mu(:)
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
     type(pair_list_t) :: pairs
     type(p_list_t), allocatable :: p_list(:)
     integer :: p_list_max
   contains
     procedure :: init
     procedure :: random_placement
     procedure :: min_dist
     procedure :: compute_force
     procedure :: srk_step
     procedure :: index_from_position
     procedure :: make_list
     procedure :: compute_force_list
  end type abp_t

  interface evec
     module procedure evec_s
     module procedure evec_1
  end interface evec

  real(kind=rk), parameter :: cut_factor = 2._rk**(1._rk/6._rk)
  real(kind=rk), parameter :: skin = 1

contains

  subroutine init(this, N, L)
    class(abp_t), intent(out) :: this
    integer, intent(in) :: N
    real(kind=rk), intent(in) :: L(2)

    allocate(this%x(2, N))
    allocate(this%x_old(2, N))
    allocate(this%v(2, N))
    allocate(this%theta(N))
    allocate(this%force(2, N))
    allocate(this%D(N))
    allocate(this%mu(N))
    allocate(this%Dr(N))
    allocate(this%v0(N))
    allocate(this%sigma(N))
    allocate(this%epsilon(N))

    this%N = N
    this%box_l = L
    this%box_l_i = 1/L

  end subroutine init

  subroutine random_placement(this)
    class(abp_t), intent(inout) :: this

    integer :: i, j
    logical :: too_close
    real(kind=rk) :: dist

    do i = 1, this%N
       too_close = .true.
       do while (too_close)
          call random_number(this%x(:,i))
          this%x(:,i) = this%x(:,i)*this%box_l
          too_close = .false.
          do j = 1, i-1
             dist = norm2(this%min_dist(this%x(:,i), this%x(:,j)))
             if (dist <= 0.9_rk*(this%sigma(i)+this%sigma(j))) then
                too_close = .true.
                exit
             end if
          end do
       end do
    end do

  call random_number(this%theta)
  this%theta = this%theta * 2 * pi

  end subroutine random_placement


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
    real(kind=rk) :: max_move

    integer :: ii, j

    if (first) then
       allocate(noise(2, this%N))
       allocate(theta_noise(this%N))
       allocate(x1(2, this%N))
       allocate(force1(2, this%N))
       first = .false.
    end if

    this%x_old = this%x
    max_move = huge(max_move)

     do ii = 1, nsteps
        ! Keep original positions for final update
        x1 = this%x

        ! Generate the noise
        call normal_distribution(noise, scale=sqrt(2*dt))
        call normal_distribution(theta_noise, scale=sqrt(2*dt))

        ! compute force at begin of timestep
        if (max_move > 0.5_rk) call this%make_list
        call this%compute_force_list
        force1 = this%force

        max_move = 0
        ! First update of the coordinates
        do j = 1, this%N
           this%x(:,j) = this%x(:,j) &
                + this%mu(j)*force1(:,j)*dt &
                + this%v0(j) * evec(this%theta(j))*dt &
                + noise(:,j)*sqrt(this%D(j))
           max_move = max(max_move, norm2(this%x(:,j)-this%x_old(:,j)))
           this%theta(j) = this%theta(j) + theta_noise(j)*sqrt(this%Dr(j))
        end do

        ! compute force at end of timestep
        if (max_move > 0.5_rk) call this%make_list
        call this%compute_force_list

        max_move = 0
        do j = 1, this%N
           this%v(:,j) = this%v0(j) * evec(this%theta(j)) &
                + this%mu(j)*(force1(:,j)+this%force(:,j))/2

           this%x(:,j) = x1(:,j) + this%v(:,j) * dt &
                + noise(:,j)*sqrt(this%D(j))
           max_move = max(max_move, norm2(this%x(:,j)-this%x_old(:,j)))

        end do

     end do

   end subroutine srk_step

   pure function index_from_position(this, i, r_max) result(idx)
     class(abp_t), intent(in) :: this
     integer, intent(in) :: i
     real(kind=rk), intent(in) :: r_max

     integer :: idx(2)

     idx = floor(modulo(this%x(:,i), this%box_l)/this%pairs%cell_l) + 1

   end function index_from_position

   subroutine make_list(this)
     class(abp_t), intent(inout) :: this

     logical, save :: first = .true.

     integer :: i, c, idx_xy(2)
     integer :: n_x, n_y
     integer :: n_max
     real(kind=rk) :: r_max, rsq, dist(2), sigma
     integer :: stencil(2, 4)
     integer :: cell_i, cell_j, n1, n2, part_1, part_2, idx_1, idx_2, count
     integer :: ncell_idx(2), ncell_i

    stencil(1,:) = [-1, 0, 1, -1]
    stencil(2,:) = [-1, -1, -1, 0]

    if (first) then
       r_max = (this%sigma(1) + this%sigma(2))*cut_factor + 1
       n_x = floor(this%box_l(1)/r_max)
       n_y = floor(this%box_l(2)/r_max)
       this%pairs%cell_l = this%box_l / [n_x, n_y]
       n_max = int(maxval(this%pairs%cell_l)**2 / (pi*minval(this%sigma)**2))
       first = .false.
       allocate(this%pairs%cell_list(n_max, n_y, n_x))
       allocate(this%pairs%cell_count(n_y, n_x))
       allocate(this%p_list(this%N))
       this%p_list_max = int(pi*((this%sigma(1) + this%sigma(2))*cut_factor + skin)**2)
       do i = 1, this%N
          allocate(this%p_list(i)%idx(this%p_list_max))
       end do
    else
       n_x = size(this%pairs%cell_count, dim=2)
       n_y = size(this%pairs%cell_count, dim=1)
    end if

    n_max = size(this%pairs%cell_list, dim=1)

    this%pairs%cell_count = 0
    do i = 1, this%N
       idx_xy = this%index_from_position(i, r_max)
       c = this%pairs%cell_count(idx_xy(2), idx_xy(1)) + 1
       if (c > n_max) stop 'exceed n_max in make_list'

       this%pairs%cell_list(c, idx_xy(2), idx_xy(1)) = i
       this%pairs%cell_count(idx_xy(2), idx_xy(1)) = c
    end do

    do i = 1, this%N
       this%p_list(i)%count = 0
    end do

    !$omp parallel
    !$omp do private(cell_i, cell_j, n1, idx_1, part_1, idx_2, part_2, &
    !$omp dist, sigma, rsq, ncell_idx, n2, c)
    do cell_i = 1, n_x
       do cell_j = 1, n_y

          ! cell self-interaction
          n1 = this%pairs%cell_count(cell_j, cell_i)

          do idx_1 = 1, n1
             part_1 = this%pairs%cell_list(idx_1, cell_j, cell_i)

             do idx_2 = idx_1 + 1, n1

                part_2 = this%pairs%cell_list(idx_2, cell_j, cell_i)

                dist = this%min_dist(this%x(:,part_1), this%x(:,part_2))
                sigma = this%sigma(part_1)+this%sigma(part_2)

                rsq = dist(1)**2 + dist(2)**2
                if (rsq < (cut_factor*sigma + skin)**2) then
                   ! add pair to list of part_1
                   c = this%p_list(part_1)%count + 1
                   if (c > this%p_list_max) stop 'count exceeds p_list_max in make_list'
                   this%p_list(part_1)%idx(c) = part_2
                   this%p_list(part_1)%count = c
                end if

             end do
          end do

          ! cell pairs
          do ncell_i = 1, size(stencil, dim=2)
             ncell_idx = [cell_i, cell_j] + stencil(:,ncell_i)
             ncell_idx = modulo(ncell_idx - 1, [n_x, n_y]) + 1

             n2 = this%pairs%cell_count(ncell_idx(2), ncell_idx(1))

             do idx_1 = 1, n1
                part_1 = this%pairs%cell_list(idx_1, cell_j, cell_i)

                do idx_2 = 1, n2

                   part_2 = this%pairs%cell_list(idx_2, ncell_idx(2), ncell_idx(1))

                   dist = this%min_dist(this%x(:,part_1), this%x(:,part_2))
                   sigma = this%sigma(part_1)+this%sigma(part_2)

                   rsq = dist(1)**2 + dist(2)**2
                   if (rsq < (cut_factor*sigma + skin)**2) then
                      ! add pair to list of part_1
                      c = this%p_list(part_1)%count + 1
                      if (c > this%p_list_max) stop 'count exceeds p_list_max in make_list'
                      this%p_list(part_1)%idx(c) = part_2
                      this%p_list(part_1)%count = c
                   end if

                end do
             end do


          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    this%x_old = this%x

   end subroutine make_list

  subroutine compute_force_list(this)
    class(abp_t), intent(inout) :: this

    real(kind=rk) :: dist(2), r, f(2)
    real(kind=rk) :: sigma, epsilon

    integer :: part_1, part_2, j
    integer :: n_in_cell

    logical, save :: first = .true.
    real(kind=rk), allocatable, save :: force(:,:)

    if (first) then
       allocate(force(2, this%N))
       first = .false.
    end if

    force = 0

    do part_1 = 1, this%N
       n_in_cell = this%p_list(part_1)%count

       do j = 1, n_in_cell

          part_2 = this%p_list(part_1)%idx(j)

          dist = this%min_dist(this%x(:,part_1), this%x(:,part_2))
          sigma = this%sigma(part_1)+this%sigma(part_2)
          epsilon = 1

          r = dist(1)**2 + dist(2)**2
          if (r < (cut_factor*sigma)**2) then
             r = sqrt(r)
             f = lj_force(dist, r, sigma, epsilon)
             force(:,part_1) = force(:,part_1) + f
             force(:,part_2) = force(:,part_2) - f
          end if

       end do
    end do

    this%force = force

  end subroutine compute_force_list


end module abp_model_m
