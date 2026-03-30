module poisson
  ! 1D Poisson solver using box-integration (finite volume) discretization.
  ! Solves: d/dz[eps(z) dPhi/dz] = -rho(z) / e0
  !
  ! Uses the Thomas algorithm (O(N)) for the resulting tridiagonal system.
  ! Supports Dirichlet-Dirichlet ('DD') and Dirichlet-Neumann ('DN') BCs.

  use definitions, only: dp, e0
  implicit none

  private
  public :: solve_poisson

contains

  subroutine solve_poisson(phi, rho, epsilon, dz, N, bc_left, bc_right, bc_type)
    ! Solve 1D Poisson equation with position-dependent dielectric.
    !
    ! Box-integration discretization (Birner et al., Acta Phys. Pol. A 110, 111 (2006)):
    !   [eps_{i+1/2}(Phi_{i+1} - Phi_i) - eps_{i-1/2}(Phi_i - Phi_{i-1})] / dz^2
    !       = -rho_i / e0
    !
    ! where eps_{i+-1/2} = (eps_i + eps_{i+-1}) / 2  (arithmetic average).
    !
    ! Inputs:
    !   rho(N)     - charge density (C/nm^3)
    !   epsilon(N) - relative dielectric constant (unitless)
    !   dz         - uniform grid spacing (nm)
    !   N          - number of grid points (>= 3)
    !   bc_left    - left boundary value (V)
    !   bc_right   - right boundary value (V), used only for DD
    !   bc_type    - 'DD' (Dirichlet-Dirichlet) or 'DN' (Dirichlet-Neumann)
    !
    ! Output:
    !   phi(N)     - electrostatic potential (V)

    real(kind=dp), intent(inout) :: phi(N)
    real(kind=dp), intent(in)    :: rho(N)
    real(kind=dp), intent(in)    :: epsilon(N)
    real(kind=dp), intent(in)    :: dz
    integer, intent(in)          :: N
    real(kind=dp), intent(in)    :: bc_left
    real(kind=dp), intent(in)    :: bc_right
    character(len=2), intent(in) :: bc_type

    ! Tridiagonal bands: a = main diagonal, b = upper, c = lower
    real(kind=dp), allocatable :: a(:), b(:), c(:), rhs(:)
    real(kind=dp) :: eps_plus, eps_minus
    integer :: i

    allocate(a(N), b(N), c(N), rhs(N))

    ! --- Interior points (i = 2 .. N-1) ---
    do i = 2, N - 1
      eps_plus  = 0.5_dp * (epsilon(i) + epsilon(i + 1))
      eps_minus = 0.5_dp * (epsilon(i) + epsilon(i - 1))

      b(i) = eps_plus  / dz**2   ! upper diagonal (coefficient of Phi_{i+1})
      c(i) = eps_minus / dz**2   ! lower diagonal (coefficient of Phi_{i-1})
      a(i) = -(b(i) + c(i))      ! main diagonal  (coefficient of Phi_i)
      rhs(i) = -rho(i) / e0      ! right-hand side
    end do

    ! --- Boundary conditions ---
    select case (trim(bc_type))
    case ('DD')
      ! Left: Dirichlet  Phi_1 = bc_left
      a(1)   = 1.0_dp
      b(1)   = 0.0_dp
      rhs(1) = bc_left

      ! Right: Dirichlet  Phi_N = bc_right
      a(N)   = 1.0_dp
      c(N)   = 0.0_dp
      rhs(N) = bc_right

    case ('DN')
      ! Left: Dirichlet  Phi_1 = bc_left
      a(1)   = 1.0_dp
      b(1)   = 0.0_dp
      rhs(1) = bc_left

      ! Right: Neumann  dPhi/dz = 0  =>  Phi_N - Phi_{N-1} = 0
      a(N)   = 1.0_dp
      c(N)   = -1.0_dp
      rhs(N) = 0.0_dp

    case default
      ! Unknown BC type — set phi to zero and return
      phi = 0.0_dp
      deallocate(a, b, c, rhs)
      return
    end select

    ! --- Solve tridiagonal system via Thomas algorithm ---
    call thomas_solve(a, b, c, rhs, N)

    ! --- Copy solution into output ---
    phi = rhs

    deallocate(a, b, c, rhs)

  end subroutine solve_poisson


  subroutine thomas_solve(a, b, c, d, N)
    ! Thomas algorithm for tridiagonal system A*x = d.
    !
    ! On entry:
    !   a(N) - main diagonal
    !   b(N) - upper diagonal (b(N) unused)
    !   c(N) - lower diagonal (c(1) unused)
    !   d(N) - right-hand side
    !
    ! On exit:
    !   d(N) - solution vector x
    !
    ! The system is:
    !   c(i)*x(i-1) + a(i)*x(i) + b(i)*x(i+1) = d(i)
    !
    ! with c(1)=0 and b(N)=0.

    real(kind=dp), intent(inout) :: a(N), b(N), c(N), d(N)
    integer, intent(in) :: N

    real(kind=dp) :: w
    integer :: i

    ! --- Forward sweep: eliminate lower diagonal ---
    do i = 2, N
      w = c(i) / a(i - 1)
      a(i) = a(i) - w * b(i - 1)
      d(i) = d(i) - w * d(i - 1)
    end do

    ! --- Back-substitution ---
    d(N) = d(N) / a(N)
    do i = N - 1, 1, -1
      d(i) = (d(i) - b(i) * d(i + 1)) / a(i)
    end do

  end subroutine thomas_solve

end module poisson
