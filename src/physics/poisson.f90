module poisson
  ! 1D and 2D Poisson solvers using box-integration discretization.
  !
  ! 1D: d/dz[eps(z) dPhi/dz] = -rho(z) / e0
  !     Solved via Thomas algorithm (O(N)) for tridiagonal systems.
  !
  ! 2D: d/dy[eps(y,z) dPhi/dy] + d/dz[eps(y,z) dPhi/dz] = -rho(y,z) / e0
  !     Solved via MKL PARDISO sparse direct solver on Ny x Nz grids.

  use definitions, only: dp, e0
  implicit none

  private
  public :: solve_poisson
  public :: solve_poisson_2d
  public :: bc_from_string
  public :: BC_DD, BC_DN

  ! Boundary condition integer constants
  integer, parameter :: BC_DD = 1   ! Dirichlet-Dirichlet
  integer, parameter :: BC_DN = 2   ! Dirichlet-Neumann

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
    !   bc_type    - BC_DD (Dirichlet-Dirichlet) or BC_DN (Dirichlet-Neumann)
    !
    ! Output:
    !   phi(N)     - electrostatic potential (V)

    integer, intent(in)          :: N
    real(kind=dp), intent(in)    :: dz
    real(kind=dp), intent(inout) :: phi(N)
    real(kind=dp), intent(in)    :: rho(N)
    real(kind=dp), intent(in)    :: epsilon(N)
    real(kind=dp), intent(in)    :: bc_left
    real(kind=dp), intent(in)    :: bc_right
    integer, intent(in)          :: bc_type

    ! Tridiagonal bands: a = main diagonal, b = upper, c = lower
    real(kind=dp), allocatable :: a(:), b(:), c(:), rhs(:)
    real(kind=dp) :: eps_plus, eps_minus
    integer :: i

    if (N < 3) then
      phi = 0.0_dp
      return
    end if

    allocate(a(N), b(N), c(N), rhs(N))
    a = 0.0_dp; b = 0.0_dp; c = 0.0_dp; rhs = 0.0_dp

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
    ! Left BC: Dirichlet Phi_1 = bc_left (common to DD and DN)
    a(1)   = 1.0_dp
    b(1)   = 0.0_dp
    rhs(1) = bc_left

    select case (bc_type)
    case (BC_DD)
      ! Right: Dirichlet  Phi_N = bc_right
      a(N)   = 1.0_dp
      c(N)   = 0.0_dp
      rhs(N) = bc_right

    case (BC_DN)
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

    integer, intent(in) :: N
    real(kind=dp), intent(inout) :: a(N), b(N), c(N), d(N)

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


  ! ------------------------------------------------------------------
  ! Convert string BC type to integer constant (boundary conversion)
  ! ------------------------------------------------------------------
  function bc_from_string(bc_str) result(bc)
    character(len=*), intent(in) :: bc_str
    integer :: bc

    select case (trim(adjustl(bc_str)))
    case ('DD')
      bc = BC_DD
    case ('DN')
      bc = BC_DN
    case default
      bc = 0
    end select

  end function bc_from_string


  ! ==================================================================
  ! 2D Poisson solver with position-dependent dielectric.
  !
  ! Solves:
  !   d/dy[eps(y,z) dPhi/dy] + d/dz[eps(y,z) dPhi/dz] = -rho(y,z) / e0
  !
  ! on an Ny x Nz rectangular grid using box-integration discretization
  ! and MKL PARDISO sparse direct solver.  Dirichlet BCs on all boundaries.
  !
  ! Flat index mapping: idx = (iy - 1)*nz + iz  (column-major).
  ! ==================================================================
  subroutine solve_poisson_2d(phi, rho, epsilon, dy, dz, ny, nz, bc_value)

    integer, intent(in)          :: ny, nz
    real(kind=dp), intent(in)    :: dy, dz
    real(kind=dp), intent(inout) :: phi(ny, nz)
    real(kind=dp), intent(in)    :: rho(ny, nz)
    real(kind=dp), intent(in)    :: epsilon(ny, nz)
    real(kind=dp), intent(in)    :: bc_value

    integer :: ntotal, nnz_est, coo_idx
    integer :: iy, iz, row, iy_n, iz_n

    ! COO assembly arrays (real-valued for Poisson)
    integer, allocatable  :: coo_rows(:), coo_cols(:)
    real(kind=dp), allocatable :: coo_vals(:)
    real(kind=dp), allocatable :: rhs(:), sol(:)

    ! PARDISO arrays
    integer(8) :: pt(64)
    integer :: iparm(64)
    integer :: maxfct, mnum, mtype, phase, nrhs, msglvl, error
    integer, allocatable :: perm(:)

    ! CSR arrays for PARDISO
    integer, allocatable :: ia(:), ja(:)
    real(kind=dp), allocatable :: a_csr(:)

    ! Discretization variables
    real(kind=dp) :: eps_yp, eps_ym, eps_zp, eps_zm
    real(kind=dp) :: dy2_inv, dz2_inv

    ! Guard: grid too small
    if (ny < 3 .or. nz < 3) then
      phi = 0.0_dp
      return
    end if

    ntotal = ny * nz
    dy2_inv = 1.0_dp / dy**2
    dz2_inv = 1.0_dp / dz**2

    ! Estimate COO size: each interior point has 5 entries, boundary has 1.
    ! Upper bound: 5 * ntotal.
    nnz_est = 5 * ntotal
    allocate(coo_rows(nnz_est), coo_cols(nnz_est), coo_vals(nnz_est))
    coo_idx = 0

    allocate(rhs(ntotal), sol(ntotal))
    rhs = 0.0_dp
    sol = 0.0_dp

    ! ================================================================
    ! Assemble sparse matrix in COO format
    ! ================================================================
    do iy = 1, ny
      do iz = 1, nz
        row = (iy - 1) * nz + iz

        if (iy == 1 .or. iy == ny .or. iz == 1 .or. iz == nz) then
          ! Boundary: Dirichlet BC
          call poisson_add_coo(coo_rows, coo_cols, coo_vals, &
            nnz_est, coo_idx, row, row, 1.0_dp)
          rhs(row) = bc_value

        else
          ! Interior point: 5-point stencil
          eps_yp = 0.5_dp * (epsilon(iy, iz) + epsilon(iy + 1, iz))
          eps_ym = 0.5_dp * (epsilon(iy, iz) + epsilon(iy - 1, iz))
          eps_zp = 0.5_dp * (epsilon(iy, iz) + epsilon(iy, iz + 1))
          eps_zm = 0.5_dp * (epsilon(iy, iz) + epsilon(iy, iz - 1))

          ! Diagonal
          call poisson_add_coo(coo_rows, coo_cols, coo_vals, &
            nnz_est, coo_idx, row, row, &
            -(eps_yp + eps_ym) * dy2_inv - (eps_zp + eps_zm) * dz2_inv)

          ! +y neighbor (iy+1, iz)
          iy_n = (iy) * nz + iz
          call poisson_add_coo(coo_rows, coo_cols, coo_vals, &
            nnz_est, coo_idx, row, iy_n, eps_yp * dy2_inv)

          ! -y neighbor (iy-1, iz)
          iy_n = (iy - 2) * nz + iz
          call poisson_add_coo(coo_rows, coo_cols, coo_vals, &
            nnz_est, coo_idx, row, iy_n, eps_ym * dy2_inv)

          ! +z neighbor (iy, iz+1)
          iz_n = row + 1
          call poisson_add_coo(coo_rows, coo_cols, coo_vals, &
            nnz_est, coo_idx, row, iz_n, eps_zp * dz2_inv)

          ! -z neighbor (iy, iz-1)
          iz_n = row - 1
          call poisson_add_coo(coo_rows, coo_cols, coo_vals, &
            nnz_est, coo_idx, row, iz_n, eps_zm * dz2_inv)

          ! RHS
          rhs(row) = -rho(iy, iz) / e0
        end if
      end do
    end do

    ! ================================================================
    ! Convert COO to CSR for PARDISO
    ! ================================================================
    call poisson_coo_to_csr(ntotal, ntotal, coo_idx, &
      coo_rows, coo_cols, coo_vals, ia, ja, a_csr)

    deallocate(coo_rows, coo_cols, coo_vals)

    ! ================================================================
    ! Solve with MKL PARDISO (real unsymmetric, mtype=11)
    ! ================================================================
    pt = 0
    iparm = 0
    maxfct = 1
    mnum = 1
    nrhs = 1
    error = 0

    mtype = 11  ! real unsymmetric
    msglvl = 0  ! no output

    iparm(1)  = 1   ! no solver default
    iparm(2)  = 2   ! fill-in reordering from METIS
    iparm(8)  = 10  ! max iterative refinement steps
    iparm(10) = 13  ! perturb pivots with 1E-13
    iparm(11) = 1   ! nonsymmetric permutation and scaling
    iparm(13) = 1   ! maximum weighted matching
    iparm(24) = 2   ! parallel factorization

    allocate(perm(ntotal))

    ! Phase 13: Analysis + Factorization + Solve
    phase = 13
    call pardiso(pt, maxfct, mnum, mtype, phase, ntotal, a_csr, ia, ja, &
      perm, nrhs, iparm, msglvl, rhs, sol, error)

    if (error /= 0) then
      print *, 'ERROR: 2D Poisson solve failed (PARDISO error=', error, '). Aborting.'
      ! Release PARDISO memory
      phase = -1
      call pardiso(pt, maxfct, mnum, mtype, phase, ntotal, a_csr, ia, ja, &
        perm, nrhs, iparm, msglvl, rhs, sol, error)
      deallocate(ia, ja, a_csr, rhs, sol, perm)
      stop 1
    end if

    ! Phase -1: Release memory
    phase = -1
    call pardiso(pt, maxfct, mnum, mtype, phase, ntotal, a_csr, ia, ja, &
      perm, nrhs, iparm, msglvl, rhs, sol, error)

    deallocate(ia, ja, a_csr, rhs, perm)

    ! ================================================================
    ! Extract solution from flat array to 2D
    ! ================================================================
    do iy = 1, ny
      do iz = 1, nz
        phi(iy, iz) = sol((iy - 1) * nz + iz)
      end do
    end do

    deallocate(sol)

  end subroutine solve_poisson_2d


  ! ------------------------------------------------------------------
  ! Add a COO entry for 2D Poisson (skip near-zero values).
  ! ------------------------------------------------------------------
  subroutine poisson_add_coo(rows, cols, vals, capacity, idx, row, col, val)
    integer, intent(inout) :: rows(:), cols(:)
    real(kind=dp), intent(inout) :: vals(:)
    integer, intent(in) :: capacity
    integer, intent(inout) :: idx
    integer, intent(in) :: row, col
    real(kind=dp), intent(in) :: val

    if (abs(val) < 1.0e-30_dp) return

    idx = idx + 1
    if (idx > capacity) then
      print *, 'ERROR: COO capacity exceeded in 2D Poisson solver'
      stop 1
    end if
    rows(idx) = row
    cols(idx) = col
    vals(idx) = val

  end subroutine poisson_add_coo


  ! ------------------------------------------------------------------
  ! Convert COO to CSR format (real-valued, 1-based indexing).
  ! Sorts by (row, col), merges duplicates, produces ia/ja/a_csr
  ! compatible with PARDISO (ia has n+1 entries, 1-based).
  ! ------------------------------------------------------------------
  subroutine poisson_coo_to_csr(nrows, ncols, nnz_in, rows, cols, vals, &
      ia, ja, a_csr)

    integer, intent(in) :: nrows, ncols, nnz_in
    integer, intent(in) :: rows(nnz_in), cols(nnz_in)
    real(kind=dp), intent(in) :: vals(nnz_in)
    integer, allocatable, intent(out) :: ia(:), ja(:)
    real(kind=dp), allocatable, intent(out) :: a_csr(:)

    integer, allocatable :: idx(:), r_sorted(:), c_sorted(:)
    real(kind=dp), allocatable :: v_sorted(:)
    integer :: i, nnz_final, row

    if (nnz_in == 0) then
      allocate(ia(nrows + 1), ja(0), a_csr(0))
      ia = 1
      return
    end if

    ! Sort COO entries by (row, col)
    allocate(idx(nnz_in))
    do i = 1, nnz_in
      idx(i) = i
    end do
    call poisson_merge_sort_coo(nnz_in, idx, rows, cols)

    ! Apply permutation and merge duplicates
    allocate(v_sorted(nnz_in), r_sorted(nnz_in), c_sorted(nnz_in))

    v_sorted(1) = vals(idx(1))
    r_sorted(1) = rows(idx(1))
    c_sorted(1) = cols(idx(1))
    nnz_final = 1

    do i = 2, nnz_in
      if (r_sorted(nnz_final) == rows(idx(i)) .and. &
          c_sorted(nnz_final) == cols(idx(i))) then
        v_sorted(nnz_final) = v_sorted(nnz_final) + vals(idx(i))
      else
        nnz_final = nnz_final + 1
        r_sorted(nnz_final) = rows(idx(i))
        c_sorted(nnz_final) = cols(idx(i))
        v_sorted(nnz_final) = vals(idx(i))
      end if
    end do

    ! Build CSR arrays
    allocate(a_csr(nnz_final), ja(nnz_final), ia(nrows + 1))
    a_csr(1:nnz_final) = v_sorted(1:nnz_final)
    ja(1:nnz_final) = c_sorted(1:nnz_final)

    ! Build rowptr
    ia = nnz_final + 1  ! sentinel
    do i = 1, nnz_final
      row = r_sorted(i)
      if (ia(row) > i) ia(row) = i
    end do
    ia(nrows + 1) = nnz_final + 1

    ! Fill empty rows
    do row = nrows, 1, -1
      if (ia(row) == nnz_final + 1) then
        ia(row) = ia(row + 1)
      end if
    end do

    deallocate(idx, v_sorted, r_sorted, c_sorted)

  end subroutine poisson_coo_to_csr


  ! ------------------------------------------------------------------
  ! Bottom-up merge sort for real COO by (row, col)
  ! ------------------------------------------------------------------
  subroutine poisson_merge_sort_coo(nnz, idx, rows, cols)
    integer, intent(in) :: nnz
    integer, intent(inout) :: idx(nnz)
    integer, intent(in) :: rows(nnz), cols(nnz)

    integer, allocatable :: work(:)
    integer :: width, start, mid, finish, i, j, k

    if (nnz <= 1) return

    allocate(work(nnz))

    width = 1
    do while (width < nnz)
      start = 1
      do while (start <= nnz)
        mid = min(start + width - 1, nnz)
        finish = min(start + 2*width - 1, nnz)

        if (mid < finish) then
          i = start
          j = mid + 1
          k = start
          do while (i <= mid .and. j <= finish)
            if (rows(idx(i)) < rows(idx(j))) then
              work(k) = idx(i); i = i + 1
            else if (rows(idx(i)) > rows(idx(j))) then
              work(k) = idx(j); j = j + 1
            else
              if (cols(idx(i)) <= cols(idx(j))) then
                work(k) = idx(i); i = i + 1
              else
                work(k) = idx(j); j = j + 1
              end if
            end if
            k = k + 1
          end do

          do while (i <= mid)
            work(k) = idx(i); i = i + 1; k = k + 1
          end do
          do while (j <= finish)
            work(k) = idx(j); j = j + 1; k = k + 1
          end do

          idx(start:finish) = work(start:finish)
        end if

        start = start + 2*width
      end do
      width = width * 2
    end do

    deallocate(work)

  end subroutine poisson_merge_sort_coo

end module poisson
