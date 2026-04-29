module finitedifferences

  use definitions

  implicit none
  private
  public :: FDmatrixDense, Identity, FDstencil, toeplitz
  public :: buildFD2ndDerivMatrix, buildFD1stDerivMatrix
  public :: FDcentralCoeffs2nd, FDcentralCoeffs1st
  public :: FDforwardCoeffs2nd, FDbackwardCoeffs2nd
  public :: FDforwardCoeffs1st, FDbackwardCoeffs1st
  public :: vandermonde_2nd_deriv, vandermonde_1st_deriv
  public :: buildStaggeredD1Inner, buildStaggeredD1Outer
  public :: interpolateToHalfPoints

  contains

    subroutine FDmatrixDense(nz, dz, order, stencil, boundary, matrix)

      integer, intent(in) :: nz, stencil, boundary, order
      real(kind=dp), intent(in) :: dz
      real(kind=dp), intent(inout), allocatable, dimension(:,:) :: matrix
      real(kind=dp), allocatable, dimension(:) :: v1


      allocate(v1(nz))
      v1 = 0

      call FDstencil(order, stencil, dz, v1)
      !boundary = 0 -> hard Wall
      !boundary = 1 -> periodic
      call toeplitz(v1, boundary, matrix)
      deallocate(v1)

    end subroutine FDmatrixDense

    subroutine FDstencil(order, stencil, d, vector)

      integer, intent(in) :: order, stencil
      real(kind=dp), intent(in) :: d
      real(kind=dp), intent(inout), allocatable, dimension(:) :: vector
      integer :: length, hf

      length = size(vector, dim=1)

      IF (ALLOCATED(vector)) DEALLOCATE(vector)
      ALLOCATE(vector(length))

      hf = int(length/2. + 0.5)

      if (order == 1) then

        select case(stencil)

        case(2)

          if (length .lt. 3) stop 'for stencil 2, n should be 3 or more'
          vector(hf-1) = -1.0_dp
          vector(hf+1) = 1.0_dp
          vector = vector/(2.0_dp*d)

        case(4)
          print *, 'case 4 not implemented yet'

        case(6)
          print *, 'case 6 not implemented yet'

        case(8)
          if (length .lt. 9) stop 'for stencil 8, n should be 9 or more'
          vector(hf-4) = 3.0_dp
          vector(hf-3) = -32.0_dp
          vector(hf-2) = 168.0_dp
          vector(hf-1) = -672.0_dp
          vector(hf) = 0.0_dp
          vector(hf+1) = 672.0_dp
          vector(hf+2) = -168.0_dp
          vector(hf+3) = 32.0_dp
          vector(hf+4) = -3.0_dp
          vector = vector/(840.0_dp*d)

        case default
          print *, 'not implemented yet'

        end select


      end if


      if (order == 2) then

        if (stencil == 2) then
          if (length .lt. 3) stop 'for stencil 2, n should be 3 or more'
          vector(hf-1) = 1.0_dp
          vector(hf) = -2.0_dp
          vector(hf+1) = 1.0_dp
        endif

        if (stencil == 4) then
          if (length .lt. 5) stop 'for stencil 4, n should be 5 or more'
          vector(hf-2) = -1.0_dp/12.0_dp
          vector(hf-1) = 4.0_dp/3.0_dp
          vector(hf) = -5.0_dp/2.0_dp
          vector(hf+1) = 4.0_dp/3.0_dp
          vector(hf+2) = -1.0_dp/12.0_dp
        endif

        if (stencil == 6) then
          if (length .lt. 7) stop 'for stencil 6, n should be 7 or more'
          vector(hf-3) = 1.0_dp/90.0_dp
          vector(hf-2) = -3.0_dp/20.0_dp
          vector(hf-1) = 3.0_dp/2.0_dp
          vector(hf) = -49.0_dp/18.0_dp
          vector(hf+1) = 3.0_dp/2.0_dp
          vector(hf+2) = -3.0_dp/20.0_dp
          vector(hf+3) = 1.0_dp/90.0_dp
        endif

        if (stencil == 8) then
          if (length .lt. 9) stop 'for stencil 8, n should be 9 or more'
          vector(hf-4) = -1.0_dp/560.0_dp
          vector(hf-3) = 8.0_dp/315.0_dp
          vector(hf-2) = -1.0_dp/5.0_dp
          vector(hf-1) = 8.0_dp/5.0_dp
          vector(hf) = -205.0_dp/72.0_dp
          vector(hf+1) = 8.0_dp/5.0_dp
          vector(hf+2) = -1.0_dp/5.0_dp
          vector(hf+3) = 8.0_dp/315.0_dp
          vector(hf+4) = -1.0_dp/560.0_dp
        endif

        if (stencil == 10) then
          if (length .lt. 11) stop 'for stencil 10, n should be 11 or more'
          vector(hf-5) = 8.0_dp/25200.0_dp
          vector(hf-4) = -125.0_dp/25200.0_dp
          vector(hf-3) = 1000.0_dp/25200.0_dp
          vector(hf-2) = -6000.0_dp/25200.0_dp
          vector(hf-1) = 42000.0_dp/25200.0_dp
          vector(hf) = -73766.0_dp/25200.0_dp
          vector(hf+1) = 42000.0_dp/25200.0_dp
          vector(hf+2) = -6000.0_dp/25200.0_dp
          vector(hf+3) = 1000.0_dp/25200.0_dp
          vector(hf+4) = -125.0_dp/25200.0_dp
          vector(hf+5) = 8.0_dp/25200.0_dp
        endif

        vector = vector/d**2

      end if

    end subroutine FDstencil

    !---------------------------------------------------------------------------
    !> Build a full NxN FD matrix for the 2nd derivative with proper boundary
    !> handling using one-sided stencils at the boundary regions.
    !>
    !> For interior points (away from boundaries by at least half_bandwidth),
    !> the central stencil is used. For points near boundaries, one-sided
    !> forward/backward stencils provide the specified accuracy.
    !>
    !> @param[in]  N        Number of grid points
    !> @param[in]  dz       Grid spacing
    !> @param[in]  FDorder  FD accuracy order (2, 4, 6, 8, or 10)
    !> @param[out] D2       N x N 2nd-derivative FD matrix
    !---------------------------------------------------------------------------
    subroutine buildFD2ndDerivMatrix(N, dz, FDorder, D2)

      integer, intent(in) :: N, FDorder
      real(kind=dp), intent(in) :: dz
      real(kind=dp), intent(inout), allocatable, dimension(:,:) :: D2

      integer :: half_bw, i, j, npts
      real(kind=dp), allocatable, dimension(:) :: coeffs

      half_bw = FDorder / 2

      if (allocated(D2)) deallocate(D2)
      allocate(D2(N, N))
      D2 = 0.0_dp

      ! Interior points: use central stencil
      do i = half_bw + 1, N - half_bw
        call FDcentralCoeffs2nd(FDorder, dz, coeffs)
        do j = -half_bw, half_bw
          D2(i, i + j) = coeffs(j + half_bw + 1)
        end do
      end do

      ! Boundary points at the left: use forward one-sided stencils
      do i = 1, min(half_bw, N)
        call FDforwardCoeffs2nd(i, FDorder, dz, coeffs)
        npts = size(coeffs)
        do j = 1, npts
          if (j <= N) D2(i, j) = coeffs(j)
        end do
      end do

      ! Boundary points at the right: use backward one-sided stencils
      do i = max(N - half_bw + 1, half_bw + 1), N
        call FDbackwardCoeffs2nd(N - i + 1, FDorder, dz, coeffs)
        npts = size(coeffs)
        do j = 1, npts
          if (N - npts + j >= 1) D2(i, N - npts + j) = coeffs(j)
        end do
      end do

      if (allocated(coeffs)) deallocate(coeffs)

    end subroutine buildFD2ndDerivMatrix

    !---------------------------------------------------------------------------
    !> Build a full NxN FD matrix for the 1st derivative with proper boundary
    !> handling using one-sided stencils at the boundary regions.
    !>
    !> @param[in]  N        Number of grid points
    !> @param[in]  dz       Grid spacing
    !> @param[in]  FDorder  FD accuracy order (2, 4, 6, 8, or 10)
    !> @param[out] D1       N x N 1st-derivative FD matrix
    !---------------------------------------------------------------------------
    subroutine buildFD1stDerivMatrix(N, dz, FDorder, D1)

      integer, intent(in) :: N, FDorder
      real(kind=dp), intent(in) :: dz
      real(kind=dp), intent(inout), allocatable, dimension(:,:) :: D1

      integer :: half_bw, i, j, npts
      real(kind=dp), allocatable, dimension(:) :: coeffs

      half_bw = FDorder / 2

      if (allocated(D1)) deallocate(D1)
      allocate(D1(N, N))
      D1 = 0.0_dp

      ! Interior points: use central stencil for 1st derivative
      do i = half_bw + 1, N - half_bw
        call FDcentralCoeffs1st(FDorder, dz, coeffs)
        do j = -half_bw, half_bw
          D1(i, i + j) = coeffs(j + half_bw + 1)
        end do
      end do

      ! Boundary points at the left: use forward one-sided stencils
      do i = 1, min(half_bw, N)
        call FDforwardCoeffs1st(i, FDorder, dz, coeffs)
        npts = size(coeffs)
        do j = 1, npts
          if (j <= N) D1(i, j) = coeffs(j)
        end do
      end do

      ! Boundary points at the right: use backward one-sided stencils
      do i = max(N - half_bw + 1, half_bw + 1), N
        call FDbackwardCoeffs1st(N - i + 1, FDorder, dz, coeffs)
        npts = size(coeffs)
        do j = 1, npts
          if (N - npts + j >= 1) D1(i, N - npts + j) = coeffs(j)
        end do
      end do

      if (allocated(coeffs)) deallocate(coeffs)

    end subroutine buildFD1stDerivMatrix

    !---------------------------------------------------------------------------
    !> Central stencil coefficients for 2nd derivative.
    !> Returns coefficients for points at offsets -half_bw, ..., 0, ..., +half_bw.
    !---------------------------------------------------------------------------
    subroutine FDcentralCoeffs2nd(FDorder, d, coeffs)

      integer, intent(in) :: FDorder
      real(kind=dp), intent(in) :: d
      real(kind=dp), allocatable, dimension(:), intent(out) :: coeffs

      integer :: npts, hf

      npts = FDorder + 1

      if (allocated(coeffs)) deallocate(coeffs)
      allocate(coeffs(npts))
      coeffs = 0.0_dp

      hf = npts / 2 + 1  ! center index (1-based)

      select case(FDorder)
      case(2)
        coeffs(hf-1) = 1.0_dp
        coeffs(hf)   = -2.0_dp
        coeffs(hf+1) = 1.0_dp
      case(4)
        coeffs(hf-2) = -1.0_dp/12.0_dp
        coeffs(hf-1) = 4.0_dp/3.0_dp
        coeffs(hf)   = -5.0_dp/2.0_dp
        coeffs(hf+1) = 4.0_dp/3.0_dp
        coeffs(hf+2) = -1.0_dp/12.0_dp
      case(6)
        coeffs(hf-3) = 1.0_dp/90.0_dp
        coeffs(hf-2) = -3.0_dp/20.0_dp
        coeffs(hf-1) = 3.0_dp/2.0_dp
        coeffs(hf)   = -49.0_dp/18.0_dp
        coeffs(hf+1) = 3.0_dp/2.0_dp
        coeffs(hf+2) = -3.0_dp/20.0_dp
        coeffs(hf+3) = 1.0_dp/90.0_dp
      case(8)
        coeffs(hf-4) = -1.0_dp/560.0_dp
        coeffs(hf-3) = 8.0_dp/315.0_dp
        coeffs(hf-2) = -1.0_dp/5.0_dp
        coeffs(hf-1) = 8.0_dp/5.0_dp
        coeffs(hf)   = -205.0_dp/72.0_dp
        coeffs(hf+1) = 8.0_dp/5.0_dp
        coeffs(hf+2) = -1.0_dp/5.0_dp
        coeffs(hf+3) = 8.0_dp/315.0_dp
        coeffs(hf+4) = -1.0_dp/560.0_dp
      case(10)
        coeffs(hf-5) = 8.0_dp/25200.0_dp
        coeffs(hf-4) = -125.0_dp/25200.0_dp
        coeffs(hf-3) = 1000.0_dp/25200.0_dp
        coeffs(hf-2) = -6000.0_dp/25200.0_dp
        coeffs(hf-1) = 42000.0_dp/25200.0_dp
        coeffs(hf)   = -73766.0_dp/25200.0_dp
        coeffs(hf+1) = 42000.0_dp/25200.0_dp
        coeffs(hf+2) = -6000.0_dp/25200.0_dp
        coeffs(hf+3) = 1000.0_dp/25200.0_dp
        coeffs(hf+4) = -125.0_dp/25200.0_dp
        coeffs(hf+5) = 8.0_dp/25200.0_dp
      case default
        print *, 'Error: unsupported FDorder:', FDorder
        stop 1
      end select

      coeffs = coeffs / d**2

    end subroutine FDcentralCoeffs2nd

    !---------------------------------------------------------------------------
    !> Central stencil coefficients for 1st derivative.
    !---------------------------------------------------------------------------
    subroutine FDcentralCoeffs1st(FDorder, d, coeffs)

      integer, intent(in) :: FDorder
      real(kind=dp), intent(in) :: d
      real(kind=dp), allocatable, dimension(:), intent(out) :: coeffs

      integer :: npts, hf

      npts = FDorder + 1

      if (allocated(coeffs)) deallocate(coeffs)
      allocate(coeffs(npts))
      coeffs = 0.0_dp

      hf = npts / 2 + 1

      select case(FDorder)
      case(2)
        coeffs(hf-1) = -1.0_dp/2.0_dp
        coeffs(hf+1) = 1.0_dp/2.0_dp
      case(4)
        coeffs(hf-2) = 1.0_dp/12.0_dp
        coeffs(hf-1) = -2.0_dp/3.0_dp
        coeffs(hf+1) = 2.0_dp/3.0_dp
        coeffs(hf+2) = -1.0_dp/12.0_dp
      case(6)
        coeffs(hf-3) = -1.0_dp/60.0_dp
        coeffs(hf-2) = 3.0_dp/20.0_dp
        coeffs(hf-1) = -3.0_dp/4.0_dp
        coeffs(hf+1) = 3.0_dp/4.0_dp
        coeffs(hf+2) = -3.0_dp/20.0_dp
        coeffs(hf+3) = 1.0_dp/60.0_dp
      case(8)
        coeffs(hf-4) = 1.0_dp/280.0_dp
        coeffs(hf-3) = -4.0_dp/105.0_dp
        coeffs(hf-2) = 1.0_dp/5.0_dp
        coeffs(hf-1) = -4.0_dp/5.0_dp
        coeffs(hf+1) = 4.0_dp/5.0_dp
        coeffs(hf+2) = -1.0_dp/5.0_dp
        coeffs(hf+3) = 4.0_dp/105.0_dp
        coeffs(hf+4) = -1.0_dp/280.0_dp
      case(10)
        coeffs(hf-5) = -1.0_dp/1260.0_dp
        coeffs(hf-4) = 5.0_dp/1008.0_dp
        coeffs(hf-3) = -5.0_dp/126.0_dp
        coeffs(hf-2) = 5.0_dp/21.0_dp
        coeffs(hf-1) = -5.0_dp/3.0_dp
        coeffs(hf+1) = 5.0_dp/3.0_dp
        coeffs(hf+2) = -5.0_dp/21.0_dp
        coeffs(hf+3) = 5.0_dp/126.0_dp
        coeffs(hf+4) = -5.0_dp/1008.0_dp
        coeffs(hf+5) = 1.0_dp/1260.0_dp
      case default
        print *, 'Error: unsupported FDorder:', FDorder
        stop 1
      end select

      coeffs = coeffs / d

    end subroutine FDcentralCoeffs1st

    !---------------------------------------------------------------------------
    !> Forward (one-sided) 2nd-derivative stencil coefficients.
    !> Evaluated at point 0, using points 0, 1, ..., npts-1.
    !> boundary_idx is the distance from the boundary (1 = boundary point).
    !> For interior accuracy at boundaries, uses fewer points when closer.
    !---------------------------------------------------------------------------
    subroutine FDforwardCoeffs2nd(boundary_idx, FDorder, d, coeffs)

      integer, intent(in) :: boundary_idx, FDorder
      real(kind=dp), intent(in) :: d
      real(kind=dp), allocatable, dimension(:), intent(out) :: coeffs

      integer :: npts, j
      real(kind=dp), allocatable :: offsets(:)

      npts = FDorder + 1
      if (allocated(coeffs)) deallocate(coeffs)
      allocate(coeffs(npts))
      coeffs = 0.0_dp

      ! Offsets from evaluation point (boundary_idx, 1-indexed) to
      ! stencil points 1..npts (1-indexed grid positions).
      allocate(offsets(npts))
      do j = 1, npts
        offsets(j) = dble(j - boundary_idx)
      end do

      call vandermonde_2nd_deriv(offsets, npts, coeffs)

      deallocate(offsets)
      coeffs = coeffs / d**2

    end subroutine FDforwardCoeffs2nd

    !---------------------------------------------------------------------------
    !> Backward (one-sided) 2nd-derivative stencil coefficients.
    !> Evaluated at point 0, using points 0, -1, ..., -(npts-1).
    !---------------------------------------------------------------------------
    subroutine FDbackwardCoeffs2nd(boundary_idx, FDorder, d, coeffs)

      integer, intent(in) :: boundary_idx, FDorder
      real(kind=dp), intent(in) :: d
      real(kind=dp), allocatable, dimension(:), intent(out) :: coeffs

      integer :: npts, j
      real(kind=dp), allocatable :: offsets(:)

      npts = FDorder + 1
      if (allocated(coeffs)) deallocate(coeffs)
      allocate(coeffs(npts))
      coeffs = 0.0_dp

      ! Backward offsets: from evaluation point (at position N-boundary_idx+1)
      ! to stencil points {N-npts+1, ..., N}.  Offsets are negative/zero.
      allocate(offsets(npts))
      do j = 1, npts
        offsets(j) = dble(-(npts - j) + (boundary_idx - 1))
      end do

      call vandermonde_2nd_deriv(offsets, npts, coeffs)

      deallocate(offsets)
      coeffs = coeffs / d**2

    end subroutine FDbackwardCoeffs2nd

    !---------------------------------------------------------------------------
    !> Solve Vandermonde system for 2nd-derivative FD coefficients.
    !> Given n offsets (in units of dz) from the evaluation point,
    !> computes coefficients c such that sum_j c_j*f(x+p_j*dz) = f''(x).
    !> Uses Gaussian elimination with partial pivoting (n <= 11).
    !---------------------------------------------------------------------------
    subroutine vandermonde_2nd_deriv(offsets, n, coeffs)
      integer, intent(in) :: n
      real(kind=dp), intent(in) :: offsets(n)
      real(kind=dp), intent(out) :: coeffs(n)

      real(kind=dp) :: V(11, 11), rhs(11), factor
      integer :: i, j, k, pivot_row

      ! Build Vandermonde matrix
      do j = 1, n
        do i = 1, n
          V(i, j) = offsets(j) ** (i - 1)
        end do
      end do
      rhs(1:n) = 0.0_dp
      rhs(3) = 2.0_dp  ! 2nd derivative: sum c_j * p_j^2 = 2

      ! Gaussian elimination with partial pivoting
      do k = 1, n
        pivot_row = k
        do i = k + 1, n
          if (abs(V(i, k)) > abs(V(pivot_row, k))) pivot_row = i
        end do
        if (pivot_row /= k) then
          do j = k, n
            factor = V(k, j); V(k, j) = V(pivot_row, j); V(pivot_row, j) = factor
          end do
          factor = rhs(k); rhs(k) = rhs(pivot_row); rhs(pivot_row) = factor
        end if
        do i = k + 1, n
          factor = V(i, k) / V(k, k)
          do j = k + 1, n
            V(i, j) = V(i, j) - factor * V(k, j)
          end do
          rhs(i) = rhs(i) - factor * rhs(k)
        end do
      end do

      ! Back-substitution
      do i = n, 1, -1
        coeffs(i) = rhs(i)
        do j = i + 1, n
          coeffs(i) = coeffs(i) - V(i, j) * coeffs(j)
        end do
        coeffs(i) = coeffs(i) / V(i, i)
      end do

    end subroutine vandermonde_2nd_deriv

    !---------------------------------------------------------------------------
    !> Solve Vandermonde system for 1st-derivative FD coefficients.
    !> Given n offsets (in units of dz) from the evaluation point,
    !> computes coefficients c such that sum_j c_j*f(x+p_j*dz) = f'(x).
    !> Uses Gaussian elimination with partial pivoting (n <= 11).
    !---------------------------------------------------------------------------
    subroutine vandermonde_1st_deriv(offsets, n, coeffs)
      integer, intent(in) :: n
      real(kind=dp), intent(in) :: offsets(n)
      real(kind=dp), intent(out) :: coeffs(n)

      real(kind=dp) :: V(11, 11), rhs(11), factor
      integer :: i, j, k, pivot_row

      ! Build Vandermonde matrix
      do j = 1, n
        do i = 1, n
          V(i, j) = offsets(j) ** (i - 1)
        end do
      end do
      rhs(1:n) = 0.0_dp
      rhs(2) = 1.0_dp  ! 1st derivative: sum c_j * p_j = 1

      ! Gaussian elimination with partial pivoting
      do k = 1, n
        pivot_row = k
        do i = k + 1, n
          if (abs(V(i, k)) > abs(V(pivot_row, k))) pivot_row = i
        end do
        if (pivot_row /= k) then
          do j = k, n
            factor = V(k, j); V(k, j) = V(pivot_row, j); V(pivot_row, j) = factor
          end do
          factor = rhs(k); rhs(k) = rhs(pivot_row); rhs(pivot_row) = factor
        end if
        do i = k + 1, n
          factor = V(i, k) / V(k, k)
          do j = k + 1, n
            V(i, j) = V(i, j) - factor * V(k, j)
          end do
          rhs(i) = rhs(i) - factor * rhs(k)
        end do
      end do

      ! Back-substitution
      do i = n, 1, -1
        coeffs(i) = rhs(i)
        do j = i + 1, n
          coeffs(i) = coeffs(i) - V(i, j) * coeffs(j)
        end do
        coeffs(i) = coeffs(i) / V(i, i)
      end do

    end subroutine vandermonde_1st_deriv

    !---------------------------------------------------------------------------
    !> Forward (one-sided) 1st-derivative stencil coefficients.
    !---------------------------------------------------------------------------
    subroutine FDforwardCoeffs1st(boundary_idx, FDorder, d, coeffs)

      integer, intent(in) :: boundary_idx, FDorder
      real(kind=dp), intent(in) :: d
      real(kind=dp), allocatable, dimension(:), intent(out) :: coeffs

      integer :: npts, j
      real(kind=dp), allocatable :: offsets(:)

      npts = FDorder + 1
      if (allocated(coeffs)) deallocate(coeffs)
      allocate(coeffs(npts))
      coeffs = 0.0_dp

      allocate(offsets(npts))
      do j = 1, npts
        offsets(j) = dble(j - boundary_idx)
      end do

      call vandermonde_1st_deriv(offsets, npts, coeffs)

      deallocate(offsets)
      coeffs = coeffs / d

    end subroutine FDforwardCoeffs1st

    !---------------------------------------------------------------------------
    !> Backward (one-sided) 1st-derivative stencil coefficients.
    !> Uses Vandermonde system for consistency with forward coefficients.
    !---------------------------------------------------------------------------
    subroutine FDbackwardCoeffs1st(boundary_idx, FDorder, d, coeffs)

      integer, intent(in) :: boundary_idx, FDorder
      real(kind=dp), intent(in) :: d
      real(kind=dp), allocatable, dimension(:), intent(out) :: coeffs

      integer :: npts, j
      real(kind=dp), allocatable :: offsets(:)

      npts = FDorder + 1

      if (allocated(coeffs)) deallocate(coeffs)
      allocate(coeffs(npts))
      coeffs = 0.0_dp

      ! Backward offsets: from evaluation point (at position N-boundary_idx+1)
      ! to stencil points {N-npts+1, ..., N}.  Offsets are negative/zero.
      allocate(offsets(npts))
      do j = 1, npts
        offsets(j) = dble(-(npts - j) + (boundary_idx - 1))
      end do

      call vandermonde_1st_deriv(offsets, npts, coeffs)

      deallocate(offsets)
      coeffs = coeffs / d

    end subroutine FDbackwardCoeffs1st

    subroutine Identity(n,matrix,factor)

      integer, intent(in) :: n
      real(kind=dp), intent(inout), allocatable, dimension(:,:) :: matrix
      real(kind=dp), intent(in), optional :: factor
      integer :: i

      IF (ALLOCATED(matrix)) DEALLOCATE(matrix)
      ALLOCATE(matrix(n,n))

      matrix = 0
      if (present(factor)) then
        do i = 1, n
          matrix(i,i) = factor
        end do
      else
        do i = 1, n
          matrix(i,i) = 1
        end do
      endif

    end subroutine

    subroutine toeplitz(vector, boundary, matrix)

      real(kind=dp), intent(in) :: vector(:)
      integer, intent(in) :: boundary
      real(kind=dp), intent(inout), allocatable :: matrix(:,:)
      integer :: length, idx, hf

      length = size(vector, dim=1)

      IF (ALLOCATED(matrix)) DEALLOCATE(matrix)
      ALLOCATE(matrix(length,length))

      hf = int(length/2. + 0.5)

      if (boundary == 0) then
        do idx = 1, length, 1
          matrix(:,idx) = eoshift(vector,shift=hf-idx)
        end do
      endif

      if (boundary == 1) then
        do idx = 1, length, 1
          matrix(:,idx) = cshift(vector,shift=hf-idx)
        end do
      end if

    end subroutine toeplitz

    !---------------------------------------------------------------------------
    !> Build (N-1) x N half-point forward first-derivative matrix.
    !> Row j computes df/dz at the half-point z_{j+1/2} = (z_j + z_{j+1})/2.
    !>
    !> For FDorder=2, this gives the standard 2-point stencil:
    !>   D_inner(j, j) = -1/dz,  D_inner(j, j+1) = 1/dz
    !>
    !> For FDorder>=4, interior half-points use FDorder-point symmetric stencils
    !> with half-integer offsets. Boundary half-points use (FDorder+1)-point
    !> one-sided stencils to maintain accuracy.
    !---------------------------------------------------------------------------
    subroutine buildStaggeredD1Inner(N, dz, FDorder, D_inner)
      integer, intent(in) :: N, FDorder
      real(kind=dp), intent(in) :: dz
      real(kind=dp), allocatable, intent(out) :: D_inner(:,:)

      integer :: half_bw, i, j, npts, col, bw_inner
      real(kind=dp) :: coeffs(11), offsets_work(11)

      half_bw = FDorder / 2

      if (allocated(D_inner)) deallocate(D_inner)
      allocate(D_inner(N-1, N))
      D_inner = 0.0_dp

      ! For FDorder=2: simple 2-point forward difference at each half-point
      if (FDorder == 2) then
        do j = 1, N - 1
          D_inner(j, j)   = -1.0_dp / dz
          D_inner(j, j+1) =  1.0_dp / dz
        end do
        return
      end if

      ! For FDorder >= 4: use (FDorder+1)-point stencils for FDorder-accurate
      ! 1st derivative at half-points.
      ! Interior: symmetric centered stencil with half-integer offsets.
      ! E.g. FDorder=4: 5-point stencil, offsets = -2.5,-1.5,-0.5,0.5,1.5
      !      covering grid points j-2,j-1,j,j+1,j+2.
      bw_inner = (FDorder + 1) / 2  ! stencil half-width for interior

      do j = 1, N - 1

        if (j >= bw_inner .and. j <= N - bw_inner) then
          ! Interior half-points: (FDorder+1)-point symmetric central stencil
          npts = FDorder + 1
          do i = 1, npts
            offsets_work(i) = dble(i - bw_inner) - 0.5_dp
          end do
          call vandermonde_1st_deriv(offsets_work(1:npts), npts, coeffs(1:npts))

          do i = 1, npts
            col = j + nint(0.5_dp + offsets_work(i))
            if (col >= 1 .and. col <= N) D_inner(j, col) = coeffs(i) / dz
          end do

        else if (j < bw_inner) then
          ! Left boundary: one-sided forward stencil
          ! Evaluation point: z_{j+0.5}.  Use grid points 1..npts.
          npts = FDorder + 1
          do i = 1, npts
            offsets_work(i) = dble(i) - dble(j) - 0.5_dp
          end do
          call vandermonde_1st_deriv(offsets_work(1:npts), npts, coeffs(1:npts))

          do i = 1, npts
            col = i
            if (col >= 1 .and. col <= N) D_inner(j, col) = coeffs(i) / dz
          end do

        else
          ! Right boundary: one-sided backward stencil
          npts = FDorder + 1
          do i = 1, npts
            offsets_work(i) = dble(N - npts + i) - dble(j) - 0.5_dp
          end do
          call vandermonde_1st_deriv(offsets_work(1:npts), npts, coeffs(1:npts))

          do i = 1, npts
            col = N - npts + i
            if (col >= 1 .and. col <= N) D_inner(j, col) = coeffs(i) / dz
          end do

        end if

      end do

    end subroutine buildStaggeredD1Inner

    !---------------------------------------------------------------------------
    !> Build N x (N-1) half-point backward first-derivative matrix.
    !> Row i computes dg/dz at grid point z_i from half-point values.
    !>
    !> Uses D_outer = -D_inner^T to ensure the product D_outer*D_inner is
    !> symmetric (self-adjoint), which is required for the conservative
    !> variable-coefficient discretization of d/dz[g(z)*d/dz].
    !---------------------------------------------------------------------------
    subroutine buildStaggeredD1Outer(N, dz, FDorder, D_outer)
      integer, intent(in) :: N, FDorder
      real(kind=dp), intent(in) :: dz
      real(kind=dp), allocatable, intent(out) :: D_outer(:,:)

      real(kind=dp), allocatable :: D_inner(:,:)
      integer :: i, j

      ! Build D_inner first, then set D_outer = -D_inner^T
      call buildStaggeredD1Inner(N, dz, FDorder, D_inner)

      if (allocated(D_outer)) deallocate(D_outer)
      allocate(D_outer(N, N-1))
      D_outer = 0.0_dp

      ! D_outer(i,j) = -D_inner(j,i)
      do j = 1, N - 1
        do i = 1, N
          D_outer(i, j) = -D_inner(j, i)
        end do
      end do

      deallocate(D_inner)

    end subroutine buildStaggeredD1Outer

    !---------------------------------------------------------------------------
    !> Interpolate profile values to half-points with order-appropriate accuracy.
    !> g_half(j) = interpolated value at z_{j+1/2} = (z_j + z_{j+1})/2.
    !>
    !> For FDorder=2:  g_half(j) = (g_j + g_{j+1}) / 2  (2nd-order)
    !> For FDorder=4:  interior uses 4th-order Lagrange interpolation at x=1/2:
    !>   g_half(j) = (-g_{j-1} + 9*g_j + 9*g_{j+1} - g_{j+2}) / 16
    !> For boundary half-points: falls back to 2nd-order averaging.
    !---------------------------------------------------------------------------
    subroutine interpolateToHalfPoints(profile_vec, N, FDorder, g_half)
      integer, intent(in) :: N, FDorder
      real(kind=dp), intent(in) :: profile_vec(N)
      real(kind=dp), intent(out) :: g_half(N-1)

      integer :: j

      ! For FDorder=2 or near boundaries: simple 2-point average
      g_half = 0.0_dp

      if (FDorder <= 2) then
        do j = 1, N - 1
          g_half(j) = 0.5_dp * (profile_vec(j) + profile_vec(j + 1))
        end do
        return
      end if

      ! FDorder >= 4: use 4th-order interior formula, 2nd-order at boundaries
      ! Boundary half-points (j=1, j=N-1): use 2-point average
      g_half(1) = 0.5_dp * (profile_vec(1) + profile_vec(2))
      g_half(N-1) = 0.5_dp * (profile_vec(N-1) + profile_vec(N))

      ! Interior half-points: 4th-order formula
      do j = 2, N - 2
        g_half(j) = (-profile_vec(j-1) + 9.0_dp*profile_vec(j) &
          &       + 9.0_dp*profile_vec(j+1) - profile_vec(j+2)) / 16.0_dp
      end do

    end subroutine interpolateToHalfPoints


end module finitedifferences
