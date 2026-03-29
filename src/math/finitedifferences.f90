module finitedifferences

  use definitions

  implicit none

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
          vector(hf-1) = 1
          vector(hf) = -2
          vector(hf+1) = 1
        endif

        if (stencil == 4) then
          if (length .lt. 5) stop 'for stencil 4, n should be 5 or more'
          vector(hf-2) = -1./12.
          vector(hf-1) = 4./3.
          vector(hf) = -5./2.
          vector(hf+1) = 4./3.
          vector(hf+2) = -1./12.
        endif

        if (stencil == 6) then
          if (length .lt. 7) stop 'for stencil 6, n should be 7 or more'
          vector(hf-3) = 1./90.
          vector(hf-2) = -3/20.
          vector(hf-1) = 3./2.
          vector(hf) = -49./18.
          vector(hf+1) = 3./2.
          vector(hf+2) = -3/20.
          vector(hf+3) = 1./90.
        endif

        if (stencil == 8) then
          if (length .lt. 9) stop 'for stencil 8, n should be 9 or more'
          vector(hf-4) = -1./560.
          vector(hf-3) = 8./315.
          vector(hf-2) = -1./5.
          vector(hf-1) = 8./5.
          vector(hf) = -205./72.
          vector(hf+1) = 8./5.
          vector(hf+2) = -1./5.
          vector(hf+3) = 8./315.
          vector(hf+4) = -1./560.
        endif

        if (stencil == 10) then
          if (length .lt. 11) stop 'for stencil 10, n should be 11 or more'
          vector(hf-5) = 8./25200.
          vector(hf-4) = -125./25200.
          vector(hf-3) = 1000./25200.
          vector(hf-2) = -6000./25200.
          vector(hf-1) = 42000./25200.
          vector(hf) = -73766./25200.
          vector(hf+1) = 42000./25200.
          vector(hf+2) = -6000./25200.
          vector(hf+3) = 1000./25200.
          vector(hf+4) = -125./25200.
          vector(hf-5) = 8./25200.
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

      integer :: npts

      ! Number of points in the one-sided stencil: need FDorder+1 points
      ! for FDorder-th order accuracy
      npts = FDorder + 1

      if (allocated(coeffs)) deallocate(coeffs)
      allocate(coeffs(npts))
      coeffs = 0.0_dp

      ! Coefficients for d^2f/dx^2 at point 0 using points 0,1,...,FDorder
      ! These give FDorder-th order accuracy
      select case(FDorder)
      case(2)
        ! 3-point forward 2nd derivative
        coeffs(1) = 1.0_dp
        coeffs(2) = -2.0_dp
        coeffs(3) = 1.0_dp
      case(4)
        ! 5-point forward 2nd derivative
        coeffs(1) = 35.0_dp/12.0_dp
        coeffs(2) = -26.0_dp/3.0_dp
        coeffs(3) = 19.0_dp/2.0_dp
        coeffs(4) = -14.0_dp/3.0_dp
        coeffs(5) = 11.0_dp/12.0_dp
      case(6)
        ! 7-point forward 2nd derivative
        coeffs(1) = 203.0_dp/45.0_dp
        coeffs(2) = -87.0_dp/5.0_dp
        coeffs(3) = 117.0_dp/4.0_dp
        coeffs(4) = -254.0_dp/9.0_dp
        coeffs(5) = 33.0_dp/2.0_dp
        coeffs(6) = -27.0_dp/5.0_dp
        coeffs(7) = 137.0_dp/180.0_dp
      case(8)
        ! 9-point forward 2nd derivative
        coeffs(1) = 29531.0_dp/5040.0_dp
        coeffs(2) = -962.0_dp/35.0_dp
        coeffs(3) = 621.0_dp/10.0_dp
        coeffs(4) = -4006.0_dp/45.0_dp
        coeffs(5) = 691.0_dp/8.0_dp
        coeffs(6) = -282.0_dp/5.0_dp
        coeffs(7) = 2143.0_dp/90.0_dp
        coeffs(8) = -206.0_dp/35.0_dp
        coeffs(9) = 363.0_dp/560.0_dp
      case(10)
        ! 11-point forward 2nd derivative
        coeffs(1)  = 177133.0_dp/25200.0_dp
        coeffs(2)  = -4861.0_dp/126.0_dp
        coeffs(3)  = 6121.0_dp/56.0_dp
        coeffs(4)  = -13082.0_dp/63.0_dp
        coeffs(5)  = 6751.0_dp/24.0_dp
        coeffs(6)  = -6877.0_dp/25.0_dp
        coeffs(7)  = 6961.0_dp/36.0_dp
        coeffs(8)  = -2006.0_dp/21.0_dp
        coeffs(9)  = 3533.0_dp/112.0_dp
        coeffs(10) = -263.0_dp/42.0_dp
        coeffs(11) = 7129.0_dp/12600.0_dp
      case default
        print *, 'Error: unsupported FDorder:', FDorder
        stop 1
      end select

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

      integer :: npts

      npts = FDorder + 1

      if (allocated(coeffs)) deallocate(coeffs)
      allocate(coeffs(npts))
      coeffs = 0.0_dp

      ! Backward stencils are the reverse of forward stencils
      select case(FDorder)
      case(2)
        coeffs(3) = 1.0_dp
        coeffs(2) = -2.0_dp
        coeffs(1) = 1.0_dp
      case(4)
        coeffs(5) = 35.0_dp/12.0_dp
        coeffs(4) = -26.0_dp/3.0_dp
        coeffs(3) = 19.0_dp/2.0_dp
        coeffs(2) = -14.0_dp/3.0_dp
        coeffs(1) = 11.0_dp/12.0_dp
      case(6)
        coeffs(7) = 203.0_dp/45.0_dp
        coeffs(6) = -87.0_dp/5.0_dp
        coeffs(5) = 117.0_dp/4.0_dp
        coeffs(4) = -254.0_dp/9.0_dp
        coeffs(3) = 33.0_dp/2.0_dp
        coeffs(2) = -27.0_dp/5.0_dp
        coeffs(1) = 137.0_dp/180.0_dp
      case(8)
        coeffs(9) = 29531.0_dp/5040.0_dp
        coeffs(8) = -962.0_dp/35.0_dp
        coeffs(7) = 621.0_dp/10.0_dp
        coeffs(6) = -4006.0_dp/45.0_dp
        coeffs(5) = 691.0_dp/8.0_dp
        coeffs(4) = -282.0_dp/5.0_dp
        coeffs(3) = 2143.0_dp/90.0_dp
        coeffs(2) = -206.0_dp/35.0_dp
        coeffs(1) = 363.0_dp/560.0_dp
      case(10)
        coeffs(11) = 177133.0_dp/25200.0_dp
        coeffs(10) = -4861.0_dp/126.0_dp
        coeffs(9)  = 6121.0_dp/56.0_dp
        coeffs(8)  = -13082.0_dp/63.0_dp
        coeffs(7)  = 6751.0_dp/24.0_dp
        coeffs(6)  = -6877.0_dp/25.0_dp
        coeffs(5)  = 6961.0_dp/36.0_dp
        coeffs(4)  = -2006.0_dp/21.0_dp
        coeffs(3)  = 3533.0_dp/112.0_dp
        coeffs(2)  = -263.0_dp/42.0_dp
        coeffs(1)  = 7129.0_dp/12600.0_dp
      case default
        print *, 'Error: unsupported FDorder:', FDorder
        stop 1
      end select

      coeffs = coeffs / d**2

    end subroutine FDbackwardCoeffs2nd

    !---------------------------------------------------------------------------
    !> Forward (one-sided) 1st-derivative stencil coefficients.
    !---------------------------------------------------------------------------
    subroutine FDforwardCoeffs1st(boundary_idx, FDorder, d, coeffs)

      integer, intent(in) :: boundary_idx, FDorder
      real(kind=dp), intent(in) :: d
      real(kind=dp), allocatable, dimension(:), intent(out) :: coeffs

      integer :: npts

      npts = FDorder + 1

      if (allocated(coeffs)) deallocate(coeffs)
      allocate(coeffs(npts))
      coeffs = 0.0_dp

      select case(FDorder)
      case(2)
        coeffs(1) = -3.0_dp/2.0_dp
        coeffs(2) = 2.0_dp
        coeffs(3) = -1.0_dp/2.0_dp
      case(4)
        coeffs(1) = -25.0_dp/12.0_dp
        coeffs(2) = 4.0_dp
        coeffs(3) = -3.0_dp
        coeffs(4) = 4.0_dp/3.0_dp
        coeffs(5) = -1.0_dp/4.0_dp
      case(6)
        coeffs(1) = -49.0_dp/20.0_dp
        coeffs(2) = 6.0_dp
        coeffs(3) = -15.0_dp/2.0_dp
        coeffs(4) = 20.0_dp/3.0_dp
        coeffs(5) = -15.0_dp/4.0_dp
        coeffs(6) = 6.0_dp/5.0_dp
        coeffs(7) = -1.0_dp/6.0_dp
      case(8)
        coeffs(1) = -761.0_dp/280.0_dp
        coeffs(2) = 8.0_dp
        coeffs(3) = -14.0_dp
        coeffs(4) = 56.0_dp/3.0_dp
        coeffs(5) = -35.0_dp/2.0_dp
        coeffs(6) = 56.0_dp/5.0_dp
        coeffs(7) = -14.0_dp/3.0_dp
        coeffs(8) = 8.0_dp/7.0_dp
        coeffs(9) = -1.0_dp/8.0_dp
      case(10)
        coeffs(1)  = -7381.0_dp/2520.0_dp
        coeffs(2)  = 10.0_dp
        coeffs(3)  = -45.0_dp/2.0_dp
        coeffs(4)  = 40.0_dp
        coeffs(5)  = -105.0_dp/2.0_dp
        coeffs(6)  = 252.0_dp/5.0_dp
        coeffs(7)  = -35.0_dp
        coeffs(8)  = 120.0_dp/7.0_dp
        coeffs(9)  = -45.0_dp/8.0_dp
        coeffs(10) = 10.0_dp/9.0_dp
        coeffs(11) = -1.0_dp/10.0_dp
      case default
        print *, 'Error: unsupported FDorder:', FDorder
        stop 1
      end select

      coeffs = coeffs / d

    end subroutine FDforwardCoeffs1st

    !---------------------------------------------------------------------------
    !> Backward (one-sided) 1st-derivative stencil coefficients.
    !---------------------------------------------------------------------------
    subroutine FDbackwardCoeffs1st(boundary_idx, FDorder, d, coeffs)

      integer, intent(in) :: boundary_idx, FDorder
      real(kind=dp), intent(in) :: d
      real(kind=dp), allocatable, dimension(:), intent(out) :: coeffs

      integer :: npts

      npts = FDorder + 1

      if (allocated(coeffs)) deallocate(coeffs)
      allocate(coeffs(npts))
      coeffs = 0.0_dp

      ! Backward = reverse sign and reverse order of forward
      select case(FDorder)
      case(2)
        coeffs(3) = 3.0_dp/2.0_dp
        coeffs(2) = -2.0_dp
        coeffs(1) = 1.0_dp/2.0_dp
      case(4)
        coeffs(5) = 25.0_dp/12.0_dp
        coeffs(4) = -4.0_dp
        coeffs(3) = 3.0_dp
        coeffs(2) = -4.0_dp/3.0_dp
        coeffs(1) = 1.0_dp/4.0_dp
      case(6)
        coeffs(7) = 49.0_dp/20.0_dp
        coeffs(6) = -6.0_dp
        coeffs(5) = 15.0_dp/2.0_dp
        coeffs(4) = -20.0_dp/3.0_dp
        coeffs(3) = 15.0_dp/4.0_dp
        coeffs(2) = -6.0_dp/5.0_dp
        coeffs(1) = 1.0_dp/6.0_dp
      case(8)
        coeffs(9) = 761.0_dp/280.0_dp
        coeffs(8) = -8.0_dp
        coeffs(7) = 14.0_dp
        coeffs(6) = -56.0_dp/3.0_dp
        coeffs(5) = 35.0_dp/2.0_dp
        coeffs(4) = -56.0_dp/5.0_dp
        coeffs(3) = 14.0_dp/3.0_dp
        coeffs(2) = -8.0_dp/7.0_dp
        coeffs(1) = 1.0_dp/8.0_dp
      case(10)
        coeffs(11) = 7381.0_dp/2520.0_dp
        coeffs(10) = -10.0_dp
        coeffs(9)  = 45.0_dp/2.0_dp
        coeffs(8)  = -40.0_dp
        coeffs(7)  = 105.0_dp/2.0_dp
        coeffs(6)  = -252.0_dp/5.0_dp
        coeffs(5)  = 35.0_dp
        coeffs(4)  = -120.0_dp/7.0_dp
        coeffs(3)  = 45.0_dp/8.0_dp
        coeffs(2)  = -10.0_dp/9.0_dp
        coeffs(1)  = 1.0_dp/10.0_dp
      case default
        print *, 'Error: unsupported FDorder:', FDorder
        stop 1
      end select

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
        forall(i = 1:n) matrix(i,i) = 1*factor
      else
        forall(i = 1:n) matrix(i,i) = 1
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



end module finitedifferences
