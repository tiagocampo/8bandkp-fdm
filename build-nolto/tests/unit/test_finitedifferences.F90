module test_finitedifferences
  use funit
  use definitions
  use finitedifferences
  implicit none

contains

  !@test
  subroutine test_fd_stencil_order2()
    ! 2nd-order stencil for d2f/dz2: should give [1, -2, 1]/dz^2
    real(kind=dp), allocatable :: v(:)
    real(kind=dp) :: dz
    dz = 0.1_dp
    allocate(v(3))
    v = 0.0_dp
    call FDstencil(2, 2, dz, v)
#line 18 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(1.0_dp/dz**2, v(1), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 18) )
  if (anyExceptions()) return
#line 19 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
#line 19 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(-2.0_dp/dz**2, v(2), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 19) )
  if (anyExceptions()) return
#line 20 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
#line 20 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(1.0_dp/dz**2, v(3), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 20) )
  if (anyExceptions()) return
#line 21 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  end subroutine test_fd_stencil_order2

  !@test
  subroutine test_fd_stencil_order4()
    ! 4th-order stencil: [-1/12, 4/3, -5/2, 4/3, -1/12] / dz^2
    real(kind=dp), allocatable :: v(:)
    real(kind=dp) :: dz
    dz = 0.1_dp
    allocate(v(5))
    v = 0.0_dp
    call FDstencil(2, 4, dz, v)
    ! Note: source uses single-precision literals for order-4 coefficients
#line 33 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(-1.0_dp/12.0_dp/dz**2, v(1), tolerance=1.0e-5_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 33) )
  if (anyExceptions()) return
#line 34 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
#line 34 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(4.0_dp/3.0_dp/dz**2, v(2), tolerance=1.0e-5_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 34) )
  if (anyExceptions()) return
#line 35 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
#line 35 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(-5.0_dp/2.0_dp/dz**2, v(3), tolerance=1.0e-5_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 35) )
  if (anyExceptions()) return
#line 36 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
#line 36 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(4.0_dp/3.0_dp/dz**2, v(4), tolerance=1.0e-5_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 36) )
  if (anyExceptions()) return
#line 37 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
#line 37 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(-1.0_dp/12.0_dp/dz**2, v(5), tolerance=1.0e-5_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 37) )
  if (anyExceptions()) return
#line 38 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  end subroutine test_fd_stencil_order4

  !@test
  subroutine test_toeplitz_hardwall()
    ! Hard-wall toeplitz from [1, -2, 1] should give tridiagonal [-2,1; 1,-2,1; ...]
    real(kind=dp), allocatable :: v(:), mat(:,:)
    allocate(v(3))
    v(1) = 1.0_dp; v(2) = -2.0_dp; v(3) = 1.0_dp
    call toeplitz(v, 0, mat)
    ! Check symmetry
#line 48 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(mat(1,2), mat(2,1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 48) )
  if (anyExceptions()) return
#line 49 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
#line 49 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(mat(2,3), mat(3,2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 49) )
  if (anyExceptions()) return
#line 50 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
    ! Check diagonal
#line 51 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(-2.0_dp, mat(1,1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 51) )
  if (anyExceptions()) return
#line 52 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
#line 52 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(-2.0_dp, mat(2,2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 52) )
  if (anyExceptions()) return
#line 53 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
    ! Check corners are zero (hard wall)
#line 54 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(0.0_dp, mat(1,3), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 54) )
  if (anyExceptions()) return
#line 55 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
#line 55 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(0.0_dp, mat(3,1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 55) )
  if (anyExceptions()) return
#line 56 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  end subroutine test_toeplitz_hardwall

  !@test
  subroutine test_toeplitz_periodic()
    ! Periodic boundary wraps the stencil around
    real(kind=dp), allocatable :: v(:), mat(:,:)
    allocate(v(3))
    v(1) = 1.0_dp; v(2) = -2.0_dp; v(3) = 1.0_dp
    call toeplitz(v, 1, mat)
    ! With periodic BC, corners should be 1 (wrapping)
#line 66 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(1.0_dp, mat(1,3), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 66) )
  if (anyExceptions()) return
#line 67 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
#line 67 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(1.0_dp, mat(3,1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 67) )
  if (anyExceptions()) return
#line 68 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  end subroutine test_toeplitz_periodic

  !@test
  subroutine test_identity()
    real(kind=dp), allocatable :: mat(:,:)
    integer :: i, n
    n = 5
    call Identity(n, mat)
    do i = 1, n
#line 77 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(1.0_dp, mat(i,i), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 77) )
  if (anyExceptions()) return
#line 78 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
    end do
    ! Off-diagonal should be zero
#line 80 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(0.0_dp, mat(1,2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 80) )
  if (anyExceptions()) return
#line 81 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
#line 81 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(0.0_dp, mat(3,5), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 81) )
  if (anyExceptions()) return
#line 82 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  end subroutine test_identity

  !@test
  subroutine test_identity_scaled()
    real(kind=dp), allocatable :: mat(:,:)
    call Identity(3, mat, factor=2.5_dp)
#line 88 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(2.5_dp, mat(1,1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 88) )
  if (anyExceptions()) return
#line 89 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
#line 89 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(2.5_dp, mat(2,2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
#line 90 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(2.5_dp, mat(3,3), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 90) )
  if (anyExceptions()) return
#line 91 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
#line 91 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(0.0_dp, mat(1,2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 91) )
  if (anyExceptions()) return
#line 92 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  end subroutine test_identity_scaled

  !@test
  subroutine test_buildFD2ndDeriv_order2()
    ! d2/dz2 of f(z) = z^2 should be exactly 2.0 at all interior points
    integer, parameter :: N = 11
    real(kind=dp) :: dz, z(N), f(N)
    real(kind=dp), allocatable :: D2(:,:)
    real(kind=dp), allocatable :: result(:)
    integer :: i

    dz = 0.1_dp
    do i = 1, N
      z(i) = (i - 1) * dz
      f(i) = z(i)**2
    end do

    call buildFD2ndDerivMatrix(N, dz, 2, D2)

    ! D2 * f should be ~2.0 at interior points
    allocate(result(N))
    result = matmul(D2, f)

    ! Check interior points (order 2: only 1 boundary point on each side)
    do i = 2, N-1
#line 117 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(2.0_dp, result(i), tolerance=1.0e-6_dp, message="Interior point d2f/dz2", &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 117) )
  if (anyExceptions()) return
#line 118 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
    end do
  end subroutine test_buildFD2ndDeriv_order2

  !@test
  subroutine test_buildFD2ndDeriv_order4()
    ! d2/dz2 of f(z) = z^3 should be exactly 6*z at all points (order 4 accurate)
    integer, parameter :: N = 21
    real(kind=dp) :: dz, z(N), f(N)
    real(kind=dp), allocatable :: D2(:,:), result(:)
    integer :: i

    dz = 0.05_dp
    do i = 1, N
      z(i) = (i - 1) * dz
      f(i) = z(i)**3
    end do

    call buildFD2ndDerivMatrix(N, dz, 4, D2)

    allocate(result(N))
    result = matmul(D2, f)

    ! Interior points: d2(z^3)/dz2 = 6*z, should be accurate to O(h^4)
    do i = 3, N-2
#line 142 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertEqual(6.0_dp*z(i), result(i), tolerance=1.0e-4_dp, message="Order-4 interior d2f/dz2", &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 142) )
  if (anyExceptions()) return
#line 143 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
    end do
  end subroutine test_buildFD2ndDeriv_order4

  !@test
  subroutine test_central_coeffs_sum_to_zero()
    ! Central stencil coefficients for d2/dz2 should sum to zero
    real(kind=dp), allocatable :: coeffs(:)
    call FDcentralCoeffs2nd(2, 1.0_dp, coeffs)
#line 151 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertTrue(abs(sum(coeffs)) < 1.0e-12_dp, message="Order-2 coeffs sum to zero", &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 151) )
  if (anyExceptions()) return
#line 152 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"

    call FDcentralCoeffs2nd(4, 1.0_dp, coeffs)
#line 154 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  call assertTrue(abs(sum(coeffs)) < 1.0e-12_dp, message="Order-4 coeffs sum to zero", &
 & location=SourceLocation( &
 & 'test_finitedifferences.pf', &
 & 154) )
  if (anyExceptions()) return
#line 155 "/data/8bandkp-fdm/tests/unit/test_finitedifferences.pf"
  end subroutine test_central_coeffs_sum_to_zero

end module test_finitedifferences

module Wraptest_finitedifferences
   use FUnit
   use test_finitedifferences
   implicit none
   private

contains


end module Wraptest_finitedifferences

function test_finitedifferences_suite() result(suite)
   use FUnit
   use test_finitedifferences
   use Wraptest_finitedifferences
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   suite = TestSuite('test_finitedifferences_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_fd_stencil_order2', &
      test_fd_stencil_order2))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_fd_stencil_order4', &
      test_fd_stencil_order4))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_toeplitz_hardwall', &
      test_toeplitz_hardwall))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_toeplitz_periodic', &
      test_toeplitz_periodic))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_identity', &
      test_identity))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_identity_scaled', &
      test_identity_scaled))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_buildFD2ndDeriv_order2', &
      test_buildFD2ndDeriv_order2))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_buildFD2ndDeriv_order4', &
      test_buildFD2ndDeriv_order4))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_central_coeffs_sum_to_zero', &
      test_central_coeffs_sum_to_zero))
   call suite%addTest(t)


end function test_finitedifferences_suite

