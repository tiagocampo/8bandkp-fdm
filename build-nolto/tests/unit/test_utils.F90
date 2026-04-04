module test_utils
  use funit
  use definitions
  use utils
  implicit none

contains

  !@test
  subroutine test_simpson_linear()
    ! Integral of f(x) = x from 0 to 1 is 0.5
    complex(kind=dp), allocatable :: f(:)
    complex(kind=dp) :: result
    real(kind=dp) :: a, b, dx
    integer :: i, N

    N = 11  ! odd number of points
    a = 0.0_dp
    b = 1.0_dp
    dx = (b - a) / (N - 1)

    allocate(f(N))
    do i = 1, N
      f(i) = cmplx(a + (i-1)*dx, 0.0_dp, kind=dp)
    end do

    result = simpson(f, a, b)
#line 28 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  call assertEqual(0.5_dp, real(result, kind=dp), tolerance=1.0e-10_dp, &
 & location=SourceLocation( &
 & 'test_utils.pf', &
 & 28) )
  if (anyExceptions()) return
#line 29 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  end subroutine test_simpson_linear

  !@test
  subroutine test_simpson_quadratic()
    ! Integral of f(x) = x^2 from 0 to 1 is 1/3
    complex(kind=dp), allocatable :: f(:)
    complex(kind=dp) :: result
    real(kind=dp) :: a, b, dx
    integer :: i, N

    N = 11
    a = 0.0_dp
    b = 1.0_dp
    dx = (b - a) / (N - 1)

    allocate(f(N))
    do i = 1, N
      f(i) = cmplx((a + (i-1)*dx)**2, 0.0_dp, kind=dp)
    end do

    result = simpson(f, a, b)
#line 50 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  call assertEqual(1.0_dp/3.0_dp, real(result, kind=dp), tolerance=1.0e-10_dp, &
 & location=SourceLocation( &
 & 'test_utils.pf', &
 & 50) )
  if (anyExceptions()) return
#line 51 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  end subroutine test_simpson_quadratic

  !@test
  subroutine test_simpson_sine()
    ! Integral of sin(x) from 0 to pi is 2.0
    complex(kind=dp), allocatable :: f(:)
    complex(kind=dp) :: result
    real(kind=dp) :: a, b, dx
    integer :: i, N

    N = 101  ! use more points for trig
    a = 0.0_dp
    b = pi_dp
    dx = (b - a) / (N - 1)

    allocate(f(N))
    do i = 1, N
      f(i) = cmplx(sin(a + (i-1)*dx), 0.0_dp, kind=dp)
    end do

    result = simpson(f, a, b)
#line 72 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  call assertEqual(2.0_dp, real(result, kind=dp), tolerance=1.0e-7_dp, &
 & location=SourceLocation( &
 & 'test_utils.pf', &
 & 72) )
  if (anyExceptions()) return
#line 73 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  end subroutine test_simpson_sine

  !@test
  subroutine test_simpson_complex()
    ! Integral of e^(i*x) from 0 to 2*pi should be 0
    complex(kind=dp), allocatable :: f(:)
    complex(kind=dp) :: result
    real(kind=dp) :: a, b, dx
    integer :: i, N

    N = 101
    a = 0.0_dp
    b = 2.0_dp * pi_dp
    dx = (b - a) / (N - 1)

    allocate(f(N))
    do i = 1, N
      f(i) = exp(cmplx(0.0_dp, a + (i-1)*dx, kind=dp))
    end do

    result = simpson(f, a, b)
#line 94 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  call assertEqual(0.0_dp, real(result, kind=dp), tolerance=1.0e-6_dp, &
 & location=SourceLocation( &
 & 'test_utils.pf', &
 & 94) )
  if (anyExceptions()) return
#line 95 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
#line 95 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  call assertEqual(0.0_dp, aimag(result), tolerance=1.0e-6_dp, &
 & location=SourceLocation( &
 & 'test_utils.pf', &
 & 95) )
  if (anyExceptions()) return
#line 96 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  end subroutine test_simpson_complex

  !@test
  subroutine test_coo_insert_and_sort()
    ! Test COO insertion and sorting
    integer(kind=4), parameter :: nnz = 10
    complex(kind=dp), allocatable :: v(:)
    integer(kind=4), allocatable :: r(:), c(:)
    integer(kind=4) :: next, nnz_count

    allocate(v(nnz), r(nnz), c(nnz))

    ! Insert entries in scrambled order
    next = 1
    call insertCOO_cmplx(v, r, c, cmplx(3.0_dp, 0.0_dp, kind=dp), 3_4, 1_4, next, nnz)
    call insertCOO_cmplx(v, r, c, cmplx(1.0_dp, 0.0_dp, kind=dp), 1_4, 1_4, next, nnz)
    call insertCOO_cmplx(v, r, c, cmplx(2.0_dp, 0.0_dp, kind=dp), 2_4, 2_4, next, nnz)

    nnz_count = next - 1
#line 115 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  call assertEqual(3, nnz_count, &
 & location=SourceLocation( &
 & 'test_utils.pf', &
 & 115) )
  if (anyExceptions()) return
#line 116 "/data/8bandkp-fdm/tests/unit/test_utils.pf"

    call finalizeCOO_cmplx(v, r, c, nnz_count)
#line 118 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  call assertEqual(3, nnz_count, &
 & location=SourceLocation( &
 & 'test_utils.pf', &
 & 118) )
  if (anyExceptions()) return
#line 119 "/data/8bandkp-fdm/tests/unit/test_utils.pf"

    ! After sort, should be ordered by (row, col)
#line 121 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  call assertEqual(1_4, r(1), &
 & location=SourceLocation( &
 & 'test_utils.pf', &
 & 121) )
  if (anyExceptions()) return
#line 122 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
#line 122 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  call assertEqual(1_4, c(1), &
 & location=SourceLocation( &
 & 'test_utils.pf', &
 & 122) )
  if (anyExceptions()) return
#line 123 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
#line 123 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  call assertEqual(1.0_dp, real(v(1), kind=dp), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_utils.pf', &
 & 123) )
  if (anyExceptions()) return
#line 124 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  end subroutine test_coo_insert_and_sort

  !@test
  subroutine test_simpson_real_quadratic()
    ! Integral of f(x) = x^2 from 0 to 1 is exactly 1/3
    real(kind=dp), allocatable :: f(:)
    real(kind=dp) :: result, a, b, dx
    integer :: i, N

    N = 11
    a = 0.0_dp
    b = 1.0_dp
    dx = (b - a) / (N - 1)

    allocate(f(N))
    do i = 1, N
      f(i) = (a + (i-1)*dx)**2
    end do

    result = simpson_real(f, a, b)
#line 144 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  call assertEqual(1.0_dp/3.0_dp, result, tolerance=1.0e-10_dp, &
 & location=SourceLocation( &
 & 'test_utils.pf', &
 & 144) )
  if (anyExceptions()) return
#line 145 "/data/8bandkp-fdm/tests/unit/test_utils.pf"

    deallocate(f)
  end subroutine test_simpson_real_quadratic


  !@test
  subroutine test_coo_duplicate_merge()
    ! Insert same (row,col) twice — should merge
    integer(kind=4), parameter :: nnz = 10
    complex(kind=dp), allocatable :: v(:)
    integer(kind=4), allocatable :: r(:), c(:)
    integer(kind=4) :: next, nnz_count

    allocate(v(nnz), r(nnz), c(nnz))

    next = 1
    call insertCOO_cmplx(v, r, c, cmplx(2.0_dp, 0.0_dp, kind=dp), 1_4, 1_4, next, nnz)
    call insertCOO_cmplx(v, r, c, cmplx(3.0_dp, 0.0_dp, kind=dp), 1_4, 1_4, next, nnz)

    nnz_count = next - 1
    call finalizeCOO_cmplx(v, r, c, nnz_count)

    ! Should merge into single entry with value 5.0
#line 168 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  call assertEqual(1, nnz_count, &
 & location=SourceLocation( &
 & 'test_utils.pf', &
 & 168) )
  if (anyExceptions()) return
#line 169 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
#line 169 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  call assertEqual(5.0_dp, real(v(1), kind=dp), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_utils.pf', &
 & 169) )
  if (anyExceptions()) return
#line 170 "/data/8bandkp-fdm/tests/unit/test_utils.pf"
  end subroutine test_coo_duplicate_merge

end module test_utils

module Wraptest_utils
   use FUnit
   use test_utils
   implicit none
   private

contains


end module Wraptest_utils

function test_utils_suite() result(suite)
   use FUnit
   use test_utils
   use Wraptest_utils
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   suite = TestSuite('test_utils_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_simpson_linear', &
      test_simpson_linear))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_simpson_quadratic', &
      test_simpson_quadratic))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_simpson_sine', &
      test_simpson_sine))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_simpson_complex', &
      test_simpson_complex))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_coo_insert_and_sort', &
      test_coo_insert_and_sort))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_simpson_real_quadratic', &
      test_simpson_real_quadratic))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_coo_duplicate_merge', &
      test_coo_duplicate_merge))
   call suite%addTest(t)


end function test_utils_suite

