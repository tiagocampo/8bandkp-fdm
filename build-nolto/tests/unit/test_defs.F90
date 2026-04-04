module test_defs
  use funit
  use definitions
  implicit none

contains

  !@test
  subroutine test_kronij_same()
#line 10 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1, kronij(3, 3), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 10) )
  if (anyExceptions()) return
#line 11 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 11 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1, kronij(1, 1), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 11) )
  if (anyExceptions()) return
#line 12 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 12 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1, kronij(0, 0), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 12) )
  if (anyExceptions()) return
#line 13 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_kronij_same

  !@test
  subroutine test_kronij_different()
#line 17 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(0, kronij(1, 2), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 17) )
  if (anyExceptions()) return
#line 18 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 18 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(0, kronij(3, 1), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 18) )
  if (anyExceptions()) return
#line 19 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 19 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(0, kronij(-1, 1), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 19) )
  if (anyExceptions()) return
#line 20 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_kronij_different

  !@test
  subroutine test_pi_constant()
    real(kind=dp) :: expected_pi
    expected_pi = acos(-1.0_dp)
#line 26 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(expected_pi, pi_dp, tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 26) )
  if (anyExceptions()) return
#line 27 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_pi_constant

  !@test
  subroutine test_imaginary_unit()
    complex(kind=dp) :: result
    result = IU * IU
#line 33 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(-1.0_dp, real(result, kind=dp), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 33) )
  if (anyExceptions()) return
#line 34 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 34 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(0.0_dp, aimag(result), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 34) )
  if (anyExceptions()) return
#line 35 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_imaginary_unit

  !@test
  subroutine test_unit_complex()
#line 39 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(0.0_dp, real(ZERO, kind=dp), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 39) )
  if (anyExceptions()) return
#line 40 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 40 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(0.0_dp, aimag(ZERO), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 40) )
  if (anyExceptions()) return
#line 41 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 41 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1.0_dp, real(UM, kind=dp), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 41) )
  if (anyExceptions()) return
#line 42 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 42 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(0.0_dp, aimag(UM), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 42) )
  if (anyExceptions()) return
#line 43 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_unit_complex

  !@test
  subroutine test_sqrt_constants()
#line 47 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(sqrt(3.0_dp), SQR3, tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 47) )
  if (anyExceptions()) return
#line 48 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 48 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(sqrt(2.0_dp), SQR2, tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 48) )
  if (anyExceptions()) return
#line 49 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 49 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1.0_dp/SQR3, RQS3, tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 49) )
  if (anyExceptions()) return
#line 50 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 50 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1.0_dp/SQR2, RQS2, tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 50) )
  if (anyExceptions()) return
#line 51 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_sqrt_constants

end module test_defs

module Wraptest_defs
   use FUnit
   use test_defs
   implicit none
   private

contains


end module Wraptest_defs

function test_defs_suite() result(suite)
   use FUnit
   use test_defs
   use Wraptest_defs
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   suite = TestSuite('test_defs_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_kronij_same', &
      test_kronij_same))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_kronij_different', &
      test_kronij_different))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_pi_constant', &
      test_pi_constant))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_imaginary_unit', &
      test_imaginary_unit))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_unit_complex', &
      test_unit_complex))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_sqrt_constants', &
      test_sqrt_constants))
   call suite%addTest(t)


end function test_defs_suite

