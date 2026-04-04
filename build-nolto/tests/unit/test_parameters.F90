module test_parameters
  use funit
  use definitions
  use parameters
  implicit none

contains

  !@test
  subroutine test_gaas_parameters()
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)

    material(1) = "GaAs"
    call paramDatabase(material, 1, params)

    ! Check key GaAs parameters
#line 18 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.067_dp, params(1)%meff, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 18) )
  if (anyExceptions()) return
#line 19 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 19 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(1.519_dp, params(1)%Eg, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 19) )
  if (anyExceptions()) return
#line 20 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 20 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.341_dp, params(1)%deltaSO, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 20) )
  if (anyExceptions()) return
#line 21 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 21 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(6.98_dp, params(1)%gamma1, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 21) )
  if (anyExceptions()) return
#line 22 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 22 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(2.06_dp, params(1)%gamma2, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 22) )
  if (anyExceptions()) return
#line 23 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 23 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(2.93_dp, params(1)%gamma3, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 23) )
  if (anyExceptions()) return
#line 24 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 24 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(28.8_dp, params(1)%EP, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 24) )
  if (anyExceptions()) return
#line 25 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_gaas_parameters

  !@test
  subroutine test_inas_parameters()
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)

    material(1) = "InAs"
    call paramDatabase(material, 1, params)

#line 35 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.417_dp, params(1)%Eg, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 35) )
  if (anyExceptions()) return
#line 36 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 36 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.39_dp, params(1)%deltaSO, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 36) )
  if (anyExceptions()) return
#line 37 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 37 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(20.0_dp, params(1)%gamma1, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 37) )
  if (anyExceptions()) return
#line 38 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
#line 38 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(8.5_dp, params(1)%gamma2, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 38) )
  if (anyExceptions()) return
#line 39 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_inas_parameters

  !@test
  subroutine test_ec_ev_consistency()
    ! EC should equal EV + Eg for all materials
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)

    material(1) = "GaAs"
    call paramDatabase(material, 1, params)
#line 49 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(params(1)%EV + params(1)%Eg, params(1)%EC, tolerance=1.0e-8_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 49) )
  if (anyExceptions()) return
#line 50 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"

    material(1) = "InAs"
    call paramDatabase(material, 1, params)
#line 53 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(params(1)%EV + params(1)%Eg, params(1)%EC, tolerance=1.0e-6_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 53) )
  if (anyExceptions()) return
#line 54 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"

    material(1) = "AlSb"
    call paramDatabase(material, 1, params)
#line 57 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(params(1)%EV + params(1)%Eg, params(1)%EC, tolerance=1.0e-6_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 57) )
  if (anyExceptions()) return
#line 58 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_ec_ev_consistency

  !@test
  subroutine test_p_from_ep()
    ! P = sqrt(EP * const) for all materials
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)

    material(1) = "GaAs"
    call paramDatabase(material, 1, params)
#line 68 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(sqrt(params(1)%EP * const), params(1)%P, tolerance=1.0e-10_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 68) )
  if (anyExceptions()) return
#line 69 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_p_from_ep

  !@test
  subroutine test_multi_material()
    type(paramStruct) :: params(3)
    character(len=255) :: material(3)

    material(1) = "AlSb"
    material(2) = "GaSb"
    material(3) = "InAs"
    call paramDatabase(material, 3, params)

    ! AlSb Eg = 2.386
#line 82 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(2.386_dp, params(1)%Eg, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 82) )
  if (anyExceptions()) return
#line 83 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
    ! GaSb Eg = 0.812
#line 84 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.812_dp, params(2)%Eg, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
    ! InAs Eg = 0.417
#line 86 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(0.417_dp, params(3)%Eg, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 86) )
  if (anyExceptions()) return
#line 87 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_multi_material

  !@test
  subroutine test_dielectric_constants()
    ! Verify eps0 values against Vurgaftman 2001 Table III
    type(paramStruct) :: params(4)
    character(len=255) :: material(4)

    material(1) = "GaAs"
    material(2) = "InAs"
    material(3) = "AlSb"
    material(4) = "InSb"
    call paramDatabase(material, 4, params)

    ! GaAs eps0 = 12.90 (Vurgaftman 2001)
#line 102 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(12.90_dp, params(1)%eps0, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 102) )
  if (anyExceptions()) return
#line 103 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
    ! InAs eps0 = 15.15 (Vurgaftman 2001)
#line 104 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(15.15_dp, params(2)%eps0, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 104) )
  if (anyExceptions()) return
#line 105 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
    ! AlSb eps0 = 12.04 (Vurgaftman 2001)
#line 106 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(12.04_dp, params(3)%eps0, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 106) )
  if (anyExceptions()) return
#line 107 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
    ! InSb eps0 = 16.80 (Vurgaftman 2001)
#line 108 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  call assertEqual(16.80_dp, params(4)%eps0, tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_parameters.pf', &
 & 108) )
  if (anyExceptions()) return
#line 109 "/data/8bandkp-fdm/tests/unit/test_parameters.pf"
  end subroutine test_dielectric_constants

end module test_parameters

module Wraptest_parameters
   use FUnit
   use test_parameters
   implicit none
   private

contains


end module Wraptest_parameters

function test_parameters_suite() result(suite)
   use FUnit
   use test_parameters
   use Wraptest_parameters
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   suite = TestSuite('test_parameters_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_gaas_parameters', &
      test_gaas_parameters))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_inas_parameters', &
      test_inas_parameters))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_ec_ev_consistency', &
      test_ec_ev_consistency))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_p_from_ep', &
      test_p_from_ep))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_multi_material', &
      test_multi_material))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_dielectric_constants', &
      test_dielectric_constants))
   call suite%addTest(t)


end function test_parameters_suite

