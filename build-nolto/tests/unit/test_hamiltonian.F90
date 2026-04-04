module test_hamiltonian
  use funit
  use definitions
  use parameters
  use hamiltonianConstructor
  implicit none

contains

  !@test
  subroutine test_bulk_hermitian()
    ! ZB8bandBulk should produce a Hermitian 8x8 matrix
    complex(kind=dp) :: HT(8,8)
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)
    type(wavevector) :: wv
    integer :: i, j
    real(kind=dp) :: max_err

    material(1) = "GaAs"
    call paramDatabase(material, 1, params)

    wv%kx = 0.05_dp
    wv%ky = 0.03_dp
    wv%kz = 0.02_dp

    HT = cmplx(0.0_dp, 0.0_dp, kind=dp)
    call ZB8bandBulk(HT, wv, params)

    ! Check H = H^dagger
    max_err = 0.0_dp
    do j = 1, 8
      do i = 1, 8
        max_err = max(max_err, abs(HT(i,j) - conjg(HT(j,i))))
      end do
    end do
#line 37 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  call assertTrue(max_err < 1.0e-12_dp, message="Bulk Hamiltonian is Hermitian", &
 & location=SourceLocation( &
 & 'test_hamiltonian.pf', &
 & 37) )
  if (anyExceptions()) return
#line 38 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  end subroutine test_bulk_hermitian

  !@test
  subroutine test_bulk_gamma_point()
    ! At k=0, the bulk Hamiltonian should be diagonal with known energies
    complex(kind=dp) :: HT(8,8)
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)
    type(wavevector) :: wv
    real(kind=dp) :: expected_CB, expected_HH, expected_SO

    material(1) = "GaAs"
    call paramDatabase(material, 1, params)

    wv%kx = 0.0_dp
    wv%ky = 0.0_dp
    wv%kz = 0.0_dp

    HT = cmplx(0.0_dp, 0.0_dp, kind=dp)
    call ZB8bandBulk(HT, wv, params)

    ! At k=0: Q=T=0, so valence diagonal is 0, SO diagonal is -deltaSO, CB is +Eg
    ! Bands 1-4 (HH,LH,LH,HH): diagonal = 0
#line 61 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  call assertEqual(0.0_dp, real(HT(1,1), kind=dp), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian.pf', &
 & 61) )
  if (anyExceptions()) return
#line 62 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
#line 62 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  call assertEqual(0.0_dp, real(HT(2,2), kind=dp), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian.pf', &
 & 62) )
  if (anyExceptions()) return
#line 63 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
    ! Bands 5-6 (SO): diagonal = -deltaSO
#line 64 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  call assertEqual(-params(1)%deltaSO, real(HT(5,5), kind=dp), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian.pf', &
 & 64) )
  if (anyExceptions()) return
#line 65 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
#line 65 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  call assertEqual(-params(1)%deltaSO, real(HT(6,6), kind=dp), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian.pf', &
 & 65) )
  if (anyExceptions()) return
#line 66 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
    ! Bands 7-8 (CB): diagonal = Eg
#line 67 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  call assertEqual(params(1)%Eg, real(HT(7,7), kind=dp), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian.pf', &
 & 67) )
  if (anyExceptions()) return
#line 68 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
#line 68 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  call assertEqual(params(1)%Eg, real(HT(8,8), kind=dp), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian.pf', &
 & 68) )
  if (anyExceptions()) return
#line 69 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"

    ! Off-diagonals should all be zero at k=0
#line 71 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  call assertEqual(0.0_dp, abs(HT(1,7)), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian.pf', &
 & 71) )
  if (anyExceptions()) return
#line 72 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
#line 72 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  call assertEqual(0.0_dp, abs(HT(1,8)), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian.pf', &
 & 72) )
  if (anyExceptions()) return
#line 73 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  end subroutine test_bulk_gamma_point

  !@test
  subroutine test_bulk_kx_only()
    ! With only kx != 0, some off-diagonal terms should be non-zero
    complex(kind=dp) :: HT(8,8)
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)
    type(wavevector) :: wv

    material(1) = "GaAs"
    call paramDatabase(material, 1, params)

    wv%kx = 0.1_dp
    wv%ky = 0.0_dp
    wv%kz = 0.0_dp

    HT = cmplx(0.0_dp, 0.0_dp, kind=dp)
    call ZB8bandBulk(HT, wv, params)

    ! PP = P * kx / sqrt(2) should give non-zero (1,7) and (7,1) elements
#line 94 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  call assertTrue(abs(HT(1,7)) > 1.0e-6_dp, message="kx produces off-diagonal CB-VB coupling", &
 & location=SourceLocation( &
 & 'test_hamiltonian.pf', &
 & 94) )
  if (anyExceptions()) return
#line 95 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
    ! Check Hermiticity of off-diagonal
#line 96 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  call assertEqual(conjg(HT(1,7)), HT(7,1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian.pf', &
 & 96) )
  if (anyExceptions()) return
#line 97 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  end subroutine test_bulk_kx_only

  !@test
  subroutine test_bulk_hermitian_inas()
    ! Also check InAs bulk Hamiltonian is Hermitian
    complex(kind=dp) :: HT(8,8)
    type(paramStruct) :: params(1)
    character(len=255) :: material(1)
    type(wavevector) :: wv
    integer :: i, j
    real(kind=dp) :: max_err

    material(1) = "InAs"
    call paramDatabase(material, 1, params)

    wv%kx = 0.08_dp
    wv%ky = 0.04_dp
    wv%kz = 0.06_dp

    HT = cmplx(0.0_dp, 0.0_dp, kind=dp)
    call ZB8bandBulk(HT, wv, params)

    max_err = 0.0_dp
    do j = 1, 8
      do i = 1, 8
        max_err = max(max_err, abs(HT(i,j) - conjg(HT(j,i))))
      end do
    end do
#line 125 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  call assertTrue(max_err < 1.0e-12_dp, message="InAs bulk Hamiltonian is Hermitian", &
 & location=SourceLocation( &
 & 'test_hamiltonian.pf', &
 & 125) )
  if (anyExceptions()) return
#line 126 "/data/8bandkp-fdm/tests/unit/test_hamiltonian.pf"
  end subroutine test_bulk_hermitian_inas

end module test_hamiltonian

module Wraptest_hamiltonian
   use FUnit
   use test_hamiltonian
   implicit none
   private

contains


end module Wraptest_hamiltonian

function test_hamiltonian_suite() result(suite)
   use FUnit
   use test_hamiltonian
   use Wraptest_hamiltonian
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   suite = TestSuite('test_hamiltonian_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_bulk_hermitian', &
      test_bulk_hermitian))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_bulk_gamma_point', &
      test_bulk_gamma_point))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_bulk_kx_only', &
      test_bulk_kx_only))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_bulk_hermitian_inas', &
      test_bulk_hermitian_inas))
   call suite%addTest(t)


end function test_hamiltonian_suite

