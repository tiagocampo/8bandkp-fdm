module test_sc_loop
  use funit
  use definitions
  use sc_loop
  implicit none

contains

  !@test
  subroutine test_linear_mix_uniform()
    ! phi_old = 1.0, phi_poisson = 2.0, alpha = 0.3
    ! phi_new = 0.7*1.0 + 0.3*2.0 = 1.3

    integer, parameter :: N = 5
    real(kind=dp) :: phi_new(N), phi_old(N), phi_poisson(N)
    real(kind=dp) :: alpha
    integer :: i

    phi_old = 1.0_dp
    phi_poisson = 2.0_dp
    alpha = 0.3_dp

    call linear_mix(phi_new, phi_old, phi_poisson, N, alpha)

    do i = 1, N
#line 26 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(1.3_dp, phi_new(i), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 26) )
  if (anyExceptions()) return
#line 27 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
    end do
  end subroutine test_linear_mix_uniform


  !@test
  subroutine test_linear_mix_alpha_zero()
    ! alpha=0 => phi_new = phi_old (no update)

    integer, parameter :: N = 3
    real(kind=dp) :: phi_new(N), phi_old(N), phi_poisson(N)

    phi_old = 5.0_dp
    phi_poisson = 100.0_dp

    call linear_mix(phi_new, phi_old, phi_poisson, N, 0.0_dp)

#line 43 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(5.0_dp, phi_new(1), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 43) )
  if (anyExceptions()) return
#line 44 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
#line 44 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(5.0_dp, phi_new(2), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 44) )
  if (anyExceptions()) return
#line 45 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
#line 45 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(5.0_dp, phi_new(3), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 45) )
  if (anyExceptions()) return
#line 46 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  end subroutine test_linear_mix_alpha_zero


  !@test
  subroutine test_linear_mix_alpha_one()
    ! alpha=1 => phi_new = phi_poisson (full replacement)

    integer, parameter :: N = 3
    real(kind=dp) :: phi_new(N), phi_old(N), phi_poisson(N)

    phi_old = 5.0_dp
    phi_poisson = 100.0_dp

    call linear_mix(phi_new, phi_old, phi_poisson, N, 1.0_dp)

#line 61 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(100.0_dp, phi_new(1), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 61) )
  if (anyExceptions()) return
#line 62 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
#line 62 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(100.0_dp, phi_new(2), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 62) )
  if (anyExceptions()) return
#line 63 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
#line 63 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(100.0_dp, phi_new(3), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 63) )
  if (anyExceptions()) return
#line 64 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  end subroutine test_linear_mix_alpha_one


  !@test
  subroutine test_diis_fallback_to_linear()
    ! DIIS with m<2 should fall back to linear mixing.
    ! iter=1, diis_len=7 => m = min(0, 7) = 0 < 2 => linear

    integer, parameter :: N = 4
    real(kind=dp) :: phi_new(N), phi_old(N), phi_poisson(N)
    real(kind=dp) :: phi_history(N, 7), res_history(N, 7)
    real(kind=dp) :: alpha
    integer :: i

    phi_old = 1.0_dp
    phi_poisson = 3.0_dp
    phi_history = 0.0_dp
    res_history = 0.0_dp
    alpha = 0.5_dp

    call diis_extrapolate(phi_new, phi_old, phi_poisson, N, &
      & phi_history, res_history, 7, alpha, 1)

    ! Linear: 0.5*1 + 0.5*3 = 2.0
    do i = 1, N
#line 89 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(2.0_dp, phi_new(i), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
    end do
  end subroutine test_diis_fallback_to_linear


  !@test
  subroutine test_find_fermi_level_simple()
    ! Two subbands, 3 k_par points, no doping.
    ! CB at E=0.5, VB at E=-0.5.
    ! With no doping, Fermi level should be between VB and CB.

    integer, parameter :: N = 3
    integer, parameter :: num_subbands = 2
    integer, parameter :: nk = 3
    integer, parameter :: nz = 3
    real(kind=dp), parameter :: dz = 1.0_dp

    real(kind=dp) :: eig_kpar(num_subbands, nk)
    complex(kind=dp) :: eigv_kpar(8*N, num_subbands, nk)
    real(kind=dp) :: kpar_grid(nk)
    type(simulation_config) :: cfg
    real(kind=dp) :: rho_doping(nz)
    real(kind=dp) :: mu

    ! CB subband at E=0.5 eV, VB at E=-0.5 eV
    eig_kpar(1, :) = -0.5_dp   ! VB
    eig_kpar(2, :) = 0.5_dp    ! CB

    ! Eigenvectors: normalized over all z-points
    ! VB (subband 1): band 1 (HH), all z-points, 1/sqrt(N) each
    ! CB (subband 2): band 7, all z-points, 1/sqrt(N) each
    eigv_kpar = cmplx(0.0_dp, 0.0_dp, kind=dp)
    ! VB: band 1
    eigv_kpar(1, 1, :) = cmplx(1.0_dp / sqrt(real(N, kind=dp)), 0.0_dp, kind=dp)
    eigv_kpar(2, 1, :) = cmplx(1.0_dp / sqrt(real(N, kind=dp)), 0.0_dp, kind=dp)
    eigv_kpar(3, 1, :) = cmplx(1.0_dp / sqrt(real(N, kind=dp)), 0.0_dp, kind=dp)
    ! CB: band 7
    eigv_kpar(6*N+1, 2, :) = cmplx(1.0_dp / sqrt(real(N, kind=dp)), 0.0_dp, kind=dp)
    eigv_kpar(6*N+2, 2, :) = cmplx(1.0_dp / sqrt(real(N, kind=dp)), 0.0_dp, kind=dp)
    eigv_kpar(6*N+3, 2, :) = cmplx(1.0_dp / sqrt(real(N, kind=dp)), 0.0_dp, kind=dp)

    kpar_grid(1) = 0.0_dp
    kpar_grid(2) = 0.1_dp
    kpar_grid(3) = 0.2_dp

    ! Minimal config
    cfg%fdStep = N
    cfg%dz = dz
    cfg%numcb = 1
    cfg%sc%temperature = 300.0_dp
    cfg%sc%fermi_mode = 0

    rho_doping = 0.0_dp

    mu = find_fermi_level(eig_kpar, eigv_kpar, kpar_grid, &
      & cfg, 8*N, num_subbands, nk, nz, dz, rho_doping)

    ! Fermi level should be between VB (-0.5) and CB (0.5)
#line 147 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertTrue(mu > -0.5_dp, message="Fermi above VB top", &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 147) )
  if (anyExceptions()) return
#line 148 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
#line 148 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertTrue(mu < 0.5_dp, message="Fermi below CB bottom", &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 148) )
  if (anyExceptions()) return
#line 149 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  end subroutine test_find_fermi_level_simple


  !@test
  subroutine test_build_epsilon()
    ! Build epsilon for a 2-layer structure with known dielectrics.

    integer, parameter :: nz = 10
    real(kind=dp) :: epsilon(nz)
    type(simulation_config) :: cfg
    integer :: iz

    ! 2 layers: layer 1 z=1..5, layer 2 z=6..10
    cfg%numLayers = 2
    allocate(cfg%intStartPos(2), cfg%intEndPos(2))
    allocate(cfg%params(2))

    cfg%intStartPos(1) = 1;  cfg%intEndPos(1) = 5
    cfg%intStartPos(2) = 6;  cfg%intEndPos(2) = 10

    cfg%params(1)%eps0 = 12.90_dp  ! GaAs
    cfg%params(2)%eps0 = 10.06_dp  ! AlAs

    call build_epsilon(epsilon, cfg, nz)

    ! Layer 1 should be GaAs epsilon
    do iz = 1, 5
#line 176 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(12.90_dp, epsilon(iz), tolerance=1.0e-10_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 176) )
  if (anyExceptions()) return
#line 177 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
    end do

    ! Layer 2 should be AlAs epsilon
    do iz = 6, 10
#line 181 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(10.06_dp, epsilon(iz), tolerance=1.0e-10_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 181) )
  if (anyExceptions()) return
#line 182 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
    end do

    deallocate(cfg%intStartPos, cfg%intEndPos, cfg%params)
  end subroutine test_build_epsilon


  !@test
  subroutine test_build_doping_charge()
    ! Layer 1: ND=1e18, NA=0 => net doping = +1e18
    ! Layer 2: ND=0, NA=5e17 => net doping = -5e17

    integer, parameter :: nz = 10
    real(kind=dp) :: rho_doping(nz)
    type(simulation_config) :: cfg
    integer :: iz

    cfg%numLayers = 2
    allocate(cfg%intStartPos(2), cfg%intEndPos(2))
    allocate(cfg%doping(2))

    cfg%intStartPos(1) = 1;  cfg%intEndPos(1) = 5
    cfg%intStartPos(2) = 6;  cfg%intEndPos(2) = 10

    cfg%doping(1)%ND = 1.0e18_dp
    cfg%doping(1)%NA = 0.0_dp
    cfg%doping(2)%ND = 0.0_dp
    cfg%doping(2)%NA = 5.0e17_dp

    call build_doping_charge(rho_doping, cfg, nz)

    do iz = 1, 5
#line 213 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(1.0e18_dp, rho_doping(iz), tolerance=1.0e-10_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 213) )
  if (anyExceptions()) return
#line 214 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
    end do
    do iz = 6, 10
#line 216 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(-5.0e17_dp, rho_doping(iz), tolerance=1.0e-10_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 216) )
  if (anyExceptions()) return
#line 217 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
    end do

    deallocate(cfg%intStartPos, cfg%intEndPos, cfg%doping)
  end subroutine test_build_doping_charge


  !@test
  subroutine test_find_fermi_level_with_doping()
    ! Same as test_find_fermi_level_simple but with n-type doping.
    ! Positive rho_doping means net donors (ND - NA > 0).
    ! With donors present, Fermi level should shift toward CB (higher than
    ! the intrinsic case where it sits between VB and CB).

    integer, parameter :: N = 3
    integer, parameter :: num_subbands = 2
    integer, parameter :: nk = 3
    integer, parameter :: nz = 3
    real(kind=dp), parameter :: dz = 1.0_dp

    real(kind=dp) :: eig_kpar(num_subbands, nk)
    complex(kind=dp) :: eigv_kpar(8*N, num_subbands, nk)
    real(kind=dp) :: kpar_grid(nk)
    type(simulation_config) :: cfg
    real(kind=dp) :: rho_doping(nz)
    real(kind=dp) :: mu, mu_intrinsic

    ! CB subband at E=0.5 eV, VB at E=-0.5 eV
    eig_kpar(1, :) = -0.5_dp   ! VB
    eig_kpar(2, :) = 0.5_dp    ! CB

    ! Eigenvectors: normalized over all z-points
    eigv_kpar = cmplx(0.0_dp, 0.0_dp, kind=dp)
    ! VB: band 1
    eigv_kpar(1, 1, :) = cmplx(1.0_dp / sqrt(real(N, kind=dp)), 0.0_dp, kind=dp)
    eigv_kpar(2, 1, :) = cmplx(1.0_dp / sqrt(real(N, kind=dp)), 0.0_dp, kind=dp)
    eigv_kpar(3, 1, :) = cmplx(1.0_dp / sqrt(real(N, kind=dp)), 0.0_dp, kind=dp)
    ! CB: band 7
    eigv_kpar(6*N+1, 2, :) = cmplx(1.0_dp / sqrt(real(N, kind=dp)), 0.0_dp, kind=dp)
    eigv_kpar(6*N+2, 2, :) = cmplx(1.0_dp / sqrt(real(N, kind=dp)), 0.0_dp, kind=dp)
    eigv_kpar(6*N+3, 2, :) = cmplx(1.0_dp / sqrt(real(N, kind=dp)), 0.0_dp, kind=dp)

    kpar_grid(1) = 0.0_dp
    kpar_grid(2) = 0.1_dp
    kpar_grid(3) = 0.2_dp

    cfg%fdStep = N
    cfg%dz = dz
    cfg%numcb = 1
    cfg%sc%temperature = 300.0_dp
    cfg%sc%fermi_mode = 0

    ! First: intrinsic (no doping) to get baseline
    rho_doping = 0.0_dp

    mu_intrinsic = find_fermi_level(eig_kpar, eigv_kpar, kpar_grid, &
      & cfg, 8*N, num_subbands, nk, nz, dz, rho_doping)

    ! Now: with n-type doping (donors push Fermi toward CB)
    rho_doping = 1.0e18_dp  ! ND - NA > 0 everywhere

    mu = find_fermi_level(eig_kpar, eigv_kpar, kpar_grid, &
      & cfg, 8*N, num_subbands, nk, nz, dz, rho_doping)

    ! With donors, Fermi level should be higher than intrinsic
#line 281 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertTrue(mu > mu_intrinsic, message="Doped Fermi above intrinsic Fermi", &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 281) )
  if (anyExceptions()) return
#line 282 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
    ! Fermi should still be in gap
#line 283 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertTrue(mu > -0.5_dp, message="Fermi above VB top", &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 283) )
  if (anyExceptions()) return
#line 284 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
#line 284 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertTrue(mu < 0.5_dp, message="Fermi below CB bottom", &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 284) )
  if (anyExceptions()) return
#line 285 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  end subroutine test_find_fermi_level_with_doping

  !@test
  subroutine test_apply_potential_to_profile()
    ! Verify profile = profile_base - phi for all 3 columns
    integer, parameter :: nz = 5
    real(kind=dp) :: profile(nz, 3), profile_base(nz, 3), phi(nz)
    integer :: iz

    do iz = 1, nz
      profile_base(iz, :) = real(iz, kind=dp) * 10.0_dp
      phi(iz) = real(iz, kind=dp)
    end do

    call apply_potential_to_profile(profile, profile_base, phi, nz)

    do iz = 1, nz
#line 302 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(profile_base(iz, 1) - phi(iz), profile(iz, 1), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 302) )
  if (anyExceptions()) return
#line 303 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
#line 303 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(profile_base(iz, 2) - phi(iz), profile(iz, 2), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 303) )
  if (anyExceptions()) return
#line 304 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
#line 304 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(profile_base(iz, 3) - phi(iz), profile(iz, 3), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 304) )
  if (anyExceptions()) return
#line 305 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
    end do
  end subroutine test_apply_potential_to_profile


  !@test
  subroutine test_map_layer_to_grid()
    ! 3-layer structure: layer 1 = z(1..4), layer 2 = z(5..7), layer 3 = z(8..10)
    ! map_layer_to_grid should return the correct layer index for each grid point.

    integer, parameter :: nz = 10
    integer, allocatable :: layer_index(:)
    type(simulation_config) :: cfg
    integer :: iz

    cfg%numLayers = 3
    allocate(cfg%intStartPos(3), cfg%intEndPos(3))
    cfg%intStartPos(1) = 1;  cfg%intEndPos(1) = 4
    cfg%intStartPos(2) = 5;  cfg%intEndPos(2) = 7
    cfg%intStartPos(3) = 8;  cfg%intEndPos(3) = 10

    call map_layer_to_grid(layer_index, cfg, nz)

    do iz = 1, 4
#line 328 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(1, layer_index(iz), message="Layer 1 mapping", &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 328) )
  if (anyExceptions()) return
#line 329 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
    end do
    do iz = 5, 7
#line 331 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(2, layer_index(iz), message="Layer 2 mapping", &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 331) )
  if (anyExceptions()) return
#line 332 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
    end do
    do iz = 8, 10
#line 334 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
  call assertEqual(3, layer_index(iz), message="Layer 3 mapping", &
 & location=SourceLocation( &
 & 'test_sc_loop.pf', &
 & 334) )
  if (anyExceptions()) return
#line 335 "/data/8bandkp-fdm/tests/unit/test_sc_loop.pf"
    end do

    deallocate(cfg%intStartPos, cfg%intEndPos)
  end subroutine test_map_layer_to_grid

end module test_sc_loop

module Wraptest_sc_loop
   use FUnit
   use test_sc_loop
   implicit none
   private

contains


end module Wraptest_sc_loop

function test_sc_loop_suite() result(suite)
   use FUnit
   use test_sc_loop
   use Wraptest_sc_loop
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   suite = TestSuite('test_sc_loop_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_linear_mix_uniform', &
      test_linear_mix_uniform))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_linear_mix_alpha_zero', &
      test_linear_mix_alpha_zero))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_linear_mix_alpha_one', &
      test_linear_mix_alpha_one))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_diis_fallback_to_linear', &
      test_diis_fallback_to_linear))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_find_fermi_level_simple', &
      test_find_fermi_level_simple))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_build_epsilon', &
      test_build_epsilon))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_build_doping_charge', &
      test_build_doping_charge))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_find_fermi_level_with_doping', &
      test_find_fermi_level_with_doping))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_apply_potential_to_profile', &
      test_apply_potential_to_profile))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_map_layer_to_grid', &
      test_map_layer_to_grid))
   call suite%addTest(t)


end function test_sc_loop_suite

