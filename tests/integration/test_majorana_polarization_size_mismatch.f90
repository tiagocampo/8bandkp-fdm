! Integration driver for Issue 04 Fix Round 1 (Important 1).
! Verifies that majorana_polarization now error-stops on size mismatch
! (replacing the previous silent zero-fill, which violated CLAUDE.md
! "No silent corrections.").
!
! This driver intentionally calls majorana_polarization with an
! undersized evec_bdg. The function must terminate the program with
! an error message containing "eigenvector size mismatch".
!
! pFUnit 4.x cannot catch Fortran `error stop` (per ADR 0002 and the
! precedent in test_parameters.pf), so this is run via shell as an
! integration test (see test_majorana_polarization_size_mismatch.sh).
program test_majorana_polarization_size_mismatch
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use topological_analysis, only: majorana_polarization, polarization_result_t
  implicit none

  integer, parameter :: n_sites = 8, half_n = 8 * n_sites
  ! Intentionally half-sized: only half_n components instead of 2*half_n.
  ! This triggers the size-mismatch branch in majorana_polarization.
  complex(kind=dp) :: evec_bdg(half_n)
  type(polarization_result_t) :: pol
  integer :: i

  do i = 1, half_n
    evec_bdg(i) = cmplx(1.0_dp / sqrt(real(half_n, kind=dp)), 0.0_dp, kind=dp)
  end do

  ! Expected: error stop with message containing "eigenvector size mismatch".
  ! If we reach the next line, the error stop did NOT fire — exit 1.
  pol = majorana_polarization(evec_bdg, n_sites)

  ! Should never reach here.
  write(*, '(a)') 'FAIL: expected error stop was not raised'
  write(*, '(a,es12.4)') 'half_wire_integral=', pol%half_wire_integral
  error stop 'majorana_polarization did not error stop on size mismatch'
end program test_majorana_polarization_size_mismatch
