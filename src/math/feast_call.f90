! Non-module wrappers for MKL FEAST calls.
! Must be outside any module to avoid gfortran allocatable descriptor issues
! when passing arrays to external routines from inside module procedures.
!
! All array arguments use explicit-shape declarations (not assumed-size *)
! to ensure gfortran passes bare data pointers with no copy-in/copy-out,
! matching the MKL C ABI that zfeast_hcsrev expects.

#ifdef USE_MKL_FEAST

subroutine feastinit_call(fpm)
  implicit none
  integer, intent(inout) :: fpm(128)
  external :: feastinit
  call feastinit(fpm)
end subroutine feastinit_call

subroutine feast_solve_hermitian_csr(uplo, n, nnz, val, rowptr, colind, &
                                      fpm, epsout, loop, emin, emax, m0, &
                                      e, x, m, res, info)
  use definitions, only: dp
  implicit none
  character(len=1), intent(in)    :: uplo
  integer, intent(in)             :: n, nnz, m0
  complex(kind=dp), intent(in)    :: val(nnz)
  integer, intent(in)             :: rowptr(n+1), colind(nnz)
  integer, intent(inout)          :: fpm(128)
  real(kind=dp), intent(out)      :: epsout
  integer, intent(out)            :: loop, m, info
  real(kind=dp), intent(in)       :: emin, emax
  real(kind=dp), intent(inout)    :: e(m0), res(m0)
  complex(kind=dp), intent(inout) :: x(n, m0)

  external :: zfeast_hcsrev

  call zfeast_hcsrev(uplo, n, val, rowptr, colind, &
                     fpm, epsout, loop, emin, emax, m0, &
                     e, x, m, res, info)
end subroutine feast_solve_hermitian_csr

#endif /* USE_MKL_FEAST */
