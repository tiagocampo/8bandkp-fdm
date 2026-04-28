module linalg
  ! Centralized explicit interface blocks for LAPACK, MKL utility, and
  ! (optionally) MKL FEAST routines.  Replaces loose `external` declarations
  ! scattered across the codebase, giving the compiler proper type-checking
  ! and eliminating gfortran strict-aliasing warnings.

  use definitions, only: dp
  implicit none
  private

  ! Standard LAPACK
  public :: zheevx
  public :: dgesv
  public :: ilaenv
  public :: dlamch

  ! MKL-specific
  public :: mkl_set_num_threads_local

  ! MKL FEAST (guarded)
#ifdef USE_MKL_FEAST
  public :: feastinit
  public :: zfeast_hcsrev
#endif

  ! ARPACK-NG (guarded)
#ifdef USE_ARPACK
  public :: znaupd
  public :: zneupd
#endif

  ! ================================================================
  ! Interface blocks
  ! ================================================================

  ! zheevx - LAPACK Hermitian eigensolver (selected eigenvalues/vectors)
  interface
    subroutine zheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
                      abstol, m, w, z, ldz, work, lwork, rwork, iwork, &
                      ifail, info)
      use definitions, only: dp
      character(len=1), intent(in) :: jobz, range, uplo
      integer, intent(in) :: n, lda, il, iu, ldz
      real(kind=dp), intent(in) :: vl, vu, abstol
      complex(kind=dp), intent(inout) :: a(lda, *)
      real(kind=dp), intent(out) :: w(*)
      complex(kind=dp), intent(out) :: z(ldz, *)
      complex(kind=dp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(kind=dp), intent(out) :: rwork(*)
      integer, intent(out) :: iwork(*), ifail(*)
      integer, intent(out) :: m, info
    end subroutine
  end interface

  ! dgesv - LAPACK general linear system (double precision)
  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      use definitions, only: dp
      integer, intent(in) :: n, nrhs, lda, ldb
      real(kind=dp), intent(inout) :: a(lda, *), b(ldb, *)
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
    end subroutine
  end interface

  ! ilaenv - LAPACK block size query
  interface
    function ilaenv(ispec, name, opts, n1, n2, n3, n4) result(value)
      character(len=*), intent(in) :: name, opts
      integer, intent(in) :: ispec, n1, n2, n3, n4
      integer :: value
    end function
  end interface

  ! dlamch - LAPACK machine constants
  interface
    function dlamch(cmach) result(value)
      use definitions, only: dp
      character(len=1), intent(in) :: cmach
      real(kind=dp) :: value
    end function
  end interface

  ! mkl_set_num_threads_local - MKL per-thread thread count
  interface
    function mkl_set_num_threads_local(nt) result(previous)
      integer, intent(in) :: nt
      integer :: previous
    end function
  end interface

#ifdef USE_MKL_FEAST
  ! feastinit - FEAST initialization
  interface
    subroutine feastinit(fpm)
      integer, intent(inout) :: fpm(128)
    end subroutine
  end interface

  ! zfeast_hcsrev - FEAST complex Hermitian CSR eigensolver
  interface
    subroutine zfeast_hcsrev(uplo, n, a, ia, ja, fpm, epsout, loop, &
                             emin, emax, m0, e, x, m, res, info)
      use definitions, only: dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, m0
      complex(kind=dp), intent(in) :: a(*)
      integer, intent(in) :: ia(*), ja(*)
      integer, intent(inout) :: fpm(128)
      real(kind=dp), intent(out) :: epsout
      integer, intent(out) :: loop, m, info
      real(kind=dp), intent(in) :: emin, emax
      real(kind=dp), intent(inout) :: e(m0), res(m0)
      complex(kind=dp), intent(inout) :: x(n, m0)
    end subroutine
  end interface
#endif

#ifdef USE_ARPACK
  ! znaupd - ARPACK-NG reverse communication for complex eigenproblems
  interface
    subroutine znaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                      iparam, ipntr, workd, workl, lworkl, rwork, info)
      use definitions, only: dp
      character(len=1), intent(in)    :: bmat
      character(len=2), intent(in)    :: which
      integer, intent(inout)          :: ido
      integer, intent(in)             :: n, nev, ncv, ldv, lworkl
      real(kind=dp), intent(in)       :: tol
      complex(kind=dp), intent(inout) :: resid(n)
      complex(kind=dp), intent(out)   :: v(ldv, ncv)
      integer, intent(inout)          :: iparam(11)
      integer, intent(out)            :: ipntr(14)
      complex(kind=dp), intent(inout) :: workd(3*n)
      complex(kind=dp), intent(inout) :: workl(lworkl)
      real(kind=dp), intent(out)      :: rwork(ncv)
      integer, intent(inout)          :: info
    end subroutine
  end interface

  ! zneupd - ARPACK-NG post-processing for complex eigenproblems
  interface
    subroutine zneupd(rvec, howmny, select, d, z, ldz, sigma, workev, &
                      bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                      iparam, ipntr, workd, workl, lworkl, rwork, info)
      use definitions, only: dp
      logical, intent(in)             :: rvec
      character(len=1), intent(in)    :: howmny, bmat
      character(len=2), intent(in)    :: which
      integer, intent(in)             :: n, nev, ncv, ldz, ldv, lworkl
      logical, intent(inout)          :: select(ncv)
      real(kind=dp), intent(in)       :: tol
      complex(kind=dp), intent(out)   :: d(nev)
      complex(kind=dp), intent(out)   :: z(ldz, nev)
      complex(kind=dp), intent(in)    :: sigma
      complex(kind=dp), intent(inout) :: workev(2*ncv)
      complex(kind=dp), intent(inout) :: resid(n)
      complex(kind=dp), intent(inout) :: v(ldv, ncv)
      integer, intent(inout)          :: iparam(11)
      integer, intent(inout)          :: ipntr(14)
      complex(kind=dp), intent(inout) :: workd(3*n)
      complex(kind=dp), intent(inout) :: workl(lworkl)
      real(kind=dp), intent(inout)    :: rwork(ncv)
      integer, intent(out)            :: info
    end subroutine
  end interface

#endif

contains
  ! (empty - this module provides only interface blocks)
end module
