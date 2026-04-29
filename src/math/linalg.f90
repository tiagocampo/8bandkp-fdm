module linalg
  ! Centralized explicit interface blocks for LAPACK, BLAS, MKL utility,
  ! MKL PARDISO, and (optionally) MKL FEAST routines.  Replaces loose
  ! `external` declarations scattered across the codebase, giving the
  ! compiler proper type-checking and eliminating gfortran strict-aliasing
  ! warnings.

  use definitions, only: dp
  use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_double, c_double_complex, c_char
  implicit none
  private

  ! Standard LAPACK
  public :: zheevx
  public :: dgesv
  public :: ilaenv
  public :: dlamch

  ! BLAS
  public :: zdotc

  ! MKL-specific
  public :: mkl_set_num_threads_local

  ! MKL PARDISO (guarded)
#ifdef USE_ARPACK
  public :: pardiso_c
#endif

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
    function mkl_set_num_threads_local(nt) result(previous) bind(C, name="MKL_Set_Num_Threads_Local")
      import :: c_int
      integer(c_int), value :: nt
      integer(c_int) :: previous
    end function
  end interface

  ! zdotc - BLAS complex dot product (conjugated first vector)
  interface
    function zdotc(n, zx, incx, zy, incy) result(res)
      import :: dp
      integer, intent(in) :: n, incx, incy
      complex(kind=dp), intent(in) :: zx(*), zy(*)
      complex(kind=dp) :: res
    end function zdotc
  end interface

#ifdef USE_ARPACK
  ! pardiso_c - MKL PARDISO direct solver (iso_c_binding)
  interface
    subroutine pardiso_c(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, &
                         nrhs, iparm, msglvl, b, x, error) bind(C, name="PARDISO")
      import :: c_int, c_int64_t, c_double_complex
      integer(c_int64_t), intent(inout) :: pt(64)
      integer(c_int), value :: maxfct, mnum, mtype, phase, n, nrhs, msglvl
      complex(c_double_complex), intent(in) :: a(*)
      integer(c_int), intent(in) :: ia(*), ja(*), perm(*)
      integer(c_int), intent(inout) :: iparm(64)
      complex(c_double_complex), intent(inout) :: b(*), x(*)
      integer(c_int), intent(out) :: error
    end subroutine pardiso_c
  end interface
#endif

#ifdef USE_MKL_FEAST
  ! feastinit - FEAST initialization (iso_c_binding)
  interface
    subroutine feastinit(fpm) bind(C, name="feastinit")
      import :: c_int
      integer(c_int), intent(inout) :: fpm(128)
    end subroutine
  end interface

  ! zfeast_hcsrev - FEAST complex Hermitian CSR eigensolver (iso_c_binding)
  interface
    subroutine zfeast_hcsrev(uplo, n, a, ia, ja, fpm, epsout, loop, &
                             emin, emax, m0, e, x, m, res, info) bind(C, name="zfeast_hcsrev")
      import :: c_int, c_double, c_double_complex, c_char
      character(c_char), intent(in) :: uplo
      integer(c_int), value :: n, m0
      complex(c_double_complex), intent(in) :: a(*)
      integer(c_int), intent(in) :: ia(*), ja(*)
      integer(c_int), intent(inout) :: fpm(128)
      real(c_double), intent(out) :: epsout
      integer(c_int), intent(out) :: loop, m, info
      real(c_double), intent(in) :: emin, emax
      real(c_double), intent(inout) :: e(m0), res(m0)
      complex(c_double_complex), intent(inout) :: x(n, m0)
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
