module green_functions

  ! ==============================================================================
  ! Green function and local density of states (LDOS) utilities.
  !
  ! LDOS(r, E) = -(1/pi) * Im[G(r, r, E + i*eta)]
  ! where G = (E + i*eta - H)^-1 is the retarded Green function.
  !
  ! Uses MKL PARDISO complex solver via pardiso_c (USE_ARPACK) to invert
  ! the shifted matrix A = E + i*eta - H for each energy point.
  ! ==============================================================================

  use definitions, only: dp, pi_dp
  use sparse_matrices
  use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t

#ifdef USE_ARPACK
  use linalg, only: pardiso_c
#endif

  implicit none
  private

#ifdef USE_ARPACK
  public :: compute_ldos_csr
contains

  ! ==============================================================================
  ! Compute LDOS at each grid point for a given energy E.
  !
  ! LDOS(r, E) = -(1/pi) * Im[G(r, r, E + i*eta)]
  ! where G = (E + i*eta - H)^-1
  !
  ! Args:
    !   H          -- CSR Hamiltonian matrix (complex)
    !   E          -- energy (real, eV)
    !   eta        -- broadening (real, eV)
    !   ldos       -- output LDOS at each grid point (real, 1/eV)
  !
  ! The shifted system is: A * x = b  with  A = E + i*eta - H
  ! For diagonal G(r,r) we solve A * e_r = e_r where e_r is unit vector at row r,
  ! then G(r,r) = x(r).  We extract diagonal by solving with unit vectors
  ! at each diagonal position.
  ! ==============================================================================
  subroutine compute_ldos_csr(H, E, eta, ldos)
    type(csr_matrix), intent(in) :: H
    real(kind=dp), intent(in) :: E
    real(kind=dp), intent(in) :: eta
    real(kind=dp), intent(out) :: ldos(:)

    integer :: N, nnz, i, r
    integer :: maxfct, mnum, mtype, phase, nrhs, msglvl, error_loc
    integer(kind=c_intptr_t) :: pt(64)
    integer :: iparm(64)
    complex(kind=dp), allocatable :: a_val(:), rhs(:), sol(:)
    integer, allocatable :: ia(:), ja(:), perm(:)
    complex(kind=dp) :: shift

    N = H%nrows

    if (N /= size(ldos)) then
      print *, 'ERROR: compute_ldos_csr: ldos size mismatch (N=', N, ', ldos=', size(ldos), ')'
      stop 1
    end if

    shift = cmplx(E, eta, kind=dp)

    ! Copy CSR structure
    nnz = H%nnz
    allocate(a_val(nnz), ia(N+1), ja(nnz))
    a_val = shift - H%values
    ia = H%rowptr
    ja = H%colind

    ! PARDISO setup
    pt = 0
    iparm = 0
    maxfct = 1
    mnum = 1
    mtype = 3  ! complex structurally symmetric
    nrhs = 1
    msglvl = 0
    error_loc = 0

    iparm(1) = 1   ! no solver default
    iparm(2) = 2   ! OpenMP nested dissection
    iparm(8) = 2   ! max iterative refinement steps
    iparm(10) = 13 ! perturb pivots with 1E-13
    iparm(11) = 1  ! scaling
    iparm(13) = 1  ! matching
    iparm(18) = -1 ! report number of nonzeros
    iparm(40) = 1  ! distributed matrix input (CSR)

    allocate(perm(N), rhs(N), sol(N))

    ! Phase 11: Symbolic factorization
    phase = 11
    call pardiso_c(pt, maxfct, mnum, mtype, phase, N, &
                   a_val, ia, ja, perm, nrhs, iparm, msglvl, &
                   rhs, sol, error_loc)
    if (error_loc /= 0) then
      print *, 'ERROR: compute_ldos_csr: symbolic factorization failed (error=', error_loc, ')'
      phase = -1
      call pardiso_c(pt, maxfct, mnum, mtype, phase, N, &
                     a_val, ia, ja, perm, nrhs, iparm, msglvl, &
                     rhs, sol, error_loc)
      stop 1
    end if

    ! Phase 22: Numeric factorization
    phase = 22
    call pardiso_c(pt, maxfct, mnum, mtype, phase, N, &
                   a_val, ia, ja, perm, nrhs, iparm, msglvl, &
                   rhs, sol, error_loc)
    if (error_loc /= 0) then
      print *, 'ERROR: compute_ldos_csr: numeric factorization failed (error=', error_loc, ')'
      phase = -1
      call pardiso_c(pt, maxfct, mnum, mtype, phase, N, &
                     a_val, ia, ja, perm, nrhs, iparm, msglvl, &
                     rhs, sol, error_loc)
      stop 1
    end if

    ! Extract diagonal of G by solving A * x = e_r for each r
    do r = 1, N
      rhs = cmplx(0.0_dp, 0.0_dp, kind=dp)
      rhs(r) = cmplx(1.0_dp, 0.0_dp, kind=dp)

      ! Phase 33: Solve
      phase = 33
      call pardiso_c(pt, maxfct, mnum, mtype, phase, N, &
                     a_val, ia, ja, perm, nrhs, iparm, msglvl, &
                     rhs, sol, error_loc)
      if (error_loc /= 0) then
        print *, 'ERROR: compute_ldos_csr: solve failed at r=', r, ' (error=', error_loc, ')'
        phase = -1
        call pardiso_c(pt, maxfct, mnum, mtype, phase, N, &
                       a_val, ia, ja, perm, nrhs, iparm, msglvl, &
                       rhs, sol, error_loc)
        stop 1
      end if

      ! LDOS = -(1/pi) * Im[G(r,r)] = -(1/pi) * Im[sol(r)]
      ldos(r) = -aimag(sol(r)) / pi_dp
    end do

    ! Phase -1: Release memory
    phase = -1
    call pardiso_c(pt, maxfct, mnum, mtype, phase, N, &
                   a_val, ia, ja, perm, nrhs, iparm, msglvl, &
                   rhs, sol, error_loc)

    deallocate(a_val, ia, ja, perm, rhs, sol)
  end subroutine compute_ldos_csr
#endif

end module green_functions