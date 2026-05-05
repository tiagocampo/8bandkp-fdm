module green_functions

  ! ==============================================================================
  ! Green function and spectral function utilities.
  !
  ! LDOS(r, E) = -(1/pi) * Im[G(r, r, E + i*eta)]  (requires PARDISO/ARPACK)
  ! Spectral function A(k, E) = sum_n Lorentzian(E - E_n, eta)
  !
  ! LDOS uses MKL PARDISO complex solver via pardiso_c (USE_ARPACK) to invert
  ! the shifted matrix A = E + i*eta - H for each energy point.
  ! ==============================================================================

  use definitions, only: dp, pi_dp, simulation_config, wavevector
  use sparse_matrices
  use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t

#ifdef USE_ARPACK
  use linalg, only: pardiso_c
#endif

  use linalg, only: zheev
  use hamiltonianConstructor, only: ZB8bandQW

  implicit none
  private

  public :: compute_spectral_function_qw

#ifdef USE_ARPACK
  public :: compute_ldos_csr
#endif

contains

  ! ==============================================================================
  ! Compute spectral function A(k, E) for a QW system.
  !
  ! A(k, E) = sum_n (1/pi) * eta / ((E - E_n(k))^2 + eta^2)
  !
  ! For each k-point the 8N x 8N QW Hamiltonian is built and diagonalized
  ! with zheev. The resulting eigenvalues are broadened with a Lorentzian
  ! of half-width eta.
  !
  ! Args:
  !   cfg      -- simulation configuration (confinement, grid, etc.)
  !   profile  -- material profile array (N x 3)
  !   kpterms  -- k.p parameter terms (N x N x 10)
  !   k_arr    -- k-points at which to compute A (nk)
  !   E_arr    -- energy grid (nE)
  !   eta      -- Lorentzian broadening half-width (eV)
  !   A_kE     -- output spectral function (nk x nE)
  ! ==============================================================================
  subroutine compute_spectral_function_qw(cfg, profile, kpterms, k_arr, E_arr, &
                                           eta, A_kE)
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: profile(:,:), kpterms(:,:,:)
    real(kind=dp), contiguous, intent(in) :: k_arr(:), E_arr(:)
    real(kind=dp), intent(in) :: eta
    real(kind=dp), allocatable, intent(out) :: A_kE(:,:)

    integer :: Ngrid, nk, nE, ik, iE, lwork, info, i
    integer :: dim
    complex(kind=dp), allocatable :: H(:,:), work(:)
    real(kind=dp), allocatable :: evals(:), rwork(:)
    real(kind=dp) :: lorntz
    type(wavevector) :: wv

    Ngrid = size(profile, 1)
    dim = 8 * Ngrid
    nk = size(k_arr)
    nE = size(E_arr)

    allocate(A_kE(nk, nE))
    A_kE = 0.0_dp

    allocate(H(dim, dim))
    allocate(evals(dim))
    allocate(rwork(max(1, 3*dim - 2)))
    allocate(work(1))

    ! Query optimal work size
    call zheev('N', 'U', dim, H, dim, evals, work, -1, rwork, info)
    lwork = nint(real(work(1)))
    deallocate(work)
    allocate(work(max(1, lwork)))

    do ik = 1, nk
      wv%kx = k_arr(ik)
      wv%ky = 0.0_dp
      wv%kz = 0.0_dp

      H = cmplx(0.0_dp, 0.0_dp, kind=dp)
      call ZB8bandQW(H, wv, profile, kpterms, cfg=cfg)

      ! Diagonalize
      call zheev('N', 'U', dim, H, dim, evals, work, size(work), rwork, info)

      ! Compute A(k, E) = sum_n delta_eta(E - E_n)
      do iE = 1, nE
        do i = 1, dim
          lorntz = eta / (pi_dp * ((E_arr(iE) - evals(i))**2 + eta**2))
          A_kE(ik, iE) = A_kE(ik, iE) + lorntz
        end do
      end do
    end do

    deallocate(H, evals, rwork, work)
  end subroutine compute_spectral_function_qw

#ifdef USE_ARPACK

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
