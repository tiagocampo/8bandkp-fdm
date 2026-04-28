module eigensolver

  use definitions, only: dp
  use sparse_matrices, only: csr_matrix
  use linalg, only: zheevx
#ifdef USE_MKL_FEAST
  use linalg, only: feastinit, zfeast_hcsrev
#endif
  implicit none

  private
  public :: eigensolver_config, eigensolver_result
  public :: feast_workspace, feast_workspace_free
  public :: solve_sparse_evp, solve_feast, solve_dense_lapack
  public :: auto_compute_energy_window, eigensolver_result_free
#ifdef USE_ARPACK
  public :: solve_arpack
#endif

  ! ------------------------------------------------------------------
  ! Configuration for the sparse eigensolver.
  ! ------------------------------------------------------------------
  type :: eigensolver_config
    character(len=10)  :: method  = 'FEAST'
    integer            :: nev     = 8
    real(kind=dp)      :: emin    = -1.0_dp
    real(kind=dp)      :: emax    =  1.0_dp
    integer            :: max_iter = 100
    real(kind=dp)      :: tol     = 1.0e-10_dp
    integer            :: feast_m0 = 0
    integer            :: ncv     = 0
    character(len=1)   :: which   = 'S'
  end type eigensolver_config

  ! ------------------------------------------------------------------
  ! Result from sparse eigensolver.
  ! ------------------------------------------------------------------
  type :: eigensolver_result
    integer                       :: nev_found = 0
    real(kind=dp), allocatable    :: eigenvalues(:)
    complex(kind=dp), allocatable :: eigenvectors(:,:)
    integer                       :: iterations = 0
    logical                       :: converged = .false.
  end type eigensolver_result

  ! ------------------------------------------------------------------
  ! Cached upper-triangle CSR structure for repeated FEAST calls.
  ! Avoids O(NNZ) scan + 3 allocations per k-point in wire kz-sweeps.
  ! ------------------------------------------------------------------
  type :: feast_workspace
    integer, allocatable          :: rowptr_loc(:)   ! (N+1) upper-triangle row pointers
    integer, allocatable          :: colind_loc(:)   ! (nnz_upper) column indices
    integer                       :: nnz_upper = 0
    integer                       :: M0 = 0
    integer                       :: N = 0
    logical                       :: initialized = .false.
  end type feast_workspace

contains

  ! ==================================================================
  ! Check whether a cached feast_workspace matches the upper-triangle
  ! sparsity pattern of the current CSR matrix.
  ! ==================================================================
  logical function feast_workspace_matches_pattern(H_csr, fw, N, M0)
    type(csr_matrix), intent(in) :: H_csr
    type(feast_workspace), intent(in) :: fw
    integer, intent(in) :: N, M0

    integer :: i, j, k

    feast_workspace_matches_pattern = .false.
    if (.not. fw%initialized) return
    if (fw%N /= N .or. fw%M0 /= M0) return
    if (.not. allocated(fw%rowptr_loc)) return
    if (.not. allocated(fw%colind_loc)) return
    if (size(fw%rowptr_loc) /= N + 1) return

    do i = 1, N
      k = fw%rowptr_loc(i)
      do j = H_csr%rowptr(i), H_csr%rowptr(i+1) - 1
        if (H_csr%colind(j) >= i) then
          if (k >= fw%rowptr_loc(i + 1)) return
          if (k > fw%nnz_upper) return
          if (fw%colind_loc(k) /= H_csr%colind(j)) return
          k = k + 1
        end if
      end do
      if (k /= fw%rowptr_loc(i + 1)) return
    end do

    feast_workspace_matches_pattern = .true.
  end function feast_workspace_matches_pattern

  ! ==================================================================
  ! Main dispatch: select FEAST or ARPACK based on config.
  ! ==================================================================
  subroutine solve_sparse_evp(H_csr, config, result, feast_ws)
    type(csr_matrix), intent(in)          :: H_csr
    type(eigensolver_config), intent(in)  :: config
    type(eigensolver_result), intent(out) :: result
    type(feast_workspace), intent(inout), optional :: feast_ws

    select case (trim(config%method))
#ifdef USE_MKL_FEAST
    case ('FEAST')
      print *, '  Eigensolver: FEAST'
      call solve_feast(H_csr, config, result, fw=feast_ws)
    case ('DENSE')
      print *, '  Eigensolver: dense LAPACK'
      call solve_dense_lapack(H_csr, config, result)
    case ('ARPACK')
      call solve_arpack_dispatch(H_csr, config, result)
#else
    case ('FEAST', 'ARPACK', 'DENSE')
      ! FEAST unavailable: fall back to dense LAPACK
      print *, 'Using dense LAPACK eigensolver.'
      call solve_dense_lapack(H_csr, config, result)
#endif
    case default
      print *, 'Error: unknown eigensolver method "', trim(config%method), '"'
      result%converged = .false.
      result%nev_found = 0
    end select
  end subroutine solve_sparse_evp

#ifdef USE_MKL_FEAST
  ! ==================================================================
  ! MKL FEAST wrapper for complex Hermitian CSR matrix.
  !
  ! Uses zfeast_hcsrev which finds all eigenvalues in [emin, emax].
  ! CSR must use 1-based indexing.  UPLO='U' (upper triangle).
  !
  ! Explicit interface from linalg module provides proper type signatures
  ! so gfortran generates correct calls without wrapper indirection.
  ! ==================================================================
  subroutine solve_feast(H_csr, config, result, fw)
    type(csr_matrix), intent(in)          :: H_csr
    type(eigensolver_config), intent(in)  :: config
    type(eigensolver_result), intent(out) :: result
    type(feast_workspace), intent(inout), optional :: fw

    integer :: N, M0, M, loop, info, i, j, k, nnz
    integer :: fpm(128)
    real(kind=dp) :: epsout
    real(kind=dp), allocatable  :: E(:), res(:)
    complex(kind=dp), allocatable :: X(:,:)
    ! Contiguous local copies for FEAST (avoids allocatable component issues)
    complex(kind=dp), allocatable :: val_loc(:)
    integer, allocatable :: rowptr_loc(:), colind_loc(:)

    N = H_csr%nrows
    if (N <= 0) then
      result%converged = .false.
      result%nev_found = 0
      return
    end if

    call feastinit(fpm)
    fpm(4) = config%max_iter  ! max FEAST refinement iterations

    if (config%tol > 0.0_dp) then
      fpm(3) = max(1, min(15, int(-log10(config%tol))))
    end if

    fpm(1) = 0  ! silent mode
    fpm(2) = 16 ! contour quadrature points (12-16 recommended for degenerate spectra)

    M0 = config%feast_m0
    if (M0 <= 0) M0 = 2 * config%nev
    M0 = max(M0, config%nev + 1)
    M0 = min(M0, N)  ! FEAST requires M0 <= N

    allocate(E(M0))
    allocate(X(N, M0))
    allocate(res(M0))
    E = 0.0_dp
    X = cmplx(0.0_dp, 0.0_dp, kind=dp)
    res = 0.0_dp

    ! Extract upper triangle only (FEAST UPLO='U' requires col >= row).
    if (present(fw) .and. feast_workspace_matches_pattern(H_csr, fw, N, M0)) then
      ! Fast path: reuse cached upper-triangle structure, just update values
      allocate(val_loc(fw%nnz_upper))
      do i = 1, N
        k = fw%rowptr_loc(i)
        do j = H_csr%rowptr(i), H_csr%rowptr(i+1) - 1
          if (H_csr%colind(j) >= i) then
            val_loc(k) = H_csr%values(j)
            k = k + 1
          end if
        end do
      end do
      ! Use cached rowptr and colind directly
      call zfeast_hcsrev('U', N, val_loc, fw%rowptr_loc, fw%colind_loc, &
                         fpm, epsout, loop, config%emin, config%emax, M0, &
                         E, X, M, res, info)
      deallocate(val_loc)
    else
      ! Free stale cache before rebuilding
      if (present(fw)) then
        if (fw%initialized) call feast_workspace_free(fw)
      end if
      ! Original path: count, allocate, fill
      allocate(rowptr_loc(N+1))
      rowptr_loc = 0
      do i = 1, N
        do k = H_csr%rowptr(i), H_csr%rowptr(i+1) - 1
          if (H_csr%colind(k) >= i) then
            rowptr_loc(i+1) = rowptr_loc(i+1) + 1
          end if
        end do
      end do
      rowptr_loc(1) = 1
      do i = 1, N
        rowptr_loc(i+1) = rowptr_loc(i) + rowptr_loc(i+1)
      end do
      nnz = rowptr_loc(N+1) - 1

      allocate(val_loc(nnz), colind_loc(nnz))
      do i = 1, N
        k = rowptr_loc(i)
        do j = H_csr%rowptr(i), H_csr%rowptr(i+1) - 1
          if (H_csr%colind(j) >= i) then
            val_loc(k) = H_csr%values(j)
            colind_loc(k) = H_csr%colind(j)
            k = k + 1
          end if
        end do
      end do

      call zfeast_hcsrev('U', N, val_loc, rowptr_loc, colind_loc, &
                         fpm, epsout, loop, config%emin, config%emax, M0, &
                         E, X, M, res, info)

      ! Cache for future calls
      if (present(fw)) then
        call move_alloc(rowptr_loc, fw%rowptr_loc)
        call move_alloc(colind_loc, fw%colind_loc)
        fw%nnz_upper = nnz
        fw%M0 = M0
        fw%N = N
        fw%initialized = .true.
      else
        deallocate(rowptr_loc, colind_loc)
      end if

      deallocate(val_loc)
    end if

    result%iterations = loop
    result%converged = (info == 0) .or. (info == 2)

    if (info < 0) then
      print *, 'FEAST error: info =', info
    else if (info == 2) then
      print *, 'FEAST warning: max iterations reached, results may be inaccurate.'
    else if (info == 3) then
      print *, '  WARNING: FEAST subspace too small (info=3). Missing eigenvalues.'
      print *, '  Increase feast_m0 or narrow the energy window.'
    end if

    if (M > 0 .and. M <= M0 .and. info >= 0) then
      result%nev_found = M
      allocate(result%eigenvalues(M))
      allocate(result%eigenvectors(N, M))
      do i = 1, M
        result%eigenvalues(i) = E(i)
        result%eigenvectors(:, i) = X(:, i)
      end do
      call sort_eigenpairs_ascending(result%eigenvalues, result%eigenvectors)
    else
      result%nev_found = 0
    end if

    deallocate(E, X, res)
  end subroutine solve_feast
#endif /* USE_MKL_FEAST */

  ! ==================================================================
  ! ARPACK-NG solver dispatch.
  !
  ! When USE_ARPACK is defined, calls the proper ARPACK solver.
  ! Otherwise falls back to dense LAPACK.
  ! ==================================================================
  subroutine solve_arpack_dispatch(H_csr, config, result)
    type(csr_matrix), intent(in)          :: H_csr
    type(eigensolver_config), intent(in)  :: config
    type(eigensolver_result), intent(out) :: result

#ifdef USE_ARPACK
    print *, '  Eigensolver: ARPACK-NG (shift-invert)'
    call solve_arpack(H_csr, config, result)
#else
    print *, '  ARPACK not available, using dense LAPACK fallback.'
    call solve_dense_lapack(H_csr, config, result)
#endif
  end subroutine solve_arpack_dispatch

#ifdef USE_ARPACK
  ! ==================================================================
  ! ARPACK-NG solver: shift-invert mode with MKL PARDISO.
  !
  ! Solves H*x = lambda*x for eigenvalues near sigma using:
  !   OP = (H - sigma*I)^{-1}
  !
  ! ARPACK finds the largest-magnitude eigenvalues of OP, which
  ! correspond to eigenvalues of H closest to sigma.
  !
  ! Phase 1: PARDISO factorisation of (H - sigma*I)
  ! Phase 2: ARPACK reverse-communication loop (OP*x via PARDISO solve)
  ! Phase 3: Extract Ritz values/vectors via zneupd
  ! ==================================================================
  subroutine solve_arpack(H_csr, config, result)
    use sparse_matrices, only: csr_spmv
    use linalg, only: znaupd, zneupd
    type(csr_matrix), intent(in)          :: H_csr
    type(eigensolver_config), intent(in)  :: config
    type(eigensolver_result), intent(out) :: result

    integer :: N, nev_want, ncv_loc, ldv, ldz, lworkl
    integer :: ido, info, i, j, nconv
    real(kind=dp) :: tol_loc, sigma_re
    complex(kind=dp) :: sigma

    ! ARPACK arrays
    complex(kind=dp), allocatable :: resid(:), v(:,:), workd(:), workl(:)
    complex(kind=dp), allocatable :: workev(:), d(:), z(:,:)
    real(kind=dp), allocatable :: rwork(:)
    logical, allocatable :: select(:)
    integer :: iparam(11), ipntr(14)

    ! PARDISO arrays for shift-invert
    integer(8) :: pt(64)
    integer :: iparm(64), maxfct, mnum, mtype, phase, nrhs, msglvl, error_loc
    complex(kind=dp), allocatable :: H_shifted_val(:)
    integer, allocatable :: H_shifted_rowptr(:), H_shifted_colind(:)
    complex(kind=dp), allocatable :: rhs(:), sol(:)
    complex(kind=dp) :: dummy_bx(1)
    external :: pardiso
    integer :: nnz, k
    integer, allocatable :: perm(:)

    N = H_csr%nrows
    if (N <= 0) then
      result%converged = .false.
      result%nev_found = 0
      return
    end if

    nev_want = min(config%nev, N - 1)
    if (nev_want <= 0) then
      result%converged = .false.
      result%nev_found = 0
      return
    end if

    ! Shift = midpoint of energy window
    sigma_re = 0.5_dp * (config%emin + config%emax)
    sigma = cmplx(sigma_re, 0.0_dp, kind=dp)

    ! --- Phase 1: Build (H - sigma*I) and factorize with PARDISO ---

    ! Copy CSR and subtract sigma from diagonal
    nnz = H_csr%rowptr(N + 1) - H_csr%rowptr(1)
    allocate(H_shifted_val(nnz), H_shifted_rowptr(N+1), H_shifted_colind(nnz))
    H_shifted_val = H_csr%values
    H_shifted_colind = H_csr%colind
    H_shifted_rowptr = H_csr%rowptr

    ! Shift diagonal
    do i = 1, N
      do k = H_shifted_rowptr(i), H_shifted_rowptr(i+1) - 1
        if (H_shifted_colind(k) == i) then
          H_shifted_val(k) = H_shifted_val(k) - sigma
          exit
        end if
      end do
    end do

    ! PARDISO init
    pt = 0
    iparm = 0
    iparm(1) = 1   ! no solver default
    iparm(2) = 2   ! OpenMP nested dissection
    iparm(3) = 1   ! reserved
    iparm(4) = 0   ! no CG iterations
    iparm(8) = 2   ! max iterative refinement steps
    iparm(10) = 13 ! perturb pivots with 1E-13
    iparm(11) = 1  ! scaling
    iparm(13) = 1  ! matching
    iparm(18) = -1 ! report number of nonzeros
    iparm(19) = -1 ! report flop count
    iparm(27) = 1  ! check matrix consistency
    iparm(40) = 1  ! distributed matrix input (CSR)
    maxfct = 1
    mnum = 1
    mtype = 3  ! complex structurally symmetric (full CSR, not triangle-only)
    nrhs = 1
    msglvl = 0  ! no output

    allocate(perm(N))
    perm = 0

    ! Symbolic factorization (phase=11)
    phase = 11
    dummy_bx = cmplx(0.0_dp, 0.0_dp, kind=dp)
    call pardiso(pt, maxfct, mnum, mtype, phase, N, &
                 H_shifted_val, H_shifted_rowptr, H_shifted_colind, &
                 perm, nrhs, iparm, msglvl, &
                 dummy_bx, dummy_bx, error_loc)
    if (error_loc /= 0) then
      print *, 'ARPACK/PARDISO: symbolic factorization error', error_loc
      result%converged = .false.
      result%nev_found = 0
      deallocate(H_shifted_val, H_shifted_rowptr, H_shifted_colind, perm)
      return
    end if

    ! Numeric factorization (phase=22)
    phase = 22
    call pardiso(pt, maxfct, mnum, mtype, phase, N, &
                 H_shifted_val, H_shifted_rowptr, H_shifted_colind, &
                 perm, nrhs, iparm, msglvl, &
                 dummy_bx, dummy_bx, error_loc)
    if (error_loc /= 0) then
      print *, 'ARPACK/PARDISO: numeric factorization error', error_loc
      result%converged = .false.
      result%nev_found = 0
      ! Cleanup PARDISO
      phase = -1
      call pardiso(pt, maxfct, mnum, mtype, phase, N, &
                   H_shifted_val, H_shifted_rowptr, H_shifted_colind, &
                   perm, nrhs, iparm, msglvl, &
                   dummy_bx, dummy_bx, error_loc)
      deallocate(H_shifted_val, H_shifted_rowptr, H_shifted_colind, perm)
      return
    end if

    ! --- Phase 2: ARPACK reverse communication loop ---

    ! Set ARPACK parameters
    ncv_loc = config%ncv
    if (ncv_loc <= 0) ncv_loc = min(max(2 * nev_want + 1, 20), N)
    ncv_loc = max(ncv_loc, nev_want + 2)

    ldv = N
    ldz = N
    lworkl = 3 * ncv_loc * ncv_loc + 5 * ncv_loc
    tol_loc = config%tol
    if (tol_loc <= 0.0_dp) tol_loc = 0.0_dp  ! machine precision

    allocate(resid(N))
    allocate(v(ldv, ncv_loc))
    allocate(workd(3 * N))
    allocate(workl(lworkl))
    allocate(rwork(ncv_loc))
    allocate(workev(2 * ncv_loc))
    allocate(select(ncv_loc))
    allocate(d(nev_want))
    allocate(z(ldz, nev_want))
    allocate(rhs(N), sol(N))

    ! Initialize
    resid = cmplx(0.0_dp, 0.0_dp, kind=dp)
    v = cmplx(0.0_dp, 0.0_dp, kind=dp)
    workd = cmplx(0.0_dp, 0.0_dp, kind=dp)
    workl = cmplx(0.0_dp, 0.0_dp, kind=dp)
    rwork = 0.0_dp
    select = .false.

    iparam = 0
    iparam(1) = 1   ! exact shifts
    iparam(3) = config%max_iter
    iparam(4) = 1   ! block size (must be 1)
    iparam(7) = 3   ! shift-invert mode: OP = (A - sigma*I)^{-1}

    ido = 0
    info = 0  ! random initial vector

    ! Reverse communication loop
    do while (ido /= 99)
      call znaupd(ido, 'I', N, 'LM', nev_want, tol_loc, resid, ncv_loc, &
                  v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info)

      select case (ido)
      case (-1, 1)
        ! Compute Y = OP * X = (H - sigma*I)^{-1} * X
        ! X = workd(ipntr(1) : ipntr(1)+N-1)
        ! Y = workd(ipntr(2) : ipntr(2)+N-1)
        rhs(1:N) = workd(ipntr(1) : ipntr(1) + N - 1)

        ! PARDISO solve (phase=33)
        phase = 33
        call pardiso(pt, maxfct, mnum, mtype, phase, N, &
                     H_shifted_val, H_shifted_rowptr, H_shifted_colind, &
                     perm, nrhs, iparm, msglvl, rhs, sol, error_loc)

        if (error_loc /= 0) then
          print *, 'ARPACK/PARDISO: solve error', error_loc
          info = -999
          exit
        end if

        workd(ipntr(2) : ipntr(2) + N - 1) = sol(1:N)

      case (2)
        ! B*x for standard problem (bmat='I'), should not happen
        ! Just copy x to y
        workd(ipntr(2) : ipntr(2) + N - 1) = &
          workd(ipntr(1) : ipntr(1) + N - 1)

      case (3)
        ! User-supplied shifts (iparam(1)=1 means ARPACK handles this)
        ! Should not reach here with default settings
        continue

      end select
    end do

    result%iterations = iparam(3)

    ! --- Phase 3: Extract Ritz values and vectors ---

    if (info == 0 .and. ido == 99) then
      call zneupd(.true., 'A', select, d, z, ldz, sigma, workev, &
                  'I', N, 'LM', nev_want, tol_loc, resid, ncv_loc, &
                  v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info)

      if (info == 0) then
        nconv = iparam(5)
        result%converged = .true.
        result%nev_found = nev_want
        allocate(result%eigenvalues(nev_want))
        allocate(result%eigenvectors(N, nev_want))
        ! For Hermitian matrix, eigenvalues are real (discard tiny imaginary parts)
        do i = 1, nev_want
          result%eigenvalues(i) = real(d(i), kind=dp)
          result%eigenvectors(:, i) = z(:, i)
        end do

        ! Sort eigenvalues (and corresponding eigenvectors) in ascending order
        call sort_eigenpairs_ascending(result%eigenvalues, result%eigenvectors)
      else
        print *, 'ARPACK zneupd error: info =', info
        result%converged = .false.
        result%nev_found = 0
      end if
    else
      if (info == 1) then
        print *, 'ARPACK: max iterations reached. Results may be inaccurate.'
      else if (info /= 0) then
        print *, 'ARPACK znaupd error: info =', info
      end if
      result%converged = .false.
      result%nev_found = 0
    end if

    ! Cleanup PARDISO
    phase = -1
    call pardiso(pt, maxfct, mnum, mtype, phase, N, &
                 H_shifted_val, H_shifted_rowptr, H_shifted_colind, &
                 perm, nrhs, iparm, msglvl, &
                 dummy_bx, dummy_bx, error_loc)

    deallocate(resid, v, workd, workl, rwork, workev, select, d, z)
    deallocate(rhs, sol)
    deallocate(H_shifted_val, H_shifted_rowptr, H_shifted_colind, perm)
  end subroutine solve_arpack
#endif /* USE_ARPACK */

  ! ==================================================================
  ! Dense fallback using LAPACK zheevx.
  !
  ! Converts CSR to dense, then calls zheevx for the N-smallest
  ! eigenvalues.  Used as ARPACK placeholder until proper ARPACK
  ! linking is set up.
  ! ==================================================================
  subroutine solve_dense_lapack(H_csr, config, result)
    type(csr_matrix), intent(in)          :: H_csr
    type(eigensolver_config), intent(in)  :: config
    type(eigensolver_result), intent(out) :: result

    integer :: N, lda, ldz, lwork, info, nev_want, nb
    complex(kind=dp), allocatable :: A(:,:), Z(:,:), work(:)
    real(kind=dp), allocatable :: rwork(:), W(:)
    integer, allocatable :: iwork(:), ifail(:)
    real(kind=dp) :: vl, vu, abstol

    N = H_csr%nrows
    if (N <= 0) then
      result%converged = .false.
      result%nev_found = 0
      return
    end if

    nev_want = min(config%nev, N)
    lda = N
    ldz = N
    abstol = 0.0_dp

    ! Convert CSR to dense
    allocate(A(N, N))
    A = cmplx(0.0_dp, 0.0_dp, kind=dp)
    call csr_to_dense_work(H_csr, A, N)

    ! Allocate output arrays — Z must hold all eigenvalues in range mode
    allocate(W(N))
    allocate(Z(N, N))
    allocate(rwork(7 * N))
    allocate(iwork(5 * N))
    allocate(ifail(N))

    ! Use range mode 'V' when emin/emax are set (non-zero), otherwise index mode 'I'
    if (config%emin /= 0.0_dp .and. config%emax /= 0.0_dp) then
      ! Range mode: return eigenvalues in [emin, emax]
      vl = config%emin
      vu = config%emax
      ! Workspace query
      allocate(work(1))
      call zheevx('V', 'V', 'U', N, A, lda, vl, vu, &
                  1, N, abstol, nb, W, Z, ldz, &
                  work, -1, rwork, iwork, ifail, info)
      lwork = max(1, nint(real(work(1))))
      deallocate(work)
      allocate(work(lwork))

      ! Solve
      call zheevx('V', 'V', 'U', N, A, lda, vl, vu, &
                  1, N, abstol, nb, W, Z, ldz, &
                  work, lwork, rwork, iwork, ifail, info)
      ! Return all eigenvalues in the energy window (no trimming)
      ! Range mode already selects the desired eigenvalue range.
    else
      ! Index mode: return nev_want smallest eigenvalues
      ! Workspace query
      allocate(work(1))
      call zheevx('V', 'I', 'U', N, A, lda, vl, vu, &
                  1, nev_want, abstol, nb, W, Z, ldz, &
                  work, -1, rwork, iwork, ifail, info)
      lwork = max(1, nint(real(work(1))))
      deallocate(work)
      allocate(work(lwork))

      ! Solve
      call zheevx('V', 'I', 'U', N, A, lda, vl, vu, &
                  1, nev_want, abstol, nb, W, Z, ldz, &
                  work, lwork, rwork, iwork, ifail, info)
    end if

    if (info == 0 .and. nb > 0) then
      result%converged = .true.
      result%nev_found = nb
      result%iterations = 1
      allocate(result%eigenvalues(nb))
      allocate(result%eigenvectors(N, nb))
      result%eigenvalues(1:nb) = W(1:nb)
      result%eigenvectors(:, 1:nb) = Z(:, 1:nb)
    else
      result%converged = .false.
      result%nev_found = 0
      if (info /= 0) then
        print *, 'Dense eigensolver error: info =', info
      end if
    end if

    deallocate(A, W, Z, work, rwork, iwork, ifail)
  end subroutine solve_dense_lapack

  ! ==================================================================
  ! Gershgorin energy window estimation.
  !
  ! Row-sum bounds: emin <= lambda_min, emax >= lambda_max.
  ! 10% margin for safety.
  ! ==================================================================
  subroutine auto_compute_energy_window(H_csr, emin, emax)
    type(csr_matrix), intent(in)  :: H_csr
    real(kind=dp), intent(out) :: emin, emax

    integer :: row, k, N
    real(kind=dp) :: row_sum, diag_val, rmin, rmax
    real(kind=dp), parameter :: margin_frac = 0.1_dp
    real(kind=dp), parameter :: margin_floor = 0.5_dp  ! eV

    N = H_csr%nrows
    if (N == 0) then
      emin = -1.0_dp
      emax = 1.0_dp
      return
    end if

    rmin = 0.0_dp
    rmax = 0.0_dp

    do row = 1, N
      diag_val = 0.0_dp
      row_sum = 0.0_dp
      do k = H_csr%rowptr(row), H_csr%rowptr(row + 1) - 1
        if (H_csr%colind(k) == row) then
          diag_val = real(H_csr%values(k), kind=dp)
        else
          row_sum = row_sum + abs(H_csr%values(k))
        end if
      end do
      if (row == 1) then
        rmin = diag_val - row_sum
        rmax = diag_val + row_sum
      else
        rmin = min(rmin, diag_val - row_sum)
        rmax = max(rmax, diag_val + row_sum)
      end if
    end do

    emin = rmin - max(margin_frac * abs(rmin), margin_floor)
    emax = rmax + max(margin_frac * abs(rmax), margin_floor)
  end subroutine auto_compute_energy_window

  ! ==================================================================
  ! Deallocate result arrays.
  ! ==================================================================
  subroutine eigensolver_result_free(result)
    type(eigensolver_result), intent(inout) :: result

    if (allocated(result%eigenvalues))  deallocate(result%eigenvalues)
    if (allocated(result%eigenvectors)) deallocate(result%eigenvectors)
    result%nev_found = 0
    result%converged = .false.
    result%iterations = 0
  end subroutine eigensolver_result_free

  ! ==================================================================
  ! Deallocate feast_workspace cached arrays.
  ! ==================================================================
  subroutine feast_workspace_free(fw)
    type(feast_workspace), intent(inout) :: fw

    if (allocated(fw%rowptr_loc)) deallocate(fw%rowptr_loc)
    if (allocated(fw%colind_loc)) deallocate(fw%colind_loc)
    fw%nnz_upper = 0
    fw%M0 = 0
    fw%N = 0
    fw%initialized = .false.
  end subroutine feast_workspace_free

  ! ==================================================================
  ! Internal: convert CSR to dense.
  ! ==================================================================
  subroutine csr_to_dense_work(H_csr, A, N)
    type(csr_matrix), intent(in)  :: H_csr
    integer, intent(in) :: N
    complex(kind=dp), intent(out) :: A(N, N)

    integer :: row, k

    A = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do row = 1, N
      do k = H_csr%rowptr(row), H_csr%rowptr(row + 1) - 1
        A(row, H_csr%colind(k)) = H_csr%values(k)
      end do
    end do
  end subroutine csr_to_dense_work


  ! ------------------------------------------------------------------
  ! Sort eigenvalues in ascending order, rearranging eigenvectors to match.
  ! Simple selection sort (sufficient for typical nev << N).
  ! ------------------------------------------------------------------
  subroutine sort_eigenpairs_ascending(evals, evecs)
    real(kind=dp), intent(inout) :: evals(:)
    complex(kind=dp), intent(inout) :: evecs(:,:)

    integer :: i, j, jmin, n
    real(kind=dp) :: tmp_eval
    complex(kind=dp), allocatable :: tmp_vec(:)

    n = size(evals)
    if (n <= 1) return

    allocate(tmp_vec(size(evecs, 1)))

    do i = 1, n - 1
      jmin = i
      do j = i + 1, n
        if (evals(j) < evals(jmin)) jmin = j
      end do
      if (jmin /= i) then
        tmp_eval = evals(i)
        evals(i) = evals(jmin)
        evals(jmin) = tmp_eval
        tmp_vec = evecs(:, i)
        evecs(:, i) = evecs(:, jmin)
        evecs(:, jmin) = tmp_vec
      end if
    end do

    deallocate(tmp_vec)
  end subroutine sort_eigenpairs_ascending

end module eigensolver
