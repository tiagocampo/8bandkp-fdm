module eigensolver

  use definitions, only: dp
  use sparse_matrices, only: csr_matrix, csr_free, csr_build_from_coo
  use linalg, only: zheevx, zheev
#ifdef USE_MKL_FEAST
  use linalg, only: feastinit, zfeast_hcsrev
#endif
  implicit none

  private
  public :: eigensolver_config, eigensolver_result
  public :: feast_workspace, feast_workspace_free
  public :: solve_feast
  public :: auto_compute_energy_window, eigensolver_result_free
  public :: eigensolver_base, dense_lapack_solver_t
#ifdef USE_MKL_FEAST
  public :: feast_solver_t
#endif
  public :: make_eigensolver
  public :: eigensolver_config_validate
  public :: EIGEN_MODE_FULL, EIGEN_MODE_INDEX, EIGEN_MODE_ENERGY

  ! ------------------------------------------------------------------
  ! Mode constants for eigensolver dispatch.
  ! ------------------------------------------------------------------
  integer, parameter :: EIGEN_MODE_FULL   = 1  ! all eigenvalues
  integer, parameter :: EIGEN_MODE_INDEX  = 2  ! eigenvalues il:iu
  integer, parameter :: EIGEN_MODE_ENERGY = 3  ! eigenvalues in [emin, emax]

  ! ------------------------------------------------------------------
  ! Configuration for the sparse eigensolver.
  ! ------------------------------------------------------------------
  type :: eigensolver_config
    character(len=10)  :: method   = 'DENSE'
    integer            :: mode     = EIGEN_MODE_FULL   ! default: FULL (valid for all methods)
    integer            :: nev      = 8
    integer            :: il       = 1
    integer            :: iu       = 8
    real(kind=dp)      :: emin     = -1.0_dp
    real(kind=dp)      :: emax     =  1.0_dp
    integer            :: max_iter = 100
    real(kind=dp)      :: tol      = 1.0e-10_dp
    integer            :: m0 = 0
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
    logical                       :: was_freed = .false.
  contains
    final :: eigensolver_result_finalize
  end type eigensolver_result

  ! ------------------------------------------------------------------
  ! Abstract base type for polymorphic eigensolver dispatch.
  ! ------------------------------------------------------------------
  type, abstract :: eigensolver_base
  contains
    procedure(solve_evp_interface), deferred :: solve
    procedure(solve_dense_interface), deferred :: solve_dense
    procedure(solve_sparse_interface), deferred :: solve_sparse
  end type eigensolver_base

  abstract interface
    subroutine solve_evp_interface(self, H_csr, config, result)
      import :: eigensolver_base, csr_matrix, eigensolver_config, eigensolver_result
      class(eigensolver_base), intent(inout) :: self
      type(csr_matrix), intent(in) :: H_csr
      type(eigensolver_config), intent(in) :: config
      type(eigensolver_result), intent(out) :: result
    end subroutine solve_evp_interface

    subroutine solve_dense_interface(self, H, config, result)
      import :: eigensolver_base, dp, eigensolver_config, eigensolver_result
      class(eigensolver_base), intent(inout) :: self
      complex(kind=dp), contiguous, intent(in) :: H(:,:)
      type(eigensolver_config), intent(in) :: config
      type(eigensolver_result), intent(out) :: result
    end subroutine solve_dense_interface

    subroutine solve_sparse_interface(self, H_csr, config, result)
      import :: eigensolver_base, csr_matrix, eigensolver_config, eigensolver_result
      class(eigensolver_base), intent(inout) :: self
      type(csr_matrix), intent(in) :: H_csr
      type(eigensolver_config), intent(in) :: config
      type(eigensolver_result), intent(out) :: result
    end subroutine solve_sparse_interface
  end interface

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
    logical                       :: was_freed = .false.
    ! Persisted subspace seed: the M0 that last produced a converged solve.
    ! A subsequent call seeds M0 from this value (capped at N, >= nev+1) so the
    ! info=3 retry tax is not re-paid on every k-point of a sweep where the
    ! eigenvalue count in the fixed energy window is ~constant. The existing
    ! info=3 retry loop remains as the backstop for any point needing more.
    integer                       :: last_successful_m0 = 0
  contains
    final :: feast_workspace_finalize
  end type feast_workspace

  ! ------------------------------------------------------------------
  ! Concrete solver types for polymorphic dispatch.
  ! ------------------------------------------------------------------
#ifdef USE_MKL_FEAST
  type, extends(eigensolver_base) :: feast_solver_t
    type(feast_workspace) :: ws
  contains
    procedure :: solve => feast_solve_dispatch        ! legacy alias -> solve_sparse
    procedure :: solve_dense => feast_solve_dense
    procedure :: solve_sparse => feast_solve_sparse_dispatch
    final :: feast_solver_finalize
  end type feast_solver_t
#endif

  type, extends(eigensolver_base) :: dense_lapack_solver_t
    ! Cached LAPACK workspace — sized on first call (or when N grows),
    ! reused thereafter. Thread-safe because each OpenMP thread allocates
    ! its own solver instance (main.f90 QW sweep).
    integer                       :: cached_n = 0
    complex(kind=dp), allocatable :: A_buf(:,:), Z_buf(:,:), work(:)
    real(kind=dp), allocatable    :: W_buf(:), rwork(:)
    integer, allocatable          :: iwork(:), ifail(:)
  contains
    procedure :: solve => dense_lapack_solve_dispatch       ! legacy alias -> solve_sparse
    procedure :: solve_dense => dense_solve_dense_dispatch
    procedure :: solve_sparse => dense_solve_sparse_dispatch
    final     :: dense_lapack_solver_finalize
  end type dense_lapack_solver_t

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
    logical :: fw_match
    ! Retry loop for info=3 (subspace too small)
    integer, parameter :: MAX_RETRY = 3
    integer :: retry, M0_initial
    logical :: cache_is_fresh

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

    M0 = config%m0
    if (M0 <= 0) M0 = 2 * config%nev
    M0 = max(M0, config%nev + 1)
    M0 = min(M0, N)  ! FEAST requires M0 <= N

    ! Seed from a previously-converged subspace size (if any) so a k-sweep does
    ! not re-pay the info=3 retry tax on every point. Never smaller than the
    ! user's setting (M0 is already >= that here), never larger than N. The
    ! persisted seed is only used when the caller passed a reusable workspace.
    if (present(fw)) then
      if (fw%last_successful_m0 > M0) then
        M0 = fw%last_successful_m0
        M0 = max(M0, config%nev + 1)
        M0 = min(M0, N)
      end if
    end if

    M0_initial = M0

    ! ------------------------------------------------------------------
    ! Retry loop: on info=3 (subspace too small), double M0 and retry.
    ! Free cached workspace on retry since the pattern no longer matches.
    ! ------------------------------------------------------------------
    cache_is_fresh = .false.
    do retry = 1, MAX_RETRY
      if (retry > 1) then
        ! Double M0, capped at N
        M0 = min(2 * M0, N)
        print *, '  FEAST info=3 retry ', retry, '/', MAX_RETRY, ': M0 increased to ', M0
        ! Invalidate workspace cache (M0 changed)
        if (present(fw)) then
          if (fw%initialized) call feast_workspace_free(fw)
        end if
        cache_is_fresh = .false.
      end if

      ! Reallocate arrays if M0 changed
      if (allocated(E)) deallocate(E)
      if (allocated(X)) deallocate(X)
      if (allocated(res)) deallocate(res)
      allocate(E(M0))
      allocate(X(N, M0))
      allocate(res(M0))
      E = 0.0_dp
      X = cmplx(0.0_dp, 0.0_dp, kind=dp)
      res = 0.0_dp

      ! Extract upper triangle only (FEAST UPLO='U' requires col >= row).
      fw_match = .false.
      if (present(fw)) fw_match = feast_workspace_matches_pattern(H_csr, fw, N, M0)
      if (fw_match) then
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
        cache_is_fresh = .true.
      else
        ! Free stale cache before rebuilding
        if (present(fw)) then
          if (fw%initialized .and. .not. cache_is_fresh) call feast_workspace_free(fw)
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
          cache_is_fresh = .true.
        else
          deallocate(rowptr_loc, colind_loc)
        end if

        deallocate(val_loc)
      end if

      ! Check convergence: exit loop unless info=3 (subspace too small)
      if (info /= 3) exit

      ! If M0 already equals N, no point retrying
      if (M0 >= N) exit
    end do

    result%iterations = loop
    result%converged = (info == 0) .or. (info == 2)

    ! Persist the subspace size that actually produced a converged solve, so
    ! the next call can seed from it and avoid re-paying the retry tax. The
    ! energy window is fixed across a k-sweep, so the needed subspace is
    ! ~constant; seeding avoids repeated info=3 retries without changing which
    ! eigenpairs FEAST returns.
    if (present(fw)) then
      if (result%converged) fw%last_successful_m0 = M0
    end if

    if (info < 0) then
      print *, 'FEAST error: info =', info
    else if (info == 2) then
      print *, 'FEAST warning: max iterations reached, results may be inaccurate.'
    else if (info == 3) then
      if (M0_initial /= M0) then
        print *, '  WARNING: FEAST subspace still too small after retries.'
        print *, '  M0 grew from ', M0_initial, ' to ', M0, '. Increase m0 or narrow energy window.'
      else
        print *, '  WARNING: FEAST subspace too small (info=3). Missing eigenvalues.'
        print *, '  Increase m0 or narrow the energy window.'
      end if
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
  ! Finalizer: automatically called when an eigensolver_result goes out of scope.
  ! Delegates to eigensolver_result_free so existing manual frees remain valid.
  ! ==================================================================
  subroutine eigensolver_result_finalize(er)
    type(eigensolver_result), intent(inout) :: er
    call eigensolver_result_free(er)
  end subroutine eigensolver_result_finalize

  ! ==================================================================
  ! Deallocate result arrays.
  ! ==================================================================
  subroutine eigensolver_result_free(result)
    type(eigensolver_result), intent(inout) :: result

    if (result%was_freed) return
    result%was_freed = .true.
    if (allocated(result%eigenvalues))  deallocate(result%eigenvalues)
    if (allocated(result%eigenvectors)) deallocate(result%eigenvectors)
    result%nev_found = 0
    result%converged = .false.
    result%iterations = 0
  end subroutine eigensolver_result_free

  ! ==================================================================
  ! Finalizer: automatically called when a feast_workspace goes out of scope.
  ! Delegates to feast_workspace_free so existing manual frees remain valid.
  ! ==================================================================
  subroutine feast_workspace_finalize(fw)
    type(feast_workspace), intent(inout) :: fw
    call feast_workspace_free(fw)
  end subroutine feast_workspace_finalize

  ! ==================================================================
  ! Deallocate feast_workspace cached arrays.
  ! ==================================================================
  subroutine feast_workspace_free(fw)
    type(feast_workspace), intent(inout) :: fw

    if (fw%was_freed) return
    fw%was_freed = .true.
    if (allocated(fw%rowptr_loc)) deallocate(fw%rowptr_loc)
    if (allocated(fw%colind_loc)) deallocate(fw%colind_loc)
    fw%nnz_upper = 0
    fw%M0 = 0
    fw%N = 0
    fw%initialized = .false.
    fw%last_successful_m0 = 0
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
    complex(kind=dp), intent(inout), contiguous :: evecs(:,:)

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

  ! ==================================================================
  ! Validate eigensolver config — rejects invalid mode/method combos.
  ! ==================================================================
  subroutine eigensolver_config_validate(config)
    type(eigensolver_config), intent(in) :: config
    character(len=256) :: msg

    if (trim(config%method) == 'FEAST' .and. config%mode == EIGEN_MODE_INDEX) then
      error stop 'eigensolver_config_validate: FEAST solver does not support INDEX mode.'
    end if
    if (config%mode == EIGEN_MODE_ENERGY) then
      if (config%emin >= config%emax) then
        write(msg, '(A,ES12.4,A,ES12.4)') &
          'emin (', config%emin, ') must be < emax (', config%emax, ')'
        error stop 'eigensolver_config_validate: ' // trim(msg)
      end if
    end if
    if (config%mode == EIGEN_MODE_INDEX) then
      if (config%il < 1 .or. config%iu < config%il) then
        error stop 'eigensolver_config_validate: invalid il/iu range for INDEX mode.'
      end if
    end if
  end subroutine eigensolver_config_validate

  ! ==================================================================
  ! Polymorphic dispatch implementations.
  ! ==================================================================

  ! ------------------------------------------------------------------
  ! FEAST solver: legacy solve -> solve_sparse
  ! ------------------------------------------------------------------
#ifdef USE_MKL_FEAST
  subroutine feast_solve_dispatch(self, H_csr, config, result)
    class(feast_solver_t), intent(inout) :: self
    type(csr_matrix), intent(in) :: H_csr
    type(eigensolver_config), intent(in) :: config
    type(eigensolver_result), intent(out) :: result

    call self%solve_sparse(H_csr, config, result)
  end subroutine feast_solve_dispatch

  ! ------------------------------------------------------------------
  ! FEAST solver: solve_sparse with mode dispatch.
  ! ------------------------------------------------------------------
  subroutine feast_solve_sparse_dispatch(self, H_csr, config, result)
    class(feast_solver_t), intent(inout) :: self
    type(csr_matrix), intent(in) :: H_csr
    type(eigensolver_config), intent(in) :: config
    type(eigensolver_result), intent(out) :: result
    type(eigensolver_config) :: cfg_local

    select case (config%mode)
    case (EIGEN_MODE_ENERGY)
      ! Native FEAST: energy window search
      call solve_feast(H_csr, config, result, fw=self%ws)
    case (EIGEN_MODE_FULL)
      ! Gershgorin-bounded window -> FEAST energy search
      cfg_local = config
      call auto_compute_energy_window(H_csr, cfg_local%emin, cfg_local%emax)
      call solve_feast(H_csr, cfg_local, result, fw=self%ws)
    case (EIGEN_MODE_INDEX)
      error stop 'feast_solve_sparse_dispatch: FEAST does not support INDEX mode.'
    case default
      error stop 'feast_solve_sparse_dispatch: unknown mode.'
    end select
  end subroutine feast_solve_sparse_dispatch

  ! ------------------------------------------------------------------
  ! FEAST solver: solve_dense converts dense -> CSR, then calls sparse.
  ! ------------------------------------------------------------------
  subroutine feast_solve_dense(self, H, config, result)
    class(feast_solver_t), intent(inout) :: self
    complex(kind=dp), contiguous, intent(in) :: H(:,:)
    type(eigensolver_config), intent(in) :: config
    type(eigensolver_result), intent(out) :: result
    type(csr_matrix) :: H_csr

    call dense_to_csr_work(H, H_csr)
    call self%solve_sparse(H_csr, config, result)
    call csr_free(H_csr)
  end subroutine feast_solve_dense

  subroutine feast_solver_finalize(self)
    type(feast_solver_t), intent(inout) :: self
    ! ws component auto-finalizes via feast_workspace_finalize
  end subroutine feast_solver_finalize
#endif

  ! ------------------------------------------------------------------
  ! Dense LAPACK solver: legacy solve -> solve_sparse
  ! ------------------------------------------------------------------
  subroutine dense_lapack_solve_dispatch(self, H_csr, config, result)
    class(dense_lapack_solver_t), intent(inout) :: self
    type(csr_matrix), intent(in) :: H_csr
    type(eigensolver_config), intent(in) :: config
    type(eigensolver_result), intent(out) :: result

    call self%solve_sparse(H_csr, config, result)
  end subroutine dense_lapack_solve_dispatch

  ! ------------------------------------------------------------------
  ! Dense LAPACK solver: solve_sparse converts CSR -> dense.
  ! Convenience path for CSR inputs; not a hot path (no current caller
  ! routes a large dense matrix through CSR). All large dense solves go
  ! via solve_dense directly; all large CSR solves use FEAST.
  ! ------------------------------------------------------------------
  subroutine dense_solve_sparse_dispatch(self, H_csr, config, result)
    class(dense_lapack_solver_t), intent(inout) :: self
    type(csr_matrix), intent(in) :: H_csr
    type(eigensolver_config), intent(in) :: config
    type(eigensolver_result), intent(out) :: result
    integer :: N
    complex(kind=dp), allocatable :: H_dense(:,:)

    N = H_csr%nrows
    if (N <= 0) then
      result%converged = .false.
      result%nev_found = 0
      return
    end if

    allocate(H_dense(N, N))
    call csr_to_dense_work(H_csr, H_dense, N)
    call self%solve_dense(H_dense, config, result)
    deallocate(H_dense)
  end subroutine dense_solve_sparse_dispatch

  ! ------------------------------------------------------------------
  ! Dense LAPACK solver: solve_dense with mode dispatch.
  ! ------------------------------------------------------------------
  subroutine dense_solve_dense_dispatch(self, H, config, result)
    class(dense_lapack_solver_t), intent(inout) :: self
    complex(kind=dp), contiguous, intent(in) :: H(:,:)
    type(eigensolver_config), intent(in) :: config
    type(eigensolver_result), intent(out) :: result

    integer :: N, lda, ldz, lwork, info, nb, il_local, iu_local
    complex(kind=dp), allocatable :: wq(:)   ! 1-element workspace-query scratch
    real(kind=dp) :: vl, vu, abstol

    N = size(H, 1)
    if (N <= 0) then
      result%converged = .false.
      result%nev_found = 0
      return
    end if

    lda = N
    ldz = N
    abstol = 0.0_dp
    nb = 0

    ! (Re)allocate N-dependent buffers only when N grows
    if (N > self%cached_n) then
      if (allocated(self%A_buf)) then
        deallocate(self%A_buf, self%Z_buf, self%W_buf, self%rwork, self%iwork, self%ifail)
      end if
      allocate(self%A_buf(N,N), self%Z_buf(N,N), self%W_buf(N))
      allocate(self%rwork(7*N), self%iwork(5*N), self%ifail(N))
      self%cached_n = N
    end if

    ! Copy H (zheev/zheevx destroys input)
    self%A_buf = H

    select case (config%mode)
    case (EIGEN_MODE_FULL)
      ! zheev: all eigenvalues. Workspace query -> reuse/grow work.
      allocate(wq(1))
      call zheev('V', 'U', N, self%A_buf, lda, self%W_buf, wq, -1, self%rwork, info)
      if (info /= 0) error stop 'dense_solve_dense_dispatch: zheev workspace query failed.'
      lwork = max(1, nint(real(wq(1))))
      deallocate(wq)
      if (.not. allocated(self%work) .or. size(self%work) < lwork) then
        if (allocated(self%work)) deallocate(self%work)
        allocate(self%work(lwork))
      end if
      call zheev('V', 'U', N, self%A_buf, lda, self%W_buf, self%work, lwork, self%rwork, info)
      if (info == 0) then
        nb = N
        self%Z_buf = self%A_buf   ! zheev returns eigenvectors in A
      end if

    case (EIGEN_MODE_INDEX)
      il_local = max(1, config%il)
      iu_local = min(N, config%iu)
      allocate(wq(1))
      call zheevx('V', 'I', 'U', N, self%A_buf, lda, vl, vu, &
                   il_local, iu_local, abstol, nb, self%W_buf, self%Z_buf, ldz, &
                   wq, -1, self%rwork, self%iwork, self%ifail, info)
      if (info /= 0) error stop 'dense_solve_dense_dispatch: zheevx(INDEX) workspace query failed.'
      lwork = max(1, nint(real(wq(1))))
      deallocate(wq)
      if (.not. allocated(self%work) .or. size(self%work) < lwork) then
        if (allocated(self%work)) deallocate(self%work)
        allocate(self%work(lwork))
      end if
      call zheevx('V', 'I', 'U', N, self%A_buf, lda, vl, vu, &
                   il_local, iu_local, abstol, nb, self%W_buf, self%Z_buf, ldz, &
                   self%work, lwork, self%rwork, self%iwork, self%ifail, info)

    case (EIGEN_MODE_ENERGY)
      vl = config%emin
      vu = config%emax
      allocate(wq(1))
      call zheevx('V', 'V', 'U', N, self%A_buf, lda, vl, vu, &
                   1, N, abstol, nb, self%W_buf, self%Z_buf, ldz, &
                   wq, -1, self%rwork, self%iwork, self%ifail, info)
      if (info /= 0) error stop 'dense_solve_dense_dispatch: zheevx(ENERGY) workspace query failed.'
      lwork = max(1, nint(real(wq(1))))
      deallocate(wq)
      if (.not. allocated(self%work) .or. size(self%work) < lwork) then
        if (allocated(self%work)) deallocate(self%work)
        allocate(self%work(lwork))
      end if
      call zheevx('V', 'V', 'U', N, self%A_buf, lda, vl, vu, &
                   1, N, abstol, nb, self%W_buf, self%Z_buf, ldz, &
                   self%work, lwork, self%rwork, self%iwork, self%ifail, info)

    case default
      error stop 'dense_solve_dense_dispatch: unknown mode.'
    end select

    if (info == 0 .and. nb > 0) then
      result%converged = .true.
      result%nev_found = nb
      result%iterations = 1
      allocate(result%eigenvalues(nb))
      allocate(result%eigenvectors(N, nb))
      result%eigenvalues(1:nb) = self%W_buf(1:nb)
      result%eigenvectors(:, 1:nb) = self%Z_buf(:, 1:nb)
    else
      result%converged = .false.
      result%nev_found = 0
      if (info /= 0) print *, 'Dense eigensolver error: info =', info
    end if
  end subroutine dense_solve_dense_dispatch

  ! ==================================================================
  ! Finalizer: frees the cached LAPACK workspace buffers owned by the
  ! solver object. Called automatically when the solver goes out of
  ! scope; manual deallocate also triggers it.
  ! ==================================================================
  subroutine dense_lapack_solver_finalize(self)
    type(dense_lapack_solver_t), intent(inout) :: self
    if (allocated(self%A_buf))  deallocate(self%A_buf)
    if (allocated(self%Z_buf))  deallocate(self%Z_buf)
    if (allocated(self%work))   deallocate(self%work)
    if (allocated(self%W_buf))  deallocate(self%W_buf)
    if (allocated(self%rwork))  deallocate(self%rwork)
    if (allocated(self%iwork))  deallocate(self%iwork)
    if (allocated(self%ifail))  deallocate(self%ifail)
    self%cached_n = 0
  end subroutine dense_lapack_solver_finalize

  ! ==================================================================
  ! Internal: convert dense matrix to CSR.
  ! ==================================================================
  subroutine dense_to_csr_work(H, H_csr)
    complex(kind=dp), contiguous, intent(in) :: H(:,:)
    type(csr_matrix), intent(out) :: H_csr

    integer :: N, nnz, i, j
    integer, allocatable :: rows(:), cols(:)
    complex(kind=dp), allocatable :: vals(:)

    N = size(H, 1)
    nnz = 0
    do j = 1, N
      do i = 1, N
        if (abs(H(i, j)) > 0.0_dp) nnz = nnz + 1
      end do
    end do

    allocate(rows(nnz), cols(nnz), vals(nnz))
    nnz = 0
    do j = 1, N
      do i = 1, N
        if (abs(H(i, j)) > 0.0_dp) then
          nnz = nnz + 1
          rows(nnz) = i
          cols(nnz) = j
          vals(nnz) = H(i, j)
        end if
      end do
    end do

    call csr_build_from_coo(H_csr, N, N, nnz, rows, cols, vals)
    deallocate(rows, cols, vals)
  end subroutine dense_to_csr_work

  function make_eigensolver(config) result(solver)
    class(eigensolver_base), allocatable :: solver
    type(eigensolver_config), intent(in) :: config

    select case (trim(config%method))
    case ('DENSE')
      allocate(dense_lapack_solver_t :: solver)
#ifdef USE_MKL_FEAST
    case ('FEAST')
      allocate(feast_solver_t :: solver)
#else
    case ('FEAST')
      print *, 'WARNING: FEAST requested but not compiled; using dense LAPACK'
      allocate(dense_lapack_solver_t :: solver)
#endif
    case default
      error stop 'make_eigensolver: Unknown eigensolver method: ' // trim(config%method)
    end select
  end function make_eigensolver

end module eigensolver
