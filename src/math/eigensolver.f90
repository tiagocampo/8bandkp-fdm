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
  public :: solve_sparse_evp, solve_feast, solve_dense_lapack
  public :: auto_compute_energy_window, eigensolver_result_free

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

contains

  ! ==================================================================
  ! Main dispatch: select FEAST or ARPACK based on config.
  ! ==================================================================
  subroutine solve_sparse_evp(H_csr, config, result)
    type(csr_matrix), intent(in)          :: H_csr
    type(eigensolver_config), intent(in)  :: config
    type(eigensolver_result), intent(out) :: result

    select case (trim(config%method))
#ifdef USE_MKL_FEAST
    case ('FEAST')
      call solve_feast(H_csr, config, result)
    case ('ARPACK')
      call solve_dense_lapack(H_csr, config, result)
#else
    case ('FEAST', 'ARPACK')
      ! FEAST unavailable: fall back to dense LAPACK
      print *, 'WARNING: FEAST/ARPACK not available. Falling back to dense LAPACK.'
      print *, '  Energy window [emin, emax] will be IGNORED.'
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
  subroutine solve_feast(H_csr, config, result)
    type(csr_matrix), intent(in)          :: H_csr
    type(eigensolver_config), intent(in)  :: config
    type(eigensolver_result), intent(out) :: result

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
    ! First pass: count upper-triangle entries per row.
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

    ! Call FEAST directly (explicit interface from linalg module)
    call zfeast_hcsrev('U', N, val_loc, rowptr_loc, colind_loc, &
                       fpm, epsout, loop, config%emin, config%emax, M0, &
                       E, X, M, res, info)

    deallocate(val_loc, rowptr_loc, colind_loc)

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
    else
      result%nev_found = 0
    end if

    deallocate(E, X, res)
  end subroutine solve_feast
#endif /* USE_MKL_FEAST */

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

    ! Allocate output arrays
    allocate(W(N))
    allocate(Z(N, nev_want))
    allocate(rwork(7 * N))
    allocate(iwork(5 * N))
    allocate(ifail(N))

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

    emin = rmin - max(0.1_dp * abs(rmin), 0.5_dp)
    emax = rmax + max(0.1_dp * abs(rmax), 0.5_dp)
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

end module eigensolver
