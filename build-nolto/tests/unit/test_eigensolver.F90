module test_eigensolver
  use funit
  use definitions
  use sparse_matrices
  use eigensolver
  implicit none

contains

  ! ==================================================================
  ! Helper: build a diagonal CSR matrix from real diagonal values.
  ! ==================================================================
  subroutine build_diagonal_csr(diag_vals, n, mat)
    integer, intent(in)             :: n
    real(kind=dp), intent(in)       :: diag_vals(n)
    type(csr_matrix), intent(out)   :: mat
    integer :: i
    integer, allocatable :: rows(:), cols(:)
    complex(kind=dp), allocatable :: vals(:)

    allocate(rows(n), cols(n), vals(n))
    do i = 1, n
      rows(i) = i
      cols(i) = i
      vals(i) = cmplx(diag_vals(i), 0.0_dp, kind=dp)
    end do
    call csr_build_from_coo(mat, n, n, n, rows, cols, vals)
    deallocate(rows, cols, vals)
  end subroutine build_diagonal_csr

  ! ==================================================================
  ! Helper: build a real symmetric tridiagonal CSR matrix.
  ! Main diagonal = dval, sub/super diagonal = oval.
  ! ==================================================================
  subroutine build_tridiag_csr(dval, oval, n, mat)
    real(kind=dp), intent(in)       :: dval, oval
    integer, intent(in)             :: n
    type(csr_matrix), intent(out)   :: mat
    integer :: i, nnz
    integer, allocatable :: rows(:), cols(:)
    complex(kind=dp), allocatable :: vals(:)

    nnz = n + 2*(n-1)
    allocate(rows(nnz), cols(nnz), vals(nnz))
    nnz = 0
    do i = 1, n
      nnz = nnz + 1
      rows(nnz) = i
      cols(nnz) = i
      vals(nnz) = cmplx(dval, 0.0_dp, kind=dp)
      if (i < n) then
        nnz = nnz + 1
        rows(nnz) = i
        cols(nnz) = i + 1
        vals(nnz) = cmplx(oval, 0.0_dp, kind=dp)
      end if
      if (i > 1) then
        nnz = nnz + 1
        rows(nnz) = i
        cols(nnz) = i - 1
        vals(nnz) = cmplx(oval, 0.0_dp, kind=dp)
      end if
    end do
    call csr_build_from_coo(mat, n, n, nnz, rows(1:nnz), &
                            cols(1:nnz), vals(1:nnz))
    deallocate(rows, cols, vals)
  end subroutine build_tridiag_csr

  ! ==================================================================
  ! Test 1: FEAST on a diagonal matrix.
  ! Eigenvalues should equal the diagonal entries within the interval.
  ! ==================================================================
  !@test
  subroutine test_feast_diagonal()
    integer, parameter :: n = 20
    real(kind=dp) :: diag(n)
    type(csr_matrix) :: H
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: res
    integer :: i

    do i = 1, n
      diag(i) = real(i, kind=dp)
    end do
    call build_diagonal_csr(diag, n, H)

    cfg%method = 'FEAST'
    cfg%emin = 0.5_dp
    cfg%emax = 5.5_dp
    cfg%nev = 5

    call solve_feast(H, cfg, res)

#line 94 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 94) )
  if (anyExceptions()) return
#line 95 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 95 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(5, res%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 95) )
  if (anyExceptions()) return
#line 96 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    ! Eigenvalues in [0.5, 5.5] are 1, 2, 3, 4, 5
    do i = 1, res%nev_found
#line 99 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(real(i, kind=dp), res%eigenvalues(i), tolerance=1.0e-8_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 99) )
  if (anyExceptions()) return
#line 100 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
    end do

    call eigensolver_result_free(res)
    call csr_free(H)
  end subroutine test_feast_diagonal

  ! ==================================================================
  ! Test 2: FEAST on a tridiagonal Laplacian (1D, hard-wall).
  !
  ! The 1D discrete Laplacian on N points with spacing h has eigenvalues:
  !   lambda_k = 2*(1 - cos(k*pi/(N+1))) / h^2,  k = 1..N
  ! ==================================================================
  !@test
  subroutine test_feast_tridiagonal_laplacian()
    integer, parameter :: n = 20
    real(kind=dp), parameter :: grid_dx = 1.0_dp
    real(kind=dp), parameter :: pi = pi_dp  ! use module constant
    type(csr_matrix) :: Hmat
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: res
    real(kind=dp) :: expected
    integer :: k

    call build_tridiag_csr(2.0_dp/grid_dx**2, -1.0_dp/grid_dx**2, n, Hmat)

    ! First verify with dense LAPACK that the matrix is correct
    cfg%method = 'ARPACK'
    cfg%nev = 3
    cfg%emin = 0.0_dp
    cfg%emax = 0.0_dp
    call solve_dense_lapack(Hmat, cfg, res)

#line 132 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 132) )
  if (anyExceptions()) return
#line 133 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 133 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(3, res%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 133) )
  if (anyExceptions()) return
#line 134 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    ! Check first 3 eigenvalues match analytical Laplacian values
    do k = 1, 3
      expected = 2.0_dp * (1.0_dp - cos(real(k, kind=dp) * pi / real(n+1, kind=dp))) / grid_dx**2
#line 138 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(expected, res%eigenvalues(k), tolerance=1.0e-6_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 138) )
  if (anyExceptions()) return
#line 139 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
    end do

    call eigensolver_result_free(res)

    ! Now test FEAST with a wider interval to verify convergence.
    ! Use emin slightly below the smallest eigenvalue and emax above the
    ! first few eigenvalues, with sufficient M0.
    cfg%method = 'FEAST'
    cfg%emin = 0.01_dp
    cfg%emax = 0.25_dp
    cfg%nev = 3
    cfg%feast_m0 = 10

    call solve_feast(Hmat, cfg, res)

#line 154 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 154) )
  if (anyExceptions()) return
#line 155 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 155 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%nev_found >= 2, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 155) )
  if (anyExceptions()) return
#line 156 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    ! Check that FEAST eigenvalues match analytical values (sorted)
    do k = 1, min(res%nev_found, 3)
      expected = 2.0_dp * (1.0_dp - cos(real(k, kind=dp) * pi / real(n+1, kind=dp))) / grid_dx**2
#line 160 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(expected, res%eigenvalues(k), tolerance=1.0e-6_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 160) )
  if (anyExceptions()) return
#line 161 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
    end do

    call eigensolver_result_free(res)
    call csr_free(Hmat)
  end subroutine test_feast_tridiagonal_laplacian

  ! ==================================================================
  ! Test 3: Eigenvector orthonormality (v^H v = I).
  ! ==================================================================
  !@test
  subroutine test_eigenvector_orthonormality()
    integer, parameter :: n = 20
    real(kind=dp) :: diag(n)
    type(csr_matrix) :: H
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: res
    complex(kind=dp) :: dot_prod
    real(kind=dp) :: ortho_err
    integer :: i, j

    do i = 1, n
      diag(i) = real(i, kind=dp)
    end do
    call build_diagonal_csr(diag, n, H)

    cfg%method = 'FEAST'
    cfg%emin = 0.5_dp
    cfg%emax = 5.5_dp
    cfg%nev = 5

    call solve_feast(H, cfg, res)

#line 193 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 193) )
  if (anyExceptions()) return
#line 194 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 194 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%nev_found >= 2, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 194) )
  if (anyExceptions()) return
#line 195 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    do i = 1, res%nev_found
      do j = 1, res%nev_found
        dot_prod = dot_product(res%eigenvectors(:,i), res%eigenvectors(:,j))
        if (i == j) then
          ortho_err = abs(dot_prod - cmplx(1.0_dp, 0.0_dp, kind=dp))
        else
          ortho_err = abs(dot_prod)
        end if
#line 204 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(ortho_err < 1.0e-8_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 204) )
  if (anyExceptions()) return
#line 205 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
      end do
    end do

    call eigensolver_result_free(res)
    call csr_free(H)
  end subroutine test_eigenvector_orthonormality

  ! ==================================================================
  ! Test 4: Identity matrix -- all eigenvalues = 1.
  ! ==================================================================
  !@test
  subroutine test_feast_identity()
    integer, parameter :: n = 10
    real(kind=dp) :: diag(n)
    type(csr_matrix) :: H
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: res
    integer :: i

    do i = 1, n
      diag(i) = 1.0_dp
    end do
    call build_diagonal_csr(diag, n, H)

    cfg%method = 'FEAST'
    cfg%emin = 0.5_dp
    cfg%emax = 1.5_dp
    cfg%nev = 10

    call solve_feast(H, cfg, res)

#line 236 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 236) )
  if (anyExceptions()) return
#line 237 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 237 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(10, res%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 237) )
  if (anyExceptions()) return
#line 238 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    do i = 1, res%nev_found
#line 240 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(1.0_dp, res%eigenvalues(i), tolerance=1.0e-8_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 240) )
  if (anyExceptions()) return
#line 241 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
    end do

    call eigensolver_result_free(res)
    call csr_free(H)
  end subroutine test_feast_identity

  ! ==================================================================
  ! Test 5: Dispatch routine selects FEAST correctly.
  ! ==================================================================
  !@test
  subroutine test_dispatch_feast()
    integer, parameter :: n = 5
    real(kind=dp) :: diag(n)
    type(csr_matrix) :: H
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: res
    integer :: i

    do i = 1, n
      diag(i) = real(i, kind=dp)
    end do
    call build_diagonal_csr(diag, n, H)

    cfg%method = 'FEAST'
    cfg%emin = 0.5_dp
    cfg%emax = 3.5_dp
    cfg%nev = 3

    call solve_sparse_evp(H, cfg, res)

#line 271 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 271) )
  if (anyExceptions()) return
#line 272 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 272 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(3, res%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 272) )
  if (anyExceptions()) return
#line 273 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    call eigensolver_result_free(res)
    call csr_free(H)
  end subroutine test_dispatch_feast

  ! ==================================================================
  ! Test 6: ARPACK fallback uses dense zheevx.
  ! ==================================================================
  !@test
  subroutine test_arpack_fallback()
    integer, parameter :: n = 5
    real(kind=dp) :: diag(n)
    type(csr_matrix) :: H
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: res
    integer :: i

    do i = 1, n
      diag(i) = real(i, kind=dp)
    end do
    call build_diagonal_csr(diag, n, H)

    cfg%method = 'ARPACK'
    cfg%nev = 3
    cfg%emin = 0.0_dp
    cfg%emax = 0.0_dp

    call solve_dense_lapack(H, cfg, res)

#line 302 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 302) )
  if (anyExceptions()) return
#line 303 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 303 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(3, res%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 303) )
  if (anyExceptions()) return
#line 304 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    do i = 1, res%nev_found
#line 306 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(real(i, kind=dp), res%eigenvalues(i), tolerance=1.0e-8_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 306) )
  if (anyExceptions()) return
#line 307 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
    end do

    call eigensolver_result_free(res)
    call csr_free(H)
  end subroutine test_arpack_fallback

  ! ==================================================================
  ! Test 7: Unknown method returns unconverged.
  ! ==================================================================
  !@test
  subroutine test_unknown_method()
    integer, parameter :: n = 4
    real(kind=dp) :: diag(n)
    type(csr_matrix) :: H
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: res
    integer :: i

    do i = 1, n
      diag(i) = real(i, kind=dp)
    end do
    call build_diagonal_csr(diag, n, H)

    cfg%method = 'UNKNOWN'
    cfg%nev = 4

    call solve_sparse_evp(H, cfg, res)

#line 335 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertFalse(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 335) )
  if (anyExceptions()) return
#line 336 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 336 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(0, res%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 336) )
  if (anyExceptions()) return
#line 337 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    call csr_free(H)
  end subroutine test_unknown_method

  ! ==================================================================
  ! Test 8: Empty matrix returns gracefully.
  ! ==================================================================
  !@test
  subroutine test_feast_empty_matrix()
    type(csr_matrix) :: H
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: res

    call csr_init(H, 0, 0)

    cfg%method = 'FEAST'
    cfg%emin = -1.0_dp
    cfg%emax = 1.0_dp
    cfg%nev = 4

    call solve_feast(H, cfg, res)

#line 359 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertFalse(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 359) )
  if (anyExceptions()) return
#line 360 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 360 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(0, res%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 360) )
  if (anyExceptions()) return
#line 361 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    call csr_free(H)
  end subroutine test_feast_empty_matrix

  ! ==================================================================
  ! Test 9: Gershgorin energy window estimation.
  ! ==================================================================
  !@test
  subroutine test_energy_window_gershgorin()
    integer, parameter :: n = 10
    real(kind=dp) :: diag(n)
    type(csr_matrix) :: H
    real(kind=dp) :: emin, emax
    integer :: i

    do i = 1, n
      diag(i) = real(i, kind=dp)  ! eigenvalues 1..10
    end do
    call build_diagonal_csr(diag, n, H)

    call auto_compute_energy_window(H, emin, emax)

    ! For diagonal matrix, Gershgorin bounds = exact diagonal values
    ! with 10% margin
#line 385 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(emin < 1.0_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 385) )
  if (anyExceptions()) return
#line 386 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 386 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(emax > 10.0_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 386) )
  if (anyExceptions()) return
#line 387 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    call csr_free(H)
  end subroutine test_energy_window_gershgorin

  ! ==================================================================
  ! Test 10: feast_workspace caches upper-triangle structure.
  ! ==================================================================
  !@test
  subroutine test_feast_workspace_cached_matches_uncached()
    ! Solve the same eigenvalue problem with and without workspace.
    ! Eigenvalues must be identical.
    integer, parameter :: n = 20
    real(kind=dp) :: diag(n)
    type(csr_matrix) :: H
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: res1, res2
    type(feast_workspace) :: fw
    integer :: i, min_nev

    do i = 1, n
      diag(i) = real(i, kind=dp)
    end do
    call build_diagonal_csr(diag, n, H)

    cfg%method = 'FEAST'
    cfg%emin = 0.5_dp
    cfg%emax = 5.5_dp
    cfg%nev = 5

    ! Solve without workspace
    call solve_feast(H, cfg, res1)

    ! Solve with workspace (first call = slow path, initializes cache)
    call solve_feast(H, cfg, res2, fw)
#line 421 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(fw%initialized, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 421) )
  if (anyExceptions()) return
#line 422 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 422 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(res1%nev_found, res2%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 422) )
  if (anyExceptions()) return
#line 423 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
    min_nev = min(res1%nev_found, res2%nev_found)
    do i = 1, min_nev
#line 425 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(res1%eigenvalues(i), res2%eigenvalues(i), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 425) )
  if (anyExceptions()) return
#line 426 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
    end do
    call eigensolver_result_free(res1)
    call eigensolver_result_free(res2)

    ! Solve again with cached workspace (tests fast path)
    call solve_feast(H, cfg, res2, fw)
#line 432 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(fw%initialized, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 432) )
  if (anyExceptions()) return
#line 433 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
    call solve_feast(H, cfg, res1)
#line 434 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(res1%nev_found, res2%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 434) )
  if (anyExceptions()) return
#line 435 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
    min_nev = min(res1%nev_found, res2%nev_found)
    do i = 1, min_nev
#line 437 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(res1%eigenvalues(i), res2%eigenvalues(i), tolerance=1.0e-12_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 437) )
  if (anyExceptions()) return
#line 438 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
    end do

    call eigensolver_result_free(res1)
    call eigensolver_result_free(res2)
    call feast_workspace_free(fw)
#line 443 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertFalse(fw%initialized, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 443) )
  if (anyExceptions()) return
#line 444 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
    call csr_free(H)
  end subroutine test_feast_workspace_cached_matches_uncached

  ! ==================================================================
  ! Test 11: FEAST workspace rebuilds on sparsity pattern change.
  !
  ! Initializes cache with diagonal upper-triangle structure (N nonzeros),
  ! then reuses workspace on same-size tridiagonal matrix with different
  ! sparsity (2N-1 upper-triangle nonzeros).  The fast path at line 148 of
  ! eigensolver.f90 only checks dimensions (N, M0) but NOT the sparsity
  ! pattern (rowptr_loc, colind_loc).  This test catches that bug: cached
  ! colind_loc has wrong entries and wrong count, producing wrong
  ! eigenvalues or a crash.
  ! ==================================================================
  !@test
  subroutine test_feast_workspace_rebuilds_on_pattern_change()
    integer, parameter :: n = 20
    real(kind=dp) :: diag(n)
    type(csr_matrix) :: H_diag, H_tri
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: cached, uncached
    type(feast_workspace) :: fw
    integer :: i, min_nev

    do i = 1, n
      diag(i) = real(i, kind=dp)
    end do
    call build_diagonal_csr(diag, n, H_diag)

    ! Tridiagonal with diagonal=3.0, off-diagonal=-0.1.  Eigenvalues cluster
    ! around 3.0 +/- 0.2*cos(k*pi/(n+1)) ~ 2.8..3.2.  Use a narrow window
    ! [2.82, 2.85] so only 1-2 eigenvalues fall inside, avoiding FEAST info=3
    ! (subspace too small) when all 20 eigenvalues sit in a wider window.
    call build_tridiag_csr(3.0_dp, -0.1_dp, n, H_tri)

    cfg%method = 'FEAST'
    cfg%emin = 2.82_dp
    cfg%emax = 2.85_dp
    cfg%nev = 3
    cfg%feast_m0 = 8

    ! Initialize cache with diagonal upper-triangle structure.
    call solve_feast(H_diag, cfg, cached, fw)
    call eigensolver_result_free(cached)
#line 488 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(fw%initialized, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 488) )
  if (anyExceptions()) return
#line 489 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    ! Reuse same workspace on same-size but different sparsity matrix.
    call solve_feast(H_tri, cfg, cached, fw)
    call solve_feast(H_tri, cfg, uncached)

#line 494 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(uncached%nev_found, cached%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 494) )
  if (anyExceptions()) return
#line 495 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
    min_nev = min(uncached%nev_found, cached%nev_found)
    do i = 1, min_nev
#line 497 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(uncached%eigenvalues(i), cached%eigenvalues(i), tolerance=1.0e-10_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 497) )
  if (anyExceptions()) return
#line 498 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
    end do

    call eigensolver_result_free(cached)
    call eigensolver_result_free(uncached)
    call feast_workspace_free(fw)
    call csr_free(H_diag)
    call csr_free(H_tri)
  end subroutine test_feast_workspace_rebuilds_on_pattern_change

end module test_eigensolver

module Wraptest_eigensolver
   use FUnit
   use test_eigensolver
   implicit none
   private

contains


end module Wraptest_eigensolver

function test_eigensolver_suite() result(suite)
   use FUnit
   use test_eigensolver
   use Wraptest_eigensolver
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   suite = TestSuite('test_eigensolver_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_feast_diagonal', &
      test_feast_diagonal))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_feast_tridiagonal_laplacian', &
      test_feast_tridiagonal_laplacian))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_eigenvector_orthonormality', &
      test_eigenvector_orthonormality))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_feast_identity', &
      test_feast_identity))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_dispatch_feast', &
      test_dispatch_feast))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_arpack_fallback', &
      test_arpack_fallback))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_unknown_method', &
      test_unknown_method))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_feast_empty_matrix', &
      test_feast_empty_matrix))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_energy_window_gershgorin', &
      test_energy_window_gershgorin))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_feast_workspace_cached_matches_uncached', &
      test_feast_workspace_cached_matches_uncached))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_feast_workspace_rebuilds_on_pattern_change', &
      test_feast_workspace_rebuilds_on_pattern_change))
   call suite%addTest(t)


end function test_eigensolver_suite

