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
    real(kind=dp), parameter :: pi = 3.14159265358979323846264338327950288_dp
    type(csr_matrix) :: Hmat
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: res
    real(kind=dp) :: expected
    integer :: k

    call build_tridiag_csr(2.0_dp/grid_dx**2, -1.0_dp/grid_dx**2, n, Hmat)

    ! First verify with ARPACK dense solver that the matrix is correct
    cfg%method = 'ARPACK'
    cfg%nev = 3
    call solve_arpack(Hmat, cfg, res)

#line 130 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 130) )
  if (anyExceptions()) return
#line 131 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 131 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(3, res%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 131) )
  if (anyExceptions()) return
#line 132 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    ! Check first 3 eigenvalues match analytical Laplacian values
    do k = 1, 3
      expected = 2.0_dp * (1.0_dp - cos(real(k, kind=dp) * pi / real(n+1, kind=dp))) / grid_dx**2
#line 136 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(expected, res%eigenvalues(k), tolerance=1.0e-6_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 136) )
  if (anyExceptions()) return
#line 137 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
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

#line 152 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 152) )
  if (anyExceptions()) return
#line 153 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 153 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%nev_found >= 2, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    ! Check that FEAST eigenvalues match analytical values (sorted)
    do k = 1, min(res%nev_found, 3)
      expected = 2.0_dp * (1.0_dp - cos(real(k, kind=dp) * pi / real(n+1, kind=dp))) / grid_dx**2
#line 158 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(expected, res%eigenvalues(k), tolerance=1.0e-6_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 158) )
  if (anyExceptions()) return
#line 159 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
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

#line 191 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 191) )
  if (anyExceptions()) return
#line 192 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 192 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%nev_found >= 2, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 192) )
  if (anyExceptions()) return
#line 193 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    do i = 1, res%nev_found
      do j = 1, res%nev_found
        dot_prod = dot_product(res%eigenvectors(:,i), res%eigenvectors(:,j))
        if (i == j) then
          ortho_err = abs(dot_prod - cmplx(1.0_dp, 0.0_dp, kind=dp))
        else
          ortho_err = abs(dot_prod)
        end if
#line 202 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(ortho_err < 1.0e-8_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 202) )
  if (anyExceptions()) return
#line 203 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
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

#line 234 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 234) )
  if (anyExceptions()) return
#line 235 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 235 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(10, res%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 235) )
  if (anyExceptions()) return
#line 236 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    do i = 1, res%nev_found
#line 238 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(1.0_dp, res%eigenvalues(i), tolerance=1.0e-8_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 238) )
  if (anyExceptions()) return
#line 239 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
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

#line 269 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 269) )
  if (anyExceptions()) return
#line 270 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 270 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(3, res%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 270) )
  if (anyExceptions()) return
#line 271 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

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

    call solve_arpack(H, cfg, res)

#line 298 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 298) )
  if (anyExceptions()) return
#line 299 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 299 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(3, res%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 299) )
  if (anyExceptions()) return
#line 300 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    do i = 1, res%nev_found
#line 302 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(real(i, kind=dp), res%eigenvalues(i), tolerance=1.0e-8_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 302) )
  if (anyExceptions()) return
#line 303 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
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

#line 331 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertFalse(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 331) )
  if (anyExceptions()) return
#line 332 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 332 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(0, res%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 332) )
  if (anyExceptions()) return
#line 333 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

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

#line 355 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertFalse(res%converged, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 355) )
  if (anyExceptions()) return
#line 356 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 356 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertEqual(0, res%nev_found, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 356) )
  if (anyExceptions()) return
#line 357 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

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
#line 381 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(emin < 1.0_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 381) )
  if (anyExceptions()) return
#line 382 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
#line 382 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"
  call assertTrue(emax > 10.0_dp, &
 & location=SourceLocation( &
 & 'test_eigensolver.pf', &
 & 382) )
  if (anyExceptions()) return
#line 383 "/data/8bandkp-fdm/tests/unit/test_eigensolver.pf"

    call csr_free(H)
  end subroutine test_energy_window_gershgorin

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


end function test_eigensolver_suite

