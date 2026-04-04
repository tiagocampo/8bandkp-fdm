module test_hamiltonian_2d
  use funit
  use definitions
  use parameters
  use sparse_matrices
  use geometry
  use finitedifferences
  use hamiltonianConstructor
  implicit none

contains

  ! ==================================================================
  ! Kronecker product identity test
  ! ==================================================================

  !@test
  subroutine test_kron_identity_size()
    ! kron(I_n, I_n) should give N^2 x N^2 identity
    integer, parameter :: n = 4
    type(csr_matrix) :: C
    complex(kind=dp) :: eye(n, n)
    integer :: i

    eye = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do i = 1, n
      eye(i, i) = cmplx(1.0_dp, 0.0_dp, kind=dp)
    end do

    call kron_dense_dense(eye, n, n, eye, n, n, C)

#line 32 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(n * n, C%nrows, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 32) )
  if (anyExceptions()) return
#line 33 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 33 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(n * n, C%ncols, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 33) )
  if (anyExceptions()) return
#line 34 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 34 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(n * n, C%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 34) )
  if (anyExceptions()) return
#line 35 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    call csr_free(C)
  end subroutine test_kron_identity_size

  ! ==================================================================
  ! 2D Laplacian sparsity and symmetry test
  ! ==================================================================

  !@test
  subroutine test_laplacian_2d_sparsity()
    ! Verify the 2D Laplacian via Kronecker products has correct
    ! dimensions and reasonable sparsity for a uniform grid.
    ! Note: boundary rows use one-sided stencils which break exact symmetry,
    ! so we check interior symmetry only.
    integer, parameter :: nx = 5, ny = 5
    real(kind=dp), parameter :: dx = 1.0_dp, dy = 1.0_dp
    integer :: ngrid
    real(kind=dp), allocatable :: D2x(:,:), D2y(:,:), Ix_mat(:,:), Iy_mat(:,:)
    complex(kind=dp), allocatable :: cD2x(:,:), cD2y(:,:), cIx(:,:), cIy(:,:)
    type(csr_matrix) :: kron_a, kron_b, lap

    ngrid = nx * ny

    call buildFD2ndDerivMatrix(nx, dx, 2, D2x)
    call buildFD2ndDerivMatrix(ny, dy, 2, D2y)
    call Identity(nx, Ix_mat)
    call Identity(ny, Iy_mat)

    allocate(cD2x(nx,nx)); cD2x = cmplx(D2x, 0.0_dp, kind=dp)
    allocate(cD2y(ny,ny)); cD2y = cmplx(D2y, 0.0_dp, kind=dp)
    allocate(cIx(nx,nx));   cIx = cmplx(Ix_mat, 0.0_dp, kind=dp)
    allocate(cIy(ny,ny));   cIy = cmplx(Iy_mat, 0.0_dp, kind=dp)

    ! Laplacian = D2x x Iy + Ix x D2y
    call kron_dense_dense(cD2x, nx, nx, cIy, ny, ny, kron_a)
    call kron_eye_dense(nx, cD2y, ny, ny, kron_b)
    call csr_add(kron_a, kron_b, lap)

#line 72 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, lap%nrows, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 72) )
  if (anyExceptions()) return
#line 73 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 73 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, lap%ncols, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 73) )
  if (anyExceptions()) return
#line 74 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 74 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(lap%nnz > 0, message="Laplacian has nonzeros", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 74) )
  if (anyExceptions()) return
#line 75 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 75 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(lap%nnz <= 2 * 3 * ngrid, message="Laplacian sparsity bounded", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 75) )
  if (anyExceptions()) return
#line 76 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Interior rows of 2D Laplacian should be symmetric
    call assert_csr_interior_symmetric(lap, nx, ny, "2D Laplacian interior")

    call csr_free(kron_a)
    call csr_free(kron_b)
    call csr_free(lap)
    deallocate(D2x, D2y, Ix_mat, Iy_mat, cD2x, cD2y, cIx, cIy)
  end subroutine test_laplacian_2d_sparsity

  ! ==================================================================
  ! Cross-derivative sparse pattern test
  ! ==================================================================

  !@test
  subroutine test_cross_derivative_pattern()
    ! d^2/dxdy = D1x x D1y.  Verify correct dimensions and nonzeros.
    integer, parameter :: nx = 4, ny = 5
    real(kind=dp), parameter :: dx = 1.0_dp, dy = 2.0_dp
    real(kind=dp), allocatable :: D1x(:,:), D1y(:,:)
    complex(kind=dp), allocatable :: cD1x(:,:), cD1y(:,:)
    type(csr_matrix) :: cross

    call buildFD1stDerivMatrix(nx, dx, 2, D1x)
    call buildFD1stDerivMatrix(ny, dy, 2, D1y)

    allocate(cD1x(nx,nx)); cD1x = cmplx(D1x, 0.0_dp, kind=dp)
    allocate(cD1y(ny,ny)); cD1y = cmplx(D1y, 0.0_dp, kind=dp)

    call kron_dense_dense(cD1x, nx, nx, cD1y, ny, ny, cross)

#line 107 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(nx * ny, cross%nrows, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 107) )
  if (anyExceptions()) return
#line 108 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 108 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(nx * ny, cross%ncols, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 108) )
  if (anyExceptions()) return
#line 109 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 109 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(cross%nnz > 0, message="Cross derivative has nonzeros", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 109) )
  if (anyExceptions()) return
#line 110 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    call csr_free(cross)
    deallocate(D1x, D1y, cD1x, cD1y)
  end subroutine test_cross_derivative_pattern

  ! ==================================================================
  ! kpterms_2d sparsity pattern test
  ! ==================================================================

  !@test
  subroutine test_kpterms_2d_sparsity()
    ! Build kpterms_2d for a simple GaAs rectangle and check sparsity.
    integer, parameter :: nx = 4, ny = 4
    real(kind=dp), parameter :: dx = 5.0_dp, dy = 5.0_dp
    integer :: ngrid, ij
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    type(region_spec), allocatable :: regions(:)
    character(len=255), allocatable :: material_names(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    integer :: mat2d(nx, ny)

    ngrid = nx * ny
    mat2d = 1  ! single material

    call grid_init_rect(grid, nx, ny, dx, dy, mat2d)

    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)

    allocate(regions(1))
    regions(1)%material = "GaAs"

    call confinementInitialization_2d(grid, params, regions, profile_2d, kpterms_2d, FDorder=2)

    ! Check array size
#line 149 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(11, size(kpterms_2d), &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 149) )
  if (anyExceptions()) return
#line 150 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Check all matrices have correct dimensions
    do ij = 1, 11
#line 153 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, kpterms_2d(ij)%nrows, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 154 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, kpterms_2d(ij)%ncols, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 154) )
  if (anyExceptions()) return
#line 155 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    end do

    ! Diagonal terms (1-4, 10) should have exactly ngrid nonzeros
#line 158 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, kpterms_2d(1)%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 158) )
  if (anyExceptions()) return
#line 159 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 159 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, kpterms_2d(2)%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 159) )
  if (anyExceptions()) return
#line 160 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 160 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, kpterms_2d(3)%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 160) )
  if (anyExceptions()) return
#line 161 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 161 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, kpterms_2d(4)%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 161) )
  if (anyExceptions()) return
#line 162 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 162 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, kpterms_2d(10)%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 162) )
  if (anyExceptions()) return
#line 163 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Laplacian-based terms (5, 7, 8) should have NNZ > ngrid
#line 165 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(kpterms_2d(5)%nnz > ngrid, message="Term 5 has off-diagonal", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 165) )
  if (anyExceptions()) return
#line 166 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 166 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(kpterms_2d(7)%nnz > ngrid, message="Term 7 has off-diagonal", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 166) )
  if (anyExceptions()) return
#line 167 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 167 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(kpterms_2d(8)%nnz > ngrid, message="Term 8 has off-diagonal", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 167) )
  if (anyExceptions()) return
#line 168 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Cross-derivative (11) should have nonzeros
#line 170 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(kpterms_2d(11)%nnz > 0, message="Term 11 has entries", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 170) )
  if (anyExceptions()) return
#line 171 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Check profile dimensions
#line 173 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, size(profile_2d, 1), &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 173) )
  if (anyExceptions()) return
#line 174 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 174 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(3, size(profile_2d, 2), &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 174) )
  if (anyExceptions()) return
#line 175 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Profile should be non-zero for active cells
#line 177 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(profile_2d(1, 1) /= 0.0_dp, message="EV non-zero", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 177) )
  if (anyExceptions()) return
#line 178 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 178 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(profile_2d(1, 3) /= 0.0_dp, message="EC non-zero", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 178) )
  if (anyExceptions()) return
#line 179 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Cleanup
    do ij = 1, 11
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_kpterms_2d_sparsity

  ! ==================================================================
  ! Hermiticity test (interior only)
  ! ==================================================================

  !@test
  subroutine test_kpterms_2d_hermitian()
    ! For a uniform material profile, Laplacian-based terms should be Hermitian
    ! (symmetric + real) in the interior.
    integer, parameter :: nnx = 5, nny = 5
    real(kind=dp), parameter :: dx = 1.0_dp, dy = 1.0_dp
    integer :: ij
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    type(region_spec), allocatable :: regions(:)
    character(len=255), allocatable :: material_names(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    integer :: mat2d(nnx, nny)

    mat2d = 1  ! single material

    call grid_init_rect(grid, nnx, nny, dx, dy, mat2d)

    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)

    allocate(regions(1))
    regions(1)%material = "GaAs"

    call confinementInitialization_2d(grid, params, regions, profile_2d, kpterms_2d, FDorder=2)

    ! Term 5: A * Laplacian -- should be Hermitian for uniform profile
    ! (interior rows only; boundary rows use one-sided stencils)
    call assert_csr_interior_hermitian(kpterms_2d(5), nnx, nny, "Term 5 (A*Laplacian)")

    ! Term 7: (gamma1-2gamma2) * Laplacian
    call assert_csr_interior_hermitian(kpterms_2d(7), nnx, nny, "Term 7 (Q*Laplacian)")

    ! Term 8: (gamma1+2gamma2) * Laplacian
    call assert_csr_interior_hermitian(kpterms_2d(8), nnx, nny, "Term 8 (T*Laplacian)")

    do ij = 1, 11
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_kpterms_2d_hermitian

  ! ==================================================================
  ! Profile construction test
  ! ==================================================================

  !@test
  subroutine test_profile_2d_band_edges()
    ! Verify profile_2d contains correct band-edge values for GaAs.
    integer, parameter :: nnx = 3, nny = 3
    real(kind=dp), parameter :: dx = 5.0_dp, dy = 5.0_dp
    integer :: ngrid, ij
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    type(region_spec), allocatable :: regions(:)
    character(len=255), allocatable :: material_names(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    integer :: mat2d(nnx, nny)
    real(kind=dp) :: expected_EV, expected_EC, expected_DeltaSO

    ngrid = nnx * nny
    mat2d = 1

    call grid_init_rect(grid, nnx, nny, dx, dy, mat2d)

    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)

    expected_EV = params(1)%EV
    expected_EC = params(1)%EC
    expected_DeltaSO = params(1)%DeltaSO

    allocate(regions(1))
    regions(1)%material = "GaAs"

    call confinementInitialization_2d(grid, params, regions, profile_2d, kpterms_2d, FDorder=2)

    ! All cells should have the same GaAs band edges
    do ij = 1, ngrid
#line 278 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(expected_EV, profile_2d(ij, 1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 278) )
  if (anyExceptions()) return
#line 279 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 279 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(expected_EV - expected_DeltaSO, profile_2d(ij, 2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 279) )
  if (anyExceptions()) return
#line 280 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 280 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(expected_EC, profile_2d(ij, 3), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 280) )
  if (anyExceptions()) return
#line 281 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    end do

    do ij = 1, 11
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_profile_2d_band_edges

  ! ==================================================================
  ! CSR add test
  ! ==================================================================

  !@test
  subroutine test_csr_add_basic()
    ! Test that csr_add correctly merges two diagonal matrices.
    integer, parameter :: n = 3
    type(csr_matrix) :: A, B, C
    complex(kind=dp) :: val_A(3), val_B(3)
    integer :: rows(3), cols(3)
    complex(kind=dp), allocatable :: dense(:,:)

    rows = [1, 2, 3]
    cols = [1, 2, 3]
    val_A = [cmplx(1.0_dp,0,dp), cmplx(2.0_dp,0,dp), cmplx(3.0_dp,0,dp)]
    val_B = [cmplx(4.0_dp,0,dp), cmplx(5.0_dp,0,dp), cmplx(6.0_dp,0,dp)]

    call csr_build_from_coo(A, n, n, 3, rows, cols, val_A)
    call csr_build_from_coo(B, n, n, 3, rows, cols, val_B)
    call csr_add(A, B, C)

#line 312 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(n, C%nrows, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 312) )
  if (anyExceptions()) return
#line 313 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 313 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(n, C%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 313) )
  if (anyExceptions()) return
#line 314 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    allocate(dense(n, n))
    call csr_to_dense(C, dense)

#line 318 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(5.0_dp,0,dp), dense(1,1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 318) )
  if (anyExceptions()) return
#line 319 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 319 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(7.0_dp,0,dp), dense(2,2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 319) )
  if (anyExceptions()) return
#line 320 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 320 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(9.0_dp,0,dp), dense(3,3), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 320) )
  if (anyExceptions()) return
#line 321 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    call csr_free(A)
    call csr_free(B)
    call csr_free(C)
    deallocate(dense)
  end subroutine test_csr_add_basic

  ! ==================================================================
  ! CSR variable coefficient test
  ! ==================================================================

  !@test
  subroutine test_csr_apply_variable_coeff()
    ! Test row scaling: result(i,j) = -profile(i) * A(i,j)
    integer, parameter :: n = 3
    type(csr_matrix) :: A, result_mat
    complex(kind=dp) :: vals(7)
    integer :: rows(7), cols(7)
    real(kind=dp) :: profile(3)
    complex(kind=dp), allocatable :: dense(:,:)

    rows = [1, 1, 2, 2, 2, 3, 3]
    cols = [1, 2, 1, 2, 3, 2, 3]
    vals = [cmplx(1.0_dp,0,dp), cmplx(2.0_dp,0,dp), cmplx(3.0_dp,0,dp), cmplx(4.0_dp,0,dp), cmplx(5.0_dp,0,dp), cmplx(6.0_dp,0,dp), cmplx(7.0_dp,0,dp)]

    call csr_build_from_coo(A, n, n, 7, rows, cols, vals)
    profile = [10.0_dp, 20.0_dp, 30.0_dp]
    call csr_apply_variable_coeff(A, profile, result_mat)

    allocate(dense(n, n))
    call csr_to_dense(result_mat, dense)

#line 353 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-10.0_dp,0,dp), dense(1,1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 353) )
  if (anyExceptions()) return
#line 354 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 354 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-20.0_dp,0,dp), dense(1,2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 354) )
  if (anyExceptions()) return
#line 355 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 355 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-60.0_dp,0,dp), dense(2,1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 355) )
  if (anyExceptions()) return
#line 356 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 356 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-80.0_dp,0,dp), dense(2,2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 356) )
  if (anyExceptions()) return
#line 357 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 357 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-100.0_dp,0,dp), dense(2,3), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 357) )
  if (anyExceptions()) return
#line 358 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 358 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-180.0_dp,0,dp), dense(3,2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 358) )
  if (anyExceptions()) return
#line 359 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 359 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-210.0_dp,0,dp), dense(3,3), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 359) )
  if (anyExceptions()) return
#line 360 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    call csr_free(A)
    call csr_free(result_mat)
    deallocate(dense)
  end subroutine test_csr_apply_variable_coeff

  ! ==================================================================
  ! 2D Laplacian eigenvalue sanity test
  ! ==================================================================

  !@test
  subroutine test_laplacian_2d_eigenvalue_sanity()
    ! Build kpterms_2d and verify term 5 diagonal is correct sign/magnitude.
    ! For uniform GaAs with A > 0 and csr_apply_variable_coeff scaling by -A,
    ! the diagonal of term 5 should be positive (negative Laplacian times -A).
    integer, parameter :: nnx = 7, nny = 7
    real(kind=dp), parameter :: dx = 1.0_dp, dy = 1.0_dp
    integer :: ngrid, ij, diag_idx
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    type(region_spec), allocatable :: regions(:)
    character(len=255), allocatable :: material_names(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    integer :: mat2d(nnx, nny)
    real(kind=dp) :: diag_val, expected_diag
    real(kind=dp) :: gaas_A

    ngrid = nnx * nny
    mat2d = 1

    call grid_init_rect(grid, nnx, nny, dx, dy, mat2d)

    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)
    gaas_A = params(1)%A

    allocate(regions(1))
    regions(1)%material = "GaAs"

    call confinementInitialization_2d(grid, params, regions, profile_2d, kpterms_2d, FDorder=2)

    ! For interior point, the 2nd derivative FD stencil gives -2/h^2 on diagonal.
    ! csr_apply_variable_coeff scales by -profile, so:
    ! diagonal = -A * (-2/dx^2 + -2/dy^2) = A * (2/dx^2 + 2/dy^2)
    ! Pick an interior point (e.g. row at center of grid)
    diag_idx = (nny/2 - 1) * nnx + nnx/2  ! interior grid point
    diag_val = 0.0_dp
    do ij = kpterms_2d(5)%rowptr(diag_idx), kpterms_2d(5)%rowptr(diag_idx+1) - 1
      if (kpterms_2d(5)%colind(ij) == diag_idx) then
        diag_val = real(kpterms_2d(5)%values(ij))
        exit
      end if
    end do

    expected_diag = gaas_A * (2.0_dp/dx**2 + 2.0_dp/dy**2)

#line 419 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(diag_val > 0.0_dp, message="Term 5 interior diagonal is positive", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 419) )
  if (anyExceptions()) return
#line 420 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 420 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(expected_diag, diag_val, tolerance=0.01_dp * expected_diag, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 420) )
  if (anyExceptions()) return
#line 421 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    do ij = 1, 11
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_laplacian_2d_eigenvalue_sanity

  ! ==================================================================
  ! Helper subroutines
  ! ==================================================================

  subroutine csr_to_dense(mat, dense)
    type(csr_matrix), intent(in) :: mat
    complex(kind=dp), intent(out) :: dense(:,:)
    integer :: row, k

    dense = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do row = 1, mat%nrows
      do k = mat%rowptr(row), mat%rowptr(row + 1) - 1
        dense(row, mat%colind(k)) = mat%values(k)
      end do
    end do
  end subroutine csr_to_dense

  ! ------------------------------------------------------------------
  ! Check symmetry for interior-interior pairs only.
  ! Interior rows/cols are those whose (ix,iy) indices are at least
  ! half_bw+1 away from boundaries.  Only pairs where BOTH ij and j
  ! are interior are checked, because boundary stencils break symmetry
  ! even when the interior row connects to a boundary column.
  ! ------------------------------------------------------------------
  subroutine assert_csr_interior_symmetric(mat, nx, ny, label)
    type(csr_matrix), intent(in) :: mat
    integer, intent(in) :: nx, ny
    character(len=*), intent(in) :: label
    complex(kind=dp), allocatable :: dense(:,:)
    real(kind=dp) :: max_err
    integer :: j, ix, iy, ij, jx, jy

    allocate(dense(mat%nrows, mat%ncols))
    call csr_to_dense(mat, dense)
    max_err = 0.0_dp

    do iy = 2, ny - 1
      do ix = 2, nx - 1
        ij = (iy - 1) * nx + ix
        do jy = 2, ny - 1
          do jx = 2, nx - 1
            j = (jy - 1) * nx + jx
            max_err = max(max_err, abs(dense(ij, j) - dense(j, ij)))
          end do
        end do
      end do
    end do

    deallocate(dense)
#line 478 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(max_err < 1.0e-10_dp, message=trim(label) // " interior symmetric", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 478) )
  if (anyExceptions()) return
#line 479 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  end subroutine assert_csr_interior_symmetric

  ! ------------------------------------------------------------------
  ! Check Hermiticity for interior-interior pairs only.
  ! ------------------------------------------------------------------
  subroutine assert_csr_interior_hermitian(mat, nx, ny, label)
    type(csr_matrix), intent(in) :: mat
    integer, intent(in) :: nx, ny
    character(len=*), intent(in) :: label
    complex(kind=dp), allocatable :: dense(:,:)
    real(kind=dp) :: max_err
    integer :: j, ix, iy, ij, jx, jy

    allocate(dense(mat%nrows, mat%ncols))
    call csr_to_dense(mat, dense)
    max_err = 0.0_dp

    do iy = 2, ny - 1
      do ix = 2, nx - 1
        ij = (iy - 1) * nx + ix
        do jy = 2, ny - 1
          do jx = 2, nx - 1
            j = (jy - 1) * nx + jx
            max_err = max(max_err, abs(dense(ij, j) - conjg(dense(j, ij))))
          end do
        end do
      end do
    end do

    deallocate(dense)
#line 509 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(max_err < 1.0e-10_dp, message=trim(label) // " interior Hermitian", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 509) )
  if (anyExceptions()) return
#line 510 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  end subroutine assert_csr_interior_hermitian

  ! ------------------------------------------------------------------
  subroutine assert_csr_hermitian(mat, label)
    type(csr_matrix), intent(in) :: mat
    character(len=*), intent(in) :: label
    complex(kind=dp), allocatable :: dense(:,:)
    real(kind=dp) :: max_err
    integer :: i, j

    allocate(dense(mat%nrows, mat%ncols))
    call csr_to_dense(mat, dense)
    max_err = 0.0_dp

    do j = 1, mat%ncols
      do i = 1, mat%nrows
        max_err = max(max_err, abs(dense(i,j) - conjg(dense(j,i))))
      end do
    end do

    deallocate(dense)
#line 531 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(max_err < 1.0e-10_dp, message=trim(label) // " Hermitian", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 531) )
  if (anyExceptions()) return
#line 532 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  end subroutine assert_csr_hermitian

end module test_hamiltonian_2d

module Wraptest_hamiltonian_2d
   use FUnit
   use test_hamiltonian_2d
   implicit none
   private

contains


end module Wraptest_hamiltonian_2d

function test_hamiltonian_2d_suite() result(suite)
   use FUnit
   use test_hamiltonian_2d
   use Wraptest_hamiltonian_2d
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   suite = TestSuite('test_hamiltonian_2d_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_kron_identity_size', &
      test_kron_identity_size))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_laplacian_2d_sparsity', &
      test_laplacian_2d_sparsity))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_cross_derivative_pattern', &
      test_cross_derivative_pattern))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_kpterms_2d_sparsity', &
      test_kpterms_2d_sparsity))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_kpterms_2d_hermitian', &
      test_kpterms_2d_hermitian))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_profile_2d_band_edges', &
      test_profile_2d_band_edges))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_csr_add_basic', &
      test_csr_add_basic))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_csr_apply_variable_coeff', &
      test_csr_apply_variable_coeff))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_laplacian_2d_eigenvalue_sanity', &
      test_laplacian_2d_eigenvalue_sanity))
   call suite%addTest(t)


end function test_hamiltonian_2d_suite

