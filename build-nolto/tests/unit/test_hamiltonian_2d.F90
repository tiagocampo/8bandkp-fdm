module test_hamiltonian_2d
  use funit
  use definitions
  use parameters
  use sparse_matrices
  use geometry
  use finitedifferences
  use hamiltonianConstructor
  use confinement_init, only: confinementInitialization_2d
  use hamiltonian_wire
  use linalg, only: zheevx, ilaenv, dlamch
  use csr_test_helpers, only: csr_to_dense, csr_hermitian_error, &
    csr_interior_symmetric_error, csr_interior_hermitian_error, &
    assert_csr_structural_invariants
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
#line 37 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(assert_csr_structural_invariants(C, .true.), message="structural invariants", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 37) )
  if (anyExceptions()) return
#line 38 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

#line 39 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(n * n, C%nrows, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 39) )
  if (anyExceptions()) return
#line 40 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 40 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(n * n, C%ncols, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 40) )
  if (anyExceptions()) return
#line 41 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 41 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(n * n, C%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 41) )
  if (anyExceptions()) return
#line 42 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
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
    type(csr_matrix) :: kron_Iy_D2x, kron_D2y_Ix, lap

    ngrid = nx * ny

    call buildFD2ndDerivMatrix(nx, dx, 2, D2x)
    call buildFD2ndDerivMatrix(ny, dy, 2, D2y)
    call Identity(nx, Ix_mat)
    call Identity(ny, Iy_mat)

    allocate(cD2x(nx,nx)); cD2x = cmplx(D2x, 0.0_dp, kind=dp)
    allocate(cD2y(ny,ny)); cD2y = cmplx(D2y, 0.0_dp, kind=dp)
    allocate(cIx(nx,nx));   cIx = cmplx(Ix_mat, 0.0_dp, kind=dp)
    allocate(cIy(ny,ny));   cIy = cmplx(Iy_mat, 0.0_dp, kind=dp)

    ! Laplacian = Iy x D2x + D2y x Ix  (column-major)
    call kron_eye_dense(ny, cD2x, nx, nx, kron_Iy_D2x)
#line 76 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(assert_csr_structural_invariants(kron_Iy_D2x), message="kron_Iy_D2x structural invariants", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 76) )
  if (anyExceptions()) return
#line 77 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    call kron_dense_dense(cD2y, ny, ny, cIx, nx, nx, kron_D2y_Ix)
#line 78 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(assert_csr_structural_invariants(kron_D2y_Ix), message="kron_D2y_Ix structural invariants", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 78) )
  if (anyExceptions()) return
#line 79 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    call csr_add(kron_Iy_D2x, kron_D2y_Ix, lap)

#line 81 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, lap%nrows, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 81) )
  if (anyExceptions()) return
#line 82 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 82 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, lap%ncols, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 82) )
  if (anyExceptions()) return
#line 83 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 83 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(lap%nnz > 0, message="Laplacian has nonzeros", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 84 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(lap%nnz <= 2 * 3 * ngrid, message="Laplacian sparsity bounded", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Interior rows of 2D Laplacian should be symmetric
#line 87 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(csr_interior_symmetric_error(lap, nx, ny) < 1.0e-10_dp, message="2D Laplacian interior symmetric", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 87) )
  if (anyExceptions()) return
#line 88 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    call csr_free(kron_Iy_D2x)
    call csr_free(kron_D2y_Ix)
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
    type(csr_matrix) :: kron_D1y_D1x

    call buildFD1stDerivMatrix(nx, dx, 2, D1x)
    call buildFD1stDerivMatrix(ny, dy, 2, D1y)

    allocate(cD1x(nx,nx)); cD1x = cmplx(D1x, 0.0_dp, kind=dp)
    allocate(cD1y(ny,ny)); cD1y = cmplx(D1y, 0.0_dp, kind=dp)

    call kron_dense_dense(cD1y, ny, ny, cD1x, nx, nx, kron_D1y_D1x)  ! cross-derivative: D1y x D1x (column-major)
#line 115 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(assert_csr_structural_invariants(kron_D1y_D1x), message="kron_D1y_D1x structural invariants", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 115) )
  if (anyExceptions()) return
#line 116 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

#line 117 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(nx * ny, kron_D1y_D1x%nrows, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 117) )
  if (anyExceptions()) return
#line 118 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 118 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(nx * ny, kron_D1y_D1x%ncols, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 118) )
  if (anyExceptions()) return
#line 119 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 119 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(kron_D1y_D1x%nnz > 0, message="Cross derivative has nonzeros", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 119) )
  if (anyExceptions()) return
#line 120 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    call csr_free(kron_D1y_D1x)
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
#line 159 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(17, size(kpterms_2d), &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 159) )
  if (anyExceptions()) return
#line 160 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Check all matrices have correct dimensions
    do ij = 1, 17
#line 163 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, kpterms_2d(ij)%nrows, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 163) )
  if (anyExceptions()) return
#line 164 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 164 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, kpterms_2d(ij)%ncols, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 164) )
  if (anyExceptions()) return
#line 165 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    end do

    ! Diagonal terms (1-4, 10) should have exactly ngrid nonzeros
#line 168 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, kpterms_2d(1)%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 168) )
  if (anyExceptions()) return
#line 169 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 169 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, kpterms_2d(2)%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 169) )
  if (anyExceptions()) return
#line 170 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 170 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, kpterms_2d(3)%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 170) )
  if (anyExceptions()) return
#line 171 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 171 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, kpterms_2d(4)%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 171) )
  if (anyExceptions()) return
#line 172 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 172 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, kpterms_2d(10)%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 172) )
  if (anyExceptions()) return
#line 173 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Laplacian-based terms (5, 7, 8) should have NNZ > ngrid
#line 175 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(kpterms_2d(5)%nnz > ngrid, message="Term 5 has off-diagonal", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 175) )
  if (anyExceptions()) return
#line 176 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 176 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(kpterms_2d(7)%nnz > ngrid, message="Term 7 has off-diagonal", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 176) )
  if (anyExceptions()) return
#line 177 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 177 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(kpterms_2d(8)%nnz > ngrid, message="Term 8 has off-diagonal", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 177) )
  if (anyExceptions()) return
#line 178 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Cross-derivative (11) should have nonzeros
#line 180 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(kpterms_2d(11)%nnz > 0, message="Term 11 has entries", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 180) )
  if (anyExceptions()) return
#line 181 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Separate x/y gradient terms (12-15) should have nonzeros
#line 183 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(kpterms_2d(12)%nnz > 0, message="Term 12 (P*d/dx) has entries", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 183) )
  if (anyExceptions()) return
#line 184 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 184 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(kpterms_2d(13)%nnz > 0, message="Term 13 (P*d/dy) has entries", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 184) )
  if (anyExceptions()) return
#line 185 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 185 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(kpterms_2d(14)%nnz > 0, message="Term 14 (g3*d/dx) has entries", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 185) )
  if (anyExceptions()) return
#line 186 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 186 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(kpterms_2d(15)%nnz > 0, message="Term 15 (g3*d/dy) has entries", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 186) )
  if (anyExceptions()) return
#line 187 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Anisotropic 2nd derivative terms (16-17) should have nonzeros
#line 189 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(kpterms_2d(16)%nnz > 0, message="Term 16 (g2*lap_diff) has entries", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 189) )
  if (anyExceptions()) return
#line 190 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Check profile dimensions
#line 192 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ngrid, size(profile_2d, 1), &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 192) )
  if (anyExceptions()) return
#line 193 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 193 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(3, size(profile_2d, 2), &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 193) )
  if (anyExceptions()) return
#line 194 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Profile should be non-zero for active cells
#line 196 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(profile_2d(1, 1) /= 0.0_dp, message="EV non-zero", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 196) )
  if (anyExceptions()) return
#line 197 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 197 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(profile_2d(1, 3) /= 0.0_dp, message="EC non-zero", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 197) )
  if (anyExceptions()) return
#line 198 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Cleanup
    do ij = 1, size(kpterms_2d)
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
#line 242 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(csr_interior_hermitian_error(kpterms_2d(5), nnx, nny) < 1.0e-10_dp, message="Term 5 (A*Laplacian) interior Hermitian", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 242) )
  if (anyExceptions()) return
#line 243 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Term 7: (gamma1-2gamma2) * Laplacian
#line 245 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(csr_interior_hermitian_error(kpterms_2d(7), nnx, nny) < 1.0e-10_dp, message="Term 7 (Q*Laplacian) interior Hermitian", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 245) )
  if (anyExceptions()) return
#line 246 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Term 8: (gamma1+2gamma2) * Laplacian
#line 248 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(csr_interior_hermitian_error(kpterms_2d(8), nnx, nny) < 1.0e-10_dp, message="Term 8 (T*Laplacian) interior Hermitian", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 248) )
  if (anyExceptions()) return
#line 249 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    do ij = 1, size(kpterms_2d)
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
#line 297 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(expected_EV, profile_2d(ij, 1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 297) )
  if (anyExceptions()) return
#line 298 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 298 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(expected_EV - expected_DeltaSO, profile_2d(ij, 2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 298) )
  if (anyExceptions()) return
#line 299 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 299 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(expected_EC, profile_2d(ij, 3), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 299) )
  if (anyExceptions()) return
#line 300 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    end do

    do ij = 1, size(kpterms_2d)
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

#line 331 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(n, C%nrows, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 331) )
  if (anyExceptions()) return
#line 332 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 332 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(n, C%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 332) )
  if (anyExceptions()) return
#line 333 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    allocate(dense(n, n))
    call csr_to_dense(C, dense)

#line 337 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(5.0_dp,0,dp), dense(1,1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 337) )
  if (anyExceptions()) return
#line 338 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 338 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(7.0_dp,0,dp), dense(2,2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 338) )
  if (anyExceptions()) return
#line 339 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 339 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(9.0_dp,0,dp), dense(3,3), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 339) )
  if (anyExceptions()) return
#line 340 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

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
    ! Test variable-coefficient scaling:
    !   diagonal:     result(i,i) = -profile(i) * A(i,i)
    !   off-diagonal: result(i,j) = -0.5*(profile(i)+profile(j)) * A(i,j)
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
#line 370 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(assert_csr_structural_invariants(result_mat), message="result_mat structural invariants", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 370) )
  if (anyExceptions()) return
#line 371 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    allocate(dense(n, n))
    call csr_to_dense(result_mat, dense)

#line 375 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-10.0_dp,0,dp), dense(1,1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 375) )
  if (anyExceptions()) return
#line 376 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 376 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-30.0_dp,0,dp), dense(1,2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 376) )
  if (anyExceptions()) return
#line 377 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 377 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-45.0_dp,0,dp), dense(2,1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 377) )
  if (anyExceptions()) return
#line 378 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 378 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-80.0_dp,0,dp), dense(2,2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 378) )
  if (anyExceptions()) return
#line 379 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 379 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-125.0_dp,0,dp), dense(2,3), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 379) )
  if (anyExceptions()) return
#line 380 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 380 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-150.0_dp,0,dp), dense(3,2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 380) )
  if (anyExceptions()) return
#line 381 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 381 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-210.0_dp,0,dp), dense(3,3), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 381) )
  if (anyExceptions()) return
#line 382 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    call csr_free(A)
    call csr_free(result_mat)
    deallocate(dense)
  end subroutine test_csr_apply_variable_coeff

  ! ==================================================================
  ! CSR variable coefficient with face fraction weighting
  ! ==================================================================

  !@test
  subroutine test_csr_apply_variable_coeff_face_frac()
    ! Build a small CSR matrix manually to test face fraction weighting.
    !
    ! 2D grid with nx=3, ny=3, ngrid=9, column-major flat index:
    !   ij = (iy-1)*nx + ix
    !   ix=1,iy=1 -> ij=1    ix=2,iy=1 -> ij=2    ix=3,iy=1 -> ij=3
    !   ix=1,iy=2 -> ij=4    ix=2,iy=2 -> ij=5    ix=3,iy=2 -> ij=6
    !   ix=1,iy=3 -> ij=7    ix=2,iy=3 -> ij=8    ix=3,iy=3 -> ij=9
    !
    ! We build a 5-entry CSR row for row 5 (center, ix=2, iy=2):
    !   (5,2): offset -3 -> y-direction
    !   (5,4): offset -1 -> x-direction
    !   (5,5): diagonal
    !   (5,6): offset +1 -> x-direction
    !   (5,8): offset +3 -> y-direction
    !
    ! With face_frac_x(5,:)=[0.5, 0.5], face_frac_y(5,:)=[0.8, 0.8]:
    !   x-avg = 0.5, y-avg = 0.8
    !   Profile = 1.0 -> scale = -1
    !   Off-diagonal: (5,2)=-0.8, (5,4)=-0.5, (5,6)=-0.5, (5,8)=-0.8
    !   Diagonal: negative sum of off-diags = -(-2.6) = 2.6 (row-sum = 0)
    !
    ! All other rows have face fractions = 1.0 (no extra scaling).
    integer, parameter :: ngrid = 9, nx = 3
    type(csr_matrix) :: A, result_mat
    complex(kind=dp) :: vals(5)
    integer :: rows(5), cols(5)
    real(kind=dp) :: profile(ngrid)
    real(kind=dp) :: ffx(ngrid, 2), ffy(ngrid, 2)
    complex(kind=dp), allocatable :: dense(:,:)

    ! Build a CSR matrix with only the 5 entries for row 5
    rows = [5, 5, 5, 5, 5]
    cols = [2, 4, 5, 6, 8]
    vals = [cmplx(1.0_dp,0,dp), cmplx(1.0_dp,0,dp), cmplx(1.0_dp,0,dp), &
            cmplx(1.0_dp,0,dp), cmplx(1.0_dp,0,dp)]

    call csr_build_from_coo(A, ngrid, ngrid, 5, rows, cols, vals)

    ! Profile = 1 -> scale = -1
    profile = 1.0_dp

    ! Face fractions: all 1.0 except row 5
    ffx = 1.0_dp
    ffy = 1.0_dp
    ffx(5, 1) = 0.5_dp;  ffx(5, 2) = 0.5_dp
    ffy(5, 1) = 0.8_dp;  ffy(5, 2) = 0.8_dp

    call csr_apply_variable_coeff(A, profile, result_mat, &
      face_frac_x=ffx, face_frac_y=ffy, grid_nx=nx)
#line 443 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(assert_csr_structural_invariants(result_mat), message="result_mat structural invariants", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 443) )
  if (anyExceptions()) return
#line 444 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    allocate(dense(ngrid, ngrid))
    dense = cmplx(0.0_dp, 0.0_dp, kind=dp)
    call csr_to_dense(result_mat, dense)

    ! Row 5 checks:
    ! col 2: y-direction, avg(ffy)=0.8 -> -1*0.8 = -0.8
#line 451 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-0.8_dp,0,dp), dense(5,2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 451) )
  if (anyExceptions()) return
#line 452 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    ! col 4: x-direction, avg(ffx)=0.5 -> -1*0.5 = -0.5
#line 453 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-0.5_dp,0,dp), dense(5,4), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 453) )
  if (anyExceptions()) return
#line 454 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    ! col 5: diagonal = negative sum of off-diagonals (row-sum = 0)
    !   off-diag sum = -0.8 + -0.5 + -0.5 + -0.8 = -2.6 -> diag = 2.6
#line 456 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(2.6_dp,0,dp), dense(5,5), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 456) )
  if (anyExceptions()) return
#line 457 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    ! col 6: x-direction -> -0.5
#line 458 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-0.5_dp,0,dp), dense(5,6), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 458) )
  if (anyExceptions()) return
#line 459 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    ! col 8: y-direction -> -0.8
#line 460 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(-0.8_dp,0,dp), dense(5,8), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 460) )
  if (anyExceptions()) return
#line 461 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! All other rows should be zero (no entries in A)
#line 463 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(0.0_dp,0,dp), dense(1,1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 463) )
  if (anyExceptions()) return
#line 464 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 464 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(cmplx(0.0_dp,0,dp), dense(3,3), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 464) )
  if (anyExceptions()) return
#line 465 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    call csr_free(A)
    call csr_free(result_mat)
    deallocate(dense)
  end subroutine test_csr_apply_variable_coeff_face_frac

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
    ! diagonal = -(A*const) * (-2/dx^2 + -2/dy^2) = A*const * (2/dx^2 + 2/dy^2)
    ! Pick an interior point (e.g. row at center of grid)
    diag_idx = (nny/2 - 1) * nnx + nnx/2  ! interior grid point
    diag_val = 0.0_dp
    do ij = kpterms_2d(5)%rowptr(diag_idx), kpterms_2d(5)%rowptr(diag_idx+1) - 1
      if (kpterms_2d(5)%colind(ij) == diag_idx) then
        diag_val = real(kpterms_2d(5)%values(ij))
        exit
      end if
    end do

    expected_diag = gaas_A * const * (2.0_dp/dx**2 + 2.0_dp/dy**2)

#line 524 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(diag_val > 0.0_dp, message="Term 5 interior diagonal is positive", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 524) )
  if (anyExceptions()) return
#line 525 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 525 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(expected_diag, diag_val, tolerance=0.01_dp * expected_diag, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 525) )
  if (anyExceptions()) return
#line 526 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    do ij = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_laplacian_2d_eigenvalue_sanity

  !@test
  subroutine test_rect_wire_keeps_dirichlet_boundary_laplacian()
    ! Rectangular wires use explicit Dirichlet boundary closures.  The
    ! variable-coefficient application must not replace those boundary
    ! diagonals with a row-sum-conserving Neumann-style diagonal.
    integer, parameter :: nnx = 5, nny = 5
    real(kind=dp), parameter :: dx = 3.0_dp, dy = 3.0_dp
    integer :: ij, corner_idx
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    type(region_spec), allocatable :: regions(:)
    character(len=255), allocatable :: material_names(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    integer :: mat2d(nnx, nny)
    real(kind=dp) :: diag_val, expected_diag
    real(kind=dp) :: gaas_A

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

    corner_idx = 1
    diag_val = 0.0_dp
    do ij = kpterms_2d(5)%rowptr(corner_idx), kpterms_2d(5)%rowptr(corner_idx + 1) - 1
      if (kpterms_2d(5)%colind(ij) == corner_idx) then
        diag_val = real(kpterms_2d(5)%values(ij))
        exit
      end if
    end do

    expected_diag = gaas_A * const * (3.0_dp/dx**2 + 3.0_dp/dy**2)
#line 576 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(expected_diag, diag_val, tolerance=1.0e-12_dp * expected_diag, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 576) )
  if (anyExceptions()) return
#line 577 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    do ij = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_rect_wire_keeps_dirichlet_boundary_laplacian

  !@test
  subroutine test_rect_wire_first_derivative_boundary_sign()
    ! First-derivative operators are not conservative Laplacians.  Their
    ! diagonal entries must retain the explicit Dirichlet closure sign.
    integer, parameter :: nnx = 5, nny = 5
    real(kind=dp), parameter :: dx = 3.0_dp, dy = 3.0_dp
    integer :: ij, corner_idx
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    type(region_spec), allocatable :: regions(:)
    character(len=255), allocatable :: material_names(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    integer :: mat2d(nnx, nny)
    real(kind=dp) :: diag_val, expected_diag
    real(kind=dp) :: gaas_P

    mat2d = 1
    call grid_init_rect(grid, nnx, nny, dx, dy, mat2d)

    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)
    gaas_P = params(1)%P

    allocate(regions(1))
    regions(1)%material = "GaAs"

    call confinementInitialization_2d(grid, params, regions, profile_2d, kpterms_2d, FDorder=2)

    corner_idx = 1
    diag_val = 0.0_dp
    do ij = kpterms_2d(12)%rowptr(corner_idx), kpterms_2d(12)%rowptr(corner_idx + 1) - 1
      if (kpterms_2d(12)%colind(ij) == corner_idx) then
        diag_val = real(kpterms_2d(12)%values(ij))
        exit
      end if
    end do

    expected_diag = -gaas_P * 0.5_dp / dx
#line 626 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(expected_diag, diag_val, tolerance=1.0e-12_dp * abs(expected_diag), &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 626) )
  if (anyExceptions()) return
#line 627 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    do ij = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_rect_wire_first_derivative_boundary_sign

  !@test
  subroutine test_inactive_wire_cells_are_hard_wall_barriers()
    ! Non-rectangular wire grids keep the rectangular storage layout, so cells
    ! outside the active cross-section must not remain zero-energy states.
    integer, parameter :: nnx = 7, nny = 7
    real(kind=dp), parameter :: dx = 1.0_dp, dy = 1.0_dp
    integer :: ij, corner_idx, cb_diag_idx, vb_diag_idx, ngrid
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    type(region_spec), allocatable :: regions(:)
    character(len=255), allocatable :: material_names(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    type(csr_matrix) :: HT_csr
    type(simulation_config) :: cfg
    integer :: mat2d(nnx, nny)
    complex(kind=dp) :: cb_diag, vb_diag

    mat2d = 1
    call grid_init_circle(grid, nnx, nny, dx, dy, 0.5_dp * nnx * dx, &
      0.5_dp * nny * dy, 2.0_dp, mat2d)

    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)

    allocate(regions(1))
    regions(1)%material = "GaAs"

    call confinementInitialization_2d(grid, params, regions, profile_2d, kpterms_2d, FDorder=2)

    cfg%grid = grid
    cfg%confinement = 'wire'
    call ZB8bandGeneralized(HT_csr, 0.0_dp, profile_2d, kpterms_2d, cfg)

    ngrid = nnx * nny
    corner_idx = 1
    vb_diag_idx = corner_idx
    cb_diag_idx = 6 * ngrid + corner_idx
    vb_diag = cmplx(0.0_dp, 0.0_dp, kind=dp)
    cb_diag = cmplx(0.0_dp, 0.0_dp, kind=dp)

    do ij = HT_csr%rowptr(vb_diag_idx), HT_csr%rowptr(vb_diag_idx + 1) - 1
      if (HT_csr%colind(ij) == vb_diag_idx) then
        vb_diag = HT_csr%values(ij)
        exit
      end if
    end do
    do ij = HT_csr%rowptr(cb_diag_idx), HT_csr%rowptr(cb_diag_idx + 1) - 1
      if (HT_csr%colind(ij) == cb_diag_idx) then
        cb_diag = HT_csr%values(ij)
        exit
      end if
    end do

#line 691 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(real(vb_diag) < -100.0_dp, message="Inactive valence cell has hard-wall barrier", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 691) )
  if (anyExceptions()) return
#line 692 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 692 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(real(cb_diag) > 100.0_dp, message="Inactive conduction cell has hard-wall barrier", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 692) )
  if (anyExceptions()) return
#line 693 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    call csr_free(HT_csr)
    do ij = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_inactive_wire_cells_are_hard_wall_barriers

  !@test
  subroutine test_heterostructure_wire_hamiltonian_is_hermitian()
    ! Material discontinuities must not make the wire Hamiltonian depend on
    ! which triangle the Hermitian eigensolver reads.
    integer, parameter :: nnx = 5, nny = 5
    real(kind=dp), parameter :: dx = 3.0_dp, dy = 3.0_dp
    integer :: ij
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    type(region_spec), allocatable :: regions(:)
    character(len=255), allocatable :: material_names(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    type(csr_matrix) :: HT_csr
    type(simulation_config) :: cfg
    integer :: mat2d(nnx, nny)

    mat2d = 1
    mat2d(3:5, :) = 2
    call grid_init_rect(grid, nnx, nny, dx, dy, mat2d)

    allocate(material_names(2))
    material_names(1) = "GaAs"
    material_names(2) = "InAs"
    allocate(params(2))
    call paramDatabase(material_names, 2, params)

    allocate(regions(2))
    regions(1)%material = "GaAs"
    regions(2)%material = "InAs"

    call confinementInitialization_2d(grid, params, regions, profile_2d, kpterms_2d, FDorder=2)

    cfg%grid = grid
    cfg%confinement = 'wire'
    call ZB8bandGeneralized(HT_csr, 0.0_dp, profile_2d, kpterms_2d, cfg)

#line 739 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(csr_hermitian_error(HT_csr) < 1.0e-10_dp, message="Heterostructure wire Hamiltonian Hermitian", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 739) )
  if (anyExceptions()) return
#line 740 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    call csr_free(HT_csr)
    do ij = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_heterostructure_wire_hamiltonian_is_hermitian

  ! ==================================================================
  ! Non-square grid full Hamiltonian test (P1-D)
  !
  ! Uses nx != ny to catch flat-index convention mismatches between
  ! geometry (flat_idx) and Kronecker product assembly.  Builds the
  ! full 8-band wire Hamiltonian via ZB8bandGeneralized and verifies:
  !   1. Correct total dimension (8 * nx * ny)
  !   2. Hermiticity of the full Hamiltonian
  !   3. Non-zero eigenvalues (not all degenerate)
  ! ==================================================================

  !@test
  subroutine test_ZB8bandGeneralized_nonsquare()
    integer, parameter :: nx = 4, ny = 6
    real(kind=dp), parameter :: dx = 2.0_dp, dy = 2.0_dp
    real(kind=dp), parameter :: kz_test = 0.1_dp
    integer :: ngrid, ntot, ij
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    type(region_spec), allocatable :: regions(:)
    character(len=255), allocatable :: material_names(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    type(csr_matrix) :: H_csr
    type(simulation_config) :: cfg
    integer :: mat2d(nx, ny)

    ngrid = nx * ny
    ntot = 8 * ngrid
    mat2d = 1  ! single material (GaAs)

    ! Build grid
    call grid_init_rect(grid, nx, ny, dx, dy, mat2d)

    ! Setup materials
    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)

    allocate(regions(1))
    regions(1)%material = "GaAs"

    ! Build kpterms_2d
    call confinementInitialization_2d(grid, params, regions, profile_2d, kpterms_2d, FDorder=2)

    ! Minimal config for ZB8bandGeneralized
    cfg%grid = grid
    cfg%confinement = 'wire'

    ! Build full Hamiltonian
    call ZB8bandGeneralized(H_csr, kz_test, profile_2d, kpterms_2d, cfg)

    ! Check dimensions
#line 803 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ntot, H_csr%nrows, message="H nrows = 8*N", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 803) )
  if (anyExceptions()) return
#line 804 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 804 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(ntot, H_csr%ncols, message="H ncols = 8*N", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 804) )
  if (anyExceptions()) return
#line 805 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 805 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(H_csr%nnz > 0, message="H has nonzeros", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 805) )
  if (anyExceptions()) return
#line 806 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! For the default second-order hard-wall wire discretization, the full
    ! Hamiltonian should now be Hermitian, including boundary rows.
    block
      complex(kind=dp), allocatable :: Hd(:,:)
      real(kind=dp) :: herr, herr_full
      integer :: i, j
      integer :: si, sj, six, siy, sjx, sjy
      logical :: i_interior, j_interior

      allocate(Hd(H_csr%nrows, H_csr%ncols))
      call csr_to_dense(H_csr, Hd)
      herr = 0.0_dp
      herr_full = 0.0_dp
      do j = 1, H_csr%ncols
        sj = mod(j-1, ngrid) + 1
        sjx = mod(sj-1, nx) + 1
        sjy = (sj-1)/nx + 1
        j_interior = (sjx > 1 .and. sjx < nx .and. sjy > 1 .and. sjy < ny)
        do i = 1, H_csr%nrows
          si = mod(i-1, ngrid) + 1
          six = mod(si-1, nx) + 1
          siy = (si-1)/nx + 1
          i_interior = (six > 1 .and. six < nx .and. siy > 1 .and. siy < ny)
          herr_full = max(herr_full, abs(Hd(i,j) - conjg(Hd(j,i))))
          if (i_interior .and. j_interior) then
            herr = max(herr, abs(Hd(i,j) - conjg(Hd(j,i))))
          end if
        end do
      end do
      print *, '  Wire full-Hermiticity max error (including boundaries): ', herr_full
      deallocate(Hd)
#line 838 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(herr < 1.0e-10_dp, message="Interior wire Hamiltonian (nx=4,ny=6) Hermitian", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 838) )
  if (anyExceptions()) return
#line 839 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 839 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(herr_full < 1.0e-10_dp, message="Full wire Hamiltonian (nx=4,ny=6) Hermitian", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 839) )
  if (anyExceptions()) return
#line 840 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    end block

    ! Cleanup
    call csr_free(H_csr)
    do ij = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_ZB8bandGeneralized_nonsquare

  ! ==================================================================
  ! Wire Hamiltonian diagonal units test
  !
  ! Verifies that the wire Hamiltonian CB diagonal includes the
  ! hbar^2/(2*m0) conversion factor.  Without it, the kinetic term
  ! is ~3.8x too small because gamma/A are dimensionless while the FD
  ! operators have units Ang^-2.
  ! ==================================================================

  !@test
  subroutine test_wire_hamiltonian_diagonal_units()
    integer, parameter :: nx = 5, ny = 5
    real(kind=dp), parameter :: dx = 3.0_dp, dy = 3.0_dp
    integer :: ngrid, ij, band7_diag_idx
    real(kind=dp) :: actual_cb_diag, expected_min

    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    character(len=255), allocatable :: material_names(:)
    type(region_spec), allocatable :: regions(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    type(csr_matrix) :: HT_csr
    type(simulation_config) :: cfg
    integer :: mat2d(nx, ny)

    ngrid = nx * ny
    mat2d = 1  ! all GaAs

    call grid_init_rect(grid, nx, ny, dx, dy, mat2d)

    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)

    allocate(regions(1))
    regions(1)%material = "GaAs"

    call confinementInitialization_2d(grid, params, regions, &
      profile_2d, kpterms_2d, FDorder=2)

    ! Minimal config for ZB8bandGeneralized
    cfg%grid = grid
    cfg%confinement = 'wire'

    call ZB8bandGeneralized(HT_csr, 0.0_dp, profile_2d, kpterms_2d, cfg)

    ! Interior grid point (ix=3, iy=3) in column-major flat index:
    !   ij = (iy-1)*nx + ix = 2*5 + 3 = 13
    ! Band 7 (CB1) offset: 6 * ngrid = 150
    ! Global index: 150 + 13 = 163
    band7_diag_idx = 6 * ngrid + 13

    ! Without const: CB_diag ~ EC + A*(2/dx^2+2/dy^2) ~ 0.72 + 6.6 ~ 7.3
    ! With const:    CB_diag ~ EC + A*const*(2/dx^2+2/dy^2) ~ 0.72 + 25.2 ~ 26.0
    ! Threshold: must exceed EC + 10 = 10.72 to confirm const is present
    expected_min = params(1)%EC + 10.0_dp

    actual_cb_diag = 0.0_dp
    do ij = HT_csr%rowptr(band7_diag_idx), HT_csr%rowptr(band7_diag_idx + 1) - 1
      if (HT_csr%colind(ij) == band7_diag_idx) then
        actual_cb_diag = real(HT_csr%values(ij), kind=dp)
        exit
      end if
    end do

#line 918 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(actual_cb_diag > expected_min, message="Wire CB diagonal includes const factor", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 918) )
  if (anyExceptions()) return
#line 919 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Also verify exact expected value: EC + A*const*(2/dx² + 2/dy²)
    block
      real(kind=dp) :: expected_exact
      expected_exact = params(1)%EC + params(1)%A * hbar2O2m0 * (2.0_dp/dx**2 + 2.0_dp/dy**2)
#line 924 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(expected_exact, actual_cb_diag, tolerance=0.5_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 924) )
  if (anyExceptions()) return
#line 925 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    end block

    call csr_free(HT_csr)
    do ij = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_wire_hamiltonian_diagonal_units

  !@test
  subroutine test_wire_kane_couplings_at_kz_zero()
    integer, parameter :: nx = 5, ny = 5
    real(kind=dp), parameter :: dx = 3.0_dp, dy = 3.0_dp
    integer :: ngrid, ij
    real(kind=dp) :: pz_27, pz_57, pp_17, pm_48

    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    character(len=255), allocatable :: material_names(:)
    type(region_spec), allocatable :: regions(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    type(csr_matrix) :: HT_csr
    complex(kind=dp), allocatable :: Hd(:,:)
    type(simulation_config) :: cfg
    integer :: mat2d(nx, ny)

    ngrid = nx * ny
    mat2d = 1

    call grid_init_rect(grid, nx, ny, dx, dy, mat2d)

    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)

    allocate(regions(1))
    regions(1)%material = "GaAs"

    call confinementInitialization_2d(grid, params, regions, &
      profile_2d, kpterms_2d, FDorder=2)

    cfg%grid = grid
    cfg%confinement = 'wire'

    call ZB8bandGeneralized(HT_csr, 0.0_dp, profile_2d, kpterms_2d, cfg)
    allocate(Hd(HT_csr%nrows, HT_csr%ncols))
    call csr_to_dense(HT_csr, Hd)

    pz_27 = maxval(abs(Hd(ngrid + 1:2*ngrid, 6*ngrid + 1:7*ngrid)))
    pz_57 = maxval(abs(Hd(4*ngrid + 1:5*ngrid, 6*ngrid + 1:7*ngrid)))
    pp_17 = maxval(abs(Hd(1:ngrid, 6*ngrid + 1:7*ngrid)))
    pm_48 = maxval(abs(Hd(3*ngrid + 1:4*ngrid, 7*ngrid + 1:8*ngrid)))
#line 980 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(pz_27 < 1.0e-12_dp, message="Wire PZ coupling vanishes at kz=0", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 980) )
  if (anyExceptions()) return
#line 981 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 981 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(pz_57 < 1.0e-12_dp, message="Wire SO-CB PZ coupling vanishes at kz=0", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 981) )
  if (anyExceptions()) return
#line 982 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 982 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(pp_17 > 1.0e-8_dp, message="Wire PP transverse coupling survives at kz=0", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 982) )
  if (anyExceptions()) return
#line 983 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 983 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(pm_48 > 1.0e-8_dp, message="Wire PM transverse coupling survives at kz=0", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 983) )
  if (anyExceptions()) return
#line 984 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    deallocate(Hd)
    call csr_free(HT_csr)
    do ij = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_wire_kane_couplings_at_kz_zero

  !@test
  subroutine test_wire_cb_block_ground_state_centered()
    integer, parameter :: nx = 11, ny = 11
    real(kind=dp), parameter :: dx = 3.0_dp, dy = 3.0_dp
    integer :: ngrid, ij, info, nb, lwork, m, il, iuu
    real(kind=dp) :: abstol, vl, vu, center_prob, corner_prob
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    character(len=255), allocatable :: material_names(:)
    type(region_spec), allocatable :: regions(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    complex(kind=dp), allocatable :: Ablk(:,:), work(:), Z(:,:)
    real(kind=dp), allocatable :: evals(:), rwork(:)
    integer, allocatable :: iwork(:), ifail(:)
    integer :: mat2d(nx, ny)

    ngrid = nx * ny
    mat2d = 1

    call grid_init_rect(grid, nx, ny, dx, dy, mat2d)

    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)

    allocate(regions(1))
    regions(1)%material = "GaAs"

    call confinementInitialization_2d(grid, params, regions, &
      profile_2d, kpterms_2d, FDorder=2)

    allocate(Ablk(ngrid, ngrid))
    Ablk = cmplx(0.0_dp, 0.0_dp, kind=dp)
    call csr_to_dense(kpterms_2d(5), Ablk)
    do ij = 1, ngrid
      Ablk(ij, ij) = Ablk(ij, ij) + cmplx(profile_2d(ij, 3), 0.0_dp, kind=dp)
    end do

    allocate(evals(ngrid), Z(ngrid, 1))
    vl = 0.0_dp
    vu = 0.0_dp
    il = 1
    iuu = 1
    abstol = 2.0_dp * dlamch('S')
    nb = ilaenv(1, 'ZHETRD', 'U', ngrid, -1, -1, -1)
    lwork = max(1, (nb + 1) * ngrid)
    allocate(work(lwork), rwork(7 * ngrid), iwork(5 * ngrid), ifail(ngrid))

    call zheevx('V', 'I', 'U', ngrid, Ablk, ngrid, vl, vu, il, iuu, abstol, &
      m, evals, Z, ngrid, work, lwork, rwork, iwork, ifail, info)

#line 1047 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(0, info, message="CB block zheevx succeeds", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1047) )
  if (anyExceptions()) return
#line 1048 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 1048 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(1, m, message="CB block returns one eigenpair", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1048) )
  if (anyExceptions()) return
#line 1049 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    corner_prob = abs(Z(1, 1))**2
    center_prob = abs(Z((ny / 2) * nx + (nx / 2) + 1, 1))**2

#line 1053 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(center_prob > corner_prob, message="Isolated CB block ground state is center-enhanced", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1053) )
  if (anyExceptions()) return
#line 1054 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    deallocate(work, rwork, iwork, ifail, evals, Z, Ablk)
    do ij = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_wire_cb_block_ground_state_centered

  ! ==================================================================
  ! wire_workspace free test
  ! ==================================================================

  !@test
  subroutine test_wire_workspace_free()
    type(wire_workspace) :: ws
    call wire_workspace_free(ws)
#line 1071 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertFalse(ws%initialized, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1071) )
  if (anyExceptions()) return
#line 1072 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  end subroutine test_wire_workspace_free

  ! ==================================================================
  ! Wire workspace: slow-path reproducibility test
  ! ==================================================================

  !@test
  subroutine test_workspace_slow_path_reproducible()
    ! Build wire Hamiltonian at one kz twice via slow path,
    ! verify identical CSR output.
    integer, parameter :: nx = 4, ny = 4
    real(kind=dp), parameter :: dx = 5.0_dp, dy = 5.0_dp
    integer :: ngrid, ij, k
    integer :: mat2d(nx, ny)
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    character(len=255), allocatable :: material_names(:)
    type(region_spec), allocatable :: regions(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    type(csr_matrix) :: H1, H2
    type(wire_coo_cache) :: coo_cache
    real(kind=dp) :: kz
    type(simulation_config) :: cfg

    ngrid = nx * ny
    mat2d = 1

    call grid_init_rect(grid, nx, ny, dx, dy, mat2d)
    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)
    allocate(regions(1))
    regions(1)%material = "GaAs"

    call confinementInitialization_2d(grid, params, regions, &
      profile_2d, kpterms_2d, FDorder=2)

    kz = 0.05_dp
    cfg%grid = grid

    ! First call
    coo_cache%initialized = .false.
    call ZB8bandGeneralized(H1, kz, profile_2d, kpterms_2d, cfg, coo_cache)
    call wire_coo_cache_free(coo_cache)

    ! Second call (same kz, fresh cache)
    coo_cache%initialized = .false.
    call ZB8bandGeneralized(H2, kz, profile_2d, kpterms_2d, cfg, coo_cache)

    ! Compare
#line 1124 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(H1%nrows, H2%nrows, message="nrows match", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1124) )
  if (anyExceptions()) return
#line 1125 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 1125 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(H1%nnz, H2%nnz, message="nnz match", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1125) )
  if (anyExceptions()) return
#line 1126 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    do k = 1, H1%nnz
#line 1127 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(H1%colind(k), H2%colind(k), message="colind match", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1127) )
  if (anyExceptions()) return
#line 1128 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    end do
    do k = 1, H1%nnz
#line 1130 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(H1%values(k), H2%values(k), tolerance=1.0e-13_dp, message="values match", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1130) )
  if (anyExceptions()) return
#line 1131 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    end do

    call csr_free(H1)
    call csr_free(H2)
    call wire_coo_cache_free(coo_cache)
    do ij = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_workspace_slow_path_reproducible

  !@test
  subroutine test_workspace_fast_path_matches_slow()
    ! Build wire Hamiltonian at kz=0.05 via slow path, then at kz=0.06
    ! via fast path (workspace initialized from first call).
    ! Verify the fast-path result matches a fresh slow-path build at kz=0.06.
    integer, parameter :: nx = 4, ny = 4
    real(kind=dp), parameter :: dx = 5.0_dp, dy = 5.0_dp
    integer :: ngrid, ij, k
    integer :: mat2d(nx, ny)
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    character(len=255), allocatable :: material_names(:)
    type(region_spec), allocatable :: regions(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    type(csr_matrix) :: H_slow, H_fast
    type(wire_coo_cache) :: coo_cache
    type(wire_workspace) :: wire_ws
    real(kind=dp) :: kz1, kz2
    type(simulation_config) :: cfg

    ngrid = nx * ny
    mat2d = 1

    call grid_init_rect(grid, nx, ny, dx, dy, mat2d)
    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)
    allocate(regions(1))
    regions(1)%material = "GaAs"

    call confinementInitialization_2d(grid, params, regions, &
      profile_2d, kpterms_2d, FDorder=2)

    kz1 = 0.05_dp
    kz2 = 0.06_dp
    cfg%grid = grid

    ! First call at kz1: initializes workspace (slow path)
    call ZB8bandGeneralized(H_slow, kz1, profile_2d, kpterms_2d, cfg, &
      ws=wire_ws)
    call csr_free(H_slow)

    ! Second call at kz2: uses fast path
    call ZB8bandGeneralized(H_fast, kz2, profile_2d, kpterms_2d, cfg, &
      ws=wire_ws)

    ! Fresh slow-path at kz2 for comparison
    coo_cache%initialized = .false.
    call ZB8bandGeneralized(H_slow, kz2, profile_2d, kpterms_2d, cfg, &
      coo_cache)

    ! Compare: must match element-by-element
#line 1197 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(H_slow%nrows, H_fast%nrows, message="nrows match", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1197) )
  if (anyExceptions()) return
#line 1198 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 1198 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(H_slow%nnz, H_fast%nnz, message="nnz match", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1198) )
  if (anyExceptions()) return
#line 1199 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    do k = 1, H_slow%nnz
#line 1200 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(H_slow%colind(k), H_fast%colind(k), message="colind match", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1200) )
  if (anyExceptions()) return
#line 1201 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    end do
    do k = 1, H_slow%nnz
#line 1203 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(H_slow%values(k), H_fast%values(k), tolerance=1.0e-13_dp, message="values match", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1203) )
  if (anyExceptions()) return
#line 1204 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    end do

    call csr_free(H_slow)
    call csr_free(H_fast)
    call wire_coo_cache_free(coo_cache)
    call wire_workspace_free(wire_ws)
    do ij = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_workspace_fast_path_matches_slow

  !@test
  subroutine test_workspace_g3_does_not_corrupt_normal_cache()
    integer, parameter :: nx = 4, ny = 4
    real(kind=dp), parameter :: dx = 5.0_dp, dy = 5.0_dp
    integer :: ij, k, normal_nnz, coo_nnz
    integer :: mat2d(nx, ny)
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    character(len=255), allocatable :: material_names(:)
    type(region_spec), allocatable :: regions(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    type(csr_matrix) :: H_fast, H_slow, H_g3
    type(wire_workspace) :: wire_ws
    type(wire_coo_cache) :: coo_cache
    type(simulation_config) :: cfg

    mat2d = 1
    call grid_init_rect(grid, nx, ny, dx, dy, mat2d)
    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)
    allocate(regions(1))
    regions(1)%material = "GaAs"

    call confinementInitialization_2d(grid, params, regions, &
      profile_2d, kpterms_2d, FDorder=2)
    cfg%grid = grid

    ! Initialize the normal workspace and its COO mapping.
    call ZB8bandGeneralized(H_fast, 0.05_dp, profile_2d, kpterms_2d, cfg, &
      ws=wire_ws)
    normal_nnz = H_fast%nnz
    coo_nnz = wire_ws%coo_nnz_in

    ! Accidental workspace use in g3 mode must not overwrite normal cache state.
    call ZB8bandGeneralized(H_g3, 1.0_dp, profile_2d, kpterms_2d, cfg, &
      ws=wire_ws, g='g3')
#line 1256 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertTrue(H_g3%nnz > 0, message="g3 matrix built", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1256) )
  if (anyExceptions()) return
#line 1257 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 1257 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(normal_nnz, H_fast%nnz, message="normal matrix retained", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1257) )
  if (anyExceptions()) return
#line 1258 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 1258 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(coo_nnz, wire_ws%coo_nnz_in, message="normal COO mapping retained", &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1258) )
  if (anyExceptions()) return
#line 1259 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"

    ! Normal fast path must still match a fresh slow-path build after g3.
    call ZB8bandGeneralized(H_fast, 0.06_dp, profile_2d, kpterms_2d, cfg, &
      ws=wire_ws)
    call ZB8bandGeneralized(H_slow, 0.06_dp, profile_2d, kpterms_2d, cfg, &
      coo_cache)

#line 1266 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(H_slow%nrows, H_fast%nrows, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1266) )
  if (anyExceptions()) return
#line 1267 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 1267 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(H_slow%nnz, H_fast%nnz, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1267) )
  if (anyExceptions()) return
#line 1268 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    do k = 1, H_slow%nnz
#line 1269 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(H_slow%colind(k), H_fast%colind(k), &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1269) )
  if (anyExceptions()) return
#line 1270 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
#line 1270 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
  call assertEqual(H_slow%values(k), H_fast%values(k), tolerance=1.0e-13_dp, &
 & location=SourceLocation( &
 & 'test_hamiltonian_2d.pf', &
 & 1270) )
  if (anyExceptions()) return
#line 1271 "/data/8bandkp-fdm/tests/unit/test_hamiltonian_2d.pf"
    end do

    call csr_free(H_fast)
    call csr_free(H_slow)
    call csr_free(H_g3)
    call wire_workspace_free(wire_ws)
    call wire_coo_cache_free(coo_cache)
    do ij = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_workspace_g3_does_not_corrupt_normal_cache

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
   allocate(t, source=TestMethod('test_csr_apply_variable_coeff_face_frac', &
      test_csr_apply_variable_coeff_face_frac))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_laplacian_2d_eigenvalue_sanity', &
      test_laplacian_2d_eigenvalue_sanity))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_rect_wire_keeps_dirichlet_boundary_laplacian', &
      test_rect_wire_keeps_dirichlet_boundary_laplacian))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_rect_wire_first_derivative_boundary_sign', &
      test_rect_wire_first_derivative_boundary_sign))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_inactive_wire_cells_are_hard_wall_barriers', &
      test_inactive_wire_cells_are_hard_wall_barriers))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_heterostructure_wire_hamiltonian_is_hermitian', &
      test_heterostructure_wire_hamiltonian_is_hermitian))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_ZB8bandGeneralized_nonsquare', &
      test_ZB8bandGeneralized_nonsquare))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_wire_hamiltonian_diagonal_units', &
      test_wire_hamiltonian_diagonal_units))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_wire_kane_couplings_at_kz_zero', &
      test_wire_kane_couplings_at_kz_zero))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_wire_cb_block_ground_state_centered', &
      test_wire_cb_block_ground_state_centered))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_wire_workspace_free', &
      test_wire_workspace_free))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_workspace_slow_path_reproducible', &
      test_workspace_slow_path_reproducible))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_workspace_fast_path_matches_slow', &
      test_workspace_fast_path_matches_slow))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_workspace_g3_does_not_corrupt_normal_cache', &
      test_workspace_g3_does_not_corrupt_normal_cache))
   call suite%addTest(t)


end function test_hamiltonian_2d_suite

