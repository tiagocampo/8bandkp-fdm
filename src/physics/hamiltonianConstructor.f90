module hamiltonianConstructor

  use definitions
  use finitedifferences
  use sparse_matrices
  use utils
  use strain_solver, only: compute_bp_scalar, bir_pikus_blocks_free

  implicit none

  interface confinementInitialization
    module procedure confinementInitialization_raw
    module procedure confinementInitialization_cfg
  end interface confinementInitialization

  ! ------------------------------------------------------------------
  ! Cache for the COO sparsity structure of the wire Hamiltonian.
  !
  ! Stores the mapping from unsorted COO indices to sorted+merged CSR
  ! positions so that subsequent k-points can skip the O(NNZ log NNZ)
  ! sort and just update values in O(NNZ) time.
  !
  ! coo_to_csr(i) gives the CSR position (1..nnz_final) that COO entry i
  ! maps to after sort and duplicate merging.  For duplicate COO entries
  ! that merge into the same CSR position, all their indices map to the
  ! same CSR position.
  !
  ! The nnz of each kp-term block is kz-independent, so this mapping is
  ! fixed across all k-points in a kz sweep.
  ! ------------------------------------------------------------------
  type :: wire_coo_cache
    integer, allocatable          :: coo_to_csr(:)   ! (coo_nnz_in) mapping
    integer                       :: coo_nnz_in = 0  ! input COO count
    logical                       :: initialized = .false.
  end type wire_coo_cache

  public :: wire_coo_cache, wire_coo_cache_free

  contains

    subroutine build_kpterm_block(kpterms, profile_vec, central, forward, &
      & backward, diag, offup, offdown, N, term_idx, scale_factor, has_diag)

      real(kind=dp), intent(inout), dimension(:,:,:) :: kpterms
      real(kind=dp), intent(in), dimension(:) :: profile_vec
      real(kind=dp), intent(in), dimension(:,:) :: central, forward, backward
      real(kind=dp), intent(inout), dimension(:) :: diag, offup, offdown
      integer, intent(in) :: N, term_idx
      real(kind=dp), intent(in) :: scale_factor
      logical, intent(in) :: has_diag

      integer :: ii

      if (has_diag) then
        call dgemv('N', N, N, 1.0_dp, central, N, profile_vec, 1, 0.0_dp, diag, 1)
      end if
      call dgemv('N', N, N, 1.0_dp, forward, N, profile_vec, 1, 0.0_dp, offup, 1)
      call dgemv('N', N, N, 1.0_dp, backward, N, profile_vec, 1, 0.0_dp, offdown, 1)

      if (has_diag) then
        forall(ii=2:N-1)
          kpterms(ii,ii,term_idx) = diag(ii)
          kpterms(ii+1,ii,term_idx) = -offup(ii)
          kpterms(ii-1,ii,term_idx) = -offdown(ii)
        end forall
        kpterms(1,1,term_idx) = diag(1)
        kpterms(N,N,term_idx) = diag(N)
        kpterms(2,1,term_idx) = -offup(1)
        kpterms(N-1,N,term_idx) = -offdown(N)
      else
        forall(ii=2:N-1)
          kpterms(ii+1,ii,term_idx) = offup(ii)
          kpterms(ii-1,ii,term_idx) = -offdown(ii)
        end forall
        kpterms(2,1,term_idx) = offup(1)
        kpterms(N-1,N,term_idx) = -offdown(N)
      end if

      kpterms(:,:,term_idx) = kpterms(:,:,term_idx) * scale_factor

    end subroutine build_kpterm_block

    subroutine externalFieldSetup_electricField(profile, Evalue, totalSize, z)

      real(kind = dp), intent(inout), allocatable, dimension(:,:) :: profile
      real(kind = dp), intent(in) :: Evalue, totalSize
      real(kind = dp), intent(in), dimension(:) :: z

      integer :: i

      do i = 1, ubound(z, dim=1), 1
          profile(i,:) = profile(i,:) - (Evalue*totalSize) * (z(i)+z(1))/(2.0_dp*z(1))
      end do

    end subroutine externalFieldSetup_electricField

    subroutine confinementInitialization_raw(z, startPos, endPos, material, nlayers,&
      & params, confDir, profile, kpterms, FDorder)

      real(kind = dp), intent(in), dimension(:) :: z
      integer, intent(in), dimension(:) :: startPos, endPos
      character(len = 255), intent(in) :: material(nlayers)
      integer, intent(in) :: nlayers
      type(paramStruct), intent(in) :: params(nlayers)
      character(len = 1), intent(in) :: confDir
      real(kind = dp), intent(inout), allocatable, dimension(:,:) :: profile
      real(kind = dp), intent(inout), dimension(:,:,:) :: kpterms
      integer, intent(in), optional :: FDorder

      real(kind=dp), allocatable, dimension(:,:) :: ScnDer, FstDer, kptermsProfile
      real(kind=dp), allocatable, dimension(:,:) :: forward, central, backward
      real(kind=dp), allocatable, dimension(:) :: diag, offup, offdown
      real(kind=dp), allocatable, dimension(:,:) :: D2, D1

      integer :: i, initIDX, endIDX, N, ii, jj
      integer :: order
      real(kind = dp) :: delta


      N = size(z, dim=1)
      delta = abs(z(2) - z(1))

      ! Resolve FD order (default 2 for backward compatibility)
      if (present(FDorder)) then
        order = FDorder
      else
        order = 2
      end if

      if (allocated(profile)) deallocate(profile)
      allocate(profile(N,3))

      if (allocated(kptermsProfile)) deallocate(kptermsProfile)
      allocate(kptermsProfile(N,5))


      if (confDir == 'z') then

        do i = 1, nlayers, 1
          if (params(i)%EV == 0.0_dp .and. params(i)%EC == 0.0_dp) then
            print *, "ERROR: Material '", trim(material(i)), &
              & "' has EV=0 and EC=0. Band offsets are required."
            print *, "  Check that the material name is correct and EV/EC are in the database."
            stop 1
          end if
        end do

        do i = 1, nlayers, 1


          profile(startPos(i):endPos(i),1) = params(i)%EV
          profile(startPos(i):endPos(i),2) = params(i)%EV - params(i)%DeltaSO
          profile(startPos(i):endPos(i),3) =  params(i)%EC

          kptermsProfile(startPos(i):endPos(i),1) = params(i)%gamma1
          kptermsProfile(startPos(i):endPos(i),2) = params(i)%gamma2
          kptermsProfile(startPos(i):endPos(i),3) = params(i)%gamma3
          kptermsProfile(startPos(i):endPos(i),4) = params(i)%A
          kptermsProfile(startPos(i):endPos(i),5) = params(i)%P

        end do

      else

        stop "confinement not implemented for directions other than z"

      end if

      forall(ii=1:N)
        kpterms(ii,ii,1) = kptermsProfile(ii,1) !gamma1
        kpterms(ii,ii,2) = kptermsProfile(ii,2) !gamma2
        kpterms(ii,ii,3) = kptermsProfile(ii,3) !gamma3
        kpterms(ii,ii,4) = kptermsProfile(ii,5) !P
        kpterms(ii,ii,10) = kptermsProfile(ii,4) !A
      end forall

      if (order == 2) then
        ! ---- Order 2: use existing tridiagonal approach (backward compatible) ----
        allocate(forward(N,N))
        allocate(backward(N,N))
        allocate(central(N,N))
        forward = 0.0_dp
        backward = 0.0_dp
        central = 0.0_dp

        forall(ii=1:N-1)
          forward(ii,ii) = 1
          forward(ii,ii+1) = 1
        end forall
        !last element done alone to avoid if inside loop
        forward(N,N) = 1
        backward = transpose(forward)
        central = backward + forward


        allocate(diag(N))
        allocate(offup(N))
        allocate(offdown(N))

        ! A*kz**2
        call build_kpterm_block(kpterms, kptermsProfile(1:N,4), central, forward, &
          & backward, diag, offup, offdown, N, 5, 1.0_dp/(2.0_dp*delta**2), .True.)


        !Q
        call build_kpterm_block(kpterms, kptermsProfile(1:N,1) - 2.0_dp*kptermsProfile(1:N,2), &
          & central, forward, backward, diag, offup, offdown, N, 7, 1.0_dp/(2.0_dp*delta**2), .True.)


        !T
        call build_kpterm_block(kpterms, kptermsProfile(1:N,1) + 2.0_dp*kptermsProfile(1:N,2), &
          & central, forward, backward, diag, offup, offdown, N, 8, 1.0_dp/(2.0_dp*delta**2), .True.)


        ! P*kz
        call build_kpterm_block(kpterms, kptermsProfile(1:N,5), central, forward, &
          & backward, diag, offup, offdown, N, 6, 1.0_dp/(4.0_dp*delta), .False.)


        ! S -> gamma3*kz
        call build_kpterm_block(kpterms, kptermsProfile(1:N,3), central, forward, &
          & backward, diag, offup, offdown, N, 9, 1.0_dp/(4.0_dp*delta), .False.)

      else
        ! ---- Higher order: use FD matrix approach ----
        ! Build 2nd-derivative and 1st-derivative FD matrices
        call buildFD2ndDerivMatrix(N, delta, order, D2)
        call buildFD1stDerivMatrix(N, delta, order, D1)

        ! A*kz**2 (term 5): profile = A(z), operator = d^2/dz^2
        call applyVariableCoeff(kpterms, kptermsProfile(1:N,4), D2, N, 5)

        ! Q (term 7): profile = gamma1 - 2*gamma2, operator = d^2/dz^2
        call applyVariableCoeff(kpterms, kptermsProfile(1:N,1) - 2.0_dp*kptermsProfile(1:N,2), &
          & D2, N, 7)

        ! T (term 8): profile = gamma1 + 2*gamma2, operator = d^2/dz^2
        call applyVariableCoeff(kpterms, kptermsProfile(1:N,1) + 2.0_dp*kptermsProfile(1:N,2), &
          & D2, N, 8)

        ! P*kz (term 6): profile = P(z), operator = d/dz
        call applyVariableCoeff(kpterms, kptermsProfile(1:N,5), D1, N, 6)

        ! S -> gamma3*kz (term 9): profile = gamma3(z), operator = d/dz
        call applyVariableCoeff(kpterms, kptermsProfile(1:N,3), D1, N, 9)

      end if

      if (allocated(kptermsProfile)) deallocate(kptermsProfile)
      if (allocated(forward)) deallocate(forward)
      if (allocated(backward)) deallocate(backward)
      if (allocated(central)) deallocate(central)
      if (allocated(diag)) deallocate(diag)
      if (allocated(offup)) deallocate(offup)
      if (allocated(offdown)) deallocate(offdown)
      if (allocated(D2)) deallocate(D2)
      if (allocated(D1)) deallocate(D1)

    end subroutine confinementInitialization_raw

    !---------------------------------------------------------------------------
    !> Cfg-based wrapper: extracts spatial_grid fields from simulation_config
    !> and delegates to confinementInitialization_raw.
    !>
    !> For QW (ndim=1), uses cfg%grid%z(:) for coordinates and
    !> grid_ngrid(cfg%grid) for the spatial DOF count.  Falls back to the
    !> legacy cfg fields (z, intStartPos, intEndPos) when the grid has not
    !> been initialised yet (grid%ndim == 0 and confinement == 1).
    !---------------------------------------------------------------------------
    subroutine confinementInitialization_cfg(cfg, profile, kpterms)

      type(simulation_config), intent(inout) :: cfg
      real(kind=dp), intent(inout), allocatable :: profile(:,:)
      real(kind=dp), intent(inout), allocatable :: kpterms(:,:,:)

      integer :: ny

      ! Determine spatial DOF count from the grid when available,
      ! otherwise from the legacy fdStep field.
      if (cfg%grid%ndim >= 1) then
        ny = grid_ngrid(cfg%grid)
      else
        ny = cfg%fdStep
      end if

      ! Allocate kpterms if not already allocated (caller may pre-allocate)
      if (.not. allocated(kpterms)) then
        allocate(kpterms(ny, ny, 10))
        kpterms = 0.0_dp
      end if

      ! Dispatch to the raw routine using the appropriate coordinate source
      if (cfg%grid%ndim >= 1 .and. allocated(cfg%grid%z)) then
        ! New path: use spatial_grid coordinates
        call confinementInitialization_raw(cfg%grid%z, cfg%intStartPos, &
          & cfg%intEndPos, cfg%materialN, cfg%numLayers, cfg%params, &
          & cfg%confDir, profile, kpterms, cfg%FDorder)
      else
        ! Legacy path: use cfg%z (backward compat before init_grid_from_config)
        call confinementInitialization_raw(cfg%z, cfg%intStartPos, &
          & cfg%intEndPos, cfg%materialN, cfg%numLayers, cfg%params, &
          & cfg%confDir, profile, kpterms, cfg%FDorder)
      end if

      ! Populate cfg%dz from the grid when available
      if (cfg%grid%ndim >= 1 .and. cfg%grid%dy > 0.0_dp) then
        cfg%dz = cfg%grid%dy
      end if

    end subroutine confinementInitialization_cfg

    !---------------------------------------------------------------------------
    !> 2D confinement initialization for quantum wire mode (ndim=2).
    !>
    !> Builds kpterms_2d (11 CSR matrices) from 1D FD operators and
    !> Kronecker products, then applies position-dependent material
    !> parameters from the spatial_grid.
    !>
    !> kpterms_2d indices:
    !>   1: gamma1 (diagonal only)
    !>   2: gamma2 (diagonal only)
    !>   3: gamma3 (diagonal only)
    !>   4: P      (diagonal only)
    !>   5: A * (d^2/dx^2 + d^2/dy^2)   -- 2nd derivative, conduction band kinetic
    !>   6: P * (d/dx + d/dy)            -- 1st derivative, momentum coupling
    !>   7: (gamma1 - 2*gamma2) * Laplacian -- Q term kinetic
    !>   8: (gamma1 + 2*gamma2) * Laplacian -- T term kinetic
    !>   9: gamma3 * (d/dx + d/dy)       -- S term, linear-k coupling
    !>  10: A (diagonal only, k^2 coefficient for CB)
    !>  11: gamma3 * d^2/dxdy             -- cross-derivative (NEW for 2D)
    !>  12: P * d/dx   (x-gradient only, for g-factor perturbation)
    !>  13: P * d/dy   (y-gradient only, for g-factor perturbation)
    !>  14: gamma3 * d/dx  (x-gradient only, for g-factor perturbation)
    !>  15: gamma3 * d/dy  (y-gradient only, for g-factor perturbation)
    !>
    !> Also builds profile_2d(Ngrid, 3) with band edges:
    !>   profile_2d(:,1) = EV
    !>   profile_2d(:,2) = EV - DeltaSO
    !>   profile_2d(:,3) = EC
    !---------------------------------------------------------------------------
    subroutine confinementInitialization_2d(grid, params, regions, &
        profile_2d, kpterms_2d, FDorder)

      type(spatial_grid), intent(in)    :: grid
      type(paramStruct), intent(in)     :: params(:)     ! size = numRegions
      type(region_spec), intent(in)     :: regions(:)    ! size = numRegions
      real(kind=dp), allocatable, intent(out) :: profile_2d(:,:)
      type(csr_matrix), allocatable, intent(out) :: kpterms_2d(:)
      integer, intent(in), optional     :: FDorder

      integer :: order, nx, ny, ngrid, ij, mid
      real(kind=dp) :: dx, dy

      ! 1D FD matrices (real, dense)
      real(kind=dp), allocatable :: D2x(:,:), D1x(:,:)
      real(kind=dp), allocatable :: D2y(:,:), D1y(:,:)
      real(kind=dp), allocatable :: Ix(:,:), Iy(:,:)

      ! Complex versions for Kron product routines
      complex(kind=dp), allocatable :: cD2x(:,:), cD1x(:,:)
      complex(kind=dp), allocatable :: cD2y(:,:), cD1y(:,:)
      complex(kind=dp), allocatable :: cIx(:,:), cIy(:,:)

      ! CSR work matrices
      type(csr_matrix) :: kron_Iy_D2x, kron_D2y_Ix, laplacian_2d
      type(csr_matrix) :: kron_Iy_D1x, kron_D1y_Ix, grad_2d
      type(csr_matrix) :: kron_D1y_D1x  ! cross-derivative
      type(csr_matrix) :: diag_csr, scaled_diag

      ! Material profiles on the 2D grid
      real(kind=dp), allocatable :: prof_gamma1(:), prof_gamma2(:)
      real(kind=dp), allocatable :: prof_gamma3(:), prof_P(:), prof_A(:)
      real(kind=dp), allocatable :: prof_gm12g2(:), prof_gp12g2(:)

      ! Resolve FD order (default 2)
      if (present(FDorder)) then
        order = FDorder
      else
        order = 2
      end if

      nx    = grid%nx
      ny    = grid%ny
      ngrid = nx * ny
      dx    = grid%dx
      dy    = grid%dy

      ! Validate grid
      if (nx < 3 .or. ny < 3) then
        print *, "ERROR: confinementInitialization_2d requires nx>=3, ny>=3"
        stop 1
      end if

      ! ====================================================================
      ! 1. Build 1D FD operators
      ! ====================================================================
      call buildFD2ndDerivMatrix(nx, dx, order, D2x)
      call buildFD1stDerivMatrix(nx, dx, order, D1x)
      call buildFD2ndDerivMatrix(ny, dy, order, D2y)
      call buildFD1stDerivMatrix(ny, dy, order, D1y)

      ! Identity matrices
      call Identity(nx, Ix)
      call Identity(ny, Iy)

      ! Convert real matrices to complex for Kron routines
      allocate(cD2x(nx, nx)); cD2x = cmplx(D2x, 0.0_dp, kind=dp)
      allocate(cD1x(nx, nx)); cD1x = cmplx(D1x, 0.0_dp, kind=dp)
      allocate(cD2y(ny, ny)); cD2y = cmplx(D2y, 0.0_dp, kind=dp)
      allocate(cD1y(ny, ny)); cD1y = cmplx(D1y, 0.0_dp, kind=dp)
      allocate(cIx(nx, nx));   cIx  = cmplx(Ix,  0.0_dp, kind=dp)
      allocate(cIy(ny, ny));   cIy  = cmplx(Iy,  0.0_dp, kind=dp)

      ! ====================================================================
      ! 2. Build 2D FD operators via Kronecker products
      ! ====================================================================

      ! Laplacian: Iy x D2x + D2y x Ix  (column-major flat_idx = (iy-1)*nx+ix)
      call kron_eye_dense(ny, cD2x, nx, nx, kron_Iy_D2x)
      call kron_dense_dense(cD2y, ny, ny, cIx, nx, nx, kron_D2y_Ix)
      call csr_add(kron_Iy_D2x, kron_D2y_Ix, laplacian_2d)

      ! Gradient: Iy x D1x + D1y x Ix  (column-major)
      call kron_eye_dense(ny, cD1x, nx, nx, kron_Iy_D1x)
      call kron_dense_dense(cD1y, ny, ny, cIx, nx, nx, kron_D1y_Ix)
      call csr_add(kron_Iy_D1x, kron_D1y_Ix, grad_2d)

      ! Cross-derivative: D1y x D1x  (column-major)
      call kron_dense_dense(cD1y, ny, ny, cD1x, nx, nx, kron_D1y_D1x)

      ! Free intermediate Kronecker products (keep kron_Iy_D1x, kron_D1y_Ix
      ! for building kpterms_2d(12-15) below)
      call csr_free(kron_Iy_D2x)
      call csr_free(kron_D2y_Ix)

      ! ====================================================================
      ! 3. Build material profiles on the 2D grid
      ! ====================================================================
      allocate(prof_gamma1(ngrid)); prof_gamma1 = 0.0_dp
      allocate(prof_gamma2(ngrid)); prof_gamma2 = 0.0_dp
      allocate(prof_gamma3(ngrid)); prof_gamma3 = 0.0_dp
      allocate(prof_P(ngrid));      prof_P      = 0.0_dp
      allocate(prof_A(ngrid));      prof_A      = 0.0_dp
      allocate(prof_gm12g2(ngrid)); prof_gm12g2 = 0.0_dp
      allocate(prof_gp12g2(ngrid)); prof_gp12g2 = 0.0_dp

      if (allocated(profile_2d)) deallocate(profile_2d)
      allocate(profile_2d(ngrid, 3))
      profile_2d = 0.0_dp

      do ij = 1, ngrid
        mid = grid%material_id(ij)
        if (mid < 1 .or. mid > size(params)) then
          ! Inactive cell or outside domain -- leave zeros
          cycle
        end if

        if (params(mid)%EV == 0.0_dp .and. params(mid)%EC == 0.0_dp) then
          print *, "ERROR: Region '", trim(regions(mid)%material), &
            "' has EV=0 and EC=0. Band offsets are required."
          stop 1
        end if

        ! Band-edge profile
        profile_2d(ij, 1) = params(mid)%EV
        profile_2d(ij, 2) = params(mid)%EV - params(mid)%DeltaSO
        profile_2d(ij, 3) = params(mid)%EC

        ! Material parameter profiles
        prof_gamma1(ij) = params(mid)%gamma1
        prof_gamma2(ij) = params(mid)%gamma2
        prof_gamma3(ij) = params(mid)%gamma3
        prof_P(ij)      = params(mid)%P
        prof_A(ij)      = params(mid)%A
        prof_gm12g2(ij) = params(mid)%gamma1 - 2.0_dp * params(mid)%gamma2
        prof_gp12g2(ij) = params(mid)%gamma1 + 2.0_dp * params(mid)%gamma2
      end do

      ! ====================================================================
      ! 4. Build kpterms_2d (15 CSR matrices)
      ! ====================================================================
      if (allocated(kpterms_2d)) deallocate(kpterms_2d)
      allocate(kpterms_2d(15))

      ! Terms 1-4: diagonal-only (gamma1, gamma2, gamma3, P)
      ! Build a diagonal CSR from the profile and store directly
      call build_diagonal_csr(ngrid, prof_gamma1, kpterms_2d(1))
      call build_diagonal_csr(ngrid, prof_gamma2, kpterms_2d(2))
      call build_diagonal_csr(ngrid, prof_gamma3, kpterms_2d(3))
      call build_diagonal_csr(ngrid, prof_P, kpterms_2d(4))

      ! Apply cut-cell face fractions to differential operators.
      ! For rectangular grids all face fractions are 1.0 (no-op).
      ! Terms 5,7,8: Laplacian;  Terms 6,9: Gradient.
      ! Term 11 (cross-derivative) mixes x/y connections so face fractions
      ! are not applied (they remain 1.0 by default).
      if (allocated(grid%face_fraction_x) .and. allocated(grid%face_fraction_y)) then
        ! Term 5: A * Laplacian
        call csr_apply_variable_coeff(laplacian_2d, prof_A, kpterms_2d(5), &
          face_frac_x=grid%face_fraction_x, &
          face_frac_y=grid%face_fraction_y, grid_nx=nx)
        ! Term 6: P * Gradient
        call csr_apply_variable_coeff(grad_2d, prof_P, kpterms_2d(6), &
          face_frac_x=grid%face_fraction_x, &
          face_frac_y=grid%face_fraction_y, grid_nx=nx)
        ! Term 7: (gamma1 - 2*gamma2) * Laplacian
        call csr_apply_variable_coeff(laplacian_2d, prof_gm12g2, kpterms_2d(7), &
          face_frac_x=grid%face_fraction_x, &
          face_frac_y=grid%face_fraction_y, grid_nx=nx)
        ! Term 8: (gamma1 + 2*gamma2) * Laplacian
        call csr_apply_variable_coeff(laplacian_2d, prof_gp12g2, kpterms_2d(8), &
          face_frac_x=grid%face_fraction_x, &
          face_frac_y=grid%face_fraction_y, grid_nx=nx)
        ! Term 9: gamma3 * Gradient
        call csr_apply_variable_coeff(grad_2d, prof_gamma3, kpterms_2d(9), &
          face_frac_x=grid%face_fraction_x, &
          face_frac_y=grid%face_fraction_y, grid_nx=nx)
      else
        ! No face fractions (rectangular grid or QW mode): standard path
        call csr_apply_variable_coeff(laplacian_2d, prof_A, kpterms_2d(5))
        call csr_apply_variable_coeff(grad_2d, prof_P, kpterms_2d(6))
        call csr_apply_variable_coeff(laplacian_2d, prof_gm12g2, kpterms_2d(7))
        call csr_apply_variable_coeff(laplacian_2d, prof_gp12g2, kpterms_2d(8))
        call csr_apply_variable_coeff(grad_2d, prof_gamma3, kpterms_2d(9))
      end if

      ! Term 10: A (diagonal only)
      call build_diagonal_csr(ngrid, prof_A, kpterms_2d(10))

      ! Term 11: gamma3 * Cross-derivative
      call csr_apply_variable_coeff(kron_D1y_D1x, prof_gamma3, kpterms_2d(11))

      ! Terms 12-15: separate x/y gradient operators for g-factor perturbation.
      ! These use the same kron products as grad_2d but split by direction.
      ! kron_Iy_D1x = Iy x D1x (x-gradient), kron_D1y_Ix = D1y x Ix (y-gradient).
      if (allocated(grid%face_fraction_x) .and. allocated(grid%face_fraction_y)) then
        ! Term 12: P * d/dx (x-gradient only)
        call csr_apply_variable_coeff(kron_Iy_D1x, prof_P, kpterms_2d(12), &
          face_frac_x=grid%face_fraction_x, &
          face_frac_y=grid%face_fraction_y, grid_nx=nx)
        ! Term 13: P * d/dy (y-gradient only)
        call csr_apply_variable_coeff(kron_D1y_Ix, prof_P, kpterms_2d(13), &
          face_frac_x=grid%face_fraction_x, &
          face_frac_y=grid%face_fraction_y, grid_nx=nx)
        ! Term 14: gamma3 * d/dx (x-gradient only)
        call csr_apply_variable_coeff(kron_Iy_D1x, prof_gamma3, kpterms_2d(14), &
          face_frac_x=grid%face_fraction_x, &
          face_frac_y=grid%face_fraction_y, grid_nx=nx)
        ! Term 15: gamma3 * d/dy (y-gradient only)
        call csr_apply_variable_coeff(kron_D1y_Ix, prof_gamma3, kpterms_2d(15), &
          face_frac_x=grid%face_fraction_x, &
          face_frac_y=grid%face_fraction_y, grid_nx=nx)
      else
        ! No face fractions (rectangular grid)
        call csr_apply_variable_coeff(kron_Iy_D1x, prof_P, kpterms_2d(12))
        call csr_apply_variable_coeff(kron_D1y_Ix, prof_P, kpterms_2d(13))
        call csr_apply_variable_coeff(kron_Iy_D1x, prof_gamma3, kpterms_2d(14))
        call csr_apply_variable_coeff(kron_D1y_Ix, prof_gamma3, kpterms_2d(15))
      end if

      ! ====================================================================
      ! 5. Cleanup
      ! ====================================================================
      call csr_free(laplacian_2d)
      call csr_free(grad_2d)
      call csr_free(kron_D1y_D1x)
      call csr_free(kron_Iy_D1x)
      call csr_free(kron_D1y_Ix)

      deallocate(D2x, D1x, D2y, D1y, Ix, Iy)
      deallocate(cD2x, cD1x, cD2y, cD1y, cIx, cIy)
      deallocate(prof_gamma1, prof_gamma2, prof_gamma3)
      deallocate(prof_P, prof_A, prof_gm12g2, prof_gp12g2)

    end subroutine confinementInitialization_2d

    !---------------------------------------------------------------------------
    !> Build a diagonal CSR matrix from a real profile vector.
    !> mat(i,i) = cmplx(profile(i), 0.0_dp).
    !---------------------------------------------------------------------------
    subroutine build_diagonal_csr(n, profile, mat)
      integer, intent(in) :: n
      real(kind=dp), intent(in) :: profile(n)
      type(csr_matrix), intent(out) :: mat

      integer :: i
      integer, allocatable :: rows(:), cols(:)
      complex(kind=dp), allocatable :: vals(:)

      allocate(rows(n), cols(n), vals(n))
      do i = 1, n
        rows(i) = i
        cols(i) = i
        vals(i) = cmplx(profile(i), 0.0_dp, kind=dp)
      end do

      call csr_build_from_coo(mat, n, n, n, rows, cols, vals)
      deallocate(rows, cols, vals)
    end subroutine build_diagonal_csr

    !---------------------------------------------------------------------------
    !> Apply variable coefficient to FD matrix using midpoint averaging.
    !> Computes: kpterms(:,:,term_idx) = -G_avg .* FD
    !> where G_avg(i,i) = profile(i) and G_avg(i,j) = (profile(i)+profile(j))/2
    !> This correctly discretizes d/dz[g(z)·d/dz] at heterointerfaces.
    !---------------------------------------------------------------------------
    subroutine applyVariableCoeff(kpterms, profile_vec, FD, N, term_idx)

      real(kind=dp), intent(inout), dimension(:,:,:) :: kpterms
      real(kind=dp), intent(in), dimension(:) :: profile_vec
      real(kind=dp), intent(in), dimension(:,:) :: FD
      integer, intent(in) :: N, term_idx

      integer :: ii, jj
      real(kind=dp) :: g_avg

      do ii = 1, N
        do jj = 1, N
          if (FD(jj, ii) == 0.0_dp) cycle
          if (jj == ii) then
            ! Diagonal: use local profile value
            kpterms(jj, ii, term_idx) = -profile_vec(jj) * FD(jj, ii)
          else
            ! Off-diagonal: midpoint average
            g_avg = 0.5_dp * (profile_vec(jj) + profile_vec(ii))
            kpterms(jj, ii, term_idx) = -g_avg * FD(jj, ii)
          end if
        end do
      end do

    end subroutine applyVariableCoeff

    subroutine ZB8bandQW(HT, wv, profile, kpterms, cfg, sparse, HT_csr, g)

      implicit none

      !input/output
      complex(kind=dp), intent(inout), dimension(:,:) :: HT
      type(csr_matrix), intent(inout), optional :: HT_csr
      type(wavevector), intent(in) :: wv
      real(kind = dp), intent(in), dimension(:,:) :: profile
      real(kind = dp), intent(in), dimension(:,:,:) :: kpterms
      type(simulation_config), intent(in), optional :: cfg
      logical, intent(in), optional:: sparse
      character(len=1), intent(in), optional :: g

      integer :: nzmax

      ! wave vector related
      real(kind=dp) :: kx2, ky2, kxky, k2
      complex(kind=dp) :: kminus, kplus

      complex(kind=dp), allocatable, dimension(:,:) :: Q, R, RC, S, SC, T, PZ, PM, PP, A

      ! constants

      integer :: i, j, N, ii, jj

      N = size(HT, dim=1)/8

      !constants


      !wave vectors
      kxky = wv%kx*wv%ky
      kx2 = wv%kx**2
      ky2 = wv%ky**2
      k2 = kx2 + ky2
      kplus   = wv%kx + IU*wv%ky
      kminus  = wv%kx - IU*wv%ky

      allocate(Q(N,N))
      allocate(T(N,N))
      allocate(S(N,N))
      allocate(SC(N,N))
      allocate(R(N,N))
      allocate(RC(N,N))
      allocate(PP(N,N))
      allocate(PM(N,N))
      allocate(PZ(N,N))
      allocate(A(N,N))
      Q = 0.0_dp
      T = 0.0_dp
      S = 0.0_dp
      SC = 0.0_dp
      R = 0.0_dp
      RC = 0.0_dp
      PP = 0.0_dp
      PM = 0.0_dp
      PZ = 0.0_dp
      A = 0.0_dp

      if (.not. present(g)) then
        forall(ii=1:N, jj=1:N)
          Q(ii,jj) = -( (kpterms(ii,jj,1) + kpterms(ii,jj,2))*k2 + kpterms(ii,jj,7) )
          T(ii,jj) = -( (kpterms(ii,jj,1) - kpterms(ii,jj,2))*k2 + kpterms(ii,jj,8) )
          S(ii,jj)  =  2.0_dp * SQR3 * kminus * kpterms(ii,jj,9) !no IU because of derivative
          SC(jj,ii) =  2.0_dp * SQR3 * kplus  * kpterms(ii,jj,9) !no IU because of derivative
          PZ(ii,jj) = kpterms(ii,jj,6) * (-IU)
          A(ii,jj)  = kpterms(ii,jj,5) + k2*kpterms(ii,jj,10)
        end forall

        forall (ii=1:N)
          R(ii,ii)  = - SQR3 * ( kpterms(ii,ii,2)*(kx2 - ky2) - 2.0_dp*IU*kpterms(ii,ii,3)*kxky )
          RC(ii,ii) = - SQR3 * ( kpterms(ii,ii,2)*(kx2 - ky2) + 2.0_dp*IU*kpterms(ii,ii,3)*kxky )
          PP(ii,ii) = kpterms(ii,ii,4) * kplus  * RQS2
          PM(ii,ii) = kpterms(ii,ii,4) * kminus * RQS2
        end forall
      else
        forall(ii=1:N, jj=1:N)
          PZ(ii,jj) = kpterms(ii,jj,4) * wv%kz
        end forall
        forall (ii=1:N)
          PP(ii,ii) = kpterms(ii,ii,4) * kplus  * RQS2
          PM(ii,ii) = kpterms(ii,ii,4) * kminus * RQS2
        end forall

      end if


      HT = 0

      ! kp matrix
      !col 1
      HT(1 + 0*N : 1*N, 1 + 0*N : 1*N) =  Q
      HT(1 + 0*N : 1*N, 1 + 1*N : 2*N) =  SC
      HT(1 + 0*N : 1*N, 1 + 2*N : 3*N) =  RC
      ! HT(1 + 0*N : 1*N, 1 + 3*N : 4*N) =  0.0_dp
      HT(1 + 0*N : 1*N, 1 + 4*N : 5*N) = -IU * RQS2 * SC
      HT(1 + 0*N : 1*N, 1 + 5*N : 6*N) =  IU * SQR2 * RC
      HT(1 + 0*N : 1*N, 1 + 6*N : 7*N) =  IU * PP
      ! HT(1 + 0*N : 1*N, 1 + 7*N : 8*N) =  0.0_dp

      !col 2
      HT(1 + 1*N : 2*N, 1 + 0*N : 1*N) =  (S)
      HT(1 + 1*N : 2*N, 1 + 1*N : 2*N) =  T
      ! HT(1 + 1*N : 2*N, 1 + 2*N : 3*N) =  0.0_dp
      HT(1 + 1*N : 2*N, 1 + 3*N : 4*N) =  RC
      HT(1 + 1*N : 2*N, 1 + 4*N : 5*N) =  IU * RQS2 * (Q - T)
      HT(1 + 1*N : 2*N, 1 + 5*N : 6*N) = -IU * SQR3 * RQS2 * SC
      HT(1 + 1*N : 2*N, 1 + 6*N : 7*N) =  SQR2 * RQS3 * PZ
      HT(1 + 1*N : 2*N, 1 + 7*N : 8*N) = -RQS3*PP

      !col 3
      HT(1 + 2*N : 3*N, 1 + 0*N : 1*N) =  R
      HT(1 + 2*N : 3*N, 1 + 2*N : 3*N) =  T
      HT(1 + 2*N : 3*N, 1 + 3*N : 4*N) = -SC
      HT(1 + 2*N : 3*N, 1 + 4*N : 5*N) =  IU * SQR3 * RQS2 * S
      HT(1 + 2*N : 3*N, 1 + 5*N : 6*N) =  IU * RQS2 * (Q - T)
      HT(1 + 2*N : 3*N, 1 + 6*N : 7*N) =  IU * RQS3 * PM
      HT(1 + 2*N : 3*N, 1 + 7*N : 8*N) =  IU * SQR2 * RQS3 * PZ

      !col 4
      HT(1 + 3*N : 4*N, 1 + 1*N : 2*N) =  R
      HT(1 + 3*N : 4*N, 1 + 2*N : 3*N) = (-S)
      HT(1 + 3*N : 4*N, 1 + 3*N : 4*N) =  Q
      HT(1 + 3*N : 4*N, 1 + 4*N : 5*N) =  IU * SQR2 * R
      HT(1 + 3*N : 4*N, 1 + 5*N : 6*N) =  IU * RQS2 * S
      ! HT(1 + 3*N : 4*N, 1 + 6*N : 7*N) =  0.0_dp
      HT(1 + 3*N : 4*N, 1 + 7*N : 8*N) = -PM

      !col 5
      HT(1 + 4*N : 5*N, 1 + 0*N : 1*N) =  IU * RQS2 * (S)
      HT(1 + 4*N : 5*N, 1 + 1*N : 2*N) = -IU * RQS2 * (Q - T)
      HT(1 + 4*N : 5*N, 1 + 2*N : 3*N) = -IU * SQR3 * RQS2 * (SC)
      HT(1 + 4*N : 5*N, 1 + 3*N : 4*N) = -IU * SQR2 * RC
      HT(1 + 4*N : 5*N, 1 + 4*N : 5*N) =  0.5_dp*(Q + T)
      ! HT(1 + 4*N : 5*N, 1 + 5*N : 6*N) =  0.0_dp
      HT(1 + 4*N : 5*N, 1 + 6*N : 7*N) =  IU * RQS3 * PZ
      HT(1 + 4*N : 5*N, 1 + 7*N : 8*N) =  IU * SQR2 * RQS3 * PP

      !col 6
      HT(1 + 5*N : 6*N, 1 + 0*N : 1*N) = -IU * SQR2 * R
      HT(1 + 5*N : 6*N, 1 + 1*N : 2*N) =  IU * SQR3 * RQS2 * (S)
      HT(1 + 5*N : 6*N, 1 + 2*N : 3*N) = -IU * RQS2 * (Q- T)
      HT(1 + 5*N : 6*N, 1 + 3*N : 4*N) = -IU * RQS2 * (SC)
      HT(1 + 5*N : 6*N, 1 + 5*N : 6*N) =  0.5_dp*(Q + T)
      HT(1 + 5*N : 6*N, 1 + 6*N : 7*N) =  SQR2 * RQS3 * PM
      HT(1 + 5*N : 6*N, 1 + 7*N : 8*N) = -RQS3 * PZ

      !col 7
      HT(1 + 6*N : 7*N, 1 + 0*N : 1*N) = -IU * PM
      HT(1 + 6*N : 7*N, 1 + 1*N : 2*N) =  SQR2 * RQS3 * (PZ)
      HT(1 + 6*N : 7*N, 1 + 2*N : 3*N) = -IU * RQS3 * PP
      HT(1 + 6*N : 7*N, 1 + 4*N : 5*N) = -IU * RQS3 * (PZ)
      HT(1 + 6*N : 7*N, 1 + 5*N : 6*N) =  SQR2 * RQS3 * PP
      HT(1 + 6*N : 7*N, 1 + 6*N : 7*N) =  A
      ! HT(7,8) =  0.0_dp

      !col 8
      HT(1 + 7*N : 8*N, 1 + 1*N : 2*N) = -RQS3 * PM
      HT(1 + 7*N : 8*N, 1 + 2*N : 3*N) = -IU * SQR2 * RQS3 * (PZ)
      HT(1 + 7*N : 8*N, 1 + 3*N : 4*N) = -PP
      HT(1 + 7*N : 8*N, 1 + 4*N : 5*N) = -IU * SQR2 * RQS3 * PM
      HT(1 + 7*N : 8*N, 1 + 5*N : 6*N) = -RQS3 * (PZ)
      HT(1 + 7*N : 8*N, 1 + 7*N : 8*N) =  A


      ! Profile
      forall(ii=1:N)
        HT(      ii,      ii) = HT(      ii,      ii) + profile(ii,1)
        HT(  N + ii,  N + ii) = HT(  N + ii,  N + ii) + profile(ii,1)
        HT(2*N + ii,2*N + ii) = HT(2*N + ii,2*N + ii) + profile(ii,1)
        HT(3*N + ii,3*N + ii) = HT(3*N + ii,3*N + ii) + profile(ii,1)
        HT(4*N + ii,4*N + ii) = HT(4*N + ii,4*N + ii) + profile(ii,2)
        HT(5*N + ii,5*N + ii) = HT(5*N + ii,5*N + ii) + profile(ii,2)
        HT(6*N + ii,6*N + ii) = HT(6*N + ii,6*N + ii) + profile(ii,3)
        HT(7*N + ii,7*N + ii) = HT(7*N + ii,7*N + ii) + profile(ii,3)
      end forall

      ! Full Bir-Pikus strain (k-independent, not in g-mode)
      if (present(cfg)) then
        if (.not. present(g) .and. allocated(cfg%strain_blocks%delta_Ec)) then
          do ii = 1, N
            call add_bp_strain_dense(HT, ii, N, cfg%strain_blocks)
          end do
        end if
      end if


      deallocate(Q)
      deallocate(T)
      deallocate(S)
      deallocate(SC)
      deallocate(R)
      deallocate(RC)
      deallocate(PP)
      deallocate(PM)
      deallocate(PZ)
      deallocate(A)

      if (present(sparse)) then
        if (sparse .eqv. .True.) then
          nzmax = (N + (N-1)*2 + 2)*64

          call dnscsr_z_mkl(nzmax, N*8, HT, HT_csr)

        end if
      end if


    end subroutine ZB8bandQW


    subroutine ZB8bandBulk(HT,wv,params,g)

      implicit none

      !input/output
      complex(kind=dp), intent(inout), dimension(:,:) :: HT
      type(wavevector), intent(in) :: wv
      type(paramStruct), intent(in) :: params(1)
      character(len=1), intent(in), optional :: g

      ! wave vector related
      real(kind=dp) :: kx2, ky2, kz2, kxky, k2
      complex(kind=dp) :: kminusz, kplusz, kminus, kplus

      !kp parameters
      real(kind=dp) :: gamma1, gamma2, gamma3, P, A
      complex(kind=dp) :: PM, PP, PZ
      complex(kind=dp) :: Q, R, RC, S, SC, T

      ! constants

      integer :: i, j


      !constants


      !wave vectors
      kx2 = wv%kx**2
      ky2 = wv%ky**2
      kz2 = wv%kz**2
      kxky = wv%kx*wv%ky
      k2 = kx2 + ky2 + kz2

      kminusz = (wv%kx - IU*wv%ky)*wv%kz
      kplusz  = (wv%kx + IU*wv%ky)*wv%kz
      kplus   = wv%kx + IU*wv%ky
      kminus  = wv%kx - IU*wv%ky

      !kp parameters
      if (present(g) .and. g == 'g') then
        gamma1 = 0!params(1)%gamma1
        gamma2 = 0!params(1)%gamma2
        gamma3 = 0!params(1)%gamma3
        P = params(1)%P
        A = 0!params(1)%A
      else
        gamma1 = params(1)%gamma1
        gamma2 = params(1)%gamma2
        gamma3 = params(1)%gamma3
        P = params(1)%P
        A = params(1)%A
      end if


      ! kp terms
      Q = -( (gamma1 + gamma2)*(kx2 + ky2) + (gamma1 - 2.0_dp*gamma2)*kz2 )
      T = -( (gamma1 - gamma2)*(kx2 + ky2) + (gamma1 + 2.0_dp*gamma2)*kz2 )

      S  =   IU * 2.0_dp * SQR3 * gamma3 * kminusz
      SC = - IU * 2.0_dp * SQR3 * gamma3 * kplusz

      R  = - SQR3 * ( gamma2*(kx2 - ky2) - 2.0_dp*IU*gamma3*kxky )
      RC = - SQR3 * ( gamma2*(kx2 - ky2) + 2.0_dp*IU*gamma3*kxky )

      PP = P * kplus  * RQS2
      PM = P * kminus * RQS2
      PZ = P * wv%kz


      HT = 0

      ! kp matrix
      !col 1
      HT(1,1) =  Q
      HT(1,2) =  SC
      HT(1,3) =  RC
      HT(1,4) =  0.0_dp
      HT(1,5) = -IU * RQS2 * SC
      HT(1,6) =  IU * SQR2 * RC
      HT(1,7) =  IU * PP
      HT(1,8) =  0.0_dp

      !col 2
      HT(2,1) =  S
      HT(2,2) =  T
      HT(2,3) =  0.0_dp
      HT(2,4) =  RC
      HT(2,5) =  IU * RQS2 * (Q - T)
      HT(2,6) = -IU * SQR3 * RQS2 * SC
      HT(2,7) =  SQR2 * RQS3 * PZ
      HT(2,8) = -RQS3*PP

      !col 3
      HT(3,1) =  R
      HT(3,3) =  T
      HT(3,4) = -SC
      HT(3,5) =  IU * SQR3 * RQS2 * S
      HT(3,6) =  IU * RQS2 * (Q - T)
      HT(3,7) =  IU * RQS3 * PM
      HT(3,8) =  IU * SQR2 * RQS3 * PZ

      !col 4
      HT(4,2) =  R
      HT(4,3) = -S
      HT(4,4) =  Q
      HT(4,5) =  IU * SQR2 * R
      HT(4,6) =  IU * RQS2 * S
      HT(4,7) =  0.0_dp
      HT(4,8) = -PM

      !col 5
      HT(5,1) =  IU * RQS2 * S
      HT(5,2) = -IU * RQS2 * (Q - T)
      HT(5,3) = -IU * SQR3 * RQS2 * SC
      HT(5,4) = -IU * SQR2 * RC
      HT(5,5) =  0.5_dp*(Q + T)
      HT(5,6) =  0.0_dp
      HT(5,7) =  IU * RQS3 * PZ
      HT(5,8) =  IU * SQR2 * RQS3 * PP

      !col 6
      HT(6,1) = -IU * SQR2 * R
      HT(6,2) =  IU * SQR3 * RQS2 * S
      HT(6,3) = -IU * RQS2 * (Q - T)
      HT(6,4) = -IU * RQS2 * SC
      HT(6,6) =  0.5_dp*(Q + T)
      HT(6,7) =  SQR2 * RQS3 * PM
      HT(6,8) = -RQS3 * PZ

      !col 7
      HT(7,1) = -IU * PM
      HT(7,2) =  SQR2 * RQS3 * PZ
      HT(7,3) = -IU * RQS3 * PP
      HT(7,5) = -IU * RQS3 * PZ
      HT(7,6) =  SQR2 * RQS3 * PP
      HT(7,7) =  A * K2
      HT(7,8) =  0.0_dp

      !col 8
      HT(8,2) = -RQS3 * PM
      HT(8,3) = -IU * SQR2 * RQS3 * PZ
      HT(8,4) = -PP
      HT(8,5) = -IU * SQR2 * RQS3 * PM
      HT(8,6) = -RQS3 * PZ
      HT(8,8) =  A * K2


      ! SOC
      HT(5,5) = HT(5,5) - params(1)%DeltaSO
      HT(6,6) = HT(6,6) - params(1)%DeltaSO

      HT(7,7) = HT(7,7) + params(1)%Eg
      HT(8,8) = HT(8,8) + params(1)%Eg

      ! ---------------------------------------------------------------
      ! Full Bir-Pikus strain Hamiltonian for bulk.
      !
      ! When strainSubstrate > 0, apply uniform biaxial [001] strain:
      !   eps_xx = eps_yy = eps_par = (a_sub - a_film) / a_film
      !   eps_zz = eps_perp = -2 C12/C11 * eps_par
      !   eps_xy = eps_xz = eps_yz = 0  (for [001] biaxial)
      !
      ! The strain Hamiltonian has the SAME matrix structure as the
      ! k-dependent terms but with eps_ij replacing k_i*k_j and
      ! deformation potentials (av, b_dp, d_dp) replacing Luttinger
      ! parameters (gamma1, gamma2, gamma3):
      !
      !   P_eps     = -av * Tr(eps)                           (hydrostatic)
      !   Q_eps     =  b_dp/2 * (eps_zz - 0.5*(eps_xx+eps_yy)) (tetragonal shear)
      !   R_eps     = -sqrt(3) * [b_dp/2*(eps_xx-eps_yy) - i*d_dp*eps_xy]
      !   S_eps     =  i*2*sqrt(3) * d_dp * (eps_xz - i*eps_yz)
      !   S_eps_bar = -i*2*sqrt(3) * d_dp * (eps_xz + i*eps_yz)
      !
      ! Diagonal:
      !   CB:  +ac * Tr(eps)
      !   HH:  +P_eps + Q_eps
      !   LH:  +P_eps - Q_eps
      !   SO:  +P_eps
      !
      ! Off-diagonal: same pattern as k-terms R, S, S_bar.
      ! For [001] biaxial: R_eps = S_eps = 0, only Q_eps survives.
      ! All terms included for physics consistency.
      ! ---------------------------------------------------------------
      if (params(1)%strainSubstrate > 0.0_dp) then
        block
          real(kind=dp) :: a_film, eps_xx, eps_zz
          type(bir_pikus_blocks) :: bp_bulk

          a_film = params(1)%a0
          if (a_film > 0.0_dp) then
            eps_xx = (params(1)%strainSubstrate - a_film) / a_film
            eps_zz = -2.0_dp * params(1)%C12 / params(1)%C11 * eps_xx

            allocate(bp_bulk%delta_Ec(1), bp_bulk%delta_EHH(1), &
              bp_bulk%delta_ELH(1), bp_bulk%delta_ESO(1), &
              bp_bulk%R_eps(1), bp_bulk%S_eps(1), bp_bulk%QT2_eps(1))

            associate(s => compute_bp_scalar(params(1), eps_xx, eps_xx, eps_zz, &
                                              0.0_dp, 0.0_dp, 0.0_dp))
              bp_bulk%delta_Ec(1)  = s%delta_Ec
              bp_bulk%delta_EHH(1) = s%delta_EHH
              bp_bulk%delta_ELH(1) = s%delta_ELH
              bp_bulk%delta_ESO(1) = s%delta_ESO
              bp_bulk%R_eps(1)     = s%R_eps
              bp_bulk%S_eps(1)     = s%S_eps
              bp_bulk%QT2_eps(1)   = s%QT2_eps
            end associate

            call add_bp_strain_dense(HT, 1, 1, bp_bulk)
            call bir_pikus_blocks_free(bp_bulk)
          end if
        end block
      end if


    end subroutine ZB8bandBulk

    !---------------------------------------------------------------------------
    !> Generalized 8-band zincblende Hamiltonian assembly for quantum wires
    !> (ndim=2).  Builds a sparse CSR Hamiltonian of size 8*Ngrid x 8*Ngrid
    !> from precomputed kpterms_2d CSR operators.
    !>
    !> For ndim < 2, delegates to the existing ZB8bandQW dense routine.
    !>
    !> Algorithm:
    !>   1. Build the 10 kp-term CSR matrices (Q, T, S, SC, R, RC, PZ, PP, PM, A)
    !>      as combinations of kpterms_2d operators and kz-dependent scalars.
    !>   2. Iterate over the 8x8 block topology, inserting all nonzeros of
    !>      each nonzero block as COO triplets with appropriate row/col offsets.
    !>   3. Add band-offset profile to diagonal blocks.
    !>   4. Finalize COO -> CSR via csr_build_from_coo.
    !>
    !> The wire has kz as free wavevector (along the wire axis), with
    !> confinement in x-y.  The kpterms_2d operators encode d/dx, d/dy,
    !> d^2/dx^2, d^2/dy^2, and cross derivatives on the 2D grid.
    !---------------------------------------------------------------------------
    subroutine ZB8bandGeneralized(HT_csr, kz, profile_2d, kpterms_2d, cfg, coo_cache, g)

      type(csr_matrix), intent(inout)         :: HT_csr
      real(kind=dp), intent(in)               :: kz
      real(kind=dp), intent(in)               :: profile_2d(:,:)
      type(csr_matrix), intent(in)            :: kpterms_2d(:)
      type(simulation_config), intent(in)     :: cfg
      type(wire_coo_cache), intent(inout), optional :: coo_cache
      character(len=2), intent(in), optional  :: g

      integer :: N, Ntot, alpha, beta, nnz_est
      real(kind=dp) :: kz2
      logical :: gmode
      integer :: gmode_dir  ! 0=none, 1=x, 2=y, 3=z

      ! CSR work matrices for the kp terms
      type(csr_matrix) :: blk_Q, blk_T, blk_S, blk_SC
      type(csr_matrix) :: blk_R, blk_RC
      type(csr_matrix) :: blk_PZ, blk_PP, blk_PM, blk_A
      type(csr_matrix) :: blk_temp, blk_temp2, blk_diff

      ! COO assembly arrays
      integer, allocatable :: coo_rows(:), coo_cols(:)
      complex(kind=dp), allocatable :: coo_vals(:)
      integer :: coo_idx, coo_capacity

      N = grid_ngrid(cfg%grid)
      Ntot = 8 * N
      kz2 = kz * kz

      gmode = .false.
      gmode_dir = 0
      if (present(g)) then
        ! Backward compat: treat 'g ' (len=1 'g' padded) as g3
        if (g == 'g1') gmode_dir = 1
        if (g == 'g2') gmode_dir = 2
        if (g == 'g3' .or. g == 'g ') gmode_dir = 3
        gmode = (gmode_dir > 0)
      end if

      if (gmode_dir == 3) then
        ! ==================================================================
        ! g='g3' mode (z direction): kz-linear perturbation terms.
        ! PP, PM, S, SC at unit kz.  All other blocks are zero.
        ! No band-offset profile is added.
        ! ==================================================================
        call build_kp_term_PP(kz, kpterms_2d, blk_PP)
        call build_kp_term_PM(kz, kpterms_2d, blk_PM)
        call build_kp_term_S(kz, kpterms_2d, blk_S)
        call build_kp_term_SC(kz, kpterms_2d, blk_SC)

      else if (gmode_dir == 1) then
        ! ==================================================================
        ! g='g1' mode (x direction): perturbation via d/dx spatial gradient.
        ! PZ uses kpterms_2d(12) = P * d/dx, S uses kpterms_2d(14) = gamma3 * d/dx.
        ! No PP/PM (those are kz-dependent; kz=0 for x perturbation).
        ! ==================================================================
        ! PZ = -IU * kpterms_2d(12)  (same prefactor as normal mode PZ)
        call csr_scale(kpterms_2d(12), blk_PZ, -IU)
        ! S = 2*sqrt(3) * kpterms_2d(14)
        call csr_scale(kpterms_2d(14), blk_S, &
          cmplx(2.0_dp * SQR3, 0.0_dp, kind=dp))
        ! SC = -S  (antisymmetric gradient)
        call csr_scale(kpterms_2d(14), blk_SC, &
          cmplx(-2.0_dp * SQR3, 0.0_dp, kind=dp))

      else if (gmode_dir == 2) then
        ! ==================================================================
        ! g='g2' mode (y direction): perturbation via d/dy spatial gradient.
        ! PZ uses kpterms_2d(13) = P * d/dy, S uses kpterms_2d(15) = gamma3 * d/dy.
        ! No PP/PM (those are kz-dependent; kz=0 for y perturbation).
        ! ==================================================================
        ! PZ = -IU * kpterms_2d(13)
        call csr_scale(kpterms_2d(13), blk_PZ, -IU)
        ! S = 2*sqrt(3) * kpterms_2d(15)
        call csr_scale(kpterms_2d(15), blk_S, &
          cmplx(2.0_dp * SQR3, 0.0_dp, kind=dp))
        ! SC = -S  (antisymmetric gradient)
        call csr_scale(kpterms_2d(15), blk_SC, &
          cmplx(-2.0_dp * SQR3, 0.0_dp, kind=dp))

      else
        ! ==================================================================
        ! Build kp-term CSR matrices (size N x N each)
        ! ==================================================================

        ! Q = -((gamma1+gamma2)*kz^2*I + kpterms_2d(7))
        call build_kp_term_Q(kz2, kpterms_2d, blk_Q)

        ! T = -(gamma1-gamma2)*kz^2*I - kpterms_2d(8)
        call build_kp_term_T(kz2, kpterms_2d, blk_T)

        ! S = 2*sqrt(3)*kz*kpterms_2d(9)  (wire: kminus = kz since kx,ky are spatial)
        call build_kp_term_S(kz, kpterms_2d, blk_S)

        ! SC = 2*sqrt(3)*kz*kpterms_2d(9)  (wire: kplus = kz since kx,ky are spatial)
        call build_kp_term_SC(kz, kpterms_2d, blk_SC)

        ! R = -sqrt(3)*gamma2*kz^2*I (diagonal only, kx^2,ky^2 absorbed into 2D ops)
        call build_kp_term_R(kz2, kpterms_2d, blk_R)

        ! RC = conjg(R) = R (real diagonal)
        call build_kp_term_RC(kz2, kpterms_2d, blk_RC)

        ! PZ = -IU * kpterms_2d(6)  (P * gradient * (-IU))
        call build_kp_term_PZ(kpterms_2d, blk_PZ)

        ! PP = kpterms_2d(4) * kz * RQS2  (diagonal)
        call build_kp_term_PP(kz, kpterms_2d, blk_PP)

        ! PM = kpterms_2d(4) * kz * RQS2  (same as PP for wire)
        call build_kp_term_PM(kz, kpterms_2d, blk_PM)

        ! A = kpterms_2d(5) + kz^2 * kpterms_2d(10)
        call build_kp_term_A(kz2, kpterms_2d, blk_A)
      end if

      ! ==================================================================
      ! COO capacity estimation and block insertion
      ! ==================================================================
      if (gmode_dir == 3) then
        ! ==================================================================
        ! g='g3' mode (z direction): PP, PM, S, SC blocks only.
        ! Same block topology as the original g='g' mode.
        ! ==================================================================
        ! PP:6, PM:5, S:6, SC:6
        nnz_est = 6*blk_PP%nnz + 5*blk_PM%nnz + 6*blk_S%nnz + 6*blk_SC%nnz
        coo_capacity = nnz_est

        allocate(coo_rows(coo_capacity))
        allocate(coo_cols(coo_capacity))
        allocate(coo_vals(coo_capacity))
        coo_idx = 0

        ! --- Row 1 (HH1) ---
        ! (1,2): SC
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 0, 1, blk_SC, N)
        ! (1,5): -IU * RQS2 * SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 0, 4, blk_SC, N, -IU * RQS2)
        ! (1,7): IU * PP
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 0, 6, blk_PP, N, IU)

        ! --- Row 2 (HH2) ---
        ! (2,1): S
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 1, 0, blk_S, N)
        ! (2,6): -IU * SQR3 * RQS2 * SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 1, 5, blk_SC, N, -IU * SQR3 * RQS2)
        ! (2,8): -RQS3 * PP
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 1, 7, blk_PP, N, cmplx(-RQS3, 0.0_dp, kind=dp))

        ! --- Row 3 (LH1) ---
        ! (3,4): -SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 2, 3, blk_SC, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))
        ! (3,5): IU * SQR3 * RQS2 * S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 2, 4, blk_S, N, IU * SQR3 * RQS2)
        ! (3,7): IU * RQS3 * PM
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 2, 6, blk_PM, N, IU * RQS3)

        ! --- Row 4 (LH2) ---
        ! (4,3): -S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 3, 2, blk_S, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))
        ! (4,6): IU * RQS2 * S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 3, 5, blk_S, N, IU * RQS2)
        ! (4,8): -PM
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 3, 7, blk_PM, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))

        ! --- Row 5 (SO1) ---
        ! (5,1): IU * RQS2 * S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 4, 0, blk_S, N, IU * RQS2)
        ! (5,3): -IU * SQR3 * RQS2 * SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 4, 2, blk_SC, N, -IU * SQR3 * RQS2)
        ! (5,8): IU * SQR2 * RQS3 * PP
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 4, 7, blk_PP, N, IU * SQR2 * RQS3)

        ! --- Row 6 (SO2) ---
        ! (6,2): IU * SQR3 * RQS2 * S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 5, 1, blk_S, N, IU * SQR3 * RQS2)
        ! (6,4): -IU * RQS2 * SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 5, 3, blk_SC, N, -IU * RQS2)
        ! (6,7): SQR2 * RQS3 * PM
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 5, 6, blk_PM, N, cmplx(SQR2 * RQS3, 0.0_dp, kind=dp))

        ! --- Row 7 (CB1) ---
        ! (7,1): -IU * PM
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 6, 0, blk_PM, N, -IU)
        ! (7,3): -IU * RQS3 * PP
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 6, 2, blk_PP, N, -IU * RQS3)
        ! (7,6): SQR2 * RQS3 * PP
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 6, 5, blk_PP, N, cmplx(SQR2 * RQS3, 0.0_dp, kind=dp))

        ! --- Row 8 (CB2) ---
        ! (8,2): -RQS3 * PM
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 7, 1, blk_PM, N, cmplx(-RQS3, 0.0_dp, kind=dp))
        ! (8,4): -PP
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 7, 3, blk_PP, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))
        ! (8,5): -IU * SQR2 * RQS3 * PM
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 7, 4, blk_PM, N, -IU * SQR2 * RQS3)

        ! No profile diagonal in g3 mode

      else if (gmode_dir == 1 .or. gmode_dir == 2) then
        ! ==================================================================
        ! g='g1' or g='g2' mode (x or y direction): PZ, S, SC blocks only.
        ! PZ uses the separate d/dx or d/dy operator (kpterms_2d(12) or (13)).
        ! S/SC use the separate d/dx or d/dy operator (kpterms_2d(14) or (15)).
        ! Block topology matches the normal mode's PZ and S/SC insertions.
        ! PZ:8, S:6, SC:6
        ! ==================================================================
        nnz_est = 8*blk_PZ%nnz + 6*blk_S%nnz + 6*blk_SC%nnz
        coo_capacity = nnz_est

        allocate(coo_rows(coo_capacity))
        allocate(coo_cols(coo_capacity))
        allocate(coo_vals(coo_capacity))
        coo_idx = 0

        ! --- S/SC blocks (same positions as normal mode S/SC) ---

        ! --- Row 1 (HH1) ---
        ! (1,2): SC
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 0, 1, blk_SC, N)
        ! (1,5): -IU * RQS2 * SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 0, 4, blk_SC, N, -IU * RQS2)

        ! --- Row 2 (HH2) ---
        ! (2,1): S
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 1, 0, blk_S, N)
        ! (2,6): -IU * SQR3 * RQS2 * SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 1, 5, blk_SC, N, -IU * SQR3 * RQS2)
        ! (2,7): SQR2 * RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 1, 6, blk_PZ, N, cmplx(SQR2 * RQS3, 0.0_dp, kind=dp))

        ! --- Row 3 (LH1) ---
        ! (3,4): -SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 2, 3, blk_SC, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))
        ! (3,5): IU * SQR3 * RQS2 * S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 2, 4, blk_S, N, IU * SQR3 * RQS2)
        ! (3,8): IU * SQR2 * RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 2, 7, blk_PZ, N, IU * SQR2 * RQS3)

        ! --- Row 4 (LH2) ---
        ! (4,3): -S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 3, 2, blk_S, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))
        ! (4,6): IU * RQS2 * S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 3, 5, blk_S, N, IU * RQS2)

        ! --- Row 5 (SO1) ---
        ! (5,1): IU * RQS2 * S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 4, 0, blk_S, N, IU * RQS2)
        ! (5,3): -IU * SQR3 * RQS2 * SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 4, 2, blk_SC, N, -IU * SQR3 * RQS2)
        ! (5,7): IU * RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 4, 6, blk_PZ, N, IU * RQS3)

        ! --- Row 6 (SO2) ---
        ! (6,2): IU * SQR3 * RQS2 * S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 5, 1, blk_S, N, IU * SQR3 * RQS2)
        ! (6,4): -IU * RQS2 * SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 5, 3, blk_SC, N, -IU * RQS2)
        ! (6,8): -RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 5, 7, blk_PZ, N, cmplx(-RQS3, 0.0_dp, kind=dp))

        ! --- Row 7 (CB1) ---
        ! (7,2): SQR2 * RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 6, 1, blk_PZ, N, cmplx(SQR2 * RQS3, 0.0_dp, kind=dp))
        ! (7,5): -IU * RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 6, 4, blk_PZ, N, -IU * RQS3)

        ! --- Row 8 (CB2) ---
        ! (8,3): -IU * SQR2 * RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 7, 2, blk_PZ, N, -IU * SQR2 * RQS3)
        ! (8,6): -RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 7, 5, blk_PZ, N, cmplx(-RQS3, 0.0_dp, kind=dp))

        ! No profile diagonal in g1/g2 mode

      else
        ! ==================================================================
        ! Exact COO capacity: sum nnz of every block insertion
        ! ==================================================================
        ! blk_diff and blk_temp are built during insert phase; their nnz is
        ! bounded by blk_Q%nnz + blk_T%nnz.  Count block multiplicities:
        ! Q:2, T:2, S:6, SC:6, R:5, RC:3, PZ:8, PP:6, PM:5, A:2,
        ! diff:4, temp:2, profile:8*N
        nnz_est = 2*blk_Q%nnz + 2*blk_T%nnz + 6*blk_S%nnz + 6*blk_SC%nnz
        nnz_est = nnz_est + 5*blk_R%nnz + 3*blk_RC%nnz + 8*blk_PZ%nnz
        nnz_est = nnz_est + 6*blk_PP%nnz + 5*blk_PM%nnz + 2*blk_A%nnz
        ! blk_diff (Q-T): at most Q%nnz + T%nnz entries, inserted 4 times
        nnz_est = nnz_est + 4*(blk_Q%nnz + blk_T%nnz)
        ! blk_temp (0.5*(Q+T)): at most Q%nnz + T%nnz entries, inserted 2 times
        nnz_est = nnz_est + 2*(blk_Q%nnz + blk_T%nnz)
        ! Profile diagonal: 8 bands * N
        nnz_est = nnz_est + 8 * N
        ! Strain: 32 entries per grid point (8 diag + 4 S_eps + 4 R_eps + 16 VB-SO)
        if (allocated(cfg%strain_blocks%delta_Ec)) then
          nnz_est = nnz_est + 32 * N
        end if
        coo_capacity = nnz_est

        allocate(coo_rows(coo_capacity))
        allocate(coo_cols(coo_capacity))
        allocate(coo_vals(coo_capacity))
        coo_idx = 0

        ! ==================================================================
        ! Insert 8x8 blocks into COO arrays
        ! ==================================================================
        ! The block topology follows ZB8bandQW exactly.
        ! Block (alpha, beta) occupies rows [(alpha-1)*N+1 : alpha*N]
        !                       cols [(beta-1)*N+1  : beta*N]

        ! --- Row 1 (HH1) ---
        ! (1,1): Q
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 0, 0, blk_Q, N)
        ! (1,2): SC
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 0, 1, blk_SC, N)
        ! (1,3): RC
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 0, 2, blk_RC, N)
        ! (1,5): -IU * RQS2 * SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 0, 4, blk_SC, N, -IU * RQS2)
        ! (1,6): IU * SQR2 * RC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 0, 5, blk_RC, N, IU * SQR2)
        ! (1,7): IU * PP
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 0, 6, blk_PP, N, IU)

        ! --- Row 2 (HH2) ---
        ! (2,1): S
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 1, 0, blk_S, N)
        ! (2,2): T
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 1, 1, blk_T, N)
        ! (2,4): RC
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 1, 3, blk_RC, N)
        ! (2,5): IU * RQS2 * (Q - T)
        call csr_add(blk_Q, blk_T, blk_diff, UM, cmplx(-1.0_dp, 0.0_dp, kind=dp))
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 1, 4, blk_diff, N, IU * RQS2)
        ! (2,6): -IU * SQR3 * RQS2 * SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 1, 5, blk_SC, N, -IU * SQR3 * RQS2)
        ! (2,7): SQR2 * RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 1, 6, blk_PZ, N, cmplx(SQR2 * RQS3, 0.0_dp, kind=dp))
        ! (2,8): -RQS3 * PP
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 1, 7, blk_PP, N, cmplx(-RQS3, 0.0_dp, kind=dp))

        ! --- Row 3 (LH1) ---
        ! (3,1): R
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 2, 0, blk_R, N)
        ! (3,3): T
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 2, 2, blk_T, N)
        ! (3,4): -SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 2, 3, blk_SC, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))
        ! (3,5): IU * SQR3 * RQS2 * S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 2, 4, blk_S, N, IU * SQR3 * RQS2)
        ! (3,6): IU * RQS2 * (Q - T)
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 2, 5, blk_diff, N, IU * RQS2)
        ! (3,7): IU * RQS3 * PM
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 2, 6, blk_PM, N, IU * RQS3)
        ! (3,8): IU * SQR2 * RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 2, 7, blk_PZ, N, IU * SQR2 * RQS3)

        ! --- Row 4 (LH2) ---
        ! (4,2): R
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 3, 1, blk_R, N)
        ! (4,3): -S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 3, 2, blk_S, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))
        ! (4,4): Q
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 3, 3, blk_Q, N)
        ! (4,5): IU * SQR2 * R
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 3, 4, blk_R, N, IU * SQR2)
        ! (4,6): IU * RQS2 * S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 3, 5, blk_S, N, IU * RQS2)
        ! (4,8): -PM
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 3, 7, blk_PM, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))

        ! --- Row 5 (SO1) ---
        ! (5,1): IU * RQS2 * S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 4, 0, blk_S, N, IU * RQS2)
        ! (5,2): -IU * RQS2 * (Q - T)
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 4, 1, blk_diff, N, -IU * RQS2)
        ! (5,3): -IU * SQR3 * RQS2 * SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 4, 2, blk_SC, N, -IU * SQR3 * RQS2)
        ! (5,4): -IU * SQR2 * RC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 4, 3, blk_RC, N, -IU * SQR2)
        ! (5,5): 0.5*(Q + T)
        call csr_add(blk_Q, blk_T, blk_temp, cmplx(0.5_dp, 0.0_dp, kind=dp), &
          cmplx(0.5_dp, 0.0_dp, kind=dp))
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 4, 4, blk_temp, N)
        call csr_free(blk_temp)
        ! (5,7): IU * RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 4, 6, blk_PZ, N, IU * RQS3)
        ! (5,8): IU * SQR2 * RQS3 * PP
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 4, 7, blk_PP, N, IU * SQR2 * RQS3)

        ! --- Row 6 (SO2) ---
        ! (6,1): -IU * SQR2 * R
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 5, 0, blk_R, N, -IU * SQR2)
        ! (6,2): IU * SQR3 * RQS2 * S
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 5, 1, blk_S, N, IU * SQR3 * RQS2)
        ! (6,3): -IU * RQS2 * (Q - T)
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 5, 2, blk_diff, N, -IU * RQS2)
        ! (6,4): -IU * RQS2 * SC
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 5, 3, blk_SC, N, -IU * RQS2)
        ! (6,6): 0.5*(Q + T)
        call csr_add(blk_Q, blk_T, blk_temp, cmplx(0.5_dp, 0.0_dp, kind=dp), &
          cmplx(0.5_dp, 0.0_dp, kind=dp))
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 5, 5, blk_temp, N)
        call csr_free(blk_temp)
        ! (6,7): SQR2 * RQS3 * PM
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 5, 6, blk_PM, N, cmplx(SQR2 * RQS3, 0.0_dp, kind=dp))
        ! (6,8): -RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 5, 7, blk_PZ, N, cmplx(-RQS3, 0.0_dp, kind=dp))

        ! --- Row 7 (CB1) ---
        ! (7,1): -IU * PM
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 6, 0, blk_PM, N, -IU)
        ! (7,2): SQR2 * RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 6, 1, blk_PZ, N, cmplx(SQR2 * RQS3, 0.0_dp, kind=dp))
        ! (7,3): -IU * RQS3 * PP
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 6, 2, blk_PP, N, -IU * RQS3)
        ! (7,5): -IU * RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 6, 4, blk_PZ, N, -IU * RQS3)
        ! (7,6): SQR2 * RQS3 * PP
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 6, 5, blk_PP, N, cmplx(SQR2 * RQS3, 0.0_dp, kind=dp))
        ! (7,7): A
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 6, 6, blk_A, N)

        ! --- Row 8 (CB2) ---
        ! (8,2): -RQS3 * PM
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 7, 1, blk_PM, N, cmplx(-RQS3, 0.0_dp, kind=dp))
        ! (8,3): -IU * SQR2 * RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 7, 2, blk_PZ, N, -IU * SQR2 * RQS3)
        ! (8,4): -PP
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 7, 3, blk_PP, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))
        ! (8,5): -IU * SQR2 * RQS3 * PM
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 7, 4, blk_PM, N, -IU * SQR2 * RQS3)
        ! (8,6): -RQS3 * PZ
        call insert_csr_block_scaled(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 7, 5, blk_PZ, N, cmplx(-RQS3, 0.0_dp, kind=dp))
        ! (8,8): A
        call insert_csr_block(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, 7, 7, blk_A, N)

        ! ==================================================================
        ! Add band-offset profile to diagonal blocks
        ! ==================================================================
        call insert_profile_diagonal(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, profile_2d, N)

        ! ==================================================================
        ! Add Bir-Pikus strain corrections (diagonal + off-diagonal)
        !
        ! Strain is k-independent so NOT added in g-mode.
        ! All strain terms are per-grid-point scalars (diagonal in spatial
        ! indices), inserted as individual COO entries.
        ! ==================================================================
        if (allocated(cfg%strain_blocks%delta_Ec)) then
          call insert_strain_coo(coo_rows, coo_cols, coo_vals, coo_capacity, &
            coo_idx, cfg%strain_blocks, N)
        end if

      end if

      ! ==================================================================
      ! Build final CSR from COO, or update values if cache is available
      ! ==================================================================
      if (present(coo_cache)) then
        if (coo_cache%initialized) then
          ! Reuse existing CSR structure: O(NNZ) value update only
          call csr_set_values_from_coo(HT_csr, coo_idx, &
            coo_cache%coo_to_csr(1:coo_idx), coo_vals(1:coo_idx))
        else
          ! First call: full build + save COO-to-CSR mapping to cache
          call csr_build_from_coo_cached(HT_csr, Ntot, Ntot, coo_idx, &
            coo_rows(1:coo_idx), coo_cols(1:coo_idx), coo_vals(1:coo_idx), &
            coo_cache%coo_to_csr)
          coo_cache%coo_nnz_in = coo_idx
          coo_cache%initialized = .true.
        end if
      else
        ! No cache: always do full build (backward compatible)
        if (coo_idx > 0) then
          call csr_build_from_coo(HT_csr, Ntot, Ntot, coo_idx, &
            coo_rows(1:coo_idx), coo_cols(1:coo_idx), coo_vals(1:coo_idx))
        else
          call csr_init(HT_csr, Ntot, Ntot)
        end if
      end if

      ! ==================================================================
      ! Cleanup
      ! ==================================================================
      deallocate(coo_rows, coo_cols, coo_vals)
      if (gmode_dir == 0) then
        ! Normal mode: free all blocks
        call csr_free(blk_Q)
        call csr_free(blk_T)
        call csr_free(blk_R)
        call csr_free(blk_RC)
        call csr_free(blk_PZ)
        call csr_free(blk_A)
        call csr_free(blk_diff)
      else if (gmode_dir == 1 .or. gmode_dir == 2) then
        ! g1/g2 mode: PZ, S, SC were built
        call csr_free(blk_PZ)
      end if
      ! S, SC, PP, PM are always built except in g1/g2 mode for PP/PM.
      ! csr_free on an uninitialized csr_matrix is safe (checks allocated).
      call csr_free(blk_S)
      call csr_free(blk_SC)
      call csr_free(blk_PP)
      call csr_free(blk_PM)

    end subroutine ZB8bandGeneralized

    ! ==================================================================
    ! Helper: Build kp-term Q for wire
    ! Q = -((gamma1+gamma2)*kz^2*I + kpterms_2d(7))
    ! kpterms_2d(1) = gamma1 diag, kpterms_2d(2) = gamma2 diag
    ! kpterms_2d(7) = -(gamma1-2gamma2)*Laplacian
    ! ==================================================================
    subroutine build_kp_term_Q(kz2, kpterms_2d, blk)
      real(kind=dp), intent(in) :: kz2
      type(csr_matrix), intent(in) :: kpterms_2d(:)
      type(csr_matrix), intent(out) :: blk

      type(csr_matrix) :: diag_sum, kz2_diag

      ! gamma1 + gamma2 diagonal, scaled by kz^2
      call csr_add(kpterms_2d(1), kpterms_2d(2), diag_sum, &
        cmplx(kz2, 0.0_dp, kind=dp), cmplx(kz2, 0.0_dp, kind=dp))
      ! Add kpterms_2d(7)
      call csr_add(diag_sum, kpterms_2d(7), blk, UM, UM)
      call csr_free(diag_sum)
      ! Negate
      call negate_csr(blk)
    end subroutine build_kp_term_Q

    ! ==================================================================
    ! Helper: Build kp-term T for wire
    ! T = -((gamma1-gamma2)*kz^2*I + kpterms_2d(8))
    ! ==================================================================
    subroutine build_kp_term_T(kz2, kpterms_2d, blk)
      real(kind=dp), intent(in) :: kz2
      type(csr_matrix), intent(in) :: kpterms_2d(:)
      type(csr_matrix), intent(out) :: blk

      type(csr_matrix) :: diag_sum

      ! (gamma1 - gamma2)*kz^2 diagonal
      call csr_add(kpterms_2d(1), kpterms_2d(2), diag_sum, &
        cmplx(kz2, 0.0_dp, kind=dp), cmplx(-kz2, 0.0_dp, kind=dp))
      ! Add kpterms_2d(8) and negate
      call csr_add(diag_sum, kpterms_2d(8), blk, UM, UM)
      call csr_free(diag_sum)
      call negate_csr(blk)
    end subroutine build_kp_term_T

    ! ==================================================================
    ! Helper: Build kp-term S for wire
    ! S = 2*sqrt(3) * kz * kpterms_2d(9)
    ! (In QW: S = 2*sqrt(3)*kminus*kpterms(9) with kminus = kx-i*ky.
    !  For wire, kx,ky are spatial so kminus = kz.)
    ! ==================================================================
    subroutine build_kp_term_S(kz, kpterms_2d, blk)
      real(kind=dp), intent(in) :: kz
      type(csr_matrix), intent(in) :: kpterms_2d(:)
      type(csr_matrix), intent(out) :: blk

      call csr_scale(kpterms_2d(9), blk, &
        cmplx(2.0_dp * SQR3 * kz, 0.0_dp, kind=dp))
    end subroutine build_kp_term_S

    ! ==================================================================
    ! Helper: Build kp-term SC for wire
    ! SC = -S = -2*sqrt(3) * kz * kpterms_2d(9)
    !
    ! In QW: SC(jj,ii) = 2*sqrt(3)*kplus*kpterms(ii,jj,9), where
    ! kplus = kx+i*ky is the conjugate of kminus.  For diagonal kpterms,
    ! SC = conjg(S).  For wire, kpterms_2d(9) = -gamma3*grad_2d is
    ! antisymmetric, so SC = S^H = S^T = -S.
    ! ==================================================================
    subroutine build_kp_term_SC(kz, kpterms_2d, blk)
      real(kind=dp), intent(in) :: kz
      type(csr_matrix), intent(in) :: kpterms_2d(:)
      type(csr_matrix), intent(out) :: blk

      call csr_scale(kpterms_2d(9), blk, &
        cmplx(-2.0_dp * SQR3 * kz, 0.0_dp, kind=dp))
    end subroutine build_kp_term_SC

    ! ==================================================================
    ! Helper: Build kp-term R for wire (diagonal only)
    ! R = -sqrt(3) * gamma2 * kz^2 * I
    ! (kx^2,ky^2 contribution absorbed into 2D operators)
    ! ==================================================================
    subroutine build_kp_term_R(kz2, kpterms_2d, blk)
      real(kind=dp), intent(in) :: kz2
      type(csr_matrix), intent(in) :: kpterms_2d(:)
      type(csr_matrix), intent(out) :: blk

      ! Scale gamma2 diagonal by -sqrt(3)*kz^2
      call csr_scale(kpterms_2d(2), blk, &
        cmplx(-SQR3 * kz2, 0.0_dp, kind=dp))
    end subroutine build_kp_term_R

    ! ==================================================================
    ! Helper: Build kp-term RC for wire
    ! RC = R (real, so conjugate is same)
    ! ==================================================================
    subroutine build_kp_term_RC(kz2, kpterms_2d, blk)
      real(kind=dp), intent(in) :: kz2
      type(csr_matrix), intent(in) :: kpterms_2d(:)
      type(csr_matrix), intent(out) :: blk

      call csr_scale(kpterms_2d(2), blk, &
        cmplx(-SQR3 * kz2, 0.0_dp, kind=dp))
    end subroutine build_kp_term_RC

    ! ==================================================================
    ! Helper: Build kp-term PZ for wire
    ! PZ = -IU * kpterms_2d(6)
    ! kpterms_2d(6) = -P*(d/dx + d/dy) (already negative from coeff application)
    ! So PZ = -IU * (-P*grad) = IU * P*grad
    ! ==================================================================
    subroutine build_kp_term_PZ(kpterms_2d, blk)
      type(csr_matrix), intent(in) :: kpterms_2d(:)
      type(csr_matrix), intent(out) :: blk

      call csr_scale(kpterms_2d(6), blk, -IU)
    end subroutine build_kp_term_PZ

    ! ==================================================================
    ! Helper: Build kp-term PP for wire (diagonal only)
    ! PP = P * kz / sqrt(2)  (kplus = kz for wire)
    ! kpterms_2d(4) = P diagonal
    ! ==================================================================
    subroutine build_kp_term_PP(kz, kpterms_2d, blk)
      real(kind=dp), intent(in) :: kz
      type(csr_matrix), intent(in) :: kpterms_2d(:)
      type(csr_matrix), intent(out) :: blk

      call csr_scale(kpterms_2d(4), blk, &
        cmplx(kz * RQS2, 0.0_dp, kind=dp))
    end subroutine build_kp_term_PP

    ! ==================================================================
    ! Helper: Build kp-term PM for wire (diagonal only)
    ! PM = P * kz / sqrt(2)  (kminus = kz for wire)
    ! ==================================================================
    subroutine build_kp_term_PM(kz, kpterms_2d, blk)
      real(kind=dp), intent(in) :: kz
      type(csr_matrix), intent(in) :: kpterms_2d(:)
      type(csr_matrix), intent(out) :: blk

      call csr_scale(kpterms_2d(4), blk, &
        cmplx(kz * RQS2, 0.0_dp, kind=dp))
    end subroutine build_kp_term_PM

    ! ==================================================================
    ! Helper: Build kp-term A for wire
    ! A = kpterms_2d(5) + kz^2 * kpterms_2d(10)
    ! kpterms_2d(5) = -A*Laplacian
    ! kpterms_2d(10) = A diagonal
    ! ==================================================================
    subroutine build_kp_term_A(kz2, kpterms_2d, blk)
      real(kind=dp), intent(in) :: kz2
      type(csr_matrix), intent(in) :: kpterms_2d(:)
      type(csr_matrix), intent(out) :: blk

      call csr_add(kpterms_2d(5), kpterms_2d(10), blk, &
        UM, cmplx(kz2, 0.0_dp, kind=dp))
    end subroutine build_kp_term_A

    ! ==================================================================
    ! Helper: Insert all nonzeros of a CSR block into COO arrays
    ! with row/col offset by (alpha_off*N, beta_off*N)
    ! ==================================================================
    subroutine insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, &
        alpha_off, beta_off, blk, N)
      integer, intent(inout) :: coo_r(:), coo_c(:)
      complex(kind=dp), intent(inout) :: coo_v(:)
      integer, intent(in) :: coo_cap
      integer, intent(inout) :: coo_idx
      integer, intent(in) :: alpha_off, beta_off, N
      type(csr_matrix), intent(in) :: blk

      integer :: row, k, g_row, g_col

      do row = 1, blk%nrows
        do k = blk%rowptr(row), blk%rowptr(row + 1) - 1
          coo_idx = coo_idx + 1
          if (coo_idx > coo_cap) then
            print *, "WARNING: COO capacity exceeded in insert_csr_block, skipping entries"
            return
          end if
          g_row = alpha_off * N + row
          g_col = beta_off * N + blk%colind(k)
          coo_r(coo_idx) = g_row
          coo_c(coo_idx) = g_col
          coo_v(coo_idx) = blk%values(k)
        end do
      end do
    end subroutine insert_csr_block

    ! ==================================================================
    ! Helper: Insert CSR block with complex scalar multiplier
    ! ==================================================================
    subroutine insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, &
        alpha_off, beta_off, blk, N, scale)
      integer, intent(inout) :: coo_r(:), coo_c(:)
      complex(kind=dp), intent(inout) :: coo_v(:)
      integer, intent(in) :: coo_cap
      integer, intent(inout) :: coo_idx
      integer, intent(in) :: alpha_off, beta_off, N
      type(csr_matrix), intent(in) :: blk
      complex(kind=dp), intent(in) :: scale

      integer :: row, k, g_row, g_col

      do row = 1, blk%nrows
        do k = blk%rowptr(row), blk%rowptr(row + 1) - 1
          coo_idx = coo_idx + 1
          if (coo_idx > coo_cap) then
            print *, "WARNING: COO capacity exceeded in insert_csr_block_scaled, skipping entries"
            return
          end if
          g_row = alpha_off * N + row
          g_col = beta_off * N + blk%colind(k)
          coo_r(coo_idx) = g_row
          coo_c(coo_idx) = g_col
          coo_v(coo_idx) = scale * blk%values(k)
        end do
      end do
    end subroutine insert_csr_block_scaled

    ! ==================================================================
    ! Helper: Insert profile diagonal into COO arrays
    ! Bands 1-4 get profile(:,1)=EV, bands 5-6 get profile(:,2)=EV-DeltaSO,
    ! bands 7-8 get profile(:,3)=EC.
    ! ==================================================================
    subroutine insert_profile_diagonal(coo_r, coo_c, coo_v, coo_cap, &
        coo_idx, profile_2d, N)
      integer, intent(inout) :: coo_r(:), coo_c(:)
      complex(kind=dp), intent(inout) :: coo_v(:)
      integer, intent(in) :: coo_cap
      integer, intent(inout) :: coo_idx
      real(kind=dp), intent(in) :: profile_2d(:,:)
      integer, intent(in) :: N

      integer :: ii, band
      integer :: band_profile_col(8)

      ! Map bands to profile columns
      band_profile_col(1:4) = 1  ! EV
      band_profile_col(5:6) = 2  ! EV - DeltaSO
      band_profile_col(7:8) = 3  ! EC

      do band = 1, 8
        do ii = 1, N
          coo_idx = coo_idx + 1
          if (coo_idx > coo_cap) then
            print *, "WARNING: COO capacity exceeded in insert_profile_diagonal, skipping entries"
            return
          end if
          coo_r(coo_idx) = (band - 1) * N + ii
          coo_c(coo_idx) = (band - 1) * N + ii
          coo_v(coo_idx) = cmplx(profile_2d(ii, band_profile_col(band)), &
            0.0_dp, kind=dp)
        end do
      end do
    end subroutine insert_profile_diagonal

    ! ==================================================================
    ! Helper: Insert Bir-Pikus strain corrections into dense 8Nx8N matrix.
    !
    ! Adds all 32 entries (8 diagonal + 4 S + 4 R + 16 VB-SO) for grid
    ! point ii. Works for both bulk (N=1, ii=1) and QW (any N).
    ! ==================================================================
    subroutine add_bp_strain_dense(HT, ii, N, bp)
      complex(kind=dp), intent(inout) :: HT(:,:)
      integer, intent(in) :: ii, N
      type(bir_pikus_blocks), intent(in) :: bp

      complex(kind=dp) :: R_eps_c, S_eps_c

      R_eps_c = conjg(bp%R_eps(ii))
      S_eps_c = conjg(bp%S_eps(ii))

      ! === Diagonal per-band ===
      HT(      ii,      ii) = HT(      ii,      ii) + bp%delta_EHH(ii)
      HT(  N + ii,  N + ii) = HT(  N + ii,  N + ii) + bp%delta_ELH(ii)
      HT(2*N + ii,2*N + ii) = HT(2*N + ii,2*N + ii) + bp%delta_ELH(ii)
      HT(3*N + ii,3*N + ii) = HT(3*N + ii,3*N + ii) + bp%delta_EHH(ii)
      HT(4*N + ii,4*N + ii) = HT(4*N + ii,4*N + ii) + bp%delta_ESO(ii)
      HT(5*N + ii,5*N + ii) = HT(5*N + ii,5*N + ii) + bp%delta_ESO(ii)
      HT(6*N + ii,6*N + ii) = HT(6*N + ii,6*N + ii) + bp%delta_Ec(ii)
      HT(7*N + ii,7*N + ii) = HT(7*N + ii,7*N + ii) + bp%delta_Ec(ii)

      ! === Off-diagonal: S_eps (HH-LH) ===
      HT(      ii,  N + ii) = HT(      ii,  N + ii) + S_eps_c
      HT(  N + ii,      ii) = HT(  N + ii,      ii) + bp%S_eps(ii)
      HT(2*N + ii,3*N + ii) = HT(2*N + ii,3*N + ii) - S_eps_c
      HT(3*N + ii,2*N + ii) = HT(3*N + ii,2*N + ii) - bp%S_eps(ii)

      ! === Off-diagonal: R_eps (HH-LH) ===
      HT(      ii,2*N + ii) = HT(      ii,2*N + ii) + R_eps_c
      HT(2*N + ii,      ii) = HT(2*N + ii,      ii) + bp%R_eps(ii)
      HT(  N + ii,3*N + ii) = HT(  N + ii,3*N + ii) + R_eps_c
      HT(3*N + ii,  N + ii) = HT(3*N + ii,  N + ii) + bp%R_eps(ii)

      ! === Off-diagonal: VB-SO coupling ===
      HT(      ii,4*N + ii) = HT(      ii,4*N + ii) - IU * RQS2 * S_eps_c
      HT(4*N + ii,      ii) = HT(4*N + ii,      ii) + IU * RQS2 * bp%S_eps(ii)
      HT(      ii,5*N + ii) = HT(      ii,5*N + ii) + IU * SQR2 * R_eps_c
      HT(5*N + ii,      ii) = HT(5*N + ii,      ii) - IU * SQR2 * bp%R_eps(ii)

      HT(  N + ii,4*N + ii) = HT(  N + ii,4*N + ii) + IU * RQS2 * bp%QT2_eps(ii)
      HT(4*N + ii,  N + ii) = HT(4*N + ii,  N + ii) - IU * RQS2 * bp%QT2_eps(ii)
      HT(  N + ii,5*N + ii) = HT(  N + ii,5*N + ii) - IU * SQR3o2 * S_eps_c
      HT(5*N + ii,  N + ii) = HT(5*N + ii,  N + ii) + IU * SQR3o2 * bp%S_eps(ii)

      HT(2*N + ii,4*N + ii) = HT(2*N + ii,4*N + ii) + IU * SQR3o2 * bp%S_eps(ii)
      HT(4*N + ii,2*N + ii) = HT(4*N + ii,2*N + ii) - IU * SQR3o2 * S_eps_c
      HT(2*N + ii,5*N + ii) = HT(2*N + ii,5*N + ii) + IU * RQS2 * bp%QT2_eps(ii)
      HT(5*N + ii,2*N + ii) = HT(5*N + ii,2*N + ii) - IU * RQS2 * bp%QT2_eps(ii)

      HT(3*N + ii,4*N + ii) = HT(3*N + ii,4*N + ii) - IU * SQR2 * bp%R_eps(ii)
      HT(4*N + ii,3*N + ii) = HT(4*N + ii,3*N + ii) + IU * SQR2 * R_eps_c
      HT(3*N + ii,5*N + ii) = HT(3*N + ii,5*N + ii) + IU * RQS2 * bp%S_eps(ii)
      HT(5*N + ii,3*N + ii) = HT(5*N + ii,3*N + ii) - IU * RQS2 * S_eps_c
    end subroutine add_bp_strain_dense

    ! ==================================================================
    ! Helper: Insert Bir-Pikus strain corrections as COO entries.
    !
    ! Adds per-band diagonal shifts and off-diagonal coupling following
    ! the same block topology as ZB8bandBulk.  All terms are per-grid-point
    ! scalars (diagonal in spatial indices).
    !
    ! Block topology (alpha_off, beta_off are 0-based band indices):
    !   Diagonal: bands 1-8, per-band shift
    !   S_eps: (0,1),(1,0),(2,3),(3,2) with sign
    !   R_eps: (0,2),(2,0),(1,3),(3,1)
    !   VB-SO: 16 entries with RQS2, SQR2, sqrt(1.5), QT2 prefactors
    ! ==================================================================
    subroutine insert_strain_coo(coo_r, coo_c, coo_v, coo_cap, &
        coo_idx, bp, N)
      integer, intent(inout) :: coo_r(:), coo_c(:)
      complex(kind=dp), intent(inout) :: coo_v(:)
      integer, intent(in) :: coo_cap
      integer, intent(inout) :: coo_idx
      type(bir_pikus_blocks), intent(in) :: bp
      integer, intent(in) :: N

      integer :: ii

      do ii = 1, N
        ! === Diagonal per-band ===
        ! band 1 (HH): delta_EHH
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = ii; coo_c(coo_idx) = ii
        coo_v(coo_idx) = cmplx(bp%delta_EHH(ii), 0.0_dp, kind=dp)

        ! band 2 (LH): delta_ELH
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = N + ii; coo_c(coo_idx) = N + ii
        coo_v(coo_idx) = cmplx(bp%delta_ELH(ii), 0.0_dp, kind=dp)

        ! band 3 (LH): delta_ELH
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 2*N + ii; coo_c(coo_idx) = 2*N + ii
        coo_v(coo_idx) = cmplx(bp%delta_ELH(ii), 0.0_dp, kind=dp)

        ! band 4 (HH): delta_EHH
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 3*N + ii; coo_c(coo_idx) = 3*N + ii
        coo_v(coo_idx) = cmplx(bp%delta_EHH(ii), 0.0_dp, kind=dp)

        ! band 5 (SO): delta_ESO
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 4*N + ii; coo_c(coo_idx) = 4*N + ii
        coo_v(coo_idx) = cmplx(bp%delta_ESO(ii), 0.0_dp, kind=dp)

        ! band 6 (SO): delta_ESO
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 5*N + ii; coo_c(coo_idx) = 5*N + ii
        coo_v(coo_idx) = cmplx(bp%delta_ESO(ii), 0.0_dp, kind=dp)

        ! band 7 (CB): delta_Ec
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 6*N + ii; coo_c(coo_idx) = 6*N + ii
        coo_v(coo_idx) = cmplx(bp%delta_Ec(ii), 0.0_dp, kind=dp)

        ! band 8 (CB): delta_Ec
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 7*N + ii; coo_c(coo_idx) = 7*N + ii
        coo_v(coo_idx) = cmplx(bp%delta_Ec(ii), 0.0_dp, kind=dp)

        ! === Off-diagonal: S_eps (HH-LH coupling) ===
        ! (1,2): S_eps_c
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = ii; coo_c(coo_idx) = N + ii
        coo_v(coo_idx) = conjg(bp%S_eps(ii))

        ! (2,1): S_eps
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = N + ii; coo_c(coo_idx) = ii
        coo_v(coo_idx) = bp%S_eps(ii)

        ! (3,4): -S_eps_c
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 2*N + ii; coo_c(coo_idx) = 3*N + ii
        coo_v(coo_idx) = -conjg(bp%S_eps(ii))

        ! (4,3): -S_eps
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 3*N + ii; coo_c(coo_idx) = 2*N + ii
        coo_v(coo_idx) = -bp%S_eps(ii)

        ! === Off-diagonal: R_eps (HH-LH coupling) ===
        ! (1,3): R_eps_c
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = ii; coo_c(coo_idx) = 2*N + ii
        coo_v(coo_idx) = conjg(bp%R_eps(ii))

        ! (3,1): R_eps
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 2*N + ii; coo_c(coo_idx) = ii
        coo_v(coo_idx) = bp%R_eps(ii)

        ! (2,4): R_eps_c
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = N + ii; coo_c(coo_idx) = 3*N + ii
        coo_v(coo_idx) = conjg(bp%R_eps(ii))

        ! (4,2): R_eps
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 3*N + ii; coo_c(coo_idx) = N + ii
        coo_v(coo_idx) = bp%R_eps(ii)

        ! === Off-diagonal: VB-SO coupling ===
        ! (1,5): -i/sqrt(2) * S_eps_c
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = ii; coo_c(coo_idx) = 4*N + ii
        coo_v(coo_idx) = -IU * RQS2 * conjg(bp%S_eps(ii))

        ! (5,1): +i/sqrt(2) * S_eps
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 4*N + ii; coo_c(coo_idx) = ii
        coo_v(coo_idx) = IU * RQS2 * bp%S_eps(ii)

        ! (1,6): +i*sqrt(2) * R_eps_c
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = ii; coo_c(coo_idx) = 5*N + ii
        coo_v(coo_idx) = IU * SQR2 * conjg(bp%R_eps(ii))

        ! (6,1): -i*sqrt(2) * R_eps
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 5*N + ii; coo_c(coo_idx) = ii
        coo_v(coo_idx) = -IU * SQR2 * bp%R_eps(ii)

        ! (2,5): +i/sqrt(2) * QT2_eps
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = N + ii; coo_c(coo_idx) = 4*N + ii
        coo_v(coo_idx) = cmplx(IU * RQS2 * bp%QT2_eps(ii), kind=dp)

        ! (5,2): -i/sqrt(2) * QT2_eps
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 4*N + ii; coo_c(coo_idx) = N + ii
        coo_v(coo_idx) = cmplx(-IU * RQS2 * bp%QT2_eps(ii), kind=dp)

        ! (2,6): -i*sqrt(3/2) * S_eps_c
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = N + ii; coo_c(coo_idx) = 5*N + ii
        coo_v(coo_idx) = -IU * SQR3o2 * conjg(bp%S_eps(ii))

        ! (6,2): +i*sqrt(3/2) * S_eps
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 5*N + ii; coo_c(coo_idx) = N + ii
        coo_v(coo_idx) = IU * SQR3o2 * bp%S_eps(ii)

        ! (3,5): +i*sqrt(3/2) * S_eps
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 2*N + ii; coo_c(coo_idx) = 4*N + ii
        coo_v(coo_idx) = IU * SQR3o2 * bp%S_eps(ii)

        ! (5,3): -i*sqrt(3/2) * S_eps_c
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 4*N + ii; coo_c(coo_idx) = 2*N + ii
        coo_v(coo_idx) = -IU * SQR3o2 * conjg(bp%S_eps(ii))

        ! (3,6): +i/sqrt(2) * QT2_eps
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 2*N + ii; coo_c(coo_idx) = 5*N + ii
        coo_v(coo_idx) = cmplx(IU * RQS2 * bp%QT2_eps(ii), kind=dp)

        ! (6,3): -i/sqrt(2) * QT2_eps
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 5*N + ii; coo_c(coo_idx) = 2*N + ii
        coo_v(coo_idx) = cmplx(-IU * RQS2 * bp%QT2_eps(ii), kind=dp)

        ! (4,5): -i*sqrt(2) * R_eps
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 3*N + ii; coo_c(coo_idx) = 4*N + ii
        coo_v(coo_idx) = -IU * SQR2 * bp%R_eps(ii)

        ! (5,4): +i*sqrt(2) * R_eps_c
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 4*N + ii; coo_c(coo_idx) = 3*N + ii
        coo_v(coo_idx) = IU * SQR2 * conjg(bp%R_eps(ii))

        ! (4,6): +i/sqrt(2) * S_eps
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 3*N + ii; coo_c(coo_idx) = 5*N + ii
        coo_v(coo_idx) = IU * RQS2 * bp%S_eps(ii)

        ! (6,4): -i/sqrt(2) * S_eps_c
        coo_idx = coo_idx + 1
        if (coo_idx > coo_cap) return
        coo_r(coo_idx) = 5*N + ii; coo_c(coo_idx) = 3*N + ii
        coo_v(coo_idx) = -IU * RQS2 * conjg(bp%S_eps(ii))
      end do
    end subroutine insert_strain_coo

    ! ==================================================================
    ! Helper: Negate a CSR matrix in-place (multiply all values by -1)
    ! ==================================================================
    subroutine negate_csr(mat)
      type(csr_matrix), intent(inout) :: mat
      integer :: k

      do k = 1, mat%nnz
        mat%values(k) = -mat%values(k)
      end do
    end subroutine negate_csr

    ! ==================================================================
    ! Free the COO cache for symbolic assembly reuse.
    ! ==================================================================
    subroutine wire_coo_cache_free(cache)
      type(wire_coo_cache), intent(inout) :: cache

      if (allocated(cache%coo_to_csr)) deallocate(cache%coo_to_csr)
      cache%coo_nnz_in = 0
      cache%initialized = .false.
    end subroutine wire_coo_cache_free


end module hamiltonianConstructor