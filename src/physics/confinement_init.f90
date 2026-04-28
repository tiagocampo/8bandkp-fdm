module confinement_init

  use definitions
  use finitedifferences
  use sparse_matrices

  implicit none

  interface confinementInitialization
    module procedure confinementInitialization_raw
    module procedure confinementInitialization_cfg
  end interface confinementInitialization

  public :: confinementInitialization
  public :: confinementInitialization_2d

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

  subroutine confinementInitialization_raw(z, startPos, endPos, material, nlayers,&
    & params, confDir, profile, kpterms, FDorder)

    real(kind = dp), intent(in), dimension(:) :: z
    integer, intent(in), dimension(:) :: startPos, endPos
    integer, intent(in) :: nlayers
    character(len = 255), intent(in) :: material(nlayers)
    type(paramStruct), intent(in) :: params(nlayers)
    character(len = 1), intent(in) :: confDir
    real(kind = dp), intent(inout), allocatable, dimension(:,:) :: profile
    real(kind = dp), intent(inout), dimension(:,:,:) :: kpterms
    integer, intent(in), optional :: FDorder

    real(kind=dp), allocatable, dimension(:,:) :: kptermsProfile
    real(kind=dp), allocatable, dimension(:,:) :: forward, central, backward
    real(kind=dp), allocatable, dimension(:) :: diag, offup, offdown
    real(kind=dp), allocatable, dimension(:,:) :: D_inner, D_outer
    real(kind=dp), allocatable, dimension(:) :: g_half

    integer :: i, N, ii, j
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

      ! Last-layer-wins assignment: later layers overwrite earlier ones.
      ! This allows 2-layer configs where barrier covers the full domain
      ! and the well layer overwrites the central region.
      do i = 1, nlayers, 1
        do j = startPos(i), endPos(i)
          profile(j, 1) = params(i)%EV
          profile(j, 2) = params(i)%EV - params(i)%DeltaSO
          profile(j, 3) = params(i)%EC
          kptermsProfile(j, 1) = params(i)%gamma1
          kptermsProfile(j, 2) = params(i)%gamma2
          kptermsProfile(j, 3) = params(i)%gamma3
          kptermsProfile(j, 4) = params(i)%A
          kptermsProfile(j, 5) = params(i)%P
        end do
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

    ! Build tridiagonal averaging matrices (shared by both FD order paths)
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
    forward(N,N) = 1
    backward = transpose(forward)
    central = backward + forward

    if (order == 2) then
      ! ---- Order 2: use existing tridiagonal approach (backward compatible) ----

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
      ! ---- Higher order: conservative variable-coefficient FD ----
      !
      ! Strategy: use the same tridiagonal structure as FDorder=2
      ! (build_kpterm_block), but with higher-order half-point
      ! interpolation of g for the 2nd-derivative terms.
      ! This gives identical results to FDorder=2 for uniform g,
      ! and improved accuracy at interfaces for FDorder>=4.
      !
      ! 2nd-derivative terms: use higher-order g_half interpolation
      ! with the staggered grid (2-point D_inner/D_outer)
      call buildStaggeredD1Inner(N, delta, 2, D_inner)
      call buildStaggeredD1Outer(N, delta, 2, D_outer)
      allocate(g_half(N - 1))

      ! A*kz**2 (term 5)
      call interpolateToHalfPoints(kptermsProfile(1:N,4), N, order, g_half)
      call applyVariableCoeffStaggered(kpterms, g_half, D_inner, D_outer, N, 5)

      ! Q (term 7): (gamma1-2*gamma2)
      call interpolateToHalfPoints(kptermsProfile(1:N,1) - 2.0_dp*kptermsProfile(1:N,2), &
        & N, order, g_half)
      call applyVariableCoeffStaggered(kpterms, g_half, D_inner, D_outer, N, 7)

      ! T (term 8): (gamma1+2*gamma2)
      call interpolateToHalfPoints(kptermsProfile(1:N,1) + 2.0_dp*kptermsProfile(1:N,2), &
        & N, order, g_half)
      call applyVariableCoeffStaggered(kpterms, g_half, D_inner, D_outer, N, 8)

      ! 1st-derivative terms: use same build_kpterm_block approach
      allocate(diag(N))
      allocate(offup(N))
      allocate(offdown(N))

      ! P*kz (term 6)
      call build_kpterm_block(kpterms, kptermsProfile(1:N,5), central, forward, &
        & backward, diag, offup, offdown, N, 6, 1.0_dp/(4.0_dp*delta), .False.)

      ! S -> gamma3*kz (term 9)
      call build_kpterm_block(kpterms, kptermsProfile(1:N,3), central, forward, &
        & backward, diag, offup, offdown, N, 9, 1.0_dp/(4.0_dp*delta), .False.)

    end if

    if (allocated(kptermsProfile)) deallocate(kptermsProfile)
    if (allocated(forward)) deallocate(forward)
    if (allocated(backward)) deallocate(backward)
    if (allocated(central)) deallocate(central)
    if (allocated(diag)) deallocate(diag)
    if (allocated(offup)) deallocate(offup)
    if (allocated(offdown)) deallocate(offdown)
    if (allocated(D_inner)) deallocate(D_inner)
    if (allocated(D_outer)) deallocate(D_outer)
    if (allocated(g_half)) deallocate(g_half)

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
  !> Builds kpterms_2d CSR matrices from 1D FD operators and
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
  !>   9: gamma3 * (d/dx + d/dy)       -- S term, linear-k coupling (LEGACY: not used in H assembly)
  !>  10: A (diagonal only, k^2 coefficient for CB)
  !>  11: gamma3 * d^2/dxdy             -- cross-derivative (for R term)
  !>  12: P * d/dx   (x-gradient only, for g-factor perturbation)
  !>  13: P * d/dy   (y-gradient only, for g-factor perturbation)
  !>  14: gamma3 * d/dx  (x-gradient only, for S term and g-factor perturbation)
  !>  15: gamma3 * d/dy  (y-gradient only, for S term and g-factor perturbation)
  !>  16: gamma2 * (d^2/dx^2 - d^2/dy^2) -- anisotropic Laplacian (for R term)
  !>  17: unused placeholder
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
    logical :: use_cut_cell_faces
    real(kind=dp), parameter :: face_tol = 1.0e-12_dp
    real(kind=dp), parameter :: inactive_barrier_ev = 1.0e3_dp

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
    use_cut_cell_faces = .false.
    if (allocated(grid%face_fraction_x) .and. allocated(grid%face_fraction_y)) then
      use_cut_cell_faces = any(abs(grid%face_fraction_x - 1.0_dp) > face_tol) .or. &
        any(abs(grid%face_fraction_y - 1.0_dp) > face_tol)
      if (allocated(grid%cell_volume)) then
        use_cut_cell_faces = use_cut_cell_faces .or. &
          any(abs(grid%cell_volume - 1.0_dp) > face_tol)
      end if
    end if

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

    ! For rectangular wires with hard-wall Dirichlet boundaries, the
    ! default one-sided closures are problematic in the multiband wire
    ! Hamiltonian. For the default second-order scheme, use cell-centered
    ! Dirichlet closures:
    !   D1: central interior stencil with ghost-cell elimination at walls
    !   D2: symmetric cell-centered Laplacian
    if (order == 2) then
      call apply_dirichlet_order2_first_derivative(D1x, dx)
      call apply_dirichlet_order2_first_derivative(D1y, dy)
      call apply_dirichlet_order2_second_derivative(D2x, dx)
      call apply_dirichlet_order2_second_derivative(D2y, dy)
    end if

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
    ! for building kpterms_2d(12-15) and kron_Iy_D2x, kron_D2y_Ix for
    ! building kpterms_2d(16-17) below)

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
        ! Inactive cells are kept in the rectangular sparse storage.  Move
        ! them out of the physical window so they do not appear as zero modes.
        profile_2d(ij, 1) = -inactive_barrier_ev
        profile_2d(ij, 2) = -inactive_barrier_ev
        profile_2d(ij, 3) =  inactive_barrier_ev
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

      ! Convention: Unlike the 1D QW path where const is applied at Hamiltonian assembly,
      ! the wire path bakes const = hbar^2/(2*m_0) into the gamma/A profiles so that
      ! kpterms_2d operators directly have energy units (eV).

      ! Material parameter profiles (scaled by const = hbar^2/(2m_0))
      ! The gamma and A parameters are dimensionless; multiplying by const
      ! converts the kpterms operators to energy units (eV).
      ! P is NOT scaled because it already includes const: P = sqrt(EP*const).
      prof_gamma1(ij) = params(mid)%gamma1 * const
      prof_gamma2(ij) = params(mid)%gamma2 * const
      prof_gamma3(ij) = params(mid)%gamma3 * const
      prof_P(ij)      = params(mid)%P
      prof_A(ij)      = params(mid)%A * const
      prof_gm12g2(ij) = (params(mid)%gamma1 - 2.0_dp * params(mid)%gamma2) * const
      prof_gp12g2(ij) = (params(mid)%gamma1 + 2.0_dp * params(mid)%gamma2) * const
    end do

    ! ====================================================================
    ! 4. Build kpterms_2d (17 CSR matrices)
    ! ====================================================================
    if (allocated(kpterms_2d)) deallocate(kpterms_2d)
    allocate(kpterms_2d(17))

    ! Terms 1-4: diagonal-only (gamma1, gamma2, gamma3, P)
    ! Build a diagonal CSR from the profile and store directly
    call build_diagonal_csr(ngrid, prof_gamma1, kpterms_2d(1))
    call build_diagonal_csr(ngrid, prof_gamma2, kpterms_2d(2))
    call build_diagonal_csr(ngrid, prof_gamma3, kpterms_2d(3))
    call build_diagonal_csr(ngrid, prof_P, kpterms_2d(4))

    ! Apply cut-cell face fractions only to second-derivative operators.
    ! Rectangular grids allocate all-one face fractions, but those are not a
    ! no-op for csr_apply_variable_coeff because it reconstructs diagonals.
    ! First-derivative operators must keep their explicit Dirichlet signs.
    ! Terms 5,7,8: Laplacian;  Terms 6,9: Gradient.
    ! Term 11 (cross-derivative) mixes x/y connections so face fractions
    ! are not applied (they remain 1.0 by default).
    if (use_cut_cell_faces) then
      ! Term 5: A * Laplacian
      call csr_apply_variable_coeff(laplacian_2d, prof_A, kpterms_2d(5), &
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
    else
      ! No face fractions (rectangular grid or QW mode): standard path
      call csr_apply_variable_coeff(laplacian_2d, prof_A, kpterms_2d(5))
      call csr_apply_variable_coeff(laplacian_2d, prof_gm12g2, kpterms_2d(7))
      call csr_apply_variable_coeff(laplacian_2d, prof_gp12g2, kpterms_2d(8))
    end if
    ! Term 6: P * Gradient
    call csr_apply_variable_coeff(grad_2d, prof_P, kpterms_2d(6))
    ! Term 9: gamma3 * Gradient
    call csr_apply_variable_coeff(grad_2d, prof_gamma3, kpterms_2d(9))

    ! Term 10: A (diagonal only)
    call build_diagonal_csr(ngrid, prof_A, kpterms_2d(10))

    ! Term 11: gamma3 * Cross-derivative
    call csr_apply_variable_coeff(kron_D1y_D1x, prof_gamma3, kpterms_2d(11))

    ! Terms 12-15: separate x/y gradient operators for g-factor perturbation.
    ! These use the same kron products as grad_2d but split by direction.
    ! kron_Iy_D1x = Iy x D1x (x-gradient), kron_D1y_Ix = D1y x Ix (y-gradient).
    ! Term 12: P * d/dx (x-gradient only)
    call csr_apply_variable_coeff(kron_Iy_D1x, prof_P, kpterms_2d(12))
    ! Term 13: P * d/dy (y-gradient only)
    call csr_apply_variable_coeff(kron_D1y_Ix, prof_P, kpterms_2d(13))
    ! Term 14: gamma3 * d/dx (x-gradient only)
    call csr_apply_variable_coeff(kron_Iy_D1x, prof_gamma3, kpterms_2d(14))
    ! Term 15: gamma3 * d/dy (y-gradient only)
    call csr_apply_variable_coeff(kron_D1y_Ix, prof_gamma3, kpterms_2d(15))

    ! Terms 16-17: separate x/y second-derivative operators for R term.
    ! Term 16: gamma2 * (d^2/dx^2 - d^2/dy^2) for R anisotropic part
    ! Term 17: placeholder (unused -- term 16 is sufficient, built from kron_Iy_D2x and kron_D2y_Ix)
    block
      type(csr_matrix) :: tmp_D2x, tmp_D2y
      if (use_cut_cell_faces) then
        call csr_apply_variable_coeff(kron_Iy_D2x, prof_gamma2, tmp_D2x, &
          face_frac_x=grid%face_fraction_x, &
          face_frac_y=grid%face_fraction_y, grid_nx=nx)
        call csr_apply_variable_coeff(kron_D2y_Ix, prof_gamma2, tmp_D2y, &
          face_frac_x=grid%face_fraction_x, &
          face_frac_y=grid%face_fraction_y, grid_nx=nx)
      else
        call csr_apply_variable_coeff(kron_Iy_D2x, prof_gamma2, tmp_D2x)
        call csr_apply_variable_coeff(kron_D2y_Ix, prof_gamma2, tmp_D2y)
      end if
      ! R_diff = gamma2*D2x - gamma2*D2y
      call csr_add(tmp_D2x, tmp_D2y, kpterms_2d(16), UM, cmplx(-1.0_dp, 0.0_dp, kind=dp))
      call csr_free(tmp_D2x)
      call csr_free(tmp_D2y)
    end block

    ! Term 17: placeholder (unused)
    call csr_init(kpterms_2d(17), ngrid, ngrid)

    ! ====================================================================
    ! 5. Cleanup
    ! ====================================================================
    call csr_free(laplacian_2d)
    call csr_free(grad_2d)
    call csr_free(kron_D1y_D1x)
    call csr_free(kron_Iy_D1x)
    call csr_free(kron_D1y_Ix)
    call csr_free(kron_Iy_D2x)
    call csr_free(kron_D2y_Ix)

    deallocate(D2x, D1x, D2y, D1y, Ix, Iy)
    deallocate(cD2x, cD1x, cD2y, cD1y, cIx, cIy)
    deallocate(prof_gamma1, prof_gamma2, prof_gamma3)
    deallocate(prof_P, prof_A, prof_gm12g2, prof_gp12g2)

  end subroutine confinementInitialization_2d

  subroutine apply_dirichlet_order2_first_derivative(D1, dz)
    real(kind=dp), intent(inout) :: D1(:,:)
    real(kind=dp), intent(in) :: dz
    integer :: n, i

    n = size(D1, 1)
    D1 = 0.0_dp
    do i = 2, n - 1
      D1(i, i - 1) = -0.5_dp / dz
      D1(i, i + 1) =  0.5_dp / dz
    end do
    ! Cell-centered Dirichlet walls with ghost elimination:
    ! u_ghost = -u_1 on the left, u_ghost = -u_n on the right.
    if (n >= 2) then
      D1(1, 1) =  0.5_dp / dz
      D1(1, 2) =  0.5_dp / dz
      D1(n, n - 1) = -0.5_dp / dz
      D1(n, n)     = -0.5_dp / dz
    end if
  end subroutine apply_dirichlet_order2_first_derivative

  subroutine apply_dirichlet_order2_second_derivative(D2, dz)
    real(kind=dp), intent(inout) :: D2(:,:)
    real(kind=dp), intent(in) :: dz
    integer :: n, i

    n = size(D2, 1)
    D2 = 0.0_dp
    do i = 1, n
      D2(i, i) = -2.0_dp / dz**2
    end do
    do i = 1, n - 1
      D2(i, i + 1) = 1.0_dp / dz**2
      D2(i + 1, i) = 1.0_dp / dz**2
    end do
    ! Cell-centered grid with Dirichlet walls at the exterior faces:
    ! the ghost-cell elimination increases the boundary diagonal from
    ! -2/h^2 to -3/h^2 at the first and last interior points.
    if (n >= 1) D2(1, 1) = -3.0_dp / dz**2
    if (n >= 2) D2(n, n) = -3.0_dp / dz**2
  end subroutine apply_dirichlet_order2_second_derivative

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
  !> Apply variable coefficient using staggered-grid conservative form.
  !> Computes: kpterms(:,:,term_idx) = -D_outer * diag(g_half) * D_inner
  !>
  !> This correctly discretizes d/dz[g(z)*d/dz] at heterointerfaces
  !> for any FD order, using the staggered-grid (half-point) approach.
  !>
  !> @param[in]  D_inner   (N-1) x N forward half-point 1st-derivative matrix
  !> @param[in]  D_outer   N x (N-1) backward half-point 1st-derivative matrix
  !> @param[in]  g_half    (N-1) coefficient values interpolated to half-points
  !---------------------------------------------------------------------------
  subroutine applyVariableCoeffStaggered(kpterms, g_half, D_inner, D_outer, N, term_idx)

    real(kind=dp), intent(inout), dimension(:,:,:) :: kpterms
    real(kind=dp), intent(in) :: g_half(:)
    real(kind=dp), intent(in) :: D_inner(:,:), D_outer(:,:)
    integer, intent(in) :: N, term_idx

    real(kind=dp), allocatable :: temp(:,:), result(:,:)
    integer :: ii, jj

    ! temp(N-1, N) = diag(g_half) * D_inner
    ! Scale each row of D_inner by g_half
    allocate(temp(N-1, N))
    temp = 0.0_dp
    do jj = 1, N
      do ii = 1, N - 1
        temp(ii, jj) = g_half(ii) * D_inner(ii, jj)
      end do
    end do

    ! result(N, N) = D_outer * temp
    allocate(result(N, N))
    result = 0.0_dp
    do ii = 1, N
      do jj = 1, N
        result(ii, jj) = dot_product(D_outer(ii, :), temp(:, jj))
      end do
    end do

    ! Store with negative sign (kinetic energy operator convention)
    do jj = 1, N
      do ii = 1, N
        kpterms(ii, jj, term_idx) = -result(ii, jj)
      end do
    end do

    deallocate(temp, result)

  end subroutine applyVariableCoeffStaggered

end module confinement_init
