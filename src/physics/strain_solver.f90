module strain_solver

  ! ==============================================================================
  ! Strain solver for semiconductor heterostructures.
  !
  ! Supports two modes:
  !   1. QW biaxial strain (confinement=1): simple algebraic calculation
  !      using lattice mismatch and elastic constants.
  !   2. Wire plane-strain PDE (confinement=2): solve Navier-Cauchy elasticity
  !      on 2D cross-section using MKL PARDISO.
  !
  ! After strain computation, apply Pikus-Bir deformation potentials to shift
  ! band edges in the profile array.
  ! ==============================================================================

  use definitions, only: dp, paramStruct, spatial_grid, strain_config, grid_ngrid

  implicit none

  ! Module-level flag to warn about invalid material_id only once
  logical, save :: warned_invalid_mat = .false.

  private
  public :: strain_result
  public :: strain_result_free
  public :: compute_strain
  public :: apply_pikus_bir

  ! ------------------------------------------------------------------
  ! Strain tensor at each grid point.
  ! For QW: arrays have length grid%ny.
  ! For wire: arrays have length grid%nx * grid%ny.
  ! ------------------------------------------------------------------
  type :: strain_result
    real(kind=dp), allocatable :: eps_xx(:)   ! (Ngrid) strain along wire axis / growth
    real(kind=dp), allocatable :: eps_yy(:)   ! (Ngrid) in-plane strain
    real(kind=dp), allocatable :: eps_zz(:)   ! (Ngrid) in-plane / growth strain
    real(kind=dp), allocatable :: eps_yz(:)   ! (Ngrid) shear strain
  end type strain_result

contains

  ! ==================================================================
  ! Top-level strain computation dispatcher.
  !
  ! Delegates to QW biaxial (ndim=1) or wire plane-strain PDE (ndim=2)
  ! based on grid%ndim.
  !
  ! If all materials are lattice-matched (a0 == a0_ref), strain is zero
  ! everywhere and no PDE solve is performed.
  ! ==================================================================
  subroutine compute_strain(grid, params, material_id, strain_cfg, a0_ref, &
      strain_out)

    type(spatial_grid), intent(in)  :: grid
    type(paramStruct), intent(in)   :: params(:)
    integer, intent(in)             :: material_id(:)
    type(strain_config), intent(in) :: strain_cfg
    real(kind=dp), intent(in)       :: a0_ref
    type(strain_result), intent(out):: strain_out

    integer :: ngrid, ij, mid
    logical :: lattice_matched

    ngrid = grid_ngrid(grid)

    ! Allocate output arrays
    allocate(strain_out%eps_xx(ngrid))
    allocate(strain_out%eps_yy(ngrid))
    allocate(strain_out%eps_zz(ngrid))
    allocate(strain_out%eps_yz(ngrid))
    strain_out%eps_xx = 0.0_dp
    strain_out%eps_yy = 0.0_dp
    strain_out%eps_zz = 0.0_dp
    strain_out%eps_yz = 0.0_dp

    ! Check if strain is enabled
    if (.not. strain_cfg%enabled) return

    ! Check if all materials are lattice-matched to reference
    lattice_matched = .true.
    do ij = 1, ngrid
      mid = material_id(ij)
      if (mid < 1 .or. mid > size(params)) cycle
      if (abs(params(mid)%a0 - a0_ref) > 1.0e-12_dp) then
        lattice_matched = .false.
        exit
      end if
    end do

    if (lattice_matched) return

    ! Dispatch based on grid dimensionality
    select case (grid%ndim)
    case (1)
      call compute_strain_qw(grid, params, material_id, a0_ref, strain_out)
    case (2)
      call compute_strain_wire(grid, params, material_id, a0_ref, strain_out)
    case default
      ! Bulk (ndim=0): no strain computation needed
      return
    end select

  end subroutine compute_strain

  ! ==================================================================
  ! QW biaxial strain (confinement=1).
  !
  ! For each layer with lattice constant a0, relative to reference a0_ref:
  !   eps_0 = (a0 - a0_ref) / a0_ref
  !   eps_xx = eps_yy = eps_0           (in-plane)
  !   eps_zz = -2*C12/C11 * eps_0       (growth direction)
  !   eps_yz = 0
  !
  ! The growth axis is z, the wire-axis convention eps_xx corresponds
  ! to the growth direction via the QW-to-wire axis mapping.
  ! ==================================================================
  subroutine compute_strain_qw(grid, params, material_id, a0_ref, strain_out)

    type(spatial_grid), intent(in)  :: grid
    type(paramStruct), intent(in)   :: params(:)
    integer, intent(in)             :: material_id(:)
    real(kind=dp), intent(in)       :: a0_ref
    type(strain_result), intent(inout) :: strain_out

    integer :: ij, mid, ngrid
    real(kind=dp) :: eps_0, a0_mat

    ngrid = grid_ngrid(grid)

    do ij = 1, ngrid
      mid = material_id(ij)
      if (mid < 1 .or. mid > size(params)) then
        if (.not. warned_invalid_mat) then
          print *, 'WARNING: invalid material_id=', mid, &
            ' at grid point', ij, '(range: 1..', size(params), '). Skipping.'
          warned_invalid_mat = .true.
        end if
        cycle
      end if

      a0_mat = params(mid)%a0
      if (a0_mat <= 0.0_dp) cycle

      eps_0 = (a0_mat - a0_ref) / a0_ref

      ! In-plane (x-y) biaxial strain
      strain_out%eps_xx(ij) = eps_0
      strain_out%eps_yy(ij) = eps_0

      ! Growth direction (z): Poisson contraction/expansion
      if (abs(params(mid)%C11) > 1.0e-12_dp) then
        strain_out%eps_zz(ij) = -2.0_dp * params(mid)%C12 / params(mid)%C11 * eps_0
      end if

      ! No shear in biaxial case
      strain_out%eps_yz(ij) = 0.0_dp
    end do

  end subroutine compute_strain_qw

  ! ==================================================================
  ! Wire plane-strain PDE solver (confinement=2).
  !
  ! Solves the Navier-Cauchy elasticity equations on the 2D cross-section
  ! for displacements u_y, u_z given fixed axial strain eps_xx at each
  ! grid point.
  !
  ! The PDE is:
  !   d(sigma_yy)/dy + d(sigma_yz)/dz = 0
  !   d(sigma_yz)/dy + d(sigma_zz)/dz = 0
  !
  ! With cubic anisotropic stress-strain:
  !   sigma_yy = C11*eps_yy + C12*(eps_zz + eps_xx)
  !   sigma_zz = C12*eps_yy + C11*eps_zz + C12*eps_xx
  !   sigma_yz = 2*C44*eps_yz
  !
  ! Discretized on a 5-point stencil, solved via MKL PARDISO.
  ! Stress-free BCs implemented via cut-cell face fractions.
  ! ==================================================================
  subroutine compute_strain_wire(grid, params, material_id, a0_ref, strain_out)

    type(spatial_grid), intent(in)  :: grid
    type(paramStruct), intent(in)   :: params(:)
    integer, intent(in)             :: material_id(:)
    real(kind=dp), intent(in)       :: a0_ref
    type(strain_result), intent(inout) :: strain_out

    integer :: nx, ny, ngrid, ndof
    integer :: ix, iy, ij, mid, comp, dof
    integer :: nnz_est, coo_idx
    real(kind=dp) :: dx, dy

    ! COO assembly arrays (real-valued for strain PDE)
    integer, allocatable  :: coo_rows(:), coo_cols(:)
    real(kind=dp), allocatable :: coo_vals(:)
    real(kind=dp), allocatable :: rhs(:), sol(:)

    ! PARDISO arrays
    integer(8) :: pt(64)
    integer :: iparm(64)
    integer :: maxfct, mnum, mtype, phase, nrhs, msglvl, error
    integer, allocatable :: perm(:)

    ! CSR arrays for PARDISO
    integer, allocatable :: ia(:), ja(:)
    real(kind=dp), allocatable :: a_csr(:)

    ! Elastic constants at each grid point
    real(kind=dp), allocatable :: C11(:), C12(:), C44(:)
    real(kind=dp), allocatable :: eps_xx_fixed(:)

    ! Face fraction averages for box-integration weighting
    real(kind=dp), allocatable :: ff_x(:), ff_y(:)

    ! Working variables for stencil
    real(kind=dp) :: c11_c, c11_e, c11_w, c11_n, c11_s
    real(kind=dp) :: c12_c, c44_c, c44_n, c44_s
    real(kind=dp) :: c12p44_c, c12p44_n, c12p44_s
    real(kind=dp) :: rhs_y, rhs_z
    real(kind=dp) :: ffx_c, ffy_c
    real(kind=dp) :: c11_avg_e, c11_avg_w, c11_avg_n, c11_avg_s
    real(kind=dp) :: c44_avg_e, c44_avg_w, c44_avg_n, c44_avg_s
    real(kind=dp) :: c12p44_avg_en, c12p44_avg_wn, c12p44_avg_es, c12p44_avg_ws

    nx    = grid%nx
    ny    = grid%ny
    ngrid = nx * ny
    ndof  = 2 * ngrid   ! two DOFs per grid point: u_y, u_z
    dx    = grid%dx
    dy    = grid%dy

    ! Allocate elastic constant arrays
    allocate(C11(ngrid), C12(ngrid), C44(ngrid))
    allocate(eps_xx_fixed(ngrid))
    C11 = 0.0_dp; C12 = 0.0_dp; C44 = 0.0_dp
    eps_xx_fixed = 0.0_dp

    ! Fill elastic constants and axial strain at each grid point
    do ij = 1, ngrid
      mid = material_id(ij)
      if (mid < 1 .or. mid > size(params)) then
        if (.not. warned_invalid_mat) then
          print *, 'WARNING: invalid material_id=', mid, &
            ' at grid point', ij, '(range: 1..', size(params), '). Skipping.'
          warned_invalid_mat = .true.
        end if
        cycle
      end if

      C11(ij) = params(mid)%C11
      C12(ij) = params(mid)%C12
      C44(ij) = params(mid)%C44

      if (params(mid)%a0 > 0.0_dp) then
        eps_xx_fixed(ij) = (params(mid)%a0 - a0_ref) / a0_ref
      end if
    end do

    ! Set axial strain in output
    strain_out%eps_xx = eps_xx_fixed

    ! Allocate face fraction averages
    allocate(ff_x(ngrid), ff_y(ngrid))
    ff_x = 1.0_dp
    ff_y = 1.0_dp

    if (allocated(grid%face_fraction_x) .and. allocated(grid%face_fraction_y)) then
      do ij = 1, ngrid
        ff_x(ij) = 0.5_dp * (grid%face_fraction_x(ij, 1) + grid%face_fraction_x(ij, 2))
        ff_y(ij) = 0.5_dp * (grid%face_fraction_y(ij, 1) + grid%face_fraction_y(ij, 2))
      end do
    end if

    ! ==================================================================
    ! Estimate COO size and allocate
    ! ==================================================================
    ! Each active grid point contributes up to 13 entries (diagonal + 4
    ! neighbors x 2 components + cross-derivative coupling).
    ! For 2 DOFs per point, each interior point gives ~13 nonzeros per
    ! row.  Total upper bound: ndof * 13.
    nnz_est = ndof * 13
    allocate(coo_rows(nnz_est), coo_cols(nnz_est), coo_vals(nnz_est))
    coo_idx = 0

    allocate(rhs(ndof), sol(ndof))
    rhs = 0.0_dp
    sol = 0.0_dp

    ! ==================================================================
    ! Assemble stiffness matrix in COO format
    ! ==================================================================
    ! DOF ordering: comp=1 is u_y, comp=2 is u_z
    ! DOF index = (comp-1)*ngrid + ij

    do iy = 1, ny
      do ix = 1, nx
        ij = wire_flat_idx(nx, ix, iy)
        mid = material_id(ij)

        ! Skip inactive cells (cell_volume == 0)
        if (allocated(grid%cell_volume)) then
          if (grid%cell_volume(ij) < 0.5_dp) cycle
        end if

        c11_c = C11(ij)
        c12_c = C12(ij)
        c44_c = C44(ij)
        c12p44_c = c12_c + c44_c
        ffx_c = ff_x(ij)
        ffy_c = ff_y(ij)

        ! Neighbor elastic constants (use current cell if at boundary)
        if (ix > 1) then
          c11_w = C11(wire_flat_idx(nx, ix-1, iy))
        else
          c11_w = c11_c
        end if
        if (ix < nx) then
          c11_e = C11(wire_flat_idx(nx, ix+1, iy))
        else
          c11_e = c11_c
        end if
        if (iy > 1) then
          c44_s = C44(wire_flat_idx(nx, ix, iy-1))
          c12p44_s = C12(wire_flat_idx(nx, ix, iy-1)) + C44(wire_flat_idx(nx, ix, iy-1))
        else
          c44_s = c44_c
          c12p44_s = c12p44_c
        end if
        if (iy < ny) then
          c11_n = C11(wire_flat_idx(nx, ix, iy+1))
          c44_n = C44(wire_flat_idx(nx, ix, iy+1))
          c12p44_n = C12(wire_flat_idx(nx, ix, iy+1)) + C44(wire_flat_idx(nx, ix, iy+1))
        else
          c11_n = c11_c
          c44_n = c44_c
          c12p44_n = c12p44_c
        end if

        ! Face-fraction-weighted averages for interface elastic constants
        ! East face average: (C_east + C_center)/2 * ff
        if (ix < nx) then
          c11_avg_e = 0.5_dp * (c11_c + c11_e) * ffx_c
          c44_avg_e = 0.5_dp * (c44_c + C44(wire_flat_idx(nx, ix+1, iy))) * ffx_c
        else
          c11_avg_e = 0.0_dp  ! boundary: no flux
          c44_avg_e = 0.0_dp
        end if
        ! West face average
        if (ix > 1) then
          c11_avg_w = 0.5_dp * (c11_c + c11_w) * ffx_c
          c44_avg_w = 0.5_dp * (c44_c + C44(wire_flat_idx(nx, ix-1, iy))) * ffx_c
        else
          c11_avg_w = 0.0_dp
          c44_avg_w = 0.0_dp
        end if
        ! North face average
        if (iy < ny) then
          c44_avg_n = 0.5_dp * (c44_c + c44_n) * ffy_c
          c11_avg_n = 0.5_dp * (c11_c + c11_n) * ffy_c
        else
          c44_avg_n = 0.0_dp
          c11_avg_n = 0.0_dp
        end if
        ! South face average
        if (iy > 1) then
          c44_avg_s = 0.5_dp * (c44_c + c44_s) * ffy_c
          c11_avg_s = 0.5_dp * (c11_c + C11(wire_flat_idx(nx, ix, iy-1))) * ffy_c
        else
          c44_avg_s = 0.0_dp
          c11_avg_s = 0.0_dp
        end if

        ! Cross-derivative coupling (C12+C44) at corners
        if (ix < nx .and. iy < ny) then
          c12p44_avg_en = 0.5_dp * (c12p44_c + c12p44_n) * &
                          sqrt(ff_x(wire_flat_idx(nx, ix, iy+1)) * ffx_c)
        else
          c12p44_avg_en = 0.0_dp
        end if
        if (ix > 1 .and. iy < ny) then
          c12p44_avg_wn = 0.5_dp * (c12p44_c + c12p44_n) * &
                          sqrt(ff_x(wire_flat_idx(nx, ix-1, iy+1)) * ffx_c)
        else
          c12p44_avg_wn = 0.0_dp
        end if
        if (ix < nx .and. iy > 1) then
          c12p44_avg_es = 0.5_dp * (c12p44_c + c12p44_s) * &
                          sqrt(ff_x(wire_flat_idx(nx, ix, iy-1)) * ffx_c)
        else
          c12p44_avg_es = 0.0_dp
        end if
        if (ix > 1 .and. iy > 1) then
          c12p44_avg_ws = 0.5_dp * (c12p44_c + c12p44_s) * &
                          sqrt(ff_x(wire_flat_idx(nx, ix-1, iy-1)) * ffx_c)
        else
          c12p44_avg_ws = 0.0_dp
        end if

        ! ============================================================
        ! Row for u_y equation (comp=1):
        !   C11*d2(uy)/dy2 + C44*d2(uy)/dz2 + (C12+C44)*d2(uz)/dydz
        !   = -dC12/dy * eps_xx  (source from axial mismatch gradient)
        !
        ! Grid-to-PDE axis mapping (wire cross-section is the x-y grid):
        !   grid x-direction (ix) -> y  in PDE
        !   grid z-direction (iy) -> z  in PDE
        ! ============================================================

        ! --- u_y equation ---
        dof = ij  ! comp=1: DOF = ij

        ! Diagonal: -(c11_avg_e + c11_avg_w)/dx^2 - (c44_avg_n + c44_avg_s)/dy^2
        call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
          dof, dof, &
          -(c11_avg_e + c11_avg_w) / dx**2 - (c44_avg_n + c44_avg_s) / dy**2)

        ! East neighbor (ix+1): u_y coupling
        if (ix < nx) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, wire_flat_idx(nx, ix+1, iy), c11_avg_e / dx**2)
        end if

        ! West neighbor (ix-1): u_y coupling
        if (ix > 1) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, wire_flat_idx(nx, ix-1, iy), c11_avg_w / dx**2)
        end if

        ! North neighbor (iy+1): u_y coupling (C44 in z-direction)
        if (iy < ny) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, wire_flat_idx(nx, ix, iy+1), c44_avg_n / dy**2)
        end if

        ! South neighbor (iy-1): u_y coupling
        if (iy > 1) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, wire_flat_idx(nx, ix, iy-1), c44_avg_s / dy**2)
        end if

        ! Cross-derivative coupling: (C12+C44)/(4*dx*dy) * (uz_{i+1,j+1} - uz_{i+1,j-1}
        !                                                 - uz_{i-1,j+1} + uz_{i-1,j-1})
        ! Northeast (ix+1, iy+1): uz coupling (+c12p44/(4*dx*dy))
        if (ix < nx .and. iy < ny) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, ngrid + wire_flat_idx(nx, ix+1, iy+1), &
            c12p44_avg_en / (4.0_dp * dx * dy))
        end if
        ! Northwest (ix-1, iy+1): uz coupling (-c12p44/(4*dx*dy))
        if (ix > 1 .and. iy < ny) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, ngrid + wire_flat_idx(nx, ix-1, iy+1), &
            -c12p44_avg_wn / (4.0_dp * dx * dy))
        end if
        ! Southeast (ix+1, iy-1): uz coupling (-c12p44/(4*dx*dy))
        if (ix < nx .and. iy > 1) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, ngrid + wire_flat_idx(nx, ix+1, iy-1), &
            -c12p44_avg_es / (4.0_dp * dx * dy))
        end if
        ! Southwest (ix-1, iy-1): uz coupling (+c12p44/(4*dx*dy))
        if (ix > 1 .and. iy > 1) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, ngrid + wire_flat_idx(nx, ix-1, iy-1), &
            c12p44_avg_ws / (4.0_dp * dx * dy))
        end if

        ! RHS for u_y equation: source from gradient of C12*eps_xx
        ! d(C12*eps_xx)/dy approximated by central differences
        rhs_y = 0.0_dp
        if (ix > 1 .and. ix < nx) then
          rhs_y = -(0.5_dp * (C12(wire_flat_idx(nx, ix+1, iy)) * eps_xx_fixed(wire_flat_idx(nx, ix+1, iy)) &
                    - C12(wire_flat_idx(nx, ix-1, iy)) * eps_xx_fixed(wire_flat_idx(nx, ix-1, iy))) / dx)
        else if (ix == 1 .and. nx > 1) then
          rhs_y = -(C12(wire_flat_idx(nx, ix+1, iy)) * eps_xx_fixed(wire_flat_idx(nx, ix+1, iy)) &
                    - C12(ij) * eps_xx_fixed(ij)) / dx
        else if (ix == nx .and. nx > 1) then
          rhs_y = -(C12(ij) * eps_xx_fixed(ij) &
                    - C12(wire_flat_idx(nx, ix-1, iy)) * eps_xx_fixed(wire_flat_idx(nx, ix-1, iy))) / dx
        end if
        rhs(dof) = rhs_y

        ! ============================================================
        ! Row for u_z equation (comp=2):
        !   C44*d2(uz)/dy2 + C11*d2(uz)/dz2 + (C12+C44)*d2(uy)/dydz
        !   = -dC12/dz * eps_xx
        ! ============================================================
        dof = ngrid + ij

        ! Diagonal
        call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
          dof, dof, &
          -(c44_avg_e + c44_avg_w) / dx**2 - (c11_avg_n + c11_avg_s) / dy**2)

        ! East neighbor (ix+1): u_z coupling (C44 in x-direction for u_z eq)
        if (ix < nx) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, ngrid + wire_flat_idx(nx, ix+1, iy), c44_avg_e / dx**2)
        end if

        ! West neighbor (ix-1): u_z coupling
        if (ix > 1) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, ngrid + wire_flat_idx(nx, ix-1, iy), c44_avg_w / dx**2)
        end if

        ! North neighbor (iy+1): u_z coupling (C11 in z-direction)
        if (iy < ny) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, ngrid + wire_flat_idx(nx, ix, iy+1), c11_avg_n / dy**2)
        end if

        ! South neighbor (iy-1): u_z coupling
        if (iy > 1) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, ngrid + wire_flat_idx(nx, ix, iy-1), c11_avg_s / dy**2)
        end if

        ! Cross-derivative coupling: same corners but coupling to u_y
        ! Northeast (ix+1, iy+1)
        if (ix < nx .and. iy < ny) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, wire_flat_idx(nx, ix+1, iy+1), &
            c12p44_avg_en / (4.0_dp * dx * dy))
        end if
        ! Northwest (ix-1, iy+1)
        if (ix > 1 .and. iy < ny) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, wire_flat_idx(nx, ix-1, iy+1), &
            -c12p44_avg_wn / (4.0_dp * dx * dy))
        end if
        ! Southeast (ix+1, iy-1)
        if (ix < nx .and. iy > 1) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, wire_flat_idx(nx, ix+1, iy-1), &
            -c12p44_avg_es / (4.0_dp * dx * dy))
        end if
        ! Southwest (ix-1, iy-1)
        if (ix > 1 .and. iy > 1) then
          call add_coo_entry(coo_rows, coo_cols, coo_vals, nnz_est, coo_idx, &
            dof, wire_flat_idx(nx, ix-1, iy-1), &
            c12p44_avg_ws / (4.0_dp * dx * dy))
        end if

        ! RHS for u_z equation: -d(C12*eps_xx)/dz
        rhs_z = 0.0_dp
        if (iy > 1 .and. iy < ny) then
          rhs_z = -(0.5_dp * (C12(wire_flat_idx(nx, ix, iy+1)) * eps_xx_fixed(wire_flat_idx(nx, ix, iy+1)) &
                    - C12(wire_flat_idx(nx, ix, iy-1)) * eps_xx_fixed(wire_flat_idx(nx, ix, iy-1))) / dy)
        else if (iy == 1 .and. ny > 1) then
          rhs_z = -(C12(wire_flat_idx(nx, ix, iy+1)) * eps_xx_fixed(wire_flat_idx(nx, ix, iy+1)) &
                    - C12(ij) * eps_xx_fixed(ij)) / dy
        else if (iy == ny .and. ny > 1) then
          rhs_z = -(C12(ij) * eps_xx_fixed(ij) &
                    - C12(wire_flat_idx(nx, ix, iy-1)) * eps_xx_fixed(wire_flat_idx(nx, ix, iy-1))) / dy
        end if
        rhs(dof) = rhs_z

      end do
    end do

    ! ==================================================================
    ! Convert COO to CSR for PARDISO
    ! ==================================================================
    call coo_to_csr_real(ndof, ndof, coo_idx, coo_rows, coo_cols, coo_vals, &
      ia, ja, a_csr)

    deallocate(coo_rows, coo_cols, coo_vals)
    deallocate(C11, C12, C44, eps_xx_fixed)
    deallocate(ff_x, ff_y)

    ! ==================================================================
    ! Solve with MKL PARDISO
    ! ==================================================================
    ! Initialize PARDISO
    pt = 0
    iparm = 0
    maxfct = 1
    mnum = 1
    nrhs = 1
    error = 0

    ! Matrix type: real symmetric indefinite
    ! The strain stiffness matrix is positive semi-definite (null space
    ! from rigid-body translation with all-Neumann BCs).
    ! Use mtype=-2 (real symmetric indefinite) with Bunch-Kaufman pivoting.
    mtype = -2
    msglvl = 0  ! no output

    ! Set default iparm values
    iparm(1)  = 1     ! no solver default
    iparm(2)  = 2     ! fill-in reordering from METIS
    iparm(8)  = 10    ! max number of iterative refinement steps
    iparm(10) = 13    ! perturb the pivot elements with 1E-13
    iparm(11) = 1     ! use nonsymmetric permutation and scaling MPS
    iparm(13) = 1     ! maximum weighted matching algorithm
    iparm(21) = 1     ! 1x1 and 2x2 Bunch-Kaufman pivoting
    iparm(24) = 2     ! parallel factorization control
    iparm(25) = 0     ! no parallel reordering

    allocate(perm(ndof))

    ! Phase 11: Analysis + Factorization + Solve
    phase = 13
    call pardiso(pt, maxfct, mnum, mtype, phase, ndof, a_csr, ia, ja, &
      perm, nrhs, iparm, msglvl, rhs, sol, error)

    if (error /= 0) then
      print *, 'ERROR: Strain PDE solve failed (PARDISO error=', error, '). Aborting.'
      ! Release PARDISO internal memory before cleanup
      phase = -1
      call pardiso(pt, maxfct, mnum, mtype, phase, ndof, a_csr, ia, ja, &
        perm, nrhs, iparm, msglvl, rhs, sol, error)
      deallocate(ia, ja, a_csr, rhs, sol, perm)
      stop 1
    end if

    ! Phase -1: Release memory
    phase = -1
    call pardiso(pt, maxfct, mnum, mtype, phase, ndof, a_csr, ia, ja, &
      perm, nrhs, iparm, msglvl, rhs, sol, error)

    deallocate(ia, ja, a_csr, rhs, perm)

    ! ==================================================================
    ! Compute strain from displacement solution
    ! ==================================================================
    ! eps_yy = du_y/dy,  eps_zz = du_z/dz,
    ! eps_yz = 0.5*(du_y/dz + du_z/dy)
    ! Using first-order forward/backward FD (central where possible).
    do iy = 1, ny
      do ix = 1, nx
        ij = wire_flat_idx(nx, ix, iy)

        ! Skip inactive cells
        if (allocated(grid%cell_volume)) then
          if (grid%cell_volume(ij) < 0.5_dp) cycle
        end if

        ! eps_yy = du_y/dx (y-direction is grid x)
        if (ix > 1 .and. ix < nx) then
          strain_out%eps_yy(ij) = &
            (sol(wire_flat_idx(nx, ix+1, iy)) - sol(wire_flat_idx(nx, ix-1, iy))) / (2.0_dp * dx)
        else if (ix == 1 .and. nx > 1) then
          strain_out%eps_yy(ij) = &
            (sol(wire_flat_idx(nx, ix+1, iy)) - sol(ij)) / dx
        else if (ix == nx .and. nx > 1) then
          strain_out%eps_yy(ij) = &
            (sol(ij) - sol(wire_flat_idx(nx, ix-1, iy))) / dx
        end if

        ! eps_zz = du_z/dz (z-direction is grid z/y)
        if (iy > 1 .and. iy < ny) then
          strain_out%eps_zz(ij) = &
            (sol(ngrid + wire_flat_idx(nx, ix, iy+1)) - sol(ngrid + wire_flat_idx(nx, ix, iy-1))) / (2.0_dp * dy)
        else if (iy == 1 .and. ny > 1) then
          strain_out%eps_zz(ij) = &
            (sol(ngrid + wire_flat_idx(nx, ix, iy+1)) - sol(ngrid + ij)) / dy
        else if (iy == ny .and. ny > 1) then
          strain_out%eps_zz(ij) = &
            (sol(ngrid + ij) - sol(ngrid + wire_flat_idx(nx, ix, iy-1))) / dy
        end if

        ! eps_yz = 0.5*(du_y/dz + du_z/dy)
        strain_out%eps_yz(ij) = 0.0_dp

        ! du_y/dz contribution
        if (iy > 1 .and. iy < ny) then
          strain_out%eps_yz(ij) = strain_out%eps_yz(ij) + &
            (sol(wire_flat_idx(nx, ix, iy+1)) - sol(wire_flat_idx(nx, ix, iy-1))) / (2.0_dp * dy)
        else if (iy == 1 .and. ny > 1) then
          strain_out%eps_yz(ij) = strain_out%eps_yz(ij) + &
            (sol(wire_flat_idx(nx, ix, iy+1)) - sol(ij)) / dy
        else if (iy == ny .and. ny > 1) then
          strain_out%eps_yz(ij) = strain_out%eps_yz(ij) + &
            (sol(ij) - sol(wire_flat_idx(nx, ix, iy-1))) / dy
        end if

        ! du_z/dy contribution
        if (ix > 1 .and. ix < nx) then
          strain_out%eps_yz(ij) = strain_out%eps_yz(ij) + &
            (sol(ngrid + wire_flat_idx(nx, ix+1, iy)) - sol(ngrid + wire_flat_idx(nx, ix-1, iy))) / (2.0_dp * dx)
        else if (ix == 1 .and. nx > 1) then
          strain_out%eps_yz(ij) = strain_out%eps_yz(ij) + &
            (sol(ngrid + wire_flat_idx(nx, ix+1, iy)) - sol(ngrid + ij)) / dx
        else if (ix == nx .and. nx > 1) then
          strain_out%eps_yz(ij) = strain_out%eps_yz(ij) + &
            (sol(ngrid + ij) - sol(ngrid + wire_flat_idx(nx, ix-1, iy))) / dx
        end if

        strain_out%eps_yz(ij) = 0.5_dp * strain_out%eps_yz(ij)

      end do
    end do

    deallocate(sol)

  end subroutine compute_strain_wire

  ! ==================================================================
  ! Apply Pikus-Bir strain Hamiltonian to profile_2d band offsets.
  !
  ! Modifies profile_2d in-place to include diagonal strain shifts:
  !   delta_Ec  = ac * Tr(eps)
  !   delta_EHH = -P_eps + Q_eps
  !   delta_ELH = -P_eps - Q_eps
  !   delta_ESO = -P_eps
  !
  ! where:
  !   Tr_eps = eps_xx + eps_yy + eps_zz
  !   P_eps  = -av * Tr_eps    (positive av convention)
  !   Q_eps  = b/2 * (eps_zz - 0.5*(eps_yy + eps_xx))
  !
  ! profile_2d convention (from hamiltonianConstructor):
  !   profile_2d(:,1) = EV          (bands 1-4: HH, LH, LH, HH)
  !   profile_2d(:,2) = EV-DeltaSO  (bands 5-6: SO)
  !   profile_2d(:,3) = EC          (bands 7-8: CB)
  !
  ! The Pikus-Bir shifts modify EV and EC:
  !   EV_new = EV + delta_EHH (for HH) or EV + delta_ELH (for LH)
  !   EC_new = EC + delta_Ec
  !
  ! Since profile_2d(:,1) is used for all 4 VB bands (HH, LH, LH, HH)
  ! we apply the HH shift. The LH correction is a sub-band splitting
  ! handled in the off-diagonal Hamiltonian terms, not in profile.
  ! For SO bands, delta_ESO is applied to profile_2d(:,2).
  ! ==================================================================
  subroutine apply_pikus_bir(strain_out, params, material_id, grid, profile_2d)

    type(strain_result), intent(in)    :: strain_out
    type(paramStruct), intent(in)      :: params(:)
    integer, intent(in)                :: material_id(:)
    type(spatial_grid), intent(in)     :: grid
    real(kind=dp), intent(inout)       :: profile_2d(:,:)

    integer :: ngrid, ij, mid
    real(kind=dp) :: Tr_eps, P_eps, Q_eps
    real(kind=dp) :: delta_Ec, delta_EHH, delta_ESO
    real(kind=dp) :: av, ac, b_dp

    ngrid = grid_ngrid(grid)

    ! Early exit if no strain data
    if (.not. allocated(strain_out%eps_xx)) return
    if (size(strain_out%eps_xx) < ngrid) return

    do ij = 1, ngrid
      mid = material_id(ij)
      if (mid < 1 .or. mid > size(params)) then
        if (.not. warned_invalid_mat) then
          print *, 'WARNING: invalid material_id=', mid, &
            ' at grid point', ij, '(range: 1..', size(params), '). Skipping.'
          warned_invalid_mat = .true.
        end if
        cycle
      end if

      Tr_eps = strain_out%eps_xx(ij) + strain_out%eps_yy(ij) + strain_out%eps_zz(ij)

      ! Check if strain is negligible
      if (abs(Tr_eps) < 1.0e-20_dp .and. &
          abs(strain_out%eps_yz(ij)) < 1.0e-20_dp) cycle

      av   = params(mid)%av
      ac   = params(mid)%ac
      b_dp = params(mid)%b_dp

      ! Pikus-Bir deformation potentials
      P_eps = -av * Tr_eps
      Q_eps = b_dp * 0.5_dp * (strain_out%eps_zz(ij) - &
        0.5_dp * (strain_out%eps_yy(ij) + strain_out%eps_xx(ij)))

      delta_Ec  = ac * Tr_eps
      delta_EHH = -P_eps + Q_eps
      delta_ESO = -P_eps

      ! Apply to profile_2d:
      ! Column 1 (EV, bands 1-4): shift by delta_EHH
      ! Column 2 (EV-DeltaSO, bands 5-6): shift by delta_ESO
      ! Column 3 (EC, bands 7-8): shift by delta_Ec
      profile_2d(ij, 1) = profile_2d(ij, 1) + delta_EHH
      profile_2d(ij, 2) = profile_2d(ij, 2) + delta_ESO
      profile_2d(ij, 3) = profile_2d(ij, 3) + delta_Ec

    end do

  end subroutine apply_pikus_bir

  ! ==================================================================
  ! Deallocate strain_result arrays
  ! ==================================================================
  subroutine strain_result_free(s)
    type(strain_result), intent(inout) :: s

    if (allocated(s%eps_xx)) deallocate(s%eps_xx)
    if (allocated(s%eps_yy)) deallocate(s%eps_yy)
    if (allocated(s%eps_zz)) deallocate(s%eps_zz)
    if (allocated(s%eps_yz)) deallocate(s%eps_yz)
  end subroutine strain_result_free

  ! ==================================================================
  ! Private helpers
  ! ==================================================================

  ! Local flat_idx to avoid circular dependency on geometry module.
  pure function wire_flat_idx(nx, ix, iy) result(ij)
    integer, intent(in) :: nx, ix, iy
    integer :: ij
    ij = (iy - 1) * nx + ix
  end function wire_flat_idx

  ! ------------------------------------------------------------------
  ! Add a COO entry (no duplicate search -- stencil produces no dupes).
  ! Duplicates are handled during COO-to-CSR merge sort.
  ! ------------------------------------------------------------------
  subroutine add_coo_entry(rows, cols, vals, capacity, idx, row, col, val)
    integer, intent(inout) :: rows(:), cols(:)
    real(kind=dp), intent(inout) :: vals(:)
    integer, intent(in) :: capacity
    integer, intent(inout) :: idx
    integer, intent(in) :: row, col
    real(kind=dp), intent(in) :: val

    if (abs(val) < 1.0e-30_dp) return

    idx = idx + 1
    if (idx > capacity) then
      print *, 'ERROR: COO capacity exceeded in strain solver'
      stop 1
    end if
    rows(idx) = row
    cols(idx) = col
    vals(idx) = val

  end subroutine add_coo_entry

  ! ------------------------------------------------------------------
  ! Convert COO to CSR format (real-valued, 1-based indexing).
  ! Sorts by (row, col), merges duplicates, produces ia/ja/a_csr
  ! compatible with PARDISO (ia has n+1 entries, 1-based).
  ! ------------------------------------------------------------------
  subroutine coo_to_csr_real(nrows, ncols, nnz_in, rows, cols, vals, &
      ia, ja, a_csr)

    integer, intent(in) :: nrows, ncols, nnz_in
    integer, intent(in) :: rows(nnz_in), cols(nnz_in)
    real(kind=dp), intent(in) :: vals(nnz_in)
    integer, allocatable, intent(out) :: ia(:), ja(:)
    real(kind=dp), allocatable, intent(out) :: a_csr(:)

    integer, allocatable :: idx(:), r_sorted(:), c_sorted(:)
    real(kind=dp), allocatable :: v_sorted(:)
    integer :: i, nnz_final, row

    if (nnz_in == 0) then
      allocate(ia(nrows + 1), ja(0), a_csr(0))
      ia = 1
      return
    end if

    ! Sort COO entries by (row, col)
    allocate(idx(nnz_in))
    do i = 1, nnz_in
      idx(i) = i
    end do
    call merge_sort_coo_real(nnz_in, idx, rows, cols)

    ! Apply permutation and merge duplicates
    allocate(v_sorted(nnz_in), r_sorted(nnz_in), c_sorted(nnz_in))

    v_sorted(1) = vals(idx(1))
    r_sorted(1) = rows(idx(1))
    c_sorted(1) = cols(idx(1))
    nnz_final = 1

    do i = 2, nnz_in
      if (r_sorted(nnz_final) == rows(idx(i)) .and. &
          c_sorted(nnz_final) == cols(idx(i))) then
        v_sorted(nnz_final) = v_sorted(nnz_final) + vals(idx(i))
      else
        nnz_final = nnz_final + 1
        r_sorted(nnz_final) = rows(idx(i))
        c_sorted(nnz_final) = cols(idx(i))
        v_sorted(nnz_final) = vals(idx(i))
      end if
    end do

    ! Build CSR arrays
    allocate(a_csr(nnz_final), ja(nnz_final), ia(nrows + 1))
    a_csr(1:nnz_final) = v_sorted(1:nnz_final)
    ja(1:nnz_final) = c_sorted(1:nnz_final)

    ! Build rowptr
    ia = nnz_final + 1  ! sentinel
    do i = 1, nnz_final
      row = r_sorted(i)
      if (ia(row) > i) ia(row) = i
    end do
    ia(nrows + 1) = nnz_final + 1

    ! Fill empty rows
    do row = nrows, 1, -1
      if (ia(row) == nnz_final + 1) then
        ia(row) = ia(row + 1)
      end if
    end do

    deallocate(idx, v_sorted, r_sorted, c_sorted)

  end subroutine coo_to_csr_real

  ! ------------------------------------------------------------------
  ! Bottom-up merge sort for real COO by (row, col)
  ! ------------------------------------------------------------------
  subroutine merge_sort_coo_real(nnz, idx, rows, cols)
    integer, intent(in) :: nnz
    integer, intent(inout) :: idx(nnz)
    integer, intent(in) :: rows(nnz), cols(nnz)

    integer, allocatable :: work(:)
    integer :: width, start, mid, finish, i, j, k

    if (nnz <= 1) return

    allocate(work(nnz))

    width = 1
    do while (width < nnz)
      start = 1
      do while (start <= nnz)
        mid = min(start + width - 1, nnz)
        finish = min(start + 2*width - 1, nnz)

        if (mid < finish) then
          i = start
          j = mid + 1
          k = start
          do while (i <= mid .and. j <= finish)
            if (rows(idx(i)) < rows(idx(j))) then
              work(k) = idx(i); i = i + 1
            else if (rows(idx(i)) > rows(idx(j))) then
              work(k) = idx(j); j = j + 1
            else
              if (cols(idx(i)) <= cols(idx(j))) then
                work(k) = idx(i); i = i + 1
              else
                work(k) = idx(j); j = j + 1
              end if
            end if
            k = k + 1
          end do

          do while (i <= mid)
            work(k) = idx(i); i = i + 1; k = k + 1
          end do
          do while (j <= finish)
            work(k) = idx(j); j = j + 1; k = k + 1
          end do

          idx(start:finish) = work(start:finish)
        end if

        start = start + 2*width
      end do
      width = width * 2
    end do

    deallocate(work)
  end subroutine merge_sort_coo_real

end module strain_solver
