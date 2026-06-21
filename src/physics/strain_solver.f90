module strain_solver

  ! ==============================================================================
  ! Strain solver for semiconductor heterostructures (Issue #06, ADR 0005).
  !
  ! Owns the Bir-Pikus deformation-potential concern and the strain
  ! block table (single source of truth for the 32-entry band-pair
  ! topology), plus the QW biaxial strain and the top-level dispatcher.
  ! The wire plane-strain Navier-Cauchy PDE was split out to the sibling
  ! module strain_pde.f90 along the concern boundary identified in
  ! ADR 0005; both modules share the strain_result type via the leaf
  ! module strain_types.f90 (no circular `use`).
  !
  ! Supports two strain modes (dispatched by compute_strain):
  !   1. QW biaxial strain (confinement='qw'): simple algebraic calculation
  !      using lattice mismatch and elastic constants (compute_strain_qw).
  !   2. Wire plane-strain PDE (confinement='wire'): Navier-Cauchy elasticity
  !      on 2D cross-section (compute_strain_wire, in strain_pde.f90).
  !
  ! After strain computation, compute_bir_pikus_blocks applies the Pikus-Bir
  ! deformation potentials to produce per-grid-point band-edge shifts that
  ! the Hamiltonian builders read via the strain table + lookup_bp_field.
  ! ==============================================================================

  use definitions, only: dp, paramStruct, spatial_grid, strain_config, &
    & grid_ngrid, SQR3, IU, RQS2, SQR2, SQR3o2, bir_pikus_blocks, bp_scalar
  use strain_types, only: strain_result, strain_result_free
  use strain_pde, only: compute_strain_wire

  implicit none

  ! Module-level flag to warn about invalid material_id only once.
  ! Encapsulated in this concern (QW biaxial + Bir-Pikus); the wire PDE
  ! concern in strain_pde.f90 has its own warn-once flag. The worst-case
  ! behavioural change vs. the pre-split single flag is that an invalid
  ! material_id in wire mode emits two warning lines instead of one.
  logical, save :: warned_invalid_mat = .false.

  private
  ! Re-export the strain_result type + manual free from strain_types so
  ! existing `use strain_solver, only: strain_result, strain_result_free`
  ! importers (main.f90, simulation_setup.f90, tests) are unchanged.
  public :: strain_result
  public :: strain_result_free
  public :: compute_strain
  public :: bir_pikus_blocks_free
  public :: compute_bir_pikus_blocks
  public :: compute_bp_scalar
  public :: strain_entry
  public :: build_strain_table
  public :: get_strain_table
  public :: init_strain_cache
  public :: lookup_bp_field

  ! ------------------------------------------------------------------
  ! Strain COO insertion table entry: describes one of the 32 block
  ! patterns that map Bir-Pikus fields into the 8x8 band structure.
  ! row_band/col_band are 0-based band offsets (0-7).
  ! field_id selects which bp component to read:
  !   1=delta_EHH, 2=delta_ELH, 3=delta_ESO, 4=delta_Ec,
  !   5=S_eps, 6=R_eps, 7=QT2_eps
  ! prefactor scales the value; use_conjg applies conjg() to the field.
  ! ------------------------------------------------------------------
  type :: strain_entry
    integer               :: row_band, col_band  ! 0-based band offsets (0-7)
    integer               :: field_id            ! 1=EHH,2=ELH,3=ESO,4=Ec,5=S,6=R,7=QT2
    complex(kind=dp)      :: prefactor           ! complex scale (e.g. IU*RQS2)
    logical               :: use_conjg           ! conjugate the field value?
  end type strain_entry

  ! Module-level cache for the strain table (32-entry compile-time constant)
  logical, save :: strain_table_cached = .false.
  type(strain_entry), save :: strain_table_cache(32)

contains

  ! ==================================================================
  ! Build the 32-entry strain insertion table.
  !
  ! Each entry describes how one band-block of the strain Hamiltonian
  ! is constructed from the Bir-Pikus fields.  The table encodes all
  ! diagonal shifts, S_eps, R_eps, and VB-SO coupling terms.
  ! ==================================================================
  function build_strain_table() result(table)
    type(strain_entry) :: table(32)

    ! --- Diagonal entries (1-8) ---
    table( 1) = strain_entry(0, 0, 1, cmplx(1.0_dp, 0.0_dp, kind=dp), .false.)
    table( 2) = strain_entry(1, 1, 2, cmplx(1.0_dp, 0.0_dp, kind=dp), .false.)
    table( 3) = strain_entry(2, 2, 2, cmplx(1.0_dp, 0.0_dp, kind=dp), .false.)
    table( 4) = strain_entry(3, 3, 1, cmplx(1.0_dp, 0.0_dp, kind=dp), .false.)
    table( 5) = strain_entry(4, 4, 3, cmplx(1.0_dp, 0.0_dp, kind=dp), .false.)
    table( 6) = strain_entry(5, 5, 3, cmplx(1.0_dp, 0.0_dp, kind=dp), .false.)
    table( 7) = strain_entry(6, 6, 4, cmplx(1.0_dp, 0.0_dp, kind=dp), .false.)
    table( 8) = strain_entry(7, 7, 4, cmplx(1.0_dp, 0.0_dp, kind=dp), .false.)

    ! --- S_eps entries (9-12) ---
    table( 9) = strain_entry(0, 1, 5, cmplx(1.0_dp, 0.0_dp, kind=dp), .true.)
    table(10) = strain_entry(1, 0, 5, cmplx(1.0_dp, 0.0_dp, kind=dp), .false.)
    table(11) = strain_entry(2, 3, 5, cmplx(-1.0_dp, 0.0_dp, kind=dp), .true.)
    table(12) = strain_entry(3, 2, 5, cmplx(-1.0_dp, 0.0_dp, kind=dp), .false.)

    ! --- R_eps entries (13-16) ---
    table(13) = strain_entry(0, 2, 6, cmplx(1.0_dp, 0.0_dp, kind=dp), .true.)
    table(14) = strain_entry(2, 0, 6, cmplx(1.0_dp, 0.0_dp, kind=dp), .false.)
    table(15) = strain_entry(1, 3, 6, cmplx(1.0_dp, 0.0_dp, kind=dp), .true.)
    table(16) = strain_entry(3, 1, 6, cmplx(1.0_dp, 0.0_dp, kind=dp), .false.)

    ! --- VB-SO coupling entries (17-32) ---
    table(17) = strain_entry(0, 4, 5, -IU * RQS2,   .true.)
    table(18) = strain_entry(4, 0, 5,  IU * RQS2,   .false.)
    table(19) = strain_entry(0, 5, 6,  IU * SQR2,   .true.)
    table(20) = strain_entry(5, 0, 6, -IU * SQR2,   .false.)
    table(21) = strain_entry(1, 4, 7,  IU * RQS2,   .false.)
    table(22) = strain_entry(4, 1, 7, -IU * RQS2,   .false.)
    table(23) = strain_entry(1, 5, 5, -IU * SQR3o2, .true.)
    table(24) = strain_entry(5, 1, 5,  IU * SQR3o2, .false.)
    table(25) = strain_entry(2, 4, 5,  IU * SQR3o2, .false.)
    table(26) = strain_entry(4, 2, 5, -IU * SQR3o2, .true.)
    table(27) = strain_entry(2, 5, 7,  IU * RQS2,   .false.)
    table(28) = strain_entry(5, 2, 7, -IU * RQS2,   .false.)
    table(29) = strain_entry(3, 4, 6, -IU * SQR2,   .false.)
    table(30) = strain_entry(4, 3, 6,  IU * SQR2,   .true.)
    table(31) = strain_entry(3, 5, 5,  IU * RQS2,   .false.)
    table(32) = strain_entry(5, 3, 5, -IU * RQS2,   .true.)

  end function build_strain_table

  function get_strain_table() result(table)
    type(strain_entry) :: table(32)
    if (.not. strain_table_cached) then
      strain_table_cache = build_strain_table()
      strain_table_cached = .true.
    end if
    table = strain_table_cache
  end function get_strain_table

  ! ==================================================================
  ! Pre-initialize the strain table cache (thread-safety).
  ! Call before OpenMP fork to avoid races on the SAVE cache.
  ! ==================================================================
  subroutine init_strain_cache()
    type(strain_entry) :: dummy(32)
    dummy = get_strain_table()
  end subroutine init_strain_cache

  ! ==================================================================
  ! Pure function: look up a Bir-Pikus field value by field_id.
  !
  ! Encapsulates the 7-way select case mapping from strain table
  ! field_id to the corresponding bp component. Used by both
  ! hamiltonianConstructor (dense) and hamiltonian_wire (COO) to
  ! avoid duplicating this dispatch logic.
  !
  ! field_id:  1=delta_EHH, 2=delta_ELH, 3=delta_ESO, 4=delta_Ec,
  !            5=S_eps, 6=R_eps, 7=QT2_eps
  ! ==================================================================
  pure function lookup_bp_field(bp, field_id, ii) result(field_val)
    type(bir_pikus_blocks), intent(in) :: bp
    integer, intent(in) :: field_id, ii
    complex(kind=dp) :: field_val

    select case (field_id)
    case (1); field_val = cmplx(bp%delta_EHH(ii), 0.0_dp, kind=dp)
    case (2); field_val = cmplx(bp%delta_ELH(ii), 0.0_dp, kind=dp)
    case (3); field_val = cmplx(bp%delta_ESO(ii), 0.0_dp, kind=dp)
    case (4); field_val = cmplx(bp%delta_Ec(ii), 0.0_dp, kind=dp)
    case (5); field_val = bp%S_eps(ii)
    case (6); field_val = bp%R_eps(ii)
    case (7); field_val = cmplx(bp%QT2_eps(ii), 0.0_dp, kind=dp)
    case default
      ! Cannot error stop in a pure function — return zero.
      ! The caller validates field_id at the call site.
      field_val = cmplx(0.0_dp, 0.0_dp, kind=dp)
    end select
  end function lookup_bp_field

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
  ! QW biaxial strain (confinement='qw').
  !
  ! For each layer with lattice constant a0, relative to reference a0_ref:
  !   eps_0 = (a0_ref - a0) / a0
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

      eps_0 = (a0_ref - a0_mat) / a0_mat

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
  ! Compute full Bir-Pikus strain Hamiltonian components per grid point.
  !
  ! Returns a bir_pikus_blocks type with per-band diagonal shifts and
  ! off-diagonal coupling terms, ready for insertion into the 8-band
  ! Hamiltonian by ZB8bandQW or ZB8bandBulk (hamiltonianConstructor.f90).
  !
  ! The diagonal shifts are:
  !   HH (bands 1,4): -P_eps + Q_eps
  !   LH (bands 2,3): -P_eps - Q_eps
  !   SO (bands 5,6): -P_eps
  !   CB (bands 7,8): ac * Tr(eps)
  !
  ! where P_eps = -av * Tr(eps), Q_eps = -(b/2) * (ezz - 0.5*(exx+eyy)).
  !
  ! Off-diagonal terms (same structure as k-dependent R, S in bulk):
  !   R_eps = -sqrt(3) * [b/2*(exx-eyy) - i*d*exy]
  !   S_eps =  i*2*sqrt(3) * d * (exz - i*eyz)
  !
  ! Note: for [001] biaxial strain, exx=eyy and exy=exz=eyz=0,
  ! so R_eps=0, S_eps=0, but QT2_eps = 2*Q_eps != 0 (LH-SO coupling).
  !
  ! eps_xz is not stored in strain_result (only xx,yyzz,yz for QW).
  ! For QW (ndim=1), eps_xz = eps_xy = 0 (pure [001] biaxial).
  ! For wire (ndim=2), we approximate eps_xz ~ 0 (plane strain in xy).
  ! ==================================================================
  subroutine compute_bir_pikus_blocks(strain_out, params, material_id, grid, bp)

    type(strain_result), intent(in)    :: strain_out
    type(paramStruct), intent(in)      :: params(:)
    integer, intent(in)                :: material_id(:)
    type(spatial_grid), intent(in)     :: grid
    type(bir_pikus_blocks), intent(out) :: bp

    integer :: ngrid, ij, mid
    real(kind=dp) :: Tr_eps
    real(kind=dp) :: eps_xx, eps_yy, eps_zz, eps_xy, eps_xz, eps_yz
    type(bp_scalar) :: s

    ngrid = grid_ngrid(grid)

    ! Early exit if no strain data
    if (.not. allocated(strain_out%eps_xx)) return
    if (size(strain_out%eps_xx) < ngrid) return

    ! Allocate all arrays
    allocate(bp%delta_Ec(ngrid))
    allocate(bp%delta_EHH(ngrid))
    allocate(bp%delta_ELH(ngrid))
    allocate(bp%delta_ESO(ngrid))
    allocate(bp%R_eps(ngrid))
    allocate(bp%S_eps(ngrid))
    allocate(bp%QT2_eps(ngrid))

    bp%delta_Ec  = 0.0_dp
    bp%delta_EHH = 0.0_dp
    bp%delta_ELH = 0.0_dp
    bp%delta_ESO = 0.0_dp
    bp%R_eps     = cmplx(0.0_dp, 0.0_dp, kind=dp)
    bp%S_eps     = cmplx(0.0_dp, 0.0_dp, kind=dp)
    bp%QT2_eps   = 0.0_dp

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

      eps_xx = strain_out%eps_xx(ij)
      eps_yy = strain_out%eps_yy(ij)
      eps_zz = strain_out%eps_zz(ij)
      eps_yz = strain_out%eps_yz(ij)
      eps_xy = 0.0_dp  ! not stored; zero for [001] biaxial
      eps_xz = 0.0_dp  ! not stored; zero for [001] biaxial

      Tr_eps = eps_xx + eps_yy + eps_zz

      ! Check if strain is negligible
      if (abs(Tr_eps) < 1.0e-20_dp .and. abs(eps_yz) < 1.0e-20_dp) cycle

      ! Single source of truth for Bir-Pikus formulas
      s = compute_bp_scalar(params(mid), eps_xx, eps_yy, eps_zz, &
                            eps_xy, eps_xz, eps_yz)
      bp%delta_Ec(ij)  = s%delta_Ec
      bp%delta_EHH(ij) = s%delta_EHH
      bp%delta_ELH(ij) = s%delta_ELH
      bp%delta_ESO(ij) = s%delta_ESO
      bp%R_eps(ij)     = s%R_eps
      bp%S_eps(ij)     = s%S_eps
      bp%QT2_eps(ij)   = s%QT2_eps
    end do

  end subroutine compute_bir_pikus_blocks

  ! ==================================================================
  ! Deallocate bir_pikus_blocks arrays
  ! ==================================================================
  subroutine bir_pikus_blocks_free(bp)
    type(bir_pikus_blocks), intent(inout) :: bp

    if (allocated(bp%delta_Ec))  deallocate(bp%delta_Ec)
    if (allocated(bp%delta_EHH)) deallocate(bp%delta_EHH)
    if (allocated(bp%delta_ELH)) deallocate(bp%delta_ELH)
    if (allocated(bp%delta_ESO)) deallocate(bp%delta_ESO)
    if (allocated(bp%R_eps))     deallocate(bp%R_eps)
    if (allocated(bp%S_eps))     deallocate(bp%S_eps)
    if (allocated(bp%QT2_eps))   deallocate(bp%QT2_eps)
  end subroutine bir_pikus_blocks_free

  ! ==================================================================
  ! Pure function: single source of truth for Bir-Pikus formulas.
  !
  ! Computes all diagonal and off-diagonal strain Hamiltonian components
  ! for a single grid point from the strain tensor and material parameters.
  ! ==================================================================
  elemental pure function compute_bp_scalar(params, eps_xx, eps_yy, eps_zz, &
      eps_xy, eps_xz, eps_yz) result(s)
    type(paramStruct), intent(in) :: params
    real(kind=dp), intent(in) :: eps_xx, eps_yy, eps_zz
    real(kind=dp), intent(in) :: eps_xy, eps_xz, eps_yz
    type(bp_scalar) :: s

    real(kind=dp) :: Tr_eps, P_eps, Q_eps, T_eps

    Tr_eps = eps_xx + eps_yy + eps_zz
    P_eps = -params%av * Tr_eps
    Q_eps = -params%b_dp * 0.5_dp * (eps_zz - 0.5_dp * (eps_yy + eps_xx))
    T_eps = -Q_eps

    s%delta_Ec  = params%ac * Tr_eps
    s%delta_EHH = -P_eps + Q_eps
    s%delta_ELH = -P_eps - Q_eps
    s%delta_ESO = -P_eps
    s%QT2_eps   = Q_eps - T_eps

    s%R_eps = -SQR3 * (params%b_dp * 0.5_dp * (eps_xx - eps_yy) &
                       - IU * params%d_dp * eps_xy)
    s%S_eps = IU * 2.0_dp * SQR3 * params%d_dp * &
              cmplx(eps_xz, -eps_yz, kind=dp)

  end function compute_bp_scalar

end module strain_solver
