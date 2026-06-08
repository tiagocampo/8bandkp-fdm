module magnetic_field

  use definitions, only: dp, e, hbar, mu_B, spatial_grid
  use sparse_matrices
  implicit none
  private

  public :: add_peierls_coo, compute_zeeman_vz, compute_gauge_shifts
  public :: zeeman_entry, get_zeeman_table, init_zeeman_cache

  ! 8x8 spin matrices in zinc-blende basis (HH↑,LH↑,LH↓,HH↓,SO↑,SO↓,CB↓,CB↑)
  ! Winkler ordering: bands 1-4 = valence (HH,LH,LH,HH), 5-6 = SO, 7-8 = CB
  ! Zeeman: Vz(b) = g_eff(b) * mu_B * B_mag, Kramers partners split opposite
  ! mu_B = e * hbar / (2 * m0) = 5.7884e-5 eV/T

  ! ------------------------------------------------------------------
  ! Zeeman diagonal table entry: describes the g_multiplier for one band.
  ! band_index is 0-based (0-7). The Zeeman energy for band b is:
  !   Vz(b) = g_factor * g_multiplier(b) * mu_B * B_mag
  ! Coefficients from compute_zeeman_vz (magnetic_field.f90):
  !   Band 0 (HH mJ=+3/2): -1.5,  Band 1 (LH mJ=+1/2): +0.5
  !   Band 2 (LH mJ=-1/2): -0.5,  Band 3 (HH mJ=-3/2): +1.5
  !   Band 4 (SO mJ=+1/2): -0.5,  Band 5 (SO mJ=-1/2): +0.5
  !   Band 6 (CB mJ=-1/2): -1.0,  Band 7 (CB mJ=+1/2): +1.0
  ! ------------------------------------------------------------------
  type :: zeeman_entry
    integer               :: band_index     ! 0-based band offset (0-7)
    real(kind=dp)         :: g_multiplier   ! spin-dependent coefficient c_b
  end type zeeman_entry

  ! Module-level cache for the Zeeman table (8-entry compile-time constant)
  logical, save :: zeeman_table_cached = .false.
  type(zeeman_entry), save :: zeeman_table_cache(8)

contains

  ! ==================================================================
  ! Build the 8-entry Zeeman diagonal table.
  !
  ! Each entry maps band_index to its g_multiplier coefficient c_b
  ! such that the Zeeman energy for band b is:
  !   Vz(b) = g_factor * c_b * mu_B * B_mag
  !
  ! Coefficients are g_J * m_J for each band in the Winkler basis
  ! (single source of truth for Zeeman diagonal g-multipliers).
  ! ==================================================================
  function build_zeeman_table() result(table)
    type(zeeman_entry) :: table(8)

    table(1) = zeeman_entry(0, -1.5_dp)  ! HH  (mJ = +3/2)
    table(2) = zeeman_entry(1,  0.5_dp)  ! LH  (mJ = +1/2)
    table(3) = zeeman_entry(2, -0.5_dp)  ! LH  (mJ = -1/2)
    table(4) = zeeman_entry(3,  1.5_dp)  ! HH  (mJ = -3/2)
    table(5) = zeeman_entry(4, -0.5_dp)  ! SO  (mJ = +1/2)
    table(6) = zeeman_entry(5,  0.5_dp)  ! SO  (mJ = -1/2)
    table(7) = zeeman_entry(6, -1.0_dp)  ! CB  (mJ = -1/2)
    table(8) = zeeman_entry(7,  1.0_dp)  ! CB  (mJ = +1/2)

  end function build_zeeman_table

  function get_zeeman_table() result(table)
    type(zeeman_entry) :: table(8)
    if (.not. zeeman_table_cached) then
      zeeman_table_cache = build_zeeman_table()
      zeeman_table_cached = .true.
    end if
    table = zeeman_table_cache
  end function get_zeeman_table

  ! ==================================================================
  ! Pre-initialize the Zeeman table cache (thread-safety).
  ! Call before OpenMP fork to avoid races on the SAVE cache.
  ! ==================================================================
  subroutine init_zeeman_cache()
    type(zeeman_entry) :: dummy(8)
    dummy = get_zeeman_table()
  end subroutine init_zeeman_cache

  subroutine add_peierls_coo(coo_vals, coo_row, coo_col, nnz_offset, &
                              grid, B_vec)
    ! Peierls substitution: k -> k - eA/hbar
    ! For Landau gauge A = (0, 0, Bx*y):
    !   kz -> kz - e*Bx*y/hbar  (phase factor e^(ieA.z/hbar))
    ! Only affects off-diagonal entries (hopping terms between sites).
    !
    ! NOTE: For multi-column wires (nx > 1) with Bx != 0, this applies
    ! phases to transverse y-direction hops which have zero z-displacement
    ! in this gauge. The correct approach is a position-dependent kz shift
    ! in Hamiltonian assembly. Single-column wires (nx = 1) are unaffected
    ! since y_i = y_j for all entries.
    !
    ! Operates on the main complex COO array so that full complex phase
    ! factors are preserved (unlike the old version which discarded the
    ! imaginary part by storing into a real array).
    complex(kind=dp), intent(inout), contiguous :: coo_vals(:)
    integer, intent(in), contiguous :: coo_row(:), coo_col(:)
    integer, intent(in) :: nnz_offset
    type(spatial_grid), intent(in) :: grid
    real(kind=dp), intent(in) :: B_vec(3)

    real(kind=dp) :: Bx, phase, y_i, y_j, dy_m, hbar_J
    integer :: idx, ngrid
    complex(kind=dp) :: exp_phase

    Bx = B_vec(1)

    ! Return early if no magnetic field in x-direction
    if (abs(Bx) < 1e-12_dp) return

    ngrid = grid%npoints()

    ! Convert dy from Angstrom to meters for SI units
    dy_m = grid%dy * 1.0e-10_dp

    ! Convert hbar from eV*s to J*s: hbar_J = hbar_eV * e_C
    ! This gives hbar in J*s for proper SI unit handling
    hbar_J = hbar * e

    ! Apply Peierls phase to all off-diagonal entries
    do idx = 1, nnz_offset
      if (coo_row(idx) /= coo_col(idx)) then
        block
          integer :: flat_i, flat_j, iy_i, iy_j

          ! Band-major basis: grid_pt = mod(row-1, ngrid) + 1
          flat_i = mod(coo_row(idx) - 1, ngrid) + 1
          flat_j = mod(coo_col(idx) - 1, ngrid) + 1

          ! Extract y-coordinate for each site. In 2D wire mode, coords(2,:)
          ! contains the y-coordinate for every flat spatial site. If coords
          ! is unavailable, grid%z has length ny, so convert flat site -> iy.
          if (grid%ndim == 2 .and. allocated(grid%coords)) then
            y_i = grid%coords(2, flat_i) * 1.0e-10_dp
            y_j = grid%coords(2, flat_j) * 1.0e-10_dp
          else
            iy_i = (flat_i - 1) / grid%nx + 1
            iy_j = (flat_j - 1) / grid%nx + 1
            y_i = grid%z(iy_i) * 1.0e-10_dp
            y_j = grid%z(iy_j) * 1.0e-10_dp
          end if
        end block

        ! Peierls phase: phi = e * Bx * (y_i - y_j) * dz / hbar
        ! Using SI units: e in C, B in T, y and dz in m, hbar in J*s
        phase = e * Bx * (y_i - y_j) * dy_m / hbar_J

        ! Phase factor: exp(-i * phase)
        exp_phase = cmplx(cos(phase), -sin(phase), kind=dp)

        ! Multiply the off-diagonal entry by the full complex phase factor
        coo_vals(idx) = coo_vals(idx) * exp_phase
      end if
    end do
  end subroutine add_peierls_coo

  pure subroutine compute_gauge_shifts(x_grid, B_vec, ky, kz, Pi_y, Pi_z)
    ! Landau gauge A = (0, Bz*x, -By*x) for B = (Bx, By, Bz).
    ! Canonical momentum shifts in 1/AA:
    !   Pi_y(i) = ky + Bz * x(i) * 1e-20 / hbar
    !   Pi_z(i) = kz - By * x(i) * 1e-20 / hbar
    real(kind=dp), intent(in), contiguous :: x_grid(:)
    real(kind=dp), intent(in) :: B_vec(3), ky, kz
    real(kind=dp), intent(out), contiguous :: Pi_y(:), Pi_z(:)

    real(kind=dp), parameter :: inv_hbar_AA = 1.0e-20_dp / hbar  ! 1/(eV*s*T*AA) -> shift factor
    integer :: i

    do i = 1, size(x_grid)
      Pi_y(i) = ky + B_vec(3) * x_grid(i) * inv_hbar_AA
      Pi_z(i) = kz - B_vec(2) * x_grid(i) * inv_hbar_AA
    end do
  end subroutine compute_gauge_shifts

  subroutine compute_zeeman_vz(g_factor, mu_B_val, B_mag, Vz)
    ! Computes the 8-component Zeeman splitting vector by reading from the
    ! Zeeman table (single source of truth). Each band gets:
    !   Vz(b) = g_factor * g_multiplier(b) * mu_B * B_mag
    ! Note: no longer pure — reads from SAVE cache via get_zeeman_table().
    ! The cache is pre-warmed by init_zeeman_cache() before OpenMP fork.
    real(kind=dp), intent(in) :: g_factor, mu_B_val, B_mag
    real(kind=dp), intent(out) :: Vz(8)
    type(zeeman_entry) :: ztable(8)
    integer :: b
    real(kind=dp) :: E0

    E0 = g_factor * mu_B_val * B_mag
    ztable = get_zeeman_table()
    do b = 1, 8
      Vz(b) = ztable(b)%g_multiplier * E0
    end do
  end subroutine compute_zeeman_vz

end module magnetic_field
