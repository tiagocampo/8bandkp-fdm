module magnetic_field

  use definitions
  use sparse_matrices
  implicit none
  private

  public :: add_peierls_coo, compute_zeeman_vz, compute_gauge_shifts

  ! 8x8 spin matrices in zinc-blende basis (HH↑,LH↑,LH↓,HH↓,SO↑,SO↓,CB↓,CB↑)
  ! Winkler ordering: bands 1-4 = valence (HH,LH,LH,HH), 5-6 = SO, 7-8 = CB
  ! Zeeman: Vz(b) = g_eff(b) * mu_B * B_mag, Kramers partners split opposite
  ! mu_B = e * hbar / (2 * m0) = 5.7884e-5 eV/T

contains

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

  pure subroutine compute_zeeman_vz(g_factor, mu_B_val, B_mag, Vz)
    ! Computes the 8-component Zeeman splitting vector:
    !   Vz(b) = g_factor * c_b * mu_B * B_mag
    ! Coefficients c_b are g_J * m_J for each band in Winkler basis:
    !   Band 1 (HH mJ=+3/2): -3/2,  Band 2 (LH mJ=+1/2): +1/2
    !   Band 3 (LH mJ=-1/2): -1/2,  Band 4 (HH mJ=-3/2): +3/2
    !   Band 5 (SO mJ=+1/2): -1/2,  Band 6 (SO mJ=-1/2): +1/2
    !   Band 7 (CB mJ=-1/2): -1,    Band 8 (CB mJ=+1/2): +1
    real(kind=dp), intent(in) :: g_factor, mu_B_val, B_mag
    real(kind=dp), intent(out) :: Vz(8)
    real(kind=dp) :: E0
    E0 = g_factor * mu_B_val * B_mag
    Vz(1) = -1.5_dp * E0  ! HH  (mJ = +3/2)
    Vz(2) =  0.5_dp * E0  ! LH  (mJ = +1/2)
    Vz(3) = -0.5_dp * E0  ! LH  (mJ = -1/2)
    Vz(4) =  1.5_dp * E0  ! HH  (mJ = -3/2)
    Vz(5) = -0.5_dp * E0  ! SO  (mJ = +1/2)
    Vz(6) =  0.5_dp * E0  ! SO  (mJ = -1/2)
    Vz(7) = -1.0_dp * E0  ! CB↓ (mJ = -1/2)
    Vz(8) =  1.0_dp * E0  ! CB↑ (mJ = +1/2)
  end subroutine compute_zeeman_vz

end module magnetic_field
