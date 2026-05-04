module magnetic_field

  use definitions
  use sparse_matrices
  implicit none
  private

  public :: add_zeeman_coo, add_peierls_coo, compute_zeeman_vz

  ! 8x8 spin matrices in zinc-blende basis (HH1,HH2,LH1,LH2,SO1,SO2,CB1,CB2)
  ! Derived from J=3/2 and J=1/2 angular momentum operators
  ! Spin-up diagonal: +mu_B * B, spin-down: -mu_B * B
  ! mu_B = e * hbar / (2 * m0) = 5.7884e-5 eV/T

contains

  subroutine add_zeeman_coo(coo_vals, coo_row, coo_col, nnz_offset, &
                             grid, B_vec, g_factor)
    ! Adds g*mu_B * B . sigma to the 8-band diagonal at each grid point.
    ! Each grid point contributes 8 diagonal COO entries (one per band).
    real(kind=dp), intent(inout) :: coo_vals(:)
    integer, intent(inout) :: coo_row(:), coo_col(:)
    integer, intent(inout) :: nnz_offset
    type(spatial_grid), intent(in) :: grid
    real(kind=dp), intent(in) :: B_vec(3), g_factor

    integer :: i, idx, n
    real(kind=dp) :: Vz(8), B_mag

    B_mag = sqrt(sum(B_vec**2))

    n = grid%npoints()
    do i = 1, n
      call compute_zeeman_vz(g_factor, mu_B, B_mag, Vz)

      do idx = 1, 8
        nnz_offset = nnz_offset + 1
        coo_row(nnz_offset) = (i-1)*8 + idx
        coo_col(nnz_offset) = (i-1)*8 + idx
        coo_vals(nnz_offset) = Vz(idx)
      end do
    end do
  end subroutine add_zeeman_coo

  subroutine add_peierls_coo(coo_vals, coo_row, coo_col, nnz_offset, &
                              grid, B_vec)
    ! Peierls substitution: k -> k - eA/hbar
    ! For Landau gauge A = (0, 0, Bx*y):
    !   kz -> kz - e*Bx*y/hbar  (phase factor e^(ieA.z/hbar))
    ! Only affects off-diagonal entries (hopping terms between sites).
    !
    ! Operates on the main complex COO array so that full complex phase
    ! factors are preserved (unlike the old version which discarded the
    ! imaginary part by storing into a real array).
    complex(kind=dp), intent(inout) :: coo_vals(:)
    integer, intent(in) :: coo_row(:), coo_col(:)
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
        ! Get y-coordinates of the two sites (z array holds y for wire mode)
        ! Grid point indices are 1-based: grid_pt = (row-1)/8 + 1
        ! For wire mode with 8 bands per grid point
        y_i = grid%z((coo_row(idx) - 1) / 8 + 1) * 1.0e-10_dp  ! meters
        y_j = grid%z((coo_col(idx) - 1) / 8 + 1) * 1.0e-10_dp  ! meters

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

  pure subroutine compute_zeeman_vz(g_factor, mu_B_val, B_mag, Vz)
    ! Computes the 8-component Zeeman splitting vector:
    !   Vz(b) = g_factor * g_J(b) * mu_B * B_mag
    ! g_J eigenvalues: HH=-3/2, LH=+1/2, SO=-1/2, CB=+-1
    real(kind=dp), intent(in) :: g_factor, mu_B_val, B_mag
    real(kind=dp), intent(out) :: Vz(8)
    Vz(1:2) = -1.5_dp * g_factor * mu_B_val * B_mag  ! HH (J_z = +/- 3/2)
    Vz(3:4) =  0.5_dp * g_factor * mu_B_val * B_mag  ! LH (J_z = +/- 1/2)
    Vz(5:6) = -0.5_dp * g_factor * mu_B_val * B_mag  ! SO (J_z = +/- 1/2)
    Vz(7)   = -1.0_dp * g_factor * mu_B_val * B_mag  ! CB1 (J_z = -1/2)
    Vz(8)   =  1.0_dp * g_factor * mu_B_val * B_mag  ! CB2 (J_z = +1/2)
  end subroutine compute_zeeman_vz

end module magnetic_field
