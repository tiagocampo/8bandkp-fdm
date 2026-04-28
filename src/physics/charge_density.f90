module charge_density
  ! Compute electron and hole charge densities from 8-band k.p eigenstates.
  !
  ! QW mode:
  !   n(z) = sum_s int_0^{k_max} |Psi_s(z,k_par)|^2 f(E_s(k_par)-mu,T)
  !          * k_par/(2*pi) dk_par   (sum over CB subbands, spin degeneracy 2)
  !
  !   |Psi_s(z)|^2 = sum_{nu=1}^{8} |Psi_s^nu(z)|^2  (all 8 components)
  !
  ! Bulk mode:
  !   n = (2/(2*pi)^3) sum_{k,s in CB} f(E_s(k)-mu,T) * (4*pi*k^2 dk)
  !   (3D DOS with spherical shell volume)

  use definitions, only: dp, pi_dp, kB_eV
  use utils, only: simpson, simpson_real
  implicit none

  private
  public :: compute_charge_density_qw
  public :: compute_charge_density_bulk
  public :: compute_charge_density_wire
  public :: fermi_dirac

contains

  ! ------------------------------------------------------------------
  ! Fermi-Dirac distribution function
  ! ------------------------------------------------------------------
  pure function fermi_dirac(energy, mu, T) result(f)
    real(kind=dp), intent(in) :: energy  ! Energy (eV)
    real(kind=dp), intent(in) :: mu      ! Fermi level (eV)
    real(kind=dp), intent(in) :: T       ! Temperature (K)
    real(kind=dp) :: f

    real(kind=dp) :: x

    if (T < 1.0e-30_dp) then
      ! T -> 0 limit: step function
      if (energy < mu) then
        f = 1.0_dp
      else
        f = 0.0_dp
      end if
    else
      x = (energy - mu) / (kB_eV * T)
      if (x > 500.0_dp) then
        f = 0.0_dp
      else if (x < -500.0_dp) then
        f = 1.0_dp
      else
        f = 1.0_dp / (1.0_dp + exp(x))
      end if
    end if
  end function fermi_dirac


  ! ------------------------------------------------------------------
  ! QW charge density
  ! ------------------------------------------------------------------
  subroutine compute_charge_density_qw(n_electron, n_hole, eigenvectors, &
      & eigenvalues_kpar, kpar_grid, fermi_level, temperature, N, num_subbands, &
      & num_kpar, numcb)
    ! Compute position-dependent electron and hole densities for a QW.
    !
    ! For each subband s and each k_par point:
    !   |Psi_s(z, k_par)|^2 = sum_{nu=1}^{8} |Psi_s^nu(z, k_par)|^2
    !
    ! Then integrate over k_par with cylindrical weight k_par/(2*pi):
    !   n(z) = 2 * sum_{s in CB} simpson( [|Psi_s|^2 * f * k_par/(2*pi)] )
    !
    ! Spin degeneracy factor 2 is included explicitly.
    !
    ! Subbands are sorted by eigenvalue at k_par=0: the highest numcb are CB,
    ! the rest are VB.
    !
    ! Output is in cm^-3 (converted from 1/nm^3 via 1e21 factor).

    integer, intent(in)          :: N
    integer, intent(in)          :: num_subbands
    integer, intent(in)          :: num_kpar
    integer, intent(in)          :: numcb
    real(kind=dp), intent(out)   :: n_electron(N)
    real(kind=dp), intent(out)   :: n_hole(N)
    complex(kind=dp), intent(in) :: eigenvectors(8*N, num_subbands, num_kpar)
    real(kind=dp), intent(in)    :: eigenvalues_kpar(num_subbands, num_kpar)
    real(kind=dp), intent(in)    :: kpar_grid(num_kpar)
    real(kind=dp), intent(in)    :: fermi_level
    real(kind=dp), intent(in)    :: temperature

    ! Local arrays
    real(kind=dp), allocatable :: psi2_z(:,:)    ! |Psi|^2 at each z, k_par
    real(kind=dp), allocatable     :: integrand(:)  ! for Simpson (real-valued)
    real(kind=dp), allocatable :: sorted_evals(:)
    integer, allocatable :: sort_idx(:)
    real(kind=dp) :: k_max
    integer :: nk, s, iz, j, band_start, idx
    integer :: numcb_local

    ! Ensure odd number of k_par points for Simpson
    nk = num_kpar
    if (mod(nk, 2) == 0) nk = nk - 1

    k_max = kpar_grid(nk)

    n_electron = 0.0_dp
    n_hole = 0.0_dp

    ! --- Classify subbands at first k_par point (k_par = 0) ---
    ! Sort eigenvalues at k_par=1 by descending value; top numcb are CB
    allocate(sorted_evals(num_subbands))
    allocate(sort_idx(num_subbands))

    do s = 1, num_subbands
      sorted_evals(s) = eigenvalues_kpar(s, 1)
      sort_idx(s) = s
    end do

    ! Simple insertion sort (descending)
    call sort_descending(sorted_evals, sort_idx, num_subbands)

    numcb_local = min(numcb, num_subbands)

    allocate(psi2_z(N, nk))
    allocate(integrand(nk))

    ! --- Electron density: CB subbands (highest numcb eigenvalues) ---
    do j = 1, numcb_local
      s = sort_idx(j)  ! original subband index

      call compute_psi2(psi2_z, eigenvectors, s, N, nk, num_subbands)

      call accumulate_band_density(n_electron, psi2_z, eigenvalues_kpar, &
        & kpar_grid, integrand, fermi_level, temperature, s, N, nk, &
        & num_subbands, k_max, is_cb=.true.)
    end do

    ! --- Hole density: VB subbands (lowest numvb eigenvalues) ---
    do j = numcb_local + 1, num_subbands
      s = sort_idx(j)  ! original subband index

      call compute_psi2(psi2_z, eigenvectors, s, N, nk, num_subbands)

      call accumulate_band_density(n_hole, psi2_z, eigenvalues_kpar, &
        & kpar_grid, integrand, fermi_level, temperature, s, N, nk, &
        & num_subbands, k_max, is_cb=.false.)
    end do

    ! --- Convert 1/AA^3 to cm^-3 ---
    n_electron = n_electron * 1.0e24_dp
    n_hole = n_hole * 1.0e24_dp

    deallocate(psi2_z, integrand, sorted_evals, sort_idx)

  end subroutine compute_charge_density_qw


  ! ------------------------------------------------------------------
  ! Bulk charge density
  ! ------------------------------------------------------------------
  subroutine compute_charge_density_bulk(n_electron, n_hole, eigenvalues_k, &
      & k_grid, fermi_level, temperature, num_subbands, num_k)
    ! Compute scalar electron and hole densities for bulk.
    !
    ! Uses 3D DOS via explicit k-sampling with spherical shell volume:
    !   n = (2/(2*pi)^3) sum_{k,s in CB} f(E_s(k)-mu, T) * 4*pi*k^2 * dk
    !
    ! Factor 2 for spin degeneracy. Output in cm^-3.

    real(kind=dp), intent(out)   :: n_electron
    real(kind=dp), intent(out)   :: n_hole
    integer, intent(in)          :: num_subbands
    integer, intent(in)          :: num_k
    real(kind=dp), intent(in)    :: eigenvalues_k(num_subbands, num_k)
    real(kind=dp), intent(in)    :: k_grid(num_k)
    real(kind=dp), intent(in)    :: fermi_level
    real(kind=dp), intent(in)    :: temperature

    integer :: s, ik
    real(kind=dp) :: dk, shell_vol, prefactor, occ

    n_electron = 0.0_dp
    n_hole = 0.0_dp

    if (num_k < 2) return

    dk = k_grid(2) - k_grid(1)
    prefactor = 2.0_dp / (2.0_dp * pi_dp)**3

    do ik = 1, num_k
      shell_vol = 4.0_dp * pi_dp * k_grid(ik)**2 * dk
      do s = 1, num_subbands
        occ = fermi_dirac(eigenvalues_k(s, ik), fermi_level, temperature)
        ! Bulk: subbands 1-6 are VB, 7-8 are CB
        ! But num_subbands is typically 8 for bulk
        if (s >= 7) then
          ! CB subbands
          n_electron = n_electron + prefactor * occ * shell_vol
        else
          ! VB subbands: holes = absence of electrons, so (1 - f)
          n_hole = n_hole + prefactor * (1.0_dp - occ) * shell_vol
        end if
      end do
    end do

    ! Convert 1/AA^3 to cm^-3
    n_electron = n_electron * 1.0e24_dp
    n_hole = n_hole * 1.0e24_dp

  end subroutine compute_charge_density_bulk


  ! ------------------------------------------------------------------
  ! Helper: insertion sort (descending)
  ! ------------------------------------------------------------------
  subroutine sort_descending(vals, idx, n)
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: vals(n)
    integer, intent(inout) :: idx(n)

    integer :: i, j, key_idx
    real(kind=dp) :: key_val

    do i = 2, n
      key_val = vals(i)
      key_idx = idx(i)
      j = i - 1
      do while (j >= 1)
        if (vals(j) >= key_val) exit
        vals(j + 1) = vals(j)
        idx(j + 1) = idx(j)
        j = j - 1
      end do
      vals(j + 1) = key_val
      idx(j + 1) = key_idx
    end do
  end subroutine sort_descending


  ! ------------------------------------------------------------------
  ! Compute |Psi_s(z, k_par)|^2 summed over all 8 band components
  ! ------------------------------------------------------------------
  subroutine compute_psi2(psi2_z, eigenvectors, s, N, nk, num_subbands)
    integer, intent(in)          :: s, N, nk, num_subbands
    real(kind=dp), intent(out)   :: psi2_z(N, nk)
    complex(kind=dp), intent(in) :: eigenvectors(8*N, num_subbands, nk)

    integer :: idx, band_start, iz

    psi2_z = 0.0_dp
    do idx = 1, nk
      do band_start = 1, 8
        do iz = 1, N
          psi2_z(iz, idx) = psi2_z(iz, idx) &
            & + abs(eigenvectors((band_start - 1)*N + iz, s, idx))**2
        end do
      end do
    end do

  end subroutine compute_psi2


  ! ------------------------------------------------------------------
  ! Accumulate band density at each z via k_par integration
  ! ------------------------------------------------------------------
  subroutine accumulate_band_density(n_density, psi2_z, eigenvalues_kpar, &
      & kpar_grid, integrand, fermi_level, temperature, s, N, nk, num_subbands, &
      & k_max, is_cb)
    integer, intent(in)          :: N, nk, num_subbands
    real(kind=dp), intent(inout) :: n_density(N)
    real(kind=dp), intent(in)    :: psi2_z(N, nk)
    real(kind=dp), intent(in)    :: eigenvalues_kpar(num_subbands, nk)
    real(kind=dp), intent(in)    :: kpar_grid(nk)
    real(kind=dp), intent(inout) :: integrand(nk)
    real(kind=dp), intent(in)    :: fermi_level, temperature, k_max
    integer, intent(in)          :: s
    logical, intent(in)          :: is_cb

    integer :: iz, idx
    real(kind=dp) :: occ(nk)

    ! Pre-compute weighted occupation (independent of z)
    do idx = 1, nk
      occ(idx) = fermi_dirac(eigenvalues_kpar(s, idx), fermi_level, temperature)
      if (.not. is_cb) occ(idx) = 1.0_dp - occ(idx)
      occ(idx) = occ(idx) * kpar_grid(idx) / (2.0_dp * pi_dp)
    end do

    do iz = 1, N
      integrand(1:nk) = psi2_z(iz, 1:nk) * occ(1:nk)
      n_density(iz) = n_density(iz) + 2.0_dp * simpson_real(integrand, 0.0_dp, k_max)
    end do

  end subroutine accumulate_band_density

  ! ------------------------------------------------------------------
  ! Wire charge density (2D confinement in y-z plane, free along x)
  ! ------------------------------------------------------------------
  subroutine compute_charge_density_wire(n_electron, n_hole, eigenvectors, &
      & eigenvalues_kx, kx_grid, fermi_level, temperature, Ny, Nz, num_subbands, &
      & num_kx, numcb)
    ! Compute 2D electron and hole densities for a quantum wire.
    !
    ! n(y,z) = 2 * sum_s int dk_x/(2*pi) * |Psi_s(y,z)|^2 * f(E_s(k_x)-EF,T)
    !
    ! s sums over CB subbands for electrons, VB subbands for holes.
    ! |Psi_s(y,z)|^2 = sum_{nu=1}^{8} |Psi_s^nu(y,z)|^2
    ! Factor 2 is spin degeneracy.
    !
    ! Output arrays are flattened (Ny*Nz), in cm^-3.

    integer, intent(in)          :: Ny, Nz, num_subbands, num_kx, numcb
    real(kind=dp), intent(out)   :: n_electron(Ny*Nz)
    real(kind=dp), intent(out)   :: n_hole(Ny*Nz)
    complex(kind=dp), intent(in) :: eigenvectors(8*Ny*Nz, num_subbands, num_kx)
    real(kind=dp), intent(in)    :: eigenvalues_kx(num_subbands, num_kx)
    real(kind=dp), intent(in)    :: kx_grid(num_kx)
    real(kind=dp), intent(in)    :: fermi_level
    real(kind=dp), intent(in)    :: temperature

    integer :: Ngrid
    real(kind=dp), allocatable :: psi2(:,:)       ! |Psi|^2 at each grid pt, kx
    real(kind=dp), allocatable :: integrand(:)     ! for Simpson
    real(kind=dp), allocatable :: sorted_evals(:)
    integer, allocatable :: sort_idx(:)
    real(kind=dp) :: kx_max
    integer :: nk, s, p, j
    integer :: numcb_local

    Ngrid = Ny * Nz

    ! Use odd number of kx points for Simpson
    nk = num_kx
    if (mod(nk, 2) == 0) nk = nk - 1

    kx_max = kx_grid(nk)

    n_electron = 0.0_dp
    n_hole = 0.0_dp

    ! --- Classify subbands at kx=0 (first point) ---
    allocate(sorted_evals(num_subbands))
    allocate(sort_idx(num_subbands))

    do s = 1, num_subbands
      sorted_evals(s) = eigenvalues_kx(s, 1)
      sort_idx(s) = s
    end do

    call sort_descending(sorted_evals, sort_idx, num_subbands)

    numcb_local = min(numcb, num_subbands)

    allocate(psi2(Ngrid, nk))
    allocate(integrand(nk))

    ! --- Electron density: CB subbands ---
    do j = 1, numcb_local
      s = sort_idx(j)

      call compute_psi2_wire(psi2, eigenvectors, s, Ngrid, nk, num_subbands)

      call accumulate_wire_band_density(n_electron, psi2, eigenvalues_kx, &
        & kx_grid, integrand, fermi_level, temperature, s, Ngrid, nk, &
        & num_subbands, kx_max, is_cb=.true.)
    end do

    ! --- Hole density: VB subbands ---
    do j = numcb_local + 1, num_subbands
      s = sort_idx(j)

      call compute_psi2_wire(psi2, eigenvectors, s, Ngrid, nk, num_subbands)

      call accumulate_wire_band_density(n_hole, psi2, eigenvalues_kx, &
        & kx_grid, integrand, fermi_level, temperature, s, Ngrid, nk, &
        & num_subbands, kx_max, is_cb=.false.)
    end do

    ! --- Convert 1/AA^3 to cm^-3 ---
    n_electron = n_electron * 1.0e24_dp
    n_hole = n_hole * 1.0e24_dp

    deallocate(psi2, integrand, sorted_evals, sort_idx)

  end subroutine compute_charge_density_wire


  ! ------------------------------------------------------------------
  ! Compute |Psi_s(y,z)|^2 summed over all 8 band components (wire)
  ! ------------------------------------------------------------------
  subroutine compute_psi2_wire(psi2, eigenvectors, s, Ngrid, nk, num_subbands)
    ! For a wire, eigenvectors are (8*Ngrid, num_subbands, nk).
    ! Component nu at flat grid point p is at index (nu-1)*Ngrid + p.

    integer, intent(in)          :: s, Ngrid, nk, num_subbands
    real(kind=dp), intent(out)   :: psi2(Ngrid, nk)
    complex(kind=dp), intent(in) :: eigenvectors(8*Ngrid, num_subbands, nk)

    integer :: idx, band_start, p

    psi2 = 0.0_dp
    do idx = 1, nk
      do band_start = 1, 8
        do p = 1, Ngrid
          psi2(p, idx) = psi2(p, idx) &
            & + abs(eigenvectors((band_start - 1)*Ngrid + p, s, idx))**2
        end do
      end do
    end do

  end subroutine compute_psi2_wire


  ! ------------------------------------------------------------------
  ! Accumulate wire band density via 1D k_x integration
  ! ------------------------------------------------------------------
  subroutine accumulate_wire_band_density(n_density, psi2, eigenvalues_kx, &
      & kx_grid, integrand, fermi_level, temperature, s, Ngrid, nk, &
      & num_subbands, kx_max, is_cb)
    ! Integrate |Psi_s|^2 * f * 1/(2*pi) over k_x using Simpson's rule.
    ! Weight is 1/(2*pi) (1D integral, no k_par factor).

    integer, intent(in)          :: Ngrid, nk, num_subbands
    real(kind=dp), intent(inout) :: n_density(Ngrid)
    real(kind=dp), intent(in)    :: psi2(Ngrid, nk)
    real(kind=dp), intent(in)    :: eigenvalues_kx(num_subbands, nk)
    real(kind=dp), intent(in)    :: kx_grid(nk)
    real(kind=dp), intent(inout) :: integrand(nk)
    real(kind=dp), intent(in)    :: fermi_level, temperature, kx_max
    integer, intent(in)          :: s
    logical, intent(in)          :: is_cb

    integer :: p, idx
    real(kind=dp) :: occ(nk)

    ! Pre-compute weighted occupation (independent of spatial point)
    do idx = 1, nk
      occ(idx) = fermi_dirac(eigenvalues_kx(s, idx), fermi_level, temperature)
      if (.not. is_cb) occ(idx) = 1.0_dp - occ(idx)
      occ(idx) = occ(idx) / (2.0_dp * pi_dp)
    end do

    do p = 1, Ngrid
      integrand(1:nk) = psi2(p, 1:nk) * occ(1:nk)
      n_density(p) = n_density(p) + 2.0_dp * simpson_real(integrand, 0.0_dp, kx_max)
    end do

  end subroutine accumulate_wire_band_density


end module charge_density
