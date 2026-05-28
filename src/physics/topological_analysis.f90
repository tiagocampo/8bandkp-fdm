module topological_analysis

  use definitions
  use sparse_matrices
  use linalg, only: zheev, zgetrf
  use hamiltonianConstructor, only: ZB8bandQW
  implicit none
  private

  public :: compute_chern_qwz
  public :: compute_berry_curvature_lattice
  public :: compute_hall_conductance
  public :: compute_conductance_kubo
  public :: compute_z2_gap
  public :: compute_z2_gap_edge
  public :: compute_z2_fukane
  public :: compute_z2_fukane_qw
  public :: compute_z2_fukane_qw_result
  public :: qw_inversion_expectation
  public :: qw_pair_inversion_sign
  public :: z2_from_trim_parities
  public :: compute_majorana_profile
  public :: build_bhz_wire_hamiltonian
  public :: bhz_wire_params
  public :: extract_edge_states_wire
  public :: fit_exponential_decay
  public :: compute_phase_diagram
  public :: compute_z2_gap_sweep
  public :: detect_z2_transitions
  public :: is_z2_transition
  public :: gap_closing_detect
  public :: bdg_zero_energy_gap

  integer, parameter :: topo_status_ok = 0
  integer, parameter :: topo_status_invalid = 1
  integer, parameter :: topo_status_asymmetric = 2
  integer, parameter :: topo_status_lapack = 3

  type :: bhz_wire_params
    real(kind=dp) :: A = 364.5_dp
    real(kind=dp) :: B = -686.0_dp
    real(kind=dp) :: D = -512.0_dp
    real(kind=dp) :: M = 10.0_dp
    real(kind=dp) :: d_wire = 58.0_dp
    integer       :: N = 200
    real(kind=dp) :: dz = 1.0_dp
  end type bhz_wire_params

contains

  function compute_chern_qwz(u, nk) result(C)
    implicit none
    real(kind=dp), intent(in) :: u
    integer, intent(in) :: nk
    integer :: C

    integer :: i, j, ip1, jp1
    complex(kind=dp), allocatable :: evecs(:,:,:)
    real(kind=dp) :: kx, ky, dk, total_flux, mz_val, E_plus, nrm
    complex(kind=dp) :: Ux, Uy, prod, off_diag, ev_tmp(2)

    dk = 2.0_dp * pi_dp / real(nk, kind=dp)
    total_flux = 0.0_dp

    allocate(evecs(nk, nk, 2))

    do j = 1, nk
      ky = -pi_dp + real(j-1, kind=dp) * dk
      do i = 1, nk
        kx = -pi_dp + real(i-1, kind=dp) * dk

        ! QWZ: H = sin(kx)*sigma_x + sin(ky)*sigma_y + mz*sigma_z
        ! H = [[mz, sin(kx)-i*sin(ky)], [sin(kx)+i*sin(ky), -mz]]
        mz_val = u + cos(kx) + cos(ky)
        off_diag = cmplx(sin(kx), -sin(ky), kind=dp)
        E_plus = sqrt(mz_val**2 + sin(kx)**2 + sin(ky)**2)

        ! Upper band eigenvector: [off_diag, E_plus - mz_val]
        ev_tmp(1) = off_diag
        ev_tmp(2) = cmplx(E_plus - mz_val, 0.0_dp, kind=dp)

        nrm = sqrt(real(conjg(ev_tmp(1))*ev_tmp(1) + conjg(ev_tmp(2))*ev_tmp(2), kind=dp))
        if (nrm > 1.0e-30_dp) then
          evecs(i,j,1) = ev_tmp(1) / nrm
          evecs(i,j,2) = ev_tmp(2) / nrm
        else
          evecs(i,j,1) = cmplx(1.0_dp, 0.0_dp, kind=dp)
          evecs(i,j,2) = cmplx(0.0_dp, 0.0_dp, kind=dp)
        end if
      end do
    end do

    ! FHS plaquette: (i,j)->(i+1,j)->(i+1,j+1)->(i,j+1)->(i,j)
    do j = 1, nk
      jp1 = mod(j, nk) + 1
      do i = 1, nk
        ip1 = mod(i, nk) + 1

        ! Ux(i,j) = <u(i,j)|u(i+1,j)>
        Ux = conjg(evecs(i,j,1))*evecs(ip1,j,1) &
           + conjg(evecs(i,j,2))*evecs(ip1,j,2)
        ! Uy(i+1,j) = <u(i+1,j)|u(i+1,j+1)>
        Uy = conjg(evecs(ip1,j,1))*evecs(ip1,jp1,1) &
           + conjg(evecs(ip1,j,2))*evecs(ip1,jp1,2)
        ! conjg(Ux(i,j+1)) = <u(i+1,j+1)|u(i,j+1)>
        prod = conjg(evecs(i,jp1,1))*evecs(ip1,jp1,1) &
             + conjg(evecs(i,jp1,2))*evecs(ip1,jp1,2)
        ! conjg(Uy(i,j)) = <u(i,j+1)|u(i,j)>
        prod = Ux * Uy * conjg(prod) &
             * conjg(conjg(evecs(i,j,1))*evecs(i,jp1,1) &
             + conjg(evecs(i,j,2))*evecs(i,jp1,2))

        total_flux = total_flux + aimag(log(prod))
      end do
    end do

    C = nint(total_flux / (2.0_dp * pi_dp))
    deallocate(evecs)
  end function compute_chern_qwz

  function compute_berry_curvature_lattice(evecs_k, kx_arr, ky_arr, n_occ) result(Omega)
    implicit none
    complex(kind=dp), contiguous, intent(in) :: evecs_k(:,:,:,:)  ! (basis, n_occ, nkx, nky)
    real(kind=dp), contiguous, intent(in) :: kx_arr(:), ky_arr(:)
    integer, intent(in) :: n_occ
    real(kind=dp), allocatable :: Omega(:,:)

    integer :: nkx, nky, i, j, ip1, jp1
    complex(kind=dp) :: M_xx, M_yx, M_xy, M_yy
    complex(kind=dp), allocatable :: overlap(:,:)
    real(kind=dp) :: dkx, dky, dA
    integer :: m, p

    nkx = size(kx_arr)
    nky = size(ky_arr)
    if (nkx < 2 .or. nky < 2 .or. n_occ < 1 .or. size(evecs_k, 2) < n_occ) then
      allocate(Omega(0, 0))
      return
    end if

    dkx = kx_arr(min(2, nkx)) - kx_arr(1)
    dky = ky_arr(min(2, nky)) - ky_arr(1)
    if (abs(dkx) <= 1.0e-14_dp .or. abs(dky) <= 1.0e-14_dp) then
      allocate(Omega(0, 0))
      return
    end if

    allocate(Omega(nkx, nky))
    Omega = 0.0_dp

    dA = dkx * dky

    allocate(overlap(n_occ, n_occ))
    do j = 1, nky
      jp1 = mod(j, nky) + 1
      do i = 1, nkx
        ip1 = mod(i, nkx) + 1

        ! U_x = det(<u_n(k)|u_m(k+x)>)
        do m = 1, n_occ
          do p = 1, n_occ
            overlap(m, p) = sum(conjg(evecs_k(:, m, i, j)) * evecs_k(:, p, ip1, j))
          end do
        end do
        M_xx = det_small(overlap, n_occ)

        ! U_y(k+x) = det(<u_n(k+x)|u_m(k+x+y)>)
        do m = 1, n_occ
          do p = 1, n_occ
            overlap(m, p) = sum(conjg(evecs_k(:, m, ip1, j)) * evecs_k(:, p, ip1, jp1))
          end do
        end do
        M_yx = det_small(overlap, n_occ)

        ! conjg(U_x(k+y)): overlap computes <u_m(k+x+y)|u_p(k+y)>
        ! which is the adjoint of <u(k+y)|u(k+x+y)>, so its det = conjg(Ux(k,y+1))
        do m = 1, n_occ
          do p = 1, n_occ
            overlap(m, p) = sum(conjg(evecs_k(:, m, ip1, jp1)) * evecs_k(:, p, i, jp1))
          end do
        end do
        M_xy = det_small(overlap, n_occ)

        ! conjg(U_y(k)): overlap computes <u_m(k+y)|u_p(k)>
        ! which is the adjoint of <u(k)|u(k+y)>, so its det = conjg(Uy(k))
        do m = 1, n_occ
          do p = 1, n_occ
            overlap(m, p) = sum(conjg(evecs_k(:, m, i, jp1)) * evecs_k(:, p, i, j))
          end do
        end do
        M_yy = det_small(overlap, n_occ)

        Omega(i, j) = aimag(log(M_xx * M_yx * M_xy * M_yy)) / dA

      end do
    end do
    deallocate(overlap)

  end function compute_berry_curvature_lattice

  function det_small(A, n) result(d)
    implicit none
    complex(kind=dp), contiguous, intent(in) :: A(:,:)
    integer, intent(in) :: n
    complex(kind=dp) :: d

    complex(kind=dp), allocatable :: LU(:,:)
    integer, allocatable :: pivot(:)
    integer :: info, i

    select case(n)
    case(1)
      d = A(1,1)
    case(2)
      d = A(1,1)*A(2,2) - A(1,2)*A(2,1)
    case default
      allocate(LU(n, n), pivot(n))
      LU = A
      call zgetrf(n, n, LU, n, pivot, info)
      d = cmplx(1.0_dp, 0.0_dp, kind=dp)
      do i = 1, n
        d = d * LU(i, i)
        if (pivot(i) /= i) d = -d
      end do
      deallocate(LU, pivot)
    end select
  end function det_small

  elemental function compute_hall_conductance(C) result(sigma_xy)
    implicit none
    integer, intent(in) :: C
    real(kind=dp) :: sigma_xy

    sigma_xy = real(C, kind=dp)
  end function compute_hall_conductance

  function compute_conductance_kubo(berry_curvature, kx_arr, ky_arr) result(sigma_xy)
    implicit none
    real(kind=dp), contiguous, intent(in) :: berry_curvature(:,:)
    real(kind=dp), contiguous, intent(in) :: kx_arr(:), ky_arr(:)
    real(kind=dp) :: sigma_xy

    real(kind=dp) :: dkx, dky

    sigma_xy = 0.0_dp
    if (size(kx_arr) < 2 .or. size(ky_arr) < 2) return
    if (size(berry_curvature, 1) /= size(kx_arr)) return
    if (size(berry_curvature, 2) /= size(ky_arr)) return

    dkx = kx_arr(2) - kx_arr(1)
    dky = ky_arr(2) - ky_arr(1)

    ! Periodic plaquette curvature: direct uniform integration in units of e^2/h.
    sigma_xy = sum(berry_curvature) * dkx * dky / (2.0_dp * pi_dp)
  end function compute_conductance_kubo

  function compute_z2_gap(eigenvalues, gap_threshold) result(z2)
    implicit none
    real(kind=dp), contiguous, intent(in) :: eigenvalues(:)
    real(kind=dp), intent(in) :: gap_threshold

    integer :: z2
    real(kind=dp) :: gap_min, gap_max
    integer :: n_in_gap, i

    z2 = 0

    if (size(eigenvalues) < 2) return

    gap_min = -gap_threshold
    gap_max = gap_threshold

    n_in_gap = 0
    do i = 1, size(eigenvalues)
      if (eigenvalues(i) > gap_min .and. eigenvalues(i) < gap_max) then
        n_in_gap = n_in_gap + 1
      end if
    end do

    ! Z2 in 1D counts edge state pairs: topological phase has exactly 2 edge
    ! states (one per end) within the bulk gap. Divide by 2 to get pair count.
    if (n_in_gap >= 2) z2 = 1

  end function compute_z2_gap

  ! ==============================================================================
  ! Z2 from gap + spatial localization (R9).
  !
  ! Like compute_z2_gap, but additionally checks that near-zero eigenvalues have
  ! eigenvectors spatially localized at the wire edges.  This distinguishes
  ! topological edge states from numerical noise.
  !
  ! The wire is assumed to have a 4-band basis per site (BHZ model).
  ! N_sites = size(eigenvectors, 1) / 4.
  !
  ! Edge localization check: compute the fraction of |psi|^2 in the first and
  ! last 10% of sites.  If this fraction exceeds edge_fraction_threshold
  ! (default 0.5), the state is edge-localized.
  ! ==============================================================================
  function compute_z2_gap_edge(eigenvalues, eigenvectors, gap_threshold, &
                                edge_fraction_threshold) result(z2)
    implicit none
    real(kind=dp), contiguous, intent(in) :: eigenvalues(:)
    complex(kind=dp), contiguous, intent(in) :: eigenvectors(:,:)
    real(kind=dp), intent(in) :: gap_threshold
    real(kind=dp), intent(in) :: edge_fraction_threshold

    integer :: z2
    real(kind=dp) :: gap_min, gap_max
    integer :: n_edge_states, i, j, k, N_sites, n_bands
    real(kind=dp) :: total_weight, edge_weight
    integer :: n_edge_sites

    z2 = 0

    if (size(eigenvalues) < 2) return
    if (size(eigenvectors, 2) < size(eigenvalues)) return

    ! Determine basis size per site
    n_bands = 4  ! BHZ 4-band model
    N_sites = size(eigenvectors, 1) / n_bands
    if (N_sites < 2) return

    gap_min = -gap_threshold
    gap_max = gap_threshold

    ! Number of sites in each edge region (~10% of wire)
    n_edge_sites = max(1, N_sites / 10)

    n_edge_states = 0
    do i = 1, size(eigenvalues)
      if (eigenvalues(i) <= gap_min .or. eigenvalues(i) >= gap_max) cycle

      ! Check spatial localization at wire edges
      total_weight = 0.0_dp
      edge_weight = 0.0_dp
      do j = 1, N_sites
        do k = 1, n_bands
          total_weight = total_weight + abs(eigenvectors((j-1)*n_bands + k, i))**2
        end do
      end do

      if (total_weight < 1.0e-30_dp) cycle

      ! Sum weight in first and last n_edge_sites
      do j = 1, n_edge_sites
        do k = 1, n_bands
          edge_weight = edge_weight + abs(eigenvectors((j-1)*n_bands + k, i))**2
        end do
      end do
      do j = N_sites - n_edge_sites + 1, N_sites
        do k = 1, n_bands
          edge_weight = edge_weight + abs(eigenvectors((j-1)*n_bands + k, i))**2
        end do
      end do

      ! State is edge-localized if edge weight fraction exceeds threshold
      if (edge_weight / total_weight > edge_fraction_threshold) then
        n_edge_states = n_edge_states + 1
      end if
    end do

    ! Z2 in 1D: topological phase has at least 2 edge states (one per end)
    if (n_edge_states >= 2) z2 = 1

  end function compute_z2_gap_edge

  function compute_z2_fukane(cfg, profile, kpterms, n_occ) result(z2)
    implicit none
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: profile(:,:), kpterms(:,:,:)
    integer, intent(in) :: n_occ
    integer :: z2

    ! Delegate to the full Fu-Kane QW implementation
    z2 = compute_z2_fukane_qw(cfg, profile, kpterms, n_occ)

  end function compute_z2_fukane

  function extract_edge_states_wire(eigenvalues, eigenvectors, grid, window) result(edge_info)
    implicit none
    real(kind=dp), contiguous, intent(in) :: eigenvalues(:)
    complex(kind=dp), contiguous, intent(in) :: eigenvectors(:,:)
    type(spatial_grid), intent(in) :: grid
    real(kind=dp), intent(in) :: window
    real(kind=dp), allocatable :: edge_info(:)

    integer :: N, Nstates, i, j, k
    real(kind=dp) :: E_min, E_max
    real(kind=dp), allocatable :: density(:,:)
    real(kind=dp), allocatable :: edge_xi_list(:)
    real(kind=dp) :: xi_min, xi_avg, xi_val
    integer :: n_edge
    logical :: success

    N = size(eigenvectors, 1) / 4
    if (N < 1 .or. size(grid%z) < N) then
      print *, 'Warning: extract_edge_states_wire: grid/eigenvector mismatch'
      print *, '  N_from_evec=', N, ' grid%npoints=', grid%npoints()
      allocate(edge_info(3)); edge_info = 0.0_dp; return
    end if
    Nstates = size(eigenvalues)

    allocate(density(N, Nstates))
    allocate(edge_xi_list(Nstates))

    do j = 1, Nstates
      do i = 1, N
        density(i, j) = 0.0_dp
        do k = 1, 4
          density(i, j) = density(i, j) + abs(eigenvectors((i-1)*4 + k, j))**2
        end do
      end do
    end do

    E_min = -window
    E_max = window

    n_edge = 0
    xi_min = 0.0_dp
    xi_avg = 0.0_dp
    edge_xi_list = 0.0_dp

    do j = 1, Nstates
      if (eigenvalues(j) >= E_min .and. eigenvalues(j) <= E_max) then
        call fit_exponential_decay(density(:, j), grid%z, xi_val, success)

        if (success .and. xi_val > 0.0_dp) then
          n_edge = n_edge + 1
          edge_xi_list(n_edge) = xi_val
          if (xi_min == 0.0_dp .or. xi_val < xi_min) then
            xi_min = xi_val
          end if
          xi_avg = xi_avg + xi_val
        end if
      end if
    end do

    if (n_edge > 0) then
      xi_avg = xi_avg / real(n_edge, kind=dp)
    else
      xi_avg = 0.0_dp
      xi_min = 0.0_dp
    end if

    allocate(edge_info(3))
    edge_info(1) = xi_min
    edge_info(2) = xi_avg
    edge_info(3) = real(n_edge, kind=dp)

    deallocate(density, edge_xi_list)

  end function extract_edge_states_wire

  ! ==============================================================================
  ! Private helper: tail-search + log-linear regression for exponential decay.
  !
  ! Finds the tail start (where density drops below threshold * peak), then
  ! performs linear regression on (position, log(density)) for nonzero points.
  !
  ! R10 fix: n_fit_actual tracks only points that contribute to the regression
  ! sums (nonzero density after log transform), not all tail positions.
  !
  ! R11: Returns converged=.false. when xi is comparable to or larger than
  ! the domain extent, with the raw xi value preserved for caller inspection.
  ! ==============================================================================
  subroutine fit_tail_exponential(positions, density, n_points, threshold, &
                                   xi, n_fit_actual, converged)
    implicit none
    real(kind=dp), contiguous, intent(in) :: positions(:)
    real(kind=dp), contiguous, intent(in) :: density(:)
    integer, intent(in) :: n_points
    real(kind=dp), intent(in) :: threshold
    real(kind=dp), intent(out) :: xi
    integer, intent(out) :: n_fit_actual
    logical, intent(out) :: converged

    integer :: peak_idx, tail_start, i
    logical :: forward_tail
    real(kind=dp) :: rho_peak, x_start, domain_extent
    real(kind=dp) :: sum_y_log, sum_x, sum_x2, sum_xy, denom, slope

    xi = 0.0_dp
    n_fit_actual = 0
    converged = .false.

    if (n_points < 3) return
    if (size(positions) < n_points .or. size(density) < n_points) return

    peak_idx = maxloc(density(1:n_points), dim=1)
    rho_peak = density(peak_idx)

    if (rho_peak <= 0.0_dp) return

    ! Search forward from peak for tail start
    tail_start = 0
    forward_tail = .true.
    do i = peak_idx + 1, n_points
      if (density(i) < threshold * rho_peak) then
        tail_start = i
        exit
      end if
    end do

    ! Search backward if no forward tail found
    if (tail_start == 0) then
      forward_tail = .false.
      do i = peak_idx - 1, 1, -1
        if (density(i) < threshold * rho_peak) then
          tail_start = i
          exit
        end if
      end do
    end if

    if (tail_start == 0) return
    if (forward_tail .and. tail_start >= n_points) return
    if (.not. forward_tail .and. tail_start <= 1) return

    x_start = positions(tail_start)

    ! Accumulate regression sums — only count nonzero-density points (R10)
    sum_y_log = 0.0_dp
    sum_x = 0.0_dp
    sum_x2 = 0.0_dp
    sum_xy = 0.0_dp
    n_fit_actual = 0

    if (forward_tail) then
      ! Tail is to the right of peak: iterate forward
      do i = tail_start, n_points
        if (density(i) > 1.0e-14_dp) then
          n_fit_actual = n_fit_actual + 1
          sum_y_log = sum_y_log + log(density(i))
          sum_x = sum_x + (positions(i) - x_start)
          sum_x2 = sum_x2 + (positions(i) - x_start)**2
          sum_xy = sum_xy + (positions(i) - x_start) * log(density(i))
        end if
      end do
    else
      ! Tail is to the left of peak: iterate backward toward index 1
      do i = tail_start, 1, -1
        if (density(i) > 1.0e-14_dp) then
          n_fit_actual = n_fit_actual + 1
          sum_y_log = sum_y_log + log(density(i))
          sum_x = sum_x + (x_start - positions(i))
          sum_x2 = sum_x2 + (x_start - positions(i))**2
          sum_xy = sum_xy + (x_start - positions(i)) * log(density(i))
        end if
      end do
    end if

    ! Fallback: if forward tail found no data, try backward direction
    if (n_fit_actual < 3 .and. forward_tail) then
      forward_tail = .false.
      do i = peak_idx - 1, 1, -1
        if (density(i) < threshold * rho_peak) then
          tail_start = i
          exit
        end if
      end do
      if (tail_start == 0 .or. tail_start <= 1) return
      x_start = positions(tail_start)
      sum_y_log = 0.0_dp
      sum_x = 0.0_dp
      sum_x2 = 0.0_dp
      sum_xy = 0.0_dp
      n_fit_actual = 0
      do i = tail_start, 1, -1
        if (density(i) > 1.0e-14_dp) then
          n_fit_actual = n_fit_actual + 1
          sum_y_log = sum_y_log + log(density(i))
          sum_x = sum_x + (x_start - positions(i))
          sum_x2 = sum_x2 + (x_start - positions(i))**2
          sum_xy = sum_xy + (x_start - positions(i)) * log(density(i))
        end if
      end do
    end if

    if (n_fit_actual < 3) return

    denom = sum_x2 - sum_x**2 / real(n_fit_actual, kind=dp)
    if (abs(denom) < 1.0e-14_dp) return

    slope = (sum_xy - sum_x * sum_y_log / real(n_fit_actual, kind=dp)) / denom
    if (abs(slope) < tiny(1.0_dp)) return

    xi = abs(-1.0_dp / slope)

    ! R11: Check near-transition regime — xi comparable to domain extent
    if (forward_tail) then
      domain_extent = abs(positions(n_points) - positions(tail_start))
    else
      domain_extent = abs(positions(tail_start) - positions(1))
    end if
    if (domain_extent > 0.0_dp .and. xi < domain_extent) then
      converged = .true.
    else if (domain_extent > 0.0_dp) then
      ! xi >= domain_extent: near-transition, return raw value unconverged
      converged = .false.
    else
      converged = .true.
    end if

  end subroutine fit_tail_exponential

  subroutine fit_exponential_decay(density, x, xi, success)
    implicit none
    real(kind=dp), contiguous, intent(in) :: density(:)
    real(kind=dp), contiguous, intent(in) :: x(:)
    real(kind=dp), intent(out) :: xi
    logical, intent(out) :: success

    integer :: n, n_fit_actual
    real(kind=dp), parameter :: tail_threshold = 0.1_dp
    logical :: converged

    n = size(density)
    success = .false.
    xi = 0.0_dp

    if (n < 3 .or. size(x) /= n) return

    call fit_tail_exponential(x, density, n, tail_threshold, &
                               xi, n_fit_actual, converged)

    success = (xi > 0.0_dp .and. n_fit_actual >= 3)

  end subroutine fit_exponential_decay

  elemental integer function band_major_row(iband, isite, nsite) result(irow)
    implicit none
    integer, intent(in) :: iband, isite, nsite
    irow = (iband - 1) * nsite + isite
  end function band_major_row

  elemental real(kind=dp) function band_inversion_parity(iband) result(p)
    implicit none
    integer, intent(in) :: iband

    if (iband <= 6) then
      p = 1.0_dp
    else
      p = -1.0_dp
    end if
  end function band_inversion_parity

  function qw_inversion_expectation(psi, nsite) result(parity)
    implicit none
    complex(kind=dp), contiguous, intent(in) :: psi(:)
    integer, intent(in) :: nsite
    real(kind=dp) :: parity
    integer :: ib, isite, imirror, irow, jrow

    parity = 0.0_dp
    if (nsite < 1 .or. size(psi) < 8 * nsite) return

    do ib = 1, 8
      do isite = 1, nsite
        imirror = nsite + 1 - isite
        irow = band_major_row(ib, isite, nsite)
        jrow = band_major_row(ib, imirror, nsite)
        parity = parity + band_inversion_parity(ib) * real(conjg(psi(irow)) * psi(jrow), kind=dp)
      end do
    end do
  end function qw_inversion_expectation

  function qw_inversion_matrix_element(psi_left, psi_right, nsite) result(element)
    implicit none
    complex(kind=dp), contiguous, intent(in) :: psi_left(:), psi_right(:)
    integer, intent(in) :: nsite
    complex(kind=dp) :: element
    integer :: ib, isite, imirror, irow, jrow

    element = cmplx(0.0_dp, 0.0_dp, kind=dp)
    if (nsite < 1) return
    if (size(psi_left) < 8 * nsite .or. size(psi_right) < 8 * nsite) return

    do ib = 1, 8
      do isite = 1, nsite
        imirror = nsite + 1 - isite
        irow = band_major_row(ib, isite, nsite)
        jrow = band_major_row(ib, imirror, nsite)
        element = element + band_inversion_parity(ib) * conjg(psi_left(irow)) * psi_right(jrow)
      end do
    end do
  end function qw_inversion_matrix_element

  function qw_pair_inversion_sign(psi_a, psi_b, nsite) result(pair_sign)
    implicit none
    complex(kind=dp), contiguous, intent(in) :: psi_a(:), psi_b(:)
    integer, intent(in) :: nsite
    real(kind=dp) :: pair_sign
    complex(kind=dp) :: p11, p22

    p11 = qw_inversion_matrix_element(psi_a, psi_a, nsite)
    p22 = qw_inversion_matrix_element(psi_b, psi_b, nsite)
    pair_sign = 0.5_dp * real(p11 + p22, kind=dp)
  end function qw_pair_inversion_sign

  integer function z2_from_trim_parities(delta_trim) result(z2)
    implicit none
    real(kind=dp), intent(in) :: delta_trim(4)
    real(kind=dp) :: trim_product

    trim_product = delta_trim(1) * delta_trim(2) * delta_trim(3) * delta_trim(4)
    if (trim_product < 0.0_dp) then
      z2 = 1
    else
      z2 = 0
    end if
  end function z2_from_trim_parities

  function compute_majorana_profile(evec_bdg, grid, energy_tol, half_n_in, profile) result(xi)
    implicit none
    complex(kind=dp), contiguous, intent(in) :: evec_bdg(:)
    type(spatial_grid), intent(in) :: grid
    real(kind=dp), intent(in) :: energy_tol
    integer, intent(in) :: half_n_in
    real(kind=dp), intent(out), optional :: profile(:)
    real(kind=dp) :: xi

    integer :: ngrid, half_n, nspatial
    real(kind=dp), allocatable :: rho(:)
    real(kind=dp), allocatable :: xx(:)
    integer :: i, ib, n_fit_actual
    real(kind=dp) :: norm_factor
    real(kind=dp), parameter :: tail_threshold = 0.1_dp
    logical :: converged

    ngrid = grid%npoints()
    half_n = half_n_in
    if (half_n * 2 > size(evec_bdg)) then
      print *, 'ERROR: compute_majorana_profile: eigenvector too small'
      print *, '  half_n=', half_n, ' size(evec)=', size(evec_bdg)
      xi = -1.0_dp
      return
    end if
    nspatial = half_n / 8

    allocate(rho(nspatial))
    do i = 1, nspatial
      rho(i) = 0.0_dp
      do ib = 1, 8
        rho(i) = rho(i) + abs(evec_bdg(band_major_row(ib, i, nspatial)))**2
        rho(i) = rho(i) + abs(evec_bdg(half_n + band_major_row(ib, i, nspatial)))**2
      end do
    end do

    norm_factor = sum(rho)
    if (norm_factor > 0.0_dp) then
      rho = rho / norm_factor
    end if

    if (present(profile)) then
      profile = 0.0_dp
      if (size(profile) >= nspatial) then
        profile(1:nspatial) = rho
      else
        profile = rho(1:size(profile))
      end if
    end if

    ! Build position array based on grid type
    if (grid%ndim == 2 .and. allocated(grid%coords)) then
      allocate(xx(nspatial))
      xx(1:nspatial) = grid%coords(2, 1:nspatial)
    else if (allocated(grid%z)) then
      allocate(xx(nspatial))
      xx = grid%z
    else
      allocate(xx(nspatial))
      do i = 1, nspatial
        xx(i) = real(i - 1, kind=dp)
      end do
    end if

    ! Delegate tail search + regression to shared helper (R10, R11, R17)
    call fit_tail_exponential(xx, rho, nspatial, tail_threshold, &
                               xi, n_fit_actual, converged)

    ! R11: near-transition regime returns raw xi with converged=.false.
    ! Return -1 only for complete failure (no fit at all)
    if (xi <= 0.0_dp .and. .not. converged) then
      xi = -1.0_dp
    end if

    deallocate(rho, xx)

  end function compute_majorana_profile

  subroutine build_bhz_wire_hamiltonian(H_csr, params)
    implicit none
    type(csr_matrix), intent(out) :: H_csr
    type(bhz_wire_params), intent(in) :: params

    integer :: N, i, row, nnz_total
    real(kind=dp) :: dz, B_plus_D_over_dz2, A_over_2dz
    complex(kind=dp), allocatable :: coo_vals(:)
    integer, allocatable :: coo_row(:), coo_col(:)
    integer :: nnz_offset

    if (params%N < 1 .or. params%dz <= 0.0_dp) then
      print *, 'ERROR: build_bhz_wire_hamiltonian: invalid parameters'
      print *, '  N=', params%N, ' dz=', params%dz
      stop 1
    end if

    N = params%N
    dz = params%dz

    B_plus_D_over_dz2 = (params%B + params%D) * 0.0001_dp / (dz * dz)
    A_over_2dz = params%A * 0.02_dp / (2.0_dp * dz)

    ! Count entries: diagonal = 4N
    ! A-fwd (i<N): 4*(N-1)
    ! A-bwd (i>1): 4*(N-1)
    ! B+D-fwd (i<N): 4*(N-1) - self-diagonal entries that will be merged
    ! B+D-bwd (i>1): 4*(N-1)
    ! Note: (row,row) entries are merged in csr_build_from_coo, so we count them
    nnz_total = 4*N + 4*(N-1) + 4*(N-1) + 4*(N-1) + 4*(N-1)
    allocate(coo_vals(nnz_total), coo_row(nnz_total), coo_col(nnz_total))

    nnz_offset = 0

    do i = 1, N
      do row = 1, 4
        nnz_offset = nnz_offset + 1
        coo_row(nnz_offset) = (i-1)*4 + row
        coo_col(nnz_offset) = (i-1)*4 + row
        if (row <= 2) then
          coo_vals(nnz_offset) = cmplx(params%M - 2.0_dp * B_plus_D_over_dz2, 0.0_dp, kind=dp)
        else
          coo_vals(nnz_offset) = cmplx(-params%M - 2.0_dp * B_plus_D_over_dz2, 0.0_dp, kind=dp)
        end if

        if (i < N) then
          nnz_offset = nnz_offset + 1
          coo_row(nnz_offset) = (i-1)*4 + row
          coo_col(nnz_offset) = i*4 + 5 - row
          if (row == 1 .or. row == 4) then
            coo_vals(nnz_offset) = cmplx(A_over_2dz, 0.0_dp, kind=dp)
          else
            coo_vals(nnz_offset) = cmplx(-A_over_2dz, 0.0_dp, kind=dp)
          end if
        end if

        if (i > 1) then
          nnz_offset = nnz_offset + 1
          coo_row(nnz_offset) = (i-1)*4 + row
          coo_col(nnz_offset) = (i-2)*4 + 5 - row
          if (row == 1 .or. row == 4) then
            coo_vals(nnz_offset) = cmplx(A_over_2dz, 0.0_dp, kind=dp)
          else
            coo_vals(nnz_offset) = cmplx(-A_over_2dz, 0.0_dp, kind=dp)
          end if
        end if

        if (i < N) then
          nnz_offset = nnz_offset + 1
          coo_row(nnz_offset) = (i-1)*4 + row
          coo_col(nnz_offset) = i*4 + mod(row, 4) + 1
          coo_vals(nnz_offset) = cmplx(B_plus_D_over_dz2, 0.0_dp, kind=dp)
        end if

        if (i > 1) then
          nnz_offset = nnz_offset + 1
          coo_row(nnz_offset) = (i-1)*4 + row
          coo_col(nnz_offset) = (i-2)*4 + mod(row + 2, 4) + 1
          coo_vals(nnz_offset) = cmplx(B_plus_D_over_dz2, 0.0_dp, kind=dp)
        end if
      end do
    end do

    call csr_build_from_coo(H_csr, 4*N, 4*N, nnz_offset, &
                            coo_row(1:nnz_offset), coo_col(1:nnz_offset), &
                            coo_vals(1:nnz_offset))

    deallocate(coo_vals, coo_row, coo_col)

  end subroutine build_bhz_wire_hamiltonian

  ! ==============================================================================
  ! BdG zero-energy gap: minimum distance of eigenvalues from zero energy.
  ! For a superconductor with particle-hole symmetry, eigenvalues come in +/- E
  ! pairs.  The superconducting gap is min|E|, NOT the minimum adjacent spacing.
  ! ==============================================================================
  pure function bdg_zero_energy_gap(eigenvalues) result(gap)
    real(kind=dp), intent(in), contiguous :: eigenvalues(:)
    real(kind=dp) :: gap

    if (size(eigenvalues) < 1) then
      gap = 0.0_dp
    else
      gap = minval(abs(eigenvalues))
    end if
  end function bdg_zero_energy_gap

  ! ==============================================================================
  ! Gap closing detection at a single (B, mu) point.
  ! Returns the minimum eigenvalue spacing (gap) around the Fermi level.
  ! ==============================================================================
  function gap_closing_detect(eigenvalues, mu, gap_threshold) result(is_gap_closing)
    implicit none
    real(kind=dp), contiguous, intent(in) :: eigenvalues(:)
    real(kind=dp), intent(in) :: mu
    real(kind=dp), intent(in) :: gap_threshold
    logical :: is_gap_closing

    real(kind=dp) :: gap_min
    integer :: i

    if (size(eigenvalues) < 2) then
      is_gap_closing = .false.
      return
    end if

    gap_min = huge(1.0_dp)
    do i = 1, size(eigenvalues) - 1
      if (abs(eigenvalues(i) - mu) < gap_threshold .or. &
          abs(eigenvalues(i+1) - mu) < gap_threshold) then
        gap_min = min(gap_min, abs(eigenvalues(i+1) - eigenvalues(i)))
      end if
    end do

    is_gap_closing = (gap_min < gap_threshold)

  end function gap_closing_detect

  ! ==============================================================================
  ! Topological phase diagram computation.
  !
  ! Sweeps B-field and chemical potential mu, tracks gap closings.
  ! Each gap closing flips the Z2 invariant.
  !
  ! Output arrays (allocated within the routine):
  !   gap_array(nB, nMu)     : minimum eigenvalue gap at each (B, mu)
  !   z2_array(nB, nMu)      : Z2 invariant (0 or 1) at each (B, mu)
  !   gap_line(n_transitions): B values where gap closes (Z2 transitions)
  ! ==============================================================================
  subroutine compute_phase_diagram(B_min, B_max, nB, mu_min, mu_max, nMu, &
                                    eigenvalues_func, gap_threshold, &
                                    gap_array, z2_array, gap_line, n_transitions)
    implicit none
    real(kind=dp), intent(in) :: B_min, B_max
    integer, intent(in) :: nB
    real(kind=dp), intent(in) :: mu_min, mu_max
    integer, intent(in) :: nMu
    interface
      function eigenvalues_func(B, mu) result(evals)
        import :: dp
        real(kind=dp), intent(in) :: B, mu
        real(kind=dp), allocatable :: evals(:)
      end function eigenvalues_func
    end interface
    real(kind=dp), intent(in) :: gap_threshold
    real(kind=dp), allocatable, intent(out) :: gap_array(:,:)
    integer, allocatable, intent(out) :: z2_array(:,:)
    real(kind=dp), allocatable, intent(out) :: gap_line(:)
    integer, intent(out) :: n_transitions

    real(kind=dp) :: B_val, mu_val, dB, dmu
    real(kind=dp), allocatable :: eigenvalues(:)
    real(kind=dp) :: gap_min
    integer :: iB, iMu, i, j
    integer :: z2_cumulative
    logical :: is_gap_closing
    real(kind=dp), allocatable :: gap_line_raw(:)
    integer :: n_raw

    if (nB < 1 .or. nMu < 1 .or. nB > 1000 .or. nMu > 1000) then
      print *, 'ERROR: compute_phase_diagram: nB, nMu must be in [1, 1000]'
      print *, '  nB=', nB, ' nMu=', nMu
      stop 1
    end if

    dB = (B_max - B_min) / real(max(1, nB - 1), kind=dp)
    dmu = (mu_max - mu_min) / real(max(1, nMu - 1), kind=dp)

    allocate(gap_array(nB, nMu))
    allocate(z2_array(nB, nMu))
    gap_array = 0.0_dp
    z2_array = 0

    ! Sweep over mu (outer) then B (inner) for cumulative Z2 tracking.
    ! Z2 starts trivial (0) at B=B_min and flips each time the gap closes.
    ! This gives physically correct phase diagram: Z2 accumulates transitions
    ! along the B-sweep rather than marking individual gap-closing points.
    do iMu = 1, nMu
      mu_val = mu_min + real(iMu - 1, kind=dp) * dmu

      ! Cumulative Z2 starts at 0 (trivial at B=B_min)
      z2_cumulative = 0

      do iB = 1, nB
        B_val = B_min + real(iB - 1, kind=dp) * dB

        eigenvalues = eigenvalues_func(B_val, mu_val)

        ! Compute minimum gap around Fermi level
        gap_min = huge(1.0_dp)
        do i = 1, size(eigenvalues) - 1
          if (abs(eigenvalues(i) - mu_val) < gap_threshold .or. &
              abs(eigenvalues(i+1) - mu_val) < gap_threshold) then
            gap_min = min(gap_min, abs(eigenvalues(i+1) - eigenvalues(i)))
          end if
        end do

        gap_array(iB, iMu) = gap_min
        is_gap_closing = (gap_min < gap_threshold)

        ! Cumulative Z2 tracking: flip when gap closes
        if (is_gap_closing) then
          z2_cumulative = 1 - z2_cumulative
        end if
        z2_array(iB, iMu) = z2_cumulative

        deallocate(eigenvalues)
      end do
    end do

    ! Track Z2 transitions along B sweep at fixed mu (use central mu)
    n_raw = 0
    allocate(gap_line_raw(nB))
    j = (nMu + 1) / 2  ! central mu index

    do iB = 1, nB - 1
      if (z2_array(iB, j) /= z2_array(iB + 1, j)) then
        n_raw = n_raw + 1
        gap_line_raw(n_raw) = B_min + real(iB - 1, kind=dp) * dB + dB / 2.0_dp
      end if
    end do

    n_transitions = n_raw
    allocate(gap_line(n_transitions))
    gap_line(1:n_transitions) = gap_line_raw(1:n_transitions)

    deallocate(gap_line_raw)

  end subroutine compute_phase_diagram

  ! ==============================================================================
  ! Fu-Kane Z2 invariant for a QW via parity eigenvalues at TRIM points.
  !
  ! Computes the Z2 topological invariant using the Fu-Kane formula:
  !   (-1)^Z2 = prod_{i=1}^{4} delta_i   (excluding Gamma for the strong index)
  ! where delta_i = prod_{m=1}^{n_occ} xi_m(Lambda_i) and xi_m is the parity
  ! eigenvalue of the m-th occupied state at TRIM point Lambda_i.
  !
  ! The 4 TRIM points in 2D are: (0,0), (pi/a,0), (0,pi/a), (pi/a,pi/a).
  !
  ! For the 8-band zinc-blende basis, inversion parity is:
  !   P_band = [+1,+1,+1,+1,+1,+1,-1,-1]
  ! (valence bands 1-6 are even under inversion, conduction bands 7-8 are odd).
  !
  ! For a QW eigenstate psi_n, the parity eigenvalue is:
  !   xi_n = sign( prod_{i,b} P_band(b)^{|psi_n(i,b)|^2} )
  ! which reduces to:
  !   xi_n = sign( sum over CB weights vs VB weights )
  !
  ! This function builds the QW Hamiltonian at each TRIM point via ZB8bandQW,
  ! diagonalizes it, and computes the parity product.
  !
  ! INPUT:
  !   cfg      - simulation configuration (grid, materials, FD order)
  !   profile  - band offset profile (N x 3)
  !   kpterms  - k.p material terms (N x N x 10)
  !   n_occ    - number of occupied bands (typically 4*N for half-filling)
  !
  ! OUTPUT:
  !   Z2 = 0 (trivial) or 1 (topological)
  ! ==============================================================================
  subroutine compute_z2_fukane_qw_result(cfg, profile, kpterms, n_occ, z2, min_gap, status)
    implicit none
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: profile(:,:), kpterms(:,:,:)
    integer, intent(in) :: n_occ
    integer, intent(out) :: z2, status
    real(kind=dp), intent(out) :: min_gap

    integer :: N, dim_H, i_trim, istate, isite, ncol
    integer :: info, lwork, n_pairs, ipair, j_best, jstate
    real(kind=dp), parameter :: parity_sign_tol = 1.0e-8_dp
    real(kind=dp) :: pair_sign, a_lat, dE_best, dE_try
    real(kind=dp) :: delta_trim(4)
    real(kind=dp) :: kx_trim(4), ky_trim(4)
    type(wavevector) :: wv

    complex(kind=dp), allocatable :: HT(:,:)
    real(kind=dp), allocatable :: evals(:), rwork(:)
    complex(kind=dp), allocatable :: work(:)

    ! Per-band parity eigenvalues and pairing workspace
    real(kind=dp), allocatable :: band_parity(:)
    integer, allocatable :: partner(:)
    logical, allocatable :: paired(:)

    z2 = 0
    min_gap = huge(1.0_dp)
    status = topo_status_invalid

    N = size(profile, 1)
    dim_H = 8 * N

    if (N < 1 .or. size(profile, 2) < 1) return
    if (size(kpterms, 1) < N .or. size(kpterms, 2) < N) return
    if (n_occ < 2 .or. mod(n_occ, 2) /= 0 .or. n_occ >= dim_H) return

    ncol = size(profile, 2)
    do isite = 1, N / 2
      if (maxval(abs(profile(isite, 1:ncol) - profile(N + 1 - isite, 1:ncol))) > 1.0e-8_dp) then
        status = topo_status_asymmetric
        return
      end if
    end do

    a_lat = qw_lattice_constant(cfg)
    if (a_lat <= 0.0_dp) return

    kx_trim = [0.0_dp, pi_dp / a_lat, 0.0_dp, pi_dp / a_lat]
    ky_trim = [0.0_dp, 0.0_dp, pi_dp / a_lat, pi_dp / a_lat]
    delta_trim = 1.0_dp

    ! Workspace for zheev: query first, then allocate
    allocate(HT(dim_H, dim_H))
    allocate(evals(dim_H))
    allocate(rwork(max(1, 3*dim_H - 2)))
    allocate(work(1))

    ! Workspace query
    call zheev('V', 'U', dim_H, HT, dim_H, evals, work, -1, rwork, info)
    if (info /= 0) then
      status = topo_status_lapack
      z2 = 0
      deallocate(HT, evals, rwork, work)
      return
    end if
    lwork = max(1, int(real(work(1), kind=dp)))
    deallocate(work)
    allocate(work(lwork))

    allocate(band_parity(n_occ))
    allocate(partner(n_occ))
    allocate(paired(n_occ))

    do i_trim = 1, 4
      ! Set up wavevector at this TRIM point
      wv%kx = kx_trim(i_trim)
      wv%ky = ky_trim(i_trim)
      wv%kz = 0.0_dp

      ! Build QW Hamiltonian at this TRIM
      HT = cmplx(0.0_dp, 0.0_dp, kind=dp)
      call ZB8bandQW(HT, wv, profile, kpterms, cfg=cfg)

      ! Diagonalize: eigenvalues in evals, eigenvectors overwrite HT
      call zheev('V', 'U', dim_H, HT, dim_H, evals, work, lwork, rwork, info)
      if (info /= 0) then
        status = topo_status_lapack
        z2 = 0
        deallocate(HT, evals, rwork, work, band_parity, partner, paired)
        return
      end if

      min_gap = min(min_gap, evals(n_occ + 1) - evals(n_occ))

      ! Step 1: Compute individual parity eigenvalue for each occupied band
      do istate = 1, n_occ
        band_parity(istate) = qw_inversion_expectation(HT(:, istate), N)
      end do

      ! Step 2: Greedy nearest-energy pairing of Kramers partners.
      ! At TRIM, Kramers partners share the same parity eigenvalue.
      ! For each unpaired band, find the closest-energy unpaired band with
      ! the same parity sign.  This is robust to zheev reordering within
      ! near-degenerate subspaces.
      paired = .false.
      partner = 0
      n_pairs = 0

      do istate = 1, n_occ
        if (paired(istate)) cycle
        j_best = 0
        dE_best = huge(1.0_dp)
        do jstate = istate + 1, n_occ
          if (paired(jstate)) cycle
          ! Kramers partners at TRIM have the same parity sign
          if (band_parity(istate) * band_parity(jstate) <= 0.0_dp) cycle
          dE_try = abs(evals(istate) - evals(jstate))
          if (dE_try < dE_best) then
            dE_best = dE_try
            j_best = jstate
          end if
        end do
        if (j_best == 0) then
          ! No valid partner found — degeneracy structure broken
          status = topo_status_invalid
          z2 = 0
          deallocate(HT, evals, rwork, work, band_parity, partner, paired)
          return
        end if
        paired(istate) = .true.
        paired(j_best) = .true.
        partner(istate) = j_best
        partner(j_best) = istate
        n_pairs = n_pairs + 1
      end do

      if (n_pairs /= n_occ / 2) then
        status = topo_status_invalid
        z2 = 0
        deallocate(HT, evals, rwork, work, band_parity, partner, paired)
        return
      end if

      ! Step 3: Compute parity product over Kramers pairs
      do ipair = 1, n_occ
        ! Only process the first member of each pair to avoid double-counting
        if (partner(ipair) <= ipair) cycle
        pair_sign = qw_pair_inversion_sign(HT(:, ipair), HT(:, partner(ipair)), N)
        if (abs(pair_sign) <= parity_sign_tol) then
          status = topo_status_invalid
          z2 = 0
          deallocate(HT, evals, rwork, work, band_parity, partner, paired)
          return
        end if
        if (pair_sign < 0.0_dp) then
          delta_trim(i_trim) = delta_trim(i_trim) * (-1.0_dp)
        end if
      end do
    end do

    z2 = z2_from_trim_parities(delta_trim)
    status = topo_status_ok

    deallocate(HT, evals, rwork, work, band_parity, partner, paired)

  end subroutine compute_z2_fukane_qw_result

  function compute_z2_fukane_qw(cfg, profile, kpterms, n_occ) result(z2)
    implicit none
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: profile(:,:), kpterms(:,:,:)
    integer, intent(in) :: n_occ
    integer :: z2
    integer :: status
    real(kind=dp) :: min_gap

    call compute_z2_fukane_qw_result(cfg, profile, kpterms, n_occ, z2, min_gap, status)
    if (status /= topo_status_ok) z2 = 0
  end function compute_z2_fukane_qw

  real(kind=dp) function qw_lattice_constant(cfg) result(a_lat)
    implicit none
    type(simulation_config), intent(in) :: cfg
    integer :: layer, i

    a_lat = 0.0_dp
    if (.not. allocated(cfg%params)) return
    if (size(cfg%params) < 1) return

    layer = max(1, (cfg%num_layers + 1) / 2)
    if (layer <= size(cfg%params)) a_lat = cfg%params(layer)%a0
    if (a_lat > 0.0_dp) return

    do i = 1, size(cfg%params)
      if (cfg%params(i)%a0 > 0.0_dp) then
        a_lat = cfg%params(i)%a0
        return
      end if
    end do
  end function qw_lattice_constant

  ! ==============================================================================
  ! Z2 phase diagram via gap sweep over (B, mu) parameter space.
  !
  ! Sweeps a linear grid of nB x nMu points, computing Z2 invariant and
  ! minimum bulk gap at each point. Detects phase transitions where Z2 flips
  ! between adjacent B values.
  !
  ! Uses evaluator selected by cfg%topo%sweep_model.
  !
  ! OUTPUT:
  !   z2_map(nMu, nB)       : Z2 invariant (0 or 1) at each grid point
  !   gap_map(nMu, nB)      : minimum gap at each grid point
  !   transitions(:, 2)     : (B, mu) transition midpoints
  ! ==============================================================================
  subroutine compute_z2_gap_sweep(cfg, B_min, B_max, nB, mu_min, mu_max, nMu, &
                                   gap_threshold, z2_map, gap_map, transitions)
    implicit none
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), intent(in) :: B_min, B_max, mu_min, mu_max, gap_threshold
    integer, intent(in) :: nB, nMu
    integer, allocatable, intent(out) :: z2_map(:,:)
    real(kind=dp), allocatable, intent(out) :: gap_map(:,:)
    real(kind=dp), allocatable, intent(out) :: transitions(:,:)

    integer :: iB, iMu, status_eval
    real(kind=dp) :: dB, dmu, B_val, mu_val

    if (nB < 1 .or. nMu < 1) then
      allocate(z2_map(0,0), gap_map(0,0), transitions(0,2))
      return
    end if

    allocate(z2_map(nMu, nB))
    allocate(gap_map(nMu, nB))
    z2_map = 0
    gap_map = 0.0_dp

    dB = 0.0_dp
    if (nB > 1) dB = (B_max - B_min) / real(nB - 1, kind=dp)

    dmu = 0.0_dp
    if (nMu > 1) dmu = (mu_max - mu_min) / real(nMu - 1, kind=dp)

    do iB = 1, nB
      B_val = B_min + real(iB - 1, kind=dp) * dB
      do iMu = 1, nMu
        mu_val = mu_min + real(iMu - 1, kind=dp) * dmu
        select case (trim(cfg%topo%sweep_model))
        case ('bhz_analytic')
          call eval_bhz_analytic(B_val, mu_val, cfg, z2_map(iMu, iB), &
            gap_map(iMu, iB), status_eval)
        case default
          print *, 'ERROR: compute_z2_gap_sweep supports sweep_model=bhz_analytic only'
          print *, '       Use topologicalAnalysis sweep mode for QW Fu-Kane or wire BdG sweeps.'
          stop 1
        end select
        if (status_eval /= topo_status_ok) then
          print *, 'ERROR: topology gap sweep evaluator failed for model ', &
            trim(cfg%topo%sweep_model)
          stop 1
        end if
      end do
    end do

    call detect_z2_transitions(z2_map, gap_map, B_min, B_max, mu_min, mu_max, &
      gap_threshold, transitions)

  end subroutine compute_z2_gap_sweep

  subroutine detect_z2_transitions(z2_map, gap_map, B_min, B_max, mu_min, mu_max, &
                                   gap_threshold, transitions)
    implicit none
    integer, contiguous, intent(in) :: z2_map(:,:)
    real(kind=dp), contiguous, intent(in) :: gap_map(:,:)
    real(kind=dp), intent(in) :: B_min, B_max, mu_min, mu_max, gap_threshold
    real(kind=dp), allocatable, intent(out) :: transitions(:,:)

    integer :: nMu, nB, iMu, iB, n_trans, max_trans
    real(kind=dp) :: dB, dmu
    real(kind=dp), allocatable :: trans_raw(:,:), trans_trim(:,:)

    nMu = size(z2_map, 1)
    nB = size(z2_map, 2)
    if (nMu < 1 .or. nB < 1) then
      allocate(transitions(0, 2))
      return
    end if

    dB = 0.0_dp
    if (nB > 1) dB = (B_max - B_min) / real(nB - 1, kind=dp)
    dmu = 0.0_dp
    if (nMu > 1) dmu = (mu_max - mu_min) / real(nMu - 1, kind=dp)

    max_trans = nMu * max(0, nB - 1) + nB * max(0, nMu - 1)
    allocate(trans_raw(max_trans, 2))
    n_trans = 0

    do iMu = 1, nMu
      do iB = 1, nB - 1
        if (is_z2_transition(z2_map(iMu, iB), z2_map(iMu, iB + 1))) then
          n_trans = n_trans + 1
          trans_raw(n_trans, 1) = B_min + (real(iB - 1, kind=dp) + 0.5_dp) * dB
          trans_raw(n_trans, 2) = mu_min + real(iMu - 1, kind=dp) * dmu
        end if
      end do
    end do

    do iMu = 1, nMu - 1
      do iB = 1, nB
        if (is_z2_transition(z2_map(iMu, iB), z2_map(iMu + 1, iB))) then
          n_trans = n_trans + 1
          trans_raw(n_trans, 1) = B_min + real(iB - 1, kind=dp) * dB
          trans_raw(n_trans, 2) = mu_min + (real(iMu - 1, kind=dp) + 0.5_dp) * dmu
        end if
      end do
    end do

    allocate(trans_trim(n_trans, 2))
    if (n_trans > 0) trans_trim = trans_raw(1:n_trans, :)
    call move_alloc(trans_trim, transitions)
    deallocate(trans_raw)

  end subroutine detect_z2_transitions

  elemental logical function is_z2_transition(z2_a, z2_b) result(is_transition)
    implicit none
    integer, intent(in) :: z2_a, z2_b

    is_transition = (z2_a /= z2_b)
  end function is_z2_transition

  subroutine eval_bhz_analytic(B_val, mu_val, cfg, z2, gap, status)
    implicit none
    real(kind=dp), intent(in) :: B_val, mu_val
    type(simulation_config), intent(in) :: cfg
    integer, intent(out) :: z2, status
    real(kind=dp), intent(out) :: gap

    real(kind=dp) :: M_eff

    M_eff = 1.0e-3_dp * cfg%topo%bhz_M + B_val - mu_val
    gap = abs(2.0_dp * M_eff)
    if (M_eff < 0.0_dp) then
      z2 = 1
    else
      z2 = 0
    end if
    status = topo_status_ok
  end subroutine eval_bhz_analytic

end module topological_analysis
