module topological_analysis

  use definitions
  use sparse_matrices
  use linalg, only: zheev, zheevd, zgetrf
  use hamiltonianConstructor, only: ZB8bandQW
  implicit none
  private

  public :: compute_chern_qwz
  public :: compute_berry_curvature
  public :: compute_berry_curvature_lattice
  public :: compute_hall_conductance
  public :: compute_z2_gap
  public :: compute_z2_fukane
  public :: compute_z2_fukane_qw
  public :: extract_edge_states
  public :: compute_majorana_profile
  public :: build_bhz_wire_hamiltonian
  public :: bhz_wire_params
  public :: extract_edge_states_wire
  public :: fit_exponential_decay
  public :: compute_phase_diagram
  public :: compute_z2_gap_sweep
  public :: gap_closing_detect

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

    dk = 2.0_dp * acos(-1.0_dp) / real(nk, kind=dp)
    total_flux = 0.0_dp

    allocate(evecs(nk, nk, 2))

    do j = 1, nk
      ky = -acos(-1.0_dp) + real(j-1, kind=dp) * dk
      do i = 1, nk
        kx = -acos(-1.0_dp) + real(i-1, kind=dp) * dk

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

    C = nint(total_flux / (2.0_dp * acos(-1.0_dp)))
    deallocate(evecs)
  end function compute_chern_qwz

  subroutine diag_2x2(H, eval, eigvec)
    implicit none
    real(kind=dp), intent(in) :: H(2,2)
    real(kind=dp), intent(out) :: eval(2), eigvec(2,2)
    real(kind=dp) :: tr, det, sqrt_term

    tr = H(1,1) + H(2,2)
    det = H(1,1)*H(2,2) - H(1,2)*H(2,1)
    sqrt_term = sqrt(max(0.0_dp, tr*tr/4.0_dp - det))
    eval(1) = tr/2.0_dp + sqrt_term; eval(2) = tr/2.0_dp - sqrt_term

    if (abs(H(1,2)) > 1e-12_dp) then
      eigvec(1,1) = H(1,2); eigvec(2,1) = eval(1) - H(1,1)
      eigvec(1,2) = H(1,2); eigvec(2,2) = eval(2) - H(1,1)
    else
      eigvec(1,1) = 1.0_dp; eigvec(2,1) = 0.0_dp
      eigvec(1,2) = 0.0_dp; eigvec(2,2) = 1.0_dp
    end if
    if (sum(eigvec(:,1)**2) > 1.0e-30_dp) then
      eigvec(:,1) = eigvec(:,1) / sqrt(sum(eigvec(:,1)**2))
    else
      eigvec(1,1) = 1.0_dp; eigvec(2,1) = 0.0_dp
    end if
    if (sum(eigvec(:,2)**2) > 1.0e-30_dp) then
      eigvec(:,2) = eigvec(:,2) / sqrt(sum(eigvec(:,2)**2))
    else
      eigvec(1,2) = 0.0_dp; eigvec(2,2) = 1.0_dp
    end if
  end subroutine diag_2x2

  function compute_berry_curvature(evecs_k, kx_arr, ky_arr, n_occ) result(Omega)
    implicit none
    complex(kind=dp), intent(in) :: evecs_k(:,:,:,:)
    real(kind=dp), intent(in) :: kx_arr(:), ky_arr(:)
    integer, intent(in) :: n_occ
    real(kind=dp), allocatable :: Omega(:,:)

    ! Delegate to lattice method
    Omega = compute_berry_curvature_lattice(evecs_k, kx_arr, ky_arr, n_occ)
  end function compute_berry_curvature

  function compute_berry_curvature_lattice(evecs_k, kx_arr, ky_arr, n_occ) result(Omega)
    implicit none
    complex(kind=dp), intent(in) :: evecs_k(:,:,:,:)  ! (basis, n_occ, nkx, nky)
    real(kind=dp), intent(in) :: kx_arr(:), ky_arr(:)
    integer, intent(in) :: n_occ
    real(kind=dp), allocatable :: Omega(:,:)

    integer :: nkx, nky, i, j, ip1, jp1
    complex(kind=dp) :: M_xx, M_yx, M_xy, M_yy
    real(kind=dp) :: dkx, dky, dA

    nkx = size(kx_arr)
    nky = size(ky_arr)
    allocate(Omega(nkx, nky))
    Omega = 0.0_dp

    dkx = kx_arr(min(2, nkx)) - kx_arr(1)
    dky = ky_arr(min(2, nky)) - ky_arr(1)
    dA = dkx * dky

    do j = 1, nky
      jp1 = mod(j, nky) + 1
      do i = 1, nkx
        ip1 = mod(i, nkx) + 1

        block
          complex(kind=dp), allocatable :: overlap(:,:)
          integer :: m, p

          allocate(overlap(n_occ, n_occ))

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

          deallocate(overlap)
        end block

      end do
    end do

  end function compute_berry_curvature_lattice

  function det_small(A, n) result(d)
    implicit none
    complex(kind=dp), intent(in) :: A(:,:)
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

  function compute_hall_conductance(C) result(sigma_xy)
    implicit none
    integer, intent(in) :: C
    real(kind=dp) :: sigma_xy

    sigma_xy = real(C, kind=dp)
  end function compute_hall_conductance

  function compute_z2_gap(N, eigenvalues, gap_threshold) result(z2)
    implicit none
    integer, intent(in) :: N
    real(kind=dp), intent(in) :: eigenvalues(:)
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

  function compute_z2_fukane(cfg, profile, kpterms, n_occ) result(z2)
    implicit none
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: profile(:,:), kpterms(:,:,:)
    integer, intent(in) :: n_occ
    integer :: z2

    ! Delegate to the full Fu-Kane QW implementation
    z2 = compute_z2_fukane_qw(cfg, profile, kpterms, n_occ)

  end function compute_z2_fukane

  function extract_edge_states(eigenvalues, eigenvectors, grid, window) result(edge_xi)
    implicit none
    real(kind=dp), intent(in) :: eigenvalues(:)
    complex(kind=dp), intent(in) :: eigenvectors(:,:)
    type(spatial_grid), intent(in) :: grid
    real(kind=dp), intent(in) :: window
    real(kind=dp) :: edge_xi

    real(kind=dp), allocatable :: edge_info(:)

    edge_info = extract_edge_states_wire(eigenvalues, eigenvectors, grid, window)
    edge_xi = edge_info(1)

  end function extract_edge_states

  function extract_edge_states_wire(eigenvalues, eigenvectors, grid, window) result(edge_info)
    implicit none
    real(kind=dp), intent(in) :: eigenvalues(:)
    complex(kind=dp), intent(in) :: eigenvectors(:,:)
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

  subroutine fit_exponential_decay(density, x, xi, success)
    implicit none
    real(kind=dp), intent(in) :: density(:)
    real(kind=dp), intent(in) :: x(:)
    real(kind=dp), intent(out) :: xi
    logical, intent(out) :: success

    integer :: n, peak_idx, tail_start, i
    real(kind=dp) :: rho_peak, tol
    real(kind=dp) :: sum_y_log, sum_x, sum_x2, sum_xy, denom, slope
    integer :: n_fit
    real(kind=dp) :: x_start

    n = size(density)
    success = .false.
    xi = 0.0_dp

    if (n < 3 .or. size(x) /= n) return

    peak_idx = maxloc(density, dim=1)

    rho_peak = density(peak_idx)
    tol = 0.1_dp
    tail_start = 0

    do i = peak_idx + 1, n
      if (density(i) < tol * rho_peak) then
        tail_start = i
        exit
      end if
    end do

    if (tail_start == 0) then
      do i = peak_idx - 1, 1, -1
        if (density(i) < tol * rho_peak) then
          tail_start = i
          exit
        end if
      end do
    end if

    if (tail_start == 0 .or. tail_start >= n) return

    x_start = x(tail_start)
    n_fit = n - tail_start + 1

    sum_y_log = 0.0_dp
    sum_x = 0.0_dp
    sum_x2 = 0.0_dp
    sum_xy = 0.0_dp

    do i = tail_start, n
      if (density(i) > 1.0e-14_dp) then
        sum_y_log = sum_y_log + log(density(i))
        sum_x = sum_x + (x(i) - x_start)
        sum_x2 = sum_x2 + (x(i) - x_start)**2
        sum_xy = sum_xy + (x(i) - x_start) * log(density(i))
      end if
    end do

    denom = sum_x2 - sum_x**2 / real(n_fit, kind=dp)
    if (abs(denom) < 1.0e-14_dp .or. n_fit < 3) return

    slope = (sum_xy - sum_x * sum_y_log / real(n_fit, kind=dp)) / denom
    if (abs(slope) < tiny(1.0_dp)) return
    xi = abs(-1.0_dp / slope)
    if (xi > 0.0_dp) then
      success = .true.
    end if

  end subroutine fit_exponential_decay

  function compute_majorana_profile(evec_bdg, grid, energy_tol, half_n_in, profile) result(xi)
    implicit none
    complex(kind=dp), intent(in) :: evec_bdg(:)
    type(spatial_grid), intent(in) :: grid
    real(kind=dp), intent(in) :: energy_tol
    integer, intent(in) :: half_n_in
    real(kind=dp), intent(out), optional :: profile(:)
    real(kind=dp) :: xi

    integer :: ngrid, half_n, nspatial
    real(kind=dp), allocatable :: rho(:)
    real(kind=dp), allocatable :: xx(:)
    integer :: i, ib, peak_idx, i_start
    real(kind=dp) :: xi_est, norm_factor
    real(kind=dp) :: sum_y_log, sum_x, sum_x2, sum_xy, denom, slope
    integer :: n_fit
    real(kind=dp) :: tol_norm

    ngrid = grid%npoints()
    half_n = half_n_in
    if (half_n * 2 > size(evec_bdg)) then
      print *, 'ERROR: compute_majorana_profile: eigenvector too small'
      print *, '  half_n=', half_n, ' size(evec)=', size(evec_bdg)
      xi = 0.0_dp
      return
    end if
    nspatial = half_n / 8

    allocate(rho(nspatial))
    do i = 1, nspatial
      rho(i) = 0.0_dp
      do ib = 1, 8
        rho(i) = rho(i) + abs(evec_bdg((i-1)*8 + ib))**2
        rho(i) = rho(i) + abs(evec_bdg(half_n + (i-1)*8 + ib))**2
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

    peak_idx = maxloc(rho, dim=1)

    if (allocated(grid%z)) then
      xx = grid%z
    else
      allocate(xx(nspatial))
      do i = 1, nspatial
        xx(i) = real(i - 1, kind=dp)
      end do
    end if

    tol_norm = 0.1_dp
    i_start = 0
    do i = peak_idx + 1, nspatial
      if (rho(i) < tol_norm * rho(peak_idx)) then
        i_start = i
        exit
      end if
    end do

    if (i_start == 0) then
      do i = peak_idx - 1, 1, -1
        if (rho(i) < tol_norm * rho(peak_idx)) then
          i_start = i
          exit
        end if
      end do
    end if

    if (i_start == 0 .or. i_start >= nspatial) then
      xi = 0.0_dp
      deallocate(rho)
      if (.not. allocated(grid%z)) deallocate(xx)
      return
    end if

    n_fit = nspatial - i_start + 1
    sum_y_log = 0.0_dp
    sum_x = 0.0_dp
    sum_x2 = 0.0_dp
    sum_xy = 0.0_dp

    do i = i_start, nspatial
      if (rho(i) > 1.0e-14_dp) then
        sum_y_log = sum_y_log + log(rho(i))
        sum_x = sum_x + (xx(i) - xx(i_start))
        sum_x2 = sum_x2 + (xx(i) - xx(i_start))**2
        sum_xy = sum_xy + (xx(i) - xx(i_start)) * log(rho(i))
      end if
    end do

    denom = sum_x2 - sum_x**2 / real(n_fit, kind=dp)
    if (abs(denom) < 1.0e-14_dp .or. n_fit < 3) then
      xi = 0.0_dp
    else
      slope = (sum_xy - sum_x * sum_y_log / real(n_fit, kind=dp)) / denom
      if (abs(slope) < tiny(1.0_dp)) then
        xi = 0.0_dp
      else
        xi_est = -1.0_dp / slope
        xi = abs(xi_est)
      end if
    end if

    deallocate(rho)
    if (.not. allocated(grid%z)) deallocate(xx)

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
          coo_col(nnz_offset) = (i-1)*4 + 5 - row
          if (row == 1 .or. row == 4) then
            coo_vals(nnz_offset) = cmplx(A_over_2dz, 0.0_dp, kind=dp)
          else
            coo_vals(nnz_offset) = cmplx(-A_over_2dz, 0.0_dp, kind=dp)
          end if
        end if

        if (i > 1) then
          nnz_offset = nnz_offset + 1
          coo_row(nnz_offset) = (i-1)*4 + row
          coo_col(nnz_offset) = (i-2)*4 + row
          if (row == 1 .or. row == 4) then
            coo_vals(nnz_offset) = cmplx(-A_over_2dz, 0.0_dp, kind=dp)
          else
            coo_vals(nnz_offset) = cmplx(A_over_2dz, 0.0_dp, kind=dp)
          end if
        end if

        if (i < N) then
          nnz_offset = nnz_offset + 1
          coo_row(nnz_offset) = (i-1)*4 + row
          coo_col(nnz_offset) = (i-1)*4 + mod(row, 4) + 1
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
  ! Gap closing detection at a single (B, mu) point.
  ! Returns the minimum eigenvalue spacing (gap) around the Fermi level.
  ! ==============================================================================
  function gap_closing_detect(eigenvalues, mu, gap_threshold) result(is_gap_closing)
    implicit none
    real(kind=dp), intent(in) :: eigenvalues(:)
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

    ! Sweep over B and mu
    do iB = 1, nB
      B_val = B_min + real(iB - 1, kind=dp) * dB
      do iMu = 1, nMu
        mu_val = mu_min + real(iMu - 1, kind=dp) * dmu

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
        if (is_gap_closing) z2_array(iB, iMu) = 1

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
  function compute_z2_fukane_qw(cfg, profile, kpterms, n_occ) result(z2)
    implicit none
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: profile(:,:), kpterms(:,:,:)
    integer, intent(in) :: n_occ
    integer :: z2

    ! Parity of each band in the 8-band zinc-blende basis
    ! Bands 1-6 (valence): even under inversion (+1)
    ! Bands 7-8 (conduction): odd under inversion (-1)
    real(kind=dp), parameter :: P_band(8) = [+1.0_dp, +1.0_dp, +1.0_dp, &
                                              +1.0_dp, +1.0_dp, +1.0_dp, &
                                             -1.0_dp, -1.0_dp]

    integer :: N, dim_H, i_trim, istate, isite, iband, irow
    integer :: info, lwork
    real(kind=dp) :: delta_trim, parity_n
    real(kind=dp) :: product_non_gamma
    real(kind=dp) :: kx_trim(4), ky_trim(4)
    type(wavevector) :: wv

    complex(kind=dp), allocatable :: HT(:,:)
    real(kind=dp), allocatable :: evals(:), rwork(:)
    complex(kind=dp), allocatable :: work(:)

    z2 = 0

    N = size(profile, 1)
    dim_H = 8 * N

    if (n_occ < 1 .or. N < 1) return
    if (n_occ > dim_H) return

    ! Lattice constant (default to 6.0 Angstrom for III-V materials)
    ! Use a default if not set in config
    block
      real(kind=dp) :: a_lat
      real(kind=dp) :: pi_val

      pi_val = acos(-1.0_dp)
      a_lat = 6.0_dp  ! Angstrom, typical III-V lattice constant

      ! 4 TRIM points in the 2D Brillouin zone
      kx_trim = [0.0_dp, pi_val / a_lat, 0.0_dp, pi_val / a_lat]
      ky_trim = [0.0_dp, 0.0_dp, pi_val / a_lat, pi_val / a_lat]
    end block

    ! Workspace for zheev: query first, then allocate
    allocate(HT(dim_H, dim_H))
    allocate(evals(dim_H))
    allocate(rwork(max(1, 3*dim_H - 2)))
    allocate(work(1))

    ! Workspace query
    call zheev('V', 'U', dim_H, HT, dim_H, evals, work, -1, rwork, info)
    if (info /= 0) then
      print *, 'WARNING: compute_z2_fukane_qw: zheev workspace query failed, info=', info
      z2 = 0
      deallocate(HT, evals, rwork, work)
      return
    end if
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))

    product_non_gamma = 1.0_dp

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
        print *, 'WARNING: compute_z2_fukane_qw: zheev failed at TRIM ', i_trim, ' info=', info
        z2 = 0
        deallocate(HT, evals, rwork, work)
        return
      end if

      ! Compute parity product over occupied states at this TRIM
      ! delta_i = prod_{n=1}^{n_occ} xi_n
      delta_trim = 1.0_dp

      do istate = 1, n_occ
        ! For each occupied eigenstate, compute parity eigenvalue.
        ! The parity of a QW eigenstate is determined by the dominant band
        ! character. In the 8-band zinc-blende basis:
        !   VB bands (1-6): even under inversion (P=+1)
        !   CB bands (7-8): odd under inversion (P=-1)
        ! For a state with mixed character, the parity eigenvalue is:
        !   xi_n = sign(prod_b P_band(b)^{weight_b})
        ! which equals (-1)^{total_CB_weight}. For a normalized state with
        ! definite parity at a TRIM, CB weight > 0.5 gives xi=-1.
        block
          real(kind=dp) :: cb_weight, vb_weight

          cb_weight = 0.0_dp
          vb_weight = 0.0_dp

          do isite = 1, N
            do iband = 1, 8
              irow = (isite - 1) * 8 + iband
              if (P_band(iband) < 0.0_dp) then
                cb_weight = cb_weight + abs(HT(irow, istate))**2
              else
                vb_weight = vb_weight + abs(HT(irow, istate))**2
              end if
            end do
          end do

          ! Parity eigenvalue: if CB character dominates, parity is odd (-1)
          if (cb_weight > vb_weight) then
            parity_n = -1.0_dp
          else
            parity_n = +1.0_dp
          end if
        end block

        delta_trim = delta_trim * parity_n
      end do

      ! The Fu-Kane formula uses product of delta_i for TRIMs excluding Gamma.
      ! Gamma is i_trim=1. Accumulate product of TRIMs 2, 3, 4.
      if (i_trim > 1) then
        product_non_gamma = product_non_gamma * delta_trim
      end if
    end do

    ! Z2 = 1 if product of delta at non-Gamma TRIMs is negative
    if (product_non_gamma < 0.0_dp) z2 = 1

    deallocate(HT, evals, rwork, work)

  end function compute_z2_fukane_qw

  ! ==============================================================================
  ! Z2 phase diagram via gap sweep over (B, mu) parameter space.
  !
  ! Sweeps a linear grid of nB x nMu points, computing Z2 invariant and
  ! minimum bulk gap at each point. Detects phase transitions where Z2 flips
  ! between adjacent B values.
  !
  ! For wire mode (confinement==2): uses BHZ M-sign heuristic (z2=1 when B>0).
  ! For other modes: z2=0 placeholder (future QW implementation).
  !
  ! OUTPUT:
  !   z2_map(nMu, nB)       : Z2 invariant (0 or 1) at each grid point
  !   gap_map(nMu, nB)      : placeholder gap values (0.0 for now)
  !   transitions(:, 2)     : (B, mu) pairs where Z2 flips, detected along
  !                           the central mu row by comparing adjacent B values
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

    integer :: iB, iMu, j_central, n_trans
    real(kind=dp) :: dB, dmu, B_val
    real(kind=dp), allocatable :: trans_raw(:,:)

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

    ! Fill Z2 map based on confinement mode
    do iB = 1, nB
      B_val = B_min + real(iB - 1, kind=dp) * dB
      do iMu = 1, nMu
        if (cfg%confinement == 2) then
          ! Wire mode: BHZ heuristic — topological when M flips sign.
          ! For the standard BHZ model, M > 0 is trivial, M < 0 is topological.
          ! In the Zeeman-dominated regime, B > 0 drives M < 0 => topological.
          if (B_val > 0.0_dp) then
            z2_map(iMu, iB) = 1
          else
            z2_map(iMu, iB) = 0
          end if
        else
          ! Bulk / QW placeholder: trivial phase
          z2_map(iMu, iB) = 0
        end if
        ! gap_map remains 0.0 (placeholder for future eigenvalue-based computation)
      end do
    end do

    ! Detect transitions along B sweep at central mu
    j_central = (nMu + 1) / 2
    allocate(trans_raw(nB, 2))
    n_trans = 0

    do iB = 1, nB - 1
      if (z2_map(j_central, iB) /= z2_map(j_central, iB + 1)) then
        n_trans = n_trans + 1
        trans_raw(n_trans, 1) = B_min + real(iB - 1, kind=dp) * dB + dB / 2.0_dp
        trans_raw(n_trans, 2) = mu_min + real(j_central - 1, kind=dp) * dmu
      end if
    end do

    ! Copy to output array
    allocate(transitions(n_trans, 2))
    if (n_trans > 0) then
      transitions(1:n_trans, :) = trans_raw(1:n_trans, :)
    end if
    deallocate(trans_raw)

  end subroutine compute_z2_gap_sweep

end module topological_analysis
