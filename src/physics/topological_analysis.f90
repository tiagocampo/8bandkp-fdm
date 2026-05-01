module topological_analysis

  use definitions
  use sparse_matrices
  use linalg
  use eigensolver
  implicit none
  private

  public :: compute_chern_qwz
  public :: compute_berry_curvature
  public :: compute_hall_conductance
  public :: compute_z2_gap
  public :: compute_z2_fukane
  public :: extract_edge_states
  public :: compute_majorana_profile
  public :: build_bhz_wire_hamiltonian
  public :: bhz_wire_params
  public :: extract_edge_states_wire
  public :: fit_exponential_decay

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
    real(kind=dp) :: kx, ky, dk, total_flux
    complex(kind=dp) :: Ux, Uy, prod
    real(kind=dp) :: H(2,2), eval(2), eigvec(2,2)
    real(kind=dp) :: sin_kx_val, sin_ky_val, cos_kx_val, cos_ky_val, mz_val
    complex(kind=dp) :: ev1_i_j, ev2_i_j, ev1_ip1_j, ev2_ip1_j, ev1_i_jp1, ev2_i_jp1

    dk = 2.0_dp * acos(-1.0_dp) / real(nk, kind=dp)
    total_flux = 0.0_dp

    allocate(evecs(nk, nk, 2))

    do j = 1, nk
      ky = -acos(-1.0_dp) + real(j-1, kind=dp) * dk
      do i = 1, nk
        kx = -acos(-1.0_dp) + real(i-1, kind=dp) * dk
        sin_kx_val = sin(kx); sin_ky_val = sin(ky)
        cos_kx_val = cos(kx); cos_ky_val = cos(ky)
        mz_val = u + cos_kx_val + cos_ky_val

        H(1,1) = mz_val; H(1,2) = sin_kx_val; H(2,1) = sin_ky_val
        H(2,2) = -mz_val

        call diag_2x2(H, eval, eigvec)
        evecs(i,j,1) = cmplx(eigvec(1,1), eigvec(1,2), kind=dp)
        evecs(i,j,2) = cmplx(eigvec(2,1), eigvec(2,2), kind=dp)
      end do
    end do

    do j = 1, nk
      jp1 = mod(j, nk) + 1
      do i = 1, nk
        ip1 = mod(i, nk) + 1

        ev1_i_j = evecs(i,j,1)
        ev2_i_j = evecs(i,j,2)
        ev1_ip1_j = evecs(ip1,j,1)
        ev2_ip1_j = evecs(ip1,j,2)
        Ux = conjg(ev1_i_j) * ev1_ip1_j + conjg(ev2_i_j) * ev2_ip1_j

        ev1_i_jp1 = evecs(i,jp1,1)
        ev2_i_jp1 = evecs(i,jp1,2)
        Uy = conjg(ev1_i_j) * ev1_i_jp1 + conjg(ev2_i_j) * ev2_i_jp1

        prod = Ux * (conjg(ev1_ip1_j) * ev1_i_jp1 + conjg(ev2_ip1_j) * ev2_i_jp1)
        prod = prod * (conjg(ev1_i_jp1) * ev1_ip1_j + conjg(ev2_i_jp1) * ev2_ip1_j)
        prod = prod * conjg(Uy)
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
    eigvec(:,1) = eigvec(:,1) / sqrt(sum(eigvec(:,1)**2))
    eigvec(:,2) = eigvec(:,2) / sqrt(sum(eigvec(:,2)**2))
  end subroutine diag_2x2

  function compute_berry_curvature(evecs_k, kx_arr, ky_arr, n_occ) result(Omega)
    implicit none
    complex(kind=dp), intent(in) :: evecs_k(:,:,:,:)
    real(kind=dp), intent(in) :: kx_arr(:), ky_arr(:)
    integer, intent(in) :: n_occ
    real(kind=dp), allocatable :: Omega(:,:)

    allocate(Omega(size(kx_arr), size(ky_arr)))
    Omega = 0.0_dp
  end function compute_berry_curvature

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

    if (n_in_gap > 0) z2 = 1

  end function compute_z2_gap

  function compute_z2_fukane(H_eigs, params) result(z2)
    implicit none
    real(kind=dp), intent(in) :: H_eigs(:)
    class(*), intent(in) :: params
    integer :: z2

    z2 = 0
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

    N = grid%npoints()
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
    real(kind=dp) :: sum_y_log, sum_x, sum_x2, sum_xy, denom
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

    xi = abs(-1.0_dp / ((sum_xy - sum_x * sum_y_log / real(n_fit, kind=dp)) / denom))
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
    real(kind=dp) :: sum_y_log, sum_x, sum_x2, sum_xy, denom
    integer :: n_fit
    real(kind=dp) :: tol_norm

    ngrid = grid%npoints()
    half_n = half_n_in
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
      profile = rho
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
      xi_est = -1.0_dp / ((sum_xy - sum_x * sum_y_log / real(n_fit, kind=dp)) / denom)
      xi = abs(xi_est)
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

    N = params%N
    dz = params%dz

    B_plus_D_over_dz2 = (params%B + params%D) / (dz * dz)
    A_over_2dz = params%A / (2.0_dp * dz)

    nnz_total = 4*N + 8*(N-1) + 8*(N-1)
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
          coo_col(nnz_offset) = (i-1)*4 + mod(row + 1, 4) + 1
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
          coo_col(nnz_offset) = i*4 + row
          coo_vals(nnz_offset) = cmplx(B_plus_D_over_dz2, 0.0_dp, kind=dp)
        end if

        if (i > 1) then
          nnz_offset = nnz_offset + 1
          coo_row(nnz_offset) = (i-1)*4 + row
          coo_col(nnz_offset) = (i-2)*4 + row
          coo_vals(nnz_offset) = cmplx(B_plus_D_over_dz2, 0.0_dp, kind=dp)
        end if
      end do
    end do

    call csr_build_from_coo(H_csr, 4*N, 4*N, nnz_offset, &
                            coo_row(1:nnz_offset), coo_col(1:nnz_offset), &
                            coo_vals(1:nnz_offset))

    deallocate(coo_vals, coo_row, coo_col)

  end subroutine build_bhz_wire_hamiltonian

end module topological_analysis
