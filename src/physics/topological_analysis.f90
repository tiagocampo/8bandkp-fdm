module topological_analysis

  use definitions
  use sparse_matrices
  use linalg
  implicit none
  private

  public :: compute_chern_qwz
  public :: compute_berry_curvature
  public :: compute_berry_curvature_lattice
  public :: compute_hall_conductance
  public :: compute_z2_gap
  public :: compute_z2_fukane
  public :: extract_edge_states
  public :: compute_majorana_profile
  public :: build_bhz_wire_hamiltonian
  public :: bhz_wire_params
  public :: extract_edge_states_wire
  public :: fit_exponential_decay
  public :: compute_phase_diagram
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

  function compute_z2_fukane(H_eigs, params) result(z2)
    implicit none
    real(kind=dp), intent(in) :: H_eigs(:)
    class(*), intent(in) :: params
    integer :: z2

    print *, 'WARNING: compute_z2_fukane is not yet implemented. Returning 0.'
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

end module topological_analysis
