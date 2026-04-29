module sc_loop
  ! Self-consistent Schrodinger-Poisson loop with DIIS acceleration.
  !
  ! The loop iterates:
  !   1. Build Hamiltonian with current profile (includes SC potential)
  !   2. Solve eigenproblem at each k_par
  !   3. Compute charge density from occupied states
  !   4. Solve Poisson equation for new potential
  !   5. Mix potential (linear + optional DIIS)
  !   6. Check convergence |Delta Phi|_inf < tolerance
  !   7. Update profile and repeat

  use definitions
  use hamiltonianConstructor
  use hamiltonian_wire, only: wire_coo_cache, ZB8bandGeneralized
  use charge_density
  use poisson
  use eigensolver, only: eigensolver_config, eigensolver_result, &
    & solve_sparse_evp, eigensolver_result_free
  use sparse_matrices, only: csr_matrix, csr_clone_structure, csr_free
  use linalg, only: zheevx, dgesv, ilaenv, dlamch
  implicit none

  private
  public :: self_consistent_loop
  public :: linear_mix
  public :: diis_extrapolate
  public :: find_fermi_level
  public :: build_epsilon
  public :: build_doping_charge
  public :: build_epsilon_2d
  public :: build_doping_charge_2d
  public :: apply_potential_to_profile
  public :: apply_potential_to_profile_2d
  public :: self_consistent_loop_wire
  public :: find_fermi_level_wire
  public :: map_layer_to_grid

  ! Unit conversion constants
  real(kind=dp), parameter :: ANGSTROM_TO_NM = 0.1_dp    ! 1 Angstrom = 0.1 nm
  real(kind=dp), parameter :: CM3_TO_PER_NM3 = 1.0e-21_dp  ! 1 cm^-3 -> 1/nm^3 (1 cm^3 = 10^21 nm^3)

  ! Fermi bisection parameters
  integer, parameter :: MAX_FERMI_BISECT = 60
  real(kind=dp), parameter :: FERMI_TOL = 1.0e-8_dp
  real(kind=dp), parameter :: FERMI_BOUND_PAD = 2.0_dp  ! eV padding for bisection bounds

contains

  ! ------------------------------------------------------------------
  ! Main self-consistent loop for QW simulations
  ! ------------------------------------------------------------------
  subroutine self_consistent_loop(profile, cfg, kpterms, HT, eig, eigv, &
      & smallk, N, il, iuu, n_electron_out, n_hole_out)

    real(kind=dp), allocatable, intent(inout) :: profile(:,:)
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), intent(in) :: kpterms(:,:,:)
    complex(kind=dp), allocatable, intent(inout) :: HT(:,:)
    real(kind=dp), allocatable, intent(inout) :: eig(:,:)
    complex(kind=dp), allocatable, intent(inout) :: eigv(:,:,:)
    type(wavevector), allocatable, intent(in) :: smallk(:)
    integer, intent(in) :: N, il, iuu
    real(kind=dp), allocatable, intent(out), optional :: n_electron_out(:)
    real(kind=dp), allocatable, intent(out), optional :: n_hole_out(:)

    ! SC loop variables
    integer :: iter, niter, info
    integer :: num_subbands, nk_actual
    real(kind=dp) :: delta_phi, fermi_level

    ! Potential arrays
    real(kind=dp), allocatable :: phi_old(:), phi_new(:), phi_poisson(:)
    real(kind=dp), allocatable :: rho(:), epsilon(:)
    real(kind=dp), allocatable :: n_electron(:), n_hole(:)

    ! k_par grid
    real(kind=dp), allocatable :: kpar_grid(:)

    ! Eigensolver workspace
    real(kind=dp), allocatable :: rwork(:), eig_kpar(:,:)
    complex(kind=dp), allocatable :: work(:), eigv_kpar(:,:,:)
    integer, allocatable :: iwork(:), ifail_arr(:)
    integer :: M_out, lwork, k_idx
    real(kind=dp) :: abstol
    integer :: NB_val

    ! DIIS history
    real(kind=dp), allocatable :: phi_history(:,:), res_history(:,:)

    ! Profile backup
    real(kind=dp), allocatable :: profile_base(:,:)

    ! Doping charge
    real(kind=dp), allocatable :: rho_doping(:)

    ! Pre-allocated work arrays for Fermi bisection
    real(kind=dp), allocatable :: fermi_ne(:), fermi_nh(:)

    ! BC type as integer
    integer :: bc_type_int

    ! Local variables
    real(kind=dp) :: dz_val, kpar_max_val
    integer :: iz, nz
    logical :: sc_converged
    type(wavevector) :: wv

    nz = cfg%fdStep
    dz_val = cfg%dz
    num_subbands = iuu - il + 1
    niter = cfg%sc%max_iterations

    ! Build k_par grid (ensure odd for Simpson)
    nk_actual = cfg%sc%num_kpar
    if (mod(nk_actual, 2) == 0) nk_actual = nk_actual - 1

    ! Auto-determine kpar_max if not set
    kpar_max_val = cfg%sc%kpar_max
    if (kpar_max_val < tolerance) then
      kpar_max_val = cfg%waveVectorMax
      if (kpar_max_val < tolerance) kpar_max_val = 0.5_dp
    end if

    allocate(kpar_grid(nk_actual))
    do iz = 1, nk_actual
      kpar_grid(iz) = real(iz - 1, kind=dp) * kpar_max_val &
        & / real(nk_actual - 1, kind=dp)
    end do

    ! Allocate working arrays
    allocate(phi_old(nz), phi_new(nz), phi_poisson(nz))
    allocate(rho(nz), epsilon(nz))
    allocate(n_electron(nz), n_hole(nz))
    allocate(profile_base(nz, 3))
    allocate(rho_doping(nz))
    allocate(phi_history(nz, cfg%sc%diis_history))
    allocate(res_history(nz, cfg%sc%diis_history))
    allocate(fermi_ne(nz), fermi_nh(nz))

    phi_history = 0.0_dp
    res_history = 0.0_dp

    ! Save initial profile as base
    profile_base = profile

    ! Build dielectric and doping arrays
    call build_epsilon(epsilon, cfg, nz)
    call build_doping_charge(rho_doping, cfg, nz)

    ! Convert BC string to integer constant once
    bc_type_int = bc_from_string(cfg%sc%bc_type)

    ! Initialize potential
    phi_old = 0.0_dp
    phi_new = 0.0_dp
    sc_converged = .false.

    ! Fermi level
    if (cfg%sc%fermi_mode == 1) then
      fermi_level = cfg%sc%fermi_level
    else
      fermi_level = 0.0_dp
    end if

    ! Eigensolver setup
    NB_val = ILAENV(1, 'ZHETRD', 'UPLO', N, N, -1, -1)
    NB_val = max(NB_val, N)
    abstol = DLAMCH('P')

    allocate(rwork(7*N))
    allocate(iwork(5*N))
    allocate(ifail_arr(N))
    allocate(eig_kpar(num_subbands, nk_actual))
    allocate(eigv_kpar(N, num_subbands, nk_actual))

    eig_kpar = 0.0_dp
    eigv_kpar = cmplx(0.0_dp, 0.0_dp, kind=dp)

    ! Initial workspace query via zheevx
    allocate(work(1))
    lwork = -1
    call zheevx('V', 'I', 'U', N, HT, N, 0.0_dp, 0.0_dp, il, iuu, abstol, M_out, &
      & eig_kpar(:,1), HT, N, work, lwork, rwork, iwork, ifail_arr, info)
    if (info /= 0) then
      print *, 'Error: zheevx workspace query failed in SC loop, info =', info
      stop 1
    end if
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))

    ! --- Main SC loop ---
    print *, '=== Self-Consistent Loop Start ==='
    print *, '  max_iterations:', niter
    print *, '  tolerance:', cfg%sc%tolerance
    print *, '  mixing_alpha:', cfg%sc%mixing_alpha
    print *, '  diis_history:', cfg%sc%diis_history
    print *, '  num_kpar:', nk_actual
    print *, '  kpar_max:', kpar_max_val
    print *, '  bc_type:', trim(cfg%sc%bc_type)
    print *, '  temperature:', cfg%sc%temperature
    print *, ''

    do iter = 1, niter

      ! Step 1: Apply current potential to profile
      call apply_potential_to_profile(profile, profile_base, phi_old, nz)

      ! Step 2: Solve eigenproblem at each k_par
      do k_idx = 1, nk_actual
        wv%kx = kpar_grid(k_idx)
        wv%ky = 0.0_dp
        wv%kz = 0.0_dp

        call ZB8bandQW(HT, wv, profile, kpterms, cfg=cfg)

        call zheevx('V', 'I', 'U', N, HT, N, 0.0_dp, 0.0_dp, il, iuu, abstol, &
          & M_out, eig_kpar(:, k_idx), HT, N, work, lwork, rwork, &
          & iwork, ifail_arr, info)
        if (info /= 0) then
          print *, 'Error: diag failed at k_par =', kpar_grid(k_idx), 'info =', info
          stop 1
        end if

        eigv_kpar(:, :, k_idx) = HT(:, 1:num_subbands)
      end do

      ! Step 3: Find Fermi level (if charge neutrality mode)
      if (cfg%sc%fermi_mode == 0) then
        fermi_level = find_fermi_level(eig_kpar, eigv_kpar, kpar_grid, &
          & cfg, N, num_subbands, nk_actual, nz, dz_val, rho_doping, &
          & fermi_ne, fermi_nh)
      end if

      ! Step 4: Compute charge density
      call compute_charge_density_qw(n_electron, n_hole, eigv_kpar, &
        & eig_kpar, kpar_grid, fermi_level, cfg%sc%temperature, &
        & nz, num_subbands, nk_actual, cfg%numcb)

      ! Step 5: Build total charge and solve Poisson
      ! Convert density from cm^-3 to C/nm^3: rho = e * n * 1e-21 (1 cm^3 = 10^21 nm^3)
      do iz = 1, nz
        rho(iz) = e * (n_hole(iz) - n_electron(iz) + rho_doping(iz)) * CM3_TO_PER_NM3
      end do

      ! dz_val in Angstrom; Poisson solver uses nm (rho in C/nm^3, e0 in C/(V*nm))
      ! 1 Angstrom = 0.1 nm
      call solve_poisson(phi_poisson, rho, epsilon, dz_val * ANGSTROM_TO_NM, nz, &
        & cfg%sc%bc_left, cfg%sc%bc_right, bc_type_int)

      ! Step 6: Mix potential
      if (iter <= cfg%sc%diis_history .or. cfg%sc%diis_history < 2) then
        call linear_mix(phi_new, phi_old, phi_poisson, nz, cfg%sc%mixing_alpha)
      else
        call diis_extrapolate(phi_new, phi_old, phi_poisson, nz, &
          & phi_history, res_history, cfg%sc%diis_history, &
          & cfg%sc%mixing_alpha, iter)
      end if

      ! Update DIIS history
      call update_diis_history(phi_history, res_history, &
        & phi_old, phi_poisson, nz, cfg%sc%diis_history, iter)

      ! Step 7: Check convergence (phi in Volts; tolerance in eV; V == eV for electron potential
      delta_phi = maxval(abs(phi_new - phi_old))

      print '(A, I4, A, ES12.4, A, F10.4, A, ES12.4)', &
        & '  iter:', iter, '  |dPhi|:', delta_phi, &
        & '  mu:', fermi_level, '  max|rho|:', maxval(abs(rho))

      if (delta_phi < cfg%sc%tolerance) then
        print *, '=== SC loop converged at iteration', iter, '==='
        phi_old = phi_new
        sc_converged = .true.
        exit
      end if

      phi_old = phi_new

    end do

    if (.not. sc_converged) then
      print *, 'Warning: SC loop did not converge after', niter, 'iterations'
      print *, '  Final |dPhi|:', delta_phi
    end if

    ! Final update: apply converged potential to profile
    call apply_potential_to_profile(profile, profile_base, phi_old, nz)

    ! Copy final eigenvalues at k_par=0 back to eig
    eig(:, 1) = eig_kpar(:, 1)

    ! --- Copy charge densities to optional outputs before cleanup ---
    if (present(n_electron_out)) then
      if (allocated(n_electron)) then
        allocate(n_electron_out(nz))
        n_electron_out = n_electron
      end if
    end if
    if (present(n_hole_out)) then
      if (allocated(n_hole)) then
        allocate(n_hole_out(nz))
        n_hole_out = n_hole
      end if
    end if

    ! --- Cleanup ---
    deallocate(kpar_grid, phi_old, phi_new, phi_poisson)
    deallocate(rho, epsilon, n_electron, n_hole)
    deallocate(profile_base, rho_doping)
    deallocate(phi_history, res_history)
    deallocate(fermi_ne, fermi_nh)
    deallocate(rwork, iwork, ifail_arr, work)
    deallocate(eig_kpar, eigv_kpar)

  end subroutine self_consistent_loop


  ! ------------------------------------------------------------------
  ! Linear mixing: phi_new = (1-alpha)*phi_old + alpha*phi_poisson
  ! ------------------------------------------------------------------
  subroutine linear_mix(phi_new, phi_old, phi_poisson, N, alpha)
    integer, intent(in) :: N
    real(kind=dp), intent(out) :: phi_new(N)
    real(kind=dp), intent(in) :: phi_old(N), phi_poisson(N), alpha

    phi_new = (1.0_dp - alpha) * phi_old + alpha * phi_poisson

  end subroutine linear_mix


  ! ------------------------------------------------------------------
  ! DIIS extrapolation with linear blend
  ! ------------------------------------------------------------------
  subroutine diis_extrapolate(phi_new, phi_old, phi_poisson, N, &
      & phi_history, res_history, diis_len, alpha, iter)

    integer, intent(in) :: N, diis_len, iter
    real(kind=dp), intent(out) :: phi_new(N)
    real(kind=dp), intent(in) :: phi_old(N), phi_poisson(N)
    real(kind=dp), intent(in) :: phi_history(N, diis_len)
    real(kind=dp), intent(in) :: res_history(N, diis_len)
    real(kind=dp), intent(in) :: alpha

    integer :: m, i, j, info
    real(kind=dp), allocatable :: B(:,:), rhs(:), coeffs(:)
    real(kind=dp), allocatable :: phi_extrap(:)
    integer, allocatable :: ipiv(:)

    m = min(iter - 1, diis_len)
    if (m < 2) then
      call linear_mix(phi_new, phi_old, phi_poisson, N, alpha)
      return
    end if

    allocate(B(m+1, m+1), rhs(m+1), coeffs(m+1), ipiv(m+1))
    B = 0.0_dp
    rhs = 0.0_dp

    do i = 1, m
      do j = 1, m
        B(i, j) = dot_product(res_history(:, i), res_history(:, j))
      end do
      B(i, m+1) = -1.0_dp
      B(m+1, i) = -1.0_dp
    end do
    B(m+1, m+1) = 0.0_dp
    rhs(m+1) = -1.0_dp

    call dgesv(m+1, 1, B, m+1, ipiv, rhs, m+1, info)

    if (info /= 0) then
      print *, '  DIIS: linear system solve failed (info=', info, '). Falling back to linear mix.'
      deallocate(B, rhs, coeffs, ipiv)
      call linear_mix(phi_new, phi_old, phi_poisson, N, alpha)
      return
    end if

    coeffs = rhs(1:m+1)

    allocate(phi_extrap(N))
    phi_extrap = 0.0_dp
    do j = 1, m
      phi_extrap = phi_extrap + coeffs(j) * phi_history(:, j)
    end do

    phi_new = (1.0_dp - alpha) * phi_old + alpha * phi_extrap

    deallocate(B, rhs, coeffs, ipiv, phi_extrap)

  end subroutine diis_extrapolate


  ! ------------------------------------------------------------------
  ! Update DIIS circular buffer, maintaining temporal order (oldest=1)
  ! ------------------------------------------------------------------
  subroutine update_diis_history(phi_history, res_history, &
      & phi_input, phi_poisson, N, diis_len, iter)

    integer, intent(in) :: N, diis_len, iter
    real(kind=dp), intent(inout) :: phi_history(N, diis_len)
    real(kind=dp), intent(inout) :: res_history(N, diis_len)
    real(kind=dp), intent(in) :: phi_input(N), phi_poisson(N)

    integer :: m, j

    m = min(iter, diis_len)

    if (iter <= diis_len) then
      ! Buffer not yet full: append at next slot
      phi_history(:, iter) = phi_input
      res_history(:, iter) = phi_poisson - phi_input
    else
      ! Buffer full: shift left (discard oldest), append at end
      do j = 1, diis_len - 1
        phi_history(:, j) = phi_history(:, j + 1)
        res_history(:, j) = res_history(:, j + 1)
      end do
      phi_history(:, diis_len) = phi_input
      res_history(:, diis_len) = phi_poisson - phi_input
    end if

  end subroutine update_diis_history


  ! ------------------------------------------------------------------
  ! Fermi level bisection for charge neutrality
  ! ------------------------------------------------------------------
  function find_fermi_level(eig_kpar, eigv_kpar, kpar_grid, &
      & cfg, N, num_subbands, nk_actual, nz, dz_val, rho_doping, &
      & work_ne, work_nh) result(mu)

    real(kind=dp) :: mu
    integer, intent(in) :: N, num_subbands, nk_actual, nz
    real(kind=dp), intent(in) :: eig_kpar(num_subbands, nk_actual)
    complex(kind=dp), intent(in) :: eigv_kpar(N, num_subbands, nk_actual)
    real(kind=dp), intent(in) :: kpar_grid(nk_actual)
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), intent(in) :: dz_val, rho_doping(nz)
    real(kind=dp), intent(inout), optional :: work_ne(nz), work_nh(nz)

    real(kind=dp), allocatable :: n_elec(:), n_hole(:)

    if (present(work_ne) .and. present(work_nh)) then
      ! Use pre-allocated work arrays from caller
      call fermi_bisect(eig_kpar, eigv_kpar, kpar_grid, cfg, N, &
        & num_subbands, nk_actual, nz, dz_val, rho_doping, &
        & work_ne, work_nh, mu)
    else
      ! Allocate locally (backward compatible)
      allocate(n_elec(nz), n_hole(nz))
      call fermi_bisect(eig_kpar, eigv_kpar, kpar_grid, cfg, N, &
        & num_subbands, nk_actual, nz, dz_val, rho_doping, &
        & n_elec, n_hole, mu)
      deallocate(n_elec, n_hole)
    end if

  end function find_fermi_level


  ! ------------------------------------------------------------------
  ! Core Fermi level bisection (uses provided work arrays)
  ! ------------------------------------------------------------------
  subroutine fermi_bisect(eig_kpar, eigv_kpar, kpar_grid, &
      & cfg, N, num_subbands, nk_actual, nz, dz_val, rho_doping, &
      & n_elec, n_hole, mu)

    integer, intent(in) :: N, num_subbands, nk_actual, nz
    real(kind=dp), intent(in) :: eig_kpar(num_subbands, nk_actual)
    complex(kind=dp), intent(in) :: eigv_kpar(N, num_subbands, nk_actual)
    real(kind=dp), intent(in) :: kpar_grid(nk_actual)
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), intent(in) :: dz_val, rho_doping(nz)
    real(kind=dp), intent(inout) :: n_elec(nz), n_hole(nz)
    real(kind=dp), intent(out) :: mu

    real(kind=dp) :: mu_lo, mu_hi, mu_mid
    real(kind=dp) :: charge_excess, target_charge
    integer :: ib, iz

    target_charge = 0.0_dp
    do iz = 1, nz
      target_charge = target_charge + rho_doping(iz) * dz_val
    end do

    mu_lo = minval(eig_kpar) - FERMI_BOUND_PAD
    mu_hi = maxval(eig_kpar) + FERMI_BOUND_PAD

    do ib = 1, MAX_FERMI_BISECT
      mu_mid = 0.5_dp * (mu_lo + mu_hi)

      call compute_charge_density_qw(n_elec, n_hole, eigv_kpar, &
        & eig_kpar, kpar_grid, mu_mid, cfg%sc%temperature, &
        & nz, num_subbands, nk_actual, cfg%numcb)

      charge_excess = 0.0_dp
      do iz = 1, nz
        charge_excess = charge_excess + (n_elec(iz) - n_hole(iz)) * dz_val
      end do
      charge_excess = charge_excess - target_charge

      if (charge_excess > 0.0_dp) then
        mu_hi = mu_mid
      else
        mu_lo = mu_mid
      end if

      if (abs(mu_hi - mu_lo) < FERMI_TOL) exit
    end do

    mu = 0.5_dp * (mu_lo + mu_hi)

  end subroutine fermi_bisect


  ! ------------------------------------------------------------------
  ! Apply electrostatic potential to all 3 profile columns
  ! ------------------------------------------------------------------
  subroutine apply_potential_to_profile(profile, profile_base, phi, nz)
    integer, intent(in)          :: nz
    real(kind=dp), intent(inout) :: profile(nz, 3)
    real(kind=dp), intent(in)    :: profile_base(nz, 3)
    real(kind=dp), intent(in)    :: phi(nz)
    integer :: iz
    do iz = 1, nz
      profile(iz, 1) = profile_base(iz, 1) - phi(iz)
      profile(iz, 2) = profile_base(iz, 2) - phi(iz)
      profile(iz, 3) = profile_base(iz, 3) - phi(iz)
    end do
  end subroutine apply_potential_to_profile


  ! ------------------------------------------------------------------
  ! Apply electrostatic potential to 2D band profile (all 3 columns)
  ! ------------------------------------------------------------------
  subroutine apply_potential_to_profile_2d(profile_2d, profile_2d_base, phi, nx, ny)
    real(kind=dp), intent(inout) :: profile_2d(:,:)
    real(kind=dp), intent(in)    :: profile_2d_base(:,:)
    real(kind=dp), intent(in)    :: phi(:,:)
    integer, intent(in)          :: nx, ny
    integer :: ix, iy, ij

    do iy = 1, ny
      do ix = 1, nx
        ij = (iy - 1) * nx + ix
        profile_2d(ij, 1) = profile_2d_base(ij, 1) - phi(ix, iy)
        profile_2d(ij, 2) = profile_2d_base(ij, 2) - phi(ix, iy)
        profile_2d(ij, 3) = profile_2d_base(ij, 3) - phi(ix, iy)
      end do
    end do
  end subroutine apply_potential_to_profile_2d


  ! ------------------------------------------------------------------
  ! Self-consistent loop for quantum wire (2D confinement)
  ! ------------------------------------------------------------------
  subroutine self_consistent_loop_wire(profile_2d, cfg, kpterms_2d, grid, &
      & coo_cache, eigen_cfg, eig_wire, eigv_wire, phi_out, n_electron_out, &
      & n_hole_out, converged_out)

    real(kind=dp), allocatable, intent(inout) :: profile_2d(:,:)
    type(simulation_config), intent(in) :: cfg
    type(csr_matrix), intent(in)        :: kpterms_2d(:)
    type(spatial_grid), intent(in)      :: grid
    type(wire_coo_cache), intent(inout) :: coo_cache
    type(eigensolver_config), intent(in) :: eigen_cfg
    real(kind=dp), allocatable, intent(inout) :: eig_wire(:,:)
    complex(kind=dp), allocatable, intent(inout) :: eigv_wire(:,:,:)
    real(kind=dp), allocatable, intent(out), optional :: phi_out(:,:)
    real(kind=dp), allocatable, intent(out), optional :: n_electron_out(:,:)
    real(kind=dp), allocatable, intent(out), optional :: n_hole_out(:,:)
    logical, intent(out), optional :: converged_out

    ! SC loop variables
    integer :: iter, niter, k_idx
    integer :: num_subbands, nk_actual
    integer :: ix, iy
    real(kind=dp) :: delta_phi, fermi_level

    ! Potential arrays
    real(kind=dp), allocatable :: phi_old(:,:), phi_new(:,:), phi_poisson(:,:)
    real(kind=dp), allocatable :: phi_old_flat(:), phi_new_flat(:), phi_poisson_flat(:)
    real(kind=dp), allocatable :: rho_2d(:,:), rho_2d_flat(:)
    real(kind=dp), allocatable :: epsilon_2d(:,:)
    real(kind=dp), allocatable :: n_electron(:), n_hole(:)

    ! kx grid for SC charge integration
    real(kind=dp), allocatable :: kx_grid(:)

    ! Eigensolver temporaries
    type(csr_matrix)          :: HT_csr_sc
    type(eigensolver_result)  :: eigen_res_sc
    integer :: Ngrid, Ntot, nev_sc, i
    integer :: num_subbands_actual

    ! DIIS history (flat arrays)
    real(kind=dp), allocatable :: phi_history(:,:), res_history(:,:)

    ! Profile backup
    real(kind=dp), allocatable :: profile_2d_base(:,:)

    ! Doping charge
    real(kind=dp), allocatable :: rho_doping_2d(:,:)

    ! Work arrays for Fermi bisection
    real(kind=dp), allocatable :: fermi_ne(:), fermi_nh(:)

    ! Grid dimensions
    integer :: nx, ny
    real(kind=dp) :: dy_val, dx_val, kpar_max_val

    logical :: sc_converged

    nx = grid%nx
    ny = grid%ny
    Ngrid = nx * ny
    Ntot = 8 * Ngrid
    dy_val = grid%dy   ! AA
    dx_val = grid%dx   ! AA — x-spacing of confinement plane
    num_subbands = eigen_cfg%nev
    niter = cfg%sc%max_iterations

    ! Build kx grid (ensure odd for Simpson)
    nk_actual = cfg%sc%num_kpar
    if (mod(nk_actual, 2) == 0) nk_actual = nk_actual - 1

    ! Auto-determine kpar_max if not set
    kpar_max_val = cfg%sc%kpar_max
    if (kpar_max_val < tolerance) then
      kpar_max_val = cfg%waveVectorMax
      if (kpar_max_val < tolerance) kpar_max_val = 0.5_dp
    end if

    allocate(kx_grid(nk_actual))
    do i = 1, nk_actual
      kx_grid(i) = real(i - 1, kind=dp) * kpar_max_val &
        & / real(nk_actual - 1, kind=dp)
    end do

    ! Allocate working arrays
    allocate(phi_old(nx, ny), phi_new(nx, ny), phi_poisson(nx, ny))
    allocate(phi_old_flat(Ngrid), phi_new_flat(Ngrid), phi_poisson_flat(Ngrid))
    allocate(rho_2d(nx, ny), rho_2d_flat(Ngrid))
    allocate(epsilon_2d(nx, ny))
    allocate(n_electron(Ngrid), n_hole(Ngrid))
    allocate(profile_2d_base(Ngrid, 3))
    allocate(phi_history(Ngrid, cfg%sc%diis_history))
    allocate(res_history(Ngrid, cfg%sc%diis_history))
    allocate(fermi_ne(Ngrid), fermi_nh(Ngrid))

    ! Allocate eigensolver arrays for SC kx sweep
    nev_sc = num_subbands
    ! Reusable eigv array for charge density: (Ntot, nev_sc, nk_actual)
    ! These are allocated externally (eigv_wire) and used for output

    phi_history = 0.0_dp
    res_history = 0.0_dp

    ! Save initial profile as base
    profile_2d_base = profile_2d

    ! Build dielectric and doping arrays
    call build_epsilon_2d(epsilon_2d, grid, cfg%params)

    ! Build doping charge (returns allocated array)
    if (allocated(cfg%doping)) then
      call build_doping_charge_2d(rho_doping_2d, grid, cfg%doping)
    else
      allocate(rho_doping_2d(nx, ny))
      rho_doping_2d = 0.0_dp
    end if

    ! Initialize potential
    phi_old = 0.0_dp
    phi_new = 0.0_dp
    sc_converged = .false.

    ! Fermi level
    if (cfg%sc%fermi_mode == 1) then
      fermi_level = cfg%sc%fermi_level
    else
      fermi_level = 0.0_dp
    end if

    ! --- Main SC loop ---
    print *, '=== Wire Self-Consistent Loop Start ==='
    print *, '  max_iterations:', niter
    print *, '  tolerance:', cfg%sc%tolerance
    print *, '  mixing_alpha:', cfg%sc%mixing_alpha
    print *, '  diis_history:', cfg%sc%diis_history
    print *, '  num_kx:', nk_actual
    print *, '  kx_max:', kpar_max_val
    print *, '  temperature:', cfg%sc%temperature
    print *, '  grid: nx=', nx, ' ny=', ny, ' Ngrid=', Ngrid
    print *, ''

    do iter = 1, niter

      ! Step 1: Apply current potential to profile
      call apply_potential_to_profile_2d(profile_2d, profile_2d_base, phi_old, nx, ny)

      ! Step 2: Solve eigenproblem at each kx
      ! Build Hamiltonian at kx=0 (first kx point) — this initializes COO cache
      call ZB8bandGeneralized(HT_csr_sc, kx_grid(1), profile_2d, &
        & kpterms_2d, cfg, coo_cache)

      do k_idx = 1, nk_actual
        if (k_idx > 1) then
          ! Rebuild with new kx (COO cache reuses CSR structure)
          call ZB8bandGeneralized(HT_csr_sc, kx_grid(k_idx), profile_2d, &
            & kpterms_2d, cfg, coo_cache)
        end if

        call solve_sparse_evp(HT_csr_sc, eigen_cfg, eigen_res_sc)

        if (.not. eigen_res_sc%converged) then
          print *, '  WARNING: eigensolver did not converge at kx=', kx_grid(k_idx)
        end if

        num_subbands_actual = min(eigen_res_sc%nev_found, nev_sc)

        ! Store eigenvalues
        if (num_subbands_actual > 0) then
          do i = 1, num_subbands_actual
            eig_wire(i, k_idx) = eigen_res_sc%eigenvalues(i)
          end do
          ! Pad with large values if fewer found
          do i = num_subbands_actual + 1, nev_sc
            eig_wire(i, k_idx) = 1.0e10_dp
          end do
          ! Store eigenvectors
          eigv_wire(:, 1:num_subbands_actual, k_idx) = &
            & eigen_res_sc%eigenvectors(:, 1:num_subbands_actual)
          ! Zero out unused subbands
          if (num_subbands_actual < nev_sc) then
            eigv_wire(:, num_subbands_actual+1:nev_sc, k_idx) = cmplx(0.0_dp, 0.0_dp, kind=dp)
          end if
        end if

        call eigensolver_result_free(eigen_res_sc)
      end do

      ! Note: HT_csr_sc is kept alive across SC iterations so the COO cache
      ! can update values in-place.  It is freed after the SC loop completes.

      ! Step 3: Find Fermi level (if charge neutrality mode)
      if (cfg%sc%fermi_mode == 0) then
        fermi_level = find_fermi_level_wire(eig_wire, eigv_wire, kx_grid, &
          & cfg, Ntot, nev_sc, nk_actual, Ngrid, dy_val * dx_val, &
          & rho_doping_2d, fermi_ne, fermi_nh, nx, ny)
      end if

      ! Step 4: Compute charge density
      call compute_charge_density_wire(n_electron, n_hole, eigv_wire, &
        & eig_wire, kx_grid, fermi_level, cfg%sc%temperature, &
        & ny, nx, nev_sc, nk_actual, cfg%numcb)
      ! Note: compute_charge_density_wire takes (Ny, Nz) but for our wire
      ! grid it's (ny, nx). The flattened array is the same size.

      ! Step 5: Build total charge and solve Poisson
      ! Convert from cm^-3 to C/nm^3: rho = e * (n_h - n_e + rho_doping) * CM3_TO_PER_NM3
      do i = 1, Ngrid
        rho_2d_flat(i) = e * (n_hole(i) - n_electron(i) &
          & + rho_doping_2d(1 + mod(i-1, nx), 1 + (i-1)/nx)) * CM3_TO_PER_NM3
      end do

      ! Reshape flat rho to 2D for solve_poisson_2d
      do i = 1, Ngrid
        rho_2d(1 + mod(i-1, nx), 1 + (i-1)/nx) = rho_2d_flat(i)
      end do

      ! solve_poisson_2d uses e0 in C/(V*nm), rho in C/nm^3.
      ! dy, dz in AA are converted to nm for consistent Poisson units.
      call solve_poisson_2d(phi_poisson, rho_2d, epsilon_2d, &
        & dy_val * ANGSTROM_TO_NM, dx_val * ANGSTROM_TO_NM, &
        & nx, ny, cfg%sc%bc_left)

      ! Step 6: Flatten phi for mixing
      do i = 1, Ngrid
        phi_old_flat(i) = phi_old(1 + mod(i-1, nx), 1 + (i-1)/nx)
        phi_poisson_flat(i) = phi_poisson(1 + mod(i-1, nx), 1 + (i-1)/nx)
      end do

      ! Mix potential
      if (iter <= cfg%sc%diis_history .or. cfg%sc%diis_history < 2) then
        call linear_mix(phi_new_flat, phi_old_flat, phi_poisson_flat, &
          & Ngrid, cfg%sc%mixing_alpha)
      else
        call diis_extrapolate(phi_new_flat, phi_old_flat, phi_poisson_flat, &
          & Ngrid, phi_history, res_history, cfg%sc%diis_history, &
          & cfg%sc%mixing_alpha, iter)
      end if

      ! Update DIIS history
      call update_diis_history(phi_history, res_history, &
        & phi_old_flat, phi_poisson_flat, Ngrid, cfg%sc%diis_history, iter)

      ! Step 7: Check convergence
      delta_phi = maxval(abs(phi_new_flat - phi_old_flat))

      print '(A, I4, A, ES12.4, A, F10.4, A, ES12.4)', &
        & '  iter:', iter, '  |dPhi|:', delta_phi, &
        & '  mu:', fermi_level, '  max|rho|:', maxval(abs(rho_2d_flat))

      if (delta_phi < cfg%sc%tolerance) then
        print *, '=== Wire SC loop converged at iteration', iter, '==='
        phi_old_flat = phi_new_flat
        sc_converged = .true.
        exit
      end if

      phi_old_flat = phi_new_flat

      ! Reshape back to 2D for next iteration
      do i = 1, Ngrid
        phi_old(1 + mod(i-1, nx), 1 + (i-1)/nx) = phi_old_flat(i)
      end do

    end do

    if (.not. sc_converged) then
      print *, 'Warning: Wire SC loop did not converge after', niter, 'iterations'
      print *, '  Final |dPhi|:', delta_phi
    end if

    ! Final update: apply converged potential to profile
    do i = 1, Ngrid
      phi_old(1 + mod(i-1, nx), 1 + (i-1)/nx) = phi_old_flat(i)
    end do
    call apply_potential_to_profile_2d(profile_2d, profile_2d_base, phi_old, nx, ny)

    ! --- Copy converged diagnostics to optional output arrays ---
    if (present(phi_out)) then
      allocate(phi_out(nx, ny))
      phi_out = phi_old
    end if
    if (present(n_electron_out)) then
      allocate(n_electron_out(nx, ny))
      do iy = 1, ny
        do ix = 1, nx
          n_electron_out(ix, iy) = n_electron((iy - 1) * nx + ix)
        end do
      end do
    end if
    if (present(n_hole_out)) then
      allocate(n_hole_out(nx, ny))
      do iy = 1, ny
        do ix = 1, nx
          n_hole_out(ix, iy) = n_hole((iy - 1) * nx + ix)
        end do
      end do
    end if

    if (present(converged_out)) converged_out = sc_converged

    ! --- Cleanup ---
    if (allocated(HT_csr_sc%values)) call csr_free(HT_csr_sc)
    deallocate(kx_grid, phi_old, phi_new, phi_poisson)
    deallocate(phi_old_flat, phi_new_flat, phi_poisson_flat)
    deallocate(rho_2d, rho_2d_flat, epsilon_2d)
    deallocate(n_electron, n_hole)
    deallocate(profile_2d_base, rho_doping_2d)
    deallocate(phi_history, res_history)
    deallocate(fermi_ne, fermi_nh)

  end subroutine self_consistent_loop_wire


  ! ------------------------------------------------------------------
  ! Fermi level bisection for wire (2D confinement)
  ! ------------------------------------------------------------------
  function find_fermi_level_wire(eig_kx, eigv_kx, kx_grid, &
      & cfg, N, num_subbands, nk_actual, Ngrid, dA, rho_doping_2d, &
      & work_ne, work_nh, nx, ny) result(mu)

    real(kind=dp) :: mu
    integer, intent(in) :: N, num_subbands, nk_actual, Ngrid, nx, ny
    real(kind=dp), intent(in) :: eig_kx(num_subbands, nk_actual)
    complex(kind=dp), intent(in) :: eigv_kx(N, num_subbands, nk_actual)
    real(kind=dp), intent(in) :: kx_grid(nk_actual)
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), intent(in) :: dA
    real(kind=dp), intent(in) :: rho_doping_2d(nx, ny)
    real(kind=dp), intent(inout) :: work_ne(Ngrid), work_nh(Ngrid)

    real(kind=dp) :: mu_lo, mu_hi, mu_mid
    real(kind=dp) :: charge_excess, target_charge
    integer :: ib, p, ix, iy, ij

    ! Target charge from doping (integrated over area, converted to linear density)
    target_charge = 0.0_dp
    do iy = 1, ny
      do ix = 1, nx
        target_charge = target_charge + rho_doping_2d(ix, iy) * dA
      end do
    end do

    mu_lo = minval(eig_kx) - FERMI_BOUND_PAD
    mu_hi = maxval(eig_kx) + FERMI_BOUND_PAD

    do ib = 1, MAX_FERMI_BISECT
      mu_mid = 0.5_dp * (mu_lo + mu_hi)

      call compute_charge_density_wire(work_ne, work_nh, eigv_kx, &
        & eig_kx, kx_grid, mu_mid, cfg%sc%temperature, &
        & ny, nx, num_subbands, nk_actual, cfg%numcb)

      charge_excess = 0.0_dp
      do p = 1, Ngrid
        charge_excess = charge_excess + (work_ne(p) - work_nh(p)) * dA
      end do
      charge_excess = charge_excess - target_charge

      if (charge_excess > 0.0_dp) then
        mu_hi = mu_mid
      else
        mu_lo = mu_mid
      end if

      if (abs(mu_hi - mu_lo) < FERMI_TOL) exit
    end do

    mu = 0.5_dp * (mu_lo + mu_hi)

  end function find_fermi_level_wire


  ! ------------------------------------------------------------------
  ! Build dielectric constant array from material parameters
  ! ------------------------------------------------------------------
  subroutine build_epsilon(epsilon, cfg, nz)
    integer, intent(in) :: nz
    real(kind=dp), intent(out) :: epsilon(nz)
    type(simulation_config), intent(in) :: cfg
    integer :: iz, ilayer
    integer, allocatable :: layer_index(:)
    call map_layer_to_grid(layer_index, cfg, nz)
    epsilon = 12.90_dp
    do iz = 1, nz
      ilayer = layer_index(iz)
      if (ilayer > 0) then
        if (cfg%params(ilayer)%eps0 > 0.0_dp) then
          epsilon(iz) = cfg%params(ilayer)%eps0
        end if
      end if
    end do
    deallocate(layer_index)
  end subroutine build_epsilon


  ! ------------------------------------------------------------------
  ! Build doping charge density array (ND - NA) in cm^-3
  ! ------------------------------------------------------------------
  subroutine build_doping_charge(rho_doping, cfg, nz)
    integer, intent(in) :: nz
    real(kind=dp), intent(out) :: rho_doping(nz)
    type(simulation_config), intent(in) :: cfg
    integer :: iz, ilayer
    integer, allocatable :: layer_index(:)
    rho_doping = 0.0_dp
    if (.not. allocated(cfg%doping)) return
    call map_layer_to_grid(layer_index, cfg, nz)
    do iz = 1, nz
      ilayer = layer_index(iz)
      if (ilayer > 0) then
        rho_doping(iz) = cfg%doping(ilayer)%ND - cfg%doping(ilayer)%NA
      end if
    end do
    deallocate(layer_index)
  end subroutine build_doping_charge


  ! ------------------------------------------------------------------
  ! Map each grid point to its layer index (1-based, 0 = unmapped)
  ! ------------------------------------------------------------------
  subroutine map_layer_to_grid(layer_index, cfg, nz)
    integer, allocatable, intent(out) :: layer_index(:)
    type(simulation_config), intent(in) :: cfg
    integer, intent(in) :: nz
    integer :: iz, ilayer
    allocate(layer_index(nz))
    layer_index = 0
    do ilayer = 1, cfg%numLayers
      do iz = cfg%intStartPos(ilayer), cfg%intEndPos(ilayer)
        if (iz >= 1 .and. iz <= nz) then
          layer_index(iz) = ilayer
        end if
      end do
    end do
  end subroutine map_layer_to_grid


  ! ------------------------------------------------------------------
  ! Build 2D dielectric constant array from wire grid material_id
  ! Returns epsilon(nx, ny) reshaped from flat grid%material_id(:).
  ! ------------------------------------------------------------------
  subroutine build_epsilon_2d(epsilon, grid, params)
    real(kind=dp), allocatable, intent(out) :: epsilon(:,:)
    type(spatial_grid), intent(in) :: grid
    type(paramStruct), intent(in)  :: params(:)

    integer :: nx, ny, ix, iy, ij, mid
    real(kind=dp) :: eps_val

    nx = grid%nx
    ny = grid%ny
    allocate(epsilon(nx, ny))

    do iy = 1, ny
      do ix = 1, nx
        ij = (iy - 1) * nx + ix
        mid = grid%material_id(ij)
        if (mid > 0 .and. mid <= size(params)) then
          eps_val = params(mid)%eps0
          if (eps_val > 0.0_dp) then
            epsilon(ix, iy) = eps_val
          else
            epsilon(ix, iy) = 12.90_dp
          end if
        else
          epsilon(ix, iy) = 12.90_dp
        end if
      end do
    end do
  end subroutine build_epsilon_2d


  ! ------------------------------------------------------------------
  ! Build 2D doping charge density (ND - NA) in cm^-3 from wire grid.
  ! Returns rho_doping(nx, ny) reshaped from flat grid%material_id(:).
  ! ------------------------------------------------------------------
  subroutine build_doping_charge_2d(rho_doping, grid, doping)
    real(kind=dp), allocatable, intent(out) :: rho_doping(:,:)
    type(spatial_grid), intent(in) :: grid
    type(doping_spec), intent(in)  :: doping(:)

    integer :: nx, ny, ix, iy, ij, mid

    nx = grid%nx
    ny = grid%ny
    allocate(rho_doping(nx, ny))
    rho_doping = 0.0_dp

    do iy = 1, ny
      do ix = 1, nx
        ij = (iy - 1) * nx + ix
        mid = grid%material_id(ij)
        if (mid > 0 .and. mid <= size(doping)) then
          rho_doping(ix, iy) = doping(mid)%ND - doping(mid)%NA
        end if
      end do
    end do
  end subroutine build_doping_charge_2d

end module sc_loop
