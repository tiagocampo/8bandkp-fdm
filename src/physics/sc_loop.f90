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
  use charge_density
  use poisson
  implicit none

  private
  public :: self_consistent_loop
  public :: linear_mix
  public :: diis_extrapolate
  public :: find_fermi_level
  public :: build_epsilon
  public :: build_doping_charge

  ! LAPACK/MKL external declarations
  integer :: ilaenv
  real(kind=dp) :: dlamch
  external :: ILAENV, DLAMCH, zheevx, dgesv

contains

  ! ------------------------------------------------------------------
  ! Main self-consistent loop for QW simulations
  ! ------------------------------------------------------------------
  subroutine self_consistent_loop(profile, cfg, kpterms, HT, eig, eigv, &
      & smallk, N, il, iuu)

    real(kind=dp), allocatable, intent(inout) :: profile(:,:)
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), intent(in) :: kpterms(:,:,:)
    complex(kind=dp), allocatable, intent(inout) :: HT(:,:)
    real(kind=dp), allocatable, intent(inout) :: eig(:,:)
    complex(kind=dp), allocatable, intent(inout) :: eigv(:,:,:)
    type(wavevector), allocatable, intent(in) :: smallk(:)
    integer, intent(in) :: N, il, iuu

    ! SC loop variables
    integer :: iter, niter, info
    integer :: num_subbands, num_kpar, nk_actual
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
    real(kind=dp) :: abstol, vl, vu
    integer :: NB_val

    ! DIIS history
    real(kind=dp), allocatable :: phi_history(:,:), res_history(:,:)

    ! Profile backup
    real(kind=dp), allocatable :: profile_base(:,:)

    ! Doping charge
    real(kind=dp), allocatable :: rho_doping(:)

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
    num_kpar = cfg%sc%num_kpar
    if (mod(num_kpar, 2) == 0) num_kpar = num_kpar - 1
    nk_actual = num_kpar

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

    phi_history = 0.0_dp
    res_history = 0.0_dp

    ! Save initial profile as base
    profile_base = profile

    ! Build dielectric and doping arrays
    call build_epsilon(epsilon, cfg, nz)
    call build_doping_charge(rho_doping, cfg, nz)

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
    vl = 0.0_dp
    vu = 0.0_dp
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
    call zheevx('V', 'I', 'U', N, HT, N, vl, vu, il, iuu, abstol, M_out, &
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
      do iz = 1, nz
        profile(iz, 1) = profile_base(iz, 1) - phi_old(iz)
        profile(iz, 2) = profile_base(iz, 2) - phi_old(iz)
        profile(iz, 3) = profile_base(iz, 3) - phi_old(iz)
      end do

      ! Step 2: Solve eigenproblem at each k_par
      do k_idx = 1, nk_actual
        wv%kx = kpar_grid(k_idx)
        wv%ky = 0.0_dp
        wv%kz = 0.0_dp

        call ZB8bandQW(HT, wv, profile, kpterms)

        call zheevx('V', 'I', 'U', N, HT, N, vl, vu, il, iuu, abstol, &
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
          & cfg, N, num_subbands, nk_actual, nz, dz_val, rho_doping)
      end if

      ! Step 4: Compute charge density
      call compute_charge_density_qw(n_electron, n_hole, eigv_kpar, &
        & eig_kpar, kpar_grid, fermi_level, cfg%sc%temperature, &
        & nz, num_subbands, nk_actual, cfg%numcb, dz_val)

      ! Step 5: Build total charge and solve Poisson
      do iz = 1, nz
        rho(iz) = e * (n_hole(iz) - n_electron(iz) + rho_doping(iz)) / 1.0e21_dp
      end do

      ! Convert dz from Angstrom to nm for Poisson solver (rho in C/nm^3, e0 in C/(V*nm))
      call solve_poisson(phi_poisson, rho, epsilon, dz_val * 0.1_dp, nz, &
        & cfg%sc%bc_left, cfg%sc%bc_right, cfg%sc%bc_type)

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
    do iz = 1, nz
      profile(iz, 1) = profile_base(iz, 1) - phi_old(iz)
      profile(iz, 2) = profile_base(iz, 2) - phi_old(iz)
      profile(iz, 3) = profile_base(iz, 3) - phi_old(iz)
    end do

    ! Copy final eigenvalues at k_par=0 back to eig
    eig(:, 1) = eig_kpar(:, 1)

    ! --- Cleanup ---
    deallocate(kpar_grid, phi_old, phi_new, phi_poisson)
    deallocate(rho, epsilon, n_electron, n_hole)
    deallocate(profile_base, rho_doping)
    deallocate(phi_history, res_history)
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
  ! Update DIIS circular buffer
  ! ------------------------------------------------------------------
  subroutine update_diis_history(phi_history, res_history, &
      & phi_input, phi_poisson, N, diis_len, iter)

    integer, intent(in) :: N, diis_len, iter
    real(kind=dp), intent(inout) :: phi_history(N, diis_len)
    real(kind=dp), intent(inout) :: res_history(N, diis_len)
    real(kind=dp), intent(in) :: phi_input(N), phi_poisson(N)

    integer :: idx

    idx = mod(iter - 1, diis_len) + 1
    phi_history(:, idx) = phi_input
    res_history(:, idx) = phi_poisson - phi_input

  end subroutine update_diis_history


  ! ------------------------------------------------------------------
  ! Fermi level bisection for charge neutrality
  ! ------------------------------------------------------------------
  function find_fermi_level(eig_kpar, eigv_kpar, kpar_grid, &
      & cfg, N, num_subbands, nk_actual, nz, dz_val, rho_doping) result(mu)

    real(kind=dp) :: mu
    real(kind=dp), intent(in) :: eig_kpar(num_subbands, nk_actual)
    complex(kind=dp), intent(in) :: eigv_kpar(N, num_subbands, nk_actual)
    real(kind=dp), intent(in) :: kpar_grid(nk_actual)
    type(simulation_config), intent(in) :: cfg
    integer, intent(in) :: N, num_subbands, nk_actual, nz
    real(kind=dp), intent(in) :: dz_val, rho_doping(nz)

    real(kind=dp) :: mu_lo, mu_hi, mu_mid
    real(kind=dp), allocatable :: n_elec(:), n_hole(:)
    real(kind=dp) :: charge_excess, target_charge
    real(kind=dp), parameter :: fermi_tol = 1.0e-8_dp
    integer, parameter :: max_bisect = 60
    integer :: ib, iz

    allocate(n_elec(nz), n_hole(nz))

    target_charge = 0.0_dp
    do iz = 1, nz
      target_charge = target_charge + rho_doping(iz) * dz_val
    end do

    mu_lo = minval(eig_kpar) - 2.0_dp
    mu_hi = maxval(eig_kpar) + 2.0_dp

    do ib = 1, max_bisect
      mu_mid = 0.5_dp * (mu_lo + mu_hi)

      call compute_charge_density_qw(n_elec, n_hole, eigv_kpar, &
        & eig_kpar, kpar_grid, mu_mid, cfg%sc%temperature, &
        & nz, num_subbands, nk_actual, cfg%numcb, dz_val)

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

      if (abs(mu_hi - mu_lo) < fermi_tol) exit
    end do

    mu = 0.5_dp * (mu_lo + mu_hi)

    deallocate(n_elec, n_hole)

  end function find_fermi_level


  ! ------------------------------------------------------------------
  ! Build dielectric constant array from material parameters
  ! ------------------------------------------------------------------
  subroutine build_epsilon(epsilon, cfg, nz)
    real(kind=dp), intent(out) :: epsilon(nz)
    type(simulation_config), intent(in) :: cfg
    integer, intent(in) :: nz

    integer :: iz, ilayer

    epsilon = 12.90_dp

    do iz = 1, nz
      do ilayer = 1, cfg%numLayers
        if (iz >= cfg%intStartPos(ilayer) .and. &
          & iz <= cfg%intEndPos(ilayer)) then
          if (cfg%params(ilayer)%eps0 > 0.0_dp) then
            epsilon(iz) = cfg%params(ilayer)%eps0
          end if
          exit
        end if
      end do
    end do

  end subroutine build_epsilon


  ! ------------------------------------------------------------------
  ! Build doping charge density array (ND - NA) in cm^-3
  ! ------------------------------------------------------------------
  subroutine build_doping_charge(rho_doping, cfg, nz)
    real(kind=dp), intent(out) :: rho_doping(nz)
    type(simulation_config), intent(in) :: cfg
    integer, intent(in) :: nz

    integer :: iz, ilayer

    rho_doping = 0.0_dp

    if (.not. allocated(cfg%doping)) return

    do iz = 1, nz
      do ilayer = 1, cfg%numLayers
        if (iz >= cfg%intStartPos(ilayer) .and. &
          & iz <= cfg%intEndPos(ilayer)) then
          rho_doping(iz) = cfg%doping(ilayer)%ND - cfg%doping(ilayer)%NA
          exit
        end if
      end do
    end do

  end subroutine build_doping_charge

end module sc_loop
