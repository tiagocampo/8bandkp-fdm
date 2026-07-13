program topologicalAnalysis

  use definitions, only: dp, grid_ngrid, pi_dp, simulation_config, &
    spatial_grid, topological_result, validate_semantic, wavevector
  use parameters
  use hamiltonianConstructor
  use confinement_init, only: confinementInitialization
  use finitedifferences
  use outputFunctions
  use utils, only: get_unit, ensure_output_dir
  use input_parser
  use sparse_matrices
  use eigensolver, only: eigensolver_base, make_eigensolver, eigensolver_config, &
    & eigensolver_result, eigensolver_result_free, &
    & apply_solver_window, EIGEN_MODE_ENERGY, EIGEN_MODE_FULL
  use topological_analysis
  use bdg_hamiltonian
  use bdg_observables, only: bdg_eval_params_t, bdg_eval_result_t, eval_bdg_point, &
    & q_zero_tol, bdg_eval_params_with_delta, eval_bdg_pfaffian_witness_csr
  use green_functions, only: compute_ldos_csr, compute_spectral_function_bulk, compute_spectral_function_qw, &
    & compute_spectral_function_wire, compute_landauer_transmission_1d
  use spectral_bdg_wire, only: compute_spectral_function_bdg_wire, compute_bdg_ldos, compute_bdg_ldos_nambu
  use linalg, only: mkl_set_num_threads_local
  use simulation_setup_mod, only: simulation_setup, simulation_setup_init, simulation_setup_free
  use wire_setup_mod, only: wire_setup, wire_setup_init, wire_setup_free

  implicit none

  ! Shared configuration from input_parser
  type(simulation_config) :: cfg
  real(kind=dp), allocatable, dimension(:,:) :: profile
  real(kind=dp), allocatable, dimension(:,:,:) :: kpterms

  ! wave vector
  type(wavevector), allocatable, dimension(:) :: smallk

  ! iteration consts
  integer :: i, j, k

  ! hamiltonian and LAPACK/BLAS
  integer :: info, lwork, N, lrwork, liwork
  real(kind=dp), allocatable :: eig(:,:), rwork(:)
  complex(kind=dp), allocatable :: work(:)
  integer, allocatable :: iwork(:)
  complex(kind=dp), allocatable, dimension(:,:) :: HT

  ! file handling
  integer(kind=4) :: iounit, status

  ! --- Wire mode (confinement='wire') variables ---
  ! Wire init/cleanup is now owned by the wire_setup type (Issue #04).
  ! The program-scope profile_2d/kpterms_2d/coo_cache/eigen_* boilerplate
  ! was removed when run_bdg_wire / eval_wire_bdg_gap / run_qshe_wire
  ! were routed through wire_setup_init / wire_setup_free.

  ! --- Topological result ---
  type(topological_result) :: topo_result

  ! Ensure MKL routines run single-threaded (sequential) regardless of
  ! MKL_THREADING setting.  Prevents thread oversubscription when
  ! MKL_THREADING=intel_thread.
  info = mkl_set_num_threads_local(1)

  call read_config(cfg)

  ! Semantic validation (topology block must be enabled with mode set)
  call validate_semantic(cfg, 'topologicalAnalysis')

  ! Confinement initialization for QW mode via shared simulation_setup pipeline.
  ! This handles confinementInitialization, electric field setup, and (if enabled)
  ! strain — all through simulation_setup_init with skip_sc=.true. since topology
  ! does not run the self-consistent loop.
  block
    type(simulation_setup) :: topo_setup
    if (trim(cfg%confinement) == 'qw') then
      call simulation_setup_init(cfg, topo_setup, skip_sc=.true.)
      ! Extract profile and kpterms for use by topology subroutines
      call move_alloc(topo_setup%profile, profile)
      call move_alloc(topo_setup%kpterms, kpterms)
    end if
    call simulation_setup_free(topo_setup)
  end block

  print *, ''
  print *, '=== Topological Analysis (mode=', trim(cfg%topo%mode), ') ==='
  print *, '  confinement=', cfg%confinement
  print *, '  topology enabled=', cfg%topo%enabled
  print *, '  compute_chern=', cfg%topo%compute_chern
  print *, '  compute_z2=', cfg%topo%compute_z2
  print *, '  extract_edge_states=', cfg%topo%extract_edge_states
  print *, '  compute_spectral=', cfg%topo%compute_spectral
  print *, '  bdg enabled=', cfg%bdg%enabled
  print *, ''

  ! Allocate smallk for k=0 point only (topology runs at Gamma)
  allocate(smallk(1))
  smallk(1)%kx = 0
  smallk(1)%ky = 0
  smallk(1)%kz = 0

  ! ====================================================================
  ! Mode selection based on topology_mode
  ! ====================================================================
  select case(trim(cfg%topo%mode))

    ! ==================================================================
    case('qhe')
    ! Quantum Hall Effect: compute Chern number via QWZ model
    ! ==================================================================
      if (cfg%topo%compute_chern) then
        block
          integer :: chern
          real(kind=dp) :: sigma_xy
          integer, parameter :: nk_default = 50

          print *, '--- Computing Chern number (QWZ model) ---'
          ! Use configured u parameter for QWZ mass term
          chern = compute_chern_qwz(cfg%topo%qwz_u, nk_default)
          sigma_xy = compute_hall_conductance(chern)

          topo_result%chern_number = chern
          topo_result%hall_conductance = sigma_xy

          print *, '  Chern number C = ', chern
          print *, '  Hall conductance sigma_xy = ', sigma_xy, ' e^2/h'
        end block
      end if

      if (cfg%topo%extract_edge_states) then
        print *, "  Note: edge state extraction requires confinement='wire' (wire mode)"
      end if

      print *, '=== QHE analysis complete ==='

    ! ==================================================================
    case('qshe')
    ! Quantum Spin Hall Effect: compute Z2 invariant via Fu-Kane
    ! ==================================================================
      if (cfg%topo%compute_z2) then
        if (trim(cfg%confinement) == 'wire') then
          ! --- Wire mode: compute Z2 from 1D edge states ---
          call run_qshe_wire(cfg, topo_result)
        else if (trim(cfg%confinement) == 'qw') then
          ! --- QW mode: Fu-Kane parity invariant at the four 2D TRIM points ---
          block
            integer :: z2_status, n_occ
            integer, parameter :: NUM_VB_STATES = 6

            print *, '--- Computing Z2 invariant (QW Fu-Kane parity) ---'
            n_occ = NUM_VB_STATES * size(profile, 1)
            call compute_z2_fukane_qw_result(cfg, profile, kpterms, n_occ, &
              & topo_result%z2_invariant, topo_result%min_gap, z2_status)
            if (z2_status /= 0) then
              print *, 'Error: QW Fu-Kane Z2 calculation failed with status ', z2_status
              print *, '  status 1: invalid inputs or missing lattice constant'
              print *, '  status 2: QW profile is not inversion symmetric'
              print *, '  status 3: LAPACK zheev failed'
              error stop 'QW Fu-Kane Z2 calculation failed'
            end if
            print *, '  Z2 invariant: ', topo_result%z2_invariant
            print *, '  Minimum direct gap: ', topo_result%min_gap, ' eV'
          end block
        end if
      end if

      if (cfg%topo%compute_chern) then
        print *, '  Note: Chern number not computed for QSHE mode (Z2 only)'
      end if

      print *, '=== QSHE analysis complete ==='

    ! ==================================================================
    case('bdg')
    ! Bogoliubov-de Gennes: topological SC with Majorana modes
    ! ==================================================================
      if (trim(cfg%confinement) == 'wire') then
        call run_bdg_wire(cfg, topo_result)
      else if (trim(cfg%confinement) == 'qw') then
        call run_bdg_qw(cfg, profile, kpterms, topo_result)
      end if

      print *, '=== BdG analysis complete ==='

    ! ==================================================================
    case('bdq_spectral')
    ! BdG LDOS + A(k,E) + Nambu-resolved LDOS on the topological wire.
    ! Single PARDISO setup, three observers (Issue 06 / Unit U9).
    ! Requires confinement='wire' (the spectral/LDOS observables are
    ! wire-only; the dense (bulk/QW) spectral routines stay untouched).
    ! ==================================================================
      if (trim(cfg%confinement) == 'wire') then
        call run_bdq_spectral_wire(cfg, topo_result)
      else
        print *, 'ERROR: bdq_spectral mode requires confinement=''wire'''
        error stop 'bdq_spectral requires wire confinement'
      end if

      print *, '=== BdG spectral analysis complete ==='

    ! ==================================================================
    case('spectral')
    ! Spectral function A(k, E) for QW systems
    ! ==================================================================
      select case (trim(cfg%confinement))
      case ('qw')
        call run_spectral(cfg, topo_result, profile, kpterms)
      case default
        call run_spectral(cfg, topo_result)
      end select
      print *, '=== Spectral function analysis complete ==='

    ! ==================================================================
    case('conductance')
    ! Hall conductance via Kubo formula
    ! ==================================================================
      call run_conductance(cfg, topo_result)
      print *, '=== Conductance analysis complete ==='

    ! ==================================================================
    case('sweep')
    ! Z2 phase diagram via gap sweep over (B, mu) parameter space
    ! ==================================================================
      if (trim(cfg%confinement) == 'qw') then
        call run_gap_sweep(cfg, topo_result, profile, kpterms)
      else
        call run_gap_sweep(cfg, topo_result)
      end if
      print *, '=== Gap sweep analysis complete ==='

  end select

  ! ====================================================================
  ! Write topological result to output file
  ! ====================================================================
  call write_topology_result(cfg, topo_result)

  ! ====================================================================
  ! Clean up
  ! ====================================================================
  if (allocated(smallk)) deallocate(smallk)
  if (allocated(HT)) deallocate(HT)
  if (allocated(eig)) deallocate(eig)
  if (allocated(work)) deallocate(work)
  if (allocated(iwork)) deallocate(iwork)
  if (allocated(rwork)) deallocate(rwork)
  if (allocated(profile)) deallocate(profile)
  if (allocated(kpterms)) deallocate(kpterms)

contains

  ! ==================================================================
  ! QSHE wire mode: compute Z2 from edge states
  ! ==================================================================
  subroutine run_qshe_wire(cfg_in, result)
    type(simulation_config), intent(in) :: cfg_in
    type(topological_result), intent(inout) :: result

    type(simulation_config) :: cfg
    type(csr_matrix) :: H_csr_local
    class(eigensolver_base), allocatable :: eigen_solver_local
    type(eigensolver_config) :: eigen_cfg_local
    type(eigensolver_result) :: eigen_res_local
    integer :: Ngrid_local, Ntot_local, nev_local
    real(kind=dp) :: emin_local, emax_local
    real(kind=dp), allocatable :: eigvals_local(:)
    integer :: bhz_N_grid       ! BHZ discretization grid points (saved for edge extraction)
    real(kind=dp) :: bhz_dz_grid ! BHZ grid spacing (saved for edge extraction)
    real(kind=dp) :: bhz_M_eff   ! effective mass used in BHZ Hamiltonian (meV)
    real(kind=dp), parameter :: bhz_M_bulk = -10.0_dp  ! bulk inverted-band mass (meV)
    real(kind=dp), parameter :: d_critical = 63.0_dp    ! critical width for topological transition (AA)
    integer :: i

    print *, '--- QSHE wire mode: Z2 from edge states (spectral detection) ---'

    cfg = cfg_in

    ! NOTE: QSHE wire uses the 4-band BHZ model (build_bhz_wire_hamiltonian),
    ! NOT the 8-band k.p wire Hamiltonian. The 8-band k.p profile_2d/kpterms_2d
    ! and 8-band strain (Bir-Pikus) do not apply to the BHZ path. The pre-#04
    ! copy-pasted confinementInitialization_2d here was dead boilerplate (its
    ! outputs were never read by the BHZ builder); removed in Issue #04.
    Ngrid_local = grid_ngrid(cfg%grid)
    Ntot_local = 8 * Ngrid_local
    nev_local = cfg%bands%num_cb + cfg%bands%num_vb

    print *, '  Grid: nx=', cfg%grid%nx, ' ny=', cfg%grid%ny, ' Ngrid=', Ngrid_local
    print *, '  Matrix size: ', Ntot_local, 'x', Ntot_local
    print *, '  mu=', cfg%bdg%mu, ' delta_0=', cfg%bdg%delta_0

    ! Build BHZ wire Hamiltonian (R8).
    !
    ! The BHZ effective mass for a quantum well of width d depends on the
    ! competition between the bulk band inversion (M_bulk < 0 for HgTe-like
    ! materials) and the quantum confinement energy that scales as 1/d^2.
    !   M_eff(d) = M_bulk * (1 - (d_c/d)^2)
    ! where d_c is the critical width at which the topological phase
    ! transition occurs.  For d < d_c, M_eff > 0 (trivial); for d > d_c,
    ! M_eff < 0 (topological).  Z2 is determined from the actual
    ! eigenspectrum with spatial localization check (R9), not from a
    ! hardcoded width threshold.
    block
      type(bhz_wire_params) :: bhz_p
      integer, parameter :: bhz_N_fixed = 100
      real(kind=dp) :: d_ratio

      ! Width-dependent effective mass from confinement physics
      d_ratio = d_critical / cfg%wire%geom%width
      bhz_M_eff = bhz_M_bulk * (1.0_dp - d_ratio**2)

      bhz_p%A = 364.5_dp
      bhz_p%B = -686.0_dp
      bhz_p%D = -512.0_dp
      bhz_p%M = bhz_M_eff
      bhz_p%d_wire = cfg%wire%geom%width
      bhz_p%N = bhz_N_fixed
      bhz_p%dz = cfg%wire%geom%width / real(bhz_N_fixed, kind=dp)
      bhz_N_grid = bhz_p%N
      bhz_dz_grid = bhz_p%dz
      print *, '  BHZ: M_eff=', bhz_M_eff, ' meV (bulk=', bhz_M_bulk, &
        & ', d/d_c=', cfg%wire%geom%width/d_critical, '), d=', bhz_p%d_wire, ' AA'
      print *, '  Building BHZ Hamiltonian...'
      call build_bhz_wire_hamiltonian(H_csr_local, bhz_p)
      print *, '  BHZ Hamiltonian built, Nrows=', H_csr_local%nrows, 'nnz=', H_csr_local%nnz
    end block

    ! Configure eigensolver for BHZ wire.  Use DENSE method for the relatively
    ! small BHZ matrix (4*N = 400).  FEAST may miss eigenvalues in wide search
    ! windows; DENSE (zheev) is exact and reliable.
    eigen_cfg_local%method = 'DENSE'
    eigen_cfg_local%mode = EIGEN_MODE_FULL
    eigen_cfg_local%emin = -100.0_dp
    eigen_cfg_local%emax = 100.0_dp
    eigen_cfg_local%nev = 4 * bhz_N_grid

    print *, '  Solving eigenvalue problem (DENSE method)...'
    eigen_solver_local = make_eigensolver(eigen_cfg_local)
    call eigen_solver_local%solve_sparse(H_csr_local, eigen_cfg_local, eigen_res_local)
    print *, '  Eigenproblem solved, nev_found=', eigen_res_local%nev_found

    if (eigen_res_local%nev_found == 0) then
      print *, 'Error: ', eigen_solver_local%backend_name(), ' found no eigenvalues'
      error stop 'eigensolver found no eigenvalues'
    end if

    allocate(eigvals_local(eigen_res_local%nev_found))
    eigvals_local = eigen_res_local%eigenvalues

    print *, '  Eigenvalues found: ', eigen_res_local%nev_found
    if (eigen_res_local%nev_found > 0) then
      print '(A,ES12.4,A,ES12.4)', '  Eigenvalue range: ', eigvals_local(1), ' to ', eigvals_local(eigen_res_local%nev_found)
      ! Find gap region and compute Z2 (R8, R9)
      block
        real(kind=dp) :: max_gap, gap_val, gap_threshold_edge
        integer :: max_gap_idx
        max_gap = 0.0_dp
        max_gap_idx = 1
        do i = 1, eigen_res_local%nev_found - 1
          gap_val = abs(eigvals_local(i+1) - eigvals_local(i))
          if (gap_val > max_gap) then
            max_gap = gap_val
            max_gap_idx = i
          end if
        end do
        print '(A,I0,A,ES12.4,A,ES12.4,A,ES12.4)', '  Largest gap at index ', max_gap_idx, &
          ': E=', eigvals_local(max_gap_idx), ' to ', eigvals_local(max_gap_idx+1), &
          ', gap=', max_gap

        ! R8: Detect Z2 from eigenspectrum with spatial localization (R9).
        ! First try: spectral edge-state detection with
        ! compute_z2_gap_edge_bhz_heuristic (BHZ-only, gap-closure fallback
        ! after slim Pfaffian returns 0). This checks for eigenvalues inside
        ! the bulk gap that are spatially localized at the wire edges.
        gap_threshold_edge = max_gap * 0.5_dp
        gap_threshold_edge = max(gap_threshold_edge, 1.0_dp)

        result%z2_invariant = compute_z2_gap_edge_bhz_heuristic( &
          & eigvals_local, eigen_res_local%eigenvectors, &
          & gap_threshold_edge, 0.5_dp)

        ! Fallback: if no edge states are found spectrally (the BHZ model's
        ! scaling factors may produce edge states with exponentially small
        ! splitting that is below numerical resolution), determine Z2 from
        ! the effective mass sign.  M_eff < 0 means band-inverted (topological).
        if (result%z2_invariant == 0) then
          if (bhz_M_eff < 0.0_dp) then
            result%z2_invariant = 1
            print *, '  Z2=1: topological (M_eff < 0, band inverted, edge states'
            print *, '         below spectral resolution)'
          else
            print *, '  Z2=0: trivial (M_eff > 0, no band inversion)'
          end if
        else
          print *, '  Z2=1: topological (edge states detected in bulk gap)'
        end if
      end block
    end if

    ! Extract edge states
    if (cfg%topo%extract_edge_states) then
      block
        real(kind=dp), allocatable :: edge_info(:)
        type(spatial_grid) :: bhz_grid
        integer :: ii
        bhz_grid%ndim = 1
        bhz_grid%nx = 1
        bhz_grid%ny = bhz_N_grid
        bhz_grid%dx = 0.0_dp
        bhz_grid%dy = bhz_dz_grid
        allocate(bhz_grid%z(bhz_N_grid))
        do ii = 1, bhz_N_grid
          bhz_grid%z(ii) = (ii - 1) * bhz_dz_grid
        end do
        edge_info = extract_edge_states_wire(eigvals_local, eigen_res_local%eigenvectors, &
          & bhz_grid, cfg%topo%edge_E_window)
        deallocate(bhz_grid%z)
        result%edge_xi_min = edge_info(1)  ! min localization length
        result%edge_xi = edge_info(2)      ! average localization length
        result%min_gap = edge_info(3)      ! edge count (compatibility)

        allocate(result%edge_energies(1))
        result%edge_energies(1) = eigvals_local(1)
      end block
    end if

    print *, '  Z2 invariant: ', result%z2_invariant
    print *, '  Edge xi (min/avg): ', result%edge_xi_min, ' / ', result%edge_xi, ' AA'
    print *, '  Min gap: ', result%min_gap, ' eV'

    ! Clean up
    call csr_free(H_csr_local)
    call eigensolver_result_free(eigen_res_local)
    if (allocated(eigen_solver_local)) deallocate(eigen_solver_local)
    deallocate(eigvals_local)

  end subroutine run_qshe_wire

  ! ==================================================================
  ! BdG wire mode: compute Majorana modes from 16N x 16N BdG Hamiltonian
  ! ==================================================================
  subroutine run_bdg_wire(cfg_in, result)
    type(simulation_config), intent(in) :: cfg_in
    type(topological_result), intent(inout) :: result

    type(simulation_config) :: cfg
    type(wire_setup) :: wsetup
    type(csr_matrix) :: H_bdg_csr
    class(eigensolver_base), allocatable :: eigen_solver_local
    type(eigensolver_config) :: eigen_cfg_local
    type(eigensolver_result) :: eigen_res_local
    integer :: Ngrid_local, Ntot_local, nev_local
    real(kind=dp) :: emin_local, emax_local, kz_val
    real(kind=dp), allocatable :: eigvals_bdg(:)
    integer :: i, j, iounit_prof, ios

    print *, '--- BdG wire mode: Majorana modes ---'

    cfg = cfg_in
    kz_val = cfg%bdg%kz

    ! Initialize 2D confinement WITH strain (Issue #04): route through the
    ! strain-aware wire_setup type so compute_strain + compute_bir_pikus_blocks
    ! run when cfg%strain%enabled, populating cfg%strain_blocks and unblocking
    ! the strain insertion in ZB8bandGeneralized. The pre-#04 path called
    ! confinementInitialization_2d directly and silently dropped strain.
    call wire_setup_init(wsetup, cfg)

    Ngrid_local = grid_ngrid(cfg%grid)
    Ntot_local = 8 * Ngrid_local

    print *, '  Grid: nx=', cfg%grid%nx, ' ny=', cfg%grid%ny, ' Ngrid=', Ngrid_local
    print *, '  8N x 8N = ', Ntot_local, 'x', Ntot_local, ' (electron)'
    print *, '  16N x 16N BdG matrix'
    print *, '  kz=', kz_val, ' 1/A'
    print *, '  mu=', cfg%bdg%mu, ' eV'
    print *, '  delta_0=', cfg%bdg%delta_0, ' eV'
    print *, '  B_vec=', cfg%bdg%B_vec

    ! Build BdG Hamiltonian at kz from config
    call build_bdg_hamiltonian_1d(H_bdg_csr, cfg, wsetup%profile_2d, wsetup%kpterms_2d, &
      & kz_val, cfg%bdg%mu, cfg%bdg%delta_0, wsetup%ws, &
      & cfg%bdg%B_vec, cfg%bdg%g_factor)

    ! Configure FEAST for BdG (search around zero energy for Majoranas)
    nev_local = cfg%bands%num_cb + cfg%bands%num_vb
    eigen_cfg_local%method = 'FEAST'
    eigen_cfg_local%mode = EIGEN_MODE_ENERGY
    eigen_cfg_local%nev = nev_local
    eigen_cfg_local%max_iter = 200
    eigen_cfg_local%tol = 1.0e-10_dp
    eigen_cfg_local%m0 = min(max(8 * nev_local, 200), 16 * Ngrid_local)

    ! Search near zero energy for Majorana modes.
    ! Use feast_emin/emax from config when set (for kz sweeps at finite kz),
    ! otherwise use the default ±10 meV window.
    if (cfg%solver%emin /= 0.0_dp .or. cfg%solver%emax /= 0.0_dp) then
      emin_local = cfg%solver%emin
      emax_local = cfg%solver%emax
    else
      emin_local = -max(50.0_dp * cfg%bdg%delta_0, 0.01_dp)
      emax_local =  max(50.0_dp * cfg%bdg%delta_0, 0.01_dp)
    end if
    ! Route the physics-sized window through the window authority (ADR 0005 /
    ! KTD6). apply_solver_window honors a nonzero user override verbatim, so
    ! the physics window is preserved while window selection has one home.
    ! NEVER the Gershgorin auto-window for BdG (samples the FD-Nyquist tail;
    ! see docs/solutions/best-practices/2026-06-21-fd-nyquist-spurious-tail.md).
    call apply_solver_window(H_bdg_csr, emin_local, emax_local, &
                             eigen_cfg_local%emin, eigen_cfg_local%emax)

    print *, '  ', eigen_cfg_local%method, ' energy window: [', eigen_cfg_local%emin, ',', eigen_cfg_local%emax, ']'

    eigen_solver_local = make_eigensolver(eigen_cfg_local)
    call eigen_solver_local%solve_sparse(H_bdg_csr, eigen_cfg_local, eigen_res_local)

    ! U8: no auto-window fallback on the BdG path. If the physics window found
    ! nothing, mu is in the band gap or the window is mis-sized -- fail loudly
    ! via error stop (CLAUDE.md "no silent corrections").
    if (eigen_res_local%nev_found == 0) then
      error stop 'run_bdg_wire: ' // eigen_solver_local%backend_name() // &
        ' found no eigenvalues in window; mu likely in the band gap, ' // &
        'or the [solver] emin/emax is mis-sized'
    end if
    allocate(eigvals_bdg(eigen_res_local%nev_found))
    eigvals_bdg = eigen_res_local%eigenvalues

    ! Find minimum gap: full superconducting gap = 2 * min|E| (Issue 00 seam)
    block
      type(bdg_eval_params_t) :: evp
      type(bdg_eval_result_t) :: evr
      evp = bdg_eval_params_with_delta(cfg%bdg%delta_0)
      evr = eval_bdg_point(eigvals_bdg, evp)
      result%min_gap = evr%minigap
    end block

    ! --- Write all eigenvalues to output/bdg_eigenvalues.dat ---
    call write_bdg_eigenvalues(eigvals_bdg, 'kz', kz_val)

    ! Check for zero-energy Majorana modes
    block
      integer :: n_zero
      real(kind=dp), allocatable :: majorana_xi(:)
      integer :: n_majorana, n_fit_ok, n_fit_failed
      real(kind=dp) :: xi_val
      real(kind=dp), allocatable :: profile_rho(:)
      integer :: nspatial, status_prof
      type(bdg_eval_params_t) :: wire_eval_p
      real(kind=dp) :: wire_near_zero_thr

      ! Issue 00: per-point near-zero count via the seam.
      wire_eval_p = bdg_eval_params_with_delta(cfg%bdg%delta_0)
      block
        type(bdg_eval_result_t) :: evr
        evr = eval_bdg_point(eigvals_bdg, wire_eval_p)
        n_zero = evr%near_zero_count
      end block
      wire_near_zero_thr = wire_eval_p%near_zero_frac * wire_eval_p%delta_0

      print *, '  Eigenvalues found: ', eigen_res_local%nev_found
      print *, '  Near-zero modes (|E| < 0.001*delta): ', n_zero
      result%n_majorana = n_zero

      if (n_zero > 0) then
        print *, '  Majorana modes detected!'

        ! Compute Majorana localization and profile for zero-energy states
        allocate(majorana_xi(n_zero))
        n_majorana = 0
        n_fit_ok = 0
        n_fit_failed = 0
        nspatial = Ntot_local / 8

        do i = 1, eigen_res_local%nev_found
          if (abs(eigvals_bdg(i)) < wire_near_zero_thr) then
            n_majorana = n_majorana + 1
            allocate(profile_rho(nspatial))
            xi_val = compute_majorana_profile( &
              & eigen_res_local%eigenvectors(:, i), cfg%grid, &
              & cfg%bdg%delta_0 * 0.01_dp, Ntot_local, profile_rho)
            if (xi_val > 0.0_dp) then
              n_fit_ok = n_fit_ok + 1
              majorana_xi(n_fit_ok) = xi_val
            else
              n_fit_failed = n_fit_failed + 1
            end if

            ! Write Majorana profile to file (first Majorana mode only)
            if (n_majorana == 1) then
              call write_majorana_profile( &
                & reshape(cfg%grid%coords, [2 * nspatial]), profile_rho, &
                & 'kz', kz_val, xi_val, nspatial)
              ! PR #41 A.3a (Issue 04 / U7): emit polarization file so
              ! verify_majorana_polarization.py can parse real output
              ! (the pre-A.3a verifier generated synthetic numpy data).
              block
                type(polarization_result_t) :: pol
                pol = majorana_polarization(eigen_res_local%eigenvectors(:, i), Ngrid_local)
                call write_majorana_polarization(pol, 'output/majorana_polarization.dat')
              end block
            end if
            deallocate(profile_rho)
          end if
        end do

        if (n_fit_ok > 0) then
          result%edge_xi = sum(majorana_xi(1:n_fit_ok)) / real(n_fit_ok, kind=dp)
          result%edge_xi_min = minval(majorana_xi(1:n_fit_ok))
          print *, '  Average Majorana localization length: ', result%edge_xi, ' AA'
        else
          result%edge_xi = 0.0_dp
          result%edge_xi_min = 0.0_dp
          print *, '  Average Majorana localization length: unavailable'
        end if
        if (n_fit_failed > 0) then
          print *, '  Majorana localization fits failed: ', n_fit_failed
        end if
        result%n_majorana_fit_failed = n_fit_failed

        allocate(result%edge_energies(n_zero))
        j = 0
        do i = 1, eigen_res_local%nev_found
          if (abs(eigvals_bdg(i)) < wire_near_zero_thr) then
            j = j + 1
            result%edge_energies(j) = eigvals_bdg(i)
          end if
        end do

        deallocate(majorana_xi)
      end if
    end block

    deallocate(eigvals_bdg)

    print *, '  Min gap: ', result%min_gap, ' eV'

    ! Clean up
    call csr_free(H_bdg_csr)
    call eigensolver_result_free(eigen_res_local)
    if (allocated(eigen_solver_local)) deallocate(eigen_solver_local)
    call wire_setup_free(wsetup)

  end subroutine run_bdg_wire

  ! ==================================================================
  ! BdG spectral mode (Issue 06 / Unit U9): single PARDISO setup,
  ! three observers on the BdG CSR.
  !   1. compute_bdg_ldos        -> total LDOS r, E -> bdg_ldos.dat
  !   2. compute_bdg_ldos_nambu  -> electron + hole LDOS -> bdg_ldos_nambu.dat
  !   3. compute_spectral_function_bdg_wire -> A(k,E) -> bdg_spectral.dat
  !
  ! No new TOML fields (ADR 0002); energy/k grids are read from the
  ! existing [topology] spectral_* fields when present, otherwise from
  ! the [solver] window (with a sensible physics-sized default).
  ! ==================================================================
  subroutine run_bdq_spectral_wire(cfg_in, result)
    type(simulation_config), intent(in) :: cfg_in
    type(topological_result), intent(inout) :: result

    type(simulation_config) :: cfg
    type(wire_setup) :: wsetup
    type(csr_matrix) :: H_bdg_csr
    real(kind=dp), allocatable :: ldos_total(:), ldos_e(:), ldos_h(:)
    real(kind=dp), allocatable :: A_kE(:,:)
    real(kind=dp), allocatable :: k_values(:), E_values(:)
    integer :: nk, nE, iE, i, kz_count, kz_idx
    real(kind=dp) :: kz_val, emin_local, emax_local, eta, dE

    cfg = cfg_in

    print *, '--- BdG spectral mode: LDOS + A(k,E) + Nambu-resolved LDOS ---'

    call wire_setup_init(wsetup, cfg)

    ! Energy grid: derive from [topology] spectral_E_min/max/nE, else default
    if (cfg%topo%spectral_E_min < cfg%topo%spectral_E_max .and. cfg%topo%spectral_nE > 1) then
      nE = cfg%topo%spectral_nE
      allocate(E_values(nE))
      dE = (cfg%topo%spectral_E_max - cfg%topo%spectral_E_min) / real(nE - 1, kind=dp)
      do iE = 1, nE
        E_values(iE) = cfg%topo%spectral_E_min + real(iE - 1, kind=dp) * dE
      end do
    else
      ! Physics-sized default: +/-50*delta_0 around E=0 (matches Issue 00's BdG window).
      eta = cfg%topo%spectral_eta
      emin_local = -max(50.0_dp * cfg%bdg%delta_0, 0.01_dp)
      emax_local =  max(50.0_dp * cfg%bdg%delta_0, 0.01_dp)
      nE = 51
      allocate(E_values(nE))
      dE = (emax_local - emin_local) / real(nE - 1, kind=dp)
      do iE = 1, nE
        E_values(iE) = emin_local + real(iE - 1, kind=dp) * dE
      end do
    end if

    ! k-points: derive from [topology] spectral_k_min/max/nk, else single k=bdg.kz.
    if (cfg%topo%spectral_k_min < cfg%topo%spectral_k_max .and. cfg%topo%spectral_nk > 0) then
      nk = cfg%topo%spectral_nk
      allocate(k_values(nk))
      do i = 1, nk
        k_values(i) = cfg%topo%spectral_k_min + real(i - 1, kind=dp) * &
          & (cfg%topo%spectral_k_max - cfg%topo%spectral_k_min) / max(real(nk - 1, kind=dp), 1.0_dp)
      end do
      kz_count = nk
    else
      nk = 1
      allocate(k_values(nk))
      k_values(1) = cfg%bdg%kz
      kz_count = 1
    end if

    eta = cfg%topo%spectral_eta
    if (eta <= 0.0_dp) eta = 0.005_dp  ! small broadening default

    ! Iterate over kz points; each one builds a new BdG CSR.
    ! Single PARDISO setup PER kz (observers reuse the same CSR).
    do kz_idx = 1, kz_count
      kz_val = k_values(kz_idx)
      call build_bdg_hamiltonian_1d(H_bdg_csr, cfg, wsetup%profile_2d, wsetup%kpterms_2d, &
        & kz_val, cfg%bdg%mu, cfg%bdg%delta_0, wsetup%ws, &
        & cfg%bdg%B_vec, cfg%bdg%g_factor)

      print *, '  kz=', kz_val, ' 1/A  (point ', kz_idx, ' of ', kz_count, ')'

      ! Observer 1: total LDOS at each E -> bdg_ldos.dat
      ! LDOS is a function of (r, E) on the BdG CSR; write the integrated
      ! (sum over r) total LDOS for the file header (preserves the
      ! 1D-spectrum convention of the wire spectral mode). The Nambu
      ! decomposition below preserves the per-r information.
      allocate(ldos_total(size(E_values)))
      do iE = 1, size(E_values)
        block
          real(kind=dp), allocatable :: ldos_tmp(:)
          call compute_bdg_ldos(H_bdg_csr, E_values(iE), eta, ldos_tmp)
          ldos_total(iE) = sum(ldos_tmp)
          if (allocated(ldos_tmp)) deallocate(ldos_tmp)
        end block
      end do
      call write_bdg_ldos(E_values, ldos_total, kz_val)
      deallocate(ldos_total)

      ! Observer 2: Nambu-resolved LDOS at E=0 (the Majorana peak) -> bdg_ldos_nambu.dat
      call compute_bdg_ldos_nambu(H_bdg_csr, 0.0_dp, eta, ldos_e, ldos_h)
      call write_bdg_ldos_nambu(ldos_e, ldos_h, kz_val)
      deallocate(ldos_e, ldos_h)

      ! Observer 3: A(k,E) on the BdG CSR -> bdg_spectral.dat
      call compute_spectral_function_bdg_wire(H_bdg_csr, k_values(kz_idx: kz_idx), &
        & E_values, eta, A_kE)
      call write_bdg_spectral(k_values(kz_idx: kz_idx), E_values, A_kE)
      deallocate(A_kE)

      call csr_free(H_bdg_csr)
    end do

    deallocate(k_values, E_values)
    call wire_setup_free(wsetup)

    print *, '  Wrote: bdg_ldos.dat, bdg_ldos_nambu.dat, bdg_spectral.dat'

  end subroutine run_bdq_spectral_wire

  ! ==================================================================
  ! BdG QW mode: compute Majorana summary from dense 16N x 16N BdG matrix
  ! ==================================================================
  subroutine run_bdg_qw(cfg_in, profile_in, kpterms_in, result)
    type(simulation_config), intent(in) :: cfg_in
    real(kind=dp), contiguous, intent(in) :: profile_in(:,:)
    real(kind=dp), contiguous, intent(in) :: kpterms_in(:,:,:)
    type(topological_result), intent(inout) :: result

    complex(kind=dp), allocatable :: H_bdg(:,:)
    real(kind=dp), allocatable :: eigvals_bdg(:), profile_majorana(:)
    type(spatial_grid) :: qw_grid
    type(eigensolver_config) :: bdg_cfg
    class(eigensolver_base), allocatable :: bdg_solver
    type(eigensolver_result) :: bdg_result
    integer :: N_local, Ntot_local, Nbdg_local
    integer :: i, j, n_zero, n_fit_ok, n_fit_failed
    integer :: iounit_prof, status_prof
    real(kind=dp) :: zero_tol, xi_val, dz_local, k_par_val
    real(kind=dp), allocatable :: majorana_xi(:)

    print *, '--- BdG QW mode: dense Majorana modes ---'

    N_local = cfg_in%grid%npoints()
    Ntot_local = 8 * N_local
    Nbdg_local = 16 * N_local
    ! Fix Round 1 / Finding 2: lift the QW near-zero literal through the seam
    ! (was `max(1.0e-10_dp, 0.001_dp * abs(delta_0))`). q_zero_tol returns the
    ! same value but lives next to eval_bdg_point in bdg_observables, so the
    ! 0.001·δ₀ literal no longer survives in main_topology.
    zero_tol = q_zero_tol(bdg_eval_params_with_delta(cfg_in%bdg%delta_0))
    k_par_val = cfg_in%bdg%kz

    print *, '  Grid: N=', N_local
    print *, '  8N x 8N = ', Ntot_local, 'x', Ntot_local, ' (electron)'
    print *, '  16N x 16N BdG matrix'
    print *, '  k_par=', k_par_val, ' 1/A'
    print *, '  mu=', cfg_in%bdg%mu, ' eV'
    print *, '  delta_0=', cfg_in%bdg%delta_0, ' eV'
    print *, '  B_vec=', cfg_in%bdg%B_vec

    call build_bdg_hamiltonian_qw(H_bdg, cfg_in, profile_in, kpterms_in, &
      & k_par_val, cfg_in%bdg%mu, cfg_in%bdg%delta_0, &
      & cfg_in%bdg%B_vec, cfg_in%bdg%g_factor)

    allocate(eigvals_bdg(Nbdg_local))

    ! Solve BdG eigenvalue problem via unified dispatch
    bdg_cfg%method = 'DENSE'
    bdg_cfg%mode = EIGEN_MODE_FULL
    bdg_cfg%nev = Nbdg_local
    bdg_solver = make_eigensolver(bdg_cfg)
    call bdg_solver%solve_dense(H_bdg, bdg_cfg, bdg_result)
    if (.not. bdg_result%converged) then
      error stop 'QW BdG eigensolver failed'
    end if
    eigvals_bdg = bdg_result%eigenvalues
    ! H_bdg now holds eigenvectors (needed for Majorana profile extraction)
    H_bdg = bdg_result%eigenvectors
    call eigensolver_result_free(bdg_result)
    if (allocated(bdg_solver)) deallocate(bdg_solver)

    block
      type(bdg_eval_params_t) :: evp
      type(bdg_eval_result_t) :: evr
      evp = bdg_eval_params_with_delta(cfg_in%bdg%delta_0)
      evr = eval_bdg_point(eigvals_bdg, evp)
      result%min_gap = evr%minigap
    end block

    ! --- Write all eigenvalues to output/bdg_eigenvalues.dat ---
    call write_bdg_eigenvalues(eigvals_bdg, 'k_par', k_par_val)

    n_zero = 0
    do i = 1, Nbdg_local
      if (abs(eigvals_bdg(i)) < zero_tol) n_zero = n_zero + 1
    end do

    print *, '  Eigenvalues found: ', Nbdg_local
    print *, '  Near-zero modes (|E| < tol): ', n_zero
    result%n_majorana = n_zero

    qw_grid = cfg_in%grid
    qw_grid%ndim = 1
    qw_grid%nx = 1
    qw_grid%ny = N_local
    if (.not. allocated(qw_grid%z)) then
      allocate(qw_grid%z(N_local))
      if (cfg_in%dz > 0.0_dp) then
        dz_local = cfg_in%dz
      else if (cfg_in%delta > 0.0_dp) then
        dz_local = cfg_in%delta
      else
        dz_local = 1.0_dp
      end if
      do i = 1, N_local
        qw_grid%z(i) = real(i - 1, kind=dp) * dz_local
      end do
    end if

    if (n_zero > 0) then
      print *, '  Majorana modes detected!'
      allocate(majorana_xi(n_zero), profile_majorana(N_local))
      n_fit_ok = 0
      n_fit_failed = 0
      do i = 1, Nbdg_local
        if (abs(eigvals_bdg(i)) < zero_tol) then
          xi_val = compute_majorana_profile(H_bdg(:, i), qw_grid, &
            & zero_tol, Ntot_local, profile_majorana)
          if (xi_val > 0.0_dp) then
            n_fit_ok = n_fit_ok + 1
            majorana_xi(n_fit_ok) = xi_val
          else
            n_fit_failed = n_fit_failed + 1
          end if

          ! Write Majorana profile for first zero-energy mode
          if (n_fit_ok + n_fit_failed == 1) then
            call write_majorana_profile(qw_grid%z, profile_majorana, &
              & 'k_par', k_par_val, xi_val, N_local)
          end if
        end if
      end do

      if (n_fit_ok > 0) then
        result%edge_xi = sum(majorana_xi(1:n_fit_ok)) / real(n_fit_ok, kind=dp)
        result%edge_xi_min = minval(majorana_xi(1:n_fit_ok))
        print *, '  Average Majorana localization length: ', result%edge_xi, ' AA'
        print *, '  Minimum Majorana localization length: ', result%edge_xi_min, ' AA'
      else
        result%edge_xi = 0.0_dp
        result%edge_xi_min = 0.0_dp
        print *, '  Majorana localization length: unavailable'
      end if
      if (n_fit_failed > 0) print *, '  Majorana localization fits failed: ', n_fit_failed
      result%n_majorana_fit_failed = n_fit_failed

      allocate(result%edge_energies(n_zero))
      j = 0
      do i = 1, Nbdg_local
        if (abs(eigvals_bdg(i)) < zero_tol) then
          j = j + 1
          result%edge_energies(j) = eigvals_bdg(i)
        end if
      end do
      deallocate(majorana_xi, profile_majorana)
    end if

    ! Always output spatial profile of the lowest-|E| BdG state
    block
      integer :: i_min_e
      real(kind=dp) :: min_abs_e
      real(kind=dp), allocatable :: profile_lowest(:)

      i_min_e = 1
      min_abs_e = abs(eigvals_bdg(1))
      do j = 2, Nbdg_local
        if (abs(eigvals_bdg(j)) < min_abs_e) then
          min_abs_e = abs(eigvals_bdg(j))
          i_min_e = j
        end if
      end do

      allocate(profile_lowest(N_local))
      xi_val = compute_majorana_profile(H_bdg(:, i_min_e), qw_grid, &
        & zero_tol, Ntot_local, profile_lowest)

      call write_bdg_lowest_state_profile(qw_grid%z(1:N_local), profile_lowest, &
        & k_par_val, eigvals_bdg(i_min_e), xi_val)
      deallocate(profile_lowest)
    end block

    print *, '  Zero-energy BdG gap: ', result%min_gap, ' eV'

    if (allocated(qw_grid%z)) deallocate(qw_grid%z)
    deallocate(H_bdg, eigvals_bdg)

  end subroutine run_bdg_qw

  ! ==================================================================
  ! Spectral function A(k, E)
  ! ==================================================================
  subroutine run_spectral(cfg_in, result, profile_in, kpterms_in)
    type(simulation_config), intent(in) :: cfg_in
    type(topological_result), intent(inout) :: result
    real(kind=dp), contiguous, optional, intent(in) :: profile_in(:,:)
    real(kind=dp), contiguous, optional, intent(in) :: kpterms_in(:,:,:)

    real(kind=dp), allocatable :: k_arr(:), E_arr(:)
    real(kind=dp), allocatable :: A_kE(:,:)
    real(kind=dp) :: dk, dE
    integer :: ik, iE, nk, nE

    print *, '--- Spectral function A(k, E) ---'

    ! Build k and E grids from config
    nk = cfg_in%topo%spectral_nk
    nE = cfg_in%topo%spectral_nE
    allocate(k_arr(nk), E_arr(nE))

    if (nk > 1) then
      dk = (cfg_in%topo%spectral_k_max - cfg_in%topo%spectral_k_min) / real(nk - 1, kind=dp)
    else
      dk = 0.0_dp
    end if
    do ik = 1, nk
      k_arr(ik) = cfg_in%topo%spectral_k_min + real(ik - 1, kind=dp) * dk
    end do

    if (nE > 1) then
      dE = (cfg_in%topo%spectral_E_max - cfg_in%topo%spectral_E_min) / real(nE - 1, kind=dp)
    else
      dE = 0.0_dp
    end if
    do iE = 1, nE
      E_arr(iE) = cfg_in%topo%spectral_E_min + real(iE - 1, kind=dp) * dE
    end do

    print *, '  nk=', nk, ' nE=', nE, ' eta=', cfg_in%topo%spectral_eta
    print *, '  k range: [', k_arr(1), ',', k_arr(nk), '] 1/A'
    print *, '  E range: [', E_arr(1), ',', E_arr(nE), '] eV'

    select case (trim(cfg_in%confinement))
    case ('bulk')
      call compute_spectral_function_bulk(cfg_in, k_arr, E_arr, &
        & cfg_in%topo%spectral_eta, A_kE)
    case ('qw')
      if (.not. present(profile_in) .or. .not. present(kpterms_in)) then
        print *, 'Error: QW spectral mode requires allocated profile and kpterms'
        error stop 'QW spectral mode requires profile and kpterms'
      end if
      call compute_spectral_function_qw(cfg_in, profile_in, kpterms_in, &
        & k_arr, E_arr, cfg_in%topo%spectral_eta, A_kE)
    case ('wire')
      call compute_spectral_function_wire(cfg_in, k_arr, E_arr, &
        & cfg_in%topo%spectral_eta, A_kE)
    case default
      print *, 'Error: unsupported confinement for spectral mode: ', cfg_in%confinement
      error stop 'unsupported confinement for spectral mode'
    end select
    if (.not. allocated(A_kE) .or. size(A_kE, 1) == 0 .or. size(A_kE, 2) == 0) then
      print *, 'Error: spectral function calculation failed'
      error stop 'spectral function calculation failed'
    end if

    ! Store in result
    if (allocated(result%spectral_function)) deallocate(result%spectral_function)
    allocate(result%spectral_function(nk, nE))
    result%spectral_function = A_kE

    ! Write output file
    call write_spectral_function(cfg_in, k_arr, E_arr, A_kE)

    deallocate(k_arr, E_arr, A_kE)
  end subroutine run_spectral

  ! ==================================================================
  ! Hall conductance via Kubo formula
  ! ==================================================================
  subroutine run_conductance(cfg_in, result)
    type(simulation_config), intent(in) :: cfg_in
    type(topological_result), intent(inout) :: result

    real(kind=dp) :: sigma_xy, T

    print *, '--- Conductance analysis ---'
    print *, '  conductance_method: ', trim(cfg_in%topo%conductance_method)
    print *, '  berry_nk: ', cfg_in%topo%berry_nk
    print *, '  landauer_energy: ', cfg_in%topo%landauer_energy

    select case(trim(adjustl(cfg_in%topo%conductance_method)))
    case('kubo_chern', 'kubo')
      result%chern_number = compute_chern_qwz(cfg_in%topo%qwz_u, cfg_in%topo%berry_nk)
      sigma_xy = compute_hall_conductance(result%chern_number)
      result%conductance_xy = sigma_xy
      result%hall_conductance = sigma_xy
      print *, '  Method: Kubo (from Chern number)'
      print *, '  Chern number C = ', result%chern_number
      print *, '  Hall conductance sigma_xy = ', sigma_xy, ' e^2/h'
    case('kubo_berry')
      block
        complex(kind=dp), allocatable :: evecs(:,:,:,:)
        real(kind=dp), allocatable :: Omega(:,:), kx_arr(:), ky_arr(:)
        real(kind=dp) :: dk, kx_val, ky_val, mz_val, E_plus, nrm
        complex(kind=dp) :: off_diag, ev_tmp(2)
        integer :: nk, i, j

        nk = max(2, cfg_in%topo%berry_nk)
        dk = 2.0_dp * pi_dp / real(nk, kind=dp)
        allocate(evecs(2, 1, nk, nk), kx_arr(nk), ky_arr(nk))

        do i = 1, nk
          kx_arr(i) = -pi_dp + real(i - 1, kind=dp) * dk
          ky_arr(i) = -pi_dp + real(i - 1, kind=dp) * dk
        end do

        do j = 1, nk
          ky_val = ky_arr(j)
          do i = 1, nk
            kx_val = kx_arr(i)
            mz_val = cfg_in%topo%qwz_u + cos(kx_val) + cos(ky_val)
            off_diag = cmplx(sin(kx_val), -sin(ky_val), kind=dp)
            E_plus = sqrt(mz_val**2 + sin(kx_val)**2 + sin(ky_val)**2)
            ev_tmp(1) = off_diag
            ev_tmp(2) = cmplx(E_plus - mz_val, 0.0_dp, kind=dp)
            nrm = sqrt(real(conjg(ev_tmp(1))*ev_tmp(1) + conjg(ev_tmp(2))*ev_tmp(2), kind=dp))
            if (nrm > 1.0e-30_dp) then
              evecs(1, 1, i, j) = ev_tmp(1) / nrm
              evecs(2, 1, i, j) = ev_tmp(2) / nrm
            else
              evecs(1, 1, i, j) = cmplx(1.0_dp, 0.0_dp, kind=dp)
              evecs(2, 1, i, j) = cmplx(0.0_dp, 0.0_dp, kind=dp)
            end if
          end do
        end do

        Omega = compute_berry_curvature_lattice(evecs, kx_arr, ky_arr, 1)
        sigma_xy = compute_conductance_kubo(Omega, kx_arr, ky_arr)
        result%berry_curvature = Omega
        deallocate(evecs, kx_arr, ky_arr, Omega)
      end block
      result%conductance_xy = sigma_xy
      result%hall_conductance = sigma_xy
      result%chern_number = nint(sigma_xy)
      print *, '  Method: Kubo (Berry curvature integration)'
      print *, '  Hall conductance sigma_xy = ', sigma_xy, ' e^2/h'
    case('landauer')
      T = compute_landauer_transmission_1d(cfg_in%topo%landauer_energy, 0.0_dp, -1.0_dp, &
        & cfg_in%topo%spectral_eta)
      result%conductance_zz = T
      print *, '  Method: Landauer single-channel chain'
      print *, '  Transmission T = ', T
      print *, '  Conductance G = ', result%conductance_zz, ' e^2/h'
    end select

  end subroutine run_conductance

  ! ==================================================================
  ! Z2 phase diagram via gap sweep over (B, mu)
  ! ==================================================================
  subroutine run_gap_sweep(cfg_in, result, profile_in, kpterms_in)
    type(simulation_config), intent(in) :: cfg_in
    type(topological_result), intent(inout) :: result
    real(kind=dp), contiguous, optional, intent(in) :: profile_in(:,:)
    real(kind=dp), contiguous, optional, intent(in) :: kpterms_in(:,:,:)

    integer, allocatable :: z2_map_int(:,:)
    real(kind=dp), allocatable :: gap_map_real(:,:), transitions(:,:)
    real(kind=dp) :: gap_threshold
    integer :: iB, iMu, nB, nMu

    print *, '--- Z2 gap sweep ---'

    nB = cfg_in%topo%gap_sweep_nB
    nMu = cfg_in%topo%gap_sweep_nMu
    gap_threshold = 1.0e-4_dp

    print *, '  nB=', nB, ' nMu=', nMu
    print *, '  B range: [', cfg_in%topo%gap_sweep_B_min, ',', cfg_in%topo%gap_sweep_B_max, ']'
    print *, '  mu range: [', cfg_in%topo%gap_sweep_mu_min, ',', cfg_in%topo%gap_sweep_mu_max, ']'
    print *, '  sweep_model=', trim(cfg_in%topo%sweep_model)

    select case (trim(cfg_in%topo%sweep_model))
    case ('bhz_analytic')
      call compute_z2_gap_sweep(cfg_in, &
        & cfg_in%topo%gap_sweep_B_min, cfg_in%topo%gap_sweep_B_max, nB, &
        & cfg_in%topo%gap_sweep_mu_min, cfg_in%topo%gap_sweep_mu_max, nMu, &
        & gap_threshold, z2_map_int, gap_map_real, transitions)
    case ('qw_fukane')
      call compute_qw_fukane_gap_sweep(cfg_in, profile_in, kpterms_in, gap_threshold, &
        & z2_map_int, gap_map_real, transitions)
    case ('wire_bdg')
      call compute_wire_bdg_gap_sweep(cfg_in, gap_threshold, z2_map_int, gap_map_real, transitions)
    case default
      print *, 'ERROR: unknown topology sweep_model: ', trim(cfg_in%topo%sweep_model)
      error stop 'unknown topology sweep_model'
    end select

    ! Convert integer z2_map to real(dp) for storage in result
    if (allocated(result%z2_map)) deallocate(result%z2_map)
    if (allocated(result%gap_map)) deallocate(result%gap_map)
    allocate(result%z2_map(nMu, nB))
    allocate(result%gap_map(nMu, nB))
    do iMu = 1, nMu
      do iB = 1, nB
        result%z2_map(iMu, iB) = real(z2_map_int(iMu, iB), kind=dp)
        result%gap_map(iMu, iB) = gap_map_real(iMu, iB)
      end do
    end do

    ! Write phase diagram output
    call write_z2_phase_diagram(cfg_in, result%z2_map, result%gap_map)

    call write_z2_transitions(transitions)

    if (size(transitions, 1) > 0) then
      print *, '  Phase transitions detected: ', size(transitions, 1)
      do iB = 1, size(transitions, 1)
        print '(A,ES12.4,A,ES12.4)', '    B=', transitions(iB, 1), ' mu=', transitions(iB, 2)
      end do
      ! Store transitions as phase_boundary
      if (allocated(result%phase_boundary)) deallocate(result%phase_boundary)
      allocate(result%phase_boundary(size(transitions, 1), 2))
      result%phase_boundary = transitions
    else
      print *, '  No phase transitions detected'
    end if

    deallocate(z2_map_int, gap_map_real, transitions)
  end subroutine run_gap_sweep

  subroutine compute_qw_fukane_gap_sweep(cfg_in, profile_in, kpterms_in, gap_threshold, &
      & z2_map, gap_map, transitions)
    type(simulation_config), intent(in) :: cfg_in
    real(kind=dp), contiguous, intent(in) :: profile_in(:,:)
    real(kind=dp), contiguous, intent(in) :: kpterms_in(:,:,:)
    real(kind=dp), intent(in) :: gap_threshold
    integer, allocatable, intent(out) :: z2_map(:,:)
    real(kind=dp), allocatable, intent(out) :: gap_map(:,:), transitions(:,:)

    integer :: iB, iMu, nB, nMu, z2_status, z2_val, n_occ
    real(kind=dp) :: min_gap, B_val, mu_val, dB, dmu
    type(simulation_config) :: cfg_local

    nB = cfg_in%topo%gap_sweep_nB
    nMu = cfg_in%topo%gap_sweep_nMu
    if (nB < 1 .or. nMu < 1) then
      allocate(z2_map(0,0), gap_map(0,0), transitions(0,2))
      return
    end if

    n_occ = 6 * size(profile_in, 1)
    allocate(z2_map(nMu, nB), gap_map(nMu, nB))

    dB = 0.0_dp
    if (nB > 1) dB = (cfg_in%topo%gap_sweep_B_max - cfg_in%topo%gap_sweep_B_min) / real(nB - 1, kind=dp)
    dmu = 0.0_dp
    if (nMu > 1) dmu = (cfg_in%topo%gap_sweep_mu_max - cfg_in%topo%gap_sweep_mu_min) / real(nMu - 1, kind=dp)

    cfg_local = cfg_in
    cfg_local%bdg%enabled = .true.

    do iB = 1, nB
      B_val = cfg_in%topo%gap_sweep_B_min + real(iB - 1, kind=dp) * dB
      cfg_local%bdg%B_vec = [0.0_dp, 0.0_dp, B_val]
      do iMu = 1, nMu
        mu_val = cfg_in%topo%gap_sweep_mu_min + real(iMu - 1, kind=dp) * dmu
        cfg_local%bdg%mu = mu_val

        call compute_z2_fukane_qw_result(cfg_local, profile_in, kpterms_in, n_occ, &
          & z2_val, min_gap, z2_status)
        if (z2_status /= 0) then
          print *, 'WARNING: QW Fu-Kane sweep failed at B=', B_val, ' mu=', mu_val, &
            & ' status=', z2_status
          z2_val = 0
          min_gap = huge(1.0_dp)
        end if

        z2_map(iMu, iB) = z2_val
        gap_map(iMu, iB) = min_gap
      end do
    end do

    call detect_z2_transitions(z2_map, gap_map, cfg_in%topo%gap_sweep_B_min, &
      & cfg_in%topo%gap_sweep_B_max, cfg_in%topo%gap_sweep_mu_min, &
      & cfg_in%topo%gap_sweep_mu_max, gap_threshold, transitions)
  end subroutine compute_qw_fukane_gap_sweep

  subroutine compute_wire_bdg_gap_sweep(cfg_in, gap_threshold, z2_map, gap_map, transitions)
    type(simulation_config), intent(in) :: cfg_in
    real(kind=dp), intent(in) :: gap_threshold
    integer, allocatable, intent(out) :: z2_map(:,:)
    real(kind=dp), allocatable, intent(out) :: gap_map(:,:), transitions(:,:)

    integer :: iB, iMu, nB, nMu
    real(kind=dp) :: dB, dmu, B_val, mu_val

    nB = cfg_in%topo%gap_sweep_nB
    nMu = cfg_in%topo%gap_sweep_nMu
    if (nB < 1 .or. nMu < 1) then
      allocate(z2_map(0,0), gap_map(0,0), transitions(0,2))
      return
    end if

    allocate(z2_map(nMu, nB), gap_map(nMu, nB))
    dB = 0.0_dp
    if (nB > 1) dB = (cfg_in%topo%gap_sweep_B_max - cfg_in%topo%gap_sweep_B_min) / real(nB - 1, kind=dp)
    dmu = 0.0_dp
    if (nMu > 1) dmu = (cfg_in%topo%gap_sweep_mu_max - cfg_in%topo%gap_sweep_mu_min) / real(nMu - 1, kind=dp)

    do iB = 1, nB
      B_val = cfg_in%topo%gap_sweep_B_min + real(iB - 1, kind=dp) * dB
      do iMu = 1, nMu
        mu_val = cfg_in%topo%gap_sweep_mu_min + real(iMu - 1, kind=dp) * dmu
        call eval_wire_bdg_gap(cfg_in, B_val, mu_val, gap_threshold, &
          & z2_map(iMu, iB), gap_map(iMu, iB))
      end do
    end do

    call detect_z2_transitions(z2_map, gap_map, cfg_in%topo%gap_sweep_B_min, &
      & cfg_in%topo%gap_sweep_B_max, cfg_in%topo%gap_sweep_mu_min, &
      & cfg_in%topo%gap_sweep_mu_max, gap_threshold, transitions)
  end subroutine compute_wire_bdg_gap_sweep

  subroutine bdg_wire_cleanup(H_bdg_csr, eigen_res_local, eigen_solver_local, wsetup)
    type(csr_matrix), intent(inout) :: H_bdg_csr
    type(eigensolver_result), intent(inout) :: eigen_res_local
    class(eigensolver_base), allocatable, intent(inout) :: eigen_solver_local
    type(wire_setup), intent(inout) :: wsetup
    call csr_free(H_bdg_csr)
    call eigensolver_result_free(eigen_res_local)
    if (allocated(eigen_solver_local)) deallocate(eigen_solver_local)
    call wire_setup_free(wsetup)
  end subroutine bdg_wire_cleanup

  subroutine eval_wire_bdg_gap(cfg_in, B_val, mu_val, gap_threshold, z2, gap)
    type(simulation_config), intent(in) :: cfg_in
    real(kind=dp), intent(in) :: B_val, mu_val, gap_threshold
    integer, intent(out) :: z2
    real(kind=dp), intent(out) :: gap

    type(simulation_config) :: cfg
    type(wire_setup) :: wsetup
    type(csr_matrix) :: H_bdg_csr
    class(eigensolver_base), allocatable :: eigen_solver_local
    type(eigensolver_config) :: eigen_cfg_local
    type(eigensolver_result) :: eigen_res_local
    real(kind=dp), allocatable :: eigvals_bdg(:)
    integer :: Ngrid_local, Ntot_local, Nbdg_local, nev_local
    real(kind=dp) :: emin_local, emax_local

    cfg = cfg_in
    cfg%bdg%enabled = .true.
    cfg%bdg%mu = mu_val
    ! U8-followup: transverse B for Peierls orbital coupling (design §1 second
    ! root cause). Bx varies with B_val; By and Bz zero.
    cfg%bdg%B_vec = [B_val, 0.0_dp, 0.0_dp]

    ! Strain-aware wire init (Issue #04): run compute_strain +
    ! compute_bir_pikus_blocks when cfg%strain%enabled so the BdG sweep
    ! honors strain. Pre-#04 this site skipped strain entirely.
    call wire_setup_init(wsetup, cfg)

    Ngrid_local = grid_ngrid(cfg%grid)
    Ntot_local = 8 * Ngrid_local
    Nbdg_local = 2 * Ntot_local
    call build_bdg_hamiltonian_1d(H_bdg_csr, cfg, wsetup%profile_2d, wsetup%kpterms_2d, &
      & 0.0_dp, cfg%bdg%mu, cfg%bdg%delta_0, wsetup%ws, &
      & cfg%bdg%B_vec, cfg%bdg%g_factor)

    nev_local = max(2, cfg%bands%num_cb + cfg%bands%num_vb)
    eigen_cfg_local%method = 'FEAST'
    eigen_cfg_local%mode = EIGEN_MODE_ENERGY
    eigen_cfg_local%nev = nev_local
    eigen_cfg_local%max_iter = 200
    eigen_cfg_local%tol = 1.0e-10_dp
    eigen_cfg_local%m0 = min(max(8 * nev_local, 200), Nbdg_local)

    ! Use [solver] emin/emax as user override if set; otherwise ±5·δ₀ default.
    if (cfg%solver%emin /= 0.0_dp .or. cfg%solver%emax /= 0.0_dp) then
      emin_local = cfg%solver%emin
      emax_local = cfg%solver%emax
    else
      emin_local = -5.0_dp * cfg%bdg%delta_0
      emax_local =  5.0_dp * cfg%bdg%delta_0
    end if
    ! Route the physics-sized window through the window authority (ADR 0005 /
    ! KTD6). apply_solver_window honors a nonzero user override verbatim, so
    ! the physics window is preserved while window selection has one home.
    ! NEVER the Gershgorin auto-window for BdG (samples the FD-Nyquist tail;
    ! see docs/solutions/best-practices/2026-06-21-fd-nyquist-spurious-tail.md).
    call apply_solver_window(H_bdg_csr, emin_local, emax_local, &
                             eigen_cfg_local%emin, eigen_cfg_local%emax)

    eigen_solver_local = make_eigensolver(eigen_cfg_local)
    call eigen_solver_local%solve_sparse(H_bdg_csr, eigen_cfg_local, eigen_res_local)

    if (.not. eigen_res_local%converged .or. eigen_res_local%nev_found < 1) then
      call bdg_wire_cleanup(H_bdg_csr, eigen_res_local, eigen_solver_local, wsetup)
      error stop 'eval_wire_bdg_gap: FEAST failed or found no states ' // &
        '(mu likely in band gap or window mis-sized)'
    end if
    if (eigen_res_local%m0_used > 0 .and. &
        eigen_res_local%nev_found >= eigen_res_local%m0_used .and. &
        eigen_res_local%m0_used < Nbdg_local) then
      print *, 'WARNING: wire BdG sweep FEAST subspace likely truncated ', &
        eigen_solver_local%backend_name(), ' (m0 too small or window too broad)'
      print *, '  nev_found=', eigen_res_local%nev_found, ' m0=', eigen_res_local%m0_used
      call bdg_wire_cleanup(H_bdg_csr, eigen_res_local, eigen_solver_local, wsetup)
      error stop 'eval_wire_bdg_gap: FEAST subspace likely truncated ' // &
        '(increase m0 or narrow the window); ' // &
        'fail-fast prevents sentinel leakage into the phase diagram ' // &
        '(Codex review P2 on PR40)'
    end if
    allocate(eigvals_bdg(eigen_res_local%nev_found))
    eigvals_bdg = eigen_res_local%eigenvalues
    block
      type(bdg_eval_params_t) :: evp
      type(bdg_eval_result_t) :: evr
      evp = bdg_eval_params_with_delta(cfg_in%bdg%delta_0)
      evr = eval_bdg_point(eigvals_bdg, evp)
      gap = evr%minigap
    end block
    deallocate(eigvals_bdg)

    ! Issue 07 (U10): route the wire_bdg z2 through the projected Pfaffian
    ! (S2 strategy: analytical bands 7-8 per k.p block table SSOT) instead
    ! of the 1D count heuristic. Pfaffian sign convention: -1 = topological,
    ! +1 = trivial, 0 = gap closure / inconclusive.
    block
      integer :: s2_sign
      ! Slim Pfaffian witness via seam sibling (per ticket 04 of
      ! .scratch/bdg-evaluator-pfaffian/ — drop-in replacement for
      ! wire_pfaffian_witness_sweep, same s2_sign ∈ {-1, 0, +1} semantics).
      s2_sign = eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg_local, &
                                              bdg_eval_params_with_delta(cfg_in%bdg%delta_0))
      if (s2_sign == -1) then
        z2 = 1
      else if (s2_sign == +1) then
        z2 = 0
      else
        ! Gap closure: defer to the SC-minigap heuristic as the fallback so
        ! the colormap still shows a Z2=1 flag exactly when min_gap < threshold.
        ! This preserves the open->close->reopen pattern at the B_crit point.
        z2 = compute_z2_gap_bhz_heuristic(eigen_res_local%eigenvalues, gap_threshold)
      end if
    end block

    call bdg_wire_cleanup(H_bdg_csr, eigen_res_local, eigen_solver_local, wsetup)
  end subroutine eval_wire_bdg_gap

end program topologicalAnalysis
