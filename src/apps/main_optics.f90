program opticalProperties

  use definitions, only: NUM_VB_STATES, dp, pi_dp, simulation_config, &
    validate_semantic, wavevector
  use parameters
  use hamiltonianConstructor
  use input_parser, only: read_config
  use simulation_setup_mod, only: simulation_setup, simulation_setup_init, &
    & simulation_setup_free, setup_build_velocity_matrices, setup_build_H, &
    & derive_eigensolver
  use optical_spectra
  use exciton_solver, only: compute_exciton_binding, apply_excitonic_corrections
  use sparse_matrices
  use eigensolver, only: eigensolver_result, eigensolver_result_free, &
    & eigensolver_config, eigensolver_base, make_eigensolver, &
    & eigensolver_config_validate, apply_solver_window, EIGEN_MODE_INDEX
  use hamiltonian_qw, only: qw_workspace, qw_workspace_free, ZB8bandQW_csr
  use linalg, only: mkl_set_num_threads_local
  use utils, only: ensure_output_dir, get_unit

  implicit none

  ! Shared configuration from input_parser
  type(simulation_config) :: cfg
  type(simulation_setup)  :: setup
  type(optics_engine)     :: oe

  ! Iteration
  integer :: i, k

  ! Hamiltonian and eigensolver
  integer :: info, N, il, iuu
  real(kind=dp) :: vl, vu
  real(kind=dp), allocatable :: eig(:,:)
  complex(kind=dp), allocatable, dimension(:,:,:) :: eigv
  complex(kind=dp), allocatable, dimension(:,:) :: HT

  ! k_par sweep
  type(wavevector), allocatable, dimension(:) :: smallk
  integer :: npts
  real(kind=dp) :: dk

  ! Simpson integration weights
  real(kind=dp), allocatable :: simpson_w(:)

  ! Wire-specific variables (confinement='wire')
  type(eigensolver_result) :: eigen_res
  integer                  :: Ngrid, Ntot, nev_wire

  ! File handling
  integer(kind=4) :: iounit

  ! ================================================================
  ! Parse input, initialize materials, external fields
  ! ================================================================
  call read_config(cfg)

  ! Semantic validation (optics block must be enabled)
  call validate_semantic(cfg, 'opticalProperties')

  ! ================================================================
  ! Branch by confinement type
  ! ================================================================
  select case (trim(cfg%confinement))

  ! ==================================================================
  ! BULK (confinement='bulk')
  ! ==================================================================
  case ('bulk')
    call simulation_setup_init(cfg, setup)
    N = setup%N   ! 8 for bulk
    il = setup%il
    iuu = setup%iuu
    if (setup%sc_was_run) cfg%sc%fermi_level = setup%fermi_level
    print '(a)', ' Bulk optics: 8x8 Hamiltonian, all 8 states'

    ! Build velocity matrices via simulation_setup
    call setup_build_velocity_matrices(setup, cfg)
    print '(a)', ' Bulk velocity matrices built successfully'

    vl = 0.0_dp
    vu = 0.0_dp

    ! Allocate eigenvalue storage
    allocate(eig(iuu - il + 1, cfg%wave_vector%nsteps))
    allocate(eigv(N, iuu - il + 1, cfg%wave_vector%nsteps))
    eig(:,:) = 0.0_dp
    eigv(:,:,:) = (0.0_dp, 0.0_dp)

    allocate(HT(N, N))
    HT = (0.0_dp, 0.0_dp)

    ! Build k-sweep array: 1D sweep from 0 to waveVectorMax
    npts = cfg%wave_vector%nsteps
    allocate(smallk(npts))
    smallk%kx = 0.0_dp
    smallk%ky = 0.0_dp
    smallk%kz = 0.0_dp
    do k = 1, npts
      smallk(k)%kz = real(k - 1, kind=dp) * cfg%wave_vector%max &
        & / real(npts - 1, kind=dp)
    end do

    ! ================================================================
    ! Initialize optics accumulation
    ! ================================================================
    cfg%optics%confinement = cfg%confinement
    call optics_init(oe, cfg%optics)

    ! ================================================================
    ! Simpson integration weights for 3D spherical k-space integration
    !
    ! For bulk, the absorption involves a 3D BZ integral.
    ! Using spherical symmetry along the z-sweep direction:
    !   int d^3k = 4*pi * int k^2 dk
    ! with additional (2*pi)^3 normalization.
    !
    ! weight = Simpson * 4*pi*k^2 / (2*pi)^3
    ! ================================================================
    if (mod(npts, 2) == 0) then
      npts = npts - 1  ! Simpson requires odd number of points
      print '(a,i0,a)', ' Warning: optics integration uses ', npts, &
        & ' k-points (Simpson requires odd count)'
    end if
    allocate(simpson_w(npts))
    dk = cfg%wave_vector%max / real(cfg%wave_vector%nsteps - 1, kind=dp)
    do i = 1, npts
      ! Base Simpson 1/3 rule weight
      if (i == 1 .or. i == npts) then
        simpson_w(i) = dk / 3.0_dp
      else if (mod(i, 2) == 0) then
        simpson_w(i) = 4.0_dp * dk / 3.0_dp
      else
        simpson_w(i) = 2.0_dp * dk / 3.0_dp
      end if
      ! Multiply by 4*pi*k^2 / (2*pi)^3 for 3D spherical BZ integration
      simpson_w(i) = simpson_w(i) * 4.0_dp * pi_dp &
        & * (real(i - 1, kind=dp) * dk)**2 / (2.0_dp * pi_dp)**3
    end do

    ! ================================================================
    ! k-sweep: diagonalize and accumulate spectra
    ! ================================================================
    info = mkl_set_num_threads_local(1)

    print '(a,i0,a)', ' Bulk optics k-sweep: ', npts, ' k-points'

    block
      type(eigensolver_config) :: bulk_cfg
      class(eigensolver_base), allocatable :: bulk_solver
      type(eigensolver_result) :: bulk_result
      complex(kind=dp), allocatable :: HT_loc(:,:)
      integer :: M_loc

      bulk_cfg%method = 'DENSE'
      bulk_cfg%mode = EIGEN_MODE_INDEX
      bulk_cfg%il = il
      bulk_cfg%iu = iuu
      bulk_cfg%nev = iuu - il + 1

      !$omp parallel firstprivate(bulk_cfg) private(k, HT_loc, bulk_solver, bulk_result, M_loc)
      allocate(HT_loc(N, N))
      HT_loc = (0.0_dp, 0.0_dp)
      bulk_solver = make_eigensolver(bulk_cfg)

      !$omp do schedule(static)
      do k = 1, npts
        ! Build 8x8 bulk Hamiltonian at this k
        call ZB8bandBulk(HT_loc, smallk(k), cfg%params(1:1), cfg=cfg)

        ! Diagonalize all 8 eigenvalues
        call bulk_solver%solve_dense(HT_loc, bulk_cfg, bulk_result)
        if (.not. bulk_result%converged) then
          !$omp critical
          print '(a,i0)', ' ERROR: diagonalization failed at k=', k
          !$omp end critical
          error stop 'eigensolver failed in bulk optics k-sweep'
        end if

        M_loc = bulk_result%nev_found
        eig(1:M_loc, k) = bulk_result%eigenvalues(1:M_loc)
        eigv(:, 1:M_loc, k) = bulk_result%eigenvectors(:, 1:M_loc)
        call eigensolver_result_free(bulk_result)
      end do
      !$omp end do

      deallocate(HT_loc)
      !$omp end parallel
    end block

    ! Accumulate optical spectra (serial, uses setup%vel)
    print '(a)', ' Accumulating optical spectra...'
    do k = 1, npts
      call optics_accumulate(oe, &
        & eig(:, k), eigv(:, :, k), simpson_w(k), &
        & setup%vel, cfg%bands%num_cb, cfg%bands%num_vb, cfg%sc%fermi_level)

      if (cfg%optics%spontaneous_enabled) then
        call optics_accumulate_spontaneous(oe, &
          & eig(:, k), eigv(:, :, k), simpson_w(k), &
          & setup%vel, cfg%bands%num_cb, cfg%bands%num_vb, cfg%sc%fermi_level)
      end if

      if (cfg%optics%gain_enabled) then
        call compute_gain_qw(oe, &
          & eig(:, k), eigv(:, :, k), simpson_w(k), &
          & setup%vel, cfg%bands%num_cb, cfg%bands%num_vb, &
          & cfg%optics%gain_carrier_density)
      end if

      if (cfg%optics%isbt_enabled) then
        call compute_isbt_absorption(oe, &
          & eig(:, k), eigv(:, :, k), &
          & setup%vel, cfg%bands%num_cb, cfg%bands%num_vb, &
          & simpson_w(k), cfg%sc%fermi_level)
      end if
    end do

    ! ================================================================
    ! Finalize: apply prefactor, write output files
    ! ================================================================
    call optics_apply_prefactor(oe)
    call optics_write_output(oe)
    call optics_free(oe)

    deallocate(simpson_w)
    call simulation_setup_free(setup)
    print '(a)', ' Bulk optical spectra written to output/'

  ! ==================================================================
  ! QUANTUM WELL (confinement='qw')
  ! ==================================================================
  case ('qw')

    call simulation_setup_init(cfg, setup)
    N = setup%N
    if (setup%sc_was_run) cfg%sc%fermi_level = setup%fermi_level
    ! Compute eigenvalue range for optics: top numvb valence + bottom numcb conduction
    il = NUM_VB_STATES*cfg%grid%npoints() - cfg%bands%num_vb + 1
    iuu = NUM_VB_STATES*cfg%grid%npoints() + cfg%bands%num_cb
    print '(a,i0,a,i0)', ' Computing states from index ', il, ' to ', iuu

    ! Build velocity matrices via simulation_setup
    call setup_build_velocity_matrices(setup, cfg)
    print '(a)', ' Velocity matrices built successfully'

    vl = 0.0_dp
    vu = 0.0_dp

    ! Allocate eigenvalue storage
    allocate(eig(iuu-il+1, cfg%wave_vector%nsteps))
    allocate(eigv(N, iuu-il+1, cfg%wave_vector%nsteps))
    eig(:,:) = 0.0_dp
    eigv(:,:,:) = (0.0_dp, 0.0_dp)

    allocate(HT(N, N))
    HT = (0.0_dp, 0.0_dp)

    ! Build k_par array: uniform sweep from 0 to waveVectorMax
    npts = cfg%wave_vector%nsteps
    allocate(smallk(npts))
    smallk%kx = 0.0_dp
    smallk%ky = 0.0_dp
    smallk%kz = 0.0_dp
    do k = 1, npts
      smallk(k)%kx = real(k - 1, kind=dp) * cfg%wave_vector%max &
        & / real(npts - 1, kind=dp)
    end do

    ! --- Write potential profile ---
    call ensure_output_dir()
    call get_unit(iounit)
    open(unit=iounit, file='output/potential_profile.dat', status='replace', action='write')
    do i = 1, cfg%grid%npoints()
      write(iounit, *) cfg%z(i), setup%profile(i,1), setup%profile(i,2), setup%profile(i,3)
    end do
    close(iounit)

    ! ================================================================
    ! Initialize optics accumulation
    ! ================================================================
    cfg%optics%confinement = cfg%confinement
    call optics_init(oe, cfg%optics)

    ! ================================================================
    ! Simpson integration weights for 2D cylindrical k_par integration
    ! int d^2k = 2*pi * int k dk  =>  weight = Simpson * 2*pi*k
    ! ================================================================
    if (mod(npts, 2) == 0) then
      npts = npts - 1  ! Simpson requires odd number of points
      print '(a,i0,a)', ' Warning: optics integration uses ', npts, &
        & ' k-points (Simpson requires odd count)'
    end if
    allocate(simpson_w(npts))
    dk = cfg%wave_vector%max / real(cfg%wave_vector%nsteps - 1, kind=dp)
    do i = 1, npts
      ! Base Simpson 1/3 rule weight
      if (i == 1 .or. i == npts) then
        simpson_w(i) = dk / 3.0_dp
      else if (mod(i, 2) == 0) then
        simpson_w(i) = 4.0_dp * dk / 3.0_dp
      else
        simpson_w(i) = 2.0_dp * dk / 3.0_dp
      end if
      ! Multiply by 2*pi*k for 2D cylindrical integration
      simpson_w(i) = simpson_w(i) * 2.0_dp * pi_dp * real(i - 1, kind=dp) * dk
    end do

    ! ================================================================
    ! k_par sweep: diagonalize and accumulate spectra
    ! ================================================================
    ! Disable MKL internal threading (serial eigensolver per k-point)
    info = mkl_set_num_threads_local(1)

    print '(a,i0,a)', ' QW optics k-sweep: ', npts, ' k-points'

    block
      type(eigensolver_config) :: qw_cfg
      class(eigensolver_base), allocatable :: qw_solver_outer
      type(eigensolver_result) :: qw_result
      complex(kind=dp), allocatable :: HT_loc(:,:)
      integer :: M_loc

      ! --- Build the k_par-sweep solver config via the single derivation
      ! seam (issue #02). This honors cfg%solver%method: AUTO resolves to
      ! DENSE+INDEX for QW (the golden default, preserved bit-identical),
      ! and an explicit method=FEAST resolves to FEAST+ENERGY. The
      ! constructed qw_solver_outer validates the config and is used only
      ! for the backend diagnostic; per-thread solvers are built inside the
      ! OpenMP region below so each thread owns its own FEAST workspace.
      call derive_eigensolver(cfg, N=N, il=il, iu=iuu, nev=iuu-il+1, &
                              eigen_cfg=qw_cfg, solver=qw_solver_outer)

      ! --- QW + FEAST: dispersion-aware sweep-envelope window (ADR 0005,
      ! issue #03). The optics k_par sweep is the same shape as the
      ! band-structure k-sweep: ONE stable window per sweep, derived by
      ! unioning the Gershgorin bounds at the two sweep endpoints
      ! (smallk(1) = k_par=0, smallk(npts) = k_par=max). DENSE ignores
      ! emin/emax, so this block is a no-op for the golden default. The
      ! authority honors a user-set [solver] window verbatim (no envelope),
      ! so the envelope is derived ONLY when the user left the window at
      ! the parser's [0,0] auto sentinel.
      if (trim(qw_cfg%method) == 'FEAST' .and. &
          cfg%solver%emin == 0.0_dp .and. cfg%solver%emax == 0.0_dp) then
        block
          type(csr_matrix) :: HT_csr_a, HT_csr_b
          type(qw_workspace) :: qw_ws_env
          ! Endpoint k-points: smallk(1) = k_par=0, smallk(npts) = k_par=max.
          call ZB8bandQW_csr(HT_csr_a, smallk(1), setup%profile, &
            setup%kpterms, cfg, ws=qw_ws_env)
          call ZB8bandQW_csr(HT_csr_b, smallk(npts), setup%profile, &
            setup%kpterms, cfg, ws=qw_ws_env)
          call apply_solver_window(HT_csr_a, HT_csr_b, &
            user_emin=cfg%solver%emin, user_emax=cfg%solver%emax, &
            emin_out=qw_cfg%emin, emax_out=qw_cfg%emax)
          call csr_free(HT_csr_a)
          call csr_free(HT_csr_b)
          call qw_workspace_free(qw_ws_env)
          print '(A,ES12.4,A,ES12.4,A)', ' QW optics auto FEAST window '// &
            '(sweep envelope): [', qw_cfg%emin, ',', qw_cfg%emax, ']'
        end block
        ! Re-validate after window update; per-thread construction below
        ! picks up the updated window.
        call eigensolver_config_validate(qw_cfg)
      end if

      print '(a,a,a)', ' QW optics eigensolver: ', &
        trim(qw_cfg%method), ' backend'

      block
        ! Per-thread solver + config (mirrors main.f90's QW DENSE path:
        ! each thread owns its own FEAST workspace via solver_loc).
        type(eigensolver_config) :: cfg_loc
        class(eigensolver_base), allocatable :: solver_loc
        integer :: nev_target, idx_lo

        ! Number of gap-straddling states the optics accumulator expects
        ! (top numvb valence + bottom numcb conduction). eig/eigv are
        ! sized to exactly this count.
        nev_target = iuu - il + 1

        !$omp parallel private(k, HT_loc, cfg_loc, solver_loc, qw_result, M_loc, idx_lo)
        allocate(HT_loc(N, N))
        HT_loc = (0.0_dp, 0.0_dp)
        cfg_loc = qw_cfg
        solver_loc = make_eigensolver(cfg_loc)

        !$omp do schedule(static)
        do k = 1, npts
          ! Build Hamiltonian for this k_par
          call ZB8bandQW(HT_loc, smallk(k), setup%profile, setup%kpterms, cfg=cfg)

          ! Diagonalize. solve_dense dispatches to LAPACK (DENSE) or to FEAST
          ! via dense->CSR conversion (feast_solve_dense, issue #01).
          call solver_loc%solve_dense(HT_loc, cfg_loc, qw_result)
          if (.not. qw_result%converged) then
            !$omp critical
            print '(a,a,a,i0)', ' ERROR: ', solver_loc%backend_name(), &
              ' diagonalization failed at k=', k
            !$omp end critical
            error stop 'eigensolver failed in QW optics k-sweep'
          end if

          M_loc = qw_result%nev_found
          if (M_loc < nev_target) then
            !$omp critical
            print '(a,a,a,i0,a,i0,a,i0,a)', ' ERROR: ', &
              solver_loc%backend_name(), ' returned only ', M_loc, &
              ' eigenvalues at k=', k, ' (need ', nev_target, &
              '); widen the energy window'
            !$omp end critical
            error stop 'insufficient eigenvalues for QW optics accumulation'
          end if

          if (M_loc == nev_target) then
            ! DENSE+INDEX: the solver already returned exactly the requested
            ! [il, iuu] slice (LAPACK zheevx range='I'); store directly.
            idx_lo = 1
          else
            ! FEAST+ENERGY: the solver returned the full in-window spectrum
            ! (sorted ascending). Extract the gap-straddling [il, iuu]
            ! slice by GLOBAL index - the same states DENSE+INDEX would
            ! return. This requires the window to cover the full spectral
            ! range below iuu (the sweep-envelope authority guarantees this).
            idx_lo = il
          end if

          eig(1:nev_target, k) = qw_result%eigenvalues(idx_lo:idx_lo+nev_target-1)
          eigv(:, 1:nev_target, k) = &
            qw_result%eigenvectors(:, idx_lo:idx_lo+nev_target-1)
          call eigensolver_result_free(qw_result)
        end do
        !$omp end do

        deallocate(HT_loc)
        !$omp end parallel
      end block
    end block

    ! Accumulate optical spectra (serial, uses setup%vel)
    print '(a)', ' Accumulating optical spectra...'
    do k = 1, npts
      call optics_accumulate(oe, &
        & eig(:, k), eigv(:, :, k), simpson_w(k), &
        & setup%vel, cfg%bands%num_cb, cfg%bands%num_vb, cfg%sc%fermi_level)

      if (cfg%optics%spontaneous_enabled) then
        call optics_accumulate_spontaneous(oe, &
          & eig(:, k), eigv(:, :, k), simpson_w(k), &
          & setup%vel, cfg%bands%num_cb, cfg%bands%num_vb, cfg%sc%fermi_level)
      end if

      if (cfg%optics%gain_enabled) then
        call compute_gain_qw(oe, &
          & eig(:, k), eigv(:, :, k), simpson_w(k), &
          & setup%vel, cfg%bands%num_cb, cfg%bands%num_vb, &
          & cfg%optics%gain_carrier_density)
      end if

      if (cfg%optics%isbt_enabled) then
        call compute_isbt_absorption(oe, &
          & eig(:, k), eigv(:, :, k), &
          & setup%vel, cfg%bands%num_cb, cfg%bands%num_vb, &
          & simpson_w(k), cfg%sc%fermi_level)
      end if
    end do

    ! ISBT dipole transition table (zone-center only)
    if (cfg%optics%isbt_enabled) then
      call compute_intersubband_transitions(eig(:, 1), eigv(:, :, 1), &
        & cfg%z, cfg%dz, cfg%bands%num_cb, cfg%bands%num_vb, cfg%grid%npoints(), &
        & "output/isbt_transitions.dat")
    end if

    ! ================================================================
    ! Finalize: apply prefactor, optional exciton, write output files
    ! ================================================================
    call optics_apply_prefactor(oe)

    ! Exciton corrections (applied after prefactor, before write_output)
    if (cfg%exciton%enabled) then
      block
        real(kind=dp) :: E_binding_ex, lambda_opt_ex, E_gap_ex
        integer :: cb_st, vb_st
        cb_st = cfg%bands%num_vb + 1
        vb_st = cfg%bands%num_vb
        E_gap_ex = eig(cb_st, 1) - eig(vb_st, 1)
        call compute_exciton_binding(eig(:, 1), eigv(:, :, 1), &
          & cfg%z, cfg%dz, cfg%num_layers, cfg%params, &
          & cfg%bands%num_cb, cfg%bands%num_vb, cfg%grid%npoints(), E_binding_ex, lambda_opt_ex, &
          & cfg%grid%material_id)
        print '(a,f8.3,a)', ' Exciton binding energy: ', E_binding_ex, ' meV'
        call apply_excitonic_corrections(oe%E_grid, oe%alpha_te, oe%alpha_tm, &
          & E_gap_ex, E_binding_ex, cfg%optics)
        print '(a)', ' Excitonic corrections applied to absorption spectra'
      end block
    end if

    call optics_write_output(oe)
    call optics_free(oe)

    deallocate(simpson_w)
    call simulation_setup_free(setup)
    print '(a)', ' Optical spectra written to output/'

  ! ==================================================================
  ! WIRE (confinement='wire')
  ! ==================================================================
  case ('wire')

    ! ----------------------------------------------------------------
    ! Initialize via simulation_setup (2D confinement, strain, SC, eigensolver)
    ! ----------------------------------------------------------------
    call simulation_setup_init(cfg, setup)
    Ntot = setup%Ntot
    nev_wire = setup%nev_wire
    if (setup%sc_was_run) cfg%sc%fermi_level = setup%fermi_level

    Ngrid = setup%Ngrid

    print '(a)', ''
    print '(a)', '=== Wire optical properties (2D confinement) ==='
    print '(a,i0,a,i0,a,i0)', '  Grid: nx=', cfg%grid%nx, ' ny=', cfg%grid%ny, &
      & ' Ngrid=', Ngrid
    print '(a,i0,a,i0)', '  Matrix size: ', Ntot, 'x', Ntot
    print '(a,i0,a,i0)', '  numcb=', cfg%bands%num_cb, ' numvb=', cfg%bands%num_vb
    print '(a)', ''

    ! ----------------------------------------------------------------
    ! Build velocity matrices via simulation_setup
    ! ----------------------------------------------------------------
    ! Need a Hamiltonian built first for velocity matrix construction
    block
      type(wavevector) :: kv_zero
      kv_zero%kx = 0.0_dp
      kv_zero%ky = 0.0_dp
      kv_zero%kz = 0.0_dp
      call setup_build_H(setup, cfg, kv_zero)
    end block

    call setup_build_velocity_matrices(setup, cfg)

    print '(a)', ' Wire velocity matrices built successfully'

    ! Initialize optics accumulation
    ! ----------------------------------------------------------------
    cfg%optics%confinement = cfg%confinement
    call optics_init(oe, cfg%optics)

    ! ----------------------------------------------------------------
    ! Simpson integration weights for 1D kz integration
    ! int dkz  =>  weight = Simpson * 1 / (2*pi)
    ! ----------------------------------------------------------------
    npts = cfg%wave_vector%nsteps
    if (mod(npts, 2) == 0) then
      npts = npts - 1  ! Simpson requires odd number of points
      print '(a,i0,a)', ' Warning: optics integration uses ', npts, &
        & ' k-points (Simpson requires odd count)'
    end if
    allocate(simpson_w(npts))
    dk = cfg%wave_vector%max / real(cfg%wave_vector%nsteps - 1, kind=dp)
    do i = 1, npts
      ! Base Simpson 1/3 rule weight
      if (i == 1 .or. i == npts) then
        simpson_w(i) = dk / 3.0_dp
      else if (mod(i, 2) == 0) then
        simpson_w(i) = 4.0_dp * dk / 3.0_dp
      else
        simpson_w(i) = 2.0_dp * dk / 3.0_dp
      end if
      ! 1D BZ: weight / (2*pi) normalization
      simpson_w(i) = simpson_w(i) / (2.0_dp * pi_dp)
    end do

    ! ----------------------------------------------------------------
    ! kz sweep: build H(kz), diagonalize with FEAST, accumulate spectra
    ! ----------------------------------------------------------------
    print '(a,i0,a)', ' Wire optics kz-sweep: ', npts, ' k-points'

    do k = 1, npts
      ! Build kz value for this sweep point
      block
        real(kind=dp) :: kz_val
        type(wavevector) :: kv_sweep

        kz_val = real(k - 1, kind=dp) * dk

        ! Build wire Hamiltonian at this kz (setup dispatches internally)
        kv_sweep%kx = 0.0_dp
        kv_sweep%ky = 0.0_dp
        kv_sweep%kz = kz_val
        call setup_build_H(setup, cfg, kv_sweep)

        ! Solve eigenvalue problem
        call setup%eigen_solver%solve_sparse(setup%HT_csr_ptr, setup%eigen_cfg, eigen_res)

        if (eigen_res%nev_found == 0) then
          print '(a,a,a,i0,a)', ' WARNING: ', setup%eigen_solver%backend_name(), &
            ' found no eigenvalues at kz-point ', k
          call eigensolver_result_free(eigen_res)
          cycle
        end if

        if (eigen_res%nev_found < cfg%bands%num_cb + cfg%bands%num_vb) then
          print '(a,a,a,i0,a,i0,a,i0,a)', ' WARNING: ', &
            setup%eigen_solver%backend_name(), ' found only ', &
            eigen_res%nev_found, ' eigenvalues at kz-point ', k, &
            ' (need ', cfg%bands%num_cb + cfg%bands%num_vb, ')'
          call eigensolver_result_free(eigen_res)
          cycle
        end if

        ! Accumulate optical spectra for this kz
        call optics_accumulate(oe, &
          & eigen_res%eigenvalues(1:nev_wire), &
          & eigen_res%eigenvectors(:, 1:nev_wire), &
          & simpson_w(k), &
          & setup%vel, cfg%bands%num_cb, cfg%bands%num_vb, cfg%sc%fermi_level)

        if (cfg%optics%spontaneous_enabled) then
          call optics_accumulate_spontaneous(oe, &
            & eigen_res%eigenvalues(1:nev_wire), &
            & eigen_res%eigenvectors(:, 1:nev_wire), &
            & simpson_w(k), &
            & setup%vel, cfg%bands%num_cb, cfg%bands%num_vb, cfg%sc%fermi_level)
        end if

        if (cfg%optics%gain_enabled) then
          call compute_gain_qw(oe, &
            & eigen_res%eigenvalues(1:nev_wire), &
            & eigen_res%eigenvectors(:, 1:nev_wire), &
            & simpson_w(k), &
            & setup%vel, cfg%bands%num_cb, cfg%bands%num_vb, &
            & cfg%optics%gain_carrier_density)
        end if

        if (cfg%optics%isbt_enabled) then
          call compute_isbt_absorption(oe, &
            & eigen_res%eigenvalues(1:nev_wire), &
            & eigen_res%eigenvectors(:, 1:nev_wire), &
            & setup%vel, cfg%bands%num_cb, cfg%bands%num_vb, &
            & simpson_w(k), cfg%sc%fermi_level)
        end if

        ! Free eigensolver result for next kz
        call eigensolver_result_free(eigen_res)
      end block
    end do

    ! ----------------------------------------------------------------
    ! Finalize: apply prefactor, write output files
    ! ----------------------------------------------------------------
    call optics_apply_prefactor(oe)
    call optics_write_output(oe)
    call optics_free(oe)

    deallocate(simpson_w)
    call simulation_setup_free(setup)
    print '(a)', ' Wire optical spectra written to output/'

  case default
    print '(a,a)', 'ERROR: unknown confinement mode ', trim(cfg%confinement)
    error stop 'unknown confinement mode'

  end select

  ! ================================================================
  ! Cleanup
  ! ================================================================
  if (allocated(smallk))  deallocate(smallk)
  if (allocated(HT))      deallocate(HT)
  if (allocated(eig))     deallocate(eig)
  if (allocated(eigv))    deallocate(eigv)

end program opticalProperties
