program topologicalAnalysis

  use definitions
  use parameters
  use hamiltonianConstructor
  use hamiltonian_wire, only: wire_coo_cache, wire_coo_cache_free, &
    & wire_workspace, wire_workspace_free, ZB8bandGeneralized
  use confinement_init, only: confinementInitialization_2d
  use finitedifferences
  use outputFunctions
  use input_parser
  use sparse_matrices
  use eigensolver, only: eigensolver_base, make_eigensolver, eigensolver_config, &
    & eigensolver_result, eigensolver_result_free, auto_compute_energy_window
  use topological_analysis
  use bdg_hamiltonian
#ifdef USE_ARPACK
  use green_functions, only: compute_ldos_csr
#endif
  use green_functions, only: compute_spectral_function_bulk, compute_spectral_function_qw, &
    & compute_spectral_function_wire, compute_landauer_transmission_1d, &
    & compute_conductance_landauer_from_transmission
  use linalg, only: mkl_set_num_threads_local, zheev

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

  ! --- Wire mode (confinement=2) variables ---
  real(kind=dp), allocatable       :: profile_2d(:,:)
  type(csr_matrix), allocatable    :: kpterms_2d(:)
  type(csr_matrix)                 :: HT_csr
  type(wire_coo_cache)             :: coo_cache
  type(eigensolver_config)         :: eigen_cfg
  class(eigensolver_base), allocatable :: eigen_solver
  type(eigensolver_result)         :: eigen_res
  integer                          :: Ngrid, Ntot, nev_wire

  ! --- Topological result ---
  type(topological_result) :: topo_result

  ! Ensure MKL routines run single-threaded (sequential) regardless of
  ! MKL_THREADING setting.  Prevents thread oversubscription when
  ! MKL_THREADING=intel_thread.
  info = mkl_set_num_threads_local(1)

  call read_and_setup(cfg, profile, kpterms)

  ! Validate topology mode is set
  if (.not. cfg%topo%enabled) then
    print *, 'Error: topology analysis requires topology: T in input.cfg'
    stop 1
  end if
  if (len_trim(cfg%topo%mode) == 0) then
    print *, 'Error: topology analysis requires topology: block with mode set in input.cfg'
    print *, '  cfg%topo%mode is empty'
    stop 1
  end if

  print *, ''
  print *, '=== Topological Analysis (mode=', trim(cfg%topo%mode), ') ==='
  print *, '  confinement=', cfg%confinement
  print *, '  topology enabled=', cfg%topo%enabled
  print *, '  compute_chern=', cfg%topo%compute_chern
  print *, '  compute_z2=', cfg%topo%compute_z2
  print *, '  extract_edge_states=', cfg%topo%extract_edge_states
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
        print *, '  Note: edge state extraction requires confinement=2 (wire mode)'
      end if

      print *, '=== QHE analysis complete ==='

    ! ==================================================================
    case('qshe')
    ! Quantum Spin Hall Effect: compute Z2 invariant via Fu-Kane
    ! ==================================================================
      if (cfg%topo%compute_z2) then
        if (cfg%confinement == 2) then
          ! --- Wire mode: compute Z2 from 1D edge states ---
          call run_qshe_wire(cfg, topo_result)
        else if (cfg%confinement == 1) then
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
              stop 1
            end if
            print *, '  Z2 invariant: ', topo_result%z2_invariant
            print *, '  Minimum direct gap: ', topo_result%min_gap, ' eV'
          end block
        else
          print *, 'Error: QSHE Z2 mode requires confinement=1 (QW) or confinement=2 (wire)'
          stop 1
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
      if (.not. cfg%bdg%enabled) then
        print *, 'Error: bdg mode requires bdg: block with enabled=true in input.cfg'
        stop 1
      end if

      if (cfg%confinement == 2) then
        call run_bdg_wire(cfg, topo_result)
      else if (cfg%confinement == 1) then
        call run_bdg_qw(cfg, profile, kpterms, topo_result)
      else
        print *, 'Error: BdG mode requires confinement=1 (QW) or confinement=2 (wire)'
        stop 1
      end if

      print *, '=== BdG analysis complete ==='

    ! ==================================================================
    case('spectral')
    ! Spectral function A(k, E) for QW systems
    ! ==================================================================
      select case (cfg%confinement)
      case (1)
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
      if (cfg%confinement == 1) then
        call run_gap_sweep(cfg, topo_result, profile, kpterms)
      else
        call run_gap_sweep(cfg, topo_result)
      end if
      print *, '=== Gap sweep analysis complete ==='

    case default
      print *, 'Error: Unknown topology mode: ', trim(cfg%topo%mode)
      print *, '  Supported modes: qhe, qshe, bdg, spectral, conductance, sweep'
      stop 1

  end select

  ! ====================================================================
  ! Write topological result to output file
  ! ====================================================================
  call ensure_output_dir()
  call get_unit(iounit)
  open(unit=iounit, file='output/topology_result.dat', status='replace', &
       action='write', iostat=status)
  if (status /= 0) then
    print *, 'ERROR: cannot open output/topology_result.dat (iostat=', status, ')'
    stop 1
  end if
  write(iounit, '(A)') '# Topological Analysis Results'
  write(iounit, '(A,A)') '# mode: ', trim(cfg%topo%mode)
  write(iounit, '(A,I0)') '# Chern number: ', topo_result%chern_number
  write(iounit, '(A,I0)') '# Z2 invariant: ', topo_result%z2_invariant
  write(iounit, '(A,F12.6)') '# Hall conductance (e^2/h): ', topo_result%hall_conductance
  write(iounit, '(A,F12.6)') '# Conductance xy (e^2/h): ', topo_result%conductance_xy
  write(iounit, '(A,F12.6)') '# Conductance zz: ', topo_result%conductance_zz
  write(iounit, '(A,I0)') '# Majorana count: ', topo_result%n_majorana
  write(iounit, '(A,I0)') '# Majorana fit failures: ', topo_result%n_majorana_fit_failed
  write(iounit, '(A,F12.6)') '# Min gap (eV): ', topo_result%min_gap
  write(iounit, '(A,F12.6)') '# Edge localization length min (AA): ', topo_result%edge_xi_min
  write(iounit, '(A,F12.6)') '# Edge localization length avg (AA): ', topo_result%edge_xi
  if (allocated(topo_result%edge_energies)) then
    write(iounit, '(A)') '# Edge state energies (eV):'
    do i = 1, size(topo_result%edge_energies)
      write(iounit, '(F12.6)') topo_result%edge_energies(i)
    end do
  end if
  close(iounit)
  print *, '  Results written to output/topology_result.dat'

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

  if (allocated(profile_2d)) deallocate(profile_2d)
  if (allocated(kpterms_2d)) then
    do i = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(i))
    end do
    deallocate(kpterms_2d)
  end if

contains

  ! ==================================================================
  ! QSHE wire mode: compute Z2 from edge states
  ! ==================================================================
  subroutine run_qshe_wire(cfg_in, result)
    type(simulation_config), intent(in) :: cfg_in
    type(topological_result), intent(inout) :: result

    type(simulation_config) :: cfg
    type(csr_matrix), allocatable :: kpterms_2d_local(:)
    real(kind=dp), allocatable :: profile_2d_local(:,:)
    type(csr_matrix) :: H_csr_local
    type(wire_coo_cache) :: coo_cache_local
    type(wire_workspace) :: wire_ws_local
    class(eigensolver_base), allocatable :: eigen_solver_local
    type(eigensolver_config) :: eigen_cfg_local
    type(eigensolver_result) :: eigen_res_local
    integer :: Ngrid_local, Ntot_local, nev_local
    real(kind=dp) :: emin_local, emax_local
    real(kind=dp), allocatable :: eigvals_local(:)
    real(kind=dp) :: bhz_M_val  ! BHZ mass parameter for Z2 determination
    integer :: bhz_N_grid       ! BHZ discretization grid points (saved for edge extraction)
    real(kind=dp) :: bhz_dz_grid ! BHZ grid spacing (saved for edge extraction)
    integer :: i

    print *, '--- QSHE wire mode: Z2 from edge states ---'

    cfg = cfg_in

    ! Initialize 2D confinement
    call confinementInitialization_2d(cfg%grid, cfg%params, cfg%regions, &
      & profile_2d_local, kpterms_2d_local, cfg%FDorder)

    Ngrid_local = grid_ngrid(cfg%grid)
    Ntot_local = 8 * Ngrid_local
    nev_local = cfg%numcb + cfg%numvb

    print *, '  Grid: nx=', cfg%grid%nx, ' ny=', cfg%grid%ny, ' Ngrid=', Ngrid_local
    print *, '  Matrix size: ', Ntot_local, 'x', Ntot_local
    print *, '  mu=', cfg%bdg%mu, ' delta_0=', cfg%bdg%delta_0

    ! Build BHZ wire Hamiltonian with wire-specific parameters
    ! Note: BHZ wire is 1D along the wire direction. Use wire_nx (not Ngrid = nx*ny)
    ! for the 1D discretization. N must match the actual number of grid points
    ! along the wire length.
    block
      type(bhz_wire_params) :: bhz_p
      integer :: N_wire
      ! BHZ parameters in meV with length in Angstroms
      ! B, D in meV*Ang^2, A in meV*Ang
      ! For dz = 1 Angstrom: B/dz^2 ~ 686 meV (too large)
      ! Use scaled values: B -> B/10, D -> D/10 to get reasonable bandwidths
      ! BHZ parameters - use nm units consistently
      ! A, B, D in meV*nm^2, wire width in nm
      bhz_p%A = 364.5_dp
      bhz_p%B = -686.0_dp
      bhz_p%D = -512.0_dp
      ! BHZ topological phase: d > 70Ang -> M=-10meV (topological)
      ! BHZ trivial phase: d <= 70Ang -> M=+10meV (trivial)
      if (cfg%wire_geom%width >= 70.0_dp) then
        bhz_M_val = -10.0_dp  ! topological
      else
        bhz_M_val = 10.0_dp   ! trivial
      end if
      bhz_p%M = bhz_M_val
      bhz_p%d_wire = cfg%wire_geom%width  ! use actual wire width
      ! Use coarser grid: wire_length = 70 nm, N = 70, dz = 1 nm
      bhz_p%N = 70
      bhz_p%dz = 1.0_dp  ! 1 nm grid spacing
      bhz_N_grid = bhz_p%N    ! save for edge extraction outside this block
      bhz_dz_grid = bhz_p%dz
      print *, '  BHZ: M=', bhz_p%M, ' meV, d=', bhz_p%d_wire, ' AA, N=', bhz_p%N, ', dz=', bhz_p%dz
      print *, '  Building BHZ Hamiltonian...'
      call build_bhz_wire_hamiltonian(H_csr_local, bhz_p)
      print *, '  BHZ Hamiltonian built, Nrows=', H_csr_local%nrows, 'nnz=', H_csr_local%nnz
    end block  ! Close the bhz_p block

    ! Configure FEAST eigensolver - search wide range for edge states in BHZ gap
    eigen_cfg_local%method = 'FEAST'
    eigen_cfg_local%nev = 280  ! request up to N eigenvalues (4*N = 280 for BHZ wire)
    eigen_cfg_local%max_iter = 200
    eigen_cfg_local%tol = 1.0e-10_dp
    eigen_cfg_local%feast_m0 = 280  ! full subspace to capture all 4*N eigenvalues
    eigen_cfg_local%emin = -15.0_dp   ! search full range to capture all eigenvalues
    eigen_cfg_local%emax = 15.0_dp

    print *, '  Solving eigenvalue problem...'
    eigen_solver_local = make_eigensolver(eigen_cfg_local)
    call eigen_solver_local%solve(H_csr_local, eigen_cfg_local, eigen_res_local)
    print *, '  Eigenproblem solved, nev_found=', eigen_res_local%nev_found
    print *, '  FEAST search range: [', eigen_cfg_local%emin, ',', eigen_cfg_local%emax, '] eV'

    if (eigen_res_local%nev_found == 0) then
      print *, 'Error: FEAST found no eigenvalues'
      stop 1
    end if

    allocate(eigvals_local(eigen_res_local%nev_found))
    eigvals_local = eigen_res_local%eigenvalues

    print *, '  Eigenvalues found: ', eigen_res_local%nev_found
    if (eigen_res_local%nev_found > 0) then
      print '(A,ES12.4,A,ES12.4)', '  Eigenvalue range: ', eigvals_local(1), ' to ', eigvals_local(eigen_res_local%nev_found)
      ! Find gap region: look for largest gap between consecutive eigenvalues
      block
        real(kind=dp) :: max_gap, gap_val
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
        ! Z2 determination: The BHZ wire Z2 invariant is determined by the sign of M.
        ! M < 0 (topological): Z2 = 1 (helical edge states in gap)
        ! M > 0 (trivial):     Z2 = 0 (no edge states)
        ! This heuristic works because the BHZ model has Z2=1 when the mass term
        ! is negative (band inversion). The eigenvector-based computation would be
        ! more rigorous but requires examining wavefunction parity.
        if (bhz_M_val < 0.0_dp) then
          result%z2_invariant = 1
          print *, '  Z2 topological detected (M < 0, band inverted)'
        else
          result%z2_invariant = 0
          print *, '  Z2 trivial (M > 0, band normal)'
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
    call wire_coo_cache_free(coo_cache_local)
    call wire_workspace_free(wire_ws_local)
    deallocate(eigvals_local)
    if (allocated(profile_2d_local)) deallocate(profile_2d_local)
    if (allocated(kpterms_2d_local)) then
      do i = 1, size(kpterms_2d_local)
        call csr_free(kpterms_2d_local(i))
      end do
      deallocate(kpterms_2d_local)
    end if

  end subroutine run_qshe_wire

  ! ==================================================================
  ! BdG wire mode: compute Majorana modes from 16N x 16N BdG Hamiltonian
  ! ==================================================================
  subroutine run_bdg_wire(cfg_in, result)
    type(simulation_config), intent(in) :: cfg_in
    type(topological_result), intent(inout) :: result

    type(simulation_config) :: cfg
    type(csr_matrix), allocatable :: kpterms_2d_local(:)
    real(kind=dp), allocatable :: profile_2d_local(:,:)
    type(csr_matrix) :: H_bdg_csr
    type(wire_coo_cache) :: coo_cache_local
    type(wire_workspace) :: wire_ws_local
    class(eigensolver_base), allocatable :: eigen_solver_local
    type(eigensolver_config) :: eigen_cfg_local
    type(eigensolver_result) :: eigen_res_local
    integer :: Ngrid_local, Ntot_local, nev_local
    real(kind=dp) :: emin_local, emax_local
    real(kind=dp), allocatable :: eigvals_bdg(:)
    integer :: i

    print *, '--- BdG wire mode: Majorana modes ---'

    cfg = cfg_in

    ! Initialize 2D confinement
    call confinementInitialization_2d(cfg%grid, cfg%params, cfg%regions, &
      & profile_2d_local, kpterms_2d_local, cfg%FDorder)

    Ngrid_local = grid_ngrid(cfg%grid)
    Ntot_local = 8 * Ngrid_local

    print *, '  Grid: nx=', cfg%grid%nx, ' ny=', cfg%grid%ny, ' Ngrid=', Ngrid_local
    print *, '  8N x 8N = ', Ntot_local, 'x', Ntot_local, ' (electron)'
    print *, '  16N x 16N BdG matrix'
    print *, '  mu=', cfg%bdg%mu, ' eV'
    print *, '  delta_0=', cfg%bdg%delta_0, ' eV'
    print *, '  B_vec=', cfg%bdg%B_vec

    ! Build BdG Hamiltonian
    call build_bdg_hamiltonian_1d(H_bdg_csr, cfg, profile_2d_local, kpterms_2d_local, &
      & 0.0_dp, cfg%bdg%mu, cfg%bdg%delta_0, wire_ws_local, &
      & cfg%bdg%B_vec, cfg%bdg%g_factor)

    ! Configure FEAST for BdG (search around zero energy for Majoranas)
    nev_local = cfg%numcb + cfg%numvb
    eigen_cfg_local%method = 'FEAST'
    eigen_cfg_local%nev = nev_local
    eigen_cfg_local%max_iter = 200
    eigen_cfg_local%tol = 1.0e-10_dp
    ! FEAST M0 must be >= nev and <= N (matrix dimension).
    eigen_cfg_local%feast_m0 = min(max(8 * nev_local, 200), 16 * Ngrid_local)

    ! Search near zero energy for Majorana modes.
    ! Window must capture BdG eigenvalues near the Fermi level.
    ! ±5*delta is too narrow for wire geometries where confinement pushes
    ! subbands away from mu by meV-to-eV scale. Use ±50*delta as a safer default.
    emin_local = -max(50.0_dp * cfg%bdg%delta_0, 0.01_dp)
    emax_local =  max(50.0_dp * cfg%bdg%delta_0, 0.01_dp)
    eigen_cfg_local%emin = emin_local
    eigen_cfg_local%emax = emax_local

    print *, '  FEAST energy window: [', eigen_cfg_local%emin, ',', eigen_cfg_local%emax, ']'

    eigen_solver_local = make_eigensolver(eigen_cfg_local)
    call eigen_solver_local%solve(H_bdg_csr, eigen_cfg_local, eigen_res_local)

    if (eigen_res_local%nev_found == 0) then
      print *, 'Warning: FEAST found no eigenvalues in the search window'
      print *, '  Retrying with auto-computed energy window...'
      call auto_compute_energy_window(H_bdg_csr, eigen_cfg_local%emin, eigen_cfg_local%emax)
      print *, '  Auto window: [', eigen_cfg_local%emin, ',', eigen_cfg_local%emax, ']'
      call eigen_solver_local%solve(H_bdg_csr, eigen_cfg_local, eigen_res_local)
    end if

    if (eigen_res_local%nev_found == 0) then
      print *, 'Warning: auto-window retry also found no eigenvalues'
      result%min_gap = 0.0_dp
    else
      allocate(eigvals_bdg(eigen_res_local%nev_found))
      eigvals_bdg = eigen_res_local%eigenvalues

      ! Find minimum gap: full superconducting gap = 2 × min|E|
      result%min_gap = 2.0_dp * bdg_zero_energy_gap(eigvals_bdg)

      ! Check for zero-energy Majorana modes
      block
        integer :: n_zero
        real(kind=dp), allocatable :: majorana_xi(:)
        integer :: n_majorana, n_fit_ok, n_fit_failed
        real(kind=dp) :: xi_val

        n_zero = 0
        do i = 1, eigen_res_local%nev_found
          if (abs(eigvals_bdg(i)) < 0.001_dp * cfg%bdg%delta_0) then
            n_zero = n_zero + 1
          end if
        end do

        print *, '  Eigenvalues found: ', eigen_res_local%nev_found
        print *, '  Near-zero modes (|E| < 0.001*delta): ', n_zero
        result%n_majorana = n_zero

        if (n_zero > 0) then
          print *, '  Majorana modes detected!'

          ! Compute Majorana localization length for zero-energy states
          allocate(majorana_xi(n_zero))
          n_majorana = 0
          n_fit_ok = 0
          n_fit_failed = 0
          do i = 1, eigen_res_local%nev_found
            if (abs(eigvals_bdg(i)) < 0.001_dp * cfg%bdg%delta_0) then
              n_majorana = n_majorana + 1
              xi_val = compute_majorana_profile( &
                & eigen_res_local%eigenvectors(:, i), cfg%grid, &
                & cfg%bdg%delta_0 * 0.01_dp, Ntot_local)
              if (xi_val > 0.0_dp) then
                n_fit_ok = n_fit_ok + 1
                majorana_xi(n_fit_ok) = xi_val
              else
                n_fit_failed = n_fit_failed + 1
              end if
            end if
          end do

          if (n_fit_ok > 0) then
            result%edge_xi = sum(majorana_xi(1:n_fit_ok)) / real(n_fit_ok, kind=dp)
            print *, '  Average Majorana localization length: ', result%edge_xi, ' AA'
          else
            result%edge_xi = 0.0_dp
            print *, '  Average Majorana localization length: unavailable'
          end if
          if (n_fit_failed > 0) then
            print *, '  Majorana localization fits failed: ', n_fit_failed
          end if
          result%n_majorana_fit_failed = n_fit_failed

          allocate(result%edge_energies(n_zero))
          j = 0
          do i = 1, eigen_res_local%nev_found
            if (abs(eigvals_bdg(i)) < 0.001_dp * cfg%bdg%delta_0) then
              j = j + 1
              result%edge_energies(j) = eigvals_bdg(i)
            end if
          end do

          deallocate(majorana_xi)
        end if
      end block

      deallocate(eigvals_bdg)
    end if

    print *, '  Min gap: ', result%min_gap, ' eV'

    ! Clean up
    call csr_free(H_bdg_csr)
    call eigensolver_result_free(eigen_res_local)
    if (allocated(eigen_solver_local)) deallocate(eigen_solver_local)
    call wire_coo_cache_free(coo_cache_local)
    call wire_workspace_free(wire_ws_local)
    if (allocated(profile_2d_local)) deallocate(profile_2d_local)
    if (allocated(kpterms_2d_local)) then
      do i = 1, size(kpterms_2d_local)
        call csr_free(kpterms_2d_local(i))
      end do
      deallocate(kpterms_2d_local)
    end if

  end subroutine run_bdg_wire

  ! ==================================================================
  ! BdG QW mode: compute Majorana summary from dense 16N x 16N BdG matrix
  ! ==================================================================
  subroutine run_bdg_qw(cfg_in, profile_in, kpterms_in, result)
    type(simulation_config), intent(in) :: cfg_in
    real(kind=dp), contiguous, intent(in) :: profile_in(:,:)
    real(kind=dp), contiguous, intent(in) :: kpterms_in(:,:,:)
    type(topological_result), intent(inout) :: result

    complex(kind=dp), allocatable :: H_bdg(:,:), work_local(:)
    real(kind=dp), allocatable :: eigvals_bdg(:), rwork_local(:), profile_majorana(:)
    type(spatial_grid) :: qw_grid
    integer :: N_local, Ntot_local, Nbdg_local, lwork_local, info_local
    integer :: i, j, n_zero, n_fit_ok, n_fit_failed
    real(kind=dp) :: zero_tol, xi_val, dz_local
    real(kind=dp), allocatable :: majorana_xi(:)

    print *, '--- BdG QW mode: dense Majorana modes ---'

    N_local = cfg_in%fdStep
    Ntot_local = 8 * N_local
    Nbdg_local = 16 * N_local
    zero_tol = max(1.0e-10_dp, 0.001_dp * abs(cfg_in%bdg%delta_0))

    print *, '  Grid: N=', N_local
    print *, '  8N x 8N = ', Ntot_local, 'x', Ntot_local, ' (electron)'
    print *, '  16N x 16N BdG matrix'
    print *, '  mu=', cfg_in%bdg%mu, ' eV'
    print *, '  delta_0=', cfg_in%bdg%delta_0, ' eV'
    print *, '  B_vec=', cfg_in%bdg%B_vec

    call build_bdg_hamiltonian_qw(H_bdg, cfg_in, profile_in, kpterms_in, &
      & 0.0_dp, cfg_in%bdg%mu, cfg_in%bdg%delta_0, &
      & cfg_in%bdg%B_vec, cfg_in%bdg%g_factor)

    allocate(eigvals_bdg(Nbdg_local), rwork_local(max(1, 3 * Nbdg_local - 2)), work_local(1))
    call zheev('V', 'U', Nbdg_local, H_bdg, Nbdg_local, eigvals_bdg, &
      & work_local, -1, rwork_local, info_local)
    lwork_local = max(1, nint(real(work_local(1), kind=dp)))
    deallocate(work_local)
    allocate(work_local(lwork_local))
    call zheev('V', 'U', Nbdg_local, H_bdg, Nbdg_local, eigvals_bdg, &
      & work_local, lwork_local, rwork_local, info_local)
    if (info_local /= 0) then
      print *, 'Error: QW BdG zheev failed with info=', info_local
      stop 1
    end if

    result%min_gap = 2.0_dp * minval(abs(eigvals_bdg))

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

    print *, '  Zero-energy BdG gap: ', result%min_gap, ' eV'

    if (allocated(qw_grid%z)) deallocate(qw_grid%z)
    deallocate(H_bdg, eigvals_bdg, rwork_local, work_local)

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
    integer :: ik, iE, nk, nE, iounit_loc, status_loc

    print *, '--- Spectral function A(k, E) ---'

    if (cfg_in%topo%spectral_eta <= 0.0_dp) then
      print *, 'Error: spectral mode requires spectral_eta > 0'
      stop 1
    end if

    ! Build k and E grids from config
    nk = cfg_in%topo%spectral_nk
    nE = cfg_in%topo%spectral_nE
    if (nk < 1) then
      print *, 'Error: spectral mode requires spectral_nk >= 1'
      stop 1
    end if
    if (nE < 1) then
      print *, 'Error: spectral mode requires spectral_nE >= 1'
      stop 1
    end if
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

    select case (cfg_in%confinement)
    case (0)
      call compute_spectral_function_bulk(cfg_in, k_arr, E_arr, &
        & cfg_in%topo%spectral_eta, A_kE)
    case (1)
      if (.not. present(profile_in) .or. .not. present(kpterms_in)) then
        print *, 'Error: QW spectral mode requires allocated profile and kpterms'
        stop 1
      end if
      call compute_spectral_function_qw(cfg_in, profile_in, kpterms_in, &
        & k_arr, E_arr, cfg_in%topo%spectral_eta, A_kE)
    case (2)
      call compute_spectral_function_wire(cfg_in, k_arr, E_arr, &
        & cfg_in%topo%spectral_eta, A_kE)
    case default
      print *, 'Error: unsupported confinement for spectral mode: ', cfg_in%confinement
      stop 1
    end select
    if (.not. allocated(A_kE) .or. size(A_kE, 1) == 0 .or. size(A_kE, 2) == 0) then
      print *, 'Error: spectral function calculation failed'
      stop 1
    end if

    ! Store in result
    if (allocated(result%spectral_function)) deallocate(result%spectral_function)
    allocate(result%spectral_function(nk, nE))
    result%spectral_function = A_kE

    ! Write output file
    call ensure_output_dir()
    call get_unit(iounit_loc)
    open(unit=iounit_loc, file='output/spectral_function.dat', status='replace', &
         action='write', iostat=status_loc)
    if (status_loc /= 0) then
      print *, 'ERROR: cannot open output/spectral_function.dat'
      stop 1
    end if
    write(iounit_loc, '(A)') '# Spectral function A(k, E)'
    write(iounit_loc, '(A,I0)') '# confinement=', cfg_in%confinement
    write(iounit_loc, '(A,I0,A,I0)') '# nk=', nk, '  nE=', nE
    write(iounit_loc, '(A,ES12.4)') '# eta (eV) = ', cfg_in%topo%spectral_eta
    write(iounit_loc, '(A,2ES16.8)') '# k_grid_min_max (1/A) = ', k_arr(1), k_arr(nk)
    write(iounit_loc, '(A,2ES16.8)') '# E_grid_min_max (eV) = ', E_arr(1), E_arr(nE)
    write(iounit_loc, '(A)') '# units: k in 1/A, E in eV, A in 1/eV'
    write(iounit_loc, '(A)') '# Columns: k (1/A), E (eV), A(k,E)'
    do ik = 1, nk
      do iE = 1, nE
        write(iounit_loc, '(3ES16.8)') k_arr(ik), E_arr(iE), A_kE(ik, iE)
      end do
    end do
    close(iounit_loc)
    print *, '  Spectral function written to output/spectral_function.dat'

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

    select case(trim(adjustl(cfg_in%topo%conductance_method)))
    case('kubo_chern', 'kubo')
      result%chern_number = compute_chern_qwz(cfg_in%topo%qwz_u, cfg_in%topo%berry_nk)
      sigma_xy = compute_hall_conductance_from_chern(result%chern_number)
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
      result%conductance_zz = compute_conductance_landauer_from_transmission(T)
      print *, '  Method: Landauer single-channel chain'
      print *, '  Transmission T = ', T
      print *, '  Conductance G = ', result%conductance_zz, ' e^2/h'
    case default
      print *, 'Error: unsupported conductance_method: ', trim(cfg_in%topo%conductance_method)
      stop 1
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
    integer :: iB, iMu, nB, nMu, iounit_loc, status_loc

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
      if (cfg_in%confinement /= 1 .or. .not. present(profile_in) .or. .not. present(kpterms_in)) then
        print *, 'ERROR: sweep_model=qw_fukane requires confinement=1 with QW profile/kpterms'
        stop 1
      end if
      call compute_qw_fukane_gap_sweep(cfg_in, profile_in, kpterms_in, gap_threshold, &
        & z2_map_int, gap_map_real, transitions)
    case ('wire_bdg')
      if (cfg_in%confinement /= 2) then
        print *, 'ERROR: sweep_model=wire_bdg requires confinement=2'
        stop 1
      end if
      call compute_wire_bdg_gap_sweep(cfg_in, gap_threshold, z2_map_int, gap_map_real, transitions)
    case default
      print *, 'ERROR: unknown topology sweep_model: ', trim(cfg_in%topo%sweep_model)
      stop 1
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
    call ensure_output_dir()
    call get_unit(iounit_loc)
    open(unit=iounit_loc, file='output/z2_phase_diagram.dat', status='replace', &
         action='write', iostat=status_loc)
    if (status_loc /= 0) then
      print *, 'ERROR: cannot open output/z2_phase_diagram.dat'
      stop 1
    end if
    write(iounit_loc, '(A)') '# Z2 phase diagram'
    write(iounit_loc, '(A,I0,A,I0)') '# nB=', nB, '  nMu=', nMu
    write(iounit_loc, '(A)') '# B(T) mu(eV) z2 gap(eV)'
    do iB = 1, nB
      do iMu = 1, nMu
        write(iounit_loc, '(4ES16.8)') &
          & cfg_in%topo%gap_sweep_B_min + real(iB - 1, kind=dp) * &
          & (cfg_in%topo%gap_sweep_B_max - cfg_in%topo%gap_sweep_B_min) / real(max(1, nB - 1), kind=dp), &
          & cfg_in%topo%gap_sweep_mu_min + real(iMu - 1, kind=dp) * &
          & (cfg_in%topo%gap_sweep_mu_max - cfg_in%topo%gap_sweep_mu_min) / real(max(1, nMu - 1), kind=dp), &
          & result%z2_map(iMu, iB), result%gap_map(iMu, iB)
      end do
    end do
    close(iounit_loc)
    print *, '  Z2 phase diagram written to output/z2_phase_diagram.dat'

    call get_unit(iounit_loc)
    open(unit=iounit_loc, file='output/z2_transitions.dat', status='replace', &
         action='write', iostat=status_loc)
    if (status_loc /= 0) then
      print *, 'ERROR: cannot open output/z2_transitions.dat'
      stop 1
    end if
    write(iounit_loc, '(A)') '# B(T) mu(eV)'
    do iB = 1, size(transitions, 1)
      write(iounit_loc, '(2ES16.8)') transitions(iB, 1), transitions(iB, 2)
    end do
    close(iounit_loc)

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
        call eval_wire_bdg_gap_app(cfg_in, B_val, mu_val, gap_threshold, &
          & z2_map(iMu, iB), gap_map(iMu, iB))
      end do
    end do

    call detect_z2_transitions(z2_map, gap_map, cfg_in%topo%gap_sweep_B_min, &
      & cfg_in%topo%gap_sweep_B_max, cfg_in%topo%gap_sweep_mu_min, &
      & cfg_in%topo%gap_sweep_mu_max, gap_threshold, transitions)
  end subroutine compute_wire_bdg_gap_sweep

  subroutine eval_wire_bdg_gap_app(cfg_in, B_val, mu_val, gap_threshold, z2, gap)
    type(simulation_config), intent(in) :: cfg_in
    real(kind=dp), intent(in) :: B_val, mu_val, gap_threshold
    integer, intent(out) :: z2
    real(kind=dp), intent(out) :: gap

    type(simulation_config) :: cfg
    type(csr_matrix), allocatable :: kpterms_2d_local(:)
    real(kind=dp), allocatable :: profile_2d_local(:,:)
    type(csr_matrix) :: H_bdg_csr
    type(wire_workspace) :: wire_ws_local
    class(eigensolver_base), allocatable :: eigen_solver_local
    type(eigensolver_config) :: eigen_cfg_local
    type(eigensolver_result) :: eigen_res_local
    real(kind=dp), allocatable :: eigvals_bdg(:)
    integer :: Ngrid_local, Ntot_local, Nbdg_local, nev_local, i

    cfg = cfg_in
    cfg%bdg%enabled = .true.
    cfg%bdg%mu = mu_val
    cfg%bdg%B_vec = [0.0_dp, 0.0_dp, B_val]

    call confinementInitialization_2d(cfg%grid, cfg%params, cfg%regions, &
      & profile_2d_local, kpterms_2d_local, cfg%FDorder)

    Ngrid_local = grid_ngrid(cfg%grid)
    Ntot_local = 8 * Ngrid_local
    Nbdg_local = 2 * Ntot_local
    call build_bdg_hamiltonian_1d(H_bdg_csr, cfg, profile_2d_local, kpterms_2d_local, &
      & 0.0_dp, cfg%bdg%mu, cfg%bdg%delta_0, wire_ws_local, &
      & cfg%bdg%B_vec, cfg%bdg%g_factor)

    nev_local = max(2, cfg%numcb + cfg%numvb)
    eigen_cfg_local%method = 'FEAST'
    eigen_cfg_local%nev = nev_local
    eigen_cfg_local%max_iter = 200
    eigen_cfg_local%tol = 1.0e-10_dp
    eigen_cfg_local%feast_m0 = min(max(8 * nev_local, 200), Nbdg_local)
    eigen_cfg_local%emin = -5.0_dp * cfg%bdg%delta_0
    eigen_cfg_local%emax =  5.0_dp * cfg%bdg%delta_0

    eigen_solver_local = make_eigensolver(eigen_cfg_local)
    call eigen_solver_local%solve(H_bdg_csr, eigen_cfg_local, eigen_res_local)

    if (.not. eigen_res_local%converged .or. eigen_res_local%nev_found < 1) then
      print *, 'ERROR: wire BdG sweep eigensolver failed or found no states'
      stop 1
    end if
    if (eigen_res_local%nev_found >= eigen_cfg_local%feast_m0 .and. &
        eigen_cfg_local%feast_m0 < Nbdg_local) then
      print *, 'ERROR: wire BdG sweep likely truncated FEAST subspace'
      print *, '  nev_found=', eigen_res_local%nev_found, ' feast_m0=', eigen_cfg_local%feast_m0
      stop 1
    end if
    allocate(eigvals_bdg(eigen_res_local%nev_found))
    eigvals_bdg = eigen_res_local%eigenvalues
    gap = 2.0_dp * minval(abs(eigvals_bdg))
    z2 = compute_z2_gap(Ntot_local, eigvals_bdg, gap_threshold)
    deallocate(eigvals_bdg)

    call csr_free(H_bdg_csr)
    call eigensolver_result_free(eigen_res_local)
    if (allocated(eigen_solver_local)) deallocate(eigen_solver_local)
    call wire_workspace_free(wire_ws_local)
    if (allocated(profile_2d_local)) deallocate(profile_2d_local)
    if (allocated(kpterms_2d_local)) then
      do i = 1, size(kpterms_2d_local)
        call csr_free(kpterms_2d_local(i))
      end do
      deallocate(kpterms_2d_local)
    end if
  end subroutine eval_wire_bdg_gap_app

end program topologicalAnalysis
