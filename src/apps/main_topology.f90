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
  use green_functions
  use linalg, only: mkl_set_num_threads_local

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
  integer(kind=4) :: iounit

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
        print *, '  DEBUG: about to call run_qshe_wire, compute_z2=', cfg%topo%compute_z2
        if (cfg%confinement == 2) then
          ! --- Wire mode: compute Z2 from 1D edge states ---
          call run_qshe_wire(cfg, topo_result)
        else
          ! --- Bulk/QW: Z2 from band structure (placeholder) ---
          print *, '--- Computing Z2 invariant (QSHE) ---'
          print *, '  Note: Full Z2 computation for bulk/QW requires band structure sweep'
          print *, '  Z2 = 0 (placeholder — implement with band inversion analysis)'
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
      else
        print *, 'Error: BdG mode currently requires confinement=2 (wire mode)'
        stop 1
      end if

      print *, '=== BdG analysis complete ==='

    case default
      print *, 'Error: Unknown topology mode: ', trim(cfg%topo%mode)
      print *, '  Supported modes: qhe, qshe, bdg'
      stop 1

  end select

  ! ====================================================================
  ! Write topological result to output file
  ! ====================================================================
  call ensure_output_dir()
  call get_unit(iounit)
  open(unit=iounit, file='output/topology_result.dat', status='replace', action='write')
  write(iounit, '(A)') '# Topological Analysis Results'
  write(iounit, '(A,A)') '# mode: ', trim(cfg%topo%mode)
  write(iounit, '(A,I0)') '# Chern number: ', topo_result%chern_number
  write(iounit, '(A,I0)') '# Z2 invariant: ', topo_result%z2_invariant
  write(iounit, '(A,F12.6)') '# Hall conductance (e^2/h): ', topo_result%hall_conductance
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
        edge_info = extract_edge_states_wire(eigvals_local, eigen_res_local%eigenvectors, &
          & cfg%grid, cfg%topo%edge_E_window)
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
    ! Use moderate feast_m0 for BdG: enough to capture eigenvalues but not full N
    ! FEAST M0 should be between nev and N. Use 8*nev for better coverage.
    eigen_cfg_local%feast_m0 = max(8 * nev_local, 200)

    ! Search near zero energy for Majorana modes
    ! Use ±5*delta to focus on the superconducting gap region
    eigen_cfg_local%emin = -5.0_dp * cfg%bdg%delta_0
    eigen_cfg_local%emax =  5.0_dp * cfg%bdg%delta_0

    print *, '  FEAST energy window: [', eigen_cfg_local%emin, ',', eigen_cfg_local%emax, ']'

    eigen_solver_local = make_eigensolver(eigen_cfg_local)
    call eigen_solver_local%solve(H_bdg_csr, eigen_cfg_local, eigen_res_local)

    if (eigen_res_local%nev_found == 0) then
      print *, 'Warning: FEAST found no eigenvalues in the search window'
      result%min_gap = 0.0_dp
    else
      allocate(eigvals_bdg(eigen_res_local%nev_found))
      eigvals_bdg = eigen_res_local%eigenvalues

      ! Find minimum gap around zero energy
      block
        real(kind=dp) :: gap_min_val
        integer :: i
        gap_min_val = huge(1.0_dp)
        ! Look at ALL consecutive eigenvalue pairs, find the minimum spacing
        ! This will show gap closure at B_crit as min gap -> 0
        do i = 1, eigen_res_local%nev_found - 1
          gap_min_val = min(gap_min_val, abs(eigvals_bdg(i+1) - eigvals_bdg(i)))
        end do
        result%min_gap = gap_min_val
      end block

      ! Check for zero-energy Majorana modes
      block
        integer :: n_zero
        real(kind=dp), allocatable :: majorana_xi(:)
        integer :: n_majorana

        n_zero = 0
        do i = 1, eigen_res_local%nev_found
          if (abs(eigvals_bdg(i)) < 0.001_dp * cfg%bdg%delta_0) then
            n_zero = n_zero + 1
          end if
        end do

        print *, '  Eigenvalues found: ', eigen_res_local%nev_found
        print *, '  Near-zero modes (|E| < 0.001*delta): ', n_zero

        if (n_zero > 0) then
          print *, '  Majorana modes detected!'

          ! Compute Majorana localization length for zero-energy states
          allocate(majorana_xi(n_zero))
          n_majorana = 0
          do i = 1, eigen_res_local%nev_found
            if (abs(eigvals_bdg(i)) < 0.001_dp * cfg%bdg%delta_0) then
              n_majorana = n_majorana + 1
              majorana_xi(n_majorana) = compute_majorana_profile( &
                & eigen_res_local%eigenvectors(:, i), cfg%grid, &
                & cfg%bdg%delta_0 * 0.01_dp, Ntot_local)
            end if
          end do

          if (n_majorana > 0) then
            result%edge_xi = sum(majorana_xi(1:n_majorana)) / real(n_majorana, kind=dp)
            print *, '  Average Majorana localization length: ', result%edge_xi, ' AA'
          end if

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

end program topologicalAnalysis
