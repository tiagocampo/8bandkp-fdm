program kpfdm

  use definitions
  use parameters
  use hamiltonianConstructor
  use finitedifferences
  use OMP_lib
  use outputFunctions
  use input_parser
  use sc_loop
  use sparse_matrices
  use eigensolver
  use strain_solver
  use linalg, only: zheevx, mkl_set_num_threads_local, ilaenv, dlamch

  implicit none

  ! Shared configuration from input_parser
  type(simulation_config) :: cfg
  real(kind=dp), allocatable, dimension(:,:) :: profile
  real(kind=dp), allocatable, dimension(:,:,:) :: kpterms

  ! wave vector
  type(wavevector), allocatable, dimension(:) :: smallk

  ! iteration consts
  integer :: i, k, ii, jj

  ! hamiltonian and LAPACK/BLAS
  integer :: info, NB, lwork, N, M, il, iuu
  real(kind=dp) :: abstol, vl, vu
  real(kind=dp), allocatable :: eig(:,:), rwork(:)
  complex(kind=dp), allocatable :: work(:)
  complex(kind=dp), allocatable, dimension(:,:,:) :: eigv
  complex(kind=dp), allocatable, dimension(:,:) :: HT, HTmp
  integer, allocatable :: iwork(:), ifail(:)

  ! file handling
  integer(kind=4) :: iounit

  ! --- Wire mode variables ---
  real(kind=dp), allocatable       :: profile_2d(:,:)
  type(csr_matrix), allocatable    :: kpterms_2d(:)
  type(csr_matrix)                 :: HT_csr
  type(wire_coo_cache)             :: coo_cache
  type(eigensolver_config)         :: eigen_cfg
  type(eigensolver_result)         :: eigen_res
  real(kind=dp), allocatable       :: eig_wire(:,:)
  integer                          :: Ngrid, Ntot, nev_wire, max_nev_found

  ! --- QW SC charge density output ---
  real(kind=dp), allocatable       :: sc_ne_qw(:), sc_nh_qw(:)

  ! Shared setup: read input, initialize materials, confinement, external field
  call read_and_setup(cfg, profile, kpterms)

  ! Build wave vector array
  allocate(smallk(cfg%waveVectorStep))
  smallk%kx = 0
  smallk%ky = 0
  smallk%kz = 0
  select case(cfg%waveVector)

    case ("kx")
      smallk%kx = [ ((i-1)*cfg%waveVectorMax/(cfg%waveVectorStep-1), i=1,cfg%waveVectorStep) ]

    case ("ky")
      smallk%ky = [ ((i-1)*cfg%waveVectorMax/(cfg%waveVectorStep-1), i=1,cfg%waveVectorStep) ]

    case ("kz")
      smallk%kz = [ ((i-1)*cfg%waveVectorMax/(cfg%waveVectorStep-1), i=1,cfg%waveVectorStep) ]

    case ("kxky")
      ! [110] direction: kx = ky = k, kz = 0
      do i = 1, cfg%waveVectorStep
        smallk(i)%kx = (i-1) * cfg%waveVectorMax / (cfg%waveVectorStep - 1)
        smallk(i)%ky = smallk(i)%kx
      end do

    case ("kxkz")
      ! [101] direction: kx = kz = k, ky = 0
      do i = 1, cfg%waveVectorStep
        smallk(i)%kx = (i-1) * cfg%waveVectorMax / (cfg%waveVectorStep - 1)
        smallk(i)%kz = smallk(i)%kx
      end do

    case ("kykz")
      ! [011] direction: ky = kz = k, kx = 0
      do i = 1, cfg%waveVectorStep
        smallk(i)%ky = (i-1) * cfg%waveVectorMax / (cfg%waveVectorStep - 1)
        smallk(i)%kz = smallk(i)%ky
      end do

    case ("k0")
      ! Already zeroed above

    case default
      stop "no such direction"

  end select

  ! ====================================================================
  ! Wire mode (confinement=2): separate sparse path using FEAST
  ! ====================================================================
  if (cfg%confinement == 2) then

    ! Initialize 2D confinement operators (sparse CSR kpterms)
    call confinementInitialization_2d(cfg%grid, cfg%params, cfg%regions, &
      & profile_2d, kpterms_2d, cfg%FDorder)

    ! --- Compute strain if enabled (wire mode) ---
    if (cfg%strain%enabled) then
      block
        type(strain_result) :: strain_out
        real(kind=dp) :: a0_ref

        ! Reference lattice constant from substrate (first material)
        a0_ref = cfg%params(1)%a0

        call compute_strain(cfg%grid, cfg%params, cfg%grid%material_id, &
          cfg%strain, a0_ref, strain_out)

        ! Write strain tensor to file (gnuplot splot format)
        ! Note: eps_xy, eps_xz are zero for plane-strain formulation
        call ensure_output_dir()
        call get_unit(iounit)
        open(unit=iounit, file='output/strain.dat', status='replace', action='write')
        write(iounit, '(A)') '# x(A) y(A) eps_xx eps_yy eps_zz eps_xy eps_xz eps_yz'
        do jj = 1, cfg%grid%ny
          do ii = 1, cfg%grid%nx
            k = (jj - 1) * cfg%grid%nx + ii
            write(unit=iounit, fmt='(8(g14.6,1x))') &
              & cfg%grid%x(ii), cfg%grid%z(jj), &
              & strain_out%eps_xx(k), strain_out%eps_yy(k), &
              & strain_out%eps_zz(k), 0.0_dp, 0.0_dp, &
              & strain_out%eps_yz(k)
          end do
          ! Blank line between y-rows for gnuplot splot
          write(iounit, '(A)') ''
        end do
        close(iounit)
        print *, '  Wire strain tensor written to output/strain.dat'

        call apply_pikus_bir(strain_out, cfg%params, cfg%grid%material_id, &
          cfg%grid, profile_2d)
        call strain_result_free(strain_out)

        print *, '  Wire strain calculation complete'
      end block
    end if

    Ngrid = grid_ngrid(cfg%grid)
    Ntot  = 8 * Ngrid

    ! Number of eigenvalues to request from FEAST
    nev_wire = cfg%numcb + cfg%numvb
    if (nev_wire > Ntot) then
      print *, "Warning: requesting more bands than matrix size. Clamping."
      nev_wire = Ntot
    end if

    print *, ''
    print *, '=== Wire mode (2D confinement) ==='
    print *, '  Grid: nx=', cfg%grid%nx, ' ny=', cfg%grid%ny, ' Ngrid=', Ngrid
    print *, '  Matrix size: ', Ntot, 'x', Ntot
    print *, '  Requesting ', nev_wire, ' eigenvalues per k-point'
    print *, '  kz sweep: ', cfg%waveVectorStep, ' points, kz_max=', cfg%waveVectorMax
    print *, ''

    ! Set up eigensolver configuration
    if (cfg%feast_m0 < 0) then
      ! Negative feast_m0 signals: use dense LAPACK instead of FEAST
      eigen_cfg%method = 'DENSE'
    else
      eigen_cfg%method = 'FEAST'
    end if
    eigen_cfg%nev      = nev_wire
    eigen_cfg%emin     = -1.0_dp   ! placeholder, set by auto_compute_energy_window
    eigen_cfg%emax     =  1.0_dp   ! placeholder
    eigen_cfg%max_iter = 100
    eigen_cfg%tol      = 1.0e-10_dp
    eigen_cfg%feast_m0 = cfg%feast_m0  ! 0 means auto: 2*nev

    ! Allocate storage for eigenvalues across k-points
    ! Use Ntot as upper bound: range mode can return many eigenvalues
    allocate(eig_wire(Ntot, cfg%waveVectorStep))
    eig_wire = 0.0_dp
    max_nev_found = 0

    ! ==================================================================
    ! Wire self-consistent Schrodinger-Poisson loop
    ! ==================================================================
    if (cfg%sc%enabled == 1) then
      block
        ! SC kx sweep arrays (separate from production kz sweep)
        integer :: nk_sc, nev_sc
        real(kind=dp), allocatable :: eig_sc(:,:)
        complex(kind=dp), allocatable :: eigv_sc(:,:,:)
        real(kind=dp), allocatable :: kx_grid_sc(:)
        integer :: ik
        ! SC diagnostics output arrays
        real(kind=dp), allocatable :: sc_phi(:,:), sc_ne(:,:), sc_nh(:,:)

        nk_sc = cfg%sc%num_kpar
        if (mod(nk_sc, 2) == 0) nk_sc = nk_sc - 1
        nev_sc = nev_wire

        allocate(eig_sc(nev_sc, nk_sc))
        allocate(eigv_sc(Ntot, nev_sc, nk_sc))
        allocate(kx_grid_sc(nk_sc))
        eig_sc = 0.0_dp
        eigv_sc = cmplx(0.0_dp, 0.0_dp, kind=dp)

        ! Build kx grid for SC charge integration
        do ik = 1, nk_sc
          kx_grid_sc(ik) = real(ik - 1, kind=dp) * cfg%sc%kpar_max &
            & / real(nk_sc - 1, kind=dp)
        end do

        print *, ''
        print *, '=== Running wire self-consistent Schrodinger-Poisson loop ==='

        call self_consistent_loop_wire(profile_2d, cfg, kpterms_2d, cfg%grid, &
          & coo_cache, eigen_cfg, eig_sc, eigv_sc, &
          & phi_out=sc_phi, n_electron_out=sc_ne, n_hole_out=sc_nh)

        print *, '  Wire SC loop complete. profile_2d updated.'
        print *, ''

        ! --- Write SC diagnostics ---
        call ensure_output_dir()

        ! sc_potential_profile.dat: band edge profile after SC
        call get_unit(iounit)
        open(unit=iounit, file='output/sc_potential_profile.dat', status='replace', action='write')
        write(iounit, '(A)') '# x(A) y(A) EV EV_DeltaSO EC'
        do jj = 1, cfg%grid%ny
          do ii = 1, cfg%grid%nx
            k = (jj - 1) * cfg%grid%nx + ii
            write(unit=iounit, fmt='(5(g14.6,1x))') &
              & cfg%grid%x(ii), cfg%grid%z(jj), &
              & profile_2d(k, 1), profile_2d(k, 2), profile_2d(k, 3)
          end do
          write(iounit, '(A)') ''
        end do
        close(iounit)
        print *, '  SC band edge profile written to output/sc_potential_profile.dat'

        ! sc_phi.dat: electrostatic potential phi(x,y)
        call get_unit(iounit)
        open(unit=iounit, file='output/sc_phi.dat', status='replace', action='write')
        write(iounit, '(A)') '# x(A) y(A) phi'
        do jj = 1, cfg%grid%ny
          do ii = 1, cfg%grid%nx
            write(unit=iounit, fmt='(3(g14.6,1x))') &
              & cfg%grid%x(ii), cfg%grid%z(jj), sc_phi(ii, jj)
          end do
          write(iounit, '(A)') ''
        end do
        close(iounit)
        print *, '  SC electrostatic potential written to output/sc_phi.dat'

        ! sc_charge.dat: electron and hole density
        call get_unit(iounit)
        open(unit=iounit, file='output/sc_charge.dat', status='replace', action='write')
        write(iounit, '(A)') '# x(A) y(A) n_e n_h'
        do jj = 1, cfg%grid%ny
          do ii = 1, cfg%grid%nx
            write(unit=iounit, fmt='(4(g14.6,1x))') &
              & cfg%grid%x(ii), cfg%grid%z(jj), sc_ne(ii, jj), sc_nh(ii, jj)
          end do
          write(iounit, '(A)') ''
        end do
        close(iounit)
        print *, '  SC charge density written to output/sc_charge.dat'

        deallocate(eig_sc, eigv_sc, kx_grid_sc)
        deallocate(sc_phi, sc_ne, sc_nh)
      end block

      ! Reset COO cache — profile_2d changed, need fresh CSR structure
      call wire_coo_cache_free(coo_cache)
    end if

    ! --- Write 2D band edge profile (after strain + SC, before k-sweep) ---
    call ensure_output_dir()
    call get_unit(iounit)
    open(unit=iounit, file='output/potential_profile.dat', status='replace', action='write')
    write(iounit, '(A)') '# x(A) y(A) EV EV_DeltaSO EC'
    do jj = 1, cfg%grid%ny
      do ii = 1, cfg%grid%nx
        k = (jj - 1) * cfg%grid%nx + ii
        write(unit=iounit, fmt='(5(g14.6,1x))') &
          & cfg%grid%x(ii), cfg%grid%z(jj), &
          & profile_2d(k, 1), profile_2d(k, 2), profile_2d(k, 3)
      end do
      ! Blank line between y-rows for gnuplot splot
      write(iounit, '(A)') ''
    end do
    close(iounit)
    print *, '  Wire band edge profile written to output/potential_profile.dat'

    ! ==================================================================
    ! Wire kz sweep: serial k=1 (build COO cache) + OpenMP k=2..N
    ! ==================================================================

    ! --- Serial k=1: build COO cache, auto energy window, solve ---
    print *, 'k-point 1/', cfg%waveVectorStep, ' kz=', smallk(1)%kz

    call ZB8bandGeneralized(HT_csr, smallk(1)%kz, profile_2d, &
      & kpterms_2d, cfg, coo_cache)

    call auto_compute_energy_window(HT_csr, eigen_cfg%emin, eigen_cfg%emax)
    ! Override with manual window if specified in config
    if (cfg%feast_emin /= 0.0_dp .or. cfg%feast_emax /= 0.0_dp) then
      eigen_cfg%emin = cfg%feast_emin
      eigen_cfg%emax = cfg%feast_emax
      print *, '  Manual energy window: [', eigen_cfg%emin, ',', eigen_cfg%emax, ']'
    else
      print *, '  Auto energy window: [', eigen_cfg%emin, ',', eigen_cfg%emax, ']'
    end if

    call solve_sparse_evp(HT_csr, eigen_cfg, eigen_res)

    if (.not. eigen_res%converged) then
      if (eigen_res%nev_found < nev_wire) then
        print *, '  ERROR: FEAST did not converge at k-point 1 and found only', &
          eigen_res%nev_found, 'eigenvalues (need', nev_wire, ')'
        stop 1
      end if
      print *, '  WARNING: FEAST subspace issue at k-point 1, but found enough eigenvalues'
    end if
    print *, '  Found ', eigen_res%nev_found, ' eigenvalues (requested ', &
      & nev_wire, ') iterations=', eigen_res%iterations

    if (eigen_res%nev_found > 0) then
      do i = 1, min(eigen_res%nev_found, Ntot)
        eig_wire(i, 1) = eigen_res%eigenvalues(i)
      end do
      max_nev_found = max(max_nev_found, eigen_res%nev_found)
    end if

    ! Write 2D eigenfunctions for k=1 (with band decomposition)
    if (eigen_res%nev_found > 0) then
      call writeEigenfunctions2d(cfg%grid, eigen_res%eigenvalues, &
        & eigen_res%eigenvectors, 1, eigen_res%nev_found, write_parts=.true.)
      print *, '  Wrote ', eigen_res%nev_found, ' 2D wavefunctions to output/eigenfunctions_k_00001_ev_*.dat'
    end if

    call eigensolver_result_free(eigen_res)

    ! --- OpenMP parallel k=2..N ---
    ! HT_csr (from k=1) provides the CSR sparsity template for each thread.
    if (cfg%waveVectorStep > 1) then
      ! Disable MKL internal threading so each OpenMP thread calls
      ! FEAST in serial — avoids oversubscription.
      info = mkl_set_num_threads_local(1)

      print '(A,I0,A)', ' Wire kz-sweep: k=2..', cfg%waveVectorStep, &
        & ' (OpenMP parallel)'

      block
        type(csr_matrix)          :: HT_csr_loc
        type(eigensolver_result)  :: eigen_res_loc

        !$omp parallel private(k, i, info, HT_csr_loc, eigen_res_loc)
        ! Each thread constrains MKL to 1 thread locally
        info = mkl_set_num_threads_local(1)

        !$omp do schedule(static)
        do k = 2, cfg%waveVectorStep
          ! Clone CSR sparsity structure for this k-point
          call csr_clone_structure(HT_csr, HT_csr_loc)
          ! Build sparse Hamiltonian (reuses read-only COO cache from k=1)
          call ZB8bandGeneralized(HT_csr_loc, smallk(k)%kz, profile_2d, &
            & kpterms_2d, cfg, coo_cache)

          ! Solve sparse eigenvalue problem
          call solve_sparse_evp(HT_csr_loc, eigen_cfg, eigen_res_loc)

          if (.not. eigen_res_loc%converged) then
            if (eigen_res_loc%nev_found < nev_wire) then
              !$omp critical
              print *, '  ERROR: FEAST did not converge at k-point', k, &
                'and found only', eigen_res_loc%nev_found, &
                'eigenvalues (need', nev_wire, ')'
              !$omp end critical
              stop 1
            end if
          end if

          ! Store eigenvalues (each thread writes to its own k column)
          if (eigen_res_loc%nev_found > 0) then
            do i = 1, min(eigen_res_loc%nev_found, Ntot)
              eig_wire(i, k) = eigen_res_loc%eigenvalues(i)
            end do
            !$omp critical
            max_nev_found = max(max_nev_found, eigen_res_loc%nev_found)
            !$omp end critical
          end if

          ! Clean up per-thread eigensolver result and CSR
          call eigensolver_result_free(eigen_res_loc)
          call csr_free(HT_csr_loc)
        end do
        !$omp end do

        !$omp end parallel
      end block
    end if

    ! --- Serial: write 2D eigenfunctions for middle and last k-points ---
    ! (HT_csr still alive as CSR sparsity template)
    block
      type(csr_matrix)         :: HT_csr_wf
      type(eigensolver_result) :: eigen_res_wf
      integer :: k_wf

      do k_wf = 1, cfg%waveVectorStep
        if (k_wf /= 1 .and. &
          & k_wf /= cfg%waveVectorStep/2 .and. &
          & k_wf /= cfg%waveVectorStep) cycle
        ! k=1 already written above
        if (k_wf == 1) cycle

        ! Clone CSR structure from k=1 template
        call csr_clone_structure(HT_csr, HT_csr_wf)
        call ZB8bandGeneralized(HT_csr_wf, smallk(k_wf)%kz, profile_2d, &
          & kpterms_2d, cfg, coo_cache)
        call solve_sparse_evp(HT_csr_wf, eigen_cfg, eigen_res_wf)

        if (eigen_res_wf%nev_found > 0) then
          call writeEigenfunctions2d(cfg%grid, eigen_res_wf%eigenvalues, &
            & eigen_res_wf%eigenvectors, k_wf, eigen_res_wf%nev_found, write_parts=.true.)
          print *, '  Wrote ', eigen_res_wf%nev_found, &
            & ' 2D wavefunctions to output/eigenfunctions_k_', k_wf, '_ev_*.dat'
        end if

        call eigensolver_result_free(eigen_res_wf)
        call csr_free(HT_csr_wf)
      end do
    end block

    ! Free Hamiltonian CSR and COO cache
    call csr_free(HT_csr)
    call wire_coo_cache_free(coo_cache)

    ! Write eigenvalues to file (only the rows actually populated)
    if (max_nev_found > 0) then
      call writeEigenvalues(smallk, eig_wire(1:max_nev_found, :), cfg%waveVectorStep, cfg)
    else
      call writeEigenvalues(smallk, eig_wire(1:nev_wire, :), cfg%waveVectorStep, cfg)
    end if
    print *, ''
    print *, 'Wire band structure written to output/eigenvalues.dat'

    ! Clean up wire mode allocations
    if (allocated(eig_wire))    deallocate(eig_wire)
    if (allocated(profile_2d))  deallocate(profile_2d)
    if (allocated(kpterms_2d)) then
      do i = 1, size(kpterms_2d)
        call csr_free(kpterms_2d(i))
      end do
      deallocate(kpterms_2d)
    end if

    ! Clean up common allocations and exit
    if (allocated(smallk))  deallocate(smallk)
    if (allocated(profile)) deallocate(profile)
    if (allocated(kpterms)) deallocate(kpterms)

    stop  ! wire mode complete

  end if

  ! ====================================================================
  ! Bulk / QW mode (confinement=0 or 1): dense LAPACK path
  ! ====================================================================

  ! Set matrix dimensions and eigenvalue range
  N = cfg%fdStep * 8
  if (cfg%confDir == 'n') then
    N = 8  ! For bulk, we only need 8x8
    if (cfg%evnum > 8) then
      print *, "Warning: requesting more bands than available in bulk (8). Using all 8 bands."
      cfg%evnum = 8
      cfg%numcb = 2
      cfg%numvb = 6
    end if
  end if

  ! For quantum well, check if requested bands are within limits
  if (cfg%confDir == 'z') then
    if (cfg%numcb > NUM_CB_STATES*cfg%fdStep) then
      print *, "Warning: requesting more conduction bands than available. Limiting to ", NUM_CB_STATES*cfg%fdStep
      cfg%numcb = NUM_CB_STATES*cfg%fdStep
    end if
    if (cfg%numvb > NUM_VB_STATES*cfg%fdStep) then
      print *, "Warning: requesting more valence bands than available. Limiting to ", NUM_VB_STATES*cfg%fdStep
      cfg%numvb = NUM_VB_STATES*cfg%fdStep
    end if
    cfg%evnum = cfg%numcb + cfg%numvb
  end if

  vl = 0.0_dp
  vu = 0.0_dp
  if (cfg%confinement == 0) then
    il = 1
    iuu = cfg%evnum
  else
    ! For quantum well, select the right range of states
    ! We want the highest numvb valence states and lowest numcb conduction states
    il = NUM_VB_STATES*cfg%fdStep - cfg%numvb + 1  ! Start from highest valence band
    iuu = NUM_VB_STATES*cfg%fdStep + cfg%numcb     ! Up to highest conduction band
    print *, "Computing states from index", il, "to", iuu
  end if

  NB = ILAENV(1, 'ZHETRD', 'UPLO', N, N, -1, -1)
  NB = MAX(NB,N)
  ABSTOL = DLAMCH('P')

  if (allocated(rwork)) deallocate(rwork)
  allocate(rwork(7*N))  ! For ZHEEVX
  if (allocated(iwork)) deallocate(iwork)
  allocate(iwork(5*N))
  if (allocated(ifail)) deallocate(ifail)
  allocate(ifail(N))
  if (allocated(eig)) deallocate(eig)
  allocate(eig(iuu-il+1,cfg%waveVectorStep))  ! Only store the states we want
  if (allocated(eigv)) deallocate(eigv)
  if (cfg%confDir == 'n') then
    allocate(eigv(8,cfg%evnum,cfg%waveVectorStep))  ! 8x8 for bulk
  else
    allocate(eigv(N,iuu-il+1,cfg%waveVectorStep))  ! Only store the states we want
  end if

  eig(:,:) = 0_dp
  eigv(:,:,:) = 0_dp

  allocate(HT(N,N))
  allocate(HTmp(8,8))  ! Temporary array for bulk diagonalization
  HT = 0.0_dp
  HTmp = 0.0_dp

  ! Initial workspace allocation
  if (allocated(work)) deallocate(work)
  allocate(work(N))  ! Will be resized after workspace query

  ! Print profile for QW mode
  if (cfg%confDir == 'z') then
    call get_unit(iounit)
    open(unit=iounit, file='output/potential_profile.dat', status="replace", action="write")
    do i = 1, cfg%fdStep, 1
      write(iounit,*) cfg%z(i), profile(i,1), profile(i,2), profile(i,3)
    end do
    close(iounit)
  end if

  ! --- Compute strain if enabled (QW mode) ---
  if (cfg%strain%enabled .and. cfg%confinement == 1) then
    block
      type(strain_result) :: strain_out_qw
      real(kind=dp) :: a0_ref

      ! Reference lattice constant from substrate (first material)
      a0_ref = cfg%params(1)%a0

      call compute_strain(cfg%grid, cfg%params, cfg%grid%material_id, &
        cfg%strain, a0_ref, strain_out_qw)
      call apply_pikus_bir(strain_out_qw, cfg%params, cfg%grid%material_id, &
        cfg%grid, profile)
      call strain_result_free(strain_out_qw)

      print *, '  QW strain calculation complete'
    end block
  end if

  ! --- Self-consistent loop (QW only) ---
  if (cfg%confDir == 'z' .and. cfg%sc%enabled == 1) then
    print *, ''
    print *, '=== Running self-consistent Schrodinger-Poisson loop ==='
    call self_consistent_loop(profile, cfg, kpterms, HT, eig, eigv, &
      & smallk, N, il, iuu, n_electron_out=sc_ne_qw, n_hole_out=sc_nh_qw)

    ! Write updated profile after SC convergence
    call get_unit(iounit)
    open(unit=iounit, file='output/sc_potential_profile.dat', status="replace", action="write")
    do i = 1, cfg%fdStep, 1
      write(iounit,*) cfg%z(i), profile(i,1), profile(i,2), profile(i,3)
    end do
    close(iounit)
    print *, 'SC potential profile written to output/sc_potential_profile.dat'

    ! Write charge density
    if (allocated(sc_ne_qw) .and. allocated(sc_nh_qw)) then
      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/sc_charge.dat', status="replace", action="write")
      write(iounit, '(A)') '# z(A) n_e(cm^-3) n_h(cm^-3)'
      do i = 1, cfg%fdStep, 1
        write(iounit, '(3(g14.6,1x))') cfg%z(i), sc_ne_qw(i), sc_nh_qw(i)
      end do
      close(iounit)
      print *, 'SC charge density written to output/sc_charge.dat'
      deallocate(sc_ne_qw, sc_nh_qw)
    end if
  end if

  ! ====================================================================
  ! Workspace query (serial, k=1) — determines lwork for zheevx
  ! ====================================================================
  if (cfg%confDir == 'n') then
    ! BULK: 8x8 workspace query
    call ZB8bandBulk(HT, smallk(1), cfg%params(1))
    HTmp = HT(1:8,1:8)
    if (allocated(work)) deallocate(work)
    allocate(work(1))
    lwork = -1
    call zheevx('V', 'I', 'U', 8, HTmp, 8, vl, vu, il, iuu, abstol, M, eig(1:cfg%evnum,1), &
               HTmp, 8, work, lwork, rwork, iwork, ifail, info)
    if (info /= 0) then
      print *, 'Error: zheevx bulk workspace query failed, info =', info
      stop 1
    end if
    lwork = int(real(work(1)))
    if (allocated(work)) deallocate(work)
    allocate(work(lwork))
  else if (cfg%confDir == 'z') then
    ! QW: NxN workspace query
    call ZB8bandQW(HT, smallk(1), profile, kpterms)
    if (allocated(work)) deallocate(work)
    allocate(work(1))
    lwork = -1
    call zheevx('V', 'I', 'U', N, HT, N, vl, vu, il, iuu, abstol, M, eig(:,1), &
               HT, N, work, lwork, rwork, iwork, ifail, info)
    if (info /= 0) then
      print *, 'Error: zheevx QW workspace query failed, info =', info
      stop 1
    end if
    lwork = int(real(work(1)))
    if (allocated(work)) deallocate(work)
    allocate(work(lwork))
  end if

  ! ====================================================================
  ! k-vector sweep: sequential for bulk, OpenMP parallel for QW
  ! ====================================================================
  if (cfg%confDir == 'n') then
    ! --- BULK (8x8, trivially fast — no parallelization needed) ---
    do k = 1, cfg%waveVectorStep
      call ZB8bandBulk(HT, smallk(k), cfg%params(1))
      HTmp = HT(1:8,1:8)

      call zheevx('V', 'I', 'U', 8, HTmp, 8, vl, vu, il, iuu, abstol, M, eig(1:cfg%evnum,k), &
                 HTmp, 8, work, lwork, rwork, iwork, ifail, info)
      if (info /= 0) then
        print *, "Diagonalization error in bulk calculation, info = ", info
        if (info < 0) print *, "Parameter ", -info, " had illegal value"
        stop "error diag"
      end if

      eigv(:,:,k) = HTmp(:,1:cfg%evnum)
    end do

    ! Write bulk eigenfunctions at start, middle, end k-points
    do k = 1, cfg%waveVectorStep
      if (k == 1 .or. k == int(cfg%waveVectorStep/2) .or. k == cfg%waveVectorStep) then
        call writeEigenfunctions(8, min(cfg%evnum,8), eigv(:,1:min(cfg%evnum,8),k), &
          & k, cfg%fdstep, cfg%z, .true., &
          & k_magnitude=sqrt(smallk(k)%kx**2 + smallk(k)%ky**2 + smallk(k)%kz**2))
      end if
    end do

  else if (cfg%confDir == 'z') then
    ! --- QUANTUM WELL (NxN, OpenMP parallel) ---
    ! Disable MKL internal threading so each OpenMP thread calls
    ! zheevx in serial — avoids oversubscription with intel_thread MKL.
    ! mkl_set_num_threads_local returns previous setting (discarded).
    info = mkl_set_num_threads_local(1)

    print '(A,I0,A)', ' QW k-sweep: ', cfg%waveVectorStep, ' k-points (OpenMP parallel)'

    block
      ! Thread-private temporaries for the parallel region
      complex(kind=dp), allocatable :: HT_loc(:,:), work_loc(:)
      real(kind=dp), allocatable    :: rwork_loc(:)
      integer, allocatable          :: iwork_loc(:), ifail_loc(:)
      integer :: info_loc, M_loc

      !$omp parallel private(k, HT_loc, work_loc, rwork_loc, iwork_loc, ifail_loc, info_loc, M_loc)
      ! Each thread allocates its own workspace
      allocate(HT_loc(N, N))
      allocate(work_loc(lwork))
      allocate(rwork_loc(7*N))
      allocate(iwork_loc(5*N))
      allocate(ifail_loc(N))
      HT_loc = (0.0_dp, 0.0_dp)

      !$omp do schedule(static)
      do k = 1, cfg%waveVectorStep
        ! Build Hamiltonian for this k-point
        call ZB8bandQW(HT_loc, smallk(k), profile, kpterms)

        ! Diagonalize
        call zheevx('V', 'I', 'U', N, HT_loc, N, vl, vu, il, iuu, abstol, M_loc, &
                   eig(:,k), HT_loc, N, work_loc, lwork, rwork_loc, iwork_loc, &
                   ifail_loc, info_loc)
        if (info_loc /= 0) then
          !$omp critical
          print *, "ERROR: diagonalization failed at k=", k, "info=", info_loc
          if (info_loc < 0) print *, "  Parameter ", -info_loc, " had illegal value"
          !$omp end critical
          stop 1
        end if

        ! Store eigenvectors (HT_loc now holds them, zheevx overwrites input)
        eigv(:,:,k) = HT_loc(:, 1:iuu-il+1)
      end do
      !$omp end do

      deallocate(HT_loc, work_loc, rwork_loc, iwork_loc, ifail_loc)
      !$omp end parallel
    end block

    ! Write QW eigenfunctions at start, middle, end k-points (serial)
    do k = 1, cfg%waveVectorStep
      if (k == 1 .or. k == int(cfg%waveVectorStep/2) .or. k == cfg%waveVectorStep) then
        call writeEigenfunctions(N, iuu-il+1, eigv(:,1:iuu-il+1,k), &
          & k, cfg%fdstep, cfg%z, .false.)
      end if
    end do

  end if
  call writeEigenvalues(smallk, eig(:,:), cfg%waveVectorStep)


  !----------------------------------------------------------------------------


  if (allocated(smallk)) deallocate(smallk)
  if (allocated(HT)) deallocate(HT)
  if (allocated(HTmp)) deallocate(HTmp)
  if (allocated(eig)) deallocate(eig)
  if (allocated(eigv)) deallocate(eigv)
  if (allocated(work)) deallocate(work)
  if (allocated(iwork)) deallocate(iwork)
  if (allocated(rwork)) deallocate(rwork)
  if (allocated(ifail)) deallocate(ifail)


end program kpfdm
