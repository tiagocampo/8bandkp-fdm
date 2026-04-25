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
  use optical_spectra
  use exciton_solver
  use scattering_solver
  use utils, only: dnscsr_z_mkl
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
  real(kind=dp), allocatable       :: prev_wire_eval(:)
  complex(kind=dp), allocatable    :: prev_wire_evec(:,:)
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

        call compute_bir_pikus_blocks(strain_out, cfg%params, cfg%grid%material_id, &
          cfg%grid, cfg%strain_blocks)
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
    if (allocated(cfg%strain_blocks%delta_Ec)) then
      write(iounit, '(A)') '# x(A) y(A) EV_top_strained ESO_strained EC_strained'
    else
      write(iounit, '(A)') '# x(A) y(A) EV EV_DeltaSO EC'
    end if
    do jj = 1, cfg%grid%ny
      do ii = 1, cfg%grid%nx
        k = (jj - 1) * cfg%grid%nx + ii
        if (allocated(cfg%strain_blocks%delta_Ec)) then
          write(unit=iounit, fmt='(5(g14.6,1x))') &
            & cfg%grid%x(ii), cfg%grid%z(jj), &
            & max(profile_2d(k, 1) + cfg%strain_blocks%delta_EHH(k), &
                  profile_2d(k, 1) + cfg%strain_blocks%delta_ELH(k)), &
            & profile_2d(k, 2) + cfg%strain_blocks%delta_ESO(k), &
            & profile_2d(k, 3) + cfg%strain_blocks%delta_Ec(k)
        else
          write(unit=iounit, fmt='(5(g14.6,1x))') &
            & cfg%grid%x(ii), cfg%grid%z(jj), &
            & profile_2d(k, 1), profile_2d(k, 2), profile_2d(k, 3)
        end if
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

    if (allocated(prev_wire_eval)) deallocate(prev_wire_eval)
    if (allocated(prev_wire_evec)) deallocate(prev_wire_evec)
    if (eigen_res%nev_found > 0) then
      allocate(prev_wire_eval(eigen_res%nev_found))
      allocate(prev_wire_evec(Ntot, eigen_res%nev_found))
      prev_wire_eval = eigen_res%eigenvalues
      prev_wire_evec = eigen_res%eigenvectors
    end if

    call eigensolver_result_free(eigen_res)

    ! --- Serial k=2..N with branch tracking against the previous k-point ---
    if (cfg%waveVectorStep > 1) then
      info = mkl_set_num_threads_local(1)

      print '(A,I0,A)', ' Wire kz-sweep: k=2..', cfg%waveVectorStep, &
        & ' (serial branch tracking)'

      do k = 2, cfg%waveVectorStep
        print *, 'k-point ', k, '/', cfg%waveVectorStep, ' kz=', smallk(k)%kz

        block
          type(csr_matrix) :: HT_csr_step

          call csr_clone_structure(HT_csr, HT_csr_step)
          call ZB8bandGeneralized(HT_csr_step, smallk(k)%kz, profile_2d, &
            & kpterms_2d, cfg, coo_cache)
          call solve_sparse_evp(HT_csr_step, eigen_cfg, eigen_res)
          call csr_free(HT_csr_step)
        end block

        if (.not. eigen_res%converged) then
          if (eigen_res%nev_found < nev_wire) then
            print *, '  ERROR: FEAST did not converge at k-point', k, &
              'and found only', eigen_res%nev_found, 'eigenvalues (need', nev_wire, ')'
            stop 1
          end if
          print *, '  WARNING: eigensolver subspace issue at k-point', k, &
            ' but found enough eigenvalues'
        end if

        if (eigen_res%nev_found > 0) then
          call reorder_wire_branches(prev_wire_eval, prev_wire_evec, &
            eigen_res%eigenvalues, eigen_res%eigenvectors)
          do i = 1, min(eigen_res%nev_found, Ntot)
            eig_wire(i, k) = eigen_res%eigenvalues(i)
          end do
          max_nev_found = max(max_nev_found, eigen_res%nev_found)

          if (allocated(prev_wire_eval)) deallocate(prev_wire_eval)
          if (allocated(prev_wire_evec)) deallocate(prev_wire_evec)
          allocate(prev_wire_eval(eigen_res%nev_found))
          allocate(prev_wire_evec(Ntot, eigen_res%nev_found))
          prev_wire_eval = eigen_res%eigenvalues
          prev_wire_evec = eigen_res%eigenvectors

          if (k == cfg%waveVectorStep/2 .or. k == cfg%waveVectorStep) then
            call writeEigenfunctions2d(cfg%grid, eigen_res%eigenvalues, &
              & eigen_res%eigenvectors, k, eigen_res%nev_found, write_parts=.false.)
            print *, '  Wrote ', eigen_res%nev_found, &
              & ' 2D wavefunctions to output/eigenfunctions_k_', k, '_ev_*.dat'
          end if
        end if

        call eigensolver_result_free(eigen_res)
      end do
    end if

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
    if (allocated(prev_wire_eval)) deallocate(prev_wire_eval)
    if (allocated(prev_wire_evec)) deallocate(prev_wire_evec)
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
    call ensure_output_dir()
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
      call compute_bir_pikus_blocks(strain_out_qw, cfg%params, cfg%grid%material_id, &
        cfg%grid, cfg%strain_blocks)
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
    call ZB8bandQW(HT, smallk(1), profile, kpterms, cfg=cfg)
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
        call ZB8bandQW(HT_loc, smallk(k), profile, kpterms, cfg=cfg)

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

  ! ====================================================================
  ! Post-sweep optical absorption (QW only, serial)
  ! ====================================================================
  if (cfg%optics%enabled .and. cfg%confDir == 'z') then
    block
      real(kind=dp), allocatable :: simpson_w(:)
      real(kind=dp) :: dk
      integer :: npts

      ! Simpson 1/3 rule weights for k_sweep integration.
      ! For a QW, the k_parallel integral over 2D BZ with cylindrical
      ! symmetry is: int d^2k = 2*pi * int k dk.  The 1D kx sweep
      ! needs an additional 2*pi*k factor per point.
      npts = cfg%waveVectorStep
      if (mod(npts, 2) == 0) then
        npts = npts - 1  ! Simpson requires odd number of points
        print '(A,I0,A)', ' Warning: optics integration uses ', npts, &
          & ' k-points (Simpson requires odd count)'
      end if
      allocate(simpson_w(npts))
      dk = cfg%waveVectorMax / real(cfg%waveVectorStep - 1, kind=dp)  ! true grid spacing
      do i = 1, npts
        ! Base Simpson weight
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

      ! Build velocity matrices for optical transitions.
      ! For a QW (1D confinement along z):
      !   vel(1) = dH/dk_x at k=0  (in-plane, Kane P matrix)
      !   vel(2) = dH/dk_y at k=0  (in-plane, Kane P matrix)
      !   vel(3) = -i [r_z, H]     (confinement, commutator on CSR)
      ! All three are k-independent, built once.
      block
        type(csr_matrix) :: vel_opt(3)
        type(csr_matrix) :: H_csr_tmp
        complex(kind=dp), allocatable :: H_k0(:,:)
        type(wavevector) :: pert_kx, pert_ky
        integer :: nzmax_tmp

        ! --- vel(3): z-direction via commutator ---
        allocate(H_k0(N, N))
        H_k0 = (0.0_dp, 0.0_dp)
        call ZB8bandQW(H_k0, smallk(1), profile, kpterms, cfg=cfg)
        nzmax_tmp = N * N
        call dnscsr_z_mkl(nzmax_tmp, N, H_k0, H_csr_tmp)
        deallocate(H_k0)

        ! Build commutator velocity (only vel(3) is non-zero for 1D QW)
        call build_velocity_matrices(H_csr_tmp, cfg%grid, vel_opt)
        call csr_free(H_csr_tmp)

        ! --- vel(1): x-direction via dH/dk_x ---
        pert_kx%kx = 1.0_dp
        pert_kx%ky = 0.0_dp
        pert_kx%kz = 0.0_dp
        allocate(H_k0(N, N))
        H_k0 = (0.0_dp, 0.0_dp)
        call ZB8bandQW(H_k0, pert_kx, profile, kpterms, cfg=cfg, g='g')
        nzmax_tmp = N * N
        call dnscsr_z_mkl(nzmax_tmp, N, H_k0, H_csr_tmp)
        deallocate(H_k0)
        call csr_clone_structure(H_csr_tmp, vel_opt(1))
        vel_opt(1)%values = H_csr_tmp%values
        call csr_free(H_csr_tmp)

        ! --- vel(2): y-direction via dH/dk_y ---
        pert_ky%kx = 0.0_dp
        pert_ky%ky = 1.0_dp
        pert_ky%kz = 0.0_dp
        allocate(H_k0(N, N))
        H_k0 = (0.0_dp, 0.0_dp)
        call ZB8bandQW(H_k0, pert_ky, profile, kpterms, cfg=cfg, g='g')
        nzmax_tmp = N * N
        call dnscsr_z_mkl(nzmax_tmp, N, H_k0, H_csr_tmp)
        deallocate(H_k0)
        call csr_clone_structure(H_csr_tmp, vel_opt(2))
        vel_opt(2)%values = H_csr_tmp%values
        call csr_free(H_csr_tmp)

        ! Initialize optics accumulation arrays
        call optics_init(cfg%optics)

        ! Reset gain quasi-Fermi state if gain is enabled
        if (cfg%optics%gain_enabled) then
          call gain_reset()
        end if

        print '(A)', ' Computing optical absorption...'
        if (cfg%optics%gain_enabled) then
          print '(A)', ' Computing gain spectrum alongside absorption...'
        end if
        do k = 1, npts
          call optics_accumulate(cfg%optics, &
            & eig(:, k), eigv(:, :, k), simpson_w(k), &
            & vel_opt, cfg%numcb, cfg%numvb, cfg%sc%fermi_level)

          ! Gain spectrum accumulation (uses separate quasi-Fermi levels)
          if (cfg%optics%gain_enabled) then
            call compute_gain_qw(cfg%optics, &
              & eigvals=eig(:, k), &
              & eigvecs=eigv(:, :, k), &
              & k_weight=simpson_w(k), &
              & nlayers=cfg%numLayers, params=cfg%params, &
              & profile=profile, kpterms=kpterms, &
              & startz=cfg%startPos(1), endz=cfg%endPos(1), dz=cfg%dz, &
              & numcb=cfg%numcb, numvb=cfg%numvb, &
              & carrier_density=cfg%optics%gain_carrier_density)
          end if

          ! ISBT accumulation
          if (cfg%optics%isbt_enabled) then
            call compute_isbt_absorption(cfg%optics, &
              & eigvals=eig(:, k), &
              & eigvecs=eigv(:, :, k), &
              & z_grid=cfg%z, dz=cfg%dz, &
              & numcb=cfg%numcb, numvb=cfg%numvb, fdstep=cfg%fdstep, &
              & k_weight=simpson_w(k), fermi_level=cfg%sc%fermi_level)
          end if
        end do

        ! Free velocity matrices
        do i = 1, 3
          call csr_free(vel_opt(i))
        end do
      end block

      ! ISBT transition table at k=0
      if (cfg%optics%isbt_enabled) then
        call compute_intersubband_transitions(eig(:, 1), eigv(:, :, 1), &
          & cfg%z, cfg%dz, cfg%numcb, cfg%numvb, cfg%fdstep, &
          & 'output/isbt_transitions.dat')
        print '(A)', ' ISBT transitions written to output/isbt_transitions.dat'
      end if

      call optics_finalize(cfg%optics)

      ! Exciton corrections (applied before cleanup to access alpha arrays)
      if (cfg%exciton%enabled) then
        block
          real(kind=dp) :: E_binding_ex, lambda_opt_ex
          real(kind=dp) :: E_gap_ex
          integer :: ie, iounit_ex, cb_st, vb_st
          cb_st = cfg%numvb + 1
          vb_st = cfg%numvb
          E_gap_ex = eig(cb_st, 1) - eig(vb_st, 1)
          call compute_exciton_binding(eig(:, 1), eigv(:, :, 1), &
            & cfg%z, cfg%dz, cfg%numLayers, cfg%params, &
            & cfg%numcb, cfg%numvb, cfg%fdstep, E_binding_ex, lambda_opt_ex, &
            & cfg%grid%material_id)
          print '(A,F8.3,A)', ' Exciton binding energy: ', E_binding_ex, ' meV'
          call apply_excitonic_corrections(E_grid, alpha_te, alpha_tm, &
            & E_gap_ex, E_binding_ex, cfg%optics)
          ! Rewrite absorption files with excitonic corrections
          call ensure_output_dir()
          call get_unit(iounit_ex)
          open(unit=iounit_ex, file='output/absorption_TE.dat', &
            & status='replace', action='write')
          write(iounit_ex, '(a)') '# TE absorption with excitonic corrections'
          write(iounit_ex, '(a)') '# E(eV)  alpha(cm^-1)'
          do ie = 1, nE
            write(iounit_ex, '(es16.8, 2x, es16.8)') E_grid(ie), alpha_te(ie)
          end do
          close(iounit_ex)
          open(unit=iounit_ex, file='output/absorption_TM.dat', &
            & status='replace', action='write')
          write(iounit_ex, '(a)') '# TM absorption with excitonic corrections'
          write(iounit_ex, '(a)') '# E(eV)  alpha(cm^-1)'
          do ie = 1, nE
            write(iounit_ex, '(es16.8, 2x, es16.8)') E_grid(ie), alpha_tm(ie)
          end do
          close(iounit_ex)
          print '(A)', ' Excitonic corrections applied to absorption spectra'
        end block
      end if

      call optics_cleanup()
      deallocate(simpson_w)
      print '(A)', ' Absorption spectra written to output/'
    end block
  end if

  ! ====================================================================
  ! Exciton binding energy (standalone, when optics is disabled)
  ! ====================================================================
  if (cfg%exciton%enabled .and. cfg%confDir == 'z' .and. .not. cfg%optics%enabled) then
    block
      real(kind=dp) :: E_binding_sa, lambda_opt_sa
      call compute_exciton_binding(eig(:, 1), eigv(:, :, 1), &
        & cfg%z, cfg%dz, cfg%numLayers, cfg%params, &
        & cfg%numcb, cfg%numvb, cfg%fdstep, E_binding_sa, lambda_opt_sa, &
        & cfg%grid%material_id)
      print '(A,F8.3,A)', ' Exciton binding energy: ', E_binding_sa, ' meV'
      print '(A,F8.2,A)', ' Variational parameter: ', lambda_opt_sa, ' AA'
    end block
  end if

  ! ====================================================================
  ! LO-phonon scattering rates (QW only, k=0)
  ! ====================================================================
  if (cfg%scattering%enabled .and. cfg%confDir == 'z') then
    call compute_phonon_scattering(cfg, eig(:, 1), eigv(:, :, 1), &
      & cfg%z, cfg%params, cfg%dz, cfg%numcb, cfg%numvb, cfg%fdstep)
    print '(A)', ' Scattering rates written to output/scattering_rates.dat'
  end if


  if (allocated(smallk)) deallocate(smallk)
  if (allocated(HT)) deallocate(HT)
  if (allocated(HTmp)) deallocate(HTmp)
  if (allocated(eig)) deallocate(eig)
  if (allocated(eigv)) deallocate(eigv)
  if (allocated(work)) deallocate(work)
  if (allocated(iwork)) deallocate(iwork)
  if (allocated(rwork)) deallocate(rwork)
  if (allocated(ifail)) deallocate(ifail)


contains

  subroutine reorder_wire_branches(prev_eval, prev_evec, curr_eval, curr_evec)
    real(kind=dp), intent(in) :: prev_eval(:)
    complex(kind=dp), intent(in) :: prev_evec(:,:)
    real(kind=dp), intent(inout) :: curr_eval(:)
    complex(kind=dp), intent(inout) :: curr_evec(:,:)

    integer :: n_prev, n_curr, n_match, i, j, best_j, next_slot
    integer :: block_start, block_end
    real(kind=dp) :: best_score, score, overlap_mag, energy_delta, parts_similarity
    logical, allocatable :: used(:)
    integer, allocatable :: perm(:)
    real(kind=dp), allocatable :: eval_tmp(:)
    complex(kind=dp), allocatable :: evec_tmp(:,:)
    real(kind=dp), allocatable :: prev_parts(:,:), curr_parts(:,:)
    real(kind=dp), parameter :: degeneracy_tol = 1.0e-8_dp

    n_prev = size(prev_eval)
    n_curr = size(curr_eval)
    if (n_prev <= 0 .or. n_curr <= 0) return

    n_match = min(n_prev, n_curr)
    allocate(used(n_curr), perm(n_curr))
    allocate(prev_parts(8, n_prev), curr_parts(8, n_curr))
    used = .false.
    perm = 0

    do i = 1, n_prev
      call compute_wire_band_parts(prev_evec(:, i), prev_parts(:, i))
    end do
    do j = 1, n_curr
      call compute_wire_band_parts(curr_evec(:, j), curr_parts(:, j))
    end do

    do i = 1, n_match
      best_j = 0
      best_score = -huge(1.0_dp)
      do j = 1, n_curr
        if (used(j)) cycle
        overlap_mag = abs(sum(conjg(prev_evec(:, i)) * curr_evec(:, j)))
        energy_delta = abs(curr_eval(j) - prev_eval(i))
        parts_similarity = sum(prev_parts(:, i) * curr_parts(:, j))
        score = overlap_mag + 0.25_dp * parts_similarity - 1.0e-3_dp * energy_delta
        if (score > best_score) then
          best_score = score
          best_j = j
        end if
      end do
      if (best_j > 0) then
        perm(i) = best_j
        used(best_j) = .true.
      end if
    end do

    next_slot = n_match + 1
    do j = 1, n_curr
      if (.not. used(j)) then
        perm(next_slot) = j
        next_slot = next_slot + 1
      end if
    end do

    ! FEAST and dense LAPACK can choose different bases inside an exactly
    ! degenerate manifold at the previous k-point. Keep the matched current
    ! states in ascending-energy order inside those blocks so the first step
    ! away from k=0 remains deterministic across solvers.
    block_start = 1
    do while (block_start <= n_match)
      block_end = block_start
      do while (block_end < n_match)
        if (abs(prev_eval(block_end + 1) - prev_eval(block_start)) > degeneracy_tol) exit
        block_end = block_end + 1
      end do
      if (block_end > block_start) then
        call sort_perm_block_by_energy(perm(block_start:block_end), curr_eval)
      end if
      block_start = block_end + 1
    end do

    allocate(eval_tmp(n_curr))
    allocate(evec_tmp(size(curr_evec, 1), n_curr))
    do i = 1, n_curr
      eval_tmp(i) = curr_eval(perm(i))
      evec_tmp(:, i) = curr_evec(:, perm(i))
    end do

    curr_eval = eval_tmp
    curr_evec = evec_tmp

    deallocate(eval_tmp, evec_tmp, prev_parts, curr_parts, perm, used)
  end subroutine reorder_wire_branches

  subroutine sort_perm_block_by_energy(block_perm, curr_eval)
    integer, intent(inout) :: block_perm(:)
    real(kind=dp), intent(in) :: curr_eval(:)

    integer :: i, j, tmp_idx

    do i = 1, size(block_perm) - 1
      do j = i + 1, size(block_perm)
        if (curr_eval(block_perm(j)) < curr_eval(block_perm(i))) then
          tmp_idx = block_perm(i)
          block_perm(i) = block_perm(j)
          block_perm(j) = tmp_idx
        end if
      end do
    end do
  end subroutine sort_perm_block_by_energy

  subroutine compute_wire_band_parts(state_vec, parts)
    complex(kind=dp), intent(in) :: state_vec(:)
    real(kind=dp), intent(out) :: parts(8)

    integer :: band, ngrid_local, start_idx, end_idx
    real(kind=dp) :: total_weight

    ngrid_local = size(state_vec) / 8
    parts = 0.0_dp
    if (ngrid_local <= 0) return

    do band = 1, 8
      start_idx = (band - 1) * ngrid_local + 1
      end_idx = band * ngrid_local
      parts(band) = real(sum(conjg(state_vec(start_idx:end_idx)) * &
        & state_vec(start_idx:end_idx)), kind=dp)
    end do

    total_weight = sum(parts)
    if (total_weight > 0.0_dp) parts = parts / total_weight
  end subroutine compute_wire_band_parts

end program kpfdm
