program opticalProperties

  use definitions
  use parameters
  use hamiltonianConstructor
  use input_parser, only: read_config
  use simulation_setup_mod, only: simulation_setup, simulation_setup_init, &
    & simulation_setup_free, setup_build_velocity_matrices, setup_build_H
  use optical_spectra
  use exciton_solver, only: compute_exciton_binding, apply_excitonic_corrections
  use sparse_matrices
  use eigensolver, only: eigensolver_result, eigensolver_result_free
  use linalg, only: zheevx, ilaenv, dlamch, mkl_set_num_threads_local
  use outputFunctions, only: ensure_output_dir, get_unit

  implicit none

  ! Shared configuration from input_parser
  type(simulation_config) :: cfg
  type(simulation_setup)  :: setup
  type(optics_engine)     :: oe

  ! Iteration
  integer :: i, k

  ! Hamiltonian and LAPACK/BLAS
  integer :: info, NB, lwork, N, M, il, iuu
  real(kind=dp) :: abstol, vl, vu
  real(kind=dp), allocatable :: eig(:,:), rwork(:)
  complex(kind=dp), allocatable :: work(:)
  complex(kind=dp), allocatable, dimension(:,:,:) :: eigv
  complex(kind=dp), allocatable, dimension(:,:) :: HT
  integer, allocatable :: iwork(:), ifail(:)

  ! k_par sweep
  type(wavevector), allocatable, dimension(:) :: smallk
  integer :: npts
  real(kind=dp) :: dk

  ! Simpson integration weights
  real(kind=dp), allocatable :: simpson_w(:)

  ! Wire-specific variables (confinement=2)
  type(eigensolver_result) :: eigen_res
  integer                  :: Ngrid, Ntot, nev_wire

  ! File handling
  integer(kind=4) :: iounit

  ! ================================================================
  ! Parse input, initialize materials, external fields
  ! ================================================================
  call read_config(cfg)

  ! Check that optics is enabled
  if (.not. cfg%optics%enabled) then
    print '(a)', 'ERROR: optics block not enabled in input.toml'
    print '(a)', '  Add an [optics] block with enabled = true to compute optical properties'
    stop 1
  end if

  ! ================================================================
  ! Branch by confinement type
  ! ================================================================
  select case (trim(cfg%confinement))

  ! ==================================================================
  ! BULK (confinement=0)
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

    ! LAPACK workspace query parameters
    vl = 0.0_dp
    vu = 0.0_dp
    NB = ILAENV(1, 'ZHETRD', 'UPLO', N, N, -1, -1)
    NB = MAX(NB, N)
    ABSTOL = DLAMCH('P')

    ! Allocate workspace
    allocate(rwork(7*N))
    allocate(iwork(5*N))
    allocate(ifail(N))
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
    ! Workspace query (k=1)
    ! ================================================================
    call ZB8bandBulk(HT, smallk(1), cfg%params(1:1), cfg=cfg)
    allocate(work(1))
    lwork = -1
    call zheevx('V', 'I', 'U', N, HT, N, vl, vu, il, iuu, abstol, M, &
      eig(:,1), HT, N, work, lwork, rwork, iwork, ifail, info)
    if (info /= 0) then
      print '(a,i0)', 'ERROR: zheevx workspace query failed, info = ', info
      stop 1
    end if
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))

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
      complex(kind=dp), allocatable :: HT_loc(:,:), work_loc(:)
      real(kind=dp), allocatable    :: rwork_loc(:)
      integer, allocatable          :: iwork_loc(:), ifail_loc(:)
      integer :: info_loc, M_loc

      !$omp parallel private(k, HT_loc, work_loc, rwork_loc, iwork_loc, ifail_loc, info_loc, M_loc)
      allocate(HT_loc(N, N))
      allocate(work_loc(lwork))
      allocate(rwork_loc(7*N))
      allocate(iwork_loc(5*N))
      allocate(ifail_loc(N))
      HT_loc = (0.0_dp, 0.0_dp)

      !$omp do schedule(static)
      do k = 1, npts
        ! Build 8x8 bulk Hamiltonian at this k
        call ZB8bandBulk(HT_loc, smallk(k), cfg%params(1:1), cfg=cfg)

        ! Diagonalize all 8 eigenvalues
        call zheevx('V', 'I', 'U', N, HT_loc, N, vl, vu, il, iuu, abstol, M_loc, &
          eig(:,k), HT_loc, N, work_loc, lwork, rwork_loc, iwork_loc, &
          ifail_loc, info_loc)
        if (info_loc /= 0) then
          !$omp critical
          print '(a,i0,a,i0)', ' ERROR: diagonalization failed at k=', k, ' info=', info_loc
          !$omp end critical
          stop 1
        end if

        ! Store eigenvectors (zheevx overwrites HT_loc with them)
        eigv(:,:,k) = HT_loc(:, 1:N)
      end do
      !$omp end do

      deallocate(HT_loc, work_loc, rwork_loc, iwork_loc, ifail_loc)
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
  ! QUANTUM WELL (confinement=1)
  ! ==================================================================
  case ('qw')
    if (cfg%conf_dir /= 'z') then
      print '(a)', 'ERROR: opticalProperties requires confDir=z for QW'
      stop 1
    end if

    call simulation_setup_init(cfg, setup)
    N = setup%N
    if (setup%sc_was_run) cfg%sc%fermi_level = setup%fermi_level
    ! Compute eigenvalue range for optics: top numvb valence + bottom numcb conduction
    il = NUM_VB_STATES*cfg%ngrid - cfg%bands%num_vb + 1
    iuu = NUM_VB_STATES*cfg%ngrid + cfg%bands%num_cb
    print '(a,i0,a,i0)', ' Computing states from index ', il, ' to ', iuu

    ! Build velocity matrices via simulation_setup
    call setup_build_velocity_matrices(setup, cfg)
    print '(a)', ' Velocity matrices built successfully'

    ! LAPACK workspace query parameters
    vl = 0.0_dp
    vu = 0.0_dp
    NB = ILAENV(1, 'ZHETRD', 'UPLO', N, N, -1, -1)
    NB = MAX(NB, N)
    ABSTOL = DLAMCH('P')

    ! Allocate workspace
    allocate(rwork(7*N))
    allocate(iwork(5*N))
    allocate(ifail(N))
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
    do i = 1, cfg%ngrid
      write(iounit, *) cfg%z(i), setup%profile(i,1), setup%profile(i,2), setup%profile(i,3)
    end do
    close(iounit)

    ! ================================================================
    ! Workspace query (k=1)
    ! ================================================================
    call ZB8bandQW(HT, smallk(1), setup%profile, setup%kpterms, cfg=cfg)
    allocate(work(1))
    lwork = -1
    call zheevx('V', 'I', 'U', N, HT, N, vl, vu, il, iuu, abstol, M, &
      eig(:,1), HT, N, work, lwork, rwork, iwork, ifail, info)
    if (info /= 0) then
      print '(a,i0)', 'ERROR: zheevx workspace query failed, info = ', info
      stop 1
    end if
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))

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
    ! Disable MKL internal threading (serial zheevx per k-point)
    info = mkl_set_num_threads_local(1)

    print '(a,i0,a)', ' QW optics k-sweep: ', npts, ' k-points'

    block
      ! Thread-private temporaries for OpenMP parallel region
      complex(kind=dp), allocatable :: HT_loc(:,:), work_loc(:)
      real(kind=dp), allocatable    :: rwork_loc(:)
      integer, allocatable          :: iwork_loc(:), ifail_loc(:)
      integer :: info_loc, M_loc

      !$omp parallel private(k, HT_loc, work_loc, rwork_loc, iwork_loc, ifail_loc, info_loc, M_loc)
      allocate(HT_loc(N, N))
      allocate(work_loc(lwork))
      allocate(rwork_loc(7*N))
      allocate(iwork_loc(5*N))
      allocate(ifail_loc(N))
      HT_loc = (0.0_dp, 0.0_dp)

      !$omp do schedule(static)
      do k = 1, npts
        ! Build Hamiltonian for this k_par
        call ZB8bandQW(HT_loc, smallk(k), setup%profile, setup%kpterms, cfg=cfg)

        ! Diagonalize
        call zheevx('V', 'I', 'U', N, HT_loc, N, vl, vu, il, iuu, abstol, M_loc, &
          eig(:,k), HT_loc, N, work_loc, lwork, rwork_loc, iwork_loc, &
          ifail_loc, info_loc)
        if (info_loc /= 0) then
          !$omp critical
          print '(a,i0,a,i0)', ' ERROR: diagonalization failed at k=', k, ' info=', info_loc
          !$omp end critical
          stop 1
        end if

        ! Store eigenvectors (zheevx overwrites HT_loc with them)
        eigv(:,:,k) = HT_loc(:, 1:iuu-il+1)
      end do
      !$omp end do

      deallocate(HT_loc, work_loc, rwork_loc, iwork_loc, ifail_loc)
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

    ! ISBT dipole transition table (zone-center only)
    if (cfg%optics%isbt_enabled) then
      call compute_intersubband_transitions(eig(:, 1), eigv(:, :, 1), &
        & cfg%z, cfg%dz, cfg%bands%num_cb, cfg%bands%num_vb, cfg%ngrid, &
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
          & cfg%bands%num_cb, cfg%bands%num_vb, cfg%ngrid, E_binding_ex, lambda_opt_ex, &
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
  ! WIRE (confinement=2)
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
        call setup%eigen_solver%solve(setup%HT_csr_ptr, setup%eigen_cfg, eigen_res)

        if (eigen_res%nev_found == 0) then
          print '(a,i0,a)', ' WARNING: FEAST found no eigenvalues at kz-point ', k
          call eigensolver_result_free(eigen_res)
          cycle
        end if

        if (eigen_res%nev_found < cfg%bands%num_cb + cfg%bands%num_vb) then
          print '(a,i0,a,i0,a,i0)', ' WARNING: FEAST found only ', &
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
    print '(a,i0)', 'ERROR: unknown confinement mode ', cfg%confinement
    stop 1

  end select

  ! ================================================================
  ! Cleanup
  ! ================================================================
  if (allocated(smallk))  deallocate(smallk)
  if (allocated(HT))      deallocate(HT)
  if (allocated(eig))     deallocate(eig)
  if (allocated(eigv))    deallocate(eigv)
  if (allocated(work))    deallocate(work)
  if (allocated(iwork))   deallocate(iwork)
  if (allocated(rwork))   deallocate(rwork)
  if (allocated(ifail))   deallocate(ifail)

end program opticalProperties
