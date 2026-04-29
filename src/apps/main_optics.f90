program opticalProperties

  use definitions
  use parameters
  use hamiltonianConstructor
  use hamiltonian_wire, only: wire_workspace, wire_workspace_free, &
    & build_velocity_matrices, ZB8bandGeneralized
  use confinement_init, only: confinementInitialization_2d
  use input_parser
  use optical_spectra
  use exciton_solver, only: compute_exciton_binding, apply_excitonic_corrections
  use sparse_matrices
  use eigensolver, only: eigensolver_base, make_eigensolver, eigensolver_config, &
    & eigensolver_result, eigensolver_result_free, &
    & auto_compute_energy_window
  use strain_solver
  use linalg, only: zheevx, ilaenv, dlamch, mkl_set_num_threads_local
  use utils, only: dnscsr_z_mkl
  use outputFunctions, only: ensure_output_dir, get_unit

  implicit none

  ! Shared configuration from input_parser
  type(simulation_config) :: cfg
  real(kind=dp), allocatable, dimension(:,:) :: profile
  real(kind=dp), allocatable, dimension(:,:,:) :: kpterms

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

  ! Velocity matrices
  type(csr_matrix) :: vel_opt(3)
  type(csr_matrix) :: H_csr_tmp
  complex(kind=dp), allocatable :: H_k0(:,:)
  type(wavevector) :: pert_kx, pert_ky
  integer :: nzmax_tmp

  ! Wire-specific variables (confinement=2)
  real(kind=dp), allocatable       :: profile_2d(:,:)
  type(csr_matrix), allocatable    :: kpterms_2d(:)
  type(csr_matrix)                 :: HT_csr
  type(wire_workspace)             :: wire_ws
  class(eigensolver_base), allocatable :: eigen_solver
  type(eigensolver_config)         :: eigen_cfg
  type(eigensolver_result)         :: eigen_res
  type(csr_matrix)                 :: vel_wire(3)
  integer                          :: Ngrid, Ntot, nev_wire

  ! File handling
  integer(kind=4) :: iounit

  ! ================================================================
  ! Parse input, initialize materials, external fields
  ! ================================================================
  call read_and_setup(cfg, profile, kpterms)

  ! Check that optics is enabled
  if (.not. cfg%optics%enabled) then
    print '(a)', 'ERROR: optics block not enabled in input.cfg'
    print '(a)', '  Add an optics: block with enabled = T to compute optical properties'
    stop 1
  end if

  ! ================================================================
  ! Branch by confinement type
  ! ================================================================
  select case (cfg%confinement)

  ! ==================================================================
  ! BULK (confinement=0)
  ! ==================================================================
  case (0)
    N = 8  ! 8x8 bulk Hamiltonian

    ! For bulk: 6 VB-like + 2 CB-like states
    il = 1
    iuu = N
    print '(a)', ' Bulk optics: 8x8 Hamiltonian, all 8 states'

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
    allocate(eig(iuu - il + 1, cfg%waveVectorStep))
    allocate(eigv(N, iuu - il + 1, cfg%waveVectorStep))
    eig(:,:) = 0.0_dp
    eigv(:,:,:) = (0.0_dp, 0.0_dp)

    allocate(HT(N, N))
    HT = (0.0_dp, 0.0_dp)

    ! Build k-sweep array: 1D sweep from 0 to waveVectorMax
    npts = cfg%waveVectorStep
    allocate(smallk(npts))
    smallk%kx = 0.0_dp
    smallk%ky = 0.0_dp
    smallk%kz = 0.0_dp
    do k = 1, npts
      smallk(k)%kz = real(k - 1, kind=dp) * cfg%waveVectorMax &
        & / real(npts - 1, kind=dp)
    end do

    ! ================================================================
    ! Workspace query (k=1)
    ! ================================================================
    call ZB8bandBulk(HT, smallk(1), cfg%params(1:1))
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
    ! Build velocity matrices at unit wavevector (k-independent)
    ! ================================================================
    ! For bulk, ZB8bandBulk(g='g') zeros gamma and A, keeping only
    ! the Kane P term. The resulting matrix is dH/dk_alpha which is
    ! k-independent (linear in k-derivative).

    ! vel(1): x-direction
    pert_kx%kx = 1.0_dp
    pert_kx%ky = 0.0_dp
    pert_kx%kz = 0.0_dp
    allocate(H_k0(N, N))
    H_k0 = (0.0_dp, 0.0_dp)
    call ZB8bandBulk(H_k0, pert_kx, cfg%params(1:1), g='g')
    nzmax_tmp = N * N
    call dnscsr_z_mkl(nzmax_tmp, N, H_k0, vel_opt(1))
    deallocate(H_k0)

    ! vel(2): y-direction
    pert_ky%kx = 0.0_dp
    pert_ky%ky = 1.0_dp
    pert_ky%kz = 0.0_dp
    allocate(H_k0(N, N))
    H_k0 = (0.0_dp, 0.0_dp)
    call ZB8bandBulk(H_k0, pert_ky, cfg%params(1:1), g='g')
    nzmax_tmp = N * N
    call dnscsr_z_mkl(nzmax_tmp, N, H_k0, vel_opt(2))
    deallocate(H_k0)

    ! vel(3): z-direction
    block
      type(wavevector) :: pert_kz
      pert_kz%kx = 0.0_dp
      pert_kz%ky = 0.0_dp
      pert_kz%kz = 1.0_dp
      allocate(H_k0(N, N))
      H_k0 = (0.0_dp, 0.0_dp)
      call ZB8bandBulk(H_k0, pert_kz, cfg%params(1:1), g='g')
      nzmax_tmp = N * N
      call dnscsr_z_mkl(nzmax_tmp, N, H_k0, vel_opt(3))
      deallocate(H_k0)
    end block

    print '(a)', ' Bulk velocity matrices built successfully'

    ! ================================================================
    ! Initialize optics accumulation
    ! ================================================================
    call optics_init(cfg%optics)

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
    dk = cfg%waveVectorMax / real(cfg%waveVectorStep - 1, kind=dp)
    do i = 1, npts
      ! Base Simpson 1/3 rule weight
      if (i == 1 .or. i == npts) then
        simpson_w(i) = dk / 3.0_dp
      else if (mod(i, 2) == 0) then
        simpson_w(i) = 4.0_dp * dk / 3.0_dp
      else
        simpson_w(i) = 2.0_dp * dk / 3.0_dp
      end if
      ! Multiply by 4*pi*k^2 / (2*pi)^3 for 3D spherical integration
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
        call ZB8bandBulk(HT_loc, smallk(k), cfg%params(1:1))

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

    ! Accumulate optical spectra (serial, uses shared vel_opt)
    print '(a)', ' Accumulating optical spectra...'
    do k = 1, npts
      call optics_accumulate(cfg%optics, &
        & eig(:, k), eigv(:, :, k), simpson_w(k), &
        & vel_opt, cfg%numcb, cfg%numvb, cfg%sc%fermi_level)

      if (cfg%optics%spontaneous_enabled) then
        call optics_accumulate_spontaneous(cfg%optics, &
          & eig(:, k), eigv(:, :, k), simpson_w(k), &
          & vel_opt, cfg%numcb, cfg%numvb, cfg%sc%fermi_level)
      end if

      if (cfg%optics%gain_enabled) then
        call compute_gain_qw(cfg%optics, &
          & eig(:, k), eigv(:, :, k), simpson_w(k), &
          & vel_opt, cfg%numcb, cfg%numvb, &
          & cfg%optics%gain_carrier_density)
      end if

      if (cfg%optics%isbt_enabled) then
        call compute_isbt_absorption(cfg%optics, &
          & eig(:, k), eigv(:, :, k), &
          & vel_opt, cfg%numcb, cfg%numvb, &
          & simpson_w(k), cfg%sc%fermi_level)
      end if
    end do

    ! ================================================================
    ! Finalize: apply prefactor, write output files
    ! ================================================================
    call optics_finalize(cfg%optics)
    call optics_cleanup()

    ! Free velocity matrices
    do i = 1, 3
      call csr_free(vel_opt(i))
    end do

    deallocate(simpson_w)
    print '(a)', ' Bulk optical spectra written to output/'

  ! ==================================================================
  ! QUANTUM WELL (confinement=1)
  ! ==================================================================
  case (1)
    if (cfg%confDir /= 'z') then
      print '(a)', 'ERROR: opticalProperties requires confDir=z for QW'
      stop 1
    end if

    N = cfg%fdStep * 8

    ! Set eigenvalue range: highest numvb valence + lowest numcb conduction
    il = NUM_VB_STATES*cfg%fdStep - cfg%numvb + 1
    iuu = NUM_VB_STATES*cfg%fdStep + cfg%numcb
    print '(a,i0,a,i0)', ' Computing states from index ', il, ' to ', iuu

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
    allocate(eig(iuu-il+1, cfg%waveVectorStep))
    allocate(eigv(N, iuu-il+1, cfg%waveVectorStep))
    eig(:,:) = 0.0_dp
    eigv(:,:,:) = (0.0_dp, 0.0_dp)

    allocate(HT(N, N))
    HT = (0.0_dp, 0.0_dp)

    ! Build k_par array: uniform sweep from 0 to waveVectorMax
    npts = cfg%waveVectorStep
    allocate(smallk(npts))
    smallk%kx = 0.0_dp
    smallk%ky = 0.0_dp
    smallk%kz = 0.0_dp
    do k = 1, npts
      smallk(k)%kx = real(k - 1, kind=dp) * cfg%waveVectorMax &
        & / real(npts - 1, kind=dp)
    end do

    ! --- Write potential profile ---
    call ensure_output_dir()
    call get_unit(iounit)
    open(unit=iounit, file='output/potential_profile.dat', status='replace', action='write')
    do i = 1, cfg%fdStep
      write(iounit, *) cfg%z(i), profile(i,1), profile(i,2), profile(i,3)
    end do
    close(iounit)

    ! ================================================================
    ! Workspace query (k=1)
    ! ================================================================
    call ZB8bandQW(HT, smallk(1), profile, kpterms, cfg=cfg)
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
    ! Build velocity matrices at k=0 (k-independent, built once)
    ! ================================================================
    ! vel(3): z-direction via commutator -i[r_z, H] on CSR
    allocate(H_k0(N, N))
    H_k0 = (0.0_dp, 0.0_dp)
    call ZB8bandQW(H_k0, smallk(1), profile, kpterms, cfg=cfg)
    nzmax_tmp = N * N
    call dnscsr_z_mkl(nzmax_tmp, N, H_k0, H_csr_tmp)
    deallocate(H_k0)

    ! build_velocity_matrices dispatches to the 1D version for QW
    call build_velocity_matrices(H_csr_tmp, cfg%grid, vel_opt)
    call csr_free(H_csr_tmp)

    ! vel(1): x-direction via dH/dk_x (Kane P matrix elements)
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

    ! vel(2): y-direction via dH/dk_y (Kane P matrix elements)
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

    print '(a)', ' Velocity matrices built successfully'

    ! ================================================================
    ! Initialize optics accumulation
    ! ================================================================
    call optics_init(cfg%optics)

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
    dk = cfg%waveVectorMax / real(cfg%waveVectorStep - 1, kind=dp)
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
        call ZB8bandQW(HT_loc, smallk(k), profile, kpterms, cfg=cfg)

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

    ! Accumulate optical spectra (serial, uses shared vel_opt)
    print '(a)', ' Accumulating optical spectra...'
    do k = 1, npts
      call optics_accumulate(cfg%optics, &
        & eig(:, k), eigv(:, :, k), simpson_w(k), &
        & vel_opt, cfg%numcb, cfg%numvb, cfg%sc%fermi_level)

      if (cfg%optics%spontaneous_enabled) then
        call optics_accumulate_spontaneous(cfg%optics, &
          & eig(:, k), eigv(:, :, k), simpson_w(k), &
          & vel_opt, cfg%numcb, cfg%numvb, cfg%sc%fermi_level)
      end if

      if (cfg%optics%gain_enabled) then
        call compute_gain_qw(cfg%optics, &
          & eig(:, k), eigv(:, :, k), simpson_w(k), &
          & vel_opt, cfg%numcb, cfg%numvb, &
          & cfg%optics%gain_carrier_density)
      end if

      if (cfg%optics%isbt_enabled) then
        call compute_isbt_absorption(cfg%optics, &
          & eig(:, k), eigv(:, :, k), &
          & vel_opt, cfg%numcb, cfg%numvb, &
          & simpson_w(k), cfg%sc%fermi_level)
      end if
    end do

    ! ================================================================
    ! Finalize: apply prefactor, write output files
    ! ================================================================
    call optics_finalize(cfg%optics)

    ! Exciton corrections (applied after finalize, before cleanup)
    if (cfg%exciton%enabled) then
      block
        real(kind=dp) :: E_binding_ex, lambda_opt_ex, E_gap_ex
        integer :: ie, iounit_ex, cb_st, vb_st
        cb_st = cfg%numvb + 1
        vb_st = cfg%numvb
        E_gap_ex = eig(cb_st, 1) - eig(vb_st, 1)
        call compute_exciton_binding(eig(:, 1), eigv(:, :, 1), &
          & cfg%z, cfg%dz, cfg%numLayers, cfg%params, &
          & cfg%numcb, cfg%numvb, cfg%fdstep, E_binding_ex, lambda_opt_ex, &
          & cfg%grid%material_id)
        print '(a,f8.3,a)', ' Exciton binding energy: ', E_binding_ex, ' meV'
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
        print '(a)', ' Excitonic corrections applied to absorption spectra'
      end block
    end if

    call optics_cleanup()

    ! Free velocity matrices
    do i = 1, 3
      call csr_free(vel_opt(i))
    end do

    deallocate(simpson_w)
    print '(a)', ' Optical spectra written to output/'

  ! ==================================================================
  ! WIRE (confinement=2)
  ! ==================================================================
  case (2)

    ! ----------------------------------------------------------------
    ! Initialize 2D confinement operators (sparse CSR kpterms)
    ! ----------------------------------------------------------------
    call confinementInitialization_2d(cfg%grid, cfg%params, cfg%regions, &
      & profile_2d, kpterms_2d, cfg%FDorder)

    ! Optional strain
    if (cfg%strain%enabled) then
      block
        type(strain_result) :: strain_out
        real(kind=dp) :: a0_ref

        a0_ref = cfg%params(1)%a0
        call compute_strain(cfg%grid, cfg%params, cfg%grid%material_id, &
          cfg%strain, a0_ref, strain_out)
        call compute_bir_pikus_blocks(strain_out, cfg%params, &
          cfg%grid%material_id, cfg%grid, cfg%strain_blocks)
        call strain_result_free(strain_out)
        print '(a)', '  Wire strain calculation complete'
      end block
    end if

    Ngrid = grid_ngrid(cfg%grid)
    Ntot  = 8 * Ngrid
    nev_wire = cfg%numcb + cfg%numvb
    if (nev_wire > Ntot) then
      print '(a,i0,a,i0)', ' Warning: requesting ', nev_wire, &
        & ' bands but matrix size is ', Ntot
      nev_wire = Ntot
    end if

    print '(a)', ''
    print '(a)', '=== Wire optical properties (2D confinement) ==='
    print '(a,i0,a,i0,a,i0)', '  Grid: nx=', cfg%grid%nx, ' ny=', cfg%grid%ny, &
      & ' Ngrid=', Ngrid
    print '(a,i0,a,i0)', '  Matrix size: ', Ntot, 'x', Ntot
    print '(a,i0,a,i0)', '  numcb=', cfg%numcb, ' numvb=', cfg%numvb
    print '(a)', ''

    ! ----------------------------------------------------------------
    ! Build wire Hamiltonian at kz=0
    ! ----------------------------------------------------------------
    call ZB8bandGeneralized(HT_csr, 0.0_dp, profile_2d, kpterms_2d, &
      & cfg, ws=wire_ws)

    ! Build commutator-based velocity matrices for x,y directions
    call build_velocity_matrices(HT_csr, cfg%grid, vel_wire(1), vel_wire(2))

    ! Build velocity matrix for z direction (free axis, uses dH/dkz)
    call ZB8bandGeneralized(vel_wire(3), 1.0_dp, profile_2d, kpterms_2d, &
      & cfg, g='g3')

    ! Keep HT_csr alive: the COO cache fast path in the kz sweep reuses its CSR structure

    print '(a)', ' Wire velocity matrices built successfully'

    ! ----------------------------------------------------------------
    ! Configure FEAST eigensolver
    ! ----------------------------------------------------------------
    eigen_cfg%method   = 'FEAST'
    eigen_cfg%nev      = nev_wire
    eigen_cfg%max_iter = 100
    eigen_cfg%tol      = 1.0e-10_dp
    eigen_cfg%feast_m0 = cfg%feast_m0

    ! Set energy window once (before the kz sweep)
    if (cfg%feast_emin /= 0.0_dp .and. cfg%feast_emax /= 0.0_dp) then
      eigen_cfg%emin = cfg%feast_emin
      eigen_cfg%emax = cfg%feast_emax
    else
      call auto_compute_energy_window(HT_csr, eigen_cfg%emin, eigen_cfg%emax)
    end if

    ! Create polymorphic eigensolver
    eigen_solver = make_eigensolver(eigen_cfg)

    ! ----------------------------------------------------------------
    ! Initialize optics accumulation
    ! ----------------------------------------------------------------
    call optics_init(cfg%optics)

    ! ----------------------------------------------------------------
    ! Simpson integration weights for 1D kz integration
    ! int dkz  =>  weight = Simpson * 1 / (2*pi)
    ! ----------------------------------------------------------------
    npts = cfg%waveVectorStep
    if (mod(npts, 2) == 0) then
      npts = npts - 1  ! Simpson requires odd number of points
      print '(a,i0,a)', ' Warning: optics integration uses ', npts, &
        & ' k-points (Simpson requires odd count)'
    end if
    allocate(simpson_w(npts))
    dk = cfg%waveVectorMax / real(cfg%waveVectorStep - 1, kind=dp)
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
    ! The COO cache preserves the CSR sparsity pattern across kz values,
    ! so only the values are updated on subsequent calls (fast path).
    ! We keep HT_csr allocated throughout the sweep.
    ! ----------------------------------------------------------------
    print '(a,i0,a)', ' Wire optics kz-sweep: ', npts, ' k-points'

    do k = 1, npts
      ! Build kz value for this sweep point
      block
        real(kind=dp) :: kz_val

        kz_val = real(k - 1, kind=dp) * dk

        ! Build wire Hamiltonian at this kz (cache reuses CSR structure)
        call ZB8bandGeneralized(HT_csr, kz_val, profile_2d, kpterms_2d, &
          & cfg, ws=wire_ws)

        ! Solve eigenvalue problem
        call eigen_solver%solve(HT_csr, eigen_cfg, eigen_res)

        if (eigen_res%nev_found == 0) then
          print '(a,i0,a)', ' WARNING: FEAST found no eigenvalues at kz-point ', k
          call eigensolver_result_free(eigen_res)
          cycle
        end if

        if (eigen_res%nev_found < cfg%numcb + cfg%numvb) then
          print '(a,i0,a,i0,a,i0)', ' WARNING: FEAST found only ', &
            eigen_res%nev_found, ' eigenvalues at kz-point ', k, &
            ' (need ', cfg%numcb + cfg%numvb, ')'
          call eigensolver_result_free(eigen_res)
          cycle
        end if

        ! Accumulate optical spectra for this kz
        call optics_accumulate(cfg%optics, &
          & eigen_res%eigenvalues(1:nev_wire), &
          & eigen_res%eigenvectors(:, 1:nev_wire), &
          & simpson_w(k), &
          & vel_wire, cfg%numcb, cfg%numvb, cfg%sc%fermi_level)

        if (cfg%optics%spontaneous_enabled) then
          call optics_accumulate_spontaneous(cfg%optics, &
            & eigen_res%eigenvalues(1:nev_wire), &
            & eigen_res%eigenvectors(:, 1:nev_wire), &
            & simpson_w(k), &
            & vel_wire, cfg%numcb, cfg%numvb, cfg%sc%fermi_level)
        end if

        if (cfg%optics%gain_enabled) then
          call compute_gain_qw(cfg%optics, &
            & eigen_res%eigenvalues(1:nev_wire), &
            & eigen_res%eigenvectors(:, 1:nev_wire), &
            & simpson_w(k), &
            & vel_wire, cfg%numcb, cfg%numvb, &
            & cfg%optics%gain_carrier_density)
        end if

        if (cfg%optics%isbt_enabled) then
          call compute_isbt_absorption(cfg%optics, &
            & eigen_res%eigenvalues(1:nev_wire), &
            & eigen_res%eigenvectors(:, 1:nev_wire), &
            & vel_wire, cfg%numcb, cfg%numvb, &
            & simpson_w(k), cfg%sc%fermi_level)
        end if

        ! Free eigensolver result for next kz
        call eigensolver_result_free(eigen_res)
      end block
    end do

    ! Free Hamiltonian CSR after sweep completes
    call csr_free(HT_csr)

    ! ----------------------------------------------------------------
    ! Finalize: apply prefactor, write output files
    ! ----------------------------------------------------------------
    call optics_finalize(cfg%optics)
    call optics_cleanup()

    ! Free velocity matrices
    do i = 1, 3
      call csr_free(vel_wire(i))
    end do

    ! Free wire-specific resources
    call wire_workspace_free(wire_ws)
    if (allocated(eigen_solver)) deallocate(eigen_solver)
    if (allocated(profile_2d)) deallocate(profile_2d)
    if (allocated(kpterms_2d)) then
      do i = 1, size(kpterms_2d)
        call csr_free(kpterms_2d(i))
      end do
      deallocate(kpterms_2d)
    end if

    deallocate(simpson_w)
    print '(a)', ' Wire optical spectra written to output/'

  case default
    print '(a,i0)', 'ERROR: unknown confinement mode ', cfg%confinement
    stop 1

  end select

  ! ================================================================
  ! Cleanup
  ! ================================================================
  if (allocated(smallk))  deallocate(smallk)
  if (allocated(profile)) deallocate(profile)
  if (allocated(kpterms)) deallocate(kpterms)
  if (allocated(HT))      deallocate(HT)
  if (allocated(eig))     deallocate(eig)
  if (allocated(eigv))    deallocate(eigv)
  if (allocated(work))    deallocate(work)
  if (allocated(iwork))   deallocate(iwork)
  if (allocated(rwork))   deallocate(rwork)
  if (allocated(ifail))   deallocate(ifail)

end program opticalProperties
