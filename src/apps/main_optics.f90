program opticalProperties

  use definitions
  use parameters
  use hamiltonianConstructor
  use input_parser
  use optical_spectra
  use sparse_matrices
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
    print '(a)', ' Optical spectra written to output/'

  ! ==================================================================
  ! WIRE (confinement=2)
  ! ==================================================================
  case (2)
    print '(a)', 'Wire optics: TODO (Phase 3)'
    stop 0

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
