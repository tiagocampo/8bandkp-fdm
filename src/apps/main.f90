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
  integer :: info, ILAENV, NB, lwork, N, M, il, iuu, vl, vu
  real(kind=dp) :: abstol, dlamch
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
  type(eigensolver_config)         :: eigen_cfg
  type(eigensolver_result)         :: eigen_res
  real(kind=dp), allocatable       :: eig_wire(:,:)
  integer                          :: Ngrid, Ntot, nev_wire

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
    print *, '  Grid: ny=', cfg%grid%ny, ' nz=', cfg%grid%nz, ' Ngrid=', Ngrid
    print *, '  Matrix size: ', Ntot, 'x', Ntot
    print *, '  Requesting ', nev_wire, ' eigenvalues per k-point'
    print *, '  kx sweep: ', cfg%waveVectorStep, ' points, kx_max=', cfg%waveVectorMax
    print *, ''

    ! Set up eigensolver configuration (FEAST with energy window)
    eigen_cfg%method   = 'FEAST'
    eigen_cfg%nev      = nev_wire
    eigen_cfg%emin     = -1.0_dp   ! placeholder, set by auto_compute_energy_window
    eigen_cfg%emax     =  1.0_dp   ! placeholder
    eigen_cfg%max_iter = 100
    eigen_cfg%tol      = 1.0e-10_dp
    eigen_cfg%feast_m0 = 0         ! auto: 2*nev

    ! Allocate storage for eigenvalues across k-points
    allocate(eig_wire(nev_wire, cfg%waveVectorStep))
    eig_wire = 0.0_dp

    ! kx sweep loop
    do k = 1, cfg%waveVectorStep

      print *, 'k-point ', k, '/', cfg%waveVectorStep, &
        & ' kx=', smallk(k)%kx

      ! Build sparse Hamiltonian at this kx
      call ZB8bandGeneralized(HT_csr, smallk(k)%kx, profile_2d, &
        & kpterms_2d, cfg)

      ! Auto-compute energy window from Gershgorin bounds on first k-point
      if (k == 1) then
        call auto_compute_energy_window(HT_csr, eigen_cfg%emin, eigen_cfg%emax)
        print *, '  Auto energy window: [', eigen_cfg%emin, ',', eigen_cfg%emax, ']'
      end if

      ! Solve sparse eigenvalue problem
      call solve_sparse_evp(HT_csr, eigen_cfg, eigen_res)

      if (.not. eigen_res%converged) then
        print *, '  WARNING: FEAST did not converge at k-point ', k
      end if

      print *, '  Found ', eigen_res%nev_found, ' eigenvalues (requested ', &
        & nev_wire, ') iterations=', eigen_res%iterations

      ! Store eigenvalues (pad with zero if fewer found than requested)
      if (eigen_res%nev_found > 0) then
        do i = 1, min(eigen_res%nev_found, nev_wire)
          eig_wire(i, k) = eigen_res%eigenvalues(i)
        end do
      end if

      ! Clean up result and Hamiltonian for this k-point
      call eigensolver_result_free(eigen_res)
      call csr_free(HT_csr)

    end do

    ! Write eigenvalues to file
    call writeEigenvalues(smallk, eig_wire, cfg%waveVectorStep)
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

  ! --- Self-consistent loop (QW only) ---
  if (cfg%confDir == 'z' .and. cfg%sc%enabled == 1) then
    print *, ''
    print *, '=== Running self-consistent Schrodinger-Poisson loop ==='
    call self_consistent_loop(profile, cfg, kpterms, HT, eig, eigv, &
      & smallk, N, il, iuu)

    ! Write updated profile after SC convergence
    call get_unit(iounit)
    open(unit=iounit, file='output/sc_potential_profile.dat', status="replace", action="write")
    do i = 1, cfg%fdStep, 1
      write(iounit,*) cfg%z(i), profile(i,1), profile(i,2), profile(i,3)
    end do
    close(iounit)
    print *, 'SC potential profile written to output/sc_potential_profile.dat'
  end if

  do k = 1, cfg%waveVectorStep, 1

    print *, smallk(k)%kx, smallk(k)%ky, smallk(k)%kz
    if (cfg%confDir == 'n') then !BULK
      call ZB8bandBulk(HT,smallk(k),cfg%params(1))
      ! Copy to temporary array
      HTmp = HT(1:8,1:8)

      ! Query optimal workspace
      if (k == 1) then
        if (allocated(work)) deallocate(work)
        allocate(work(1))
        lwork = -1
        call zheevx('V', 'I', 'U', 8, HTmp, 8, vl, vu, il, iuu, abstol, M, eig(1:cfg%evnum,k), &
                   HTmp, 8, work, lwork, rwork, iwork, ifail, info)
        if (info /= 0) then
          print *, 'Error: zheevx bulk workspace query failed, info =', info
          stop 1
        end if
        lwork = int(real(work(1)))
        if (allocated(work)) deallocate(work)
        allocate(work(lwork))
      end if

      ! Actual diagonalization - use zheevx for selected eigenvalues
      call zheevx('V', 'I', 'U', 8, HTmp, 8, vl, vu, il, iuu, abstol, M, eig(1:cfg%evnum,k), &
                 HTmp, 8, work, lwork, rwork, iwork, ifail, info)
      if (info /= 0) then
        print *, "Diagonalization error in bulk calculation, info = ", info
        if (info < 0) then
          print *, "Parameter ", -info, " had illegal value"
        end if
        stop "error diag"
      end if

      ! For bulk, just copy the selected eigenvectors
      eigv(:,:,k) = HTmp(:,1:cfg%evnum)

    else if (cfg%confDir == 'z') then ! QUANTUM WELL
      call ZB8bandQW(HT, smallk(k), profile, kpterms)

      ! Query optimal workspace for full matrix
      if (k == 1) then
        if (allocated(work)) deallocate(work)
        allocate(work(1))
        lwork = -1
        call zheevx('V', 'I', 'U', N, HT, N, vl, vu, il, iuu, abstol, M, eig(:,k), &
                   HT, N, work, lwork, rwork, iwork, ifail, info)
        if (info /= 0) then
          print *, 'Error: zheevx QW workspace query failed, info =', info
          stop 1
        end if
        lwork = int(real(work(1)))
        if (allocated(work)) deallocate(work)
        allocate(work(lwork))
      end if

      ! Actual diagonalization of full matrix using zheevx
      call zheevx('V', 'I', 'U', N, HT, N, vl, vu, il, iuu, abstol, M, eig(:,k), &
                 HT, N, work, lwork, rwork, iwork, ifail, info)
      if (info /= 0) then
        print *, "Diagonalization error in quantum well calculation, info = ", info
        if (info < 0) then
          print *, "Parameter ", -info, " had illegal value"
        end if
        stop "error diag"
      end if

      ! For quantum well, copy only the selected eigenvectors
      eigv(:,:,k) = HT(:,1:iuu-il+1)
    else
      stop "verify confinement direction or something else"
    end if

    ! Write eigenfunctions only at k=0 or specific k-points
    if (k == 1 .or. k == int(cfg%waveVectorStep/2) .or. k == cfg%waveVectorStep) then  ! Write at start, middle, and end k-points
      if (cfg%confDir == 'n') then
        call writeEigenfunctions(8, min(cfg%evnum,8), HTmp(:,1:min(cfg%evnum,8)), k, cfg%fdstep, cfg%z, .true.)
      else
        ! For quantum well, write only the requested number of states
        call writeEigenfunctions(N, iuu-il+1, HT(:,1:iuu-il+1), k, cfg%fdstep, cfg%z, .false.)
      end if
    end if

  end do
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
