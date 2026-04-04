program gfactor

  use definitions
  use parameters
  use hamiltonianConstructor
  use finitedifferences
  use OMP_lib
  use outputFunctions
  use gfactorFunctions
  use utils
  use input_parser
  use sparse_matrices
  use eigensolver
  use strain_solver

  implicit NONE

  ! Shared configuration from input_parser
  type(simulation_config) :: cfg
  real(kind=dp), allocatable, dimension(:,:) :: profile
  real(kind=dp), allocatable, dimension(:,:,:) :: kpterms

  ! wave vector
  type(wavevector), allocatable, dimension(:) :: smallk

  ! iteration consts
  integer :: i, j, k

  ! hamiltonian and LAPACK/BLAS
  integer :: info, ILAENV, NB, lwork, N, M, il, iuu, vl, vu, lrwork, liwork
  real(kind=dp) :: abstol, dlamch
  real(kind=dp), allocatable :: eig(:,:), rwork(:)
  complex(kind=dp), allocatable :: work(:)
  complex(kind=dp), allocatable, dimension(:,:,:) :: eigv
  integer, allocatable :: iwork(:), ifail(:)
  complex(kind=dp), allocatable, dimension(:,:) :: HT

  ! gfactor
  complex(kind=dp), allocatable, dimension(:,:) :: cb_state, vb_state
  real(kind=dp), allocatable, dimension(:) :: cb_value, vb_value
  integer :: whichBand, bandIdx
  complex(kind=dp), allocatable, dimension(:,:,:) :: tensor
  complex(kind=dp) :: aa, bb, cc, dd
  real(kind=dp) :: gfac(2,3)

  ! file handling
  integer(kind=4) :: iounit

  ! --- Wire mode (confinement=2) variables ---
  real(kind=dp), allocatable       :: profile_2d(:,:)
  type(csr_matrix), allocatable    :: kpterms_2d(:)
  type(csr_matrix)                 :: HT_csr
  type(wire_coo_cache)             :: coo_cache
  type(eigensolver_config)         :: eigen_cfg
  type(eigensolver_result)         :: eigen_res
  integer                          :: Ngrid, Ntot, nev_wire

  ! Shared setup: read input, initialize materials, confinement, external field
  call read_and_setup(cfg, profile, kpterms)

  ! g-factor specific validation
  if (cfg%waveVectorStep /= 0 .and. cfg%waveVector /= 'k0') stop 'g-factor calculation requires only k=0'
  if (cfg%waveVectorStep /= 0) then
    print *, 'Warning: Setting wvStep to 0 for g-factor calculation'
    cfg%waveVectorStep = 0
  endif

  ! Initialize k for workspace query
  k = 1

  ! Allocate arrays for single k-point
  if (allocated(smallk)) deallocate(smallk)
  allocate(smallk(1))
  smallk(1)%kx = 0
  smallk(1)%ky = 0
  smallk(1)%kz = 0

  select case(cfg%waveVector)
    case ("k0")
      ! Already set to zero
    case default
      stop "no such direction"
  end select

  whichBand = cfg%whichBand
  bandIdx = cfg%bandIdx

  ! ====================================================================
  ! Branch: wire (confinement=2) vs bulk/QW (confinement=0,1)
  ! ====================================================================
  if (cfg%confinement == 2) then

    ! ----------------------------------------------------------------
    ! Wire mode: sparse Hamiltonian + FEAST eigensolver
    ! ----------------------------------------------------------------

    ! For wire mode, numcb/numvb are used directly from input (not * fdStep)
    cfg%evnum = cfg%numcb + cfg%numvb

    ! Initialize 2D confinement operators (sparse CSR kpterms)
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
        call apply_pikus_bir(strain_out, cfg%params, cfg%grid%material_id, &
          cfg%grid, profile_2d)
        call strain_result_free(strain_out)

        print *, '  Wire strain calculation complete'
      end block
    end if

    Ngrid = grid_ngrid(cfg%grid)
    Ntot  = 8 * Ngrid
    nev_wire = cfg%evnum
    if (nev_wire > Ntot) then
      print *, "Warning: requesting more bands than matrix size. Clamping."
      nev_wire = Ntot
    end if

    print *, ''
    print *, '=== Wire g-factor (2D confinement) ==='
    print *, '  Grid: nx=', cfg%grid%nx, ' ny=', cfg%grid%ny, ' Ngrid=', Ngrid
    print *, '  Matrix size: ', Ntot, 'x', Ntot
    print *, '  numcb=', cfg%numcb, ' numvb=', cfg%numvb
    print *, ''

    ! Build wire Hamiltonian at kz=0
    call ZB8bandGeneralized(HT_csr, 0.0_dp, profile_2d, kpterms_2d, cfg, coo_cache)

    ! Configure FEAST with auto energy window
    eigen_cfg%method   = 'FEAST'
    eigen_cfg%nev      = nev_wire
    eigen_cfg%max_iter = 100
    eigen_cfg%tol      = 1.0e-10_dp
    eigen_cfg%feast_m0 = 0  ! auto: 2*nev

    call auto_compute_energy_window(HT_csr, eigen_cfg%emin, eigen_cfg%emax)

    ! Solve eigenvalue problem
    call solve_sparse_evp(HT_csr, eigen_cfg, eigen_res)

    if (eigen_res%nev_found == 0) then
      print *, 'Error: FEAST found no eigenvalues. Check energy window.'
      print *, '  Window: [', eigen_cfg%emin, ',', eigen_cfg%emax, ']'
      stop 1
    end if

    if (.not. eigen_res%converged) then
      print *, 'Warning: FEAST did not converge for g-factor calculation'
    end if

    if (eigen_res%nev_found < cfg%numcb + cfg%numvb) then
      print *, 'Error: FEAST found', eigen_res%nev_found, &
        ' eigenvalues but need', cfg%numcb + cfg%numvb
      print *, '  numcb=', cfg%numcb, ' numvb=', cfg%numvb
      stop 1
    end if

    if (bandIdx < 1 .or. bandIdx + 1 > cfg%numcb) then
      print *, 'Error: bandIdx=', bandIdx, ' out of range for numcb=', cfg%numcb
      print *, '  Need bandIdx >= 1 and bandIdx+1 <= numcb'
      stop 1
    end if

    N = Ntot

    ! Write wire eigenfunctions at kz=0
    call writeEigenfunctions2d(cfg%grid, eigen_res%eigenvalues, &
      & eigen_res%eigenvectors, 1, eigen_res%nev_found, write_parts=.true.)

    ! Extract CB/VB states from sorted eigenvalues
    ! Eigenvalues are sorted ascending: VB first (lower energy), then CB (higher)
    allocate(cb_value(cfg%numcb))
    allocate(vb_value(cfg%numvb))
    allocate(cb_state(N, cfg%numcb))
    allocate(vb_state(N, cfg%numvb))

    cb_value(:) = eigen_res%eigenvalues(cfg%numvb+1:cfg%numvb+cfg%numcb)
    vb_value(:) = eigen_res%eigenvalues(cfg%numvb:1:-1)

    cb_state(:,:) = eigen_res%eigenvectors(:, cfg%numvb+1:cfg%numvb+cfg%numcb)
    vb_state(:,:) = eigen_res%eigenvectors(:, cfg%numvb:1:-1)

    ! Compute g-factor tensor
    allocate(tensor(2,2,3))
    tensor = 0

    call gfactorCalculation_wire(tensor, whichBand, bandIdx, cfg%numcb, &
      & cfg%numvb, cb_state, vb_state, cb_value, vb_value, cfg, &
      & profile_2d, kpterms_2d)

    ! Compute optical transitions
    block
      type(optical_transition), allocatable :: transitions(:)
      integer :: num_trans, it

      call compute_optical_matrix_wire(transitions, num_trans, &
        cb_state, vb_state, cb_value, vb_value, cfg%numcb, cfg%numvb, &
        profile_2d, kpterms_2d, cfg)

      ! Write to file
      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/optical_transitions.dat', status='replace', action='write')
      write(iounit, '(A)') '# CB VB dE(eV) |px|^2 |py|^2 |pz|^2 f_osc'
      do it = 1, num_trans
        write(iounit, '(2(I4,1x),5(g14.6,1x))') &
          transitions(it)%cb_idx, transitions(it)%vb_idx, &
          transitions(it)%energy, transitions(it)%px, &
          transitions(it)%py, transitions(it)%pz, &
          transitions(it)%oscillator_strength
      end do
      close(iounit)
      print *, '  Optical transitions written to output/optical_transitions.dat'

      deallocate(transitions)
    end block

    ! Free wire-specific resources
    call csr_free(HT_csr)
    call wire_coo_cache_free(coo_cache)
    call eigensolver_result_free(eigen_res)
    if (allocated(profile_2d)) deallocate(profile_2d)
    if (allocated(kpterms_2d)) then
      do i = 1, size(kpterms_2d)
        call csr_free(kpterms_2d(i))
      end do
      deallocate(kpterms_2d)
    end if

  else

    ! ----------------------------------------------------------------
    ! Bulk / QW mode: dense Hamiltonian + LAPACK diagonalization
    ! ----------------------------------------------------------------

    print *, 'adjusting numcb and numvb to get all states'
    cfg%numcb = NUM_CB_STATES*cfg%fdStep
    cfg%numvb = NUM_VB_STATES*cfg%fdStep

    cfg%evnum = cfg%numcb + cfg%numvb

    N = cfg%fdStep*8

    if (cfg%evnum /= N) stop 'evnum not equal to total matrix size'

    NB = ILAENV(1, 'ZHETRD', 'UPLO', N, N, -1, -1)
    NB = MAX(NB,N)
    ABSTOL = DLAMCH('P')

    ! Allocate arrays
    if (allocated(eig)) deallocate(eig)
    allocate(eig(N,1))  ! Only one k-point for g-factor
    eig(:,:) = 0_dp

    allocate(HT(N,N))
    HT = 0.0_dp

    ! Workspace query for LAPACK
    if (allocated(work)) deallocate(work)
    allocate(work(1))
    if (allocated(rwork)) deallocate(rwork)
    allocate(rwork(1))
    if (allocated(iwork)) deallocate(iwork)
    allocate(iwork(1))

    ! Workspace query
    if (cfg%numLayers == 1) then
      call zheev('V', 'U', N, HT, N, eig(:,1), work, -1, rwork, info)
      if (info /= 0) then
        print *, 'Error: zheev workspace query failed, info =', info
        stop 1
      end if
      lwork = int(work(1))
      if (allocated(work)) deallocate(work)
      allocate(work(lwork))
      if (allocated(rwork)) deallocate(rwork)
      allocate(rwork(max(1,3*N-2)))
    else
      call zheevd('V', 'U', N, HT, N, eig(:,1), work, -1, rwork, -1, iwork, -1, info)
      if (info /= 0) then
        print *, 'Error: zheevd workspace query failed, info =', info
        stop 1
      end if
      lwork = int(work(1))
      lrwork = int(rwork(1))
      liwork = iwork(1)
      if (allocated(work)) deallocate(work)
      allocate(work(lwork))
      if (allocated(rwork)) deallocate(rwork)
      allocate(rwork(lrwork))
      if (allocated(iwork)) deallocate(iwork)
      allocate(iwork(liwork))
    end if

    ! Print profile for QW mode
    if (cfg%confDir == 'z') then
      do i = 1, cfg%fdStep, 1
        write(101,*) cfg%z(i), profile(i,1), profile(i,2), profile(i,3)
      end do
    end if

    !only Gamma points
    k = 1

    if (cfg%confDir == 'n') then !BULK
      call ZB8bandBulk(HT,smallk(k),cfg%params(1))
    else if (cfg%confDir == 'z') then ! QUANTUM WELL
      call ZB8bandQW(HT, smallk(k), profile, kpterms)
    else
      stop "verify confinement direction or something else"
    end if


    ! full diagonalization

    if (cfg%numLayers == 1) call zheev('V', 'L', N, HT, N, eig(:,k), work, lwork, rwork, info)

    if (cfg%numLayers > 1) call zheevd('V', 'U', N, HT, N, eig(:,k), work, lwork, rwork, lrwork, iwork, liwork, info)

    if (info /= 0) then
      print *, "Diagonalization error in g-factor calculation, info = ", info
      if (info < 0) print *, "Parameter ", -info, " had illegal value"
      stop 1
    end if


     if (cfg%numLayers > 1 ) call writeEigenfunctions(N, N, HT, k, cfg%fdstep, cfg%z, cfg%numLayers==1)


    allocate(cb_state(N,cfg%numcb))
    allocate(vb_state(N,cfg%numvb))

    allocate(cb_value(cfg%numcb))
    allocate(vb_value(cfg%numvb))

    cb_state(:,:) = HT(:, cfg%numvb+1:cfg%numvb+cfg%numcb)
    vb_state(:,:) = HT(:, cfg%numvb:1:-1)

    cb_value(:) = eig(cfg%numvb+1:cfg%numvb+cfg%numcb,1)
    vb_value(:) = eig(cfg%numvb:1:-1,1)


    allocate(tensor(2,2,3))
    tensor = 0

    if (cfg%numLayers == 1) call gfactorCalculation(tensor, whichBand, bandIdx, cfg%numcb, &
    & cfg%numvb, cb_state, vb_state, cb_value, vb_value, cfg%numLayers, cfg%params, &
    & cfg%startPos(1), cfg%endPos(1), dz=cfg%dz)
    if (cfg%numLayers > 1)  call gfactorCalculation(tensor, whichBand, bandIdx, cfg%numcb, &
    & cfg%numvb, cb_state, vb_state, cb_value, vb_value, cfg%numLayers, cfg%params, &
    & cfg%startPos(1), cfg%endPos(1), profile=profile, kpterms=kpterms, dz=cfg%dz)

  end if  ! confinement /= 2

  print *, 'tensor'
  do i = 1, 2
    write(*,'(2(f9.4,1x,f9.4))') (tensor(i,j,1), j=1,2)
  end do
  print *, ' '
  do i = 1, 2
    write(*,'(2(f9.4,1x,f9.4))') (tensor(i,j,2), j=1,2)
  end do
  print *, ' '
  do i = 1, 2
    write(*,'(2(f9.4,1x,f9.4))') (tensor(i,j,3), j=1,2)
  end do

  if (allocated(rwork)) deallocate(rwork)
  allocate(rwork(3*2-2))
  NB = ILAENV(1, 'ZHETRD', 'UPLO', 2, 2, -1, -1)
  NB = MAX(NB,2)
  lwork = 2*(NB+1)*N
  if (allocated(work)) deallocate(work)
  allocate(work(lwork))

  print *, 'gx'
  aa = tensor(1,1,1)
  bb = tensor(1,2,1)
  cc = tensor(2,1,1)
  dd = tensor(2,2,1)
  print *, (0.5 * ( -sqrt(aa**2 -2*aa*dd+4*bb*cc + dd**2) + aa + dd ) - 0.5 * ( sqrt(aa**2 -2*aa*dd+4*bb*cc + dd**2) + aa + dd ))
  call zheev('N', 'U', 2, tensor(:,:,1), 2, gfac(:,1), work, lwork, rwork, info)
  print *, 2*gfac(1,1) !+ 2.00231

  print *, 'gy'
  aa = tensor(1,1,2)
  bb = tensor(1,2,2)
  cc = tensor(2,1,2)
  dd = tensor(2,2,2)
  print *, (0.5 * ( -sqrt(aa**2 -2*aa*dd+4*bb*cc + dd**2) + aa + dd ) - 0.5 * ( sqrt(aa**2 -2*aa*dd+4*bb*cc + dd**2) + aa + dd ))
  call zheev('N', 'U', 2, tensor(:,:,2), 2, gfac(:,2), work, lwork, rwork, info)
  print *, 2*gfac(1,2) !+ 2.00231

  print *, 'gz'
  aa = tensor(1,1,3)
  bb = tensor(1,2,3)
  cc = tensor(2,1,3)
  dd = tensor(2,2,3)
  print *, (0.5 * ( -sqrt(aa**2 -2*aa*dd+4*bb*cc + dd**2) + aa + dd ) - 0.5 * ( sqrt(aa**2 -2*aa*dd+4*bb*cc + dd**2) + aa + dd ))
  call zheev('N', 'U', 2, tensor(:,:,3), 2, gfac(:,3), work, lwork, rwork, info)
  print *, 2*gfac(1,3) !+ 2.00231

  call get_unit(iounit)
  open(unit=iounit, file='output/gfactor.dat', status="replace", action="write")
  write(iounit,*) 2*gfac(1,1), 2*gfac(1,2), 2*gfac(1,3)
  close(iounit)

  !----------------------------------------------------------------------------


  if (allocated(smallk)) deallocate(smallk)
  if (allocated(HT)) deallocate(HT)
  if (allocated(eig)) deallocate(eig)
  if (allocated(work)) deallocate(work)
  if (allocated(iwork)) deallocate(iwork)
  if (allocated(rwork)) deallocate(rwork)
  if (allocated(cb_state)) deallocate(cb_state)
  if (allocated(vb_state)) deallocate(vb_state)
  if (allocated(cb_value)) deallocate(cb_value)
  if (allocated(vb_value)) deallocate(vb_value)
  if (allocated(tensor)) deallocate(tensor)


end program gfactor
