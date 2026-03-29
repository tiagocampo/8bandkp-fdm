program kpfdm

  use definitions
  use parameters
  use hamiltonianConstructor
  use finitedifferences
  use OMP_lib
  use outputFunctions
  use input_parser

  implicit none

  ! Shared configuration from input_parser
  type(simulation_config) :: cfg
  real(kind=dp), allocatable, dimension(:,:) :: profile
  real(kind=dp), allocatable, dimension(:,:,:) :: kpterms

  ! wave vector
  type(wavevector), allocatable, dimension(:) :: smallk

  ! iteration consts
  integer :: i, k, ii, jj, lin, col

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
    if (cfg%numcb > 2*cfg%fdStep) then
      print *, "Warning: requesting more conduction bands than available. Limiting to ", 2*cfg%fdStep
      cfg%numcb = 2*cfg%fdStep
    end if
    if (cfg%numvb > 6*cfg%fdStep) then
      print *, "Warning: requesting more valence bands than available. Limiting to ", 6*cfg%fdStep
      cfg%numvb = 6*cfg%fdStep
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
    il = 6*cfg%fdStep - cfg%numvb + 1  ! Start from highest valence band
    iuu = 6*cfg%fdStep + cfg%numcb     ! Up to highest conduction band
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
  if (allocated(work)) deallocate(work)
  if (allocated(iwork)) deallocate(iwork)
  if (allocated(rwork)) deallocate(rwork)
  if (allocated(ifail)) deallocate(ifail)


end program kpfdm
