program kpfdm

  use definitions
  use parameters
  use hamiltonianConstructor
  use finitedifferences
  use OMP_lib
  use outputFunctions

  implicit none

  ! configuration / initialization
  integer :: confinement, nlayers
  character (len = 255), allocatable, dimension(:) :: material
  type(paramStruct), allocatable, dimension(:) :: params
  real(kind=dp), allocatable, dimension(:) :: z, startPos, endPos, bshift
  integer, allocatable, dimension(:) :: intStartPos, intEndPos
  integer :: fdStep, evnum, numcb, numvb
  character (len = 1) :: confDir
  real(kind=dp) :: totalSize, delta
  integer :: externalField
  character(len = 2) :: EFtype
  real(kind=dp) :: Evalue

  ! wave vector
  character (len = 2) :: wvDir
  real(kind=sp) :: wvMax ![1/AA]
  integer :: wvStep
  type(wavevector), allocatable, dimension(:) :: smallk

  ! iteration consts
  integer :: i, k, ii ,jj, lin, col

  ! hamiltonian and LAPACK/BLAS
  integer :: info, ILAENV, NB,  lwork, N, M, il, iuu, vl, vu
  real(kind=dp) :: abstol, dlamch
  real(kind=dp), allocatable :: eig(:,:), rwork(:)
  complex(kind=dp), allocatable :: work(:)
  complex(kind=dp), allocatable, dimension(:,:,:) :: eigv
  complex(kind=dp), allocatable, dimension(:,:) :: HT, HTmp
  integer, allocatable :: iwork(:), ifail(:)
  real(kind = dp), allocatable, dimension(:,:,:) :: kpterms
  real(kind = dp), allocatable, dimension(:,:) :: profile

  ! file handling
  integer ( kind = 4 ) :: data_unit, iounit
  integer :: status
  character ( len = 255 ) :: data_filename, label

  !reading input/configuration file
  call get_unit ( data_unit )
  data_filename = 'input.cfg'
  open( data_unit ,file=data_filename ,status="old")


  read(data_unit, *) label, wvDir
  print *, trim(label), wvDir
  read(data_unit, *) label, wvMax
  print *, trim(label), wvMax
  read(data_unit, *) label, wvStep
  print *, trim(label), wvStep
  read(data_unit, *) label, confinement
  print *, trim(label), confinement
  read(data_unit, *) label, fdStep
  print *, trim(label), fdStep
  read(data_unit, *) label, nlayers
  print *, trim(label), nlayers

  N = fdStep*8

  allocate(params(nlayers))
  allocate(material(nlayers))

  allocate(z(fdStep))

  if (confinement == 0 .and. nlayers == 1)  then
    read(data_unit, *) label, material(1)
    print *, trim(label), material(1)
    material(1) = trim(material(1))
    confDir = 'n'
  end if

  if (confinement == 1) then

    allocate(startPos(nlayers))
    allocate(endPos(nlayers))
    allocate(intStartPos(nlayers))
    allocate(intEndPos(nlayers))
    allocate(bshift(nlayers))

    do i = 1, nlayers, 1
      read(data_unit, *) label, material(i), startPos(i), endPos(i)
      print *, trim(label), trim(material(i)), startPos(i), endPos(i)
    end do

    totalSize = endPos(1)-startPos(1)
    delta = totalSize/(fdStep-1)

    z = [ (startPos(1) + (i-1)*delta, i=1,fdStep) ]

    confDir = 'z'


    do i = 1, nlayers, 1
      intStartPos(i) = int(abs( (startPos(1)-startPos(i))/(totalSize/(fdStep-1))) )
      intEndPos(i) = int (intStartPos(1) + ((endPos(1)+endPos(i))/delta) )
    end do
    intStartPos = intStartPos + 1
    intEndPos = intEndPos + 1

  end if

  allocate(smallk(wvStep))
  smallk%kx = 0
  smallk%ky = 0
  smallk%kz = 0
  select case(wvDir)

    case ("kx")
      smallk%kx = [ ((i-1)*wvMax/(wvStep-1), i=1,wvStep) ]![ (i*wvMax/(wvStep-1), i=-int(wvStep/2),int(wvStep/2)) ]

    case ("ky")
      smallk%ky = [ ((i-1)*wvMax/(wvStep-1), i=1,wvStep) ]

    case ("kz")
      smallk%kz = [ ((i-1)*wvMax/(wvStep-1), i=1,wvStep) ]

    case ("k0")
      smallk%kx = 0
      smallk%ky = 0
      smallk%kz = 0

    case default
      stop "no such direction"

  end select

  read(data_unit, *) label, numcb
  print *, trim(label), numcb
  read(data_unit, *) label, numvb
  print *, trim(label), numvb
  evnum = numcb + numvb
  
  ! Set matrix dimensions and eigenvalue range
  N = fdStep * 8
  if (confDir == 'n') then
    N = 8  ! For bulk, we only need 8x8
    if (evnum > 8) then
      print *, "Warning: requesting more bands than available in bulk (8). Using all 8 bands."
      evnum = 8
      numcb = 2
      numvb = 6
    end if
  end if

  ! For quantum well, check if requested bands are within limits
  if (confDir == 'z') then
    if (numcb > 2*fdStep) then
      print *, "Warning: requesting more conduction bands than available. Limiting to ", 2*fdStep
      numcb = 2*fdStep
    end if
    if (numvb > 6*fdStep) then
      print *, "Warning: requesting more valence bands than available. Limiting to ", 6*fdStep
      numvb = 6*fdStep
    end if
    evnum = numcb + numvb
  end if

  vl = 0.0_dp
  vu = 0.0_dp
  if (confinement == 0) then
    il = 1
    iuu = evnum
  else
    ! For quantum well, select the right range of states
    ! We want the highest numvb valence states and lowest numcb conduction states
    il = 6*fdStep - numvb + 1  ! Start from highest valence band
    iuu = 6*fdStep + numcb     ! Up to highest conduction band
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
  allocate(eig(iuu-il+1,wvStep))  ! Only store the states we want
  if (allocated(eigv)) deallocate(eigv)
  if (confDir == 'n') then
    allocate(eigv(8,evnum,wvStep))  ! 8x8 for bulk
  else
    allocate(eigv(N,iuu-il+1,wvStep))  ! Only store the states we want
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

  read(data_unit, *) label, externalField, EFtype
  print *, trim(label), externalField, EFtype
  if (EFtype == "EF") then
    read(data_unit, *) label, Evalue
    print *, trim(label), Evalue
  else
    stop "Type of external field not implemented"
  end if

  print *, ' '

  !----------------------------------------------------------------------------
  call paramDatabase(material,nlayers,params)

  if (confDir == 'z') then
    allocate(kpterms(fdStep,fdStep,10))
    kpterms = 0.0_dp
    call confinementInitialization(z, intStartPos, intEndPos, material, &
    & nlayers, params, bshift, confDir, profile, kpterms)
    if (externalField == 1 .and. EFType == "EF") then
      call externalFiledSetup_electricField(profile, Evalue, totalSize, z)
    end if
    !printing profile
    call get_unit(iounit)
    open(unit=iounit, file='output/potential_profile.dat', status="replace", action="write")
    do i = 1, fdStep, 1
      write(iounit,*) z(i), profile(i,1), profile(i,2), profile(i,3)
    end do
    close(iounit)
  end if


  do k = 1, wvStep, 1

    print *, smallk(k)%kx, smallk(k)%ky, smallk(k)%kz
    if (confDir == 'n') then !BULK
      call ZB8bandBulk(HT,smallk(k),params(1))
      ! Copy to temporary array
      HTmp = HT(1:8,1:8)

      ! Query optimal workspace
      if (k == 1) then
        if (allocated(work)) deallocate(work)
        allocate(work(1))
        lwork = -1
        call zheevx('V', 'I', 'U', 8, HTmp, 8, vl, vu, il, iuu, abstol, M, eig(1:evnum,k), &
                   HTmp, 8, work, lwork, rwork, iwork, ifail, info)
        lwork = int(real(work(1)))
        if (allocated(work)) deallocate(work)
        allocate(work(lwork))
      end if

      ! Actual diagonalization - use zheevx for selected eigenvalues
      call zheevx('V', 'I', 'U', 8, HTmp, 8, vl, vu, il, iuu, abstol, M, eig(1:evnum,k), &
                 HTmp, 8, work, lwork, rwork, iwork, ifail, info)
      if (info /= 0) then
        print *, "Diagonalization error in bulk calculation, info = ", info
        if (info < 0) then
          print *, "Parameter ", -info, " had illegal value"
        end if
        stop "error diag"
      end if

      ! For bulk, just copy the selected eigenvectors
      eigv(:,:,k) = HTmp(:,1:evnum)

    else if (confDir == 'z') then ! QUANTUM WELL
      call ZB8bandQW(HT, smallk(k), profile, kpterms)

      ! Query optimal workspace for full matrix
      if (k == 1) then
        if (allocated(work)) deallocate(work)
        allocate(work(1))
        lwork = -1
        call zheevx('V', 'I', 'U', N, HT, N, vl, vu, il, iuu, abstol, M, eig(:,k), &
                   HT, N, work, lwork, rwork, iwork, ifail, info)
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
    if (k == 1 .or. k == int(wvStep/2) .or. k == wvStep) then  ! Write at start, middle, and end k-points
      if (confDir == 'n') then
        call writeEigenfunctions(8, min(evnum,8), HTmp(:,1:min(evnum,8)), k, fdstep, z, .true.)
      else
        ! For quantum well, write only the requested number of states
        call writeEigenfunctions(N, iuu-il+1, HT(:,1:iuu-il+1), k, fdstep, z, .false.)
      end if
    end if

  end do
  call writeEigenvalues(smallk, eig(:,:), wvStep)
  ! call writeEigenvalues(smallk, eig(il:iu,:), wvStep)


  !----------------------------------------------------------------------------


  if (allocated(smallk)) deallocate(smallk)
  if (allocated(HT)) deallocate(HT)
  if (allocated(HTmp)) deallocate(HTmp)
  if (allocated(work)) deallocate(work)
  if (allocated(iwork)) deallocate(iwork)
  if (allocated(rwork)) deallocate(rwork)
  if (allocated(ifail)) deallocate(ifail)
  if (allocated(params)) deallocate(params)
  if (allocated(material)) deallocate(material)
  if (allocated(z)) deallocate(z)
  if (allocated(bshift)) deallocate(bshift)
  if (allocated(profile)) deallocate(profile)
  if (allocated(kpterms)) deallocate(kpterms)


end program kpfdm