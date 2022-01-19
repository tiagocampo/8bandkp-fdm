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
  integer, allocatable :: iwork(:), ifail(:)
  complex(kind=dp), allocatable, dimension(:,:) :: HT
  real(kind = dp), allocatable, dimension(:,:,:) :: kpterms
  real(kind = dp), allocatable, dimension(:,:) :: profile

  ! file handling
  integer ( kind = 4 ) :: data_unit
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
  evnum = numcb+numvb
  vl = 0
  vu = 0
  il = 6*fdStep - numvb + 1
  iuu = 6*fdStep + numcb

  NB = ILAENV(1, 'ZHETRD', 'UPLO', N, N, -1, -1)
  NB = MAX(NB,N)
  ABSTOL = DLAMCH('P')

  if (allocated(rwork)) deallocate(rwork)
  allocate(rwork(7*N))
  if (allocated(iwork)) deallocate(iwork)
  allocate(iwork(5*N))
  if (allocated(ifail)) deallocate(ifail)
  allocate(ifail(N))
  if (allocated(eig)) deallocate(eig)
  allocate(eig(N,wvStep))
  if (allocated(eigv)) deallocate(eigv)
  allocate(eigv(N,evnum,wvStep))

  eig(:,:) = 0_dp
  eigv(:,:,:) = 0_dp

  ! lwork = 2*(NB+1)*N
  ! if (allocated(work)) deallocate(work)
  ! allocate(work(lwork))

  allocate(HT(N,N))
  HT = 0.0_dp


  if (allocated(work)) deallocate(work)
  allocate(work(1))
  lwork = -1
  call zheevx('V','I','U',N,HT,N,vl,vu,il,iuu,ABSTOL,M,eig(:,k),eigv(:,:,k),&
  & N,work,lwork,rwork,iwork,ifail,info)
  lwork = work(1)
  if (allocated(work)) deallocate(work)
  allocate(work(lwork))

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
    do i = 1, fdStep, 1
      write(101,*) z(i), profile(i,1), profile(i,2), profile(i,3)
    end do
  end if


  do k = 1, wvStep, 1

    print *, smallk(k)%kx, smallk(k)%ky, smallk(k)%kz
    if (confDir == 'n') then !BULK

      call ZB8bandBulk(HT,smallk(k),params(1))

    else if (confDir == 'z') then ! QUANTUM WELL

      call ZB8bandQW(HT, smallk(k), profile, kpterms)

      ! do lin=1,N
      !  do col=lin,N
      !     if ( zabs(HT(lin,col) - dconjg(HT(col,lin))) .gt. 1e-4) then
      !        print *, lin, col, HT(lin,col), HT(col,lin), zabs(HT(lin,col)-dconjg(HT(col,lin)))
      !     end if
      !  end do
      ! end do
      !
      ! do ii = 1, N, 1
      !   write(103,*) (real(HT(ii,jj)), jj=1, N)
      ! end do
      ! stop

    else

      stop "verify confinement direction or something else"

    end if

    ! full diagonalization
    ! call zheev('V', 'U', N, HT, N, eig(:,k), work, lwork, rwork, info)
    ! if (info /= 0) stop "error diag"
    !print *, "k = ", smallk(k)%kx, eig(:,k)
    ! write(100,*) smallk(k)%kx, eig(il:iu,k)

    ! print *, 'before diag'
    call zheevx('V', 'I', 'U', N, HT, N, vl, vu, il, iuu, ABSTOL, M, &
    & eig(:,k), eigv(1:N,:,k), N, work, lwork, rwork, iwork, ifail, info)
    if (info /= 0) stop "error diag"
    ! print *, 'after diag'
    ! write(100,*) smallk(k)%kx, eig(1:evnum,k)

    if (smallk(k)%kx==0) call writeEigenfunctions(N, evnum, &
    & eigv(1:N,1:evnum,k), k, fdstep, z)
    ! call writeEigenfunctions(N, evnum, HT(1:N,il:iu), k, fdstep)

  end do
  call writeEigenvalues(smallk, eig(1:evnum,:), wvStep)
  ! call writeEigenvalues(smallk, eig(il:iu,:), wvStep)


  !----------------------------------------------------------------------------


  if (allocated(smallk)) deallocate(smallk)
  if (allocated(HT)) deallocate(HT)
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
