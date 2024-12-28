program gfactor

  use definitions
  use parameters
  use hamiltonianConstructor
  use finitedifferences
  use OMP_lib
  use outputFunctions
  use gfactorFunctions
  !use mkl_spblas
  use utils


  implicit NONE

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
  integer :: i, k, j, ii ,jj, lin, col

  ! hamiltonian and LAPACK/BLAS
  integer :: info, ILAENV, NB,  lwork, N, M, il, iuu, vl, vu, lrwork, liwork
  real(kind=dp) :: abstol, dlamch
  real(kind=dp), allocatable :: eig(:,:), rwork(:)
  complex(kind=dp), allocatable :: work(:)
  complex(kind=dp), allocatable, dimension(:,:,:) :: eigv
  integer, allocatable :: iwork(:), ifail(:)
  complex(kind=dp), allocatable, dimension(:,:) :: HT
  real(kind = dp), allocatable, dimension(:,:,:) :: kpterms
  real(kind = dp), allocatable, dimension(:,:) :: profile

  ! file handling
  integer ( kind = 4 ) :: data_unit, iounit
  integer :: status
  character ( len = 255 ) :: data_filename, label

  ! gfactor
  complex(kind=dp), allocatable, dimension(:,:) :: cb_state, vb_state
  real(kind=dp), allocatable, dimension(:) :: cb_value, vb_value
  integer :: whichBand, bandIdx, whichG
  complex(kind=dp), allocatable, dimension(:,:,:) :: tensor
  complex(kind=dp) :: aa, bb, cc, dd
  real(kind=dp) :: gfac(2,3)

  complex(kind=dp) :: res
  real(kind=dp) :: lim

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
    print *, trim(label), trim(material(1))
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

  if (wvStep /= 0 .and. wvDir /= 'k0') stop 'g-factor calculation requires only k=0'
  if (wvStep /= 0) then
    print *, 'Warning: Setting wvStep to 0 for g-factor calculation'
    wvStep = 0
  endif

  ! Initialize k for workspace query
  k = 1

  ! Allocate arrays for single k-point
  if (allocated(smallk)) deallocate(smallk)
  allocate(smallk(1))
  smallk(1)%kx = 0
  smallk(1)%ky = 0
  smallk(1)%kz = 0

  select case(wvDir)
    case ("k0")
      ! Already set to zero
    case default
      stop "no such direction"
  end select

  read(data_unit, *) label, numcb
  print *, trim(label), numcb
  read(data_unit, *) label, numvb
  print *, trim(label), numvb

  print *, 'adjusting numcb and numvb to get all states'
  numcb = 2!2*fdStep
  numvb = 6*fdStep

  evnum = numcb+numvb

  if (evnum /= N) stop 'evnum not equal to total matrix size'

  NB = ILAENV(1, 'ZHETRD', 'UPLO', N, N, -1, -1)
  NB = MAX(NB,N)
  ABSTOL = DLAMCH('P')

  ! For g-factor calculation, we only need one k-point
  if (wvDir /= 'k0') stop 'g-factor calculation requires only k=0'

  ! Initialize k for workspace query
  k = 1

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
  if (nlayers == 1) then
    call zheev('V', 'U', N, HT, N, eig(:,1), work, -1, rwork, info)
    lwork = int(work(1))
    if (allocated(work)) deallocate(work)
    allocate(work(lwork))
    if (allocated(rwork)) deallocate(rwork)
    allocate(rwork(max(1,3*N-2)))
  else
    call zheevd('V', 'U', N, HT, N, eig(:,1), work, -1, rwork, -1, iwork, -1, info)
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

  ! Actual diagonalization
  if (nlayers == 1) then
    call zheev('V', 'U', N, HT, N, eig(:,1), work, lwork, rwork, info)
  else
    call zheevd('V', 'U', N, HT, N, eig(:,1), work, lwork, rwork, lrwork, iwork, liwork, info)
  end if

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

  !only Gamma points
  k = 1

  if (confDir == 'n') then !BULK
    call ZB8bandBulk(HT,smallk(k),params(1))
  else if (confDir == 'z') then ! QUANTUM WELL
    call ZB8bandQW(HT, smallk(k), profile, kpterms)
  else
    stop "verify confinement direction or something else"
  end if

  ! do lin=1,N
  !  do col=lin,N
  !     if ( zabs(HT(lin,col) - dconjg(HT(col,lin))) .gt. 1e-4) then
  !        print *, lin, col, HT(lin,col), HT(col,lin), zabs(HT(lin,col)-dconjg(HT(col,lin)))
  !     end if
  !  end do
  ! end do

  ! LIM = 100.0 * MAXVAL(ABS(SUM(HT(:,:),2))) * EPSILON(ABS(HT(1,1)))
  ! DO I=1,SIZE(HT,2)
  !    WHERE(ABS(REAL(HT(:,I)))  < LIM) HT(:,I)=dCMPLX(0.0,AIMAG(HT(:,I)))
  !    WHERE(ABS(AIMAG(HT(:,I))) < LIM) HT(:,I)=dCMPLX(REAL(HT(:,I)),0.0)
  ! ENDDO

  ! full diagonalization

  if (nlayers == 1) call zheev('V', 'L', N, HT, N, eig(:,k), work, lwork, rwork, info)
  ! if (info /= 0) stop "error diag"
  ! print *, "k = ", smallk(k)%kx, eig(:,k)

  if (nlayers > 1) call zheevd('V', 'U', N, HT, N, eig(:,k), work, lwork, rwork, lrwork, iwork, liwork, info)

  if (info /= 0) stop "error diag"
  ! stop
  ! print *, eig(:,1)

  ! print *, HT

  ! if (smallk(k)%kx==0)
   if (nlayers > 1 ) call writeEigenfunctions(N, N, HT, k, fdstep, z, nlayers==1)

  ! call writeEigenfunctions(N, evnum, HT(1:N,il:iuu), k, fdstep)

  ! call writeEigenvalues(smallk, eig(1:evnum,:), wvStep)
  ! call writeEigenvalues(smallk, eig(il:iuu ,:), wvStep)

  whichBand = 0 !cb
  bandIdx = 1 !first band
  numcb = 2!2*fdStep
  numvb = 6*fdStep

  allocate(cb_state(N,numcb))
  allocate(vb_state(N,numvb))

  allocate(cb_value(numcb))
  allocate(vb_value(numvb))

  do i = 1, N
    cb_state(i,1:numcb) = HT(i,numvb+1:numvb+numcb)
    vb_state(i,1:numvb) = HT(i,numvb:1:-1)
  end do

  cb_value(:) = eig(numvb+1:numvb+numcb,1)
  vb_value(:) = eig(numvb:1:-1,1)


  allocate(tensor(2,2,3))
  tensor = 0

  if (nlayers == 1) call gfactorCalculation(tensor, whichBand, bandIdx, numcb, &
  & numvb, cb_state, vb_state, cb_value, vb_value, nlayers, params, &
  & startPos(1), endPos(1))
  if (nlayers > 1)  call gfactorCalculation(tensor, whichBand, bandIdx, numcb, &
  & numvb, cb_state, vb_state, cb_value, vb_value, nlayers, params, &
  & startPos(1), endPos(1), profile=profile, kpterms=kpterms, dz=(z(2)-z(1)))

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
  if (allocated(cb_state)) deallocate(cb_state)
  if (allocated(vb_state)) deallocate(vb_state)
  if (allocated(cb_value)) deallocate(cb_value)
  if (allocated(vb_value)) deallocate(vb_value)


end program gfactor