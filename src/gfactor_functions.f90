module gfactorFunctions

  use definitions
  use hamiltonianConstructor
  use mkl_spblas
  use OMP_lib

  implicit none

  contains

complex(kind=dp) function sigmaElem(state1, state2, dir, fdstep, startz, endz, dz, nlayers)

  complex(kind=dp), intent(in), dimension(:) :: state1, state2
  character(len=1), intent(in) :: dir
  real(kind=dp), intent(in) :: startz, endz, dz
  integer, intent(in) :: fdstep, nlayers

  complex(kind=dp) :: aux(8,8), aux_v(8), aux1_v(8), Y(8)
  complex(kind=dp), allocatable, dimension(:) :: v
  integer :: i, j

  complex(kind=dp) :: SIGMA_X(8,8), SIGMA_Y(8,8), SIGMA_Z(8,8), alpha, beta
  complex(kind=dp) :: auxx(8,8), auxy(8,8), auxz(8,8), zdotc

  integer :: lin, col

  SIGMA_X(1,:) = (/ ZERO, -IU*RQS3, ZERO, ZERO, UM*SQR2o3, ZERO, ZERO, ZERO /)
  SIGMA_X(2,:) = (/ IU*RQS3, ZERO, -IU*2.0_dp/3.0_dp, ZERO, ZERO, -UM*SQR2/3.0_dp, ZERO, ZERO /)
  SIGMA_X(3,:) = (/ ZERO, 2.0_dp*IU/3.0_dp, ZERO, -IU*RQS3, UM*SQR2/3.0_dp, ZERO, ZERO, ZERO /)
  SIGMA_X(4,:) = (/ ZERO, ZERO, IU*RQS3, ZERO, ZERO, -UM*SQR2o3, ZERO, ZERO /)
  SIGMA_X(5,:) = (/ UM*SQR2o3, ZERO, UM*SQR2/3.0_dp, ZERO, ZERO, -IU/3.0_dp, ZERO, ZERO /)
  SIGMA_X(6,:) = (/ ZERO, -UM*SQR2/3.0_dp, ZERO, -UM*SQR2o3, IU/3.0_dp, ZERO, ZERO, ZERO /)
  SIGMA_X(7,:) = (/ ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, UM /)
  SIGMA_X(8,:) = (/ ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, UM, ZERO /)

  ! do lin=1,8
  !  do col=lin,8
  !     if ( zabs(sigma_x(lin,col) - dconjg(sigma_x(col,lin))) .gt. 1e-4) then
  !        print *, lin, col, sigma_x(lin,col), sigma_x(col,lin), zabs(sigma_x(lin,col)-dconjg(sigma_x(col,lin)))
  !     end if
  !  end do
  ! end do

  SIGMA_Y(1,:) = (/ ZERO, UM*RQS3, ZERO, ZERO, IU*SQR2o3, ZERO, ZERO, ZERO /)
  SIGMA_Y(2,:) = (/ UM*RQS3, ZERO, UM*2.0_dp/3.0_dp, ZERO, ZERO, -IU*SQR2/3.0_dp, ZERO, ZERO /)
  SIGMA_Y(3,:) = (/ ZERO, UM*2.0_dp/3.0_dp, ZERO, UM*RQS3, -IU*SQR2/3.0_dp, ZERO, ZERO, ZERO /)
  SIGMA_Y(4,:) = (/ ZERO, ZERO, UM*RQS3, ZERO, ZERO, IU*SQR2o3, ZERO, ZERO /)
  SIGMA_Y(5,:) = (/ -IU*SQR2o3, ZERO, IU*SQR2/3.0_dp, ZERO, ZERO, UM/3.0_dp, ZERO, ZERO /)
  SIGMA_Y(6,:) = (/ ZERO, IU*SQR2/3.0_dp, ZERO, -IU*SQR2o3, UM/3.0_dp, ZERO, ZERO, ZERO /)
  SIGMA_Y(7,:) = (/ ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, IU /)
  SIGMA_Y(8,:) = (/ ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, -IU, ZERO /)

  ! do lin=1,8
  !  do col=lin,8
  !     if ( zabs(sigma_y(lin,col) - dconjg(sigma_y(col,lin))) .gt. 1e-4) then
  !        print *, lin, col, sigma_y(lin,col), sigma_y(col,lin), zabs(sigma_y(lin,col)-dconjg(sigma_y(col,lin)))
  !     end if
  !  end do
  ! end do

  SIGMA_Z(1,:) = (/ UM, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO /)
  SIGMA_Z(2,:) = (/ ZERO, UM/3.0_dp, ZERO, ZERO, -2.0_dp*IU*SQR2/3.0_dp, ZERO, ZERO, ZERO /)
  SIGMA_Z(3,:) = (/ ZERO, ZERO, -UM/3.0_dp, ZERO, ZERO, 2.0_dp*IU*SQR2/3.0_dp, ZERO, ZERO /)
  SIGMA_Z(4,:) = (/ ZERO, ZERO, ZERO, -UM, ZERO, ZERO, ZERO, ZERO /)
  SIGMA_Z(5,:) = (/ ZERO, 2.0_dp*IU*SQR2/3.0_dp, ZERO, ZERO, -UM/3.0_dp, ZERO, ZERO, ZERO /)
  SIGMA_Z(6,:) = (/ ZERO, ZERO, -2.0_dp*IU*SQR2/3.0_dp, ZERO, ZERO, UM/3.0_dp, ZERO, ZERO /)
  SIGMA_Z(7,:) = (/ ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, UM, ZERO /)
  SIGMA_Z(8,:) = (/ ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, -UM /)

  ! do lin=1,8
  !  do col=lin,8
  !     if ( zabs(sigma_z(lin,col) - dconjg(sigma_z(col,lin))) .gt. 1e-4) then
  !        print *, lin, col, sigma_z(lin,col), sigma_z(col,lin), zabs(sigma_z(lin,col)-dconjg(sigma_z(col,lin)))
  !     end if
  !  end do
  ! end do

  ! do i = 1, 8
  !   write(200,'(8(f9.4,1x,f9.4))') (SIGMA_X(i,j), j=1,8)
  ! end do
  ! print *, ' '
  ! do i = 1, 8
  !   write(201,'(8(f9.4,1x,f9.4))') (SIGMA_Y(i,j), j=1,8)
  ! end do
  ! print *, ' '
  ! do i = 1, 8
  !   write(202,'(8(f9.4,1x,f9.4))') (SIGMA_Z(i,j), j=1,8)
  ! end do
  ! stop

  allocate(v(fdstep))
  v = 0
  alpha = dcmplx(1.,0.)
  beta = dcmplx(0.,0.)
  sigmaElem = 0
  aux_v = 0
  aux1_v = 0


  if (nlayers > 1) then

    Y = 0

    if ( dir == 'x' ) then
      do i = 1, fdstep
        aux_v = [( state2(i + (j-1)*fdstep), j=1,8 )]
        aux1_v = [( state1(i + (j-1)*fdstep), j=1,8 )]
        call zgemv('N',8,8,alpha,SIGMA_X,8,aux_v,1,beta,Y,1)
        v(i) = zdotc(8,aux1_v,1,Y(:),1)
      end do
    elseif ( dir == 'y' ) then
      do i = 1, fdstep
        aux_v = [( state2(i + (j-1)*fdstep), j=1,8 )]
        aux1_v = [( state1(i + (j-1)*fdstep), j=1,8 )]
        call zgemv('N',8,8,alpha,SIGMA_Y,8,aux_v,1,beta,Y,1)
        v(i) = zdotc(8,aux1_v,1,Y(:),1)
      end do
    elseif ( dir == 'z' ) then
      do i = 1, fdstep
        aux_v = [( state2(i + (j-1)*fdstep), j=1,8 )]
        aux1_v = [( state1(i + (j-1)*fdstep), j=1,8 )]
        call zgemv('N',8,8,alpha,SIGMA_Z,8,aux_v,1,beta,Y,1)
        v(i) = zdotc(8,aux1_v,1,Y(:),1)
      end do
    end if


    sigmaElem = simpson(v,startz,endz)/dz
  else
    Y = 0
    ! print *, state2
    ! print *, state1
    if ( dir == 'x' ) then
      call zgemv('N',8,8,alpha,SIGMA_X,8,state2,1,beta,Y,1)
    elseif ( dir == 'y' ) then
      call zgemv('N',8,8,alpha,SIGMA_Y,8,state2,1,beta,Y,1)
    elseif ( dir == 'z' ) then
      call zgemv('N',8,8,alpha,SIGMA_Z,8,state2,1,beta,Y,1)
    end if

    sigmaElem = zdotc(8,state1,1,Y(:),1)

  end if

  deallocate(v)

end function


subroutine pMatrixEleCalc(Pele,d,state1,state2,nlayers,params, Ppartial, &
  & profile, kpterms, HT_csr, hkp_out, startz, endz, dz)

  complex(kind=dp), intent(inout) :: Pele
  integer, intent(in) :: d
  complex(kind=dp), intent(in), dimension(:) :: state1, state2
  integer, intent(in) :: nlayers
  type(paramStruct), intent(in) :: params(nlayers)
  complex(kind=dp), intent(inout), optional, dimension(:) :: Ppartial
  real(kind = dp), intent(in), optional, dimension(:,:) :: profile
  real(kind = dp), intent(in), optional, dimension(:,:,:) :: kpterms
  type(sparse_matrix_T), intent(in), optional :: HT_csr
  real(kind=dp), intent(in), optional :: startz, endz, dz
  complex(kind=dp), intent(in), optional, dimension(:,:) :: hkp_out

  type(wavevector), allocatable, dimension(:) :: smallk
  integer :: dimax
  complex(kind=dp), allocatable, dimension(:,:) :: hkp

  type(matrix_descr) :: descrA
  integer :: pm(128), info

  complex(kind=dp) :: alpha, beta
  complex(kind=dp) :: zdotc

  integer :: i, j, m, n, fdstep
  complex(kind=dp) :: aux(8,8), aux_v(8), aux1_v(8), Y(8)
  complex(kind=dp), allocatable, dimension(:) :: v, Y1

  dimax = size(state1,1)
  fdstep = dimax/8

  allocate(smallk(1))
  if (.not. present(HT_csr)) then
    allocate(hkp(dimax,dimax))
    hkp(:,:) = cmplx(0,0)
  end if

  if (present(hkp_out)) hkp = hkp_out

  alpha = cmplx(1.,0.)
  beta = cmplx(0.,0.)

  if ( d==1 ) then
    smallk%kx = 1
    smallk%ky = 0
    smallk%kz = 0
  end if
  if ( d==2 ) then
    smallk%kx = 0
    smallk%ky = 1
    smallk%kz = 0
  end if
  if ( d==3 ) then
    smallk%kx = 0
    smallk%ky = 0
    smallk%kz = 1
  end if
  if (nlayers == 1) then
    call ZB8bandBulk(hkp,smallk(1),params(1),g='g')
  end if
  if (nlayers == 2 .and. .not. present(HT_csr) .and. .not. present(hkp_out)) then
    ! print *, 'here'
    call ZB8bandQW(hkp, smallk(1), profile, kpterms,g='g')
  end if


  if(present(Ppartial)) then
     call zgemv('N',dimax,dimax,alpha,hkp(:,:),dimax,state2(:),1,beta,Ppartial,1)
  elseif ( present(HT_csr) ) then


    allocate(Y1(dimax))
    !create matrix descriptor
    descrA%TYPE =SPARSE_MATRIX_TYPE_GENERAL

    info = mkl_sparse_z_mv (SPARSE_OPERATION_NON_TRANSPOSE, alpha, HT_csr, descrA, state2(:), beta, Y1)

    Pele = zdotc(dimax,state1(:),1,Y1(:),1)

    deallocate(Y1)

  else
    ! print *, 'here'
    allocate(v(fdstep))
    v = 0
    aux_v = 0
    aux1_v = 0

    !$OMP PARALLEL DO private(i, n, m, aux, aux_v, aux1_v, Y) shared(v, hkp, fdstep, state1, state2, alpha, beta)
    do i = 1, fdstep

      do n = 1, 8, 1
        do m = 1, 8, 1
          ! print *, i + (m-1)*fdstep, i + (n-1)*fdstep
          aux(m,n) = hkp(i + (m-1)*fdstep,i + (n-1)*fdstep)
        end do
      end do

      aux_v = [( state2(i + (j-1)*fdstep), j=1,8 )]
      aux1_v = [( state1(i + (j-1)*fdstep), j=1,8 )]

      call zgemv('N',8,8,alpha,aux,8,aux_v,1,beta,Y,1)
      v(i) = zdotc(8,aux1_v,1,Y(:),1)
    end do
    !$OMP END PARALLEL DO
    if (nlayers == 1) Pele = v(1)
    if (nlayers > 1)  Pele = simpson(v,startz,endz)/dz

  end if
  ! stop

  if (allocated(hkp)) deallocate(hkp)

end subroutine pMatrixEleCalc

subroutine gfactorCalculation(tensor, whichBand, bandIdx, numcb, numvb, &
  & cb_state, vb_state, cb_value, vb_value, nlayers, params, startz, endz, &
  & profile, kpterms, dz)

  implicit none


  complex(kind=dp), intent(inout), dimension(:,:,:) :: tensor
  integer, intent(in) :: whichBand !=0 cb, =1 vb
  integer, intent(in) :: bandIdx !idx of band to compute. 1, 2, 3 ... (mostly for confined systems)
  integer, intent(in) :: numcb, numvb ! total number of cb and vb bands
  complex(kind=dp), intent(in), dimension(:,:) :: cb_state, vb_state !env func at k=0
  real(kind=dp), intent(in), dimension(:) :: cb_value, vb_value ! energies at k=0
  integer, intent(in) :: nlayers
  type(paramStruct), intent(in) :: params(nlayers)
  real(kind=dp), intent(in) :: startz, endz

  real(kind = dp), intent(in), optional, dimension(:,:) :: profile
  real(kind = dp), intent(in), optional, dimension(:,:,:) :: kpterms
  real(kind=dp), intent(in), optional :: dz

  integer :: d, init, endit !loop variables
  integer :: mod1, mod2 !which matrix elements to compute
  integer :: i, j, ii, jj, n, m, l!aux loop variables
  integer :: dimax, fdStep

  complex(kind=dp) :: Pele1, Pele2, Pele3, Pele4
  complex(kind=dp) :: sigma(2,2,3)

  real(kind=dp), allocatable, dimension(:,:) :: derivative
  complex(kind=dp), allocatable, dimension(:,:) :: der

  complex(kind=dp) :: alpha, beta
  complex(kind=dp), allocatable, dimension(:) :: Y
  complex(kind=dp) :: zdotc

  type(wavevector), allocatable, dimension(:) :: smallk
  type(sparse_matrix_T) :: HT_csr_mod1,HT_csr_mod2, der_csr
  complex(kind=dp), allocatable, dimension(:,:) :: hkp
  logical :: sparse
  integer :: nzmax, info
  type(matrix_descr) :: descrA

  dimax = size(cb_state,1)
  fdStep = int(dimax/8)

  Pele1 = 0.
  Pele2 = 0.
  Pele3 = 0.
  Pele4 = 0.

  if (present(dz)) then
    ! call FDmatrixDense(fdstep, dz, 1, 2, 0, derivative)
    ! ! der = dcmplx(derivative,0.0_dp)
    ! nzmax = (dimax+(dimax-1)*2)
    ! ! print *, nzmax
    ! call dnscsr_z_mkl(nzmax, dimax, dcmplx(derivative, 0.0_dp), der_csr)
    ! allocate(Y(dimax))
    ! deallocate(derivative)
    ! descrA%TYPE = SPARSE_MATRIX_TYPE_GENERAL
  end if

  sparse = .True.

  sigma = 0

  if (whichBand == 0) then

    print *, 'gfactor for conduction band'

    ii = 0  ! Initialize counter
    do n = bandIdx,bandIdx+1
      ii = ii+1
      jj = 0
      do m = bandIdx,bandIdx+1
        jj = jj+1
        sigma(ii,jj,1) = sigmaElem(cb_state(:,n), cb_state(:,m),'x',fdStep, startz,endz,dz,nlayers)
        sigma(ii,jj,2) = sigmaElem(cb_state(:,n), cb_state(:,m),'y',fdStep, startz,endz,dz,nlayers)
        sigma(ii,jj,3) = sigmaElem(cb_state(:,n), cb_state(:,m),'z',fdStep, startz,endz,dz,nlayers)
      end do
    end do
    print *, 'sigma'
    do i = 1, 2
      write(*,'(2(f9.4,1x,f9.4))') (sigma(i,j,1), j=1,2)
    end do
    print *, ' '
    do i = 1, 2
      write(*,'(2(f9.4,1x,f9.4))') (sigma(i,j,2), j=1,2)
    end do
    print *, ' '
    do i = 1, 2
      write(*,'(2(f9.4,1x,f9.4))') (sigma(i,j,3), j=1,2)
    end do
    ! stop

    do d = 1,3
      !compute sigma element, spin
      ii = 0
      jj = 0

      ! select which matrix element to compute
      if(d==1) then
        mod1 = 2
        mod2 = 3
      elseif ( d==2 ) then
        mod1 = 3
        mod2 = 1
      elseif ( d==3 ) then
        mod1 = 1
        mod2 = 2
      end if

      if (nlayers == 2 .and. sparse) then
        ! print *, 'sparse'
        if (allocated(smallk)) deallocate(smallk)
        allocate(smallk(1))
        allocate(hkp(dimax,dimax))

        if ( mod1==1 ) then
          smallk%kx = 1
          smallk%ky = 0
          smallk%kz = 0
        elseif ( mod1==2 ) then
          smallk%kx = 0
          smallk%ky = 1
          smallk%kz = 0
        elseif ( mod1==3 ) then
          smallk%kx = 0
          smallk%ky = 0
          smallk%kz = 1
        end if

        call ZB8bandQW(hkp, smallk(1), profile, kpterms, sparse=.True., HT_csr=HT_csr_mod1,g='g')

        if ( mod2==1 ) then
          smallk%kx = 1
          smallk%ky = 0
          smallk%kz = 0
        elseif ( mod2==2 ) then
          smallk%kx = 0
          smallk%ky = 1
          smallk%kz = 0
        elseif ( mod2==3 ) then
          smallk%kx = 0
          smallk%ky = 0
          smallk%kz = 1
        end if
        call ZB8bandQW(hkp, smallk(1), profile, kpterms, sparse=.True., HT_csr=HT_csr_mod2,g='g')

        if (allocated(hkp)) deallocate(hkp)

      end if

      if (nlayers == 2 .and. .not. sparse) then
        ! print *, 'sparse'
        if (allocated(smallk)) deallocate(smallk)
        allocate(smallk(1))
        if (allocated(hkp)) deallocate(hkp)
        allocate(hkp(dimax,dimax))

        if ( mod1==1 ) then
          smallk%kx = 1
          smallk%ky = 0
          smallk%kz = 0
        elseif ( mod1==2 ) then
          smallk%kx = 0
          smallk%ky = 1
          smallk%kz = 0
        elseif ( mod1==3 ) then
          smallk%kx = 0
          smallk%ky = 0
          smallk%kz = 1
        end if

        call ZB8bandQW(hkp, smallk(1), profile, kpterms, sparse=.False.,g='g')

        if ( mod2==1 ) then
          smallk%kx = 1
          smallk%ky = 0
          smallk%kz = 0
        elseif ( mod2==2 ) then
          smallk%kx = 0
          smallk%ky = 1
          smallk%kz = 0
        elseif ( mod2==3 ) then
          smallk%kx = 0
          smallk%ky = 0
          smallk%kz = 1
        end if
        call ZB8bandQW(hkp, smallk(1), profile, kpterms, sparse=.False.,g='g')

      end if


      ! tensor loop
      ii = 0
      do n=bandIdx,bandIdx+1
        ii = ii + 1
        jj = 0
        ! print *, 'n', n
        do m=bandIdx,bandIdx+1
          jj = jj+1

          ! print *, 'm', m

          !now run over numvb valence bands
          do l=1,numvb

            Pele1 = 0
            Pele2 = 0
            Pele3 = 0
            Pele4 = 0

            !have to compute <cb_n|p|vb_l> and <vb_l|p|cb_m> for mod1 and mod2
            !========Pele1========
            if (nlayers == 1) then
              call pMatrixEleCalc(Pele1, mod1, cb_state(:,n), vb_state(:,l), nlayers, params)
            else
              if(sparse) then
                call pMatrixEleCalc(Pele1, mod1, cb_state(:,n), vb_state(:,l), nlayers, params, profile=profile, kpterms=kpterms, HT_csr=HT_csr_mod1)
              else
                call pMatrixEleCalc(Pele1, mod1, cb_state(:,n), vb_state(:,l), nlayers, params, hkp_out=hkp, profile=profile, kpterms=kpterms, startz=startz, endz=endz,dz=dz)
              end if

              ! if ( mod1==3 ) then
              !   Y = 0
              !   info = mkl_sparse_z_mv (SPARSE_OPERATION_NON_TRANSPOSE, alpha, der_csr, descrA, vb_state(:,l), beta, Y)
              !   Pele1 = Pele1 -cmplx(0.,1.)*hbar*zdotc(dimax,cb_state(:,n),1,Y(:),1)
              ! end if
            end if
            !========Pele1========

            !========Pele2========
            if (nlayers == 1) then
              call pMatrixEleCalc(Pele2, mod2, vb_state(:,l), cb_state(:,m), nlayers, params)
            else
              if(sparse) then
                call pMatrixEleCalc(Pele2, mod2, vb_state(:,l), cb_state(:,m), nlayers, params, profile=profile, kpterms=kpterms, HT_csr=HT_csr_mod2)
              else
                call pMatrixEleCalc(Pele2, mod2, vb_state(:,l), cb_state(:,m), nlayers, params, hkp_out=hkp, profile=profile, kpterms=kpterms, startz=startz, endz=endz,dz=dz)
              end if

              ! if ( mod2==3 ) then
              !   Y = 0
              !   info = mkl_sparse_z_mv (SPARSE_OPERATION_NON_TRANSPOSE, alpha, der_csr, descrA, cb_state(:,m), beta, Y)
              !   Pele2 = Pele2 -cmplx(0.,1.)*hbar*zdotc(dimax,vb_state(:,l),1,Y(:),1)
              ! end if
            end if
            !========Pele2========

            !========Pele3========
            if (nlayers == 1) then
              call pMatrixEleCalc(Pele3, mod2, cb_state(:,n), vb_state(:,l), nlayers, params)
            else
              if(sparse) then
                call pMatrixEleCalc(Pele3, mod2, cb_state(:,n), vb_state(:,l), nlayers, params, hkp_out=hkp, profile=profile, kpterms=kpterms, HT_csr=HT_csr_mod2)
              else
                call pMatrixEleCalc(Pele3, mod2, cb_state(:,n), vb_state(:,l), nlayers, params, profile=profile, kpterms=kpterms, startz=startz, endz=endz,dz=dz)
              end if

              ! if ( mod2==3 ) then
              !   Y = 0
              !   info = mkl_sparse_z_mv (SPARSE_OPERATION_NON_TRANSPOSE, alpha, der_csr, descrA, vb_state(:,l), beta, Y)
              !   Pele3 = Pele3 -cmplx(0.,1.)*hbar*zdotc(dimax,cb_state(:,n),1,Y(:),1)
              ! end if
            end if
            !========Pele3========

            !========Pele4========
            if (nlayers == 1) then
              call pMatrixEleCalc(Pele4, mod1, vb_state(:,l), cb_state(:,m), nlayers, params)
            else
              if(sparse) then
                call pMatrixEleCalc(Pele4, mod1, vb_state(:,l), cb_state(:,m), nlayers, params, hkp_out=hkp, profile=profile, kpterms=kpterms, HT_csr=HT_csr_mod1)
              else
                call pMatrixEleCalc(Pele4, mod1, vb_state(:,l), cb_state(:,m), nlayers, params, profile=profile, kpterms=kpterms, startz=startz, endz=endz,dz=dz)
              end if

              ! if ( mod1==3 ) then
              !   Y = 0
              !   info = mkl_sparse_z_mv (SPARSE_OPERATION_NON_TRANSPOSE, alpha, der_csr, descrA, cb_state(:,m), beta, Y)
              !   Pele4 = Pele4 -cmplx(0.,1.)*hbar*zdotc(dimax,vb_state(:,l),1,Y(:),1)
              ! end if
            end if

            !========Pele4========


            tensor(ii,jj,d) = tensor(ii,jj,d) + ( (Pele1*Pele2-Pele3*Pele4) / ( ( cb_value(n)-vb_value(l) ) + ( cb_value(m)-vb_value(l) ) ) )

          end do

        end do
      end do

      if ( numcb > 2 ) then
        ii = 0
        do n=bandidx,bandidx+1
          ! print *, 'n', n
          ii = ii + 1
          jj = 0
          do m=bandidx,bandidx+1
            jj = jj + 1
            ! print *, 'm', m
            do l=1,numcb

              Pele1 = 0
              Pele2 = 0
              Pele3 = 0
              Pele4 = 0

              if ( ( cb_value(n)-cb_value(l) ) /= 0. .and. ( cb_value(m)-cb_value(l) ) /= 0. ) then

                !========Pele1========
                if (nlayers == 1) then
                  call pMatrixEleCalc(Pele1, mod1, cb_state(:,n), cb_state(:,l), nlayers, params)
                else
                  if(sparse) then
                    call pMatrixEleCalc(Pele1, mod1, cb_state(:,n), cb_state(:,l), nlayers, params, profile=profile, kpterms=kpterms, HT_csr=HT_csr_mod1)
                  else
                    call pMatrixEleCalc(Pele1, mod1, cb_state(:,n), cb_state(:,l), nlayers, params, profile=profile, kpterms=kpterms, startz=startz, endz=endz,dz=dz)
                  end if

                  ! if ( mod1==3 ) then
                  !   Y = 0
                  !   info = mkl_sparse_z_mv (SPARSE_OPERATION_NON_TRANSPOSE, alpha, der_csr, descrA, cb_state(:,l), beta, Y)
                  !   Pele1 = Pele1 -cmplx(0.,1.)*hbar*zdotc(dimax,cb_state(:,n),1,Y(:),1)
                  ! end if

                end if

                !========Pele1========

                !========Pele2========
                if (nlayers == 1) then
                  call pMatrixEleCalc(Pele2, mod2, cb_state(:,l), cb_state(:,m), nlayers, params)
                else
                  if(sparse) then
                    call pMatrixEleCalc(Pele2, mod2, cb_state(:,l), cb_state(:,m), nlayers, params, profile=profile, kpterms=kpterms, HT_csr=HT_csr_mod2)
                  else
                    call pMatrixEleCalc(Pele2, mod2, cb_state(:,l), cb_state(:,m), nlayers, params, profile=profile, kpterms=kpterms, startz=startz, endz=endz,dz=dz)
                  end if

                  ! if ( mod2==3 ) then
                  !   Y = 0
                  !   info = mkl_sparse_z_mv (SPARSE_OPERATION_NON_TRANSPOSE, alpha, der_csr, descrA, cb_state(:,m), beta, Y)
                  !   Pele2 = Pele2 -cmplx(0.,1.)*hbar*zdotc(dimax,cb_state(:,l),1,Y(:),1)
                  ! end if
                end if
                !========Pele2========

                !========Pele3========
                if (nlayers == 1) then
                  call pMatrixEleCalc(Pele3, mod2, cb_state(:,n), cb_state(:,l), nlayers, params)
                else
                  if(sparse) then
                    call pMatrixEleCalc(Pele3, mod2, cb_state(:,n), cb_state(:,l), nlayers, params, profile=profile, kpterms=kpterms, HT_csr=HT_csr_mod2)
                  else
                    call pMatrixEleCalc(Pele3, mod2, cb_state(:,n), cb_state(:,l), nlayers, params, profile=profile, kpterms=kpterms, startz=startz, endz=endz,dz=dz)
                  end if
                  !
                  ! if ( mod2==3 ) then
                  !   Y = 0
                  !   info = mkl_sparse_z_mv (SPARSE_OPERATION_NON_TRANSPOSE, alpha, der_csr, descrA, cb_state(:,l), beta, Y)
                  !   Pele3 = Pele3 -cmplx(0.,1.)*hbar*zdotc(dimax,cb_state(:,n),1,Y(:),1)
                  ! end if
                end if
                !========Pele3========

                !========Pele4========
                if (nlayers == 1) then
                  call pMatrixEleCalc(Pele4, mod1, cb_state(:,l), cb_state(:,m), nlayers, params)
                else
                  if(sparse) then
                    call pMatrixEleCalc(Pele4, mod1, cb_state(:,l), cb_state(:,m), nlayers, params, profile=profile, kpterms=kpterms, HT_csr=HT_csr_mod1)
                  else
                    call pMatrixEleCalc(Pele4, mod1, cb_state(:,l), cb_state(:,m), nlayers, params, profile=profile, kpterms=kpterms, startz=startz, endz=endz,dz=dz)
                  end if

                  ! if ( mod1==3 ) then
                  !   Y = 0
                  !   info = mkl_sparse_z_mv (SPARSE_OPERATION_NON_TRANSPOSE, alpha, der_csr, descrA, cb_state(:,m), beta, Y)
                  !   Pele4 = Pele4 -cmplx(0.,1.)*hbar*zdotc(dimax,cb_state(:,l),1,Y(:),1)
                  ! end if
                end if
                !========Pele4========

                tensor(ii,jj,d) = tensor(ii,jj,d) + ( (Pele1*Pele2-Pele3*Pele4) / ( ( cb_value(n)-cb_value(l) ) + ( cb_value(m)-cb_value(l) ) ) )

              end if
            end do
          end do
        end do
      end if
    end do ! end main loop
  end if

  tensor(1:2,1:2,1:3) = -cmplx(0.,1.)*tensor(1:2,1:2,1:3)/hbar2O2m0
  tensor(1:2,1:2,1:3) = tensor(1:2,1:2,1:3) - (2.00231_dp/2.0_dp)*sigma(1:2,1:2,1:3)


  if (allocated(Y)) deallocate(Y)

end subroutine gfactorCalculation



end module gfactorFunctions
