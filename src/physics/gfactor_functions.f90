module gfactorFunctions

  use definitions
  use hamiltonianConstructor
  use mkl_spblas
  use OMP_lib

  implicit none

  contains

subroutine set_perturbation_direction(d, smallk)
  integer, intent(in) :: d
  type(wavevector), intent(inout) :: smallk

  smallk%kx = 0
  smallk%ky = 0
  smallk%kz = 0
  select case(d)
  case(1); smallk%kx = 1
  case(2); smallk%ky = 1
  case(3); smallk%kz = 1
  case default
    print *, 'Error: perturbation direction must be 1, 2, or 3, got:', d
    stop 1
  end select
end subroutine set_perturbation_direction

subroutine compute_pele(Pele, direction, state_a, state_b, nlayers, params, &
  & sparse_mode, profile, kpterms, HT_csr_dir, hkp_mat, startz, endz, dz)

  complex(kind=dp), intent(inout) :: Pele
  integer, intent(in) :: direction
  complex(kind=dp), intent(in), dimension(:) :: state_a, state_b
  integer, intent(in) :: nlayers
  type(paramStruct), intent(in) :: params(nlayers)
  logical, intent(in) :: sparse_mode
  real(kind=dp), intent(in), optional, dimension(:,:) :: profile
  real(kind=dp), intent(in), optional, dimension(:,:,:) :: kpterms
  type(sparse_matrix_T), intent(in), optional :: HT_csr_dir
  complex(kind=dp), intent(in), optional, dimension(:,:) :: hkp_mat
  real(kind=dp), intent(in), optional :: startz, endz, dz

  if (nlayers == 1) then
    call pMatrixEleCalc(Pele, direction, state_a, state_b, nlayers, params)
  else
    if (sparse_mode) then
      call pMatrixEleCalc(Pele, direction, state_a, state_b, nlayers, params, &
        & profile=profile, kpterms=kpterms, HT_csr=HT_csr_dir)
    else
      call pMatrixEleCalc(Pele, direction, state_a, state_b, nlayers, params, &
        & hkp_out=hkp_mat, profile=profile, kpterms=kpterms, &
        & startz=startz, endz=endz, dz=dz)
    end if
  end if

end subroutine compute_pele

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

  ! Spin matrices in the 8-band zincblende basis.
  ! Basis ordering: 1=|3/2,+3/2> (HH), 2=|3/2,+1/2> (LH), 3=|3/2,-1/2> (LH),
  !   4=|3/2,-3/2> (HH), 5=|1/2,+1/2> (SO), 6=|1/2,-1/2> (SO),
  !   7=|S,+1/2> (CB), 8=|S,-1/2> (CB)
  ! Phase convention matches ZB8bandBulk in hamiltonianConstructor.f90.
  ! CB block (bands 7-8) gives standard Pauli matrices.
  ! VB-VB and VB-SO off-diagonal elements may differ from the Chuang & Chang
  ! convention by factors of +/-i (Winkler, Table 2.3).
  SIGMA_X(1,:) = (/ ZERO, -IU*RQS3, ZERO, ZERO, UM*SQR2o3, ZERO, ZERO, ZERO /)
  SIGMA_X(2,:) = (/ IU*RQS3, ZERO, -IU*2.0_dp/3.0_dp, ZERO, ZERO, -UM*SQR2/3.0_dp, ZERO, ZERO /)
  SIGMA_X(3,:) = (/ ZERO, 2.0_dp*IU/3.0_dp, ZERO, -IU*RQS3, UM*SQR2/3.0_dp, ZERO, ZERO, ZERO /)
  SIGMA_X(4,:) = (/ ZERO, ZERO, IU*RQS3, ZERO, ZERO, -UM*SQR2o3, ZERO, ZERO /)
  SIGMA_X(5,:) = (/ UM*SQR2o3, ZERO, UM*SQR2/3.0_dp, ZERO, ZERO, -IU/3.0_dp, ZERO, ZERO /)
  SIGMA_X(6,:) = (/ ZERO, -UM*SQR2/3.0_dp, ZERO, -UM*SQR2o3, IU/3.0_dp, ZERO, ZERO, ZERO /)
  SIGMA_X(7,:) = (/ ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, UM /)
  SIGMA_X(8,:) = (/ ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, UM, ZERO /)


  SIGMA_Y(1,:) = (/ ZERO, UM*RQS3, ZERO, ZERO, IU*SQR2o3, ZERO, ZERO, ZERO /)
  SIGMA_Y(2,:) = (/ UM*RQS3, ZERO, UM*2.0_dp/3.0_dp, ZERO, ZERO, -IU*SQR2/3.0_dp, ZERO, ZERO /)
  SIGMA_Y(3,:) = (/ ZERO, UM*2.0_dp/3.0_dp, ZERO, UM*RQS3, -IU*SQR2/3.0_dp, ZERO, ZERO, ZERO /)
  SIGMA_Y(4,:) = (/ ZERO, ZERO, UM*RQS3, ZERO, ZERO, IU*SQR2o3, ZERO, ZERO /)
  SIGMA_Y(5,:) = (/ -IU*SQR2o3, ZERO, IU*SQR2/3.0_dp, ZERO, ZERO, UM/3.0_dp, ZERO, ZERO /)
  SIGMA_Y(6,:) = (/ ZERO, IU*SQR2/3.0_dp, ZERO, -IU*SQR2o3, UM/3.0_dp, ZERO, ZERO, ZERO /)
  SIGMA_Y(7,:) = (/ ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, IU /)
  SIGMA_Y(8,:) = (/ ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, -IU, ZERO /)


  SIGMA_Z(1,:) = (/ UM, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO /)
  SIGMA_Z(2,:) = (/ ZERO, UM/3.0_dp, ZERO, ZERO, -2.0_dp*IU*SQR2/3.0_dp, ZERO, ZERO, ZERO /)
  SIGMA_Z(3,:) = (/ ZERO, ZERO, -UM/3.0_dp, ZERO, ZERO, 2.0_dp*IU*SQR2/3.0_dp, ZERO, ZERO /)
  SIGMA_Z(4,:) = (/ ZERO, ZERO, ZERO, -UM, ZERO, ZERO, ZERO, ZERO /)
  SIGMA_Z(5,:) = (/ ZERO, 2.0_dp*IU*SQR2/3.0_dp, ZERO, ZERO, -UM/3.0_dp, ZERO, ZERO, ZERO /)
  SIGMA_Z(6,:) = (/ ZERO, ZERO, -2.0_dp*IU*SQR2/3.0_dp, ZERO, ZERO, UM/3.0_dp, ZERO, ZERO /)
  SIGMA_Z(7,:) = (/ ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, UM, ZERO /)
  SIGMA_Z(8,:) = (/ ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, -UM /)



  allocate(v(fdstep))
  v = 0
  alpha = cmplx(1.0_dp, 0.0_dp, kind=dp)
  beta = cmplx(0.0_dp, 0.0_dp, kind=dp)
  sigmaElem = 0
  aux_v = 0
  aux1_v = 0

  if (nlayers > 1) then

    if (abs(dz) < tolerance) then
      print *, 'Error: sigmaElem called with dz=0. Grid spacing is zero.'
      stop 1
    end if

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
  integer :: info

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
    hkp(:,:) = cmplx(0, 0, kind=dp)
  end if

  if (present(hkp_out)) hkp = hkp_out

  alpha = cmplx(1.0_dp, 0.0_dp, kind=dp)
  beta = cmplx(0.0_dp, 0.0_dp, kind=dp)

  call set_perturbation_direction(d, smallk(1))
  if (nlayers == 1) then
    call ZB8bandBulk(hkp,smallk(1),params(1),g='g')
  end if
  if (nlayers > 1 .and. .not. present(HT_csr) .and. .not. present(hkp_out)) then
    call ZB8bandQW(hkp, smallk(1), profile, kpterms,g='g')
  end if


  if(present(Ppartial)) then
     call zgemv('N',dimax,dimax,alpha,hkp(:,:),dimax,state2(:),1,beta,Ppartial,1)
  elseif ( present(HT_csr) ) then


    allocate(Y1(dimax))
    !create matrix descriptor
    descrA%TYPE =SPARSE_MATRIX_TYPE_GENERAL

    info = mkl_sparse_z_mv (SPARSE_OPERATION_NON_TRANSPOSE, alpha, HT_csr, descrA, state2(:), beta, Y1)

    if (info /= 0) then
      print *, 'Error: mkl_sparse_z_mv failed with info =', info
      stop 1
    end if

    Pele = zdotc(dimax,state1(:),1,Y1(:),1)

    deallocate(Y1)

  else
    allocate(v(fdstep))
    v = 0
    aux_v = 0
    aux1_v = 0

    !$OMP PARALLEL DO private(i, n, m, aux, aux_v, aux1_v, Y) shared(v, hkp, fdstep, state1, state2, alpha, beta)
    do i = 1, fdstep

      do n = 1, 8, 1
        do m = 1, 8, 1
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
  integer :: skip_count
  integer :: mod1, mod2 !which matrix elements to compute
  integer :: i, j, ii, jj, n, m, l!aux loop variables
  integer :: dimax, fdStep
  real(kind=dp) :: denom

  complex(kind=dp) :: Pele1, Pele2, Pele3, Pele4
  complex(kind=dp) :: sigma(2,2,3)


  complex(kind=dp) :: alpha, beta
  complex(kind=dp), allocatable, dimension(:) :: Y
  complex(kind=dp) :: zdotc

  type(wavevector), allocatable, dimension(:) :: smallk
  type(sparse_matrix_T) :: HT_csr_mod1,HT_csr_mod2
  complex(kind=dp), allocatable, dimension(:,:) :: hkp
  logical :: sparse
  integer :: info
  type(matrix_descr) :: descrA

  dimax = size(cb_state,1)
  fdStep = int(dimax/8)

  Pele1 = 0.
  Pele2 = 0.
  Pele3 = 0.
  Pele4 = 0.

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

    skip_count = 0
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

      if (nlayers > 1 .and. sparse) then
        if (allocated(smallk)) deallocate(smallk)
        allocate(smallk(1))
        allocate(hkp(dimax,dimax))

        call set_perturbation_direction(mod1, smallk(1))
        call ZB8bandQW(hkp, smallk(1), profile, kpterms, sparse=.True., HT_csr=HT_csr_mod1,g='g')

        call set_perturbation_direction(mod2, smallk(1))
        call ZB8bandQW(hkp, smallk(1), profile, kpterms, sparse=.True., HT_csr=HT_csr_mod2,g='g')

        if (allocated(hkp)) deallocate(hkp)

      end if

      ! tensor loop
      ii = 0
      do n=bandIdx,bandIdx+1
        ii = ii + 1
        jj = 0
        do m=bandIdx,bandIdx+1
          jj = jj+1


          !now run over numvb valence bands
          do l=1,numvb

            Pele1 = 0
            Pele2 = 0
            Pele3 = 0
            Pele4 = 0

            !have to compute <cb_n|p|vb_l> and <vb_l|p|cb_m> for mod1 and mod2
            call compute_pele(Pele1, mod1, cb_state(:,n), vb_state(:,l), nlayers, params, &
              & sparse, profile, kpterms, HT_csr_mod1, hkp, startz, endz, dz)

            call compute_pele(Pele2, mod2, vb_state(:,l), cb_state(:,m), nlayers, params, &
              & sparse, profile, kpterms, HT_csr_mod2, hkp, startz, endz, dz)

            call compute_pele(Pele3, mod2, cb_state(:,n), vb_state(:,l), nlayers, params, &
              & sparse, profile, kpterms, HT_csr_mod2, hkp, startz, endz, dz)

            call compute_pele(Pele4, mod1, vb_state(:,l), cb_state(:,m), nlayers, params, &
              & sparse, profile, kpterms, HT_csr_mod1, hkp, startz, endz, dz)


            denom = (cb_value(n) - vb_value(l)) + (cb_value(m) - vb_value(l))
            if (abs(denom) > tolerance) then
              tensor(ii, jj, d) = tensor(ii, jj, d) + (Pele1*Pele2 - Pele3*Pele4) / denom
            else
              skip_count = skip_count + 1
              if (skip_count <= 5) then
                print *, 'WARNING: near-zero energy denominator, n=', n, 'm=', m, 'l=', l, 'denom=', denom
              end if
            end if

          end do

        end do
      end do

      ii = 0
      do n=bandidx,bandidx+1
        ii = ii + 1
        jj = 0
        do m=bandidx,bandidx+1
          jj = jj + 1
          do l=1,numcb

              ! Skip self-interaction (l==n or l==m gives zero contribution)
              if (l == n .or. l == m) cycle

              Pele1 = 0
              Pele2 = 0
              Pele3 = 0
              Pele4 = 0

                call compute_pele(Pele1, mod1, cb_state(:,n), cb_state(:,l), nlayers, params, &
                  & sparse, profile, kpterms, HT_csr_mod1, hkp, startz, endz, dz)

                call compute_pele(Pele2, mod2, cb_state(:,l), cb_state(:,m), nlayers, params, &
                  & sparse, profile, kpterms, HT_csr_mod2, hkp, startz, endz, dz)

                call compute_pele(Pele3, mod2, cb_state(:,n), cb_state(:,l), nlayers, params, &
                  & sparse, profile, kpterms, HT_csr_mod2, hkp, startz, endz, dz)

                call compute_pele(Pele4, mod1, cb_state(:,l), cb_state(:,m), nlayers, params, &
                  & sparse, profile, kpterms, HT_csr_mod1, hkp, startz, endz, dz)

                denom = (cb_value(n) - cb_value(l)) + (cb_value(m) - cb_value(l))
                if (abs(denom) > tolerance) then
                  tensor(ii, jj, d) = tensor(ii, jj, d) + (Pele1*Pele2 - Pele3*Pele4) / denom
                else
                  skip_count = skip_count + 1
                  if (skip_count <= 5) then
                    print *, 'WARNING: near-zero energy denominator, n=', n, 'm=', m, 'l=', l, 'denom=', denom
                  end if
                end if

            end do
          end do
        end do
    if (skip_count > 0) print *, 'Total skipped contributions:', skip_count
    end do ! end main loop

  else if (whichBand == 1) then

    print *, 'gfactor for valence band'

    ! Sigma matrix between VB doublet (Kramers partners)
    ii = 0
    do n = bandIdx,bandIdx+1
      ii = ii+1
      jj = 0
      do m = bandIdx,bandIdx+1
        jj = jj+1
        sigma(ii,jj,1) = sigmaElem(vb_state(:,n), vb_state(:,m),'x',fdStep, startz,endz,dz,nlayers)
        sigma(ii,jj,2) = sigmaElem(vb_state(:,n), vb_state(:,m),'y',fdStep, startz,endz,dz,nlayers)
        sigma(ii,jj,3) = sigmaElem(vb_state(:,n), vb_state(:,m),'z',fdStep, startz,endz,dz,nlayers)
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

    skip_count = 0
    do d = 1,3
      ii = 0
      jj = 0

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

      if (nlayers > 1 .and. sparse) then
        if (allocated(smallk)) deallocate(smallk)
        allocate(smallk(1))
        allocate(hkp(dimax,dimax))

        call set_perturbation_direction(mod1, smallk(1))
        call ZB8bandQW(hkp, smallk(1), profile, kpterms, sparse=.True., HT_csr=HT_csr_mod1,g='g')

        call set_perturbation_direction(mod2, smallk(1))
        call ZB8bandQW(hkp, smallk(1), profile, kpterms, sparse=.True., HT_csr=HT_csr_mod2,g='g')

        if (allocated(hkp)) deallocate(hkp)
      end if

      ! VB doublet (n,m) with CB intermediates (l)
      ii = 0
      do n=bandIdx,bandIdx+1
        ii = ii + 1
        jj = 0
        do m=bandIdx,bandIdx+1
          jj = jj+1

          do l=1,numcb

            Pele1 = 0
            Pele2 = 0
            Pele3 = 0
            Pele4 = 0

            call compute_pele(Pele1, mod1, vb_state(:,n), cb_state(:,l), nlayers, params, &
              & sparse, profile, kpterms, HT_csr_mod1, hkp, startz, endz, dz)

            call compute_pele(Pele2, mod2, cb_state(:,l), vb_state(:,m), nlayers, params, &
              & sparse, profile, kpterms, HT_csr_mod2, hkp, startz, endz, dz)

            call compute_pele(Pele3, mod2, vb_state(:,n), cb_state(:,l), nlayers, params, &
              & sparse, profile, kpterms, HT_csr_mod2, hkp, startz, endz, dz)

            call compute_pele(Pele4, mod1, cb_state(:,l), vb_state(:,m), nlayers, params, &
              & sparse, profile, kpterms, HT_csr_mod1, hkp, startz, endz, dz)

            denom = (vb_value(n) - cb_value(l)) + (vb_value(m) - cb_value(l))
            if (abs(denom) > tolerance) then
              tensor(ii, jj, d) = tensor(ii, jj, d) + (Pele1*Pele2 - Pele3*Pele4) / denom
            else
              skip_count = skip_count + 1
              if (skip_count <= 5) then
                print *, 'WARNING: near-zero energy denominator, n=', n, 'm=', m, 'l=', l, 'denom=', denom
              end if
            end if

          end do

        end do
      end do

      ! VB doublet (n,m) with VB intermediates (l), skip self-term
      ii = 0
      do n=bandIdx,bandIdx+1
        ii = ii + 1
        jj = 0
        do m=bandIdx,bandIdx+1
          jj = jj+1
          do l=1,numvb

            if (l == n .or. l == m) cycle

            Pele1 = 0
            Pele2 = 0
            Pele3 = 0
            Pele4 = 0

            call compute_pele(Pele1, mod1, vb_state(:,n), vb_state(:,l), nlayers, params, &
              & sparse, profile, kpterms, HT_csr_mod1, hkp, startz, endz, dz)

            call compute_pele(Pele2, mod2, vb_state(:,l), vb_state(:,m), nlayers, params, &
              & sparse, profile, kpterms, HT_csr_mod2, hkp, startz, endz, dz)

            call compute_pele(Pele3, mod2, vb_state(:,n), vb_state(:,l), nlayers, params, &
              & sparse, profile, kpterms, HT_csr_mod2, hkp, startz, endz, dz)

            call compute_pele(Pele4, mod1, vb_state(:,l), vb_state(:,m), nlayers, params, &
              & sparse, profile, kpterms, HT_csr_mod1, hkp, startz, endz, dz)

            denom = (vb_value(n) - vb_value(l)) + (vb_value(m) - vb_value(l))
            if (abs(denom) > tolerance) then
              tensor(ii, jj, d) = tensor(ii, jj, d) + (Pele1*Pele2 - Pele3*Pele4) / denom
            else
              skip_count = skip_count + 1
              if (skip_count <= 5) then
                print *, 'WARNING: near-zero energy denominator, n=', n, 'm=', m, 'l=', l, 'denom=', denom
              end if
            end if

          end do
        end do
      end do
    if (skip_count > 0) print *, 'Total skipped contributions:', skip_count
    end do ! end main loop
  end if

  tensor(1:2,1:2,1:3) = -cmplx(0.0_dp, 1.0_dp, kind=dp)*tensor(1:2,1:2,1:3)/hbar2O2m0
  tensor(1:2,1:2,1:3) = tensor(1:2,1:2,1:3) - (2.00231_dp/2.0_dp)*sigma(1:2,1:2,1:3)


end subroutine gfactorCalculation


end module gfactorFunctions