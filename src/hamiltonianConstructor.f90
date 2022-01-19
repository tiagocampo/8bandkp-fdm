module hamiltonianConstructor

  use definitions
  use finitedifferences
  use utils
  use mkl_spblas

  implicit none

  contains

    subroutine externalFiledSetup_electricField(profile, Evalue, totalSize, z)

      real(kind = dp), intent(inout), allocatable, dimension(:,:) :: profile
      real(kind = dp), intent(in) :: Evalue, totalSize
      real(kind = dp), intent(in), dimension(:) :: z

      integer :: i

      do i = 1, ubound(z, dim=1), 1
          profile(i,:) = profile(i,:) - (Evalue*totalSize) * (z(i)+z(1))/(2.0_dp*z(1))
      end do

    end subroutine externalFiledSetup_electricField

    subroutine confinementInitialization(z, startPos, endPos, material, nlayers,&
      & params, bshift, confDir, profile, kpterms)

      real(kind = dp), intent(in), dimension(:) :: z, bshift
      integer, intent(in), dimension(:) :: startPos, endPos
      character(len = 255), intent(in) :: material(nlayers)
      integer, intent(in) :: nlayers
      type(paramStruct), intent(in) :: params(nlayers)
      character(len = 1), intent(in) :: confDir
      real(kind = dp), intent(inout), allocatable, dimension(:,:) :: profile
      real(kind = dp), intent(inout), dimension(:,:,:) :: kpterms

      real(kind=dp), allocatable, dimension(:,:) :: ScnDer, FstDer, kptermsProfile
      real(kind=dp), allocatable, dimension(:,:) :: forward, central, backward
      real(kind=dp), allocatable, dimension(:) :: diag, offup, offdown

      integer :: i, initIDX, endIDX, N, ii, jj
      real(kind = dp) :: delta

      ! complex(kind=dp) :: zgemv

      N = size(z, dim=1)
      delta = abs(z(2) - z(1))
      print *, delta

      if (allocated(profile)) deallocate(profile)
      allocate(profile(N,3))

      if (allocated(kptermsProfile)) deallocate(kptermsProfile)
      allocate(kptermsProfile(N,5))

      ! allocate(ScnDer(N,N))
      ! allocate(Fstder(N,N))
      ! call FDmatrixDense(N, delta, 2, 2, 1, ScnDer)
      ! call FDmatrixDense(N, delta, 1, 2, 1, FstDer)
      !call Identity(N, Idn)

      if (confDir == 'z') then

        do i = 1, nlayers, 1

          ! profile(startPos(i):endPos(i),1) = (params(1)%Eg - params(i)%Eg)*bshift(i)
          ! profile(startPos(i):endPos(i),2) = (params(1)%Eg - params(i)%Eg)*bshift(i) - params(i)%DeltaSO
          ! profile(startPos(i):endPos(i),3) =  params(1)%Eg + (params(i)%Eg - params(1)%Eg)*(1.0_dp - bshift(i))

          profile(startPos(i):endPos(i),1) = params(i)%EV
          profile(startPos(i):endPos(i),2) = params(i)%EV - params(i)%DeltaSO
          profile(startPos(i):endPos(i),3) =  params(i)%EC

          kptermsProfile(startPos(i):endPos(i),1) = params(i)%gamma1
          kptermsProfile(startPos(i):endPos(i),2) = params(i)%gamma2
          kptermsProfile(startPos(i):endPos(i),3) = params(i)%gamma3
          kptermsProfile(startPos(i):endPos(i),4) = params(i)%A
          kptermsProfile(startPos(i):endPos(i),5) = params(i)%P

        end do

      else

        stop "confinement not implemented for directions other than z"

      end if

      forall(ii=1:N)
        kpterms(ii,ii,1) = kptermsProfile(ii,1) !gamma1
        kpterms(ii,ii,2) = kptermsProfile(ii,2) !gamma2
        kpterms(ii,ii,3) = kptermsProfile(ii,3) !gamma3
        kpterms(ii,ii,4) = kptermsProfile(ii,5) !P
        kpterms(ii,ii,10) = kptermsProfile(ii,4) !A
      end forall

      allocate(forward(N,N))
      allocate(backward(N,N))
      allocate(central(N,N))
      forward = 0.0_dp
      backward = 0.0_dp
      central = 0.0_dp

      forall(ii=1:N-1)
        forward(ii,ii) = 1
        forward(ii,ii+1) = 1
      end forall
      !last element done alone to avoid if inside loop
      forward(N,N) = 1
      backward = transpose(forward)
      central = backward + forward

      ! do ii = 1, N, 1
      !   write(104,*) (forward(ii,jj), jj=1,N)
      !   write(105,*) (backward(ii,jj), jj=1,N)
      !   write(106,*) (central(ii,jj), jj=1,N)
      ! end do
      ! stop

      allocate(diag(N))
      allocate(offup(N))
      allocate(offdown(N))

      ! A*kz**2
      call dgemv('N', N, N, 1.0_dp, central, N, kptermsProfile(1:N,4), 1, 0.0_dp, diag, 1)
      call dgemv('N', N, N, 1.0_dp, forward, N, kptermsProfile(1:N,4), 1, 0.0_dp, offup, 1)
      call dgemv('N', N, N, 1.0_dp, backward, N, kptermsProfile(1:N,4), 1, 0.0_dp, offdown, 1)
      forall(ii=2:N-1)
        kpterms(ii,ii,5) = diag(ii)
        kpterms(ii+1,ii,5) = -offup(ii)
        kpterms(ii-1,ii,5) = -offdown(ii)
      end forall
      kpterms(1,1,5) = diag(1)
      kpterms(N,N,5) = diag(N)
      kpterms(2,1,5) = -offup(1)
      kpterms(N-1,N,5) = -offdown(N)
      kpterms(:,:,5) = kpterms(:,:,5)*(1.0_dp/(2.0_dp*delta**2))

      ! do ii = 1, N, 1
      !   write(104,*) (kpterms(ii,jj,5), jj=1,N)
      ! end do

      !Q
      call dgemv('N', N, N, 1.0_dp, central, N, kptermsProfile(1:N,1) - 2.0_dp*kptermsProfile(1:N,2), 1, 0.0_dp, diag, 1)
      call dgemv('N', N, N, 1.0_dp, forward, N, kptermsProfile(1:N,1) - 2.0_dp*kptermsProfile(1:N,2), 1, 0.0_dp, offup, 1)
      call dgemv('N', N, N, 1.0_dp, backward, N, kptermsProfile(1:N,1) - 2.0_dp*kptermsProfile(1:N,2), 1, 0.0_dp, offdown, 1)
      forall(ii=2:N-1)
        kpterms(ii,ii,7) = diag(ii)
        kpterms(ii+1,ii,7) = -offup(ii)
        kpterms(ii-1,ii,7) = -offdown(ii)
      end forall
      kpterms(1,1,7) = diag(1)
      kpterms(N,N,7) = diag(N)
      kpterms(2,1,7) = -offup(1)
      kpterms(N-1,N,7) = -offdown(N)
      kpterms(:,:,7) = kpterms(:,:,7)*(1.0_dp/(2.0_dp*delta**2))

      ! do ii = 1, N, 1
      !   write(105,*) (kpterms(ii,jj,7), jj=1,N)
      ! end do

      !T
      call dgemv('N', N, N, 1.0_dp, central, N, kptermsProfile(1:N,1) + 2.0_dp*kptermsProfile(1:N,2), 1, 0.0_dp, diag, 1)
      call dgemv('N', N, N, 1.0_dp, forward, N, kptermsProfile(1:N,1) + 2.0_dp*kptermsProfile(1:N,2), 1, 0.0_dp, offup, 1)
      call dgemv('N', N, N, 1.0_dp, backward, N, kptermsProfile(1:N,1) + 2.0_dp*kptermsProfile(1:N,2), 1, 0.0_dp, offdown, 1)
      forall(ii=2:N-1)
        kpterms(ii,ii,8) = diag(ii)
        kpterms(ii+1,ii,8) = -offup(ii)
        kpterms(ii-1,ii,8) = -offdown(ii)
      end forall
      kpterms(1,1,8) = diag(1)
      kpterms(N,N,8) = diag(N)
      kpterms(2,1,8) = -offup(1)
      kpterms(N-1,N,8) = -offdown(N)
      kpterms(:,:,8) = kpterms(:,:,8)*(1.0_dp/(2.0_dp*delta**2))

      ! do ii = 1, N, 1
      !   write(106,*) (kpterms(ii,jj,8), jj=1,N)
      ! end do


      ! P*kz
      call dgemv('N', N, N, 1.0_dp, forward, N, kptermsProfile(1:N,5), 1, 0.0_dp, offup, 1)
      call dgemv('N', N, N, 1.0_dp, backward, N, kptermsProfile(1:N,5), 1, 0.0_dp, offdown, 1)
      forall(ii=2:N-1)
        kpterms(ii+1,ii,6) = offup(ii)
        kpterms(ii-1,ii,6) = -offdown(ii)
      end forall
      kpterms(2,1,6) = offup(1)
      kpterms(N-1,N,6) = -offdown(N)
      kpterms(:,:,6) = kpterms(:,:,6)*(1.0_dp/(4.0_dp*delta))

      ! do ii = 1, N, 1
      !   write(107,*) (kpterms(ii,jj,6), jj=1,N)
      ! end do

      ! S -> gamma3*kz
      call dgemv('N', N, N, 1.0_dp, forward, N, kptermsProfile(1:N,3), 1, 0.0_dp, offup, 1)
      call dgemv('N', N, N, 1.0_dp, backward, N, kptermsProfile(1:N,3), 1, 0.0_dp, offdown, 1)
      forall(ii=2:N-1)
        kpterms(ii+1,ii,9) = offup(ii)
        kpterms(ii-1,ii,9) = -offdown(ii)
      end forall
      kpterms(2,1,9) = offup(1)
      kpterms(N-1,N,9) = -offdown(N)
      kpterms(:,:,9) = kpterms(:,:,9)*(1.0_dp/(4.0_dp*delta))

      ! do ii = 1, N, 1
      !   write(108,*) (kpterms(ii,jj,9), jj=1,N)
      ! end do
      !
      ! stop

      ! forall(ii=1:N, jj=1:N)
      !   kpterms(ii,jj,5) = -ScnDer(ii,jj)*kptermsProfile(ii,4) !A*kz**2
      !   kpterms(ii,jj,6) = -FstDer(ii,jj)*kptermsProfile(ii,5) !Pz
      !   kpterms(ii,jj,7) = -ScnDer(ii,jj)*(kptermsProfile(ii,1) - 2.0_dp*kptermsProfile(ii,1)) !Q
      !   kpterms(ii,jj,8) = -ScnDer(ii,jj)*(kptermsProfile(ii,1) + 2.0_dp*kptermsProfile(ii,1)) !T
      !   kpterms(ii,jj,9) = -FstDer(ii,jj)*kptermsProfile(ii,3) !S
      ! end forall

      if (allocated(ScnDer)) deallocate(ScnDer)
      if (allocated(FstDer)) deallocate(FstDer)
      if (allocated(kptermsProfile)) deallocate(kptermsProfile)
      if (allocated(forward)) deallocate(forward)
      if (allocated(backward)) deallocate(backward)
      if (allocated(central)) deallocate(central)
      if (allocated(diag)) deallocate(diag)
      if (allocated(offup)) deallocate(offup)
      if (allocated(offdown)) deallocate(offdown)

    end subroutine confinementInitialization

    subroutine ZB8bandQW(HT, wv, profile, kpterms, sparse, HT_csr, g)

      implicit none

      !input/output
      complex(kind=dp), intent(inout), dimension(:,:) :: HT
      type(sparse_matrix_T), intent(inout), optional :: HT_csr
      type(wavevector), intent(in) :: wv
      real(kind = dp), intent(in), dimension(:,:) :: profile
      real(kind = dp), intent(in), dimension(:,:,:) :: kpterms
      logical, intent(in), optional:: sparse
      character(len=1), intent(in), optional :: g

      integer :: nzmax

      ! wave vector related
      real(kind=dp) :: kx2, ky2, kxky, k2
      complex(kind=dp) :: kminus, kplus

      complex(kind=dp), allocatable, dimension(:,:) :: Q, R, RC, S, SC, T, PZ, PM, PP, A

      ! constants
      ! complex(kind=dp) :: IU
      ! real(kind=dp) :: SQR3, RQS2, SQR2, RQS3

      integer :: i, j, N, ii, jj

      N = size(HT, dim=1)/8

      !constants
      ! IU = dcmplx(0.0_dp,1.0_dp)
      ! SQR3 = dsqrt(3.0_dp)
      ! SQR2 = dsqrt(2.0_dp)
      ! RQS2 = 1.0_dp/dsqrt(2.0_dp)
      ! RQS3 = 1.0_dp/dsqrt(3.0_dp)


      !wave vectors
      kxky = wv%kx*wv%ky
      kx2 = wv%kx**2
      ky2 = wv%ky**2
      k2 = kx2 + ky2
      kplus   = wv%kx + IU*wv%ky
      kminus  = wv%kx - IU*wv%ky

      allocate(Q(N,N))
      allocate(T(N,N))
      allocate(S(N,N))
      allocate(SC(N,N))
      allocate(R(N,N))
      allocate(RC(N,N))
      allocate(PP(N,N))
      allocate(PM(N,N))
      allocate(PZ(N,N))
      allocate(A(N,N))
      Q = 0.0_dp
      T = 0.0_dp
      S = 0.0_dp
      SC = 0.0_dp
      R = 0.0_dp
      RC = 0.0_dp
      PP = 0.0_dp
      PM = 0.0_dp
      PZ = 0.0_dp
      A = 0.0_dp

      if (.not. present(g)) then
        forall(ii=1:N, jj=1:N)
          Q(ii,jj) = -( (kpterms(ii,jj,1) + kpterms(ii,jj,2))*k2 + kpterms(ii,jj,7) )
          T(ii,jj) = -( (kpterms(ii,jj,1) - kpterms(ii,jj,2))*k2 + kpterms(ii,jj,8) )
          S(ii,jj)  =  2.0_dp * SQR3 * kminus * kpterms(ii,jj,9) !no IU because of derivative
          SC(jj,ii) =  2.0_dp * SQR3 * kplus  * kpterms(ii,jj,9) !no IU because of derivative
          PZ(ii,jj) = kpterms(ii,jj,6) * (-IU)
          A(ii,jj)  = kpterms(ii,jj,5) + k2*kpterms(ii,jj,10)
        end forall

        forall (ii=1:N)
          R(ii,ii)  = - SQR3 * ( kpterms(ii,ii,2)*(kx2 - ky2) - 2.0_dp*IU*kpterms(ii,ii,3)*kxky )
          RC(ii,ii) = - SQR3 * ( kpterms(ii,ii,2)*(kx2 - ky2) + 2.0_dp*IU*kpterms(ii,ii,3)*kxky )
          PP(ii,ii) = kpterms(ii,ii,4) * kplus  * RQS2
          PM(ii,ii) = kpterms(ii,ii,4) * kminus * RQS2
        end forall
      else
        ! print *, wv%kz, kplus, kminus
        forall(ii=1:N, jj=1:N)
          PZ(ii,jj) = kpterms(ii,jj,4) * wv%kz
        end forall
        forall (ii=1:N)
          PP(ii,ii) = kpterms(ii,ii,4) * kplus  * RQS2
          PM(ii,ii) = kpterms(ii,ii,4) * kminus * RQS2
        end forall

      end if


      HT = 0

      ! kp matrix
      !col 1
      HT(1 + 0*N : 1*N, 1 + 0*N : 1*N) =  Q
      HT(1 + 0*N : 1*N, 1 + 1*N : 2*N) =  SC
      HT(1 + 0*N : 1*N, 1 + 2*N : 3*N) =  RC
      ! HT(1 + 0*N : 1*N, 1 + 3*N : 4*N) =  0.0_dp
      HT(1 + 0*N : 1*N, 1 + 4*N : 5*N) = -IU * RQS2 * SC
      HT(1 + 0*N : 1*N, 1 + 5*N : 6*N) =  IU * SQR2 * RC
      HT(1 + 0*N : 1*N, 1 + 6*N : 7*N) =  IU * PP
      ! HT(1 + 0*N : 1*N, 1 + 7*N : 8*N) =  0.0_dp

      !col 2
      HT(1 + 1*N : 2*N, 1 + 0*N : 1*N) =  (S)
      HT(1 + 1*N : 2*N, 1 + 1*N : 2*N) =  T
      ! HT(1 + 1*N : 2*N, 1 + 2*N : 3*N) =  0.0_dp
      HT(1 + 1*N : 2*N, 1 + 3*N : 4*N) =  RC
      HT(1 + 1*N : 2*N, 1 + 4*N : 5*N) =  IU * RQS2 * (Q - T)
      HT(1 + 1*N : 2*N, 1 + 5*N : 6*N) = -IU * SQR3 * RQS2 * SC
      HT(1 + 1*N : 2*N, 1 + 6*N : 7*N) =  SQR2 * RQS3 * PZ
      HT(1 + 1*N : 2*N, 1 + 7*N : 8*N) = -RQS3*PP

      !col 3
      HT(1 + 2*N : 3*N, 1 + 0*N : 1*N) =  R
      !HT(3,2) =  0.0_dp
      HT(1 + 2*N : 3*N, 1 + 2*N : 3*N) =  T
      HT(1 + 2*N : 3*N, 1 + 3*N : 4*N) = -SC
      HT(1 + 2*N : 3*N, 1 + 4*N : 5*N) =  IU * SQR3 * RQS2 * S
      HT(1 + 2*N : 3*N, 1 + 5*N : 6*N) =  IU * RQS2 * (Q - T)
      HT(1 + 2*N : 3*N, 1 + 6*N : 7*N) =  IU * RQS3 * PM
      HT(1 + 2*N : 3*N, 1 + 7*N : 8*N) =  IU * SQR2 * RQS3 * PZ

      !col 4
      !HT(4,1) =  0.0_dp
      HT(1 + 3*N : 4*N, 1 + 1*N : 2*N) =  R
      HT(1 + 3*N : 4*N, 1 + 2*N : 3*N) = (-S)
      HT(1 + 3*N : 4*N, 1 + 3*N : 4*N) =  Q
      HT(1 + 3*N : 4*N, 1 + 4*N : 5*N) =  IU * SQR2 * R
      HT(1 + 3*N : 4*N, 1 + 5*N : 6*N) =  IU * RQS2 * S
      ! HT(1 + 3*N : 4*N, 1 + 6*N : 7*N) =  0.0_dp
      HT(1 + 3*N : 4*N, 1 + 7*N : 8*N) = -PM

      !col 5
      HT(1 + 4*N : 5*N, 1 + 0*N : 1*N) =  IU * RQS2 * (S)
      HT(1 + 4*N : 5*N, 1 + 1*N : 2*N) = -IU * RQS2 * (Q - T)
      HT(1 + 4*N : 5*N, 1 + 2*N : 3*N) = -IU * SQR3 * RQS2 * (SC)
      HT(1 + 4*N : 5*N, 1 + 3*N : 4*N) = -IU * SQR2 * RC
      HT(1 + 4*N : 5*N, 1 + 4*N : 5*N) =  0.5_dp*(Q + T)
      ! HT(1 + 4*N : 5*N, 1 + 5*N : 6*N) =  0.0_dp
      HT(1 + 4*N : 5*N, 1 + 6*N : 7*N) =  IU * RQS3 * PZ
      HT(1 + 4*N : 5*N, 1 + 7*N : 8*N) =  IU * SQR2 * RQS3 * PP

      !col 6
      HT(1 + 5*N : 6*N, 1 + 0*N : 1*N) = -IU * SQR2 * R
      HT(1 + 5*N : 6*N, 1 + 1*N : 2*N) =  IU * SQR3 * RQS2 * (S)
      HT(1 + 5*N : 6*N, 1 + 2*N : 3*N) = -IU * RQS2 * (Q- T)
      HT(1 + 5*N : 6*N, 1 + 3*N : 4*N) = -IU * RQS2 * (SC)
      !HT(6,5) =  0.0_dp
      HT(1 + 5*N : 6*N, 1 + 5*N : 6*N) =  0.5_dp*(Q + T)
      HT(1 + 5*N : 6*N, 1 + 6*N : 7*N) =  SQR2 * RQS3 * PM
      HT(1 + 5*N : 6*N, 1 + 7*N : 8*N) = -RQS3 * PZ

      !col 7
      HT(1 + 6*N : 7*N, 1 + 0*N : 1*N) = -IU * PM
      HT(1 + 6*N : 7*N, 1 + 1*N : 2*N) =  SQR2 * RQS3 * (PZ)
      HT(1 + 6*N : 7*N, 1 + 2*N : 3*N) = -IU * RQS3 * PP
      !HT(7,4) =  0.0_dp
      HT(1 + 6*N : 7*N, 1 + 4*N : 5*N) = -IU * RQS3 * (PZ)
      HT(1 + 6*N : 7*N, 1 + 5*N : 6*N) =  SQR2 * RQS3 * PP
      HT(1 + 6*N : 7*N, 1 + 6*N : 7*N) =  A
      ! HT(7,8) =  0.0_dp

      !col 8
      !HT(8,1) =  0.0_dp
      HT(1 + 7*N : 8*N, 1 + 1*N : 2*N) = -RQS3 * PM
      HT(1 + 7*N : 8*N, 1 + 2*N : 3*N) = -IU * SQR2 * RQS3 * (PZ)
      HT(1 + 7*N : 8*N, 1 + 3*N : 4*N) = -PP
      HT(1 + 7*N : 8*N, 1 + 4*N : 5*N) = -IU * SQR2 * RQS3 * PM
      HT(1 + 7*N : 8*N, 1 + 5*N : 6*N) = -RQS3 * (PZ)
      !HT(8,7) =  0.0_dp
      HT(1 + 7*N : 8*N, 1 + 7*N : 8*N) =  A


      ! Profile
      forall(ii=1:N)
        HT(      ii,      ii) = HT(      ii,      ii) + profile(ii,1)
        HT(  N + ii,  N + ii) = HT(  N + ii,  N + ii) + profile(ii,1)
        HT(2*N + ii,2*N + ii) = HT(2*N + ii,2*N + ii) + profile(ii,1)
        HT(3*N + ii,3*N + ii) = HT(3*N + ii,3*N + ii) + profile(ii,1)
        HT(4*N + ii,4*N + ii) = HT(4*N + ii,4*N + ii) + profile(ii,2)
        HT(5*N + ii,5*N + ii) = HT(5*N + ii,5*N + ii) + profile(ii,2)
        HT(6*N + ii,6*N + ii) = HT(6*N + ii,6*N + ii) + profile(ii,3)
        HT(7*N + ii,7*N + ii) = HT(7*N + ii,7*N + ii) + profile(ii,3)
      end forall


      ! do i = 1, 8, 1
      !   write(102,*) (real(HT(i,j)), j=1,8)
      !   write(103,*) (aimag(HT(i,j)), j=1,8)
      ! end do
      ! stop

      ! HT = transpose(HT)

      deallocate(Q)
      deallocate(T)
      deallocate(S)
      deallocate(SC)
      deallocate(R)
      deallocate(RC)
      deallocate(PP)
      deallocate(PM)
      deallocate(PZ)
      deallocate(A)

      if (present(sparse)) then
        if (sparse .eqv. .True.) then
          nzmax = (N + (N-1)*2 + 2)*64
          ! print *, nzmax

          call dnscsr_z_mkl(nzmax, N*8, HT, HT_csr)

          ! allocate(Aa_a(nzmax))
          ! allocate(Aa_ia(nzmax))
          ! allocate(Aa_ja(nzmax))
          !
          ! next = 1
          ! do i = 1, N*8
          !   do j = 1, N*8
          !     call insertCOO_cmplx(Aa_a, Aa_ia, Aa_ja, HT(i,j), i, j, next, nzmax)
          !   end do
          ! end do
          ! ! print *, next-1
          !
          ! allocate(A_a(0))
          ! allocate(A_ia(0))
          ! allocate(A_ja(0))
          !
          ! A_a = [ Aa_a(1:next-1) ]
          ! A_ia = [ Aa_ia(1:next-1) ]
          ! A_ja = [ Aa_ja(1:next-1) ]
          !
          ! deallocate(Aa_a, Aa_ia, Aa_ja)
          !
          ! nzmax = ubound(A_a, 1)
          !
          ! info = mkl_sparse_z_create_coo (HT_coo, SPARSE_INDEX_BASE_ONE, N*8, &
          ! & N*8, ubound(A_a, dim=1), A_ia, A_ja, A_a)
          !
          ! if (info /= 0) stop 'error creating HT_coo'
          !
          ! info = mkl_sparse_convert_csr (HT_coo, SPARSE_OPERATION_NON_TRANSPOSE, HT_csr)
          !
          ! if (info /= 0) stop 'error converting HT_coo to HT_csr'
          !
          ! deallocate(A_a, A_ia, A_ja)
          ! deallocate (HT)

        end if
      end if



    end subroutine ZB8bandQW


    subroutine ZB8bandBulk(HT,wv,params,g)

      implicit none

      !input/output
      complex(kind=dp), intent(inout), dimension(:,:) :: HT
      type(wavevector), intent(in) :: wv
      type(paramStruct), intent(in) :: params(1)
      character(len=1), intent(in), optional :: g

      ! wave vector related
      real(kind=dp) :: kx2, ky2, kz2, kxky, k2
      complex(kind=dp) :: kminusz, kplusz, kminus, kplus

      !kp parameters
      real(kind=dp) :: gamma1, gamma2, gamma3, P, A
      complex(kind=dp) :: PM, PP, PZ
      complex(kind=dp) :: Q, R, RC, S, SC, T

      ! constants
      ! complex(kind=dp) :: IU
      ! real(kind=dp) :: SQR3, RQS2, SQR2, RQS3

      integer :: i, j


      !constants
      ! IU = dcmplx(0.0_dp,1.0_dp)
      ! SQR3 = dsqrt(3.0_dp)
      ! SQR2 = dsqrt(2.0_dp)
      ! RQS2 = 1.0_dp/dsqrt(2.0_dp)
      ! RQS3 = 1.0_dp/dsqrt(3.0_dp)


      !wave vectors
      kx2 = wv%kx**2
      ky2 = wv%ky**2
      kz2 = wv%kz**2
      kxky = wv%kx*wv%ky
      k2 = kx2 + ky2 + kz2

      kminusz = (wv%kx - IU*wv%ky)*wv%kz
      kplusz  = (wv%kx + IU*wv%ky)*wv%kz
      kplus   = wv%kx + IU*wv%ky
      kminus  = wv%kx - IU*wv%ky

      !kp parameters
      if (present(g) .and. g == 'g') then
        gamma1 = 0!params(1)%gamma1
        gamma2 = 0!params(1)%gamma2
        gamma3 = 0!params(1)%gamma3
        P = params(1)%P
        A = 0!params(1)%A
      else
        gamma1 = params(1)%gamma1
        gamma2 = params(1)%gamma2
        gamma3 = params(1)%gamma3
        P = params(1)%P
        A = params(1)%A
      end if


      ! kp terms
      Q = -( (gamma1 + gamma2)*(kx2 + ky2) + (gamma1 - 2.0_dp*gamma2)*kz2 )
      T = -( (gamma1 - gamma2)*(kx2 + ky2) + (gamma1 + 2.0_dp*gamma2)*kz2 )

      S  =   IU * 2.0_dp * SQR3 * gamma3 * kminusz
      SC = - IU * 2.0_dp * SQR3 * gamma3 * kplusz

      R  = - SQR3 * ( gamma2*(kx2 - ky2) - 2.0_dp*IU*gamma3*kxky )
      RC = - SQR3 * ( gamma2*(kx2 - ky2) + 2.0_dp*IU*gamma3*kxky )

      PP = P * kplus  * RQS2
      PM = P * kminus * RQS2
      PZ = P * wv%kz


      HT = 0

      ! kp matrix
      !col 1
      HT(1,1) =  Q
      HT(1,2) =  SC
      HT(1,3) =  RC
      HT(1,4) =  0.0_dp
      HT(1,5) = -IU * RQS2 * SC
      HT(1,6) =  IU * SQR2 * RC
      HT(1,7) =  IU * PP
      HT(1,8) =  0.0_dp

      !col 2
      HT(2,1) =  S
      HT(2,2) =  T
      HT(2,3) =  0.0_dp
      HT(2,4) =  RC
      HT(2,5) =  IU * RQS2 * (Q - T)
      HT(2,6) = -IU * SQR3 * RQS2 * SC
      HT(2,7) =  SQR2 * RQS3 * PZ
      HT(2,8) = -RQS3*PP

      !col 3
      HT(3,1) =  R
      !HT(3,2) =  0.0_dp
      HT(3,3) =  T
      HT(3,4) = -SC
      HT(3,5) =  IU * SQR3 * RQS2 * S
      HT(3,6) =  IU * RQS2 * (Q - T)
      HT(3,7) =  IU * RQS3 * PM
      HT(3,8) =  IU * SQR2 * RQS3 * PZ

      !col 4
      !HT(4,1) =  0.0_dp
      HT(4,2) =  R
      HT(4,3) = -S
      HT(4,4) =  Q
      HT(4,5) =  IU * SQR2 * R
      HT(4,6) =  IU * RQS2 * S
      HT(4,7) =  0.0_dp
      HT(4,8) = -PM

      !col 5
      HT(5,1) =  IU * RQS2 * S
      HT(5,2) = -IU * RQS2 * (Q - T)
      HT(5,3) = -IU * SQR3 * RQS2 * SC
      HT(5,4) = -IU * SQR2 * RC
      HT(5,5) =  0.5_dp*(Q + T)
      HT(5,6) =  0.0_dp
      HT(5,7) =  IU * RQS3 * PZ
      HT(5,8) =  IU * SQR2 * RQS3 * PP

      !col 6
      HT(6,1) = -IU * SQR2 * R
      HT(6,2) =  IU * SQR3 * RQS2 * S
      HT(6,3) = -IU * RQS2 * (Q - T)
      HT(6,4) = -IU * RQS2 * SC
      !HT(6,5) =  0.0_dp
      HT(6,6) =  0.5_dp*(Q + T)
      HT(6,7) =  SQR2 * RQS3 * PM
      HT(6,8) = -RQS3 * PZ

      !col 7
      HT(7,1) = -IU * PM
      HT(7,2) =  SQR2 * RQS3 * PZ
      HT(7,3) = -IU * RQS3 * PP
      !HT(7,4) =  0.0_dp
      HT(7,5) = -IU * RQS3 * PZ
      HT(7,6) =  SQR2 * RQS3 * PP
      HT(7,7) =  A * K2
      HT(7,8) =  0.0_dp

      !col 8
      !HT(8,1) =  0.0_dp
      HT(8,2) = -RQS3 * PM
      HT(8,3) = -IU * SQR2 * RQS3 * PZ
      HT(8,4) = -PP
      HT(8,5) = -IU * SQR2 * RQS3 * PM
      HT(8,6) = -RQS3 * PZ
      !HT(8,7) =  0.0_dp
      HT(8,8) =  A * K2


      ! SOC
      ! HT(1:6,1:6) = HT(1:6,1:6) + params(1)%EV
      HT(5,5) = HT(5,5) - params(1)%DeltaSO
      HT(6,6) = HT(6,6) - params(1)%DeltaSO

      HT(7,7) = HT(7,7) + params(1)%Eg
      HT(8,8) = HT(8,8) + params(1)%Eg
      ! HT(7,7) = HT(7,7) + params(1)%Ec
      ! HT(8,8) = HT(8,8) + params(1)%Ec

      ! do i = 1, 8, 1
      !   write(102,*) (real(HT(i,j)), j=1,8)
      !   write(103,*) (aimag(HT(i,j)), j=1,8)
      ! end do
      ! stop

    end subroutine ZB8bandBulk







end module hamiltonianConstructor
