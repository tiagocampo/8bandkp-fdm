module hamiltonianConstructor

  use definitions
  use finitedifferences
  use utils
  use mkl_spblas

  implicit none

  contains

    subroutine build_kpterm_block(kpterms, profile_vec, central, forward, &
      & backward, diag, offup, offdown, N, term_idx, scale_factor, has_diag)

      real(kind=dp), intent(inout), dimension(:,:,:) :: kpterms
      real(kind=dp), intent(in), dimension(:) :: profile_vec
      real(kind=dp), intent(in), dimension(:,:) :: central, forward, backward
      real(kind=dp), intent(inout), dimension(:) :: diag, offup, offdown
      integer, intent(in) :: N, term_idx
      real(kind=dp), intent(in) :: scale_factor
      logical, intent(in) :: has_diag

      integer :: ii

      if (has_diag) then
        call dgemv('N', N, N, 1.0_dp, central, N, profile_vec, 1, 0.0_dp, diag, 1)
      end if
      call dgemv('N', N, N, 1.0_dp, forward, N, profile_vec, 1, 0.0_dp, offup, 1)
      call dgemv('N', N, N, 1.0_dp, backward, N, profile_vec, 1, 0.0_dp, offdown, 1)

      if (has_diag) then
        forall(ii=2:N-1)
          kpterms(ii,ii,term_idx) = diag(ii)
          kpterms(ii+1,ii,term_idx) = -offup(ii)
          kpterms(ii-1,ii,term_idx) = -offdown(ii)
        end forall
        kpterms(1,1,term_idx) = diag(1)
        kpterms(N,N,term_idx) = diag(N)
        kpterms(2,1,term_idx) = -offup(1)
        kpterms(N-1,N,term_idx) = -offdown(N)
      else
        forall(ii=2:N-1)
          kpterms(ii+1,ii,term_idx) = offup(ii)
          kpterms(ii-1,ii,term_idx) = -offdown(ii)
        end forall
        kpterms(2,1,term_idx) = offup(1)
        kpterms(N-1,N,term_idx) = -offdown(N)
      end if

      kpterms(:,:,term_idx) = kpterms(:,:,term_idx) * scale_factor

    end subroutine build_kpterm_block

    subroutine externalFieldSetup_electricField(profile, Evalue, totalSize, z)

      real(kind = dp), intent(inout), allocatable, dimension(:,:) :: profile
      real(kind = dp), intent(in) :: Evalue, totalSize
      real(kind = dp), intent(in), dimension(:) :: z

      integer :: i

      do i = 1, ubound(z, dim=1), 1
          profile(i,:) = profile(i,:) - (Evalue*totalSize) * (z(i)+z(1))/(2.0_dp*z(1))
      end do

    end subroutine externalFieldSetup_electricField

    subroutine confinementInitialization(z, startPos, endPos, material, nlayers,&
      & params, confDir, profile, kpterms, FDorder)

      real(kind = dp), intent(in), dimension(:) :: z
      integer, intent(in), dimension(:) :: startPos, endPos
      character(len = 255), intent(in) :: material(nlayers)
      integer, intent(in) :: nlayers
      type(paramStruct), intent(in) :: params(nlayers)
      character(len = 1), intent(in) :: confDir
      real(kind = dp), intent(inout), allocatable, dimension(:,:) :: profile
      real(kind = dp), intent(inout), dimension(:,:,:) :: kpterms
      integer, intent(in), optional :: FDorder

      real(kind=dp), allocatable, dimension(:,:) :: ScnDer, FstDer, kptermsProfile
      real(kind=dp), allocatable, dimension(:,:) :: forward, central, backward
      real(kind=dp), allocatable, dimension(:) :: diag, offup, offdown
      real(kind=dp), allocatable, dimension(:,:) :: D2, D1

      integer :: i, initIDX, endIDX, N, ii, jj
      integer :: order
      real(kind = dp) :: delta


      N = size(z, dim=1)
      delta = abs(z(2) - z(1))
      print *, delta

      ! Resolve FD order (default 2 for backward compatibility)
      if (present(FDorder)) then
        order = FDorder
      else
        order = 2
      end if

      if (allocated(profile)) deallocate(profile)
      allocate(profile(N,3))

      if (allocated(kptermsProfile)) deallocate(kptermsProfile)
      allocate(kptermsProfile(N,5))


      if (confDir == 'z') then

        do i = 1, nlayers, 1
          if (params(i)%EV == 0.0_dp .and. params(i)%EC == 0.0_dp) then
            print *, "ERROR: Material '", trim(material(i)), &
              & "' has EV=0 and EC=0. Band offsets are required."
            print *, "  Check that the material name is correct and EV/EC are in the database."
            stop 1
          end if
        end do

        do i = 1, nlayers, 1


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

      if (order == 2) then
        ! ---- Order 2: use existing tridiagonal approach (backward compatible) ----
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


        allocate(diag(N))
        allocate(offup(N))
        allocate(offdown(N))

        ! A*kz**2
        call build_kpterm_block(kpterms, kptermsProfile(1:N,4), central, forward, &
          & backward, diag, offup, offdown, N, 5, 1.0_dp/(2.0_dp*delta**2), .True.)


        !Q
        call build_kpterm_block(kpterms, kptermsProfile(1:N,1) - 2.0_dp*kptermsProfile(1:N,2), &
          & central, forward, backward, diag, offup, offdown, N, 7, 1.0_dp/(2.0_dp*delta**2), .True.)


        !T
        call build_kpterm_block(kpterms, kptermsProfile(1:N,1) + 2.0_dp*kptermsProfile(1:N,2), &
          & central, forward, backward, diag, offup, offdown, N, 8, 1.0_dp/(2.0_dp*delta**2), .True.)


        ! P*kz
        call build_kpterm_block(kpterms, kptermsProfile(1:N,5), central, forward, &
          & backward, diag, offup, offdown, N, 6, 1.0_dp/(4.0_dp*delta), .False.)


        ! S -> gamma3*kz
        call build_kpterm_block(kpterms, kptermsProfile(1:N,3), central, forward, &
          & backward, diag, offup, offdown, N, 9, 1.0_dp/(4.0_dp*delta), .False.)

      else
        ! ---- Higher order: use FD matrix approach ----
        ! Build 2nd-derivative and 1st-derivative FD matrices
        call buildFD2ndDerivMatrix(N, delta, order, D2)
        call buildFD1stDerivMatrix(N, delta, order, D1)

        ! A*kz**2 (term 5): profile = A(z), operator = d^2/dz^2
        call applyVariableCoeff(kpterms, kptermsProfile(1:N,4), D2, N, 5)

        ! Q (term 7): profile = gamma1 - 2*gamma2, operator = d^2/dz^2
        call applyVariableCoeff(kpterms, kptermsProfile(1:N,1) - 2.0_dp*kptermsProfile(1:N,2), &
          & D2, N, 7)

        ! T (term 8): profile = gamma1 + 2*gamma2, operator = d^2/dz^2
        call applyVariableCoeff(kpterms, kptermsProfile(1:N,1) + 2.0_dp*kptermsProfile(1:N,2), &
          & D2, N, 8)

        ! P*kz (term 6): profile = P(z), operator = d/dz
        call applyVariableCoeff(kpterms, kptermsProfile(1:N,5), D1, N, 6)

        ! S -> gamma3*kz (term 9): profile = gamma3(z), operator = d/dz
        call applyVariableCoeff(kpterms, kptermsProfile(1:N,3), D1, N, 9)

      end if

      if (allocated(kptermsProfile)) deallocate(kptermsProfile)
      if (allocated(forward)) deallocate(forward)
      if (allocated(backward)) deallocate(backward)
      if (allocated(central)) deallocate(central)
      if (allocated(diag)) deallocate(diag)
      if (allocated(offup)) deallocate(offup)
      if (allocated(offdown)) deallocate(offdown)
      if (allocated(D2)) deallocate(D2)
      if (allocated(D1)) deallocate(D1)

    end subroutine confinementInitialization

    !---------------------------------------------------------------------------
    !> Apply variable coefficient to FD matrix and store in kpterms.
    !> Computes: kpterms(:,:,term_idx) = -diag(profile) @ FD
    !> This gives the correct sign for both 1st and 2nd derivative terms
    !> in the k.p Hamiltonian.
    !---------------------------------------------------------------------------
    subroutine applyVariableCoeff(kpterms, profile_vec, FD, N, term_idx)

      real(kind=dp), intent(inout), dimension(:,:,:) :: kpterms
      real(kind=dp), intent(in), dimension(:) :: profile_vec
      real(kind=dp), intent(in), dimension(:,:) :: FD
      integer, intent(in) :: N, term_idx

      integer :: ii, jj

      ! kpterms(j,i) = -profile_vec(j) * FD(j,i)
      forall(ii=1:N, jj=1:N)
        kpterms(jj, ii, term_idx) = -profile_vec(jj) * FD(jj, ii)
      end forall

    end subroutine applyVariableCoeff

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

      integer :: i, j, N, ii, jj

      N = size(HT, dim=1)/8

      !constants


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
      HT(1 + 2*N : 3*N, 1 + 2*N : 3*N) =  T
      HT(1 + 2*N : 3*N, 1 + 3*N : 4*N) = -SC
      HT(1 + 2*N : 3*N, 1 + 4*N : 5*N) =  IU * SQR3 * RQS2 * S
      HT(1 + 2*N : 3*N, 1 + 5*N : 6*N) =  IU * RQS2 * (Q - T)
      HT(1 + 2*N : 3*N, 1 + 6*N : 7*N) =  IU * RQS3 * PM
      HT(1 + 2*N : 3*N, 1 + 7*N : 8*N) =  IU * SQR2 * RQS3 * PZ

      !col 4
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
      HT(1 + 5*N : 6*N, 1 + 5*N : 6*N) =  0.5_dp*(Q + T)
      HT(1 + 5*N : 6*N, 1 + 6*N : 7*N) =  SQR2 * RQS3 * PM
      HT(1 + 5*N : 6*N, 1 + 7*N : 8*N) = -RQS3 * PZ

      !col 7
      HT(1 + 6*N : 7*N, 1 + 0*N : 1*N) = -IU * PM
      HT(1 + 6*N : 7*N, 1 + 1*N : 2*N) =  SQR2 * RQS3 * (PZ)
      HT(1 + 6*N : 7*N, 1 + 2*N : 3*N) = -IU * RQS3 * PP
      HT(1 + 6*N : 7*N, 1 + 4*N : 5*N) = -IU * RQS3 * (PZ)
      HT(1 + 6*N : 7*N, 1 + 5*N : 6*N) =  SQR2 * RQS3 * PP
      HT(1 + 6*N : 7*N, 1 + 6*N : 7*N) =  A
      ! HT(7,8) =  0.0_dp

      !col 8
      HT(1 + 7*N : 8*N, 1 + 1*N : 2*N) = -RQS3 * PM
      HT(1 + 7*N : 8*N, 1 + 2*N : 3*N) = -IU * SQR2 * RQS3 * (PZ)
      HT(1 + 7*N : 8*N, 1 + 3*N : 4*N) = -PP
      HT(1 + 7*N : 8*N, 1 + 4*N : 5*N) = -IU * SQR2 * RQS3 * PM
      HT(1 + 7*N : 8*N, 1 + 5*N : 6*N) = -RQS3 * (PZ)
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

          call dnscsr_z_mkl(nzmax, N*8, HT, HT_csr)

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

      integer :: i, j


      !constants


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
      HT(3,3) =  T
      HT(3,4) = -SC
      HT(3,5) =  IU * SQR3 * RQS2 * S
      HT(3,6) =  IU * RQS2 * (Q - T)
      HT(3,7) =  IU * RQS3 * PM
      HT(3,8) =  IU * SQR2 * RQS3 * PZ

      !col 4
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
      HT(6,6) =  0.5_dp*(Q + T)
      HT(6,7) =  SQR2 * RQS3 * PM
      HT(6,8) = -RQS3 * PZ

      !col 7
      HT(7,1) = -IU * PM
      HT(7,2) =  SQR2 * RQS3 * PZ
      HT(7,3) = -IU * RQS3 * PP
      HT(7,5) = -IU * RQS3 * PZ
      HT(7,6) =  SQR2 * RQS3 * PP
      HT(7,7) =  A * K2
      HT(7,8) =  0.0_dp

      !col 8
      HT(8,2) = -RQS3 * PM
      HT(8,3) = -IU * SQR2 * RQS3 * PZ
      HT(8,4) = -PP
      HT(8,5) = -IU * SQR2 * RQS3 * PM
      HT(8,6) = -RQS3 * PZ
      HT(8,8) =  A * K2


      ! SOC
      HT(5,5) = HT(5,5) - params(1)%DeltaSO
      HT(6,6) = HT(6,6) - params(1)%DeltaSO

      HT(7,7) = HT(7,7) + params(1)%Eg
      HT(8,8) = HT(8,8) + params(1)%Eg


    end subroutine ZB8bandBulk


end module hamiltonianConstructor