module hamiltonianConstructor

  use definitions
  use finitedifferences
  use sparse_matrices
  use utils
  use strain_solver, only: compute_bp_scalar, bir_pikus_blocks_free
  use confinement_init
  use hamiltonian_wire

  implicit none

  private

  public :: externalFieldSetup_electricField
  public :: ZB8bandQW, ZB8bandBulk
  public :: add_bp_strain_dense

  contains

    subroutine externalFieldSetup_electricField(profile, Evalue, totalSize, z)

      real(kind = dp), intent(inout), allocatable, dimension(:,:) :: profile
      real(kind = dp), intent(in) :: Evalue, totalSize
      real(kind = dp), intent(in), dimension(:) :: z

      integer :: i

      do i = 1, ubound(z, dim=1), 1
          profile(i,:) = profile(i,:) - (Evalue*totalSize) * (z(i)+z(1))/(2.0_dp*z(1))
      end do

    end subroutine externalFieldSetup_electricField
    subroutine ZB8bandQW(HT, wv, profile, kpterms, cfg, sparse, HT_csr, g)

      implicit none

      !input/output
      complex(kind=dp), intent(inout), dimension(:,:) :: HT
      type(csr_matrix), intent(inout), optional :: HT_csr
      type(wavevector), intent(in) :: wv
      real(kind = dp), intent(in), dimension(:,:) :: profile
      real(kind = dp), intent(in), dimension(:,:,:) :: kpterms
      type(simulation_config), intent(in), optional :: cfg
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
        do jj = 1, N
          do ii = 1, N
            Q(ii,jj) = -( (kpterms(ii,jj,1) + kpterms(ii,jj,2))*k2 + kpterms(ii,jj,7) )
            T(ii,jj) = -( (kpterms(ii,jj,1) - kpterms(ii,jj,2))*k2 + kpterms(ii,jj,8) )
            S(ii,jj)  =  2.0_dp * SQR3 * kminus * kpterms(ii,jj,9) !no IU because of derivative
            SC(jj,ii) =  2.0_dp * SQR3 * kplus  * kpterms(ii,jj,9) !no IU because of derivative
            PZ(ii,jj) = kpterms(ii,jj,6) * (-IU)
            A(ii,jj)  = kpterms(ii,jj,5) + k2*kpterms(ii,jj,10)
          end do
        end do

        do ii = 1, N
          R(ii,ii)  = - SQR3 * ( kpterms(ii,ii,2)*(kx2 - ky2) - 2.0_dp*IU*kpterms(ii,ii,3)*kxky )
          RC(ii,ii) = - SQR3 * ( kpterms(ii,ii,2)*(kx2 - ky2) + 2.0_dp*IU*kpterms(ii,ii,3)*kxky )
          PP(ii,ii) = kpterms(ii,ii,4) * kplus  * RQS2
          PM(ii,ii) = kpterms(ii,ii,4) * kminus * RQS2
        end do
      else
        do jj = 1, N
          do ii = 1, N
            PZ(ii,jj) = kpterms(ii,jj,4) * wv%kz
          end do
        end do
        do ii = 1, N
          PP(ii,ii) = kpterms(ii,ii,4) * kplus  * RQS2
          PM(ii,ii) = kpterms(ii,ii,4) * kminus * RQS2
        end do

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
      do ii = 1, N
        HT(      ii,      ii) = HT(      ii,      ii) + profile(ii,1)
        HT(  N + ii,  N + ii) = HT(  N + ii,  N + ii) + profile(ii,1)
        HT(2*N + ii,2*N + ii) = HT(2*N + ii,2*N + ii) + profile(ii,1)
        HT(3*N + ii,3*N + ii) = HT(3*N + ii,3*N + ii) + profile(ii,1)
        HT(4*N + ii,4*N + ii) = HT(4*N + ii,4*N + ii) + profile(ii,2)
        HT(5*N + ii,5*N + ii) = HT(5*N + ii,5*N + ii) + profile(ii,2)
        HT(6*N + ii,6*N + ii) = HT(6*N + ii,6*N + ii) + profile(ii,3)
        HT(7*N + ii,7*N + ii) = HT(7*N + ii,7*N + ii) + profile(ii,3)
      end do

      ! Full Bir-Pikus strain (k-independent, not in g-mode)
      if (present(cfg)) then
        if (.not. present(g) .and. allocated(cfg%strain_blocks%delta_Ec)) then
          do ii = 1, N
            call add_bp_strain_dense(HT, ii, N, cfg%strain_blocks)
          end do
        end if
      end if


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
          ! 8x8 band blocks, each with (N + 2*(N-1) + 2) potential entries
          nzmax = (N + (N-1)*2 + 2) * 8 * 8

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
        gamma1 = 0.0_dp
        gamma2 = 0.0_dp
        gamma3 = 0.0_dp
        P = params(1)%P
        A = 0.0_dp
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
      HT(7,7) =  A * k2
      HT(7,8) =  0.0_dp

      !col 8
      HT(8,2) = -RQS3 * PM
      HT(8,3) = -IU * SQR2 * RQS3 * PZ
      HT(8,4) = -PP
      HT(8,5) = -IU * SQR2 * RQS3 * PM
      HT(8,6) = -RQS3 * PZ
      HT(8,8) =  A * k2


      ! SOC
      HT(5,5) = HT(5,5) - params(1)%DeltaSO
      HT(6,6) = HT(6,6) - params(1)%DeltaSO

      HT(7,7) = HT(7,7) + params(1)%Eg
      HT(8,8) = HT(8,8) + params(1)%Eg

      ! ---------------------------------------------------------------
      ! Full Bir-Pikus strain Hamiltonian for bulk.
      !
      ! When strainSubstrate > 0, apply uniform biaxial [001] strain:
      !   eps_xx = eps_yy = eps_par = (a_sub - a_film) / a_film
      !   eps_zz = eps_perp = -2 C12/C11 * eps_par
      !   eps_xy = eps_xz = eps_yz = 0  (for [001] biaxial)
      !
      ! The strain Hamiltonian has the SAME matrix structure as the
      ! k-dependent terms but with eps_ij replacing k_i*k_j and
      ! deformation potentials (av, b_dp, d_dp) replacing Luttinger
      ! parameters (gamma1, gamma2, gamma3):
      !
      !   P_eps     = -av * Tr(eps)                           (hydrostatic)
      !   Q_eps     =  b_dp/2 * (eps_zz - 0.5*(eps_xx+eps_yy)) (tetragonal shear)
      !   R_eps     = -sqrt(3) * [b_dp/2*(eps_xx-eps_yy) - i*d_dp*eps_xy]
      !   S_eps     =  i*2*sqrt(3) * d_dp * (eps_xz - i*eps_yz)
      !   S_eps_bar = -i*2*sqrt(3) * d_dp * (eps_xz + i*eps_yz)
      !
      ! Diagonal:
      !   CB:  +ac * Tr(eps)
      !   HH:  -P_eps + Q_eps
      !   LH:  -P_eps - Q_eps
      !   SO:  -P_eps
      !
      ! Off-diagonal: same pattern as k-terms R, S, S_bar.
      ! For [001] biaxial: R_eps = S_eps = 0, only Q_eps survives.
      ! All terms included for physics consistency.
      ! ---------------------------------------------------------------
      if (params(1)%strainSubstrate > 0.0_dp) then
        block
          real(kind=dp) :: a_film, eps_xx, eps_zz
          type(bir_pikus_blocks) :: bp_bulk

          a_film = params(1)%a0
          if (a_film > 0.0_dp) then
            eps_xx = (params(1)%strainSubstrate - a_film) / a_film
            eps_zz = -2.0_dp * params(1)%C12 / params(1)%C11 * eps_xx

            allocate(bp_bulk%delta_Ec(1), bp_bulk%delta_EHH(1), &
              bp_bulk%delta_ELH(1), bp_bulk%delta_ESO(1), &
              bp_bulk%R_eps(1), bp_bulk%S_eps(1), bp_bulk%QT2_eps(1))

            associate(s => compute_bp_scalar(params(1), eps_xx, eps_xx, eps_zz, &
                                              0.0_dp, 0.0_dp, 0.0_dp))
              bp_bulk%delta_Ec(1)  = s%delta_Ec
              bp_bulk%delta_EHH(1) = s%delta_EHH
              bp_bulk%delta_ELH(1) = s%delta_ELH
              bp_bulk%delta_ESO(1) = s%delta_ESO
              bp_bulk%R_eps(1)     = s%R_eps
              bp_bulk%S_eps(1)     = s%S_eps
              bp_bulk%QT2_eps(1)   = s%QT2_eps
            end associate

            call add_bp_strain_dense(HT, 1, 1, bp_bulk)
            call bir_pikus_blocks_free(bp_bulk)
          end if
        end block
      end if


    end subroutine ZB8bandBulk
    subroutine add_bp_strain_dense(HT, ii, N, bp)
      complex(kind=dp), intent(inout), contiguous :: HT(:,:)
      integer, intent(in) :: ii, N
      type(bir_pikus_blocks), intent(in) :: bp

      complex(kind=dp) :: R_eps_c, S_eps_c

      R_eps_c = conjg(bp%R_eps(ii))
      S_eps_c = conjg(bp%S_eps(ii))

      ! === Diagonal per-band ===
      HT(      ii,      ii) = HT(      ii,      ii) + bp%delta_EHH(ii)
      HT(  N + ii,  N + ii) = HT(  N + ii,  N + ii) + bp%delta_ELH(ii)
      HT(2*N + ii,2*N + ii) = HT(2*N + ii,2*N + ii) + bp%delta_ELH(ii)
      HT(3*N + ii,3*N + ii) = HT(3*N + ii,3*N + ii) + bp%delta_EHH(ii)
      HT(4*N + ii,4*N + ii) = HT(4*N + ii,4*N + ii) + bp%delta_ESO(ii)
      HT(5*N + ii,5*N + ii) = HT(5*N + ii,5*N + ii) + bp%delta_ESO(ii)
      HT(6*N + ii,6*N + ii) = HT(6*N + ii,6*N + ii) + bp%delta_Ec(ii)
      HT(7*N + ii,7*N + ii) = HT(7*N + ii,7*N + ii) + bp%delta_Ec(ii)

      ! === Off-diagonal: S_eps (HH-LH) ===
      HT(      ii,  N + ii) = HT(      ii,  N + ii) + S_eps_c
      HT(  N + ii,      ii) = HT(  N + ii,      ii) + bp%S_eps(ii)
      HT(2*N + ii,3*N + ii) = HT(2*N + ii,3*N + ii) - S_eps_c
      HT(3*N + ii,2*N + ii) = HT(3*N + ii,2*N + ii) - bp%S_eps(ii)

      ! === Off-diagonal: R_eps (HH-LH) ===
      HT(      ii,2*N + ii) = HT(      ii,2*N + ii) + R_eps_c
      HT(2*N + ii,      ii) = HT(2*N + ii,      ii) + bp%R_eps(ii)
      HT(  N + ii,3*N + ii) = HT(  N + ii,3*N + ii) + R_eps_c
      HT(3*N + ii,  N + ii) = HT(3*N + ii,  N + ii) + bp%R_eps(ii)

      ! === Off-diagonal: VB-SO coupling ===
      HT(      ii,4*N + ii) = HT(      ii,4*N + ii) - IU * RQS2 * S_eps_c
      HT(4*N + ii,      ii) = HT(4*N + ii,      ii) + IU * RQS2 * bp%S_eps(ii)
      HT(      ii,5*N + ii) = HT(      ii,5*N + ii) + IU * SQR2 * R_eps_c
      HT(5*N + ii,      ii) = HT(5*N + ii,      ii) - IU * SQR2 * bp%R_eps(ii)

      HT(  N + ii,4*N + ii) = HT(  N + ii,4*N + ii) + IU * RQS2 * bp%QT2_eps(ii)
      HT(4*N + ii,  N + ii) = HT(4*N + ii,  N + ii) - IU * RQS2 * bp%QT2_eps(ii)
      HT(  N + ii,5*N + ii) = HT(  N + ii,5*N + ii) - IU * SQR3o2 * S_eps_c
      HT(5*N + ii,  N + ii) = HT(5*N + ii,  N + ii) + IU * SQR3o2 * bp%S_eps(ii)

      HT(2*N + ii,4*N + ii) = HT(2*N + ii,4*N + ii) + IU * SQR3o2 * bp%S_eps(ii)
      HT(4*N + ii,2*N + ii) = HT(4*N + ii,2*N + ii) - IU * SQR3o2 * S_eps_c
      HT(2*N + ii,5*N + ii) = HT(2*N + ii,5*N + ii) + IU * RQS2 * bp%QT2_eps(ii)
      HT(5*N + ii,2*N + ii) = HT(5*N + ii,2*N + ii) - IU * RQS2 * bp%QT2_eps(ii)

      HT(3*N + ii,4*N + ii) = HT(3*N + ii,4*N + ii) - IU * SQR2 * bp%R_eps(ii)
      HT(4*N + ii,3*N + ii) = HT(4*N + ii,3*N + ii) + IU * SQR2 * R_eps_c
      HT(3*N + ii,5*N + ii) = HT(3*N + ii,5*N + ii) + IU * RQS2 * bp%S_eps(ii)
      HT(5*N + ii,3*N + ii) = HT(5*N + ii,3*N + ii) - IU * RQS2 * S_eps_c
    end subroutine add_bp_strain_dense

end module hamiltonianConstructor
