module hamiltonianConstructor

  use definitions
  use finitedifferences
  use sparse_matrices
  use utils
  use strain_solver, only: compute_bp_scalar, bir_pikus_blocks_free, &
    get_strain_table, strain_entry, zeeman_entry, get_zeeman_table
  use confinement_init
  use hamiltonian_wire
  use magnetic_field, only: compute_gauge_shifts

  implicit none

  private

  public :: externalFieldSetup_electricField
  public :: ZB8bandQW, ZB8bandBulk, ZB8bandLandau

  contains

    subroutine externalFieldSetup_electricField(profile, Evalue, totalSize, z)

      real(kind = dp), intent(inout), allocatable, dimension(:,:) :: profile
      real(kind = dp), intent(in) :: Evalue, totalSize
      real(kind = dp), intent(in), dimension(:) :: z

      integer :: i

      if (size(z) > 0 .and. abs(z(1)) < tolerance) then
        print *, "ERROR: externalFieldSetup_electricField called with z(1)=0"
        print *, "  This causes division by zero in the linear potential."
        print *, "  Ensure z-coordinates start at a non-zero value."
        stop 1
      end if

      do i = 1, ubound(z, dim=1), 1
          profile(i,:) = profile(i,:) - (Evalue*totalSize) * (z(i)+z(1))/(2.0_dp*z(1))
      end do

    end subroutine externalFieldSetup_electricField
    subroutine ZB8bandQW(HT, wv, profile, kpterms, cfg, sparse, HT_csr, g)

      implicit none

      !input/output
      complex(kind=dp), intent(inout), contiguous :: HT(:,:)
      type(csr_matrix), intent(inout), optional :: HT_csr
      type(wavevector), intent(in) :: wv
      real(kind = dp), intent(in), contiguous :: profile(:,:)
      real(kind = dp), intent(in), contiguous :: kpterms(:,:,:)
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

      ! Full Bir-Pikus strain (k-independent, not in g-mode).
      ! Table-driven via get_strain_table() from strain_solver.
      if (present(cfg)) then
        if (.not. present(g) .and. allocated(cfg%strain_blocks%delta_Ec)) then
          call apply_strain_table_dense(HT, N, cfg%strain_blocks)
        end if
      end if

      ! ---------------------------------------------------------------
      ! Zeeman splitting: g*mu_B * |B| * sigma
      ! Applied after profile and strain so all terms are present.
      ! Uses full B-field magnitude (not just B_z). This is physically
      ! correct: spin splitting depends on total |B|, regardless of
      ! direction. Orbital effects (Landau levels) from B_perp are
      ! handled separately via Peierls substitution in wire mode.
      ! Uses the data-driven Zeeman table from strain_solver.
      ! ---------------------------------------------------------------
      if (present(cfg)) then
        if (.not. present(g) .and. cfg%bdg%enabled) then
          block
            real(kind=dp) :: B_mag, E0
            type(zeeman_entry) :: ztable(8)
            integer :: b
            B_mag = sqrt(sum(cfg%bdg%B_vec**2))
            E0 = cfg%bdg%g_factor * mu_B * B_mag
            ztable = get_zeeman_table()
            do ii = 1, N
              do b = 1, 8
                HT(ztable(b)%band_index*N + ii, &
                   ztable(b)%band_index*N + ii) = &
                   HT(ztable(b)%band_index*N + ii, &
                      ztable(b)%band_index*N + ii) &
                   + ztable(b)%g_multiplier * E0
              end do
            end do
          end block
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

    subroutine ZB8bandLandau(HT, wv, profile, kpterms, x_grid, cfg)
      ! Build the 8N x 8N Landau Hamiltonian on a 1D x-grid.
      !
      ! Landau gauge: A = (0, Bz*x, -By*x)
      !   Pi_y(i) = ky + Bz * x(i) * 1e-20 / hbar   [1/AA]
      !   Pi_z(i) = kz - By * x(i) * 1e-20 / hbar   [1/AA]
      !
      ! B_x contributes only Zeeman splitting (parallel to discretized direction).
      ! The block structure is identical to ZB8bandQW.
      implicit none

      complex(kind=dp), intent(inout), contiguous :: HT(:,:)
      type(wavevector), intent(in) :: wv
      real(kind=dp), intent(in), contiguous :: profile(:,:)
      real(kind=dp), intent(in), contiguous :: kpterms(:,:,:)
      real(kind=dp), intent(in), contiguous :: x_grid(:)
      type(simulation_config), intent(in), optional :: cfg

      integer :: N, ii, jj
      real(kind=dp) :: By, Bz, B_mag
      real(kind=dp), allocatable :: Pi_y(:), Pi_z(:)
      real(kind=dp), allocatable :: kx2_arr(:), ky2_arr(:), k2_arr(:), kxky_arr(:)
      complex(kind=dp), allocatable :: kplus_arr(:), kminus_arr(:)

      complex(kind=dp), allocatable :: Q(:,:), R(:,:), RC(:,:), S(:,:), SC(:,:)
      complex(kind=dp), allocatable :: T(:,:), PZ(:,:), PP(:,:), PM(:,:), Amat(:,:)

      N = size(HT, dim=1) / 8

      ! Default B-field values
      By = 0.0_dp
      Bz = 0.0_dp
      B_mag = 0.0_dp
      if (present(cfg)) then
        By = cfg%bdg%B_vec(2)
        Bz = cfg%bdg%B_vec(3)
        B_mag = sqrt(sum(cfg%bdg%B_vec**2))
      end if

      ! Compute position-dependent gauge-shifted k-values
      allocate(Pi_y(N), Pi_z(N))
      call compute_gauge_shifts(x_grid, [0.0_dp, By, Bz], wv%ky, wv%kz, Pi_y, Pi_z)

      ! Compute position-dependent k-values
      allocate(kx2_arr(N), ky2_arr(N), k2_arr(N), kxky_arr(N))
      allocate(kplus_arr(N), kminus_arr(N))

      do ii = 1, N
        kx2_arr(ii) = Pi_y(ii)**2
        ky2_arr(ii) = Pi_z(ii)**2
        k2_arr(ii) = kx2_arr(ii) + ky2_arr(ii)
        kxky_arr(ii) = Pi_y(ii) * Pi_z(ii)
        kplus_arr(ii) = Pi_y(ii) + IU * Pi_z(ii)
        kminus_arr(ii) = Pi_y(ii) - IU * Pi_z(ii)
      end do

      ! Compute midpoint values for off-diagonal FD hopping
      ! (used conceptually: S/SC hermiticity guaranteed via SC = conjg(S^T))

      ! Allocate kp matrices
      allocate(Q(N,N), T(N,N), S(N,N), SC(N,N))
      allocate(R(N,N), RC(N,N), PP(N,N), PM(N,N), PZ(N,N), Amat(N,N))
      Q = ZERO; T = ZERO; S = ZERO; SC = ZERO
      R = ZERO; RC = ZERO; PP = ZERO; PM = ZERO; PZ = ZERO; Amat = ZERO

      ! Build Q, T, A matrices using kpterms FD stencils
      ! Q = -(kpterms1 + kpterms2)*k2 + kpterms7) with position-dependent k2
      ! But kpterms already encodes the FD stencil for the x-direction kinetic
      ! energy. Here k2 = kx2 + ky2 is the in-plane (Pi_y, Pi_z) contribution.
      ! We need to add kpterms (x-derivative part) + in-plane k2 terms.
      !
      ! In ZB8bandQW: Q = -(kpterms(1)+kpterms(2))*k2 + kpterms(7))
      !   kpterms(1) = gamma1 (diagonal), kpterms(2) = gamma2 (diagonal)
      !   kpterms(7) = FD stencil for (gamma1 - 2*gamma2)*d2/dx2
      ! Here we separate: kpterms(7) is the x-FD part, and we add position-dependent
      ! in-plane k2 contribution to the diagonal.
      !
      ! Actually, looking at confinementInitialization_landau:
      !   kpterms(ii,ii,1) = gamma1  (diagonal material param)
      !   kpterms(ii,ii,2) = gamma2  (diagonal material param)
      !   kpterms(ii,ii,3) = gamma3  (diagonal material param)
      !   kpterms(ii,ii,4) = P       (diagonal material param)
      !   kpterms(:,:,5)   = A * FD stencil (x-kinetic energy for CB)
      !   kpterms(:,:,6)   = P * FD stencil (x-momentum)
      !   kpterms(:,:,7)   = (gamma1-2gamma2) * FD stencil (x-Q term)
      !   kpterms(:,:,8)   = (gamma1+2gamma2) * FD stencil (x-T term)
      !   kpterms(:,:,9)   = gamma3 * FD stencil (x-S term)
      !   kpterms(ii,ii,10)= A       (diagonal material param)
      !
      ! In ZB8bandQW (with uniform k):
      !   Q(ii,jj) = -((kpterms1+kpterms2)*k2 + kpterms7)
      ! Since kpterms1,2 are diagonal, and kpterms7 has FD off-diagonal,
      ! this is: Q = -(gamma1+gamma2)*k2 * I - FD_x_term
      !
      ! For Landau, k2 is position-dependent, so we need per-point scaling.

      ! Diagonal elements: use position-dependent k2 at each grid point
      do jj = 1, N
        do ii = 1, N
          ! Q = -(kpterms1 + kpterms2) * k2(ii) + kpterms7
          ! For diagonal terms, kpterms1,2 are nonzero only on diagonal
          Q(ii,jj) = -((kpterms(ii,jj,1) + kpterms(ii,jj,2)) * k2_arr(jj) &
            & + kpterms(ii,jj,7))
          ! T = -(kpterms1 - kpterms2) * k2(ii) + kpterms8
          T(ii,jj) = -((kpterms(ii,jj,1) - kpterms(ii,jj,2)) * k2_arr(jj) &
            & + kpterms(ii,jj,8))
          ! A = kpterms5 + kpterms10 * k2(ii)
          Amat(ii,jj) = kpterms(ii,jj,5) + kpterms(ii,jj,10) * k2_arr(jj)
          ! PZ = kpterms6 * (-IU)  [P * d/dx FD stencil]
          PZ(ii,jj) = kpterms(ii,jj,6) * (-IU)
        end do
      end do

      ! S and SC: position-dependent kminus/kplus with FD stencil
      ! S = 2 * sqrt(3) * kminus(i) * kpterms(i,j,9)
      ! In QW (uniform k): SC(j,i) = 2*SQR3 * kplus * kpterms(i,j,9) = conjg(S(i,j))
      ! For Landau: use per-element kminus for S, then SC = conjg(transpose(S))
      ! to guarantee hermiticity by construction.
      do jj = 1, N
        do ii = 1, N
          if (abs(kpterms(ii,jj,9)) > 0.0_dp) then
            S(ii,jj) = 2.0_dp * SQR3 * kminus_arr(ii) * kpterms(ii,jj,9)
          end if
        end do
      end do
      ! SC = S^H (guarantees hermiticity)
      do jj = 1, N
        do ii = 1, N
          SC(ii,jj) = conjg(S(jj,ii))
        end do
      end do

      ! R and RC: purely diagonal (like ZB8bandQW)
      ! R = -sqrt(3) * (gamma2*(kx2-ky2) - 2*i*gamma3*kxky)
      ! With position-dependent k-values:
      do ii = 1, N
        R(ii,ii) = -SQR3 * (kpterms(ii,ii,2) * (kx2_arr(ii) - ky2_arr(ii)) &
          & - 2.0_dp * IU * kpterms(ii,ii,3) * kxky_arr(ii))
        RC(ii,ii) = -SQR3 * (kpterms(ii,ii,2) * (kx2_arr(ii) - ky2_arr(ii)) &
          & + 2.0_dp * IU * kpterms(ii,ii,3) * kxky_arr(ii))
        ! PP and PM: position-dependent
        PP(ii,ii) = kpterms(ii,ii,4) * kplus_arr(ii) * RQS2
        PM(ii,ii) = kpterms(ii,ii,4) * kminus_arr(ii) * RQS2
      end do

      ! Initialize HT
      HT = ZERO

      ! Assemble 8Nx8N Hamiltonian (identical block structure to ZB8bandQW)
      ! col 1
      HT(1 + 0*N : 1*N, 1 + 0*N : 1*N) =  Q
      HT(1 + 0*N : 1*N, 1 + 1*N : 2*N) =  SC
      HT(1 + 0*N : 1*N, 1 + 2*N : 3*N) =  RC
      HT(1 + 0*N : 1*N, 1 + 4*N : 5*N) = -IU * RQS2 * SC
      HT(1 + 0*N : 1*N, 1 + 5*N : 6*N) =  IU * SQR2 * RC
      HT(1 + 0*N : 1*N, 1 + 6*N : 7*N) =  IU * PP

      ! col 2
      HT(1 + 1*N : 2*N, 1 + 0*N : 1*N) =  S
      HT(1 + 1*N : 2*N, 1 + 1*N : 2*N) =  T
      HT(1 + 1*N : 2*N, 1 + 3*N : 4*N) =  RC
      HT(1 + 1*N : 2*N, 1 + 4*N : 5*N) =  IU * RQS2 * (Q - T)
      HT(1 + 1*N : 2*N, 1 + 5*N : 6*N) = -IU * SQR3 * RQS2 * SC
      HT(1 + 1*N : 2*N, 1 + 6*N : 7*N) =  SQR2 * RQS3 * PZ
      HT(1 + 1*N : 2*N, 1 + 7*N : 8*N) = -RQS3 * PP

      ! col 3
      HT(1 + 2*N : 3*N, 1 + 0*N : 1*N) =  R
      HT(1 + 2*N : 3*N, 1 + 2*N : 3*N) =  T
      HT(1 + 2*N : 3*N, 1 + 3*N : 4*N) = -SC
      HT(1 + 2*N : 3*N, 1 + 4*N : 5*N) =  IU * SQR3 * RQS2 * S
      HT(1 + 2*N : 3*N, 1 + 5*N : 6*N) =  IU * RQS2 * (Q - T)
      HT(1 + 2*N : 3*N, 1 + 6*N : 7*N) =  IU * RQS3 * PM
      HT(1 + 2*N : 3*N, 1 + 7*N : 8*N) =  IU * SQR2 * RQS3 * PZ

      ! col 4
      HT(1 + 3*N : 4*N, 1 + 1*N : 2*N) =  R
      HT(1 + 3*N : 4*N, 1 + 2*N : 3*N) = -S
      HT(1 + 3*N : 4*N, 1 + 3*N : 4*N) =  Q
      HT(1 + 3*N : 4*N, 1 + 4*N : 5*N) =  IU * SQR2 * R
      HT(1 + 3*N : 4*N, 1 + 5*N : 6*N) =  IU * RQS2 * S
      HT(1 + 3*N : 4*N, 1 + 7*N : 8*N) = -PM

      ! col 5
      HT(1 + 4*N : 5*N, 1 + 0*N : 1*N) =  IU * RQS2 * S
      HT(1 + 4*N : 5*N, 1 + 1*N : 2*N) = -IU * RQS2 * (Q - T)
      HT(1 + 4*N : 5*N, 1 + 2*N : 3*N) = -IU * SQR3 * RQS2 * SC
      HT(1 + 4*N : 5*N, 1 + 3*N : 4*N) = -IU * SQR2 * RC
      HT(1 + 4*N : 5*N, 1 + 4*N : 5*N) =  0.5_dp * (Q + T)
      HT(1 + 4*N : 5*N, 1 + 6*N : 7*N) =  IU * RQS3 * PZ
      HT(1 + 4*N : 5*N, 1 + 7*N : 8*N) =  IU * SQR2 * RQS3 * PP

      ! col 6
      HT(1 + 5*N : 6*N, 1 + 0*N : 1*N) = -IU * SQR2 * R
      HT(1 + 5*N : 6*N, 1 + 1*N : 2*N) =  IU * SQR3 * RQS2 * S
      HT(1 + 5*N : 6*N, 1 + 2*N : 3*N) = -IU * RQS2 * (Q - T)
      HT(1 + 5*N : 6*N, 1 + 3*N : 4*N) = -IU * RQS2 * SC
      HT(1 + 5*N : 6*N, 1 + 5*N : 6*N) =  0.5_dp * (Q + T)
      HT(1 + 5*N : 6*N, 1 + 6*N : 7*N) =  SQR2 * RQS3 * PM
      HT(1 + 5*N : 6*N, 1 + 7*N : 8*N) = -RQS3 * PZ

      ! col 7
      HT(1 + 6*N : 7*N, 1 + 0*N : 1*N) = -IU * PM
      HT(1 + 6*N : 7*N, 1 + 1*N : 2*N) =  SQR2 * RQS3 * PZ
      HT(1 + 6*N : 7*N, 1 + 2*N : 3*N) = -IU * RQS3 * PP
      HT(1 + 6*N : 7*N, 1 + 4*N : 5*N) = -IU * RQS3 * PZ
      HT(1 + 6*N : 7*N, 1 + 5*N : 6*N) =  SQR2 * RQS3 * PP
      HT(1 + 6*N : 7*N, 1 + 6*N : 7*N) =  Amat

      ! col 8
      HT(1 + 7*N : 8*N, 1 + 1*N : 2*N) = -RQS3 * PM
      HT(1 + 7*N : 8*N, 1 + 2*N : 3*N) = -IU * SQR2 * RQS3 * PZ
      HT(1 + 7*N : 8*N, 1 + 3*N : 4*N) = -PP
      HT(1 + 7*N : 8*N, 1 + 4*N : 5*N) = -IU * SQR2 * RQS3 * PM
      HT(1 + 7*N : 8*N, 1 + 5*N : 6*N) = -RQS3 * PZ
      HT(1 + 7*N : 8*N, 1 + 7*N : 8*N) =  Amat

      ! Add profile (band edges)
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

      ! Zeeman splitting (full |B| magnitude)
      ! Uses the data-driven Zeeman table from strain_solver.
      if (present(cfg)) then
        if (B_mag > 1.0e-12_dp) then
          block
            real(kind=dp) :: E0
            type(zeeman_entry) :: ztable(8)
            integer :: b
            E0 = cfg%bdg%g_factor * mu_B * B_mag
            ztable = get_zeeman_table()
            do ii = 1, N
              do b = 1, 8
                HT(ztable(b)%band_index*N + ii, &
                   ztable(b)%band_index*N + ii) = &
                   HT(ztable(b)%band_index*N + ii, &
                      ztable(b)%band_index*N + ii) &
                   + ztable(b)%g_multiplier * E0
              end do
            end do
          end block
        end if
      end if

      ! Clean up
      deallocate(Q, T, S, SC, R, RC, PP, PM, PZ, Amat)
      deallocate(Pi_y, Pi_z, kx2_arr, ky2_arr, k2_arr, kxky_arr)
      deallocate(kplus_arr, kminus_arr)

    end subroutine ZB8bandLandau

    subroutine ZB8bandBulk(HT,wv,params,cfg,g)

      implicit none

      !input/output
      complex(kind=dp), intent(inout), dimension(:,:) :: HT
      type(wavevector), intent(in) :: wv
      type(paramStruct), intent(in) :: params(1)
      type(simulation_config), intent(in), optional :: cfg
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
      if (present(g)) then
        if (g == 'g') then
          gamma1 = 0.0_dp
          gamma2 = 0.0_dp
          gamma3 = 0.0_dp
          P = params(1)%P
          A = 0.0_dp
        else
          gamma1 = params(1)%gamma1 * const
          gamma2 = params(1)%gamma2 * const
          gamma3 = params(1)%gamma3 * const
          P = params(1)%P
          A = params(1)%A * const
        end if
      else
        ! Scale gamma and A by const = hbar^2/(2*m_0) so that k-dependent
        ! diagonal terms have energy units (eV). P already includes const
        ! via P = sqrt(EP*const), so it is NOT scaled here.
        gamma1 = params(1)%gamma1 * const
        gamma2 = params(1)%gamma2 * const
        gamma3 = params(1)%gamma3 * const
        P = params(1)%P
        A = params(1)%A * const
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

      ! External potential shifts (EF + SC)
      if (present(cfg)) then
        if (cfg%ExternalField == 1 .and. cfg%EFtype == 'EF') then
          do i = 1, 8
            HT(i,i) = HT(i,i) + cmplx(cfg%Evalue, 0.0_dp, kind=dp)
          end do
        end if
        if (cfg%sc_potential_shift /= 0.0_dp) then
          do i = 1, 8
            HT(i,i) = HT(i,i) + cmplx(cfg%sc_potential_shift, 0.0_dp, kind=dp)
          end do
        end if
      end if

      ! ---------------------------------------------------------------
      ! Zeeman splitting: g*mu_B * B . sigma
      ! Applied after SOC and band-edge shifts so all terms are present.
      ! Note: This is SPIN Zeeman only. Orbital Landau level quantization
      ! requires Peierls substitution (kz -> kz - eBx*y/hbar) which needs
      ! spatial y-discretization. Bulk 8x8 at k=0 has no y information
      ! so only Zeeman (spin) is possible here. For orbital physics, use
      ! wire mode (confinement=2) which has 2D grid and can do Peierls.
      ! Uses the data-driven Zeeman table from strain_solver.
      ! ---------------------------------------------------------------
      if (present(cfg)) then
        if (cfg%bdg%enabled) then
          block
            real(kind=dp) :: B_mag, E0
            type(zeeman_entry) :: ztable(8)
            integer :: b
            B_mag = sqrt(sum(cfg%bdg%B_vec**2))
            E0 = cfg%bdg%g_factor * mu_B * B_mag
            ztable = get_zeeman_table()
            do b = 1, 8
              HT(ztable(b)%band_index + 1, ztable(b)%band_index + 1) = &
                HT(ztable(b)%band_index + 1, ztable(b)%band_index + 1) &
                + ztable(b)%g_multiplier * E0
            end do
          end block
        end if
      end if

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
      !   Q_eps     = -(b_dp/2) * (eps_zz - 0.5*(eps_xx+eps_yy)) (tetragonal shear)
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

            call apply_strain_table_dense(HT, 1, bp_bulk)
            call bir_pikus_blocks_free(bp_bulk)
          end if
        end block
      end if


    end subroutine ZB8bandBulk

    ! ==================================================================
    ! Table-driven Bir-Pikus strain insertion for dense Hamiltonian.
    !
    ! Reads the 32-entry strain table from get_strain_table() and applies
    ! all strain terms to the dense 8N x 8N Hamiltonian.  Uses scalar
    ! field lookups from bir_pikus_blocks (no array temporaries) to
    ! avoid -O3 + OpenMP stack corruption.
    !
    ! Works for both QW (N > 1) and bulk (N = 1).
    ! ==================================================================
    subroutine apply_strain_table_dense(HT, N, bp)
      complex(kind=dp), intent(inout), contiguous :: HT(:,:)
      integer, intent(in) :: N
      type(bir_pikus_blocks), intent(in) :: bp

      type(strain_entry) :: table(32)
      integer :: ii, e, row, col
      complex(kind=dp) :: field_val

      table = get_strain_table()

      do ii = 1, N
        do e = 1, 32
          ! Scalar field lookup from bir_pikus_blocks — no array temporaries
          select case (table(e)%field_id)
          case (1); field_val = cmplx(bp%delta_EHH(ii), 0.0_dp, kind=dp)
          case (2); field_val = cmplx(bp%delta_ELH(ii), 0.0_dp, kind=dp)
          case (3); field_val = cmplx(bp%delta_ESO(ii), 0.0_dp, kind=dp)
          case (4); field_val = cmplx(bp%delta_Ec(ii), 0.0_dp, kind=dp)
          case (5); field_val = bp%S_eps(ii)
          case (6); field_val = bp%R_eps(ii)
          case (7); field_val = cmplx(bp%QT2_eps(ii), 0.0_dp, kind=dp)
          end select

          if (table(e)%use_conjg) field_val = conjg(field_val)

          row = table(e)%row_band * N + ii
          col = table(e)%col_band * N + ii

          HT(row, col) = HT(row, col) + table(e)%prefactor * field_val
        end do
      end do
    end subroutine apply_strain_table_dense

end module hamiltonianConstructor
