module hamiltonian_blocks

  ! ==============================================================================
  ! 8-band zinc-blende k.p Hamiltonian block structure as DATA.
  !
  ! Defines a table of 52 entries mapping (band_pair, kp_term, prefactor).
  ! Each entry describes one block in the 8x8 band space: which two bands
  ! couple, which k.p term (Q, R, S, T, P, etc.) is involved, and what
  ! complex prefactor multiplies it.
  !
  ! This is the foundation for Phase C where both dense and COO builders
  ! import this table instead of hard-coding the block structure.
  ! ==============================================================================

  use definitions, only: dp, IU, SQR2, SQR3, RQS2, RQS3, SQR3o2
  implicit none

  private

  ! ------------------------------------------------------------------
  ! Named constants for k.p terms.
  ! Each constant identifies one of the matrix blocks computed from
  ! material parameters and wavevector components:
  !   Q, T    : diagonal HH/LH kinetic energy blocks
  !   S, SC   : off-diagonal linear-k shear coupling
  !   R, RC   : off-diagonal quadratic-k coupling (conjugate pair)
  !   PP, PM  : Kane momentum * kplus/kminus
  !   PZ      : Kane momentum * kz
  !   A       : conduction band kinetic energy
  !   KP_DIFF     : Q - T  (LH-T mixing combination)
  !   KP_HALF_SUM : 0.5*(Q + T)  (SO diagonal combination)
  ! ------------------------------------------------------------------
  integer, parameter, public :: KP_Q         = 1
  integer, parameter, public :: KP_T         = 2
  integer, parameter, public :: KP_S         = 3
  integer, parameter, public :: KP_SC        = 4
  integer, parameter, public :: KP_R         = 5
  integer, parameter, public :: KP_RC        = 6
  integer, parameter, public :: KP_PP        = 7
  integer, parameter, public :: KP_PM        = 8
  integer, parameter, public :: KP_PZ        = 9
  integer, parameter, public :: KP_A         = 10
  integer, parameter, public :: KP_DIFF      = 11   ! Q - T
  integer, parameter, public :: KP_HALF_SUM  = 12   ! 0.5*(Q + T)

  ! ------------------------------------------------------------------
  ! Table entry type: one block in the 8x8 band Hamiltonian.
  !
  ! row_band, col_band : 0-based band indices (0-7)
  !   Bands 0,3 = HH, 1,2 = LH, 4,5 = SO, 6,7 = CB
  ! kp_term : which k.p matrix block to use (KP_Q, KP_R, ...)
  ! prefactor : complex scalar multiplier
  ! use_conjg : if true, apply conjg() to the k.p matrix block before
  !   multiplying by prefactor (reserved for future use; currently .false.
  !   for all k.p entries since conjugation is encoded in the prefactor)
  ! ------------------------------------------------------------------
  type, public :: kp_entry
    integer          :: row_band = 0    ! 0-based band index (0-7)
    integer          :: col_band = 0    ! 0-based band index (0-7)
    integer          :: kp_term  = 0    ! KP_Q, KP_R, etc.
    complex(kind=dp) :: prefactor = cmplx(0.0_dp, 0.0_dp, kind=dp)
    logical          :: use_conjg = .false.
  end type kp_entry

  public :: get_kp_block_table

  ! Module-level cache
  logical, save :: kp_table_cached = .false.
  type(kp_entry), save :: kp_table_cache(52)

contains

  ! ==================================================================
  ! Return the 52-entry k.p block table (cached).
  !
  ! Each entry is extracted from the HT() assignments in ZB8bandQW
  ! (hamiltonianConstructor.f90).  The 8x8 block structure is:
  !
  !   Band layout (0-based):
  !     0,3 = HH (heavy hole)
  !     1,2 = LH (light hole)
  !     4,5 = SO (split-off)
  !     6,7 = CB (conduction band)
  !
  !   Matrix blocks:
  !     Q  = -((gamma1+gamma2)*k2 + kpterms7)     (HH diagonal)
  !     T  = -((gamma1-gamma2)*k2 + kpterms8)     (LH diagonal)
  !     S  = 2*sqrt(3)*kminus*kpterms9            (shear)
  !     SC = 2*sqrt(3)*kplus*kpterms9             (shear conjugate)
  !     R  = -sqrt(3)*(gamma2*(kx2-ky2) - 2i*gamma3*kxky)
  !     RC = -sqrt(3)*(gamma2*(kx2-ky2) + 2i*gamma3*kxky)
  !     PP = P * kplus * RQS2
  !     PM = P * kminus * RQS2
  !     PZ = P * kz * (-IU)   (or P * kz in g-mode)
  !     A  = kpterms5 + k2*kpterms10              (CB diagonal)
  !     DIFF = Q - T
  !     HALF_SUM = 0.5*(Q + T)
  ! ==================================================================
  function get_kp_block_table() result(table)
    type(kp_entry) :: table(52)

    if (.not. kp_table_cached) then
      kp_table_cache = build_kp_table()
      kp_table_cached = .true.
    end if
    table = kp_table_cache
  end function get_kp_block_table

  function build_kp_table() result(table)
    type(kp_entry) :: table(52)

    ! ==================================================================
    ! Row 0 (HH1): 6 entries
    ! HT(0,0) = Q
    ! HT(0,1) = SC
    ! HT(0,2) = RC
    ! HT(0,4) = -IU*RQS2*SC
    ! HT(0,5) = IU*SQR2*RC
    ! HT(0,6) = IU*PP
    ! ==================================================================
    table( 1) = kp_entry(0, 0, KP_Q,  cmplx(1.0_dp, 0.0_dp, kind=dp))
    table( 2) = kp_entry(0, 1, KP_SC, cmplx(1.0_dp, 0.0_dp, kind=dp))
    table( 3) = kp_entry(0, 2, KP_RC, cmplx(1.0_dp, 0.0_dp, kind=dp))
    table( 4) = kp_entry(0, 4, KP_SC, -IU*RQS2)
    table( 5) = kp_entry(0, 5, KP_RC, IU*SQR2)
    table( 6) = kp_entry(0, 6, KP_PP, IU)

    ! ==================================================================
    ! Row 1 (HH2/LH1): 7 entries
    ! HT(1,0) = S
    ! HT(1,1) = T
    ! HT(1,3) = RC
    ! HT(1,4) = IU*RQS2*(Q-T)
    ! HT(1,5) = -IU*SQR3*RQS2*SC
    ! HT(1,6) = SQR2*RQS3*PZ
    ! HT(1,7) = -RQS3*PP
    ! ==================================================================
    table( 7) = kp_entry(1, 0, KP_S,  cmplx(1.0_dp, 0.0_dp, kind=dp))
    table( 8) = kp_entry(1, 1, KP_T,  cmplx(1.0_dp, 0.0_dp, kind=dp))
    table( 9) = kp_entry(1, 3, KP_RC, cmplx(1.0_dp, 0.0_dp, kind=dp))
    table(10) = kp_entry(1, 4, KP_DIFF, IU*RQS2)
    table(11) = kp_entry(1, 5, KP_SC, -IU*SQR3*RQS2)
    table(12) = kp_entry(1, 6, KP_PZ, cmplx(SQR2*RQS3, 0.0_dp, kind=dp))
    table(13) = kp_entry(1, 7, KP_PP, cmplx(-RQS3, 0.0_dp, kind=dp))

    ! ==================================================================
    ! Row 2 (LH1/LH2): 7 entries
    ! HT(2,0) = R
    ! HT(2,2) = T
    ! HT(2,3) = -SC
    ! HT(2,4) = IU*SQR3*RQS2*S
    ! HT(2,5) = IU*RQS2*(Q-T)
    ! HT(2,6) = IU*RQS3*PM
    ! HT(2,7) = IU*SQR2*RQS3*PZ
    ! ==================================================================
    table(14) = kp_entry(2, 0, KP_R,  cmplx(1.0_dp, 0.0_dp, kind=dp))
    table(15) = kp_entry(2, 2, KP_T,  cmplx(1.0_dp, 0.0_dp, kind=dp))
    table(16) = kp_entry(2, 3, KP_SC, cmplx(-1.0_dp, 0.0_dp, kind=dp))
    table(17) = kp_entry(2, 4, KP_S,  IU*SQR3*RQS2)
    table(18) = kp_entry(2, 5, KP_DIFF, IU*RQS2)
    table(19) = kp_entry(2, 6, KP_PM, IU*RQS3)
    table(20) = kp_entry(2, 7, KP_PZ, IU*SQR2*RQS3)

    ! ==================================================================
    ! Row 3 (HH2): 6 entries
    ! HT(3,1) = R
    ! HT(3,2) = -S
    ! HT(3,3) = Q
    ! HT(3,4) = IU*SQR2*R
    ! HT(3,5) = IU*RQS2*S
    ! HT(3,7) = -PM
    ! ==================================================================
    table(21) = kp_entry(3, 1, KP_R,  cmplx(1.0_dp, 0.0_dp, kind=dp))
    table(22) = kp_entry(3, 2, KP_S,  cmplx(-1.0_dp, 0.0_dp, kind=dp))
    table(23) = kp_entry(3, 3, KP_Q,  cmplx(1.0_dp, 0.0_dp, kind=dp))
    table(24) = kp_entry(3, 4, KP_R,  IU*SQR2)
    table(25) = kp_entry(3, 5, KP_S,  IU*RQS2)
    table(26) = kp_entry(3, 7, KP_PM, cmplx(-1.0_dp, 0.0_dp, kind=dp))

    ! ==================================================================
    ! Row 4 (SO1): 7 entries
    ! HT(4,0) = IU*RQS2*S
    ! HT(4,1) = -IU*RQS2*(Q-T)
    ! HT(4,2) = -IU*SQR3*RQS2*SC
    ! HT(4,3) = -IU*SQR2*RC
    ! HT(4,4) = 0.5*(Q+T)
    ! HT(4,6) = IU*RQS3*PZ
    ! HT(4,7) = IU*SQR2*RQS3*PP
    ! ==================================================================
    table(27) = kp_entry(4, 0, KP_S,  IU*RQS2)
    table(28) = kp_entry(4, 1, KP_DIFF, -IU*RQS2)
    table(29) = kp_entry(4, 2, KP_SC, -IU*SQR3*RQS2)
    table(30) = kp_entry(4, 3, KP_RC, -IU*SQR2)
    table(31) = kp_entry(4, 4, KP_HALF_SUM, cmplx(1.0_dp, 0.0_dp, kind=dp))
    table(32) = kp_entry(4, 6, KP_PZ, IU*RQS3)
    table(33) = kp_entry(4, 7, KP_PP, IU*SQR2*RQS3)

    ! ==================================================================
    ! Row 5 (SO2): 7 entries
    ! HT(5,0) = -IU*SQR2*R
    ! HT(5,1) = IU*SQR3*RQS2*S
    ! HT(5,2) = -IU*RQS2*(Q-T)
    ! HT(5,3) = -IU*RQS2*SC
    ! HT(5,5) = 0.5*(Q+T)
    ! HT(5,6) = SQR2*RQS3*PM
    ! HT(5,7) = -RQS3*PZ
    ! ==================================================================
    table(34) = kp_entry(5, 0, KP_R,  -IU*SQR2)
    table(35) = kp_entry(5, 1, KP_S,  IU*SQR3*RQS2)
    table(36) = kp_entry(5, 2, KP_DIFF, -IU*RQS2)
    table(37) = kp_entry(5, 3, KP_SC, -IU*RQS2)
    table(38) = kp_entry(5, 5, KP_HALF_SUM, cmplx(1.0_dp, 0.0_dp, kind=dp))
    table(39) = kp_entry(5, 6, KP_PM, cmplx(SQR2*RQS3, 0.0_dp, kind=dp))
    table(40) = kp_entry(5, 7, KP_PZ, cmplx(-RQS3, 0.0_dp, kind=dp))

    ! ==================================================================
    ! Row 6 (CB1): 6 entries
    ! HT(6,0) = -IU*PM
    ! HT(6,1) = SQR2*RQS3*PZ
    ! HT(6,2) = -IU*RQS3*PP
    ! HT(6,4) = -IU*RQS3*PZ
    ! HT(6,5) = SQR2*RQS3*PP
    ! HT(6,6) = A
    ! ==================================================================
    table(41) = kp_entry(6, 0, KP_PM, -IU)
    table(42) = kp_entry(6, 1, KP_PZ, cmplx(SQR2*RQS3, 0.0_dp, kind=dp))
    table(43) = kp_entry(6, 2, KP_PP, -IU*RQS3)
    table(44) = kp_entry(6, 4, KP_PZ, -IU*RQS3)
    table(45) = kp_entry(6, 5, KP_PP, cmplx(SQR2*RQS3, 0.0_dp, kind=dp))
    table(46) = kp_entry(6, 6, KP_A,  cmplx(1.0_dp, 0.0_dp, kind=dp))

    ! ==================================================================
    ! Row 7 (CB2): 6 entries
    ! HT(7,1) = -RQS3*PM
    ! HT(7,2) = -IU*SQR2*RQS3*PZ
    ! HT(7,3) = -PP
    ! HT(7,4) = -IU*SQR2*RQS3*PM
    ! HT(7,5) = -RQS3*PZ
    ! HT(7,7) = A
    ! ==================================================================
    table(47) = kp_entry(7, 1, KP_PM, cmplx(-RQS3, 0.0_dp, kind=dp))
    table(48) = kp_entry(7, 2, KP_PZ, -IU*SQR2*RQS3)
    table(49) = kp_entry(7, 3, KP_PP, cmplx(-1.0_dp, 0.0_dp, kind=dp))
    table(50) = kp_entry(7, 4, KP_PM, -IU*SQR2*RQS3)
    table(51) = kp_entry(7, 5, KP_PZ, cmplx(-RQS3, 0.0_dp, kind=dp))
    table(52) = kp_entry(7, 7, KP_A,  cmplx(1.0_dp, 0.0_dp, kind=dp))

  end function build_kp_table

end module hamiltonian_blocks
