module definitions

  implicit none

  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: qp = selected_real_kind(33, 4931)
  integer, parameter :: iknd = selected_int_kind(8)

  real(kind=sp), parameter :: pi_sp = 3.14159265358979323846264338327950288_sp
  real(kind=dp), parameter :: pi_dp = 3.14159265358979323846264338327950288_dp
  real(kind=qp), parameter :: pi_qp = 3.14159265358979323846264338327950288_qp

  real(kind=dp), parameter :: a0 = 0.5291772_dp !bohr radius in AA
  real(kind=dp), parameter :: hbar = 6.582119514_dp*1E-16! eV s
  real(kind=dp), parameter :: const = 3.8099820104685644_dp !ev AA**2
  real(kind=dp), parameter :: c = 2.99792458_dp*1E8*1E10 ! AA / s
  real(kind=dp), parameter :: m0 = 0.5109989461_dp*1E6/c**2 ! eV / c^2
  real(kind=dp), parameter :: e = 1.602176565_dp*1E-19
  real(kind=dp), parameter :: hbar2O2m0 = (hbar)**2 / (2.*m0) !eV AA**2
  real(kind=dp), parameter :: e0 = 8.854187817_dp*1E-12*(1E-9) ! C V-1 nm-1
  real(kind=dp), parameter :: tolerance=1e-7
  logical, parameter :: renormalization = .False.

  complex(kind=dp), parameter :: IU = cmplx(0.0_dp, 1.0_dp, kind=dp)
  real(kind=dp), parameter :: SQR3 = dsqrt(3.0_dp)
  real(kind=dp), parameter :: SQR2 = dsqrt(2.0_dp)
  real(kind=dp), parameter :: SQR2o3 = dsqrt(2.0_dp/3.0_dp)
  real(kind=dp), parameter :: RQS3 = 1.0_dp/SQR3
  real(kind=dp), parameter :: RQS2 = 1.0_dp/SQR2
  complex(kind=dp), parameter :: ZERO = cmplx(0.0_dp, 0.0_dp, kind=dp)
  complex(kind=dp), parameter :: UM = cmplx(1.0_dp, 0.0_dp, kind=dp)


  type wavevector

    real(kind=dp) :: kx, ky, kz

  end type wavevector

  type paramStruct

    real(kind=dp) :: meff, gamma1, gamma2, gamma3, P, A, deltaSO, EP, Eg, EV, EC

  end type paramStruct

  type simulation_config
    integer :: confinement = 0
    integer :: fdStep = 1
    integer :: numLayers = 1
    integer :: numcb = 2
    integer :: numvb = 6
    integer :: evnum = 8
    integer :: waveVectorStep = 100
    integer :: ExternalField = 0
    character(len=2) :: waveVector = 'k0'
    character(len=1) :: confDir = 'n'
    character(len=2) :: EFtype = '  '
    real(kind=dp) :: waveVectorMax = 0.0_dp
    real(kind=dp) :: Evalue = 0.0_dp
    real(kind=dp) :: totalSize = 0.0_dp
    real(kind=dp) :: delta = 0.0_dp
    real(kind=dp) :: dz = 0.0_dp
    real(kind=dp), allocatable :: startPos(:)
    real(kind=dp), allocatable :: endPos(:)
    real(kind=dp), allocatable :: z(:)
    integer, allocatable :: intStartPos(:)
    integer, allocatable :: intEndPos(:)
    character(len=255), allocatable :: materialN(:)
    type(paramStruct), allocatable :: params(:)
  end type simulation_config

  type group
      integer :: order    ! original order of unsorted data
      real(kind=dp) :: value       ! values to be sorted by
  end type group


  contains

  pure function kronij(i,j)
    integer, intent(in) :: i,j
    integer :: kronij
    kronij = merge(1, 0, i == j)
  end function

  subroutine tick(t)
      integer, intent(OUT) :: t

      call system_clock(t)
  end subroutine tick

  ! returns time in seconds from now to time described by t
  real function tock(t)
      integer, intent(in) :: t
      integer :: now, clock_rate

      call system_clock(now,clock_rate)

      tock = real(now - t)/real(clock_rate)
  end function tock


end module
