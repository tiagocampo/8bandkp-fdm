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

  ! Numerical perturbation constants
  real(kind=dp), parameter :: mu_B = 5.7883818012e-5_dp  ! Bohr magneton in eV/T
  real(kind=dp), parameter :: g_e = 2.0023193043622_dp    ! Free electron g-factor
  real(kind=dp), parameter :: perturbation_tolerance_default = 1e-12_dp  ! Default convergence tolerance
  real(kind=dp), parameter :: min_magnetic_field = 1e-8_dp      ! Minimum magnetic field [T]
  real(kind=dp), parameter :: max_magnetic_field = 10.0_dp     ! Maximum magnetic field [T]
  real(kind=dp), parameter :: default_step_size = 1e-4_dp      ! Default perturbation step [T]
  real(kind=dp), parameter :: min_step_size = 1e-8_dp          ! Minimum step size [T]
  real(kind=dp), parameter :: max_step_size = 1e-2_dp          ! Maximum step size [T]
  real(kind=dp), parameter :: degeneracy_threshold = 1e-6_dp  ! Energy difference for degenerate states [eV]
  real(kind=dp), parameter :: small_bandgap_threshold = 0.01_dp ! Small band gap threshold [eV]
  real(kind=dp), parameter :: hermitian_threshold = 1e-12_dp  ! Hermitian detection threshold

  complex(kind=dp), parameter :: IU = dcmplx(0.0_dp, 1.0_dp)
  real(kind=dp), parameter :: SQR3 = dsqrt(3.0_dp)
  real(kind=dp), parameter :: SQR2 = dsqrt(2.0_dp)
  real(kind=dp), parameter :: SQR2o3 = dsqrt(2.0_dp/3.0_dp)
  real(kind=dp), parameter :: RQS3 = 1.0_dp/SQR3
  real(kind=dp), parameter :: RQS2 = 1.0_dp/SQR2
  complex(kind=dp), parameter :: ZERO = dcmplx(0.0_dp,0.0_dp)
  complex(kind=dp), parameter :: UM = dcmplx(1.0_dp,0.0_dp)


  type wavevector

    real(kind=dp) :: kx, ky, kz

  end type wavevector

  type paramStruct

    real(kind=dp) :: meff, gamma1, gamma2, gamma3, P, A, deltaSO, EP, Eg, EV, EC, kappa, q

  end type paramStruct

  type group
      integer :: order    ! original order of unsorted data
      real(kind=dp) :: value       ! values to be sorted by
  end type group

  type hilbertspace
    integer :: globalIndex
    integer :: localIndex(2)
    !integer :: qn_nx1, qn_ny1, qn_nx2, qn_ny2
    integer :: basisSpinIndex(2)
    real(kind=dp) :: energies
  end type hilbertspace

  contains

  pure function kronij(i,j)
    integer, intent(in) :: i,j
    integer :: kronij
    kronij = int((float((i+j)-abs(i-j)))/(float((i+j)+abs(i-j))))
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

  recursive function gcd_rec(u, v) result(gcd)
      integer             :: gcd
      integer, intent(in) :: u, v

      if (mod(u, v) /= 0) then
          gcd = gcd_rec(v, mod(u, v))
      else
          gcd = v
      end if
  end function gcd_rec

end module
