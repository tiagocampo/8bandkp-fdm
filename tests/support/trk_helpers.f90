module trk_helpers

  use definitions, only: dp, wavevector, paramStruct, spatial_grid
  use hamiltonianConstructor, only: ZB8bandBulk
  use sparse_matrices, only: csr_matrix, csr_spmv, csr_free
  use linalg, only: zheevd
  implicit none

  private
  public :: compute_double_commutator_bulk
  public :: build_bulk_velocity_matrix
  public :: compute_band_curvature_bulk
  public :: compute_trk_sum_explicit
  public :: compute_double_commutator_csr
  public :: compute_trk_sum_explicit_csr

contains

  ! ------------------------------------------------------------------
  ! Compute d2H/dk2 for the bulk 8x8 Hamiltonian using central finite
  ! differences:  d2H/dk2 = (H(+d) - 2H(0) + H(-d)) / d^2
  !
  ! Arguments:
  !   d2Hdk2    -- output 8x8 matrix (second derivative)
  !   params    -- material parameters (size 1)
  !   direction -- 1=x, 2=y, 3=z
  ! ------------------------------------------------------------------
  subroutine compute_double_commutator_bulk(d2Hdk2, params, direction)
    complex(kind=dp), intent(out) :: d2Hdk2(8, 8)
    type(paramStruct), intent(in) :: params(1)
    integer, intent(in) :: direction  ! 1=x, 2=y, 3=z

    real(kind=dp), parameter :: delta = 1.0e-4_dp  ! Angstrom^-1

    complex(kind=dp) :: H_plus(8, 8), H_zero(8, 8), H_minus(8, 8)
    type(wavevector) :: wv_plus, wv_zero, wv_minus

    ! k=0 wavevector
    wv_zero%kx = 0.0_dp
    wv_zero%ky = 0.0_dp
    wv_zero%kz = 0.0_dp

    ! +/-d wavevectors
    wv_plus  = wv_zero
    wv_minus = wv_zero

    select case(direction)
    case(1)
      wv_plus%kx  = delta
      wv_minus%kx = -delta
    case(2)
      wv_plus%ky  = delta
      wv_minus%ky = -delta
    case(3)
      wv_plus%kz  = delta
      wv_minus%kz = -delta
    end select

    ! Build Hamiltonians at three k-points
    H_zero = cmplx(0.0_dp, 0.0_dp, kind=dp)
    H_plus = cmplx(0.0_dp, 0.0_dp, kind=dp)
    H_minus = cmplx(0.0_dp, 0.0_dp, kind=dp)

    call ZB8bandBulk(H_zero, wv_zero, params)
    call ZB8bandBulk(H_plus, wv_plus, params)
    call ZB8bandBulk(H_minus, wv_minus, params)

    ! Central finite difference for second derivative
    d2Hdk2 = (H_plus - 2.0_dp * H_zero + H_minus) / (delta * delta)
  end subroutine compute_double_commutator_bulk

  ! ------------------------------------------------------------------
  ! Build the velocity matrix v_alpha = dH/dk_alpha for bulk at k=0
  ! using central finite differences: dH/dk = (H(+d) - H(-d)) / (2d)
  !
  ! Arguments:
  !   vel       -- output 8x8 velocity matrix
  !   params    -- material parameters (size 1)
  !   direction -- 1=x, 2=y, 3=z
  ! ------------------------------------------------------------------
  subroutine build_bulk_velocity_matrix(vel, params, direction)
    complex(kind=dp), intent(out) :: vel(8, 8)
    type(paramStruct), intent(in) :: params(1)
    integer, intent(in) :: direction

    real(kind=dp), parameter :: delta = 1.0e-4_dp

    complex(kind=dp) :: H_plus(8, 8), H_minus(8, 8)
    type(wavevector) :: wv_plus, wv_minus

    wv_plus%kx = 0.0_dp
    wv_plus%ky = 0.0_dp
    wv_plus%kz = 0.0_dp
    wv_minus = wv_plus

    select case(direction)
    case(1)
      wv_plus%kx  = delta
      wv_minus%kx = -delta
    case(2)
      wv_plus%ky  = delta
      wv_minus%ky = -delta
    case(3)
      wv_plus%kz  = delta
      wv_minus%kz = -delta
    end select

    H_plus  = cmplx(0.0_dp, 0.0_dp, kind=dp)
    H_minus = cmplx(0.0_dp, 0.0_dp, kind=dp)

    call ZB8bandBulk(H_plus,  wv_plus,  params)
    call ZB8bandBulk(H_minus, wv_minus, params)

    vel = (H_plus - H_minus) / (2.0_dp * delta)
  end subroutine build_bulk_velocity_matrix

  ! ------------------------------------------------------------------
  ! Compute band curvature d2E_n/dk2 using finite differences on
  ! eigenvalues: d2E/dk2 = (E(+d) - 2E(0) + E(-d)) / d^2
  !
  ! Arguments:
  !   params      -- material parameters (size 1)
  !   direction   -- 1=x, 2=y, 3=z
  !   band_idx    -- which band (1-8, in ascending eigenvalue order)
  !   curvature   -- output: d2E/dk2 for the specified band
  ! ------------------------------------------------------------------
  subroutine compute_band_curvature_bulk(params, direction, band_idx, curvature)
    type(paramStruct), intent(in) :: params(1)
    integer, intent(in) :: direction, band_idx
    real(kind=dp), intent(out) :: curvature

    real(kind=dp), parameter :: delta = 1.0e-4_dp

    complex(kind=dp) :: HT(8, 8)
    type(wavevector) :: wv
    real(kind=dp) :: evals_zero(8), evals_plus(8), evals_minus(8)
    complex(kind=dp), allocatable :: work(:)
    real(kind=dp), allocatable :: rwork(:)
    integer, allocatable :: iwork(:)
    integer :: lwork, lrwork, liwork, info

    ! Diagonalize at k=0
    wv%kx = 0.0_dp; wv%ky = 0.0_dp; wv%kz = 0.0_dp
    HT = cmplx(0.0_dp, 0.0_dp, kind=dp)
    call ZB8bandBulk(HT, wv, params)
    allocate(work(1), rwork(1), iwork(1))
    call zheevd('N', 'U', 8, HT, 8, evals_zero, work, -1, rwork, -1, iwork, -1, info)
    lwork = int(real(work(1))); lrwork = int(rwork(1)); liwork = iwork(1)
    deallocate(work, rwork, iwork)
    allocate(work(lwork), rwork(lrwork), iwork(liwork))
    call zheevd('N', 'U', 8, HT, 8, evals_zero, work, lwork, rwork, lrwork, &
      iwork, liwork, info)
    deallocate(work, rwork, iwork)

    ! Diagonalize at k=+d
    select case(direction)
    case(1); wv%kx = delta
    case(2); wv%ky = delta
    case(3); wv%kz = delta
    end select
    HT = cmplx(0.0_dp, 0.0_dp, kind=dp)
    call ZB8bandBulk(HT, wv, params)
    allocate(work(1), rwork(1), iwork(1))
    call zheevd('N', 'U', 8, HT, 8, evals_plus, work, -1, rwork, -1, iwork, -1, info)
    lwork = int(real(work(1))); lrwork = int(rwork(1)); liwork = iwork(1)
    deallocate(work, rwork, iwork)
    allocate(work(lwork), rwork(lrwork), iwork(liwork))
    call zheevd('N', 'U', 8, HT, 8, evals_plus, work, lwork, rwork, lrwork, &
      iwork, liwork, info)
    deallocate(work, rwork, iwork)

    ! Diagonalize at k=-d
    select case(direction)
    case(1); wv%kx = -delta
    case(2); wv%ky = -delta
    case(3); wv%kz = -delta
    end select
    HT = cmplx(0.0_dp, 0.0_dp, kind=dp)
    call ZB8bandBulk(HT, wv, params)
    allocate(work(1), rwork(1), iwork(1))
    call zheevd('N', 'U', 8, HT, 8, evals_minus, work, -1, rwork, -1, iwork, -1, info)
    lwork = int(real(work(1))); lrwork = int(rwork(1)); liwork = iwork(1)
    deallocate(work, rwork, iwork)
    allocate(work(lwork), rwork(lrwork), iwork(liwork))
    call zheevd('N', 'U', 8, HT, 8, evals_minus, work, lwork, rwork, lrwork, &
      iwork, liwork, info)
    deallocate(work, rwork, iwork)

    ! Central finite difference
    curvature = (evals_plus(band_idx) - 2.0_dp * evals_zero(band_idx) + &
      evals_minus(band_idx)) / (delta * delta)
  end subroutine compute_band_curvature_bulk

  ! ------------------------------------------------------------------
  ! Compute the explicit TRK sum-over-states for the 8-band bulk system:
  !
  !   trk_sum = sum_{m != n} |<m|v|n>|^2 / (E_n - E_m)
  !
  ! This is the sum-over-states contribution to the effective mass
  ! from 2nd-order k.p perturbation theory.
  !
  ! Arguments:
  !   evals      -- eigenvalues (8), ascending order
  !   evecs      -- eigenvectors as columns (8 x 8)
  !   vel        -- velocity matrix (8 x 8)
  !   state_idx  -- index of initial state (1-based)
  !   trk_sum    -- output: the sum-over-states value
  ! ------------------------------------------------------------------
  subroutine compute_trk_sum_explicit(evals, evecs, vel, state_idx, trk_sum)
    real(kind=dp), intent(in) :: evals(8)
    complex(kind=dp), intent(in) :: evecs(8, 8)
    complex(kind=dp), intent(in) :: vel(8, 8)
    integer, intent(in) :: state_idx
    real(kind=dp), intent(out) :: trk_sum

    integer :: m, k
    real(kind=dp) :: E0
    complex(kind=dp) :: mel

    E0 = evals(state_idx)
    trk_sum = 0.0_dp

    do m = 1, 8
      if (abs(evals(m) - E0) < 1.0e-10_dp) cycle  ! skip degenerate states
      ! <m|v|state_idx> = sum_k evecs(k,m)^* * sum_j vel(k,j) * evecs(j,state_idx)
      mel = cmplx(0.0_dp, 0.0_dp, kind=dp)
      do k = 1, 8
        mel = mel + conjg(evecs(k, m)) * sum(vel(k, :) * evecs(:, state_idx))
      end do
      trk_sum = trk_sum + abs(mel)**2 / (E0 - evals(m))
    end do
  end subroutine compute_trk_sum_explicit

  ! ------------------------------------------------------------------
  ! Compute the TRK sum-over-states for CSR velocity matrices:
  !
  !   trk_sum = sum_{m != n} |<m|v|n>|^2 / (E_n - E_m)
  !
  ! Uses CSR SpMV for the matrix elements: <m|v|n> = <m| v|n>.
  !
  ! Arguments:
  !   evals      -- eigenvalues (n), ascending order
  !   evecs      -- eigenvectors as columns (n x n)
  !   vel_csr    -- CSR velocity matrix
  !   state_idx  -- index of initial state (1-based)
  !   trk_sum    -- output: the sum-over-states value
  ! ------------------------------------------------------------------
  subroutine compute_trk_sum_explicit_csr(evals, evecs, vel_csr, state_idx, trk_sum)
    real(kind=dp), intent(in), contiguous    :: evals(:)
    complex(kind=dp), intent(in), contiguous :: evecs(:,:)
    type(csr_matrix), intent(in)    :: vel_csr
    integer, intent(in)             :: state_idx
    real(kind=dp), intent(out)      :: trk_sum

    integer :: m, n
    real(kind=dp) :: E0
    complex(kind=dp) :: mel
    complex(kind=dp), allocatable :: vn(:), mvn(:)

    n = size(evals)
    E0 = evals(state_idx)

    allocate(vn(n), mvn(n))
    vn = cmplx(0.0_dp, 0.0_dp, kind=dp)

    ! Compute v|n>
    call csr_spmv(vel_csr, evecs(:, state_idx), vn, &
      cmplx(1.0_dp, 0.0_dp, kind=dp), cmplx(0.0_dp, 0.0_dp, kind=dp))

    trk_sum = 0.0_dp

    do m = 1, n
      if (m == state_idx) cycle
      if (abs(evals(m) - E0) < 1.0e-10_dp) cycle  ! skip degenerate states

      ! <m|v|n> = evecs(:,m)^H . vn
      mel = cmplx(0.0_dp, 0.0_dp, kind=dp)
      do concurrent (integer :: k = 1:n)
        mvn(k) = conjg(evecs(k, m)) * vn(k)
      end do
      ! Sum-reduction: must use regular do (not do concurrent) since iterations
      ! are not independent -- they accumulate into mel.
      do k = 1, n
        mel = mel + mvn(k)
      end do

      trk_sum = trk_sum + abs(mel)**2 / (E0 - evals(m))
    end do

    deallocate(vn, mvn)
  end subroutine compute_trk_sum_explicit_csr

  ! ------------------------------------------------------------------
  ! Compute double commutator from CSR Hamiltonian:
  !
  !   dc = 1/2 * <n| [r, [r, H]] |n>
  !      = 1/2 * sum_{i,j} (x_i - x_j)^2 * H_ij * psi0*(i) * psi0(j)
  !
  ! Loops over all CSR entries. The double commutator [r,[r,H]] has
  ! matrix elements (z_i - z_j)^2 * H_ij.  For Hermitian H, summing
  ! over all entries and taking the real part gives the correct result.
  !
  ! The TRK sum rule states:
  !   sum_{m!=n} |<m|v|n>|^2 / (E_n - E_m) = 1/2 * <n| [r, [r, H]] |n>
  ! so dc_value should equal the explicit sum-over-states.
  !
  ! Arguments:
  !   H_csr      -- CSR Hamiltonian (Hermitian)
  !   grid       -- spatial grid with positions
  !   psi0       -- ground-state eigenvector (8*N elements)
  !   direction  -- 'z' for QW z-velocity, 'x' or 'y' for wire
  !   dc_value   -- output: the double commutator value
  ! ------------------------------------------------------------------
  subroutine compute_double_commutator_csr(H_csr, grid, psi0, direction, dc_value)
    type(csr_matrix), intent(in)    :: H_csr
    type(spatial_grid), intent(in)  :: grid
    complex(kind=dp), intent(in), contiguous :: psi0(:)
    character(len=1), intent(in)    :: direction  ! 'x', 'y', 'z'
    real(kind=dp), intent(out)      :: dc_value

    integer :: k, row, col, Ngrid, sp_row, sp_col
    real(kind=dp) :: pos_diff

    Ngrid = grid%npoints()
    dc_value = 0.0_dp

    do row = 1, H_csr%nrows
      do k = H_csr%rowptr(row), H_csr%rowptr(row + 1) - 1
        col = H_csr%colind(k)

        sp_row = mod(row - 1, Ngrid) + 1
        sp_col = mod(col - 1, Ngrid) + 1

        select case(direction)
        case('z')
          pos_diff = grid%z(sp_row) - grid%z(sp_col)
        case('x')
          pos_diff = grid%coords(1, sp_row) - grid%coords(1, sp_col)
        case('y')
          pos_diff = grid%coords(2, sp_row) - grid%coords(2, sp_col)
        case default
          pos_diff = 0.0_dp
        end select

        if (abs(pos_diff) < 1.0e-15_dp) cycle

        ! sum_{ij} (x_i-x_j)^2 * H_ij * psi0*(i)*psi0(j)
        ! For Hermitian H: each (i,j) pair appears as H_ij and H_ji = H_ij*,
        ! and psi0*(i)*psi0(j) and psi0*(j)*psi0(i) = conjg of above.
        ! So summing all entries and taking real part gives correct answer.
        dc_value = dc_value + pos_diff * pos_diff * &
          real(H_csr%values(k) * conjg(psi0(row)) * psi0(col), kind=dp)
      end do
    end do

    dc_value = 0.5_dp * dc_value
  end subroutine compute_double_commutator_csr

end module trk_helpers
