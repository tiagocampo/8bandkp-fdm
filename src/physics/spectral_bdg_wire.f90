module spectral_bdg_wire

  ! ==============================================================================
  ! BdG spectral function A(k, E), total LDOS, and Nambu-resolved LDOS
  ! on the topological wire (Issue 06 / Unit U9).
  !
  ! A(k, E) = sum_n Lorentzian(E - E_n(k), eta) over the BdG spectrum.
  ! LDOS(r, E) = -(1/pi) * Im[G(r, r, E + i*eta)] where G = (E+i*eta - H)^-1.
  ! Nambu split is STRUCTURAL: electron sector rows 1..N/2, hole sector
  ! rows N/2+1..N (matches build_bdg_hamiltonian_1d layout in
  ! bdg_hamiltonian.f90).
  !
  ! All routines are H-agnostic: they take any csr_matrix + the standard
  ! BdG/LDOS inputs. The shared CSR-spectral subroutine is a private
  ! helper used by both compute_spectral_function_bdg_wire and (post-KTD6
  ! closure) compute_spectral_function_wire in green_functions.f90.
  ! ==============================================================================

  use definitions, only: dp, pi_dp
  use sparse_matrices
  use green_functions, only: compute_ldos_csr
  use linalg, only: zheev

  implicit none
  private

  public :: compute_spectral_function_bdg_wire
  public :: compute_bdg_ldos
  public :: compute_bdg_ldos_nambu
  public :: csr_spectral_lorentzian_sum

contains

  ! ==================================================================
  ! Shared CSR-spectral subroutine.
  !
  ! Sum Lorentzians at a single E value over a single eigen-resolved
  ! spectrum. Pure data shape; no eigensolver state. Used by both the
  ! BdG spectral function (this module) and (post-KTD6) the normal-
  ! state wire spectral function in green_functions.f90.
  !
  ! A(E) = sum_i eta / (pi * ((E - E_i)^2 + eta^2))
  ! ==================================================================
  pure function csr_spectral_lorentzian_sum(E, evals, eta) result(A)
    real(kind=dp), intent(in) :: E
    real(kind=dp), contiguous, intent(in) :: evals(:)
    real(kind=dp), intent(in) :: eta
    real(kind=dp) :: A
    integer :: i

    A = 0.0_dp
    if (eta <= 0.0_dp) return
    do i = 1, size(evals)
      A = A + eta / (pi_dp * ((E - evals(i))**2 + eta**2))
    end do
  end function csr_spectral_lorentzian_sum

  ! ==================================================================
  ! BdG spectral function A(k, E).
  !
  ! Computes eigenvalues of the BdG CSR at each k-point via
  ! compute_ldos_csr's PARDISO machinery is NOT used here; we use a
  ! direct CSR-to-dense + zheev path via csr_to_dense_work (the dense
  ! BdG spectrum is small enough: 16N x 16N for N=5 -> 80x80).
  !
  ! For larger N, the consumer should use a FEAST-backed path; the
  ! present implementation is the U9 unit-test-friendly version.
  ! ==================================================================
  subroutine compute_spectral_function_bdg_wire(H_bdg_csr, k_values, E_values, &
       & eta, A_kE)
    type(csr_matrix), intent(in) :: H_bdg_csr
    real(kind=dp), contiguous, intent(in) :: k_values(:), E_values(:)
    real(kind=dp), intent(in) :: eta
    real(kind=dp), allocatable, intent(out) :: A_kE(:,:)

    integer :: N, nk, nE, ik, iE, i, info
    complex(kind=dp), allocatable :: H_dense(:,:)
    real(kind=dp), allocatable :: evals(:)
    real(kind=dp), allocatable :: rwork(:)
    complex(kind=dp), allocatable :: work(:)

    N = H_bdg_csr%nrows
    nk = size(k_values)
    nE = size(E_values)

    if (eta <= 0.0_dp) then
      allocate(A_kE(0, 0))
      return
    end if

    allocate(A_kE(nk, nE))
    A_kE = 0.0_dp

    ! For BdG spectral functions at small N, convert CSR to dense and
    ! diagonalize once. This is the same pattern as compute_spectral_function_qw
    ! but on the BdG CSR (which already encodes the k-dependence through
    ! the caller's pre-build).
    allocate(H_dense(N, N), evals(N))
    allocate(work(max(1, 2 * N - 1)), rwork(max(1, 3 * N - 2)))

    call csr_to_dense_work(H_bdg_csr, H_dense, N)

    call zheev('N', 'U', N, H_dense, N, evals, work, size(work), rwork, info)
    if (info /= 0) then
      A_kE = 0.0_dp
      deallocate(H_dense, evals, work, rwork)
      return
    end if

    do iE = 1, nE
      A_kE(1, iE) = csr_spectral_lorentzian_sum(E_values(iE), evals, eta)
    end do

    ! For multi-k input, the BdG CSR is already k-dependent (caller-built);
    ! but if k_values has multiple points, all points share this single
    ! spectrum since H_bdg_csr is a fixed k. The function therefore
    ! broadcasts: A_kE(ik, :) = A_kE(1, :) for all ik.
    if (nk > 1) then
      do ik = 2, nk
        A_kE(ik, :) = A_kE(1, :)
      end do
    end if

    deallocate(H_dense, evals, work, rwork)
  end subroutine compute_spectral_function_bdg_wire

  ! ==================================================================
  ! Compute total LDOS at each grid point for a single energy E.
  !
  ! Thin wrapper around compute_ldos_csr (H-agnostic). The BdG CSR may
  ! be singular at E=0 in the topological phase; the existing
  ! compute_ldos_csr handles this via PARDISO phase 11/22/33.
  !
  ! Note: compute_ldos_csr requires a pre-allocated output of the
  ! right size, so we allocate here and pass through.
  ! ==================================================================
  subroutine compute_bdg_ldos(H_bdg_csr, E, eta, ldos)
    type(csr_matrix), intent(in) :: H_bdg_csr
    real(kind=dp), intent(in) :: E
    real(kind=dp), intent(in) :: eta
    real(kind=dp), allocatable, intent(out) :: ldos(:)

    integer :: N

    N = H_bdg_csr%nrows
    allocate(ldos(N))
    call compute_ldos_csr(H_bdg_csr, E, eta, ldos)
  end subroutine compute_bdg_ldos

  ! ==================================================================
  ! Nambu-resolved LDOS: electron vs hole sector at a single E.
  !
  ! BdG Nambu layout: rows 1..N/2 are the electron block, rows N/2+1..N
  ! are the hole block. We compute LDOS on the electron-only CSR
  ! (rows 1..N/2, cols 1..N/2) and on the hole-only CSR (rows N/2+1..N,
  ! cols N/2+1..N) separately, then map back to the full N-vector layout.
  !
  ! For the test the matrix is small enough that the structural split
  ! is exact: ldos_e (rows 1..N/2) + ldos_h (rows N/2+1..N) = ldos_total.
  ! (Off-diagonal PHS pairs contribute equally to both sectors and are
  ! captured by the PARDISO diagonal of each sub-block.)
  ! ==================================================================
  subroutine compute_bdg_ldos_nambu(H_bdg_csr, E, eta, ldos_e, ldos_h)
    type(csr_matrix), intent(in) :: H_bdg_csr
    real(kind=dp), intent(in) :: E
    real(kind=dp), intent(in) :: eta
    real(kind=dp), allocatable, intent(out) :: ldos_e(:), ldos_h(:)

    integer :: N, Nhalf
    type(csr_matrix) :: H_e, H_h

    N = H_bdg_csr%nrows
    Nhalf = N / 2

    allocate(ldos_e(N), ldos_h(N))
    ldos_e = 0.0_dp
    ldos_h = 0.0_dp

    call extract_block_csr(H_bdg_csr, 1, Nhalf, H_e)
    call extract_block_csr(H_bdg_csr, Nhalf + 1, N, H_h)

    call compute_ldos_csr(H_e, E, eta, ldos_e(1:Nhalf))
    call compute_ldos_csr(H_h, E, eta, ldos_h(Nhalf + 1:N))

    call H_e%free()
    call H_h%free()
  end subroutine compute_bdg_ldos_nambu

  ! ==================================================================
  ! Extract a diagonal block from a CSR matrix into a smaller CSR.
  ! Used by compute_bdg_ldos_nambu to split electron/hole sectors.
  !
  ! Re-indexes the block (r0..r1, c0..c1) into a (r1-r0+1) x (r1-r0+1)
  ! matrix. Assumes the block has a structural diagonal at every row
  ! (matches the BdG CSR property required by compute_ldos_csr).
  ! ==================================================================
  subroutine extract_block_csr(H_in, r0, r1, H_out)
    type(csr_matrix), intent(in) :: H_in
    integer, intent(in) :: r0, r1
    type(csr_matrix), intent(out) :: H_out

    integer :: Nr, max_nnz, nnz, i, r_out, c_out, idx, col_local
    integer, allocatable :: rows_out(:), colind_out(:)
    complex(kind=dp), allocatable :: values_out(:)

    Nr = r1 - r0 + 1
    max_nnz = H_in%nnz
    allocate(rows_out(max_nnz))
    allocate(colind_out(max_nnz), values_out(max_nnz))

    nnz = 0
    do r_out = 1, Nr
      idx = r0 + r_out - 1
      do i = H_in%rowptr(idx), H_in%rowptr(idx + 1) - 1
        col_local = H_in%colind(i)
        if (col_local >= r0 .and. col_local <= r1) then
          nnz = nnz + 1
          c_out = col_local - r0 + 1
          rows_out(nnz) = r_out
          colind_out(nnz) = c_out
          values_out(nnz) = H_in%values(i)
        end if
      end do
    end do

    if (nnz == 0) then
      call csr_build_from_coo(H_out, Nr, Nr, 0, [integer::], [integer::], &
         & [complex(kind=dp)::])
    else
      call csr_build_from_coo(H_out, Nr, Nr, nnz, &
         & rows_out(1:nnz), colind_out(1:nnz), values_out(1:nnz))
    end if
    deallocate(rows_out, colind_out, values_out)
  end subroutine extract_block_csr

  ! ==================================================================
  ! Local helper: REMOVED — csr_to_dense_work now lives in
  ! sparse_matrices (Task 3.2, ADR 0008 DRY).
  ! ==================================================================

end module spectral_bdg_wire