module krylov_helpers

  use definitions, only: dp
  use sparse_matrices, only: csr_matrix, csr_spmv
  implicit none

  private
  public :: krylov_chain
  public :: krylov_compare

contains

  ! ------------------------------------------------------------------
  ! Apply k SpMV iterations from a seed vector, storing all intermediate
  ! vectors.  vectors(:, 1) = seed, vectors(:, i+1) = A * vectors(:, i).
  !
  ! This is a fixed-seed power iteration chain (not a Krylov subspace
  ! method).  Used for structural testing: the chain propagates errors
  ! across connected components, catching CSR index-arithmetic bugs
  ! that preserve eigenvalues but corrupt structure.
  ! ------------------------------------------------------------------
  subroutine krylov_chain(A, seed, k, vectors)
    type(csr_matrix), intent(in) :: A
    complex(kind=dp), intent(in) :: seed(:)
    integer, intent(in) :: k
    complex(kind=dp), allocatable, intent(out) :: vectors(:,:)

    integer :: n, i
    complex(kind=dp), parameter :: ONE = cmplx(1.0_dp, 0.0_dp, dp)
    complex(kind=dp), parameter :: ZERO = cmplx(0.0_dp, 0.0_dp, dp)

    n = A%nrows
    allocate(vectors(n, k + 1))
    vectors(:, 1) = seed

    do i = 1, k
      vectors(:, i + 1) = ZERO
      call csr_spmv(A, vectors(:, i), vectors(:, i + 1), ONE, ZERO)
    end do
  end subroutine krylov_chain

  ! ------------------------------------------------------------------
  ! Compare computed Krylov vectors against reference.
  ! Returns .true. if all vectors match within tolerance.
  ! On failure, reports which iteration diverged and the magnitude.
  ! ------------------------------------------------------------------
  function krylov_compare(vectors, reference, tolerance, failing_iteration, &
      max_divergence) result(pass)
    complex(kind=dp), intent(in) :: vectors(:,:)
    complex(kind=dp), intent(in) :: reference(:,:)
    real(kind=dp), intent(in) :: tolerance
    integer, intent(out), optional :: failing_iteration
    real(kind=dp), intent(out), optional :: max_divergence
    logical :: pass

    integer :: i, n, k
    real(kind=dp) :: err, max_err, scale

    n = size(vectors, 1)
    k = size(vectors, 2)

    if (size(reference, 1) /= n .or. size(reference, 2) /= k) then
      pass = .false.
      if (present(failing_iteration)) failing_iteration = 0
      if (present(max_divergence)) max_divergence = huge(1.0_dp)
      return
    end if

    pass = .true.
    max_err = 0.0_dp

    do i = 1, k
      err = maxval(abs(vectors(:, i) - reference(:, i)))
      scale = max(maxval(abs(reference(:, i))), 1.0_dp)
      max_err = max(max_err, err / scale)
      if (err / scale > tolerance) then
        pass = .false.
        if (present(failing_iteration)) failing_iteration = i
        if (present(max_divergence)) max_divergence = err / scale
        return
      end if
    end do

    if (present(failing_iteration)) failing_iteration = 0
    if (present(max_divergence)) max_divergence = max_err
  end function krylov_compare

end module krylov_helpers
