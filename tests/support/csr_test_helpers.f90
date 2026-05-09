module csr_test_helpers

  use definitions, only: dp
  use sparse_matrices, only: csr_matrix, csr_free
  implicit none

  private
  public :: assert_csr_structural_invariants
  public :: check_rowptr_monotonic
  public :: check_colind_sorted
  public :: check_colind_bounds
  public :: check_nnz_consistency
  public :: check_diagonal_present
  public :: csr_to_dense
  public :: csr_hermitian_error
  public :: csr_interior_symmetric_error
  public :: csr_interior_hermitian_error

contains

  ! ------------------------------------------------------------------
  ! Primary fixture: validate all CSR structural invariants.
  ! If any invariant is violated, prints diagnostic and stops.
  !
  ! Invariants checked:
  !   1. rowptr is monotonically non-decreasing, 1-based
  !   2. colind sorted within each row
  !   3. colind entries in [1, ncols]
  !   4. nnz = rowptr(nrows+1) - 1
  !   5. (optional) diagonal entry present in every row
  ! ------------------------------------------------------------------
  subroutine assert_csr_structural_invariants(mat, require_diagonal)
    type(csr_matrix), intent(in) :: mat
    logical, intent(in), optional :: require_diagonal

    if (.not. check_rowptr_monotonic(mat)) stop 1
    if (.not. check_colind_sorted(mat))    stop 1
    if (.not. check_colind_bounds(mat))    stop 1
    if (.not. check_nnz_consistency(mat))  stop 1

    if (present(require_diagonal)) then
      if (require_diagonal) then
        if (.not. check_diagonal_present(mat)) stop 1
      end if
    end if
  end subroutine assert_csr_structural_invariants

  ! ------------------------------------------------------------------
  ! Check rowptr is monotonically non-decreasing and 1-based.
  ! rowptr(1) must be 1, rowptr(i+1) >= rowptr(i).
  ! ------------------------------------------------------------------
  function check_rowptr_monotonic(mat) result(ok)
    type(csr_matrix), intent(in) :: mat
    logical :: ok
    integer :: i

    ok = .true.
    if (mat%nrows < 1) return

    if (mat%rowptr(1) /= 1) then
      print *, 'FAIL: rowptr(1) = ', mat%rowptr(1), ' expected 1'
      ok = .false.
      return
    end if

    do i = 1, mat%nrows
      if (mat%rowptr(i + 1) < mat%rowptr(i)) then
        print *, 'FAIL: rowptr not monotonic at row ', i, &
          ' rowptr(i)=', mat%rowptr(i), ' rowptr(i+1)=', mat%rowptr(i + 1)
        ok = .false.
        return
      end if
    end do
  end function check_rowptr_monotonic

  ! ------------------------------------------------------------------
  ! Check colind is sorted in ascending order within each row.
  ! ------------------------------------------------------------------
  function check_colind_sorted(mat) result(ok)
    type(csr_matrix), intent(in) :: mat
    logical :: ok
    integer :: i, k

    ok = .true.
    do i = 1, mat%nrows
      do k = mat%rowptr(i) + 1, mat%rowptr(i + 1) - 1
        if (mat%colind(k) < mat%colind(k - 1)) then
          print *, 'FAIL: colind not sorted at row ', i, &
            ' position ', k, ' colind(k-1)=', mat%colind(k - 1), &
            ' colind(k)=', mat%colind(k)
          ok = .false.
          return
        end if
      end do
    end do
  end function check_colind_sorted

  ! ------------------------------------------------------------------
  ! Check all colind entries are in [1, ncols].
  ! ------------------------------------------------------------------
  function check_colind_bounds(mat) result(ok)
    type(csr_matrix), intent(in) :: mat
    logical :: ok
    integer :: k

    ok = .true.
    do k = 1, mat%nnz
      if (mat%colind(k) < 1 .or. mat%colind(k) > mat%ncols) then
        print *, 'FAIL: colind(', k, ') = ', mat%colind(k), &
          ' out of range [1, ', mat%ncols, ']'
        ok = .false.
        return
      end if
    end do
  end function check_colind_bounds

  ! ------------------------------------------------------------------
  ! Check nnz = rowptr(nrows+1) - 1.
  ! ------------------------------------------------------------------
  function check_nnz_consistency(mat) result(ok)
    type(csr_matrix), intent(in) :: mat
    logical :: ok
    integer :: expected_nnz

    ok = .true.
    expected_nnz = mat%rowptr(mat%nrows + 1) - 1
    if (mat%nnz /= expected_nnz) then
      print *, 'FAIL: nnz = ', mat%nnz, ' but rowptr(n+1)-1 = ', expected_nnz
      ok = .false.
    end if
  end function check_nnz_consistency

  ! ------------------------------------------------------------------
  ! Check diagonal entry present in every row (square matrices only).
  ! ------------------------------------------------------------------
  function check_diagonal_present(mat) result(ok)
    type(csr_matrix), intent(in) :: mat
    logical :: ok
    integer :: i, k
    logical :: found

    ok = .true.
    do i = 1, mat%nrows
      found = .false.
      do k = mat%rowptr(i), mat%rowptr(i + 1) - 1
        if (mat%colind(k) == i) then
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        print *, 'FAIL: no diagonal entry in row ', i
        ok = .false.
        return
      end if
    end do
  end function check_diagonal_present

  ! ------------------------------------------------------------------
  ! Convert CSR matrix to dense (nrows x ncols) array.
  ! ------------------------------------------------------------------
  subroutine csr_to_dense(mat, dense)
    type(csr_matrix), intent(in) :: mat
    complex(kind=dp), intent(out) :: dense(:,:)
    integer :: row, k

    dense = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do row = 1, mat%nrows
      do k = mat%rowptr(row), mat%rowptr(row + 1) - 1
        dense(row, mat%colind(k)) = mat%values(k)
      end do
    end do
  end subroutine csr_to_dense

  ! ------------------------------------------------------------------
  ! Compute max Hermiticity error: max |dense(i,j) - conjg(dense(j,i))|.
  ! ------------------------------------------------------------------
  function csr_hermitian_error(mat) result(max_err)
    type(csr_matrix), intent(in) :: mat
    real(kind=dp) :: max_err
    complex(kind=dp), allocatable :: dense(:,:)
    integer :: i, j

    allocate(dense(mat%nrows, mat%ncols))
    call csr_to_dense(mat, dense)
    max_err = 0.0_dp

    do j = 1, mat%ncols
      do i = 1, mat%nrows
        max_err = max(max_err, abs(dense(i,j) - conjg(dense(j,i))))
      end do
    end do

    deallocate(dense)
  end function csr_hermitian_error

  ! ------------------------------------------------------------------
  ! Compute max symmetry error for interior-interior pairs.
  ! ------------------------------------------------------------------
  function csr_interior_symmetric_error(mat, nx, ny) result(max_err)
    type(csr_matrix), intent(in) :: mat
    integer, intent(in) :: nx, ny
    real(kind=dp) :: max_err
    complex(kind=dp), allocatable :: dense(:,:)
    integer :: j, ix, iy, ij, jx, jy

    allocate(dense(mat%nrows, mat%ncols))
    call csr_to_dense(mat, dense)
    max_err = 0.0_dp

    do iy = 2, ny - 1
      do ix = 2, nx - 1
        ij = (iy - 1) * nx + ix
        do jy = 2, ny - 1
          do jx = 2, nx - 1
            j = (jy - 1) * nx + jx
            max_err = max(max_err, abs(dense(ij, j) - dense(j, ij)))
          end do
        end do
      end do
    end do

    deallocate(dense)
  end function csr_interior_symmetric_error

  ! ------------------------------------------------------------------
  ! Compute max Hermiticity error for interior-interior pairs.
  ! ------------------------------------------------------------------
  function csr_interior_hermitian_error(mat, nx, ny) result(max_err)
    type(csr_matrix), intent(in) :: mat
    integer, intent(in) :: nx, ny
    real(kind=dp) :: max_err
    complex(kind=dp), allocatable :: dense(:,:)
    integer :: j, ix, iy, ij, jx, jy

    allocate(dense(mat%nrows, mat%ncols))
    call csr_to_dense(mat, dense)
    max_err = 0.0_dp

    do iy = 2, ny - 1
      do ix = 2, nx - 1
        ij = (iy - 1) * nx + ix
        do jy = 2, ny - 1
          do jx = 2, nx - 1
            j = (jy - 1) * nx + jx
            max_err = max(max_err, abs(dense(ij, j) - conjg(dense(j, ij))))
          end do
        end do
      end do
    end do

    deallocate(dense)
  end function csr_interior_hermitian_error

end module csr_test_helpers
