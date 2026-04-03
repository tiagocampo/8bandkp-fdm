module sparse_matrices

  ! ==============================================================================
  ! CSR sparse matrix type and utility routines for complex-valued sparse
  ! matrices.  Used for 2D quantum-wire Hamiltonian assembly (Phase 0b) and
  ! Kronecker product construction (Phase 0c).
  !
  ! Build strategy: caller collects COO triplets externally, then passes them
  ! to csr_build_from_coo which sorts, merges duplicates, and produces CSR.
  !
  ! All routines are library code -- no I/O.
  ! ==============================================================================

  use definitions, only: dp, iknd

  implicit none

  private
  public :: csr_matrix
  public :: csr_init, csr_free, csr_build_from_coo
  public :: kron_dense_dense, kron_dense_eye, kron_eye_dense, kron_dense_dense_1d
  public :: csr_set_values_from_coo, csr_build_from_coo_cached
  public :: csr_add, csr_scale, csr_apply_variable_coeff

  ! ------------------------------------------------------------------
  ! Compressed Sparse Row (CSR) matrix type for complex-valued
  ! sparse matrices.  1-based indexing throughout (compatible with
  ! MKL and LAPACK conventions used in this codebase).
  !
  ! rowptr(i) .. rowptr(i+1)-1  gives the range of entries in
  ! colind(:) and values(:) that belong to row i.
  ! ------------------------------------------------------------------
  type :: csr_matrix
    integer :: nrows = 0
    integer :: ncols = 0
    integer :: nnz   = 0
    complex(kind=dp), allocatable :: values(:)   ! (nnz)
    integer, allocatable          :: colind(:)   ! (nnz)
    integer, allocatable          :: rowptr(:)   ! (nrows+1)
  end type csr_matrix

contains

  ! ==================================================================
  ! CSR basic utilities (Phase 0b)
  ! ==================================================================

  ! ------------------------------------------------------------------
  ! Initialize an empty CSR matrix with the given dimensions.
  ! Allocates rowptr(nrows+1) set to 1 (no entries), and zero-length
  ! values/colind arrays.
  ! ------------------------------------------------------------------
  subroutine csr_init(mat, nrows, ncols)
    type(csr_matrix), intent(out) :: mat
    integer, intent(in) :: nrows, ncols

    mat%nrows = nrows
    mat%ncols = ncols
    mat%nnz   = 0

    allocate(mat%rowptr(nrows + 1))
    mat%rowptr = 1   ! empty: all rows start at index 1

    allocate(mat%values(0))
    allocate(mat%colind(0))
  end subroutine csr_init

  ! ------------------------------------------------------------------
  ! Deallocate all arrays in a CSR matrix.
  ! ------------------------------------------------------------------
  subroutine csr_free(mat)
    type(csr_matrix), intent(inout) :: mat

    if (allocated(mat%values))  deallocate(mat%values)
    if (allocated(mat%colind))  deallocate(mat%colind)
    if (allocated(mat%rowptr))  deallocate(mat%rowptr)
    mat%nrows = 0
    mat%ncols = 0
    mat%nnz   = 0
  end subroutine csr_free

  ! ------------------------------------------------------------------
  ! Build a CSR matrix from COO (coordinate) triplets.
  !
  ! Input:
  !   nrows, ncols  -- matrix dimensions
  !   nnz_in        -- number of COO triplets supplied
  !   rows(1:nnz_in), cols(1:nnz_in) -- 1-based row/col indices
  !   vals(1:nnz_in)                  -- complex values
  !
  ! Output:
  !   mat -- CSR matrix with entries sorted by (row, col) and
  !          duplicate entries merged by summation.
  !
  ! Uses bottom-up merge sort O(N log N) to sort COO triplets by
  ! (row, col), then merges duplicates by summation.
  ! ------------------------------------------------------------------
  subroutine csr_build_from_coo(mat, nrows, ncols, nnz_in, rows, cols, vals)
    type(csr_matrix), intent(out) :: mat
    integer, intent(in) :: nrows, ncols, nnz_in
    integer, intent(in) :: rows(nnz_in), cols(nnz_in)
    complex(kind=dp), intent(in) :: vals(nnz_in)

    integer, allocatable :: idx(:), r_sorted(:), c_sorted(:)
    complex(kind=dp), allocatable :: v_sorted(:)
    integer :: i, nnz_final, row

    mat%nrows = nrows
    mat%ncols = ncols

    if (nnz_in == 0) then
      mat%nnz = 0
      allocate(mat%rowptr(nrows + 1))
      mat%rowptr = 1
      allocate(mat%values(0))
      allocate(mat%colind(0))
      return
    end if

      ! Build index array and sort by (row, col)
    allocate(idx(nnz_in))
    do i = 1, nnz_in
      idx(i) = i
    end do

    ! Merge sort of idx by (rows(idx), cols(idx))
    call merge_sort_coo(nnz_in, idx, rows, cols)

    ! Apply permutation and merge duplicates
    allocate(v_sorted(nnz_in))
    allocate(r_sorted(nnz_in))
    allocate(c_sorted(nnz_in))

    v_sorted(1) = vals(idx(1))
    r_sorted(1) = rows(idx(1))
    c_sorted(1) = cols(idx(1))
    nnz_final = 1

    do i = 2, nnz_in
      if (r_sorted(nnz_final) == rows(idx(i)) .and. &
          c_sorted(nnz_final) == cols(idx(i))) then
        ! Merge duplicate
        v_sorted(nnz_final) = v_sorted(nnz_final) + vals(idx(i))
      else
        nnz_final = nnz_final + 1
        r_sorted(nnz_final) = rows(idx(i))
        c_sorted(nnz_final) = cols(idx(i))
        v_sorted(nnz_final) = vals(idx(i))
      end if
    end do

    ! Build CSR arrays
    mat%nnz = nnz_final
    allocate(mat%values(nnz_final))
    allocate(mat%colind(nnz_final))
    allocate(mat%rowptr(nrows + 1))

    mat%values(1:nnz_final) = v_sorted(1:nnz_final)
    mat%colind(1:nnz_final) = c_sorted(1:nnz_final)

    ! Build rowptr: scan sorted row indices
    mat%rowptr = nnz_final + 1  ! sentinel
    do i = 1, nnz_final
      row = r_sorted(i)
      if (mat%rowptr(row) > i) mat%rowptr(row) = i
    end do
    mat%rowptr(nrows + 1) = nnz_final + 1

    ! Fill any empty rows: rowptr(row) should equal rowptr(row+1)
    do row = nrows, 1, -1
      if (mat%rowptr(row) == nnz_final + 1) then
        mat%rowptr(row) = mat%rowptr(row + 1)
      end if
    end do

    deallocate(idx, v_sorted, r_sorted, c_sorted)
  end subroutine csr_build_from_coo

  ! ------------------------------------------------------------------
  ! Build a CSR matrix from COO triplets AND save the COO-to-CSR
  ! mapping for subsequent fast value updates.
  !
  ! Same as csr_build_from_coo but additionally produces coo_to_csr
  ! array: for each input COO index i (1..nnz_in), coo_to_csr(i) gives
  ! the final CSR position after sort and duplicate merging.  Multiple
  ! COO entries that merge into the same CSR position share the same
  ! mapping value.
  !
  ! Use the saved mapping with csr_set_values_from_coo for O(NNZ) value
  ! updates on subsequent k-points (no sort needed).
  ! ------------------------------------------------------------------
  subroutine csr_build_from_coo_cached(mat, nrows, ncols, nnz_in, rows, cols, &
      vals, coo_to_csr)
    type(csr_matrix), intent(out) :: mat
    integer, intent(in) :: nrows, ncols, nnz_in
    integer, intent(in) :: rows(nnz_in), cols(nnz_in)
    complex(kind=dp), intent(in) :: vals(nnz_in)
    integer, intent(out), allocatable :: coo_to_csr(:)

    integer, allocatable :: idx(:), r_sorted(:), c_sorted(:)
    complex(kind=dp), allocatable :: v_sorted(:)
    integer :: i, nnz_final, row

    mat%nrows = nrows
    mat%ncols = ncols

    if (nnz_in == 0) then
      mat%nnz = 0
      allocate(mat%rowptr(nrows + 1))
      mat%rowptr = 1
      allocate(mat%values(0))
      allocate(mat%colind(0))
      allocate(coo_to_csr(0))
      return
    end if

    ! Build index array and sort by (row, col)
    allocate(idx(nnz_in))
    do i = 1, nnz_in
      idx(i) = i
    end do
    call merge_sort_coo(nnz_in, idx, rows, cols)

    ! Apply permutation and merge duplicates, tracking the mapping
    allocate(v_sorted(nnz_in))
    allocate(r_sorted(nnz_in))
    allocate(c_sorted(nnz_in))
    allocate(coo_to_csr(nnz_in))

    v_sorted(1) = vals(idx(1))
    r_sorted(1) = rows(idx(1))
    c_sorted(1) = cols(idx(1))
    coo_to_csr(idx(1)) = 1
    nnz_final = 1

    do i = 2, nnz_in
      if (r_sorted(nnz_final) == rows(idx(i)) .and. &
          c_sorted(nnz_final) == cols(idx(i))) then
        ! Merge duplicate -- map to same CSR position
        v_sorted(nnz_final) = v_sorted(nnz_final) + vals(idx(i))
        coo_to_csr(idx(i)) = nnz_final
      else
        nnz_final = nnz_final + 1
        r_sorted(nnz_final) = rows(idx(i))
        c_sorted(nnz_final) = cols(idx(i))
        v_sorted(nnz_final) = vals(idx(i))
        coo_to_csr(idx(i)) = nnz_final
      end if
    end do

    ! Build CSR arrays
    mat%nnz = nnz_final
    allocate(mat%values(nnz_final))
    allocate(mat%colind(nnz_final))
    allocate(mat%rowptr(nrows + 1))

    mat%values(1:nnz_final) = v_sorted(1:nnz_final)
    mat%colind(1:nnz_final) = c_sorted(1:nnz_final)

    ! Build rowptr
    mat%rowptr = nnz_final + 1  ! sentinel
    do i = 1, nnz_final
      row = r_sorted(i)
      if (mat%rowptr(row) > i) mat%rowptr(row) = i
    end do
    mat%rowptr(nrows + 1) = nnz_final + 1

    ! Fill any empty rows
    do row = nrows, 1, -1
      if (mat%rowptr(row) == nnz_final + 1) then
        mat%rowptr(row) = mat%rowptr(row + 1)
      end if
    end do

    deallocate(idx, v_sorted, r_sorted, c_sorted)
  end subroutine csr_build_from_coo_cached

  ! ------------------------------------------------------------------
  ! Fast value update of an existing CSR matrix using a cached mapping.
  !
  ! Given new COO values and the coo_to_csr mapping saved by
  ! csr_build_from_coo_cached, this routine accumulates values into
  ! mat%values in O(NNZ_in) time with no sorting.
  !
  ! The CSR structure (rowptr, colind) is NOT modified.  The caller
  ! must ensure the same COO sparsity pattern as the original build.
  ! ------------------------------------------------------------------
  subroutine csr_set_values_from_coo(mat, nnz_in, coo_to_csr, vals)
    type(csr_matrix), intent(inout) :: mat
    integer, intent(in) :: nnz_in
    integer, intent(in) :: coo_to_csr(nnz_in)
    complex(kind=dp), intent(in) :: vals(nnz_in)

    integer :: i, csr_pos

    ! Zero out all values first (they will be re-accumulated)
    mat%values = cmplx(0.0_dp, 0.0_dp, kind=dp)

    ! Accumulate each COO entry into its mapped CSR position
    do i = 1, nnz_in
      csr_pos = coo_to_csr(i)
      if (csr_pos >= 1 .and. csr_pos <= mat%nnz) then
        mat%values(csr_pos) = mat%values(csr_pos) + vals(i)
      end if
    end do
  end subroutine csr_set_values_from_coo

  ! ==================================================================
  ! Private sort helper
  ! ==================================================================

  ! ------------------------------------------------------------------
  ! Bottom-up merge sort of index array idx by key (rows(idx), cols(idx)).
  !
  ! Stable sort: equal keys preserve their original relative order,
  ! which is required for correct duplicate merging afterwards.
  !
  ! Uses a workspace array of the same size as idx for the merge step.
  ! Complexity: O(N log N) where N = nnz.
  ! ------------------------------------------------------------------
  subroutine merge_sort_coo(nnz, idx, rows, cols)
    integer, intent(in)    :: nnz
    integer, intent(inout) :: idx(nnz)
    integer, intent(in)    :: rows(nnz)
    integer, intent(in)    :: cols(nnz)

    integer, allocatable :: work(:)
    integer :: width, start, mid, finish, i, j, k

    if (nnz <= 1) return

    allocate(work(nnz))

    ! Bottom-up merge sort: merge sub-arrays of doubling width
    width = 1
    do while (width < nnz)
      ! Merge pairs of adjacent runs [start:mid] and [mid+1:finish]
      start = 1
      do while (start <= nnz)
        mid = min(start + width - 1, nnz)
        finish = min(start + 2*width - 1, nnz)

        if (mid < finish) then
          ! Merge idx(start:mid) and idx(mid+1:finish) into work
          i = start
          j = mid + 1
          k = start
          do while (i <= mid .and. j <= finish)
            if (rows(idx(i)) < rows(idx(j))) then
              work(k) = idx(i); i = i + 1
            else if (rows(idx(i)) > rows(idx(j))) then
              work(k) = idx(j); j = j + 1
            else
              ! Same row -- compare columns
              if (cols(idx(i)) <= cols(idx(j))) then
                work(k) = idx(i); i = i + 1
              else
                work(k) = idx(j); j = j + 1
              end if
            end if
            k = k + 1
          end do

          ! Copy remaining elements from left run
          do while (i <= mid)
            work(k) = idx(i); i = i + 1; k = k + 1
          end do

          ! Copy remaining elements from right run
          do while (j <= finish)
            work(k) = idx(j); j = j + 1; k = k + 1
          end do

          ! Copy back from work to idx
          idx(start:finish) = work(start:finish)
        end if

        start = start + 2*width
      end do

      width = width * 2
    end do

    deallocate(work)
  end subroutine merge_sort_coo

  ! ==================================================================
  ! Kronecker product routines (Phase 0c)
  ! ==================================================================

  ! ------------------------------------------------------------------
  ! Kronecker product of two dense complex matrices, result in CSR.
  !
  ! C = A (x) B   where A is (na1 x na2) and B is (nb1 x nb2).
  !
  ! Result C is (na1*nb1 x na2*nb2).  Only entries with |value| >=
  ! tol are stored.  Memory is allocated in a single pass after
  ! counting nonzeros.
  !
  ! Algorithm: for each (i,j) in A and (k,l) in B, the output entry
  !   C((i-1)*nb1+k, (j-1)*nb2+l) = A(i,j) * B(k,l)
  ! is stored if |A(i,j) * B(k,l)| >= tol.
  ! ------------------------------------------------------------------
  subroutine kron_dense_dense(A, na1, na2, B, nb1, nb2, C, tol)
    complex(kind=dp), intent(in) :: A(na1, na2)
    integer, intent(in) :: na1, na2
    complex(kind=dp), intent(in) :: B(nb1, nb2)
    integer, intent(in) :: nb1, nb2
    type(csr_matrix), intent(out) :: C
    real(kind=dp), intent(in), optional :: tol

    real(kind=dp) :: threshold
    integer :: nc1, nc2, i, j, k, l, nnz_est, idx_c
    integer :: row_c, col_c
    complex(kind=dp) :: val_prod
    real(kind=dp) :: abs_prod

    integer, allocatable :: rows_coo(:), cols_coo(:)
    complex(kind=dp), allocatable :: vals_coo(:)

    threshold = 1.0e-15_dp
    if (present(tol)) threshold = tol

    nc1 = na1 * nb1
    nc2 = na2 * nb2

    ! First pass: count nonzeros to allocate exactly
    nnz_est = 0
    do j = 1, na2
      do i = 1, na1
        if (abs(A(i, j)) < threshold) cycle
        do l = 1, nb2
          do k = 1, nb1
            if (abs(B(k, l)) < threshold) cycle
            abs_prod = abs(A(i, j)) * abs(B(k, l))
            if (abs_prod >= threshold) nnz_est = nnz_est + 1
          end do
        end do
      end do
    end do

    if (nnz_est == 0) then
      call csr_init(C, nc1, nc2)
      return
    end if

    ! Allocate COO arrays
    allocate(rows_coo(nnz_est))
    allocate(cols_coo(nnz_est))
    allocate(vals_coo(nnz_est))

    ! Second pass: fill COO triplets
    idx_c = 0
    do j = 1, na2
      do i = 1, na1
        if (abs(A(i, j)) < threshold) cycle
        do l = 1, nb2
          do k = 1, nb1
            if (abs(B(k, l)) < threshold) cycle
            val_prod = A(i, j) * B(k, l)
            if (abs(val_prod) < threshold) cycle
            idx_c = idx_c + 1
            row_c = (i - 1) * nb1 + k
            col_c = (j - 1) * nb2 + l
            rows_coo(idx_c) = row_c
            cols_coo(idx_c) = col_c
            vals_coo(idx_c) = val_prod
          end do
        end do
      end do
    end do

    ! Build CSR from COO (sorts and merges duplicates)
    call csr_build_from_coo(C, nc1, nc2, idx_c, rows_coo(1:idx_c), &
                            cols_coo(1:idx_c), vals_coo(1:idx_c))

    deallocate(rows_coo, cols_coo, vals_coo)
  end subroutine kron_dense_dense

  ! ------------------------------------------------------------------
  ! Kronecker product A (x) I_n  where I_n is the n x n identity.
  !
  ! A is a dense (na x na) matrix.  Result C is (na*n x na*n) in CSR.
  !
  ! This is efficient: only nnz(A) blocks of size (n x n) are produced,
  ! but since I_n is diagonal, each block contributes exactly one entry
  ! per A nonzero per identity row.  Total nnz = nnz(A) * n_eye
  ! (entries above threshold).
  !
  ! C((i-1)*n+k, (j-1)*n+k) = A(i,j)  for k = 1..n
  ! ------------------------------------------------------------------
  subroutine kron_dense_eye(A, na, n_eye, C, tol)
    complex(kind=dp), intent(in) :: A(na, na)
    integer, intent(in) :: na, n_eye
    type(csr_matrix), intent(out) :: C
    real(kind=dp), intent(in), optional :: tol

    real(kind=dp) :: threshold
    integer :: nc, i, j, k, idx_c
    integer :: row_c, col_c

    integer, allocatable :: rows_coo(:), cols_coo(:)
    complex(kind=dp), allocatable :: vals_coo(:)

    threshold = 1.0e-15_dp
    if (present(tol)) threshold = tol

    nc = na * n_eye

    ! Count nonzeros
    idx_c = 0
    do j = 1, na
      do i = 1, na
        if (abs(A(i, j)) >= threshold) then
          idx_c = idx_c + n_eye
        end if
      end do
    end do

    if (idx_c == 0) then
      call csr_init(C, nc, nc)
      return
    end if

    ! Fill COO
    allocate(rows_coo(idx_c))
    allocate(cols_coo(idx_c))
    allocate(vals_coo(idx_c))

    idx_c = 0
    do j = 1, na
      do i = 1, na
        if (abs(A(i, j)) < threshold) cycle
        do k = 1, n_eye
          idx_c = idx_c + 1
          row_c = (i - 1) * n_eye + k
          col_c = (j - 1) * n_eye + k
          rows_coo(idx_c) = row_c
          cols_coo(idx_c) = col_c
          vals_coo(idx_c) = A(i, j)
        end do
      end do
    end do

    call csr_build_from_coo(C, nc, nc, idx_c, rows_coo, cols_coo, vals_coo)

    deallocate(rows_coo, cols_coo, vals_coo)
  end subroutine kron_dense_eye

  ! ------------------------------------------------------------------
  ! Kronecker product I_n (x) B  where I_n is the n x n identity.
  !
  ! B is a dense (nb1 x nb2) matrix.  Result C is (n*nb1 x n*nb2) in CSR.
  !
  ! This places a copy of B at each diagonal block:
  !   C((s-1)*nb1+i, (s-1)*nb2+j) = B(i,j)  for s = 1..n
  !
  ! Total nnz = n * nnz(B) (entries above threshold).
  ! ------------------------------------------------------------------
  subroutine kron_eye_dense(n_eye, B, nb1, nb2, C, tol)
    integer, intent(in) :: n_eye
    complex(kind=dp), intent(in) :: B(nb1, nb2)
    integer, intent(in) :: nb1, nb2
    type(csr_matrix), intent(out) :: C
    real(kind=dp), intent(in), optional :: tol

    real(kind=dp) :: threshold
    integer :: nc1, nc2, i, j, s, idx_c, nnz_B
    integer :: row_c, col_c

    integer, allocatable :: rows_coo(:), cols_coo(:)
    complex(kind=dp), allocatable :: vals_coo(:)

    threshold = 1.0e-15_dp
    if (present(tol)) threshold = tol

    nc1 = n_eye * nb1
    nc2 = n_eye * nb2

    ! Count nonzeros in B first
    nnz_B = 0
    do j = 1, nb2
      do i = 1, nb1
        if (abs(B(i, j)) >= threshold) nnz_B = nnz_B + 1
      end do
    end do

    idx_c = n_eye * nnz_B

    if (idx_c == 0) then
      call csr_init(C, nc1, nc2)
      return
    end if

    ! Fill COO: for each diagonal block s, copy nonzeros of B
    allocate(rows_coo(idx_c))
    allocate(cols_coo(idx_c))
    allocate(vals_coo(idx_c))

    idx_c = 0
    do s = 1, n_eye
      do j = 1, nb2
        do i = 1, nb1
          if (abs(B(i, j)) < threshold) cycle
          idx_c = idx_c + 1
          row_c = (s - 1) * nb1 + i
          col_c = (s - 1) * nb2 + j
          rows_coo(idx_c) = row_c
          cols_coo(idx_c) = col_c
          vals_coo(idx_c) = B(i, j)
        end do
      end do
    end do

    call csr_build_from_coo(C, nc1, nc2, idx_c, rows_coo, cols_coo, vals_coo)

    deallocate(rows_coo, cols_coo, vals_coo)
  end subroutine kron_eye_dense

  ! ------------------------------------------------------------------
  ! Kronecker product of two 1D dense arrays (treated as vectors)
  ! producing a 2D CSR matrix.
  !
  ! A is a complex array of length na, B is a complex array of length nb.
  ! The result C is (na x nb) where C(i, j) = A(i) * B(j), stored in
  ! CSR format.  Used for building outer products from 1D FD operators.
  !
  ! This is the outer product: C = A * B^T  (rank-1 Kronecker).
  ! ------------------------------------------------------------------
  subroutine kron_dense_dense_1d(A, na, B, nb, C, tol)
    complex(kind=dp), intent(in) :: A(na)
    integer, intent(in) :: na
    complex(kind=dp), intent(in) :: B(nb)
    integer, intent(in) :: nb
    type(csr_matrix), intent(out) :: C
    real(kind=dp), intent(in), optional :: tol

    real(kind=dp) :: threshold
    integer :: i, j, idx_c
    complex(kind=dp) :: val_prod

    integer, allocatable :: rows_coo(:), cols_coo(:)
    complex(kind=dp), allocatable :: vals_coo(:)

    threshold = 1.0e-15_dp
    if (present(tol)) threshold = tol

    ! Count nonzeros
    idx_c = 0
    do j = 1, nb
      if (abs(B(j)) < threshold) cycle
      do i = 1, na
        if (abs(A(i)) < threshold) cycle
        val_prod = A(i) * B(j)
        if (abs(val_prod) >= threshold) idx_c = idx_c + 1
      end do
    end do

    if (idx_c == 0) then
      call csr_init(C, na, nb)
      return
    end if

    ! Fill COO
    allocate(rows_coo(idx_c))
    allocate(cols_coo(idx_c))
    allocate(vals_coo(idx_c))

    idx_c = 0
    do j = 1, nb
      if (abs(B(j)) < threshold) cycle
      do i = 1, na
        if (abs(A(i)) < threshold) cycle
        val_prod = A(i) * B(j)
        if (abs(val_prod) < threshold) cycle
        idx_c = idx_c + 1
        rows_coo(idx_c) = i
        cols_coo(idx_c) = j
        vals_coo(idx_c) = val_prod
      end do
    end do

    call csr_build_from_coo(C, na, nb, idx_c, rows_coo(1:idx_c), &
                            cols_coo(1:idx_c), vals_coo(1:idx_c))

    deallocate(rows_coo, cols_coo, vals_coo)
  end subroutine kron_dense_dense_1d

  ! ==================================================================
  ! CSR arithmetic routines (Phase 1b)
  ! ==================================================================

  ! ------------------------------------------------------------------
  ! Add two CSR matrices: C = alpha*A + beta*B.
  ! A and B must have the same dimensions.  C is newly allocated.
  ! Uses COO accumulation: sparsifies union of sparsity patterns.
  ! ------------------------------------------------------------------
  subroutine csr_add(A, B, C, alpha, beta)
    type(csr_matrix), intent(in)  :: A, B
    type(csr_matrix), intent(out) :: C
    complex(kind=dp), intent(in), optional :: alpha, beta
    complex(kind=dp) :: sa, sb

    integer :: nnz_est, idx_c, i, row, k
    integer, allocatable :: rows_coo(:), cols_coo(:)
    complex(kind=dp), allocatable :: vals_coo(:)

    ! Guard against self-aliasing: C must not be the same as A or B
    if (loc(C) == loc(A) .or. loc(C) == loc(B)) then
      print *, 'ERROR: csr_add self-aliasing detected'
      stop 1
    end if

    sa = cmplx(1.0_dp, 0.0_dp, kind=dp)
    sb = cmplx(1.0_dp, 0.0_dp, kind=dp)
    if (present(alpha)) sa = alpha
    if (present(beta))  sb = beta

    nnz_est = A%nnz + B%nnz

    if (nnz_est == 0) then
      call csr_init(C, A%nrows, A%ncols)
      return
    end if

    allocate(rows_coo(nnz_est))
    allocate(cols_coo(nnz_est))
    allocate(vals_coo(nnz_est))

    idx_c = 0

    ! Entries from A
    do row = 1, A%nrows
      do k = A%rowptr(row), A%rowptr(row + 1) - 1
        idx_c = idx_c + 1
        rows_coo(idx_c) = row
        cols_coo(idx_c) = A%colind(k)
        vals_coo(idx_c) = sa * A%values(k)
      end do
    end do

    ! Entries from B
    do row = 1, B%nrows
      do k = B%rowptr(row), B%rowptr(row + 1) - 1
        idx_c = idx_c + 1
        rows_coo(idx_c) = row
        cols_coo(idx_c) = B%colind(k)
        vals_coo(idx_c) = sb * B%values(k)
      end do
    end do

    call csr_build_from_coo(C, A%nrows, A%ncols, idx_c, &
      rows_coo(1:idx_c), cols_coo(1:idx_c), vals_coo(1:idx_c))

    deallocate(rows_coo, cols_coo, vals_coo)
  end subroutine csr_add

  ! ------------------------------------------------------------------
  ! Scale a CSR matrix by a complex scalar: result = scale * A.
  ! Copies A's structure and multiplies all values by scale.
  ! ------------------------------------------------------------------
  subroutine csr_scale(A, result_mat, scale)
    type(csr_matrix), intent(in)  :: A
    type(csr_matrix), intent(out) :: result_mat
    complex(kind=dp), intent(in)  :: scale

    result_mat%nrows = A%nrows
    result_mat%ncols = A%ncols
    result_mat%nnz   = A%nnz

    allocate(result_mat%rowptr(A%nrows + 1))
    allocate(result_mat%colind(A%nnz))
    allocate(result_mat%values(A%nnz))

    result_mat%rowptr = A%rowptr
    result_mat%colind = A%colind
    result_mat%values = A%values * scale
  end subroutine csr_scale

  ! ------------------------------------------------------------------
  ! Apply variable coefficient (element-wise row scaling) to a CSR
  ! matrix: result(i,j) = -profile(i) * A(i,j).
  !
  ! The profile array has length A%nrows.  Each row i of A is scaled
  ! by -profile(i).  The negative sign follows the convention in the
  ! 1D applyVariableCoeff (kpterms = -profile * FD).
  !
  ! Optionally, cut-cell face fractions can modify boundary stencil
  ! weights.  When face_frac_y and face_frac_z are present, each
  ! nonzero in row i is additionally scaled by the geometric mean of
  ! the face fractions for the corresponding direction:
  !   d^2/dx^2 terms: multiply by (face_frac_x(i,1) + face_frac_x(i,2)) / 2
  !   d^2/dy^2 terms: multiply by (face_frac_y(i,1) + face_frac_y(i,2)) / 2
  !   d/dx terms:     multiply by (face_frac_x(i,2) - face_frac_x(i,1)) / 2
  !   d/dy terms:     multiply by (face_frac_y(i,2) - face_frac_y(i,1)) / 2
  !
  ! For now (Phase 1b), the face fractions are applied as a simple
  ! row volume scaling: result(i,j) = -profile(i) * vol(i) * A(i,j),
  ! where vol(i) = cell_volume(i).  This is the correct box-integration
  ! weight for a 2D Laplacian with cut cells.
  ! ------------------------------------------------------------------
  subroutine csr_apply_variable_coeff(A, profile, result_mat, cell_volume)
    type(csr_matrix), intent(in)  :: A
    real(kind=dp), intent(in)     :: profile(:)     ! (A%nrows)
    type(csr_matrix), intent(out) :: result_mat
    real(kind=dp), intent(in), optional :: cell_volume(:)  ! (A%nrows)

    integer :: i, k, row, nnz_out
    real(kind=dp) :: scale
    integer, allocatable :: rows_coo(:), cols_coo(:)
    complex(kind=dp), allocatable :: vals_coo(:)

    nnz_out = A%nnz

    if (nnz_out == 0) then
      call csr_init(result_mat, A%nrows, A%ncols)
      return
    end if

    allocate(rows_coo(nnz_out))
    allocate(cols_coo(nnz_out))
    allocate(vals_coo(nnz_out))

    ! Build COO from A with row scaling
    k = 0
    do row = 1, A%nrows
      scale = -profile(row)
      if (present(cell_volume)) then
        if (cell_volume(row) > 0.0_dp) then
          scale = scale * cell_volume(row)
        else
          scale = 0.0_dp
        end if
      end if
      do i = A%rowptr(row), A%rowptr(row + 1) - 1
        k = k + 1
        rows_coo(k) = row
        cols_coo(k) = A%colind(i)
        vals_coo(k) = scale * A%values(i)
      end do
    end do

    call csr_build_from_coo(result_mat, A%nrows, A%ncols, k, &
      rows_coo(1:k), cols_coo(1:k), vals_coo(1:k))

    deallocate(rows_coo, cols_coo, vals_coo)
  end subroutine csr_apply_variable_coeff

end module sparse_matrices
