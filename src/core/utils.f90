module utils

use definitions
use mkl_spblas

implicit none

contains

  subroutine dnscsr_z_mkl(nzmax, N, dns, csr)


    type(sparse_matrix_T), intent(out) :: csr
    integer, intent(inout) :: nzmax
    integer, intent(in) :: N
    complex(kind=dp), intent(in), dimension(:,:) :: dns

    integer :: ierr, info, next, i, j, nnz_count
    complex(kind=dp), allocatable, dimension(:) :: Aa_a, A_a
    integer, allocatable, dimension(:) :: Aa_ia, Aa_ja, A_ia, A_ja
    type(sparse_matrix_T) :: HT_coo


    allocate(Aa_a(nzmax))
    allocate(Aa_ia(nzmax))
    allocate(Aa_ja(nzmax))

    next = 1
    do i = 1, N
      do j = 1, N
        call insertCOO_cmplx(Aa_a, Aa_ia, Aa_ja, dns(i,j), i, j, next, nzmax)
      end do
    end do

    ! Finalize: sort by (row, col) and merge duplicates
    if (next > 1) then
      nnz_count = next - 1
      call finalizeCOO_cmplx(Aa_a, Aa_ia, Aa_ja, nnz_count)
    else
      nnz_count = 0
    end if

    allocate(A_a(0))
    allocate(A_ia(0))
    allocate(A_ja(0))

    A_a = [ Aa_a(1:nnz_count) ]
    A_ia = [ Aa_ia(1:nnz_count) ]
    A_ja = [ Aa_ja(1:nnz_count) ]

    deallocate(Aa_a, Aa_ia, Aa_ja)

    nzmax = ubound(A_a, 1)

    info = mkl_sparse_z_create_coo (HT_coo, SPARSE_INDEX_BASE_ONE, N, &
    & N, ubound(A_a, dim=1), A_ia, A_ja, A_a)

    if (info /= 0) stop 'error creating HT_coo'

    info = mkl_sparse_convert_csr (HT_coo, SPARSE_OPERATION_NON_TRANSPOSE, csr)

    if (info /= 0) stop 'error converting HT_coo to HT_csr'

    deallocate(A_a, A_ia, A_ja)

  end subroutine dnscsr_z_mkl

  subroutine insertCOO_cmplx(v, r, c, vvalue, i, j, next, nnz)

      complex ( kind = dp ), intent(inout) :: v(:)
      integer ( kind = 4 ), intent(inout) :: r(:), c(:), next
      complex ( kind = dp ), intent(in) :: vvalue
      integer ( kind = 4 ), intent(in) :: i, j, nnz

      if ( vvalue /= cmplx(0.0_dp, 0.0_dp, kind=dp) .and. abs(vvalue) >= 10E-10 ) then

          if (nnz < next ) THEN
            print *, nnz, next
            stop 'next greater than nnz'
          endif

          v(next) = vvalue
          r(next) = i
          c(next) = j
          next = next + 1

      end if

  end subroutine insertCOO_cmplx

  subroutine finalizeCOO_cmplx(v, r, c, nnz)
      complex ( kind = dp ), intent(inout), allocatable :: v(:)
      integer ( kind = 4 ), intent(inout), allocatable :: r(:), c(:)
      integer ( kind = 4 ), intent(inout) :: nnz

      integer ( kind = 4 ), allocatable :: idx(:)
      integer ( kind = 4 ) :: i, nnz_final
      complex ( kind = dp ), allocatable :: v_sorted(:)
      integer ( kind = 4 ), allocatable :: r_sorted(:), c_sorted(:)

      ! Nothing to do for empty or single-element arrays
      if (nnz <= 1) return

      ! Build index array and sort by (row, col)
      allocate(idx(nnz))
      do i = 1, nnz
        idx(i) = i
      end do

      ! Index sort: sort idx by (r, c) pairs using insertion sort
      call index_sort_by_rc(idx, r, c, nnz)

      ! Apply permutation into temporary arrays
      allocate(v_sorted(nnz))
      allocate(r_sorted(nnz))
      allocate(c_sorted(nnz))
      do i = 1, nnz
        v_sorted(i) = v(idx(i))
        r_sorted(i) = r(idx(i))
        c_sorted(i) = c(idx(i))
      end do

      ! Merge duplicates in-place in the sorted arrays
      nnz_final = 1
      do i = 2, nnz
        if (r_sorted(nnz_final) == r_sorted(i) .and. &
          & c_sorted(nnz_final) == c_sorted(i)) then
          ! Merge: add value to previous unique entry
          v_sorted(nnz_final) = v_sorted(nnz_final) + v_sorted(i)
        else
          ! New unique entry
          nnz_final = nnz_final + 1
          v_sorted(nnz_final) = v_sorted(i)
          r_sorted(nnz_final) = r_sorted(i)
          c_sorted(nnz_final) = c_sorted(i)
        end if
      end do

      ! Copy merged results back
      v(1:nnz_final) = v_sorted(1:nnz_final)
      r(1:nnz_final) = r_sorted(1:nnz_final)
      c(1:nnz_final) = c_sorted(1:nnz_final)
      nnz = nnz_final

      deallocate(idx, v_sorted, r_sorted, c_sorted)

  end subroutine finalizeCOO_cmplx

  subroutine index_sort_by_rc(idx, r, c, n)
      integer ( kind = 4 ), intent(inout) :: idx(:)
      integer ( kind = 4 ), intent(in) :: r(:), c(:)
      integer ( kind = 4 ), intent(in) :: n

      integer ( kind = 4 ) :: i, j, key_idx
      integer ( kind = 4 ) :: key_r, key_c

      ! Insertion sort of idx by (r(idx), c(idx)) pairs
      do i = 2, n
        key_idx = idx(i)
        key_r = r(key_idx)
        key_c = c(key_idx)
        j = i - 1
        do while (j >= 1)
          if (r(idx(j)) < key_r) exit
          if (r(idx(j)) == key_r .and. c(idx(j)) <= key_c) exit
          idx(j+1) = idx(j)
          j = j - 1
        end do
        idx(j+1) = key_idx
      end do

  end subroutine index_sort_by_rc


complex(kind=dp) function simpson(f,a,b)

  complex(kind=dp), intent(in) :: f(:)
  real(kind=dp), intent(in) :: a, b

  real(kind=dp) :: h
  integer :: num, i, j, N
  real(kind=dp), allocatable :: sc(:)

  if (mod(size(f), 2) == 0) then
    print *, 'Error: Simpson integration requires odd number of points.'
    stop
  end if

  num = size(f)
  h = (b-a)/(num-1)
  h = h/3.0_dp

  N = num
  allocate(sc(num))

  sc = 1.0_dp
  forall(i=2:N-1:2) sc(i) = 2.0_dp
  forall(i=3:N-1:2) sc(i) = 4.0_dp

  simpson= h*sum(sc*f)

  deallocate(sc)

end function simpson


end module utils
