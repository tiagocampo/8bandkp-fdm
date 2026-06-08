module utils

use definitions, only: dp
use sparse_matrices, only: csr_matrix, csr_build_from_coo

implicit none
private
public :: dnscsr_z_mkl, simpson, simpson_real
public :: insertCOO_cmplx, finalizeCOO_cmplx
public :: tick, tock
public :: get_unit, ensure_output_dir

contains

  subroutine dnscsr_z_mkl(nzmax, N, dns, csr)


    type(csr_matrix), intent(out) :: csr
    integer, intent(inout) :: nzmax
    integer, intent(in) :: N
    complex(kind=dp), intent(in), contiguous, dimension(:,:) :: dns

    integer :: next, i, j, nnz_count
    complex(kind=dp), allocatable, dimension(:) :: Aa_a, A_a
    integer, allocatable, dimension(:) :: Aa_ia, Aa_ja, A_ia, A_ja


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

    call csr_build_from_coo(csr, N, N, nzmax, A_ia, A_ja, A_a)

    deallocate(A_a, A_ia, A_ja)

  end subroutine dnscsr_z_mkl

  subroutine insertCOO_cmplx(v, r, c, vvalue, i, j, next, nnz)

      complex ( kind = dp ), intent(inout) :: v(:)
      integer, intent(inout) :: r(:), c(:), next
      complex ( kind = dp ), intent(in) :: vvalue
      integer, intent(in) :: i, j, nnz

      if ( vvalue /= cmplx(0.0_dp, 0.0_dp, kind=dp) .and. abs(vvalue) >= 1.0e-10_dp ) then

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
      integer, intent(inout), allocatable :: r(:), c(:)
      integer, intent(inout) :: nnz

      integer, allocatable :: idx(:)
      integer :: i, nnz_final
      complex ( kind = dp ), allocatable :: v_sorted(:)
      integer, allocatable :: r_sorted(:), c_sorted(:)

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
      integer, intent(inout) :: idx(:)
      integer, intent(in) :: r(:), c(:)
      integer, intent(in) :: n

      integer :: i, j, key_idx
      integer :: key_r, key_c

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
  integer :: num, i, N
  real(kind=dp), allocatable :: sc(:)

  if (size(f) < 3) then
    print *, 'Error: Simpson integration requires at least 3 points, got', size(f)
    stop
  end if

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
  ! Simpson's 1/3 rule: even-index interior points get 4, odd-index get 2
  do i = 2, N - 1, 2
    sc(i) = 4.0_dp
  end do
  do i = 3, N - 1, 2
    sc(i) = 2.0_dp
  end do

  simpson= h*sum(sc*f)

  deallocate(sc)

end function simpson


real(kind=dp) function simpson_real(f, a, b)
  ! Real-valued Simpson 1/3 composite rule (avoids complex arithmetic overhead).

  real(kind=dp), intent(in) :: f(:)
  real(kind=dp), intent(in) :: a, b

  real(kind=dp) :: h, weighted_sum
  integer :: num, i

  if (size(f) < 3) then
    print *, 'Error: simpson_real requires at least 3 points, got', size(f)
    stop
  end if

  if (mod(size(f), 2) == 0) then
    print *, 'Error: simpson_real requires odd number of points.'
    stop
  end if

  num = size(f)
  h = (b - a) / (num - 1) / 3.0_dp

  weighted_sum = f(1) + f(num)
  do i = 2, num - 1, 2
    weighted_sum = weighted_sum + 4.0_dp * f(i)
  end do
  do i = 3, num - 1, 2
    weighted_sum = weighted_sum + 2.0_dp * f(i)
  end do

  simpson_real = h * weighted_sum

end function simpson_real


  subroutine tick(t)
      integer, intent(out) :: t

      call system_clock(t)
  end subroutine tick

  ! returns time in seconds from now to time described by t
  real(kind=dp) function tock(t)
      integer, intent(in) :: t
      integer :: now, clock_rate

      call system_clock(now,clock_rate)

      tock = real(now - t, kind=dp) / real(clock_rate, kind=dp)
  end function tock

  ! ==================================================================
  ! Ensure the output directory exists (creates it if needed).
  ! ==================================================================
  subroutine ensure_output_dir()
    character(len=*), parameter :: OUTPUT_DIR = 'output'
    logical :: dir_exists

    inquire(file=OUTPUT_DIR//'/.', exist=dir_exists)
    if (.not. dir_exists) then
      call execute_command_line('mkdir -p '//OUTPUT_DIR)
    end if
  end subroutine ensure_output_dir

  ! ==================================================================
  ! GET_UNIT returns a free FORTRAN unit number.
  !
  ! A "free" FORTRAN unit number is a value between 1 and 99 which
  ! is not currently associated with an I/O device.
  !
  ! If IUNIT = 0, then no free FORTRAN unit could be found, although
  ! all 99 units were checked (except for units 5, 6 and 9, which
  ! are commonly reserved for console I/O).
  !
  ! Otherwise, IUNIT is a value between 1 and 99, representing a
  ! free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
  ! are special, and will never return those values.
  ! ==================================================================
  subroutine get_unit(iunit)
    integer(kind=4), intent(out) :: iunit

    integer(kind=4) :: i, ios
    logical(kind=4) :: lopen

    iunit = 0

    do i = 1, 99
      if (i /= 5 .and. i /= 6 .and. i /= 9) then
        inquire(unit=i, opened=lopen, iostat=ios)
        if (ios == 0) then
          if (.not. lopen) then
            iunit = i
            return
          end if
        end if
      end if
    end do
  end subroutine get_unit


end module utils
