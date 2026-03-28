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

    integer :: ierr, info, next, i, j
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
    ! print *, next-1

    allocate(A_a(0))
    allocate(A_ia(0))
    allocate(A_ja(0))

    A_a = [ Aa_a(1:next-1) ]
    A_ia = [ Aa_ia(1:next-1) ]
    A_ja = [ Aa_ja(1:next-1) ]

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

      integer ( kind = 4 ) :: exists, ii

      exists = 0
      if ( vvalue /= cmplx(0.0d0,0.0d0) .and. abs(vvalue) >= 10E-10 ) then
      !if ( vvalue /= 0.0d0) then
          !print *, vvalue
          if (nnz < next ) THEN
            print *, nnz, next
            stop 'next greater than nnz'
          endif

          if (next > 1) then
              do ii = 1, next-1
                  if (r(ii) == i .and. c(ii) == j .and. exists /= 1) then
                      v(ii) = v(ii) + vvalue
                      exists = 1
                  end if
              end do
              if (exists == 0) then
                  v(next) = vvalue
                  r(next) = i
                  c(next) = j
                  next = next + 1
                  exists = 0
              end if
          else
              v(next) = vvalue
              r(next) = i
              c(next) = j
              next = next + 1
              exists = 0
          end if

      end if

  end subroutine insertCOO_cmplx

  subroutine dnscsr ( nrow, ncol, nzmax, dns, ndns, a, ja, ia, ierr )

!*****************************************************************************80
!
!! DNSCSR converts Dense to Compressed Row Sparse format.
!
!  Discussion:
!
!    This routine converts a densely stored matrix into a row orientied
!    compactly sparse matrix.  It is the reverse of CSRDNS.
!
!    This routine does not check whether an element is small.  It considers
!    that A(I,J) is zero only if it is exactly equal to zero.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix.
!
!    Input, integer ( kind = 4 ) NZMAX, the maximum number of nonzero elements
!    allowed.  This should be set to be the lengths of the arrays A and JA.
!
!    Input, real DNS(NDNS,NCOL), an NROW by NCOL dense matrix.
!
!    Input, integer ( kind = 4 ) NDNS, the first dimension of DNS, which must be
!    at least NROW.
!
!    Output, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, integer ( kind = 4 ) IERR, error indicator.
!    0 means normal return;
!    I, means that the the code stopped while processing row I, because
!       there was no space left in A and JA, as defined by NZMAX.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) ndns
  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) dns(ndns,ncol)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nzmax

  ierr = 0
  next = 1
  ia(1) = 1

  do i = 1, nrow

    do j = 1, ncol

      if ( dns(i,j) /= 0.0D+00 ) then

        if ( nzmax < next ) then
          ierr = i
          return
        end if

        ja(next) = j
        a(next) = dns(i,j)
        next = next + 1

      end if

    end do

    ia(i+1) = next

  end do

  return
end

complex(kind=dp) function simpson(f,a,b)

  complex(kind=dp), intent(in) :: f(:)
  real(kind=dp), intent(in) :: a, b

  real(kind=dp) :: h
  integer :: num, i, j, N
  real(kind=dp), allocatable :: sc(:)

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
