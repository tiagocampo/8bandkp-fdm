module finitedifferences

  use definitions

  implicit none

  contains

    subroutine FDmatrixDense(nz, dz, order, stencil, boundary, matrix)

      integer, intent(in) :: nz, stencil, boundary, order
      real(kind=dp), intent(in) :: dz
      real(kind=dp), intent(inout), allocatable, dimension(:,:) :: matrix
      real(kind=dp), allocatable, dimension(:) :: v1


      allocate(v1(nz))
      v1 = 0

      call FDstencil(order, stencil, dz, v1)
      !boundary = 0 -> hard Wall
      !boundary = 1 -> periodic
      call toeplitz(v1, boundary, matrix)
      deallocate(v1)

    end subroutine FDmatrixDense

    subroutine FDstencil(order, stencil, d, vector, type)

      integer, intent(in) :: order, stencil
      real(kind=dp), intent(in) :: d
      real(kind=dp), intent(inout), allocatable, dimension(:) :: vector
      integer :: length, hf
      character(len = 1), optional :: type

      length = size(vector, dim=1)

      IF (ALLOCATED(vector)) DEALLOCATE(vector)
      ALLOCATE(vector(length))

      hf = int(length/2. + 0.5)

      if (order == 1) then

        select case(stencil)

        case(2)

          if (length .lt. 3) stop 'for stencil 2, n should be 3 or more'
          vector(hf-1) = -1.0_dp
          vector(hf+1) = 1.0_dp
          vector = vector/(2.0_dp*d)

        case(4)
          print *, 'case 4 not implemented yet'

        case(6)
          print *, 'case 6 not implemented yet'

        case(8)
          if (length .lt. 9) stop 'for stencil 8, n should be 9 or more'
          vector(hf-4) = 3.0_dp
          vector(hf-3) = -32.0_dp
          vector(hf-2) = 168.0_dp
          vector(hf-1) = -672.0_dp
          vector(hf) = 0.0_dp
          vector(hf+1) = 672.0_dp
          vector(hf+2) = -168.0_dp
          vector(hf+3) = 32.0_dp
          vector(hf+4) = -3.0_dp
          vector = vector/(840.0_dp*d)

        case default
          print *, 'not implemented yet'

        end select


      end if


      if (order == 2) then

        if (stencil == 2) then
          if (length .lt. 3) stop 'for stencil 2, n should be 3 or more'
          vector(hf-1) = 1
          vector(hf) = -2
          vector(hf+1) = 1
        endif

        if (stencil == 4) then
          if (length .lt. 5) stop 'for stencil 4, n should be 5 or more'
          vector(hf-2) = -1./12.
          vector(hf-1) = 4./3.
          vector(hf) = -5./2.
          vector(hf+1) = 4./3.
          vector(hf+2) = -1./12.
        endif

        if (stencil == 6) then
          if (length .lt. 7) stop 'for stencil 6, n should be 7 or more'
          vector(hf-3) = 1./90.
          vector(hf-2) = -3/20.
          vector(hf-1) = 3./2.
          vector(hf) = -49./18.
          vector(hf+1) = 3./2.
          vector(hf+2) = -3/20.
          vector(hf+3) = 1./90.
        endif

        if (stencil == 8) then
          if (length .lt. 9) stop 'for stencil 8, n should be 9 or more'
          vector(hf-4) = -1./560.
          vector(hf-3) = 8./315.
          vector(hf-2) = -1./5.
          vector(hf-1) = 8./5.
          vector(hf) = -205./72.
          vector(hf+1) = 8./5.
          vector(hf+2) = -1./5.
          vector(hf+3) = 8./315.
          vector(hf+4) = -1./560.
        endif

        if (stencil == 10) then
          if (length .lt. 11) stop 'for stencil 10, n should be 11 or more'
          vector(hf-5) = 8./25200.
          vector(hf-4) = -125./25200.
          vector(hf-3) = 1000./25200.
          vector(hf-2) = -6000./25200.
          vector(hf-1) = 42000./25200.
          vector(hf) = -73766./25200.
          vector(hf+1) = 42000./25200.
          vector(hf+2) = -6000./25200.
          vector(hf+3) = 1000./25200.
          vector(hf+4) = -125./25200.
          vector(hf-5) = 8./25200.
        endif

        vector = vector/d**2

      end if

    end subroutine FDstencil

    subroutine Identity(n,matrix,factor)

      integer, intent(in) :: n
      real(kind=dp), intent(inout), allocatable, dimension(:,:) :: matrix
      real(kind=dp), intent(in), optional :: factor
      integer :: i

      IF (ALLOCATED(matrix)) DEALLOCATE(matrix)
      ALLOCATE(matrix(n,n))

      matrix = 0
      if (present(factor)) then
        forall(i = 1:n) matrix(i,i) = 1*factor
      else
        forall(i = 1:n) matrix(i,i) = 1
      endif

    end subroutine

    subroutine toeplitz(vector, boundary, matrix)

      real(kind=dp), intent(in) :: vector(:)
      integer, intent(in) :: boundary
      real(kind=dp), intent(inout), allocatable :: matrix(:,:)
      integer :: length, idx, hf

      length = size(vector, dim=1)

      IF (ALLOCATED(matrix)) DEALLOCATE(matrix)
      ALLOCATE(matrix(length,length))

      hf = int(length/2. + 0.5)

      if (boundary == 0) then
        do idx = 1, length, 1
          matrix(:,idx) = eoshift(vector,shift=hf-idx)
        end do
      endif

      if (boundary == 1) then
        do idx = 1, length, 1
          matrix(:,idx) = cshift(vector,shift=hf-idx)
        end do
      end if

    end subroutine toeplitz


    ! from https://rosettacode.org/wiki/Kronecker_product#Fortran
    SUBROUTINE KPRODUCT(A,B,AB)	!AB = Kronecker product of A and B, both two-dimensional arrays.
  !Considers the arrays to be addressed as A(row,column), despite any storage order arrangements.        .
  !Creating array AB to fit here, adjusting the caller's array AB, may not work on some compilers.
     real(kind=dp) ::  A(:,:),B(:,:)		!Two-dimensional arrays, lower bound one.
     real(kind=dp), ALLOCATABLE :: AB(:,:)	!To be created to fit.
     INTEGER R,RA,RB,C,CA,CB,I,J	!Assistants.
      RA = UBOUND(A,DIM = 1)	!Ascertain the upper bounds of the incoming arrays.
      CA = UBOUND(A,DIM = 2)	!Their lower bounds will be deemed one,
      RB = UBOUND(B,DIM = 1)	!And the upper bound as reported will correspond.
      CB = UBOUND(B,DIM = 2)	!UBOUND(A) would give an array of two values, RA and CA, more for higher dimensionality.
      WRITE (6,1) "A",RA,CA,"B",RB,CB,"A.k.B",RA*RB,CA*CB	!Announce.
  1     FORMAT (3(A," is ",I0,"x",I0,1X))	!Three sets of sizes.
      !IF (ALLOCATED(AB)) DEALLOCATE(AB)	!Discard any lingering storage.
      !ALLOCATE (AB(RA*RB,CA*CB))		!Obtain the exact desired size.
      R = 0		!Syncopation: start the row offset.

      DO I = 1,RA	!Step down the rows of A.
        C = 0		!For each row, start the column offset.
        DO J = 1,CA		!Step along the columns of A.
          AB(R + 1:R + RB,C + 1:C + CB) = AB(R + 1:R + RB,C + 1:C + CB) + A(I,J)*B	!Place a block of B values.
          C = C + CB		!Advance a block of columns.
        END DO		!On to the next column of A.
        R = R + RB		!Advance a block of rows.

      END DO	!On to the next row of A.

    END SUBROUTINE KPRODUCT	!No tests for bad parameters, or lack of storage...

end module finitedifferences
