module outputFunctions

  use definitions
  use utils

  implicit NONE

  CONTAINS


    subroutine writeEigenfunctions(N, evnum, A, k, fdstep, z, is_bulk)

      integer, intent(in) :: N, evnum, k, fdstep
      real(kind=dp), intent(in) :: z(:)
      complex(kind=dp), intent(in) :: A(N, evnum)
      logical, intent(in) :: is_bulk

      integer (kind=4) :: iounit, iounit2, ios
      character (len = 255) :: filename, filename2, x1, x2
      integer :: i, j, ncomp

      real(kind=dp), allocatable :: R(:,:)
      complex(kind=dp), allocatable :: C(:,:)
      real(kind=dp) :: lim
      real(kind=dp) :: parts(evnum,8)

      character(len=8) :: fmt ! format descriptor
      fmt = '(I5.5)' ! an integer of width 5 with zeros at the left

      write (x1,fmt) k ! converting integer to string using a 'internal file'

      if (is_bulk) then
        ! For bulk, we have 8 components per state
        ncomp = 8
        allocate(C(evnum, ncomp))
        ! For bulk, preserve complex components
        do i = 1, evnum
          do j = 1, ncomp
            C(i,j) = A(j,i)
          end do
        end do
      else
        ! For quantum well, we have fdstep points for each of the 8 components
        ncomp = N
        allocate(R(evnum, ncomp))
        ! For quantum well, calculate probability density
        R = transpose(real(conjg(A)*A))
      end if

      if (.not. is_bulk) then
        LIM = 100.0 * MAXVAL(ABS(SUM(R(:,:),2))) * EPSILON(ABS(R(1,1)))
        DO I=1,SIZE(R,2)
           WHERE(ABS(R(:,I)) < LIM) R(:,I)=0.0_dp
        ENDDO
      end if

      do j = 1, evnum
        write (x2,fmt) j ! converting integer to string using a 'internal file'
        filename = 'eigenfunctions_k_'//trim(x1)//'_ev_'//trim(x2)//'.dat'
        call get_unit(iounit)

        open(unit=iounit, file=filename, iostat=ios, status="replace", action="write")
        if ( ios /= 0 ) stop "Error opening file "

        if (is_bulk) then
          ! For bulk, write real and imaginary parts of each component
          write(unit=iounit, fmt="(2(g14.6))",iostat=ios) &
          & (real(C(j,i)), aimag(C(j,i)), i=1,8)
        else
          ! For quantum well, write z and all components
          write(unit=iounit, fmt="(9(g14.6))",iostat=ios) &
          & (z(i), R(j,i), R(j,fdstep+i), &
          & R(j,2*fdstep+i), R(j,3*fdstep+i), &
          & R(j,4*fdstep+i), R(j,5*fdstep+i), &
          & R(j,6*fdstep+i), R(j,7*fdstep+i), i=1,fdstep)
        end if
        if ( ios /= 0 ) stop "Write error in file unit "
      end do

      filename2 = 'parts.dat'
      call get_unit(iounit2)
      open(unit=iounit2, file=filename2, iostat=ios, status="replace", action="write")
      if ( ios /= 0 ) stop "Error opening file "

      if (is_bulk) then
        ! For bulk, calculate probability density for parts
        do i = 1, 8
          do j = 1, evnum
            parts(j,i) = real(conjg(C(j,i))*C(j,i))
          end do
        end do
        if (allocated(C)) deallocate(C)
      else
        ! For quantum well, integrate over z
        do i = 1, 8
          do j = 1, evnum
            parts(j,i) = simpson(dcmplx(R(j,(i-1)*fdstep+1:i*fdstep),0.0_dp),z(1),z(fdstep))
          end do
        end do
        parts = parts/(z(2)-z(1))
        if (allocated(R)) deallocate(R)
      end if

      write(unit=iounit2, fmt="(8(g14.6))",iostat=ios) &
      & (parts(j,:), j=1,evnum)
      if ( ios /= 0 ) stop "Write error in file unit "

      close ( unit = iounit )
      close ( unit = iounit2 )

    end subroutine

    subroutine writeEigenvalues(smallk, eig, wvStep)

      real(kind=dp), intent(in) :: eig(:,:)
      type(wavevector), intent(in) :: smallk(:)
      integer, intent(in) :: wvStep

      integer (kind=4) :: iounit, ios
      character (len = 255) :: filename
      integer :: i

      filename = 'eigenvalues.dat'
      call get_unit(iounit)

      open(unit=iounit, file=filename, iostat=ios, status="replace", action="write")
      if ( ios /= 0 ) stop "Error opening file "

      write(iounit, *) '#k, values'
      do i = 1, wvStep, 1
        write(unit=iounit, fmt="(1000(1x,g14.6))", iostat=ios) &
        & smallk(i)%kx, eig(:,i)
        if ( ios /= 0 ) stop "Write error in file unit "
      end do

      close ( unit = iounit )

    end subroutine

    subroutine get_unit ( iunit )

    !*****************************************************************************80
    !
    !! GET_UNIT returns a free FORTRAN unit number.
    !
    !  Discussion:
    !
    !    A "free" FORTRAN unit number is a value between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5, 6 and 9, which
    !    are commonly reserved for console I/O).
    !
    !    Otherwise, IUNIT is a value between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 October 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) IUNIT, the free unit number.
    !
      implicit none

      integer ( kind = 4 ) i
      integer ( kind = 4 ) ios
      integer ( kind = 4 ) iunit
      logical ( kind = 4 ) lopen

      iunit = 0

      do i = 1, 99

        if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

          inquire ( unit = i, opened = lopen, iostat = ios )

          if ( ios == 0 ) then
            if ( .not. lopen ) then
              iunit = i
              return
            end if
          end if

        end if

      end do

      return
    end



end module outputFunctions