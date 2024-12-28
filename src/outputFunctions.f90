module outputFunctions

  use definitions
  use utils

  implicit NONE

  private
  public :: writeEigenfunctions, writeEigenvalues, get_unit

  character(len=*), parameter :: OUTPUT_DIR = 'output'

  CONTAINS

    subroutine ensure_output_dir()
      ! Local variables
      logical :: dir_exists
      integer :: status

      ! Check if directory exists
      inquire(file=OUTPUT_DIR//'/.', exist=dir_exists)
      
      ! Create directory if it doesn't exist
      if (.not. dir_exists) then
        call system('mkdir -p '//OUTPUT_DIR)
      end if
    end subroutine

    subroutine writeEigenfunctions(N, evnum, A, k, fdstep, z, is_bulk)

      integer, intent(in) :: N, evnum, k, fdstep
      real(kind=dp), intent(in) :: z(:)
      complex(kind=dp), intent(in) :: A(:,:)
      logical, intent(in) :: is_bulk

      integer (kind=4) :: iounit, iounit2, ios
      character (len = 255) :: filename
      integer :: i, j, m
      real(kind=dp), allocatable :: parts(:,:)
      real(kind=dp), allocatable :: eigv_abs(:,:)
      real(kind=dp), allocatable :: component(:)

      ! Ensure output directory exists
      call ensure_output_dir()

      ! Allocate arrays
      allocate(parts(evnum,8))
      allocate(component(fdstep))
      parts = 0.0_dp
      if (.not. is_bulk) then
        allocate(eigv_abs(fdstep,8))
        eigv_abs = 0.0_dp
      end if

      ! Write eigenfunctions
      do j = 1, evnum
        write(filename,'(a,i0.5,a,i0.5,a)') OUTPUT_DIR//'/eigenfunctions_k_', k, '_ev_', j, '.dat'
        call get_unit(iounit)
        open(unit=iounit, file=filename, status='replace', action='write')

        if (is_bulk) then
          ! For bulk, write only the components
          write(unit=iounit, fmt="(8(g14.6))", iostat=ios) &
          & (abs(A(i,j)), i=1,8)
        else
          ! For quantum well, get each component safely
          do m = 1, 8
            call get_eigenvector_component(A, j, m, fdstep, N, component)
            eigv_abs(:,m) = abs(component)
          end do

          ! Write z and all components
          do i = 1, fdstep
            write(unit=iounit, fmt="(9(g14.6))", iostat=ios) &
            & z(i), (eigv_abs(i,m), m=1,8)
          end do
        end if

        close(iounit)
      end do

      ! Calculate and write parts
      call get_unit(iounit2)
      open(unit=iounit2, file=OUTPUT_DIR//'/parts.dat', status='replace', action='write')

      if (is_bulk) then
        ! For bulk, calculate probability density for parts
        do j = 1, evnum
          do i = 1, 8
            parts(j,i) = abs(A(i,j))**2
          end do
        end do
      else
        ! For quantum well, integrate over z
        do j = 1, evnum
          ! Initialize parts for this eigenstate
          parts(j,:) = 0.0_dp

          ! Get each component safely and calculate its contribution
          do m = 1, 8
            call get_eigenvector_component(A, j, m, fdstep, N, component)
            eigv_abs(:,m) = abs(component)
            ! Calculate parts by integrating over z
            parts(j,m) = sum(eigv_abs(:,m)**2) * (z(2)-z(1))
          end do
        end do
      end if

      ! Write parts
      do j = 1, evnum
        write(unit=iounit2, fmt="(8(g14.6))", iostat=ios) (parts(j,m), m=1,8)
      end do

      close(iounit2)
      if (.not. is_bulk) deallocate(eigv_abs)
      deallocate(parts)
      deallocate(component)

    end subroutine writeEigenfunctions

    subroutine writeEigenvalues(smallk, eig, wvStep)

      real(kind=dp), intent(in) :: eig(:,:)
      type(wavevector), intent(in) :: smallk(:)
      integer, intent(in) :: wvStep

      integer (kind=4) :: iounit, ios
      character (len = 255) :: filename
      integer :: i

      ! Ensure output directory exists
      call ensure_output_dir()

      filename = OUTPUT_DIR//'/eigenvalues.dat'
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

    subroutine get_eigenvector_component(eigenvectors, eigenstate_index, band_index, fdstep, total_points, component)
      implicit none
      complex(kind=dp), intent(in) :: eigenvectors(:,:)
      integer, intent(in) :: eigenstate_index, band_index, fdstep, total_points
      real(kind=dp), intent(out), dimension(fdstep) :: component
      integer :: i, base_idx

      component = 0.0_dp  ! Initialize all components to zero
      if (band_index <= 8) then  ! Only process valid band indices
        do i = 1, fdstep
          base_idx = (band_index - 1) * fdstep + i
          if (base_idx <= total_points) then
            component(i) = abs(eigenvectors(base_idx, eigenstate_index))
          end if
        end do
      end if
    end subroutine

end module outputFunctions