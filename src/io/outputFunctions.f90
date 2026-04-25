module outputFunctions

  use definitions
  use utils

  implicit NONE

  private
  public :: writeEigenfunctions, writeEigenfunctions2d, writeParts2d, writeEigenvalues, get_unit, ensure_output_dir

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

    subroutine writeEigenfunctions(N, evnum, A, k, fdstep, z, is_bulk, k_magnitude)

      integer, intent(in) :: N, evnum, k, fdstep
      real(kind=dp), intent(in) :: z(:)
      complex(kind=dp), intent(in) :: A(:,:)
      logical, intent(in) :: is_bulk
      real(kind=dp), intent(in), optional :: k_magnitude  ! |k| for gnuplot header

      integer (kind=4) :: iounit, iounit2, ios
      character (len = 255) :: filename
      integer :: i, j, m
      real(kind=dp) :: kmag
      real(kind=dp), allocatable :: parts(:,:)
      real(kind=dp), allocatable :: eigv_abs(:,:)
      real(kind=dp), allocatable :: component(:)

      kmag = 0.0_dp
      if (present(k_magnitude)) kmag = k_magnitude

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

      if (is_bulk) then
        ! For bulk, calculate probability density for parts
        do j = 1, evnum
          do i = 1, 8
            parts(j,i) = abs(A(i,j))**2
          end do
        end do

        ! Write parts in multi-block gnuplot format (append for k > 1)
        if (k <= 1) then
          open(unit=iounit2, file=OUTPUT_DIR//'/parts.dat', status='replace', action='write')
        else
          open(unit=iounit2, file=OUTPUT_DIR//'/parts.dat', position='append', action='write')
        end if

        ! k-point header (gnuplot index separator)
        write(unit=iounit2, fmt='("# k = ", g14.6)', iostat=ios) kmag

        ! Write parts for this k-point
        do j = 1, evnum
          write(unit=iounit2, fmt="(8(g14.6))", iostat=ios) (parts(j,m), m=1,8)
        end do

        ! Blank line as gnuplot block separator
        write(unit=iounit2, fmt='("")', iostat=ios)

        close(iounit2)
      else
        ! For quantum well, integrate over z
        open(unit=iounit2, file=OUTPUT_DIR//'/parts.dat', status='replace', action='write')

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

        ! Write parts
        do j = 1, evnum
          write(unit=iounit2, fmt="(8(g14.6))", iostat=ios) (parts(j,m), m=1,8)
        end do

        close(iounit2)
      end if
      if (.not. is_bulk) deallocate(eigv_abs)
      deallocate(parts)
      deallocate(component)

    end subroutine writeEigenfunctions

    subroutine writeEigenvalues(smallk, eig, wvStep, cfg)

      real(kind=dp), intent(in) :: eig(:,:)
      type(wavevector), intent(in) :: smallk(:)
      integer, intent(in) :: wvStep
      type(simulation_config), intent(in), optional :: cfg

      integer (kind=4) :: iounit, ios
      character (len = 255) :: filename
      integer :: i

      ! Ensure output directory exists
      call ensure_output_dir()

      filename = OUTPUT_DIR//'/eigenvalues.dat'
      call get_unit(iounit)

      open(unit=iounit, file=filename, iostat=ios, status="replace", action="write")
      if ( ios /= 0 ) stop "Error opening file "

      ! Write metadata header for wire mode
      if (present(cfg)) then
        if (cfg%confinement == 2) then
          write(iounit, '(A,I0,A,I0)', advance='no') &
            & '# confinement: ', cfg%confinement, &
            & ' nx: ', cfg%grid%nx
          write(iounit, '(A,I0,A,g14.6)', advance='no') &
            & ' ny: ', cfg%grid%ny, ' dx: ', cfg%grid%dx
          write(iounit, '(A,g14.6,A,A)') &
            & ' dy: ', cfg%grid%dy, ' shape: ', trim(cfg%wire_geom%shape)
        end if
      end if

      write(iounit, *) '#k, values'
      do i = 1, wvStep, 1
        write(unit=iounit, fmt="(*(1x,g14.6))", iostat=ios) &
        & sqrt(smallk(i)%kx**2 + smallk(i)%ky**2 + smallk(i)%kz**2), eig(:,i)
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

    ! ==================================================================
    ! Write 2D probability density |psi(x,y)|^2 for wire eigenstates.
    !
    ! Each eigenvector has 8*Ngrid complex entries.  The probability
    ! density at grid point (ix, iy) is summed over all 8 bands:
    !   |psi(x,y)|^2 = sum_{band=1}^{8} |vec((band-1)*Ngrid + (iy-1)*grid%nx + ix)|^2
    !
    ! Output format (one file per eigenstate per k-point):
    !   x  y  |psi|^2     (three columns, nx*ny rows)
    ! Directly plottable with gnuplot: splot 'file' using 1:2:3
    ! ==================================================================
    subroutine writeEigenfunctions2d(grid, eigenvalues, eigenvectors, k_index, nev, write_parts)

      type(spatial_grid), intent(in)  :: grid
      real(kind=dp), intent(in)       :: eigenvalues(:)
      complex(kind=dp), intent(in)    :: eigenvectors(:,:)
      integer, intent(in)             :: k_index
      integer, intent(in)             :: nev
      logical, intent(in), optional   :: write_parts

      integer(kind=4) :: iounit, ios
      character(len=255) :: filename
      integer :: n, band, ix, iy, Ngrid, flat_idx, nev_actual
      real(kind=dp) :: prob
      logical :: do_parts

      ! Ensure output directory exists
      call ensure_output_dir()

      Ngrid = grid%nx * grid%ny
      if (Ngrid == 0) return

      nev_actual = min(nev, size(eigenvectors, 2))

      do_parts = .false.
      if (present(write_parts)) do_parts = write_parts

      ! Compute band decomposition if requested (delegate to writeParts2d)
      if (do_parts) then
        call writeParts2d(grid, eigenvectors, nev_actual)
      end if

      do n = 1, nev_actual
        write(filename,'(a,i0.5,a,i0.5,a)') OUTPUT_DIR//'/eigenfunctions_k_', k_index, '_ev_', n, '.dat'
        call get_unit(iounit)
        open(unit=iounit, file=filename, status='replace', action='write')

        ! Header with eigenvalue
        write(iounit, '(a,g14.6)') '# E = ', eigenvalues(n)

        ! Write probability density grid: x  y  |psi|^2
        do iy = 1, grid%ny
          do ix = 1, grid%nx
            prob = 0.0_dp
            do band = 1, 8
              flat_idx = (band - 1) * Ngrid + (iy - 1) * grid%nx + ix
              prob = prob + abs(eigenvectors(flat_idx, n))**2
            end do
            write(unit=iounit, fmt='(3(g14.6,1x))', iostat=ios) &
              & grid%x(ix), grid%z(iy), prob
          end do
          ! Blank line between y-rows for gnuplot splot
          write(iounit, '(a)', iostat=ios) ''
        end do

        close(iounit)
      end do

    end subroutine writeEigenfunctions2d

    ! ==================================================================
    ! Write 8-band decomposition (parts.dat) for wire eigenstates.
    !
    ! Decomposes each eigenvector into contributions from each of the
    ! 8 basis bands by integrating |psi_b(x,y)|^2 over the 2D grid.
    ! The flat index layout matches the Hamiltonian construction:
    !   flat_idx = (band-1)*Ngrid + (iy-1)*nx + ix
    !
    ! Output: output/parts.dat with 8 columns (one row per eigenstate).
    ! Values are band fractions that sum to 1.0 for each eigenstate.
    ! ==================================================================
    subroutine writeParts2d(grid, eigenvectors, nev)
      type(spatial_grid), intent(in)  :: grid
      complex(kind=dp), intent(in)    :: eigenvectors(:,:)
      integer, intent(in)             :: nev

      integer(kind=4) :: iounit, ios
      integer :: n, band, ix, iy, Ngrid, flat_idx, nev_actual
      real(kind=dp) :: dA
      real(kind=dp), allocatable :: parts(:,:)

      ! Ensure output directory exists
      call ensure_output_dir()

      Ngrid = grid%nx * grid%ny
      if (Ngrid == 0) return

      nev_actual = min(nev, size(eigenvectors, 2))
      dA = grid%dx * grid%dy

      allocate(parts(nev_actual, 8))
      parts = 0.0_dp

      do n = 1, nev_actual
        do band = 1, 8
          do iy = 1, grid%ny
            do ix = 1, grid%nx
              flat_idx = (band - 1) * Ngrid + (iy - 1) * grid%nx + ix
              parts(n, band) = parts(n, band) + abs(eigenvectors(flat_idx, n))**2 * dA
            end do
          end do
        end do
        ! Normalize so band fractions sum to 1
        if (sum(parts(n, :)) > 0.0_dp) then
          parts(n, :) = parts(n, :) / sum(parts(n, :))
        end if
      end do

      ! Write parts.dat
      call get_unit(iounit)
      open(unit=iounit, file=OUTPUT_DIR//'/parts.dat', status='replace', action='write')
      write(iounit, '(A)') '# Band decomposition (wire eigenstates)'
      write(iounit, '(A)') '# HH1  HH2  LH1  LH2  SO1  SO2  CB1  CB2'
      do n = 1, nev_actual
        write(unit=iounit, fmt='(8(g14.6))', iostat=ios) (parts(n, band), band=1, 8)
      end do
      close(iounit)

      deallocate(parts)

    end subroutine writeParts2d

end module outputFunctions