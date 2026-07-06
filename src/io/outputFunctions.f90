module outputFunctions

  use definitions, only: dp, simulation_config, spatial_grid, &
    wavevector, optical_transition, topological_result
  use utils
  use topological_analysis, only: polarization_result_t

  implicit NONE

  private
  public :: writeEigenfunctions, writeEigenfunctions2d, writeParts2d, writeEigenvalues
  public :: write_bdg_eigenvalues, write_optical_transitions, write_profile_1d
  public :: write_bdg_ldos, write_bdg_ldos_nambu, write_bdg_spectral
  public :: write_topology_result, write_majorana_profile, write_bdg_lowest_state_profile, &
    & write_spectral_function, write_z2_phase_diagram, write_z2_transitions, &
    & write_majorana_polarization

  character(len=*), parameter :: OUTPUT_DIR = 'output'

  CONTAINS

    subroutine writeEigenfunctions(N, evnum, A, k, fdstep, z, is_bulk, k_magnitude)

      integer, intent(in) :: N, evnum, k, fdstep
      real(kind=dp), intent(in) :: z(:)
      complex(kind=dp), intent(in), contiguous :: A(:,:)
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
        if (trim(cfg%confinement) == 'wire') then
          write(iounit, '(A,A,A,I0)', advance='no') &
            & '# confinement: ', trim(cfg%confinement), &
            & ' nx: ', cfg%grid%nx
          write(iounit, '(A,I0,A,g14.6)', advance='no') &
            & ' ny: ', cfg%grid%ny, ' dx: ', cfg%grid%dx
          write(iounit, '(A,g14.6,A,A)') &
            & ' dy: ', cfg%grid%dy, ' shape: ', trim(cfg%wire%geom%shape)
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

    subroutine get_eigenvector_component(eigenvectors, eigenstate_index, band_index, fdstep, total_points, component)
      implicit none
      complex(kind=dp), intent(in), contiguous :: eigenvectors(:,:)
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
      complex(kind=dp), intent(in), contiguous    :: eigenvectors(:,:)
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

      Ngrid = grid%npoints()
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
      complex(kind=dp), intent(in), contiguous    :: eigenvectors(:,:)
      integer, intent(in)             :: nev

      integer(kind=4) :: iounit, ios
      integer :: n, band, ix, iy, Ngrid, flat_idx, nev_actual
      real(kind=dp) :: dA
      real(kind=dp), allocatable :: parts(:,:)

      ! Ensure output directory exists
      call ensure_output_dir()

      Ngrid = grid%npoints()
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

    ! ==================================================================
    ! Write BdG eigenvalues to output/bdg_eigenvalues.dat.
    !
    ! Consolidates the two near-identical write blocks formerly inlined
    ! in main_topology.f90 (BdG-wire and BdG-QW paths).  The two paths
    ! differ only in the second header line label ("kz" vs "k_par") and
    ! the scalar printed there; the remaining bytes (header text, count
    ! line, column-order line, data format '(I6,ES20.12)', confirmation
    ! message) are identical.  axis_label supplies the discriminator.
    ! ==================================================================
    subroutine write_bdg_eigenvalues(eigvals, axis_label, kval)

      real(kind=dp), intent(in), contiguous :: eigvals(:)
      character(len=*), intent(in) :: axis_label  ! e.g. 'kz' or 'k_par'
      real(kind=dp), intent(in) :: kval

      integer(kind=4) :: iounit, ios, status
      integer :: i, nev

      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/bdg_eigenvalues.dat', status='replace', &
           action='write', iostat=status)
      if (status /= 0) return

      nev = size(eigvals)
      write(iounit, '(A)') '# BdG eigenvalues (eV)'
      write(iounit, '(A,ES16.8)') '# ' // axis_label // ' (1/A) = ', kval
      write(iounit, '(A,I0)') '# n_eigenvalues = ', nev
      write(iounit, '(A)') '# Columns: index, energy (eV)'
      do i = 1, nev
        write(iounit, '(I6,ES20.12)') i, eigvals(i)
      end do
      close(iounit)
      print *, '  Eigenvalues written to output/bdg_eigenvalues.dat'

    end subroutine write_bdg_eigenvalues

    ! ==================================================================
    ! Write optical transition table to output/optical_transitions.dat.
    !
    ! Consolidates the two byte-identical write blocks formerly inlined
    ! in main_gfactor.f90 (wire path and QW path).  Format preserved
    ! exactly: header comment, then '(2(I4,1x),5(g14.6,1x))' rows.
    ! ==================================================================
    subroutine write_optical_transitions(transitions)

      type(optical_transition), intent(in), contiguous :: transitions(:)
      integer(kind=4) :: iounit, ios
      integer :: it, num_trans

      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/optical_transitions.dat', status='replace', action='write')
      write(iounit, '(A)') '# CB VB dE(eV) |px|^2 |py|^2 |pz|^2 f_osc'
      num_trans = size(transitions)
      do it = 1, num_trans
        write(iounit, '(2(I4,1x),5(g14.6,1x))') &
          transitions(it)%cb_idx, transitions(it)%vb_idx, &
          transitions(it)%energy, transitions(it)%px, &
          transitions(it)%py, transitions(it)%pz, &
          transitions(it)%oscillator_strength
      end do
      close(iounit)
      print *, '  Optical transitions written to output/optical_transitions.dat'

    end subroutine write_optical_transitions

    ! ==================================================================
    ! Write a 1D band-edge profile to output/potential_profile.dat.
    !
    ! Consolidates the byte-identical list-directed writes formerly
    ! inlined in main.f90 (QW and Landau), main_gfactor.f90 (QW) and
    ! main_optics.f90 (QW): one coordinate column + three profile
    ! columns (EV, EV_DeltaSO, EC), written with 'write(unit,*)'.
    ! The coordinate array (cfg%z for QW, cfg%grid%x for Landau) is
    ! supplied by the caller.
    ! ==================================================================
    subroutine write_profile_1d(coord, profile)

      real(kind=dp), intent(in), contiguous :: coord(:)
      real(kind=dp), intent(in), contiguous :: profile(:,:)
      integer(kind=4) :: iounit
      integer :: i, n

      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/potential_profile.dat', status='replace', action='write')
      n = size(coord)
      do i = 1, n
        write(iounit, *) coord(i), profile(i,1), profile(i,2), profile(i,3)
      end do
      close(iounit)

    end subroutine write_profile_1d

    ! ==================================================================
    ! Issue 06 / Unit U9: BdG LDOS at the specified energy grid.
    !
    ! Writes output/bdg_ldos.dat with three columns: r, E, ldos(r,E).
    ! Header records the kz point. Append mode across multiple kz points
    ! keeps a single multi-row file (the dispatch runs one kz per call).
    ! ==================================================================
    subroutine write_bdg_ldos(E_values, ldos, kz_val)

      real(kind=dp), intent(in), contiguous :: E_values(:)
      real(kind=dp), intent(in), contiguous :: ldos(:)
      real(kind=dp), intent(in) :: kz_val
      integer(kind=4) :: iounit, status
      integer :: iE, n

      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/bdg_ldos.dat', status='replace', action='write', &
           iostat=status)
      if (status /= 0) return

      n = size(E_values)
      write(iounit, '(A)') '# BdG LDOS (1/eV)'
      write(iounit, '(A,ES16.8)') '# kz (1/A) = ', kz_val
      write(iounit, '(A)') '# Columns: r, E (eV), LDOS (1/eV)'
      do iE = 1, n
        write(iounit, '(I6,1X,ES16.8,1X,ES16.8)') 1, E_values(iE), ldos(iE)
      end do
      close(iounit)

    end subroutine write_bdg_ldos

    ! ==================================================================
    ! Issue 06 / Unit U9: Nambu-resolved LDOS at E=0 (Majorana peak).
    !
    ! Writes output/bdg_ldos_nambu.dat with three columns:
    !   r, LDOS_electron(r, E=0), LDOS_hole(r, E=0)
    ! Header records the kz point. Single-row per r for the Majorana
    ! zero-energy slice.
    ! ==================================================================
    subroutine write_bdg_ldos_nambu(ldos_e, ldos_h, kz_val)

      real(kind=dp), intent(in), contiguous :: ldos_e(:)
      real(kind=dp), intent(in), contiguous :: ldos_h(:)
      real(kind=dp), intent(in) :: kz_val
      integer(kind=4) :: iounit, status
      integer :: r, N

      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/bdg_ldos_nambu.dat', status='replace', action='write', &
           iostat=status)
      if (status /= 0) return

      N = size(ldos_e)
      write(iounit, '(A)') '# BdG Nambu-resolved LDOS at E=0 (1/eV)'
      write(iounit, '(A,ES16.8)') '# kz (1/A) = ', kz_val
      write(iounit, '(A)') '# Columns: r, LDOS_electron, LDOS_hole'
      do r = 1, N
        write(iounit, '(I6,1X,ES16.8,1X,ES16.8)') r, ldos_e(r), ldos_h(r)
      end do
      close(iounit)

    end subroutine write_bdg_ldos_nambu

    ! ==================================================================
    ! Issue 06 / Unit U9: BdG spectral function A(k,E).
    !
    ! Writes output/bdg_spectral.dat as a 2D-grid: rows are k-points,
    ! columns are E-points. Header documents the k and E ranges.
    ! Format: 'i, k_i, j, E_j, A_ij' per row.
    ! ==================================================================
    subroutine write_bdg_spectral(k_values, E_values, A_kE)

      real(kind=dp), intent(in), contiguous :: k_values(:)
      real(kind=dp), intent(in), contiguous :: E_values(:)
      real(kind=dp), intent(in), contiguous :: A_kE(:,:)
      integer(kind=4) :: iounit, status
      integer :: ik, iE, nk, nE

      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/bdg_spectral.dat', status='replace', action='write', &
           iostat=status)
      if (status /= 0) return

      nk = size(k_values)
      nE = size(E_values)
      write(iounit, '(A)') '# BdG spectral function A(k, E) (1/eV)'
      write(iounit, '(A,I0)') '# nk = ', nk
      write(iounit, '(A,I0)') '# nE = ', nE
      write(iounit, '(A)') '# Columns: ik, k (1/A), iE, E (eV), A(k,E)'
      do ik = 1, nk
        do iE = 1, nE
          write(iounit, '(I6,1X,ES16.8,1X,I6,1X,ES16.8,1X,ES16.8)') &
            ik, k_values(ik), iE, E_values(iE), A_kE(ik, iE)
        end do
      end do
      close(iounit)

    end subroutine write_bdg_spectral

    ! ==================================================================
    ! Task 3.3: Write topological result summary to output/topology_result.dat.
    !
    ! Replaces the inline block formerly at main_topology.f90:236-261.
    ! Consolidates header + diagnostics + (optional) edge-state energies.
    ! Error path mirrors original (ERROR print + error stop on open failure).
    ! ==================================================================
    subroutine write_topology_result(cfg, topo_result)

      type(simulation_config), intent(in) :: cfg
      type(topological_result), intent(in) :: topo_result

      integer(kind=4) :: iounit, status
      integer :: i

      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/topology_result.dat', status='replace', &
           action='write', iostat=status)
      if (status /= 0) then
        print *, 'ERROR: cannot open output/topology_result.dat (iostat=', status, ')'
        error stop 'cannot open topology_result.dat'
      end if
      write(iounit, '(A)') '# Topological Analysis Results'
      write(iounit, '(A,A)') '# mode: ', trim(cfg%topo%mode)
      write(iounit, '(A,I0)') '# Chern number: ', topo_result%chern_number
      write(iounit, '(A,I0)') '# Z2 invariant: ', topo_result%z2_invariant
      write(iounit, '(A,F12.6)') '# Hall conductance (e^2/h): ', topo_result%hall_conductance
      write(iounit, '(A,F12.6)') '# Conductance xy (e^2/h): ', topo_result%conductance_xy
      write(iounit, '(A,F12.6)') '# Conductance zz: ', topo_result%conductance_zz
      write(iounit, '(A,I0)') '# Majorana count: ', topo_result%n_majorana
      write(iounit, '(A,I0)') '# Majorana fit failures: ', topo_result%n_majorana_fit_failed
      write(iounit, '(A,F12.6)') '# Min gap (eV): ', topo_result%min_gap
      write(iounit, '(A,F12.6)') '# Edge localization length min (AA): ', topo_result%edge_xi_min
      write(iounit, '(A,F12.6)') '# Edge localization length avg (AA): ', topo_result%edge_xi
      if (allocated(topo_result%edge_energies)) then
        write(iounit, '(A)') '# Edge state energies (eV):'
        do i = 1, size(topo_result%edge_energies)
          write(iounit, '(F12.6)') topo_result%edge_energies(i)
        end do
      end if
      close(iounit)
      print *, '  Results written to output/topology_result.dat'

    end subroutine write_topology_result

    ! ==================================================================
    ! Task 3.3: Write Majorana wavefunction spatial profile.
    !
    ! Replaces the inline blocks formerly at main_topology.f90:614-628 (wire,
    ! 2D coords) and :910-923 (QW, 1D z-grid). The axis_label and coord data
    ! are passed by the caller so the same writer serves both geometries:
    !   - wire: coord is cfg%grid%coords (2, nspatial) -> index, x, y, rho
    !   - QW:   coord is qw_grid%z (nspatial)          -> z, rho
    ! axis_label distinguishes 'kz' (wire) from 'k_par' (QW) in the header.
    ! ==================================================================
    subroutine write_majorana_profile(coord, profile, axis_label, axis_val, xi_val, n)

      real(kind=dp), intent(in), contiguous :: coord(:)
      real(kind=dp), intent(in), contiguous :: profile(:)
      character(len=*), intent(in) :: axis_label
      real(kind=dp), intent(in) :: axis_val, xi_val
      integer, intent(in) :: n

      integer(kind=4) :: iounit, status
      integer :: j
      logical :: is_2d

      is_2d = (size(coord) == 2 * n)

      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/majorana_profile.dat', &
           status='replace', action='write', iostat=status)
      if (status /= 0) return

      write(iounit, '(A)') '# Majorana wavefunction spatial profile'
      write(iounit, '(A,ES16.8)') '# ' // axis_label // ' (1/A) = ', axis_val
      write(iounit, '(A,ES16.8)') '# xi (AA) = ', xi_val
      write(iounit, '(A,I0)') '# n_spatial = ', n
      if (is_2d) then
        write(iounit, '(A)') '# Columns: index, x (AA), y (AA), rho'
        do j = 1, n
          write(iounit, '(I6,2ES16.8,ES20.12)') j, coord(j), coord(n + j), profile(j)
        end do
      else
        write(iounit, '(A)') '# Columns: z (AA), rho'
        do j = 1, n
          write(iounit, '(ES16.8,ES20.12)') coord(j), profile(j)
        end do
      end if
      close(iounit)
      print *, '  Majorana profile written to output/majorana_profile.dat'

    end subroutine write_majorana_profile

    ! ==================================================================
    ! PR #41 A.3a (Issue 04 / U7): Write Majorana polarization profile
    ! (Sticlet P_M and <tau_z>) to output/majorana_polarization.dat.
    !
    ! Consumes a polarization_result_t produced by
    ! topological_analysis::majorana_polarization applied to one BdG
    ! eigenvector. The file format is fixed-width, parseable by
    ! tests/integration/verify_majorana_polarization.py (A.3b): one row
    ! per site with columns (index, P_M, tau_z, half_wire_integral).
    ! half_wire_integral is a scalar so it is repeated on every row for
    ! convenience (verifier reads row 1 for the integral).
    ! ==================================================================
    subroutine write_majorana_polarization(pol, filename)

      type(polarization_result_t), intent(in) :: pol
      character(len=*), intent(in) :: filename

      integer(kind=4) :: iounit, status
      integer :: i, n

      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file=filename, status='replace', action='write', iostat=status)
      if (status /= 0) then
        print *, 'ERROR: cannot open ', filename
        error stop 'cannot open majorana_polarization.dat'
      end if

      n = size(pol%P_M)
      write(iounit, '(A)') '# Majorana polarization profile (Sticlet P_M, <tau_z>)'
      write(iounit, '(A,I0)') '# n_sites = ', n
      write(iounit, '(A,ES15.7)') '# half_wire_integral = ', pol%half_wire_integral
      write(iounit, '(A)') '# Columns: index, P_M, tau_z, half_wire_integral'
      do i = 1, n
        write(iounit, '(I6,3(1X,ES15.7))') i, pol%P_M(i), pol%tau_z(i), pol%half_wire_integral
      end do
      close(iounit)
      print *, '  Majorana polarization written to ', filename

    end subroutine write_majorana_polarization

    ! ==================================================================
    ! Task 3.3: Write the lowest-|E| BdG state spatial profile.
    !
    ! Replaces the inline block formerly at main_topology.f90:972-986.
    ! Header records k_par, the minimum-energy BdG eigenvalue, and xi.
    ! Format: 'z (AA), rho' (one row per spatial point on the QW z-grid).
    ! ==================================================================
    subroutine write_bdg_lowest_state_profile(z_grid, profile, k_par_val, energy, xi_val)

      real(kind=dp), intent(in), contiguous :: z_grid(:)
      real(kind=dp), intent(in), contiguous :: profile(:)
      real(kind=dp), intent(in) :: k_par_val, energy, xi_val

      integer(kind=4) :: iounit, status
      integer :: j, n

      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/bdg_lowest_state_profile.dat', &
           status='replace', action='write', iostat=status)
      if (status /= 0) return

      n = size(z_grid)
      write(iounit, '(A)') '# Lowest-|E| BdG state spatial profile'
      write(iounit, '(A,ES16.8)') '# k_par (1/A) = ', k_par_val
      write(iounit, '(A,ES16.8)') '# Energy (eV) = ', energy
      write(iounit, '(A,ES16.8)') '# xi (AA) = ', xi_val
      write(iounit, '(A)') '# Columns: z (AA), rho'
      do j = 1, n
        write(iounit, '(ES16.8,ES20.12)') z_grid(j), profile(j)
      end do
      close(iounit)
      print *, '  Lowest-state profile written to output/bdg_lowest_state_profile.dat'

    end subroutine write_bdg_lowest_state_profile

    ! ==================================================================
    ! Task 3.3: Write spectral function A(k, E) to output/spectral_function.dat.
    !
    ! Replaces the inline block formerly at main_topology.f90:1070-1091.
    ! The writer reads only the data the header needs (confinement, eta, k/E
    ! ranges) plus the A(k,E) array itself, so the call site stays compact.
    ! Format: 'k (1/A), E (eV), A(k,E)' (one row per (ik, iE) pair).
    ! ==================================================================
    subroutine write_spectral_function(cfg, k_arr, E_arr, A_kE)

      type(simulation_config), intent(in) :: cfg
      real(kind=dp), intent(in), contiguous :: k_arr(:)
      real(kind=dp), intent(in), contiguous :: E_arr(:)
      real(kind=dp), intent(in), contiguous :: A_kE(:,:)

      integer(kind=4) :: iounit, status
      integer :: ik, iE, nk, nE

      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/spectral_function.dat', status='replace', &
           action='write', iostat=status)
      if (status /= 0) then
        print *, 'ERROR: cannot open output/spectral_function.dat'
        error stop 'cannot open spectral_function.dat'
      end if

      nk = size(k_arr)
      nE = size(E_arr)
      write(iounit, '(A)') '# Spectral function A(k, E)'
      write(iounit, '(A,A)') '# confinement=', trim(cfg%confinement)
      write(iounit, '(A,I0,A,I0)') '# nk=', nk, '  nE=', nE
      write(iounit, '(A,ES12.4)') '# eta (eV) = ', cfg%topo%spectral_eta
      write(iounit, '(A,2ES16.8)') '# k_grid_min_max (1/A) = ', k_arr(1), k_arr(nk)
      write(iounit, '(A,2ES16.8)') '# E_grid_min_max (eV) = ', E_arr(1), E_arr(nE)
      write(iounit, '(A)') '# units: k in 1/A, E in eV, A in 1/eV'
      write(iounit, '(A)') '# Columns: k (1/A), E (eV), A(k,E)'
      do ik = 1, nk
        do iE = 1, nE
          write(iounit, '(3ES16.8)') k_arr(ik), E_arr(iE), A_kE(ik, iE)
        end do
      end do
      close(iounit)
      print *, '  Spectral function written to output/spectral_function.dat'

    end subroutine write_spectral_function

    ! ==================================================================
    ! Task 3.3: Write Z2 phase diagram to output/z2_phase_diagram.dat.
    !
    ! Replaces the inline block formerly at main_topology.f90:1232-1253.
    ! B and mu grids are reconstructed from cfg%topo sweep ranges so the
    ! caller only passes the result grids (z2_map, gap_map). The (iB, iMu)
    ! loop order matches the original: rows iterate over B, columns over mu.
    ! ==================================================================
    subroutine write_z2_phase_diagram(cfg, z2_map, gap_map)

      type(simulation_config), intent(in) :: cfg
      real(kind=dp), intent(in), contiguous :: z2_map(:,:)
      real(kind=dp), intent(in), contiguous :: gap_map(:,:)

      integer(kind=4) :: iounit, status
      integer :: iB, iMu, nB, nMu

      nB = size(z2_map, 2)
      nMu = size(z2_map, 1)

      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/z2_phase_diagram.dat', status='replace', &
           action='write', iostat=status)
      if (status /= 0) then
        print *, 'ERROR: cannot open output/z2_phase_diagram.dat'
        error stop 'cannot open z2_phase_diagram.dat'
      end if
      write(iounit, '(A)') '# Z2 phase diagram'
      write(iounit, '(A,I0,A,I0)') '# nB=', nB, '  nMu=', nMu
      write(iounit, '(A)') '# B(T) mu(eV) z2 gap(eV)'
      do iB = 1, nB
        do iMu = 1, nMu
          write(iounit, '(4ES16.8)') &
            & cfg%topo%gap_sweep_B_min + real(iB - 1, kind=dp) * &
            & (cfg%topo%gap_sweep_B_max - cfg%topo%gap_sweep_B_min) / real(max(1, nB - 1), kind=dp), &
            & cfg%topo%gap_sweep_mu_min + real(iMu - 1, kind=dp) * &
            & (cfg%topo%gap_sweep_mu_max - cfg%topo%gap_sweep_mu_min) / real(max(1, nMu - 1), kind=dp), &
            & z2_map(iMu, iB), gap_map(iMu, iB)
        end do
      end do
      close(iounit)
      print *, '  Z2 phase diagram written to output/z2_phase_diagram.dat'

    end subroutine write_z2_phase_diagram

    ! ==================================================================
    ! Task 3.3: Write Z2 phase transition list to output/z2_transitions.dat.
    !
    ! Replaces the inline block formerly at main_topology.f90:1255-1266.
    ! One row per detected transition: 'B (T), mu (eV)'. Silent on empty.
    ! ==================================================================
    subroutine write_z2_transitions(transitions)

      real(kind=dp), intent(in), contiguous :: transitions(:,:)
      integer(kind=4) :: iounit, status
      integer :: iB, n

      n = size(transitions, 1)

      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/z2_transitions.dat', status='replace', &
           action='write', iostat=status)
      if (status /= 0) then
        print *, 'ERROR: cannot open output/z2_transitions.dat'
        error stop 'cannot open z2_transitions.dat'
      end if
      write(iounit, '(A)') '# B(T) mu(eV)'
      do iB = 1, n
        write(iounit, '(2ES16.8)') transitions(iB, 1), transitions(iB, 2)
      end do
      close(iounit)

    end subroutine write_z2_transitions

end module outputFunctions