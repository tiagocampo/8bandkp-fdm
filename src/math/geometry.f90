module geometry

  ! ==============================================================================
  ! Cut-cell immersed boundary geometry for 2D quantum wire cross-sections.
  !
  ! Computes cell-volume and face-fraction arrays for arbitrary wire shapes
  ! (rectangle, circle, hexagon, polygon) on a Cartesian FD grid.
  ! All shapes are defined in the y-z plane.  Grid coordinates are cell-centered:
  !   y(i) = (i - 0.5) * dy,   z(j) = (j - 0.5) * dz
  !
  ! Uses analytical intersection wherever possible; falls back to marching-squares
  ! sub-sampling for general polygons.
  !
  ! Public procedures:
  !   grid_init_rect    -- rectangle (trivial, all fractions = 1)
  !   grid_init_circle  -- circular wire (analytical chord intersection)
  !   grid_init_hexagon -- regular hexagon (6 half-plane intersections)
  !   grid_init_polygon -- arbitrary polygon (Sutherland-Hodgman clip)
  !   grid_build_coords -- build flattened coords(:,:) from y(:), z(:)
  !   grid_free         -- deallocate all allocatable fields
  ! ==============================================================================

  use definitions, only: dp, spatial_grid, pi_dp
  use sparse_matrices

  implicit none

  private
  public :: grid_init_rect, grid_init_circle, grid_init_hexagon
  public :: grid_init_polygon, grid_build_coords, grid_free

  ! Number of sub-sample points per cell edge for marching-squares integration
  integer, parameter :: MS_N = 16

contains

  ! ==================================================================
  ! Internal helpers
  ! ==================================================================

  ! ------------------------------------------------------------------
  ! Map 2D (iy, iz) to flattened column-major index: ij = (iz-1)*ny + iy
  ! ------------------------------------------------------------------
  pure function flat_idx(ny, iy, iz) result(ij)
    integer, intent(in) :: ny, iy, iz
    integer :: ij
    ij = (iz - 1) * ny + iy
  end function flat_idx

  ! ------------------------------------------------------------------
  ! Sutherland-Hodgman clip of a convex polygon against a half-plane.
  ! The half-plane is defined by:  n . (p - p0) >= 0
  !
  ! Input:  vx, vz  -- vertex coordinates (nvert points)
  ! Output: clipped polygon appended to cx, cz; nc updated.
  !
  ! Uses work arrays cx, cz which must be large enough.
  ! ------------------------------------------------------------------
  pure subroutine clip_halfplane(vx, vz, nvert, cx, cz, nc, nx, nz, px, pz)
    real(kind=dp), intent(in)  :: vx(nvert), vz(nvert)
    integer, intent(in)        :: nvert
    real(kind=dp), intent(out) :: cx(nvert + 1), cz(nvert + 1)
    integer, intent(out)       :: nc
    real(kind=dp), intent(in)  :: nx, nz     ! outward normal
    real(kind=dp), intent(in)  :: px, pz     ! point on boundary

    integer       :: i, j
    real(kind=dp) :: dx_i, dz_i, dx_j, dz_j, ti, tj, t

    nc = 0
    j = nvert
    dx_j = vx(j) - px
    dz_j = vz(j) - pz
    tj = nx * dx_j + nz * dz_j

    do i = 1, nvert
      dx_i = vx(i) - px
      dz_i = vz(i) - pz
      ti = nx * dx_i + nz * dz_i

      if (ti >= 0.0_dp) then
        ! Current vertex is inside
        if (tj < 0.0_dp) then
          ! Edge crosses from outside to inside -- add intersection
          t = tj / (tj - ti)
          nc = nc + 1
          cx(nc) = vx(j) + t * (vx(i) - vx(j))
          cz(nc) = vz(j) + t * (vz(i) - vz(j))
        end if
        nc = nc + 1
        cx(nc) = vx(i)
        cz(nc) = vz(i)
      else if (tj >= 0.0_dp) then
        ! Edge crosses from inside to outside -- add intersection
        t = tj / (tj - ti)
        nc = nc + 1
        cx(nc) = vx(j) + t * (vx(i) - vx(j))
        cz(nc) = vz(j) + t * (vz(i) - vz(j))
      end if

      j  = i
      tj = ti
    end do
  end subroutine clip_halfplane

  ! ------------------------------------------------------------------
  ! Area of a simple polygon (shoelace formula).
  ! Vertices in order (CCW or CW).  Returns positive area.
  ! ------------------------------------------------------------------
  pure function polygon_area(vx, vz, n) result(area)
    integer, intent(in)        :: n
    real(kind=dp), intent(in)  :: vx(n), vz(n)
    real(kind=dp) :: area
    integer       :: i, j

    area = 0.0_dp
    do i = 1, n
      j = modulo(i, n) + 1
      area = area + vx(i) * vz(j) - vx(j) * vz(i)
    end do
    area = abs(area) * 0.5_dp
  end function polygon_area

  ! ------------------------------------------------------------------
  ! Clip a rectangle (cell) against a convex polygon using successive
  ! half-plane clips (Sutherland-Hodgman).  Returns the clipped area.
  !
  ! The polygon is defined by edges with outward normals computed
  ! from consecutive vertices.  Each edge defines a half-plane that
  ! the cell polygon must be clipped against.
  ! ------------------------------------------------------------------
  pure function clip_cell_against_polygon(ylo, yhi, zlo, zhi, &
      pvx, pvz, npoly) result(area)
    real(kind=dp), intent(in)  :: ylo, yhi, zlo, zhi
    integer, intent(in)        :: npoly
    real(kind=dp), intent(in)  :: pvx(npoly), pvz(npoly)
    real(kind=dp) :: area

    real(kind=dp) :: wx(8), wz(8)
    real(kind=dp) :: ox(8), oz(8)
    integer       :: nw, no, ie, ip
    real(kind=dp) :: ex, ez, nx, nz, len

    ! Start with the cell rectangle (CCW)
    wx(1) = ylo; wz(1) = zlo
    wx(2) = yhi; wz(2) = zlo
    wx(3) = yhi; wz(3) = zhi
    wx(4) = ylo; wz(4) = zhi
    nw = 4

    ! Clip against each edge of the polygon
    do ie = 1, npoly
      ip = modulo(ie, npoly) + 1
      ! Edge from vertex ie to ip
      ex = pvx(ip) - pvx(ie)
      ez = pvz(ip) - pvz(ie)
      len = sqrt(ex * ex + ez * ez)
      if (len < 1.0e-30_dp) cycle
      ! Inward normal (left-hand normal for CCW polygon points inward)
      nx = -ez / len
      nz =  ex / len

      if (nw == 0) exit

      call clip_halfplane(wx(1:nw), wz(1:nw), nw, ox, oz, no, nx, nz, &
                          pvx(ie), pvz(ie))
      if (no == 0) then
        nw = 0
        exit
      end if
      wx(1:no) = ox(1:no)
      wz(1:no) = oz(1:no)
      nw = no
    end do

    if (nw < 3) then
      area = 0.0_dp
    else
      area = polygon_area(wx, wz, nw)
    end if
  end function clip_cell_against_polygon

  ! ------------------------------------------------------------------
  ! Fraction of a line segment [a, b] inside a circle centered at
  ! (cx, cz) with given radius.  Returns fraction in [0, 1].
  !
  ! The segment is horizontal (constant z) for y-face computation, or
  ! vertical (constant y) for z-face computation.  Uses the analytical
  ! chord-length formula.
  ! ------------------------------------------------------------------
  pure function segment_circle_fraction(a, b, coord, is_horizontal, &
      cx, cz, radius) result(frac)
    real(kind=dp), intent(in)  :: a, b          ! segment endpoints along main axis
    real(kind=dp), intent(in)  :: coord          ! perpendicular coordinate
    logical, intent(in)        :: is_horizontal  ! T: segment along y at z=coord
    real(kind=dp), intent(in)  :: cx, cz, radius
    real(kind=dp) :: frac

    real(kind=dp) :: d2, r2, half_chord, t1, t2, s_lo, s_hi
    real(kind=dp) :: dx, dy

    r2 = radius * radius

    if (is_horizontal) then
      ! Segment from (a, coord) to (b, coord); distance to center
      dy = coord - cz
      d2 = dy * dy
    else
      ! Segment from (coord, a) to (coord, b)
      dx = coord - cx
      d2 = dx * dx
    end if

    if (d2 >= r2) then
      frac = 0.0_dp
      return
    end if

    half_chord = sqrt(r2 - d2)

    if (is_horizontal) then
      t1 = cx - half_chord
      t2 = cx + half_chord
    else
      t1 = cz - half_chord
      t2 = cz + half_chord
    end if

    ! Intersection of [a,b] with [t1,t2]
    s_lo = max(a, t1)
    s_hi = min(b, t2)

    if (s_hi <= s_lo) then
      frac = 0.0_dp
    else
      frac = (s_hi - s_lo) / (b - a)
    end if
  end function segment_circle_fraction

  ! ------------------------------------------------------------------
  ! Build y(:) and z(:) cell-centered coordinate arrays.
  ! ------------------------------------------------------------------
  subroutine build_coord_arrays(grid)
    type(spatial_grid), intent(inout) :: grid
    integer :: i

    if (.not. allocated(grid%y)) allocate(grid%y(grid%ny))
    if (.not. allocated(grid%z)) allocate(grid%z(grid%nz))

    do i = 1, grid%ny
      grid%y(i) = (i - 0.5_dp) * grid%dy
    end do
    do i = 1, grid%nz
      grid%z(i) = (i - 0.5_dp) * grid%dz
    end do
  end subroutine build_coord_arrays

  ! ------------------------------------------------------------------
  ! Build the ghost_map for a 2D grid.  For each cell, find the
  ! nearest active neighbor in each of the 4 directions (N=+z, S=-z,
  ! W=-y, E=+y).  For inactive cells (cell_volume==0), search outward
  ! until an active cell is found.  Boundary cells map to themselves.
  !
  ! Directions: 1=N(+z), 2=S(-z), 3=W(-y), 4=E(+y)
  ! ------------------------------------------------------------------
  subroutine build_ghost_map(grid)
    type(spatial_grid), intent(inout) :: grid
    integer :: iy, iz, ij, dy_idx, dz_idx, ny_idx, nz_idx, ngrid
    integer :: search, sy, sz, ty, tz, tij
    logical :: found

    ny_idx = grid%ny
    nz_idx = grid%nz
    ngrid  = ny_idx * nz_idx

    if (.not. allocated(grid%ghost_map)) then
      allocate(grid%ghost_map(ngrid, 4))
    end if

    do iz = 1, nz_idx
      do iy = 1, ny_idx
        ij = flat_idx(ny_idx, iy, iz)

        ! Direction 1: N (+z)
        if (iz < nz_idx) then
          tij = flat_idx(ny_idx, iy, iz + 1)
          grid%ghost_map(ij, 1) = tij
        else
          grid%ghost_map(ij, 1) = ij
        end if

        ! Direction 2: S (-z)
        if (iz > 1) then
          tij = flat_idx(ny_idx, iy, iz - 1)
          grid%ghost_map(ij, 2) = tij
        else
          grid%ghost_map(ij, 2) = ij
        end if

        ! Direction 3: W (-y)
        if (iy > 1) then
          tij = flat_idx(ny_idx, iy - 1, iz)
          grid%ghost_map(ij, 3) = tij
        else
          grid%ghost_map(ij, 3) = ij
        end if

        ! Direction 4: E (+y)
        if (iy < ny_idx) then
          tij = flat_idx(ny_idx, iy + 1, iz)
          grid%ghost_map(ij, 4) = tij
        else
          grid%ghost_map(ij, 4) = ij
        end if

        ! For inactive cells, search outward for nearest active neighbor
        if (grid%cell_volume(ij) < 0.5_dp) then
          do search = 1, max(ny_idx, nz_idx)
            found = .false.
            ! N
            if (.not. found) then
              do sz = 1, min(search, nz_idx - iz)
                tij = flat_idx(ny_idx, iy, iz + sz)
                if (grid%cell_volume(tij) > 0.5_dp) then
                  grid%ghost_map(ij, 1) = tij
                  found = .true.
                  exit
                end if
              end do
            end if
            ! S
            found = .false.
            do sz = 1, min(search, iz - 1)
              tij = flat_idx(ny_idx, iy, iz - sz)
              if (grid%cell_volume(tij) > 0.5_dp) then
                grid%ghost_map(ij, 2) = tij
                found = .true.
                exit
              end if
            end do
            ! W
            found = .false.
            do sy = 1, min(search, iy - 1)
              tij = flat_idx(ny_idx, iy - sy, iz)
              if (grid%cell_volume(tij) > 0.5_dp) then
                grid%ghost_map(ij, 3) = tij
                found = .true.
                exit
              end if
            end do
            ! E
            found = .false.
            do sy = 1, min(search, ny_idx - iy)
              tij = flat_idx(ny_idx, iy + sy, iz)
              if (grid%cell_volume(tij) > 0.5_dp) then
                grid%ghost_map(ij, 4) = tij
                found = .true.
                exit
              end if
            end do
          end do
          ! If still no active neighbor found, map to self
          if (grid%cell_volume(grid%ghost_map(ij, 1)) < 0.5_dp) &
            grid%ghost_map(ij, 1) = ij
          if (grid%cell_volume(grid%ghost_map(ij, 2)) < 0.5_dp) &
            grid%ghost_map(ij, 2) = ij
          if (grid%cell_volume(grid%ghost_map(ij, 3)) < 0.5_dp) &
            grid%ghost_map(ij, 3) = ij
          if (grid%cell_volume(grid%ghost_map(ij, 4)) < 0.5_dp) &
            grid%ghost_map(ij, 4) = ij
        end if
      end do
    end do
  end subroutine build_ghost_map

  ! ==================================================================
  ! Public grid initialization routines
  ! ==================================================================

  ! ------------------------------------------------------------------
  ! Rectangle wire.  All cells are fully active (volume=1, faces=1).
  ! material_id_2d(ny,nz) maps each 2D cell to a material layer.
  ! ------------------------------------------------------------------
  subroutine grid_init_rect(grid, ny, nz, dy, dz, material_id_2d)
    type(spatial_grid), intent(inout) :: grid
    integer, intent(in)               :: ny, nz
    real(kind=dp), intent(in)         :: dy, dz
    integer, intent(in)               :: material_id_2d(ny, nz)

    integer :: iy, iz, ij, ngrid

    grid%ndim = 2
    grid%ny   = ny
    grid%nz   = nz
    grid%dy   = dy
    grid%dz   = dz
    ngrid     = ny * nz

    call build_coord_arrays(grid)

    ! Allocate cut-cell fields
    allocate(grid%material_id(ngrid))
    allocate(grid%cell_volume(ngrid))
    allocate(grid%face_fraction_y(ngrid, 2))
    allocate(grid%face_fraction_z(ngrid, 2))

    do iz = 1, nz
      do iy = 1, ny
        ij = flat_idx(ny, iy, iz)
        grid%material_id(ij)      = material_id_2d(iy, iz)
        grid%cell_volume(ij)      = 1.0_dp
        grid%face_fraction_y(ij, 1) = 1.0_dp
        grid%face_fraction_y(ij, 2) = 1.0_dp
        grid%face_fraction_z(ij, 1) = 1.0_dp
        grid%face_fraction_z(ij, 2) = 1.0_dp
      end do
    end do

    call build_ghost_map(grid)
  end subroutine grid_init_rect

  ! ------------------------------------------------------------------
  ! Circular wire centered at (cx, cz) with given radius.
  ! Uses analytical chord intersection for face fractions and
  ! cell-volume computation via marching-squares sub-sampling.
  ! ------------------------------------------------------------------
  subroutine grid_init_circle(grid, ny, nz, dy, dz, cx, cz, radius, &
      material_id_2d)
    type(spatial_grid), intent(inout) :: grid
    integer, intent(in)               :: ny, nz
    real(kind=dp), intent(in)         :: dy, dz, cx, cz, radius
    integer, intent(in)               :: material_id_2d(ny, nz)

    integer       :: iy, iz, ij, ngrid
    real(kind=dp) :: ylo, yhi, zlo, zhi, r2, dy2, dz2
    real(kind=dp) :: dist_sq
    integer       :: ky, kz, n_inside, n_total

    grid%ndim = 2
    grid%ny   = ny
    grid%nz   = nz
    grid%dy   = dy
    grid%dz   = dz
    ngrid     = ny * nz

    call build_coord_arrays(grid)

    allocate(grid%material_id(ngrid))
    allocate(grid%cell_volume(ngrid))
    allocate(grid%face_fraction_y(ngrid, 2))
    allocate(grid%face_fraction_z(ngrid, 2))

    r2 = radius * radius

    do iz = 1, nz
      do iy = 1, ny
        ij = flat_idx(ny, iy, iz)

        ylo = (iy - 1) * dy
        yhi = iy * dy
        zlo = (iz - 1) * dz
        zhi = iz * dz

        ! Quick rejection/acceptance
        dy2 = max(0.0_dp, ylo - cx)**2
        dy2 = dy2 + max(0.0_dp, cx - yhi)**2
        dz2 = max(0.0_dp, zlo - cz)**2
        dz2 = dz2 + max(0.0_dp, cz - zhi)**2
        dist_sq = dy2 + dz2

        if (dist_sq >= r2) then
          ! Cell entirely outside
          grid%material_id(ij)      = 0
          grid%cell_volume(ij)      = 0.0_dp
          grid%face_fraction_y(ij,:) = 0.0_dp
          grid%face_fraction_z(ij,:) = 0.0_dp
          cycle
        end if

        ! Check if fully inside: all corners within circle
        if ((ylo-cx)**2+(zlo-cz)**2 < r2 .and. &
            (yhi-cx)**2+(zlo-cz)**2 < r2 .and. &
            (ylo-cx)**2+(zhi-cz)**2 < r2 .and. &
            (yhi-cx)**2+(zhi-cz)**2 < r2) then
          ! Fully inside
          grid%material_id(ij)      = material_id_2d(iy, iz)
          grid%cell_volume(ij)      = 1.0_dp
          grid%face_fraction_y(ij,:) = 1.0_dp
          grid%face_fraction_z(ij,:) = 1.0_dp
          cycle
        end if

        ! Partial cell: use sub-sampling for volume
        n_total = MS_N * MS_N
        n_inside = 0
        do kz = 1, MS_N
          do ky = 1, MS_N
            dy2 = (ylo + (ky - 0.5_dp) * dy / MS_N - cx)**2
            dz2 = (zlo + (kz - 0.5_dp) * dz / MS_N - cz)**2
            if (dy2 + dz2 < r2) n_inside = n_inside + 1
          end do
        end do

        grid%cell_volume(ij) = real(n_inside, kind=dp) / real(n_total, kind=dp)

        if (grid%cell_volume(ij) > 0.0_dp) then
          grid%material_id(ij) = material_id_2d(iy, iz)
        else
          grid%material_id(ij) = 0
        end if

        ! Face fractions: left y-face (y = ylo), right y-face (y = yhi)
        grid%face_fraction_y(ij, 1) = segment_circle_fraction( &
          zlo, zhi, ylo, .false., cx, cz, radius)
        grid%face_fraction_y(ij, 2) = segment_circle_fraction( &
          zlo, zhi, yhi, .false., cx, cz, radius)

        ! Face fractions: bottom z-face (z = zlo), top z-face (z = zhi)
        grid%face_fraction_z(ij, 1) = segment_circle_fraction( &
          ylo, yhi, zlo, .true., cx, cz, radius)
        grid%face_fraction_z(ij, 2) = segment_circle_fraction( &
          ylo, yhi, zhi, .true., cx, cz, radius)
      end do
    end do

    call build_ghost_map(grid)
  end subroutine grid_init_circle

  ! ------------------------------------------------------------------
  ! Regular hexagon wire centered at (cx, cz) with given side length.
  ! The hexagon has flat-top orientation: vertices at angles 0, 60, ...,
  ! 300 degrees.  Uses Sutherland-Hodgman polygon clipping for cells.
  ! ------------------------------------------------------------------
  subroutine grid_init_hexagon(grid, ny, nz, dy, dz, cx, cz, side_length, &
      material_id_2d)
    type(spatial_grid), intent(inout) :: grid
    integer, intent(in)               :: ny, nz
    real(kind=dp), intent(in)         :: dy, dz, cx, cz, side_length
    integer, intent(in)               :: material_id_2d(ny, nz)

    real(kind=dp) :: hx(6), hz(6)
    integer       :: i

    ! Build hexagon vertices (flat-top: first vertex to the right)
    do i = 1, 6
      hx(i) = cx + side_length * cos(real(i - 1, kind=dp) * pi_dp / 3.0_dp)
      hz(i) = cz + side_length * sin(real(i - 1, kind=dp) * pi_dp / 3.0_dp)
    end do

    call grid_init_polygon(grid, ny, nz, dy, dz, 6, hx, hz, material_id_2d)
  end subroutine grid_init_hexagon

  ! ------------------------------------------------------------------
  ! Arbitrary polygon wire from vertex list.
  ! Uses Sutherland-Hodgman clip for each cell to compute cut-cell
  ! fractions.
  ! ------------------------------------------------------------------
  subroutine grid_init_polygon(grid, ny, nz, dy, dz, nvert, vx, vz, &
      material_id_2d)
    type(spatial_grid), intent(inout) :: grid
    integer, intent(in)               :: ny, nz, nvert
    real(kind=dp), intent(in)         :: dy, dz
    real(kind=dp), intent(in)         :: vx(nvert), vz(nvert)
    integer, intent(in)               :: material_id_2d(ny, nz)

    integer       :: iy, iz, ij, ngrid, ic
    real(kind=dp) :: ylo, yhi, zlo, zhi, cell_area, full_area
    real(kind=dp) :: ry(4), rz(4)
    logical       :: all_inside

    grid%ndim = 2
    grid%ny   = ny
    grid%nz   = nz
    grid%dy   = dy
    grid%dz   = dz
    ngrid     = ny * nz

    call build_coord_arrays(grid)

    allocate(grid%material_id(ngrid))
    allocate(grid%cell_volume(ngrid))
    allocate(grid%face_fraction_y(ngrid, 2))
    allocate(grid%face_fraction_z(ngrid, 2))

    full_area = dy * dz

    do iz = 1, nz
      do iy = 1, ny
        ij = flat_idx(ny, iy, iz)

        ylo = (iy - 1) * dy
        yhi = iy * dy
        zlo = (iz - 1) * dz
        zhi = iz * dz

        ! Quick test: check if all 4 corners are inside the polygon
        ! using a point-in-polygon test
        all_inside = .true.
        ry(1) = ylo; rz(1) = zlo
        ry(2) = yhi; rz(2) = zlo
        ry(3) = yhi; rz(3) = zhi
        ry(4) = ylo; rz(4) = zhi
        do ic = 1, 4
          if (.not. point_in_polygon(ry(ic), rz(ic), vx, vz, nvert)) then
            all_inside = .false.
            exit
          end if
        end do

        if (all_inside) then
          grid%material_id(ij)      = material_id_2d(iy, iz)
          grid%cell_volume(ij)      = 1.0_dp
          grid%face_fraction_y(ij,:) = 1.0_dp
          grid%face_fraction_z(ij,:) = 1.0_dp
          cycle
        end if

        ! Clip cell rectangle against polygon
        cell_area = clip_cell_against_polygon(ylo, yhi, zlo, zhi, &
          vx, vz, nvert)

        grid%cell_volume(ij) = cell_area / full_area

        if (grid%cell_volume(ij) > 0.0_dp) then
          grid%material_id(ij) = material_id_2d(iy, iz)
        else
          grid%material_id(ij) = 0
        end if

        ! Face fractions: clip line segment against polygon
        ! Left y-face:  segment from (ylo, zlo) to (ylo, zhi) -- vertical at y=ylo
        grid%face_fraction_y(ij, 1) = segment_polygon_fraction( &
          zlo, zhi, ylo, .false., vx, vz, nvert)
        ! Right y-face: segment at y=yhi
        grid%face_fraction_y(ij, 2) = segment_polygon_fraction( &
          zlo, zhi, yhi, .false., vx, vz, nvert)
        ! Bottom z-face: segment from (ylo, zlo) to (yhi, zlo) -- horizontal at z=zlo
        grid%face_fraction_z(ij, 1) = segment_polygon_fraction( &
          ylo, yhi, zlo, .true., vx, vz, nvert)
        ! Top z-face: segment at z=zhi
        grid%face_fraction_z(ij, 2) = segment_polygon_fraction( &
          ylo, yhi, zhi, .true., vx, vz, nvert)
      end do
    end do

    call build_ghost_map(grid)
  end subroutine grid_init_polygon

  ! ------------------------------------------------------------------
  ! Point-in-polygon test (ray casting algorithm).
  ! Returns .true. if point (py, pz) is inside the polygon defined by
  ! vertices (vx, vz).
  ! ------------------------------------------------------------------
  pure function point_in_polygon(py, pz, vx, vz, n) result(inside)
    real(kind=dp), intent(in)  :: py, pz
    integer, intent(in)        :: n
    real(kind=dp), intent(in)  :: vx(n), vz(n)
    logical :: inside

    integer       :: i, j
    logical       :: c

    c = .false.
    j = n
    do i = 1, n
      if (((vz(i) > pz) .neqv. (vz(j) > pz)) .and. &
          (py < (vx(j) - vx(i)) * (pz - vz(i)) / (vz(j) - vz(i)) + vx(i))) then
        c = .not. c
      end if
      j = i
    end do
    inside = c
  end function point_in_polygon

  ! ------------------------------------------------------------------
  ! Fraction of a line segment inside a polygon.
  !
  ! For a horizontal segment at z=coord from a to b:
  !   finds all intersection t-parameters, sorts them, and sums up
  !   the portions inside the polygon.
  ! For a vertical segment at y=coord from a to b:
  !   similar treatment.
  ! ------------------------------------------------------------------
  pure function segment_polygon_fraction(a, b, coord, is_horizontal, &
      vx, vz, n) result(frac)
    real(kind=dp), intent(in)  :: a, b
    real(kind=dp), intent(in)  :: coord
    logical, intent(in)        :: is_horizontal
    integer, intent(in)        :: n
    real(kind=dp), intent(in)  :: vx(n), vz(n)
    real(kind=dp) :: frac

    ! We use a sub-sampling approach for robustness.
    ! Number of sample points along the segment.
    integer, parameter :: ns = 32
    integer       :: k, n_inside
    real(kind=dp) :: t, py, pz, dt

    dt = (b - a) / real(ns, kind=dp)
    n_inside = 0
    do k = 1, ns
      t = a + (k - 0.5_dp) * dt
      if (is_horizontal) then
        py = t
        pz = coord
      else
        py = coord
        pz = t
      end if
      if (point_in_polygon(py, pz, vx, vz, n)) n_inside = n_inside + 1
    end do
    frac = real(n_inside, kind=dp) / real(ns, kind=dp)
  end function segment_polygon_fraction

  ! ==================================================================
  ! Utility routines
  ! ==================================================================

  ! ------------------------------------------------------------------
  ! Build flattened coords(:,:) from y(:) and z(:).
  ! coords(1, ij) = y(iy), coords(2, ij) = z(iz)
  ! with column-major ordering: ij = (iz-1)*ny + iy
  ! ------------------------------------------------------------------
  subroutine grid_build_coords(grid)
    type(spatial_grid), intent(inout) :: grid

    integer :: iy, iz, ij, ngrid

    ngrid = grid%ny * grid%nz

    if (.not. allocated(grid%y) .or. .not. allocated(grid%z)) return

    if (allocated(grid%coords)) deallocate(grid%coords)
    allocate(grid%coords(2, ngrid))

    do iz = 1, grid%nz
      do iy = 1, grid%ny
        ij = flat_idx(grid%ny, iy, iz)
        grid%coords(1, ij) = grid%y(iy)
        grid%coords(2, ij) = grid%z(iz)
      end do
    end do
  end subroutine grid_build_coords

  ! ------------------------------------------------------------------
  ! Deallocate all allocatable fields in a spatial_grid.
  ! ------------------------------------------------------------------
  subroutine grid_free(grid)
    type(spatial_grid), intent(inout) :: grid

    if (allocated(grid%y))              deallocate(grid%y)
    if (allocated(grid%z))              deallocate(grid%z)
    if (allocated(grid%coords))         deallocate(grid%coords)
    if (allocated(grid%material_id))    deallocate(grid%material_id)
    if (allocated(grid%cell_volume))    deallocate(grid%cell_volume)
    if (allocated(grid%face_fraction_y)) deallocate(grid%face_fraction_y)
    if (allocated(grid%face_fraction_z)) deallocate(grid%face_fraction_z)
    if (allocated(grid%ghost_map))      deallocate(grid%ghost_map)

    grid%ndim = 0
    grid%ny   = 1
    grid%nz   = 1
    grid%dy   = 0.0_dp
    grid%dz   = 0.0_dp
  end subroutine grid_free

end module geometry
