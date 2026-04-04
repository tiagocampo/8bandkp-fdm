module geometry

  ! ==============================================================================
  ! Cut-cell immersed boundary geometry for 2D quantum wire cross-sections.
  !
  ! Computes cell-volume and face-fraction arrays for arbitrary wire shapes
  ! (rectangle, circle, hexagon, polygon) on a Cartesian FD grid.
  ! All shapes are defined in the x-y plane.  Grid coordinates are cell-centered:
  !   x(i) = (i - 0.5) * dx,   y(j) = (j - 0.5) * dy
  !
  ! Uses analytical intersection wherever possible; falls back to marching-squares
  ! sub-sampling for general polygons.
  !
  ! Public procedures:
  !   grid_init_rect    -- rectangle (trivial, all fractions = 1)
  !   grid_init_circle  -- circular wire (analytical chord intersection)
  !   grid_init_hexagon -- regular hexagon (6 half-plane intersections)
  !   grid_init_polygon -- arbitrary polygon (Sutherland-Hodgman clip)
  !   grid_build_coords -- build flattened coords(:,:) from x(:), z(:)
  !   grid_free         -- deallocate all allocatable fields
  ! ==============================================================================

  use definitions, only: dp, spatial_grid, pi_dp

  implicit none

  private
  public :: grid_init_rect, grid_init_circle, grid_init_hexagon
  public :: grid_init_polygon, grid_build_coords, grid_free
  public :: init_wire_from_config

  ! Number of sub-sample points per cell edge for marching-squares integration
  integer, parameter :: MS_N = 16

contains

  ! ==================================================================
  ! Internal helpers
  ! ==================================================================

  ! ------------------------------------------------------------------
  ! Map 2D (ix, iy) to flattened column-major index: ij = (iy-1)*nx + ix
  ! ------------------------------------------------------------------
  pure function flat_idx(nx, ix, iy) result(ij)
    integer, intent(in) :: nx, ix, iy
    integer :: ij
    ij = (iy - 1) * nx + ix
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
  pure function clip_cell_against_polygon(xlo, xhi, ylo, yhi, &
      pvx, pvy, npoly) result(area)
    real(kind=dp), intent(in)  :: xlo, xhi, ylo, yhi
    integer, intent(in)        :: npoly
    real(kind=dp), intent(in)  :: pvx(npoly), pvy(npoly)
    real(kind=dp) :: area

    real(kind=dp) :: wx(8), wy(8)
    real(kind=dp) :: ox(8), oy(8)
    integer       :: nw, no, ie, ip
    real(kind=dp) :: ex, ey, nx, ny, len

    ! Start with the cell rectangle (CCW)
    wx(1) = xlo; wy(1) = ylo
    wx(2) = xhi; wy(2) = ylo
    wx(3) = xhi; wy(3) = yhi
    wx(4) = xlo; wy(4) = yhi
    nw = 4

    ! Clip against each edge of the polygon
    do ie = 1, npoly
      ip = modulo(ie, npoly) + 1
      ! Edge from vertex ie to ip
      ex = pvx(ip) - pvx(ie)
      ey = pvy(ip) - pvy(ie)
      len = sqrt(ex * ex + ey * ey)
      if (len < 1.0e-30_dp) cycle
      ! Inward normal (left-hand normal for CCW polygon points inward)
      nx = -ey / len
      ny =  ex / len

      if (nw == 0) exit

      call clip_halfplane(wx(1:nw), wy(1:nw), nw, ox, oy, no, nx, ny, &
                          pvx(ie), pvy(ie))
      if (no == 0) then
        nw = 0
        exit
      end if
      wx(1:no) = ox(1:no)
      wy(1:no) = oy(1:no)
      nw = no
    end do

    if (nw < 3) then
      area = 0.0_dp
    else
      area = polygon_area(wx, wy, nw)
    end if
  end function clip_cell_against_polygon

  ! ------------------------------------------------------------------
  ! Fraction of a line segment [a, b] inside a circle centered at
  ! (cx, cy) with given radius.  Returns fraction in [0, 1].
  !
  ! The segment is horizontal (constant y) for x-face computation, or
  ! vertical (constant x) for y-face computation.  Uses the analytical
  ! chord-length formula.
  ! ------------------------------------------------------------------
  pure function segment_circle_fraction(a, b, coord, is_horizontal, &
      cx, cy, radius) result(frac)
    real(kind=dp), intent(in)  :: a, b          ! segment endpoints along main axis
    real(kind=dp), intent(in)  :: coord          ! perpendicular coordinate
    logical, intent(in)        :: is_horizontal  ! T: segment along x at y=coord
    real(kind=dp), intent(in)  :: cx, cy, radius
    real(kind=dp) :: frac

    real(kind=dp) :: d2, r2, half_chord, t1, t2, s_lo, s_hi
    real(kind=dp) :: ddx, ddy

    r2 = radius * radius

    if (is_horizontal) then
      ! Segment from (a, coord) to (b, coord); distance to center
      ddy = coord - cy
      d2 = ddy * ddy
    else
      ! Segment from (coord, a) to (coord, b)
      ddx = coord - cx
      d2 = ddx * ddx
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
      t1 = cy - half_chord
      t2 = cy + half_chord
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
  ! Build x(:) and z(:) cell-centered coordinate arrays.
  ! ------------------------------------------------------------------
  subroutine build_coord_arrays(grid)
    type(spatial_grid), intent(inout) :: grid
    integer :: i

    if (.not. allocated(grid%x)) allocate(grid%x(grid%nx))
    if (.not. allocated(grid%z)) allocate(grid%z(grid%ny))

    do i = 1, grid%nx
      grid%x(i) = (i - 0.5_dp) * grid%dx
    end do
    do i = 1, grid%ny
      grid%z(i) = (i - 0.5_dp) * grid%dy
    end do
  end subroutine build_coord_arrays

  ! ------------------------------------------------------------------
  ! Build the ghost_map for a 2D grid.  For each cell, find the
  ! nearest active neighbor in each of the 4 directions (N=+y, S=-y,
  ! W=-x, E=+x).  For inactive cells (cell_volume==0), search outward
  ! until an active cell is found.  Boundary cells map to themselves.
  !
  ! Directions: 1=N(+y), 2=S(-y), 3=W(-x), 4=E(+x)
  ! ------------------------------------------------------------------
  subroutine build_ghost_map(grid)
    type(spatial_grid), intent(inout) :: grid
    integer :: ix, iy, ij, ngrid
    integer :: nx_idx, ny_idx
    integer :: search, sx, sy, tx, ty, tij
    logical :: found

    nx_idx = grid%nx
    ny_idx = grid%ny
    ngrid  = nx_idx * ny_idx

    if (.not. allocated(grid%ghost_map)) then
      allocate(grid%ghost_map(ngrid, 4))
    end if

    do iy = 1, ny_idx
      do ix = 1, nx_idx
        ij = flat_idx(nx_idx, ix, iy)

        ! Direction 1: N (+y)
        if (iy < ny_idx) then
          tij = flat_idx(nx_idx, ix, iy + 1)
          grid%ghost_map(ij, 1) = tij
        else
          grid%ghost_map(ij, 1) = ij
        end if

        ! Direction 2: S (-y)
        if (iy > 1) then
          tij = flat_idx(nx_idx, ix, iy - 1)
          grid%ghost_map(ij, 2) = tij
        else
          grid%ghost_map(ij, 2) = ij
        end if

        ! Direction 3: W (-x)
        if (ix > 1) then
          tij = flat_idx(nx_idx, ix - 1, iy)
          grid%ghost_map(ij, 3) = tij
        else
          grid%ghost_map(ij, 3) = ij
        end if

        ! Direction 4: E (+x)
        if (ix < nx_idx) then
          tij = flat_idx(nx_idx, ix + 1, iy)
          grid%ghost_map(ij, 4) = tij
        else
          grid%ghost_map(ij, 4) = ij
        end if

        ! For inactive cells, search outward for nearest active neighbor
        if (grid%cell_volume(ij) < 0.5_dp) then
          do search = 1, max(nx_idx, ny_idx)
            found = .false.
            ! N
            if (.not. found) then
              do sy = 1, min(search, ny_idx - iy)
                tij = flat_idx(nx_idx, ix, iy + sy)
                if (grid%cell_volume(tij) > 0.5_dp) then
                  grid%ghost_map(ij, 1) = tij
                  found = .true.
                  exit
                end if
              end do
            end if
            ! S
            found = .false.
            do sy = 1, min(search, iy - 1)
              tij = flat_idx(nx_idx, ix, iy - sy)
              if (grid%cell_volume(tij) > 0.5_dp) then
                grid%ghost_map(ij, 2) = tij
                found = .true.
                exit
              end if
            end do
            ! W
            found = .false.
            do sx = 1, min(search, ix - 1)
              tij = flat_idx(nx_idx, ix - sx, iy)
              if (grid%cell_volume(tij) > 0.5_dp) then
                grid%ghost_map(ij, 3) = tij
                found = .true.
                exit
              end if
            end do
            ! E
            found = .false.
            do sx = 1, min(search, nx_idx - ix)
              tij = flat_idx(nx_idx, ix + sx, iy)
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
  ! material_id_2d(nx,ny) maps each 2D cell to a material layer.
  ! ------------------------------------------------------------------
  subroutine grid_init_rect(grid, nx, ny, dx, dy, material_id_2d)
    type(spatial_grid), intent(inout) :: grid
    integer, intent(in)               :: nx, ny
    real(kind=dp), intent(in)         :: dx, dy
    integer, intent(in)               :: material_id_2d(nx, ny)

    integer :: ix, iy, ij, ngrid

    grid%ndim = 2
    grid%nx   = nx
    grid%ny   = ny
    grid%dx   = dx
    grid%dy   = dy
    ngrid     = nx * ny

    call build_coord_arrays(grid)

    ! Allocate cut-cell fields
    allocate(grid%material_id(ngrid))
    allocate(grid%cell_volume(ngrid))
    allocate(grid%face_fraction_x(ngrid, 2))
    allocate(grid%face_fraction_y(ngrid, 2))

    do iy = 1, ny
      do ix = 1, nx
        ij = flat_idx(nx, ix, iy)
        grid%material_id(ij)        = material_id_2d(ix, iy)
        grid%cell_volume(ij)        = 1.0_dp
        grid%face_fraction_x(ij, 1) = 1.0_dp
        grid%face_fraction_x(ij, 2) = 1.0_dp
        grid%face_fraction_y(ij, 1) = 1.0_dp
        grid%face_fraction_y(ij, 2) = 1.0_dp
      end do
    end do

    call build_ghost_map(grid)
  end subroutine grid_init_rect

  ! ------------------------------------------------------------------
  ! Circular wire centered at (cx, cy) with given radius.
  ! Uses analytical chord intersection for face fractions and
  ! cell-volume computation via marching-squares sub-sampling.
  ! ------------------------------------------------------------------
  subroutine grid_init_circle(grid, nx, ny, dx, dy, cx, cy, radius, &
      material_id_2d)
    type(spatial_grid), intent(inout) :: grid
    integer, intent(in)               :: nx, ny
    real(kind=dp), intent(in)         :: dx, dy, cx, cy, radius
    integer, intent(in)               :: material_id_2d(nx, ny)

    integer       :: ix, iy, ij, ngrid
    real(kind=dp) :: xlo, xhi, ylo, yhi, r2, dx2, dy2
    real(kind=dp) :: dist_sq
    integer       :: kx, ky, n_inside, n_total

    grid%ndim = 2
    grid%nx   = nx
    grid%ny   = ny
    grid%dx   = dx
    grid%dy   = dy
    ngrid     = nx * ny

    call build_coord_arrays(grid)

    allocate(grid%material_id(ngrid))
    allocate(grid%cell_volume(ngrid))
    allocate(grid%face_fraction_x(ngrid, 2))
    allocate(grid%face_fraction_y(ngrid, 2))

    r2 = radius * radius

    do iy = 1, ny
      do ix = 1, nx
        ij = flat_idx(nx, ix, iy)

        xlo = (ix - 1) * dx
        xhi = ix * dx
        ylo = (iy - 1) * dy
        yhi = iy * dy

        ! Quick rejection/acceptance
        dx2 = max(0.0_dp, xlo - cx)**2
        dx2 = dx2 + max(0.0_dp, cx - xhi)**2
        dy2 = max(0.0_dp, ylo - cy)**2
        dy2 = dy2 + max(0.0_dp, cy - yhi)**2
        dist_sq = dx2 + dy2

        if (dist_sq >= r2) then
          ! Cell entirely outside
          grid%material_id(ij)      = 0
          grid%cell_volume(ij)      = 0.0_dp
          grid%face_fraction_x(ij,:) = 0.0_dp
          grid%face_fraction_y(ij,:) = 0.0_dp
          cycle
        end if

        ! Check if fully inside: all corners within circle
        if ((xlo-cx)**2+(ylo-cy)**2 < r2 .and. &
            (xhi-cx)**2+(ylo-cy)**2 < r2 .and. &
            (xlo-cx)**2+(yhi-cy)**2 < r2 .and. &
            (xhi-cx)**2+(yhi-cy)**2 < r2) then
          ! Fully inside
          grid%material_id(ij)      = material_id_2d(ix, iy)
          grid%cell_volume(ij)      = 1.0_dp
          grid%face_fraction_x(ij,:) = 1.0_dp
          grid%face_fraction_y(ij,:) = 1.0_dp
          cycle
        end if

        ! Partial cell: use sub-sampling for volume
        n_total = MS_N * MS_N
        n_inside = 0
        do ky = 1, MS_N
          do kx = 1, MS_N
            dx2 = (xlo + (kx - 0.5_dp) * dx / MS_N - cx)**2
            dy2 = (ylo + (ky - 0.5_dp) * dy / MS_N - cy)**2
            if (dx2 + dy2 < r2) n_inside = n_inside + 1
          end do
        end do

        grid%cell_volume(ij) = real(n_inside, kind=dp) / real(n_total, kind=dp)

        if (grid%cell_volume(ij) > 0.0_dp) then
          grid%material_id(ij) = material_id_2d(ix, iy)
        else
          grid%material_id(ij) = 0
        end if

        ! Face fractions: left x-face (x = xlo), right x-face (x = xhi)
        grid%face_fraction_x(ij, 1) = segment_circle_fraction( &
          ylo, yhi, xlo, .false., cx, cy, radius)
        grid%face_fraction_x(ij, 2) = segment_circle_fraction( &
          ylo, yhi, xhi, .false., cx, cy, radius)

        ! Face fractions: bottom y-face (y = ylo), top y-face (y = yhi)
        grid%face_fraction_y(ij, 1) = segment_circle_fraction( &
          xlo, xhi, ylo, .true., cx, cy, radius)
        grid%face_fraction_y(ij, 2) = segment_circle_fraction( &
          xlo, xhi, yhi, .true., cx, cy, radius)
      end do
    end do

    call build_ghost_map(grid)
  end subroutine grid_init_circle

  ! ------------------------------------------------------------------
  ! Regular hexagon wire centered at (cx, cy) with given side length.
  ! The hexagon has flat-top orientation: vertices at angles 0, 60, ...,
  ! 300 degrees.  Uses Sutherland-Hodgman polygon clipping for cells.
  ! ------------------------------------------------------------------
  subroutine grid_init_hexagon(grid, nx, ny, dx, dy, cx, cy, side_length, &
      material_id_2d)
    type(spatial_grid), intent(inout) :: grid
    integer, intent(in)               :: nx, ny
    real(kind=dp), intent(in)         :: dx, dy, cx, cy, side_length
    integer, intent(in)               :: material_id_2d(nx, ny)

    real(kind=dp) :: hx(6), hy(6)
    integer       :: i

    ! Build hexagon vertices (flat-top: first vertex to the right)
    do i = 1, 6
      hx(i) = cx + side_length * cos(real(i - 1, kind=dp) * pi_dp / 3.0_dp)
      hy(i) = cy + side_length * sin(real(i - 1, kind=dp) * pi_dp / 3.0_dp)
    end do

    call grid_init_polygon(grid, nx, ny, dx, dy, 6, hx, hy, material_id_2d)
  end subroutine grid_init_hexagon

  ! ------------------------------------------------------------------
  ! Arbitrary polygon wire from vertex list.
  ! Uses Sutherland-Hodgman clip for each cell to compute cut-cell
  ! fractions.
  ! ------------------------------------------------------------------
  subroutine grid_init_polygon(grid, nx, ny, dx, dy, nvert, vx, vy, &
      material_id_2d)
    type(spatial_grid), intent(inout) :: grid
    integer, intent(in)               :: nx, ny, nvert
    real(kind=dp), intent(in)         :: dx, dy
    real(kind=dp), intent(in)         :: vx(nvert), vy(nvert)
    integer, intent(in)               :: material_id_2d(nx, ny)

    integer       :: ix, iy, ij, ngrid, ic
    real(kind=dp) :: xlo, xhi, ylo, yhi, cell_area, full_area
    real(kind=dp) :: rx(4), ry(4)
    logical       :: all_inside

    grid%ndim = 2
    grid%nx   = nx
    grid%ny   = ny
    grid%dx   = dx
    grid%dy   = dy
    ngrid     = nx * ny

    call build_coord_arrays(grid)

    allocate(grid%material_id(ngrid))
    allocate(grid%cell_volume(ngrid))
    allocate(grid%face_fraction_x(ngrid, 2))
    allocate(grid%face_fraction_y(ngrid, 2))

    full_area = dx * dy

    do iy = 1, ny
      do ix = 1, nx
        ij = flat_idx(nx, ix, iy)

        xlo = (ix - 1) * dx
        xhi = ix * dx
        ylo = (iy - 1) * dy
        yhi = iy * dy

        ! Quick test: check if all 4 corners are inside the polygon
        ! using a point-in-polygon test
        all_inside = .true.
        rx(1) = xlo; ry(1) = ylo
        rx(2) = xhi; ry(2) = ylo
        rx(3) = xhi; ry(3) = yhi
        rx(4) = xlo; ry(4) = yhi
        do ic = 1, 4
          if (.not. point_in_polygon(rx(ic), ry(ic), vx, vy, nvert)) then
            all_inside = .false.
            exit
          end if
        end do

        if (all_inside) then
          grid%material_id(ij)      = material_id_2d(ix, iy)
          grid%cell_volume(ij)      = 1.0_dp
          grid%face_fraction_x(ij,:) = 1.0_dp
          grid%face_fraction_y(ij,:) = 1.0_dp
          cycle
        end if

        ! Clip cell rectangle against polygon
        cell_area = clip_cell_against_polygon(xlo, xhi, ylo, yhi, &
          vx, vy, nvert)

        grid%cell_volume(ij) = cell_area / full_area

        if (grid%cell_volume(ij) > 0.0_dp) then
          grid%material_id(ij) = material_id_2d(ix, iy)
        else
          grid%material_id(ij) = 0
        end if

        ! Face fractions: clip line segment against polygon
        ! Left x-face:  segment from (xlo, ylo) to (xlo, yhi) -- vertical at x=xlo
        grid%face_fraction_x(ij, 1) = segment_polygon_fraction( &
          ylo, yhi, xlo, .false., vx, vy, nvert)
        ! Right x-face: segment at x=xhi
        grid%face_fraction_x(ij, 2) = segment_polygon_fraction( &
          ylo, yhi, xhi, .false., vx, vy, nvert)
        ! Bottom y-face: segment from (xlo, ylo) to (xhi, ylo) -- horizontal at y=ylo
        grid%face_fraction_y(ij, 1) = segment_polygon_fraction( &
          xlo, xhi, ylo, .true., vx, vy, nvert)
        ! Top y-face: segment at y=yhi
        grid%face_fraction_y(ij, 2) = segment_polygon_fraction( &
          xlo, xhi, yhi, .true., vx, vy, nvert)
      end do
    end do

    call build_ghost_map(grid)
  end subroutine grid_init_polygon

  ! ------------------------------------------------------------------
  ! Point-in-polygon test (ray casting algorithm).
  ! Returns .true. if point (px, py) is inside the polygon defined by
  ! vertices (vx, vy).
  ! ------------------------------------------------------------------
  pure function point_in_polygon(px, py, vx, vy, n) result(inside)
    real(kind=dp), intent(in)  :: px, py
    integer, intent(in)        :: n
    real(kind=dp), intent(in)  :: vx(n), vy(n)
    logical :: inside

    integer       :: i, j
    logical       :: c

    c = .false.
    j = n
    do i = 1, n
      if (((vy(i) > py) .neqv. (vy(j) > py)) .and. &
          (px < (vx(j) - vx(i)) * (py - vy(i)) / (vy(j) - vy(i)) + vx(i))) then
        c = .not. c
      end if
      j = i
    end do
    inside = c
  end function point_in_polygon

  ! ------------------------------------------------------------------
  ! Fraction of a line segment inside a polygon.
  !
  ! For a horizontal segment at y=coord from a to b:
  !   finds all intersection t-parameters, sorts them, and sums up
  !   the portions inside the polygon.
  ! For a vertical segment at x=coord from a to b:
  !   similar treatment.
  ! ------------------------------------------------------------------
  pure function segment_polygon_fraction(a, b, coord, is_horizontal, &
      vx, vy, n) result(frac)
    real(kind=dp), intent(in)  :: a, b
    real(kind=dp), intent(in)  :: coord
    logical, intent(in)        :: is_horizontal
    integer, intent(in)        :: n
    real(kind=dp), intent(in)  :: vx(n), vy(n)
    real(kind=dp) :: frac

    ! We use a sub-sampling approach for robustness.
    ! Number of sample points along the segment.
    integer, parameter :: ns = 32
    integer       :: k, n_inside
    real(kind=dp) :: t, px, py, dt

    dt = (b - a) / real(ns, kind=dp)
    n_inside = 0
    do k = 1, ns
      t = a + (k - 0.5_dp) * dt
      if (is_horizontal) then
        px = t
        py = coord
      else
        px = coord
        py = t
      end if
      if (point_in_polygon(px, py, vx, vy, n)) n_inside = n_inside + 1
    end do
    frac = real(n_inside, kind=dp) / real(ns, kind=dp)
  end function segment_polygon_fraction

  ! ==================================================================
  ! High-level wire initialization from simulation_config
  ! ==================================================================

  ! ------------------------------------------------------------------
  ! Initialize wire geometry from a simulation_config.
  !
  ! Expects cfg%grid%ndim/nx/ny/dx/dy to already be set (done by
  ! init_grid_from_config).  This routine:
  !   1. Computes the grid center (midpoint of the x-y domain)
  !   2. Builds material_id_2d(nx, ny) from regions using 2D distance
  !   3. Dispatches to the appropriate grid_init_* routine
  !   4. Builds the flattened coords(:,:) array
  ! ------------------------------------------------------------------
  subroutine init_wire_from_config(cfg)
    use definitions, only: simulation_config, region_spec, dp
    type(simulation_config), intent(inout) :: cfg

    integer :: nx, ny, ix, iy, i, ngrid
    real(kind=dp) :: dx, dy, cx, cy, dist
    integer, allocatable :: mat2d(:,:)

    nx = cfg%grid%nx
    ny = cfg%grid%ny
    dx = cfg%grid%dx
    dy = cfg%grid%dy
    ngrid = nx * ny

    ! Grid center (center of the domain [0, nx*dx] x [0, ny*dy])
    cx = 0.5_dp * real(nx, kind=dp) * dx
    cy = 0.5_dp * real(ny, kind=dp) * dy

    ! Build material_id_2d from regions using 2D distance from center
    allocate(mat2d(nx, ny))
    mat2d = 0
    if (cfg%numRegions > 0) then
      do iy = 1, ny
        do ix = 1, nx
          ! Cell-centered coordinate
          dist = sqrt(((ix - 0.5_dp) * dx - cx)**2 + &
                      ((iy - 0.5_dp) * dy - cy)**2)
          do i = 1, cfg%numRegions
            if (cfg%regions(i)%inner <= dist .and. dist <= cfg%regions(i)%outer) then
              mat2d(ix, iy) = i
              exit
            end if
          end do
        end do
      end do
    end if

    ! Dispatch to the appropriate geometry initializer
    select case (trim(cfg%wire_geom%shape))
    case ('rectangle')
      call grid_init_rect(cfg%grid, nx, ny, dx, dy, mat2d)

    case ('circle')
      call grid_init_circle(cfg%grid, nx, ny, dx, dy, &
        cx, cy, cfg%wire_geom%radius, mat2d)

    case ('hexagon')
      call grid_init_hexagon(cfg%grid, nx, ny, dx, dy, &
        cx, cy, cfg%wire_geom%radius, mat2d)

    case ('polygon')
      call grid_init_polygon(cfg%grid, nx, ny, dx, dy, &
        cfg%wire_geom%nverts, &
        cfg%wire_geom%verts(1, 1:cfg%wire_geom%nverts), &
        cfg%wire_geom%verts(2, 1:cfg%wire_geom%nverts), mat2d)

    case default
      print *, 'Error: Unknown wire shape in init_wire_from_config: ', &
        trim(cfg%wire_geom%shape)
      stop 1
    end select

    ! Build the flattened coords(:,:) array for downstream use
    call grid_build_coords(cfg%grid)

    deallocate(mat2d)
  end subroutine init_wire_from_config

  ! ==================================================================
  ! Utility routines
  ! ==================================================================

  ! ------------------------------------------------------------------
  ! Build flattened coords(:,:) from x(:) and z(:).
  ! coords(1, ij) = x(ix), coords(2, ij) = y(iy) (stored in z(:))
  ! with column-major ordering: ij = (iy-1)*nx + ix
  ! ------------------------------------------------------------------
  subroutine grid_build_coords(grid)
    type(spatial_grid), intent(inout) :: grid

    integer :: ix, iy, ij, ngrid

    ngrid = grid%nx * grid%ny

    if (.not. allocated(grid%x) .or. .not. allocated(grid%z)) return

    if (allocated(grid%coords)) deallocate(grid%coords)
    allocate(grid%coords(2, ngrid))

    do iy = 1, grid%ny
      do ix = 1, grid%nx
        ij = flat_idx(grid%nx, ix, iy)
        grid%coords(1, ij) = grid%x(ix)
        grid%coords(2, ij) = grid%z(iy)
      end do
    end do
  end subroutine grid_build_coords

  ! ------------------------------------------------------------------
  ! Deallocate all allocatable fields in a spatial_grid.
  ! ------------------------------------------------------------------
  subroutine grid_free(grid)
    type(spatial_grid), intent(inout) :: grid

    if (allocated(grid%x))              deallocate(grid%x)
    if (allocated(grid%z))              deallocate(grid%z)
    if (allocated(grid%coords))         deallocate(grid%coords)
    if (allocated(grid%material_id))    deallocate(grid%material_id)
    if (allocated(grid%cell_volume))    deallocate(grid%cell_volume)
    if (allocated(grid%face_fraction_x)) deallocate(grid%face_fraction_x)
    if (allocated(grid%face_fraction_y)) deallocate(grid%face_fraction_y)
    if (allocated(grid%ghost_map))      deallocate(grid%ghost_map)

    grid%ndim = 0
    grid%nx   = 1
    grid%ny   = 1
    grid%dx   = 0.0_dp
    grid%dy   = 0.0_dp
  end subroutine grid_free

end module geometry
