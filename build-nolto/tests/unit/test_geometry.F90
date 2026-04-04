module test_geometry
  use funit
  use definitions
  use geometry
  implicit none

contains

  ! ==================================================================
  ! Rectangle tests
  ! ==================================================================

  !@test
  subroutine test_rect_all_active()
    type(spatial_grid) :: grid
    integer, parameter :: nx = 4, ny = 5
    real(kind=dp), parameter :: dx = 1.0_dp, dy = 2.0_dp
    integer :: mat2d(nx, ny)
    integer :: ix, iy, ij

    mat2d = 1
    call grid_init_rect(grid, nx, ny, dx, dy, mat2d)

#line 24 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(2, grid%ndim, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 24) )
  if (anyExceptions()) return
#line 25 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 25 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(nx, grid%nx, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 25) )
  if (anyExceptions()) return
#line 26 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 26 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(ny, grid%ny, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 26) )
  if (anyExceptions()) return
#line 27 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 27 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertTrue(allocated(grid%cell_volume), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 27) )
  if (anyExceptions()) return
#line 28 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 28 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertTrue(allocated(grid%face_fraction_x), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 28) )
  if (anyExceptions()) return
#line 29 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 29 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertTrue(allocated(grid%face_fraction_y), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 29) )
  if (anyExceptions()) return
#line 30 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 30 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertTrue(allocated(grid%ghost_map), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 30) )
  if (anyExceptions()) return
#line 31 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"

    do iy = 1, ny
      do ix = 1, nx
        ij = (iy - 1) * nx + ix
#line 35 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1.0_dp, grid%cell_volume(ij), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 35) )
  if (anyExceptions()) return
#line 36 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 36 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1.0_dp, grid%face_fraction_x(ij, 1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 36) )
  if (anyExceptions()) return
#line 37 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 37 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1.0_dp, grid%face_fraction_x(ij, 2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 37) )
  if (anyExceptions()) return
#line 38 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 38 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1.0_dp, grid%face_fraction_y(ij, 1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 38) )
  if (anyExceptions()) return
#line 39 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 39 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1.0_dp, grid%face_fraction_y(ij, 2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 39) )
  if (anyExceptions()) return
#line 40 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 40 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1, grid%material_id(ij), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 40) )
  if (anyExceptions()) return
#line 41 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
      end do
    end do

    call grid_free(grid)
#line 45 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertFalse(allocated(grid%cell_volume), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 45) )
  if (anyExceptions()) return
#line 46 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  end subroutine test_rect_all_active

  !@test
  subroutine test_rect_material_id()
    type(spatial_grid) :: grid
    integer, parameter :: nx = 3, ny = 3
    integer :: mat2d(nx, ny)
    integer :: ij

    mat2d = reshape([1,1,1, 2,2,2, 3,3,3], [nx, ny])
    call grid_init_rect(grid, nx, ny, 1.0_dp, 1.0_dp, mat2d)

    ij = (1-1)*nx + 1
#line 59 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1, grid%material_id(ij), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 59) )
  if (anyExceptions()) return
#line 60 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
    ij = (2-1)*nx + 1
#line 61 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(2, grid%material_id(ij), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 61) )
  if (anyExceptions()) return
#line 62 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
    ij = (3-1)*nx + 3
#line 63 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(3, grid%material_id(ij), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 63) )
  if (anyExceptions()) return
#line 64 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"

    call grid_free(grid)
  end subroutine test_rect_material_id

  ! ==================================================================
  ! Circle tests
  ! ==================================================================

  !@test
  subroutine test_circle_total_area()
    type(spatial_grid) :: grid
    integer, parameter :: nx = 100, ny = 100
    real(kind=dp), parameter :: dx = 0.1_dp, dy = 0.1_dp
    real(kind=dp), parameter :: cx = 5.0_dp, cy = 5.0_dp, radius = 4.0_dp
    real(kind=dp) :: total_area, expected_area, rel_err
    logical :: ok
    integer :: mat2d(nx, ny), ij

    mat2d = 1
    call grid_init_circle(grid, nx, ny, dx, dy, cx, cy, radius, mat2d)

    total_area = 0.0_dp
    do ij = 1, nx * ny
      total_area = total_area + grid%cell_volume(ij)
    end do
    total_area = total_area * dx * dy

    expected_area = pi_dp * radius**2
    rel_err = abs(total_area - expected_area) / expected_area
    ok = (rel_err < 0.01_dp)
#line 94 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertTrue(ok, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 94) )
  if (anyExceptions()) return
#line 95 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"

    call grid_free(grid)
  end subroutine test_circle_total_area

  !@test
  subroutine test_circle_outside_zero()
    type(spatial_grid) :: grid
    integer, parameter :: nx = 20, ny = 20
    real(kind=dp), parameter :: dx = 1.0_dp, dy = 1.0_dp
    real(kind=dp), parameter :: cx = 10.0_dp, cy = 10.0_dp, radius = 3.0_dp
    integer :: mat2d(nx, ny)
    integer :: ij

    mat2d = 1
    call grid_init_circle(grid, nx, ny, dx, dy, cx, cy, radius, mat2d)

    ij = (1-1)*nx + 1
#line 112 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(0.0_dp, grid%cell_volume(ij), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 112) )
  if (anyExceptions()) return
#line 113 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 113 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(0, grid%material_id(ij), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 113) )
  if (anyExceptions()) return
#line 114 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 114 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(0.0_dp, grid%face_fraction_x(ij, 1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 114) )
  if (anyExceptions()) return
#line 115 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 115 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(0.0_dp, grid%face_fraction_y(ij, 1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 115) )
  if (anyExceptions()) return
#line 116 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"

    call grid_free(grid)
  end subroutine test_circle_outside_zero

  !@test
  subroutine test_circle_center_inside()
    type(spatial_grid) :: grid
    integer, parameter :: nx = 20, ny = 20
    real(kind=dp), parameter :: dx = 1.0_dp, dy = 1.0_dp
    real(kind=dp), parameter :: cx = 10.0_dp, cy = 10.0_dp, radius = 8.0_dp
    integer :: mat2d(nx, ny)
    integer :: ij

    mat2d = 1
    call grid_init_circle(grid, nx, ny, dx, dy, cx, cy, radius, mat2d)

    ! Cell ix=10,iy=10 spans x=[9,10], y=[9,10]; all corners within radius 8
    ij = (10-1)*nx + 10
#line 134 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1.0_dp, grid%cell_volume(ij), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 134) )
  if (anyExceptions()) return
#line 135 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 135 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1, grid%material_id(ij), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 135) )
  if (anyExceptions()) return
#line 136 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"

    call grid_free(grid)
  end subroutine test_circle_center_inside

  ! ==================================================================
  ! Hexagon tests
  ! ==================================================================

  !@test
  subroutine test_hexagon_total_area()
    type(spatial_grid) :: grid
    integer, parameter :: nx = 80, ny = 80
    real(kind=dp), parameter :: dx = 0.1_dp, dy = 0.1_dp
    real(kind=dp), parameter :: cx = 4.0_dp, cy = 4.0_dp
    real(kind=dp), parameter :: side = 3.0_dp
    real(kind=dp) :: total_area, expected_area, rel_err
    logical :: ok
    integer :: mat2d(nx, ny), ij

    mat2d = 1
    call grid_init_hexagon(grid, nx, ny, dx, dy, cx, cy, side, mat2d)

    total_area = 0.0_dp
    do ij = 1, nx * ny
      total_area = total_area + grid%cell_volume(ij)
    end do
    total_area = total_area * dx * dy

    expected_area = 1.5_dp * sqrt(3.0_dp) * side**2
    rel_err = abs(total_area - expected_area) / expected_area
    ok = (rel_err < 0.02_dp)
#line 167 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertTrue(ok, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 167) )
  if (anyExceptions()) return
#line 168 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"

    call grid_free(grid)
  end subroutine test_hexagon_total_area

  ! ==================================================================
  ! Polygon tests
  ! ==================================================================

  !@test
  subroutine test_polygon_triangle_area()
    type(spatial_grid) :: grid
    integer, parameter :: nx = 100, ny = 100
    real(kind=dp), parameter :: dx = 0.1_dp, dy = 0.1_dp
    real(kind=dp) :: vx(3), vy(3)
    real(kind=dp) :: total_area, expected_area, rel_err
    logical :: ok
    integer :: mat2d(nx, ny), ij

    vx = [2.0_dp, 8.0_dp, 2.0_dp]
    vy = [2.0_dp, 2.0_dp, 8.0_dp]
    expected_area = 0.5_dp * 6.0_dp * 6.0_dp

    mat2d = 1
    call grid_init_polygon(grid, nx, ny, dx, dy, 3, vx, vy, mat2d)

    total_area = 0.0_dp
    do ij = 1, nx * ny
      total_area = total_area + grid%cell_volume(ij)
    end do
    total_area = total_area * dx * dy

    rel_err = abs(total_area - expected_area) / expected_area
    ok = (rel_err < 0.02_dp)
#line 201 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertTrue(ok, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 201) )
  if (anyExceptions()) return
#line 202 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"

    call grid_free(grid)
  end subroutine test_polygon_triangle_area

  ! ==================================================================
  ! Coordinate consistency tests
  ! ==================================================================

  !@test
  subroutine test_grid_coords_consistency()
    type(spatial_grid) :: grid
    integer, parameter :: nx = 5, ny = 7
    real(kind=dp), parameter :: dx = 0.5_dp, dy = 1.0_dp
    integer :: mat2d(nx, ny)
    integer :: ix, iy, ij

    mat2d = 1
    call grid_init_rect(grid, nx, ny, dx, dy, mat2d)
    call grid_build_coords(grid)

#line 222 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertTrue(allocated(grid%coords), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 222) )
  if (anyExceptions()) return
#line 223 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 223 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(2, size(grid%coords, 1), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 223) )
  if (anyExceptions()) return
#line 224 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 224 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(nx * ny, size(grid%coords, 2), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 224) )
  if (anyExceptions()) return
#line 225 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"

    do iy = 1, ny
      do ix = 1, nx
        ij = (iy - 1) * nx + ix
#line 229 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual((ix - 0.5_dp) * dx, grid%coords(1, ij), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 229) )
  if (anyExceptions()) return
#line 230 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 230 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual((iy - 0.5_dp) * dy, grid%coords(2, ij), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 230) )
  if (anyExceptions()) return
#line 231 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
      end do
    end do

    call grid_free(grid)
  end subroutine test_grid_coords_consistency

  !@test
  subroutine test_grid_coord_arrays()
    type(spatial_grid) :: grid
    integer, parameter :: nx = 4, ny = 3
    real(kind=dp), parameter :: dx = 2.0_dp, dy = 3.0_dp
    integer :: mat2d(nx, ny)
    integer :: i

    mat2d = 1
    call grid_init_rect(grid, nx, ny, dx, dy, mat2d)

#line 248 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertTrue(allocated(grid%x), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 248) )
  if (anyExceptions()) return
#line 249 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 249 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertTrue(allocated(grid%z), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 249) )
  if (anyExceptions()) return
#line 250 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 250 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(nx, size(grid%x), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 250) )
  if (anyExceptions()) return
#line 251 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 251 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(ny, size(grid%z), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 251) )
  if (anyExceptions()) return
#line 252 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"

    do i = 1, nx
#line 254 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual((i - 0.5_dp) * dx, grid%x(i), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 254) )
  if (anyExceptions()) return
#line 255 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
    end do
    do i = 1, ny
#line 257 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual((i - 0.5_dp) * dy, grid%z(i), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 257) )
  if (anyExceptions()) return
#line 258 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
    end do

    call grid_free(grid)
  end subroutine test_grid_coord_arrays

  ! ==================================================================
  ! Ghost map tests
  ! ==================================================================

  !@test
  subroutine test_rect_ghost_map()
    type(spatial_grid) :: grid
    integer, parameter :: nx = 3, ny = 3
    integer :: mat2d(nx, ny)
    integer :: ij

    mat2d = 1
    call grid_init_rect(grid, nx, ny, 1.0_dp, 1.0_dp, mat2d)

    ! Interior cell (2,2): neighbors should be immediate
    ij = (2-1)*nx + 2
    ! N (+y): (2,3) -> ij = (3-1)*3 + 2 = 8
#line 280 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(8, grid%ghost_map(ij, 1), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 280) )
  if (anyExceptions()) return
#line 281 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
    ! S (-y): (2,1) -> ij = (1-1)*3 + 2 = 2
#line 282 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(2, grid%ghost_map(ij, 2), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 282) )
  if (anyExceptions()) return
#line 283 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
    ! W (-x): (1,2) -> ij = (2-1)*3 + 1 = 4
#line 284 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(4, grid%ghost_map(ij, 3), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 284) )
  if (anyExceptions()) return
#line 285 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
    ! E (+x): (3,2) -> ij = (2-1)*3 + 3 = 6
#line 286 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(6, grid%ghost_map(ij, 4), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 286) )
  if (anyExceptions()) return
#line 287 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"

    ! Boundary cell (1,1): S and W should map to self
    ij = (1-1)*nx + 1
#line 290 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(ij, grid%ghost_map(ij, 2), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 290) )
  if (anyExceptions()) return
#line 291 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 291 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(ij, grid%ghost_map(ij, 3), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 291) )
  if (anyExceptions()) return
#line 292 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"

    call grid_free(grid)
  end subroutine test_rect_ghost_map

  ! ==================================================================
  ! Grid free tests
  ! ==================================================================

  !@test
  subroutine test_grid_free_resets()
    type(spatial_grid) :: grid
    integer :: mat2d(2, 2)

    mat2d = 1
    call grid_init_rect(grid, 2, 2, 1.0_dp, 1.0_dp, mat2d)
    call grid_free(grid)

#line 309 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(0, grid%ndim, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 309) )
  if (anyExceptions()) return
#line 310 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 310 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1, grid%nx, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 310) )
  if (anyExceptions()) return
#line 311 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 311 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1, grid%ny, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 311) )
  if (anyExceptions()) return
#line 312 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 312 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(0.0_dp, grid%dx, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 312) )
  if (anyExceptions()) return
#line 313 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 313 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(0.0_dp, grid%dy, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 313) )
  if (anyExceptions()) return
#line 314 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 314 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertFalse(allocated(grid%x), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 314) )
  if (anyExceptions()) return
#line 315 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 315 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertFalse(allocated(grid%z), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 315) )
  if (anyExceptions()) return
#line 316 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 316 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertFalse(allocated(grid%coords), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 316) )
  if (anyExceptions()) return
#line 317 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 317 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertFalse(allocated(grid%cell_volume), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 317) )
  if (anyExceptions()) return
#line 318 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 318 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertFalse(allocated(grid%ghost_map), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 318) )
  if (anyExceptions()) return
#line 319 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  end subroutine test_grid_free_resets

  ! ==================================================================
  ! init_wire_from_config integration tests
  ! ==================================================================

  !@test
  subroutine test_wire_rect_single_material()
    ! Rectangle wire with single region: all cells active, material_id=1.
    type(simulation_config) :: cfg
    integer :: ij, ngrid
    real(kind=dp) :: xval, yval

    cfg%confinement = 2
    cfg%wire_nx = 5
    cfg%wire_ny = 5
    cfg%wire_dx = 2.0_dp
    cfg%wire_dy = 2.0_dp
    cfg%FDorder = 2

    ! Set grid dimensions (as init_grid_from_config would do)
    cfg%grid%ndim = 2
    cfg%grid%nx   = 5
    cfg%grid%ny   = 5
    cfg%grid%dx   = 2.0_dp
    cfg%grid%dy   = 2.0_dp

    cfg%wire_geom%shape = 'rectangle'

    ! Single region covering entire wire
    cfg%numRegions = 1
    allocate(cfg%regions(1))
    cfg%regions(1)%material = 'GaAs'
    cfg%regions(1)%inner = 0.0_dp
    cfg%regions(1)%outer = 100.0_dp  ! larger than any distance

    call init_wire_from_config(cfg)

    ngrid = 5 * 5

    ! All cells should be active (rectangle, fully inside)
#line 360 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertTrue(allocated(cfg%grid%cell_volume), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 360) )
  if (anyExceptions()) return
#line 361 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 361 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertTrue(allocated(cfg%grid%material_id), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 361) )
  if (anyExceptions()) return
#line 362 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 362 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertTrue(allocated(cfg%grid%ghost_map), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 362) )
  if (anyExceptions()) return
#line 363 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 363 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertTrue(allocated(cfg%grid%coords), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 363) )
  if (anyExceptions()) return
#line 364 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"

    do ij = 1, ngrid
#line 366 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1.0_dp, cfg%grid%cell_volume(ij), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 366) )
  if (anyExceptions()) return
#line 367 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 367 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1.0_dp, cfg%grid%face_fraction_x(ij, 1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 367) )
  if (anyExceptions()) return
#line 368 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 368 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1.0_dp, cfg%grid%face_fraction_x(ij, 2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 368) )
  if (anyExceptions()) return
#line 369 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 369 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1.0_dp, cfg%grid%face_fraction_y(ij, 1), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 369) )
  if (anyExceptions()) return
#line 370 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 370 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1.0_dp, cfg%grid%face_fraction_y(ij, 2), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 370) )
  if (anyExceptions()) return
#line 371 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 371 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1, cfg%grid%material_id(ij), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 371) )
  if (anyExceptions()) return
#line 372 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
    end do

    ! Coords should be cell-centered: x(1) = 0.5*dx = 1.0, y(1) = 0.5*dy = 1.0
    xval = cfg%grid%coords(1, 1)
    yval = cfg%grid%coords(2, 1)
#line 377 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1.0_dp, xval, tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 377) )
  if (anyExceptions()) return
#line 378 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
#line 378 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1.0_dp, yval, tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 378) )
  if (anyExceptions()) return
#line 379 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"

    call grid_free(cfg%grid)
    deallocate(cfg%regions)
  end subroutine test_wire_rect_single_material

  !@test
  subroutine test_wire_rect_two_regions()
    ! Rectangle wire with core (region 1, r<=4) and shell (region 2, r>4).
    ! Grid: 10x10, dx=dy=1. Grid center at (5,5).
    type(simulation_config) :: cfg
    integer :: ix, iy, ij
    real(kind=dp) :: cx, cy, dist

    cfg%confinement = 2
    cfg%wire_nx = 10
    cfg%wire_ny = 10
    cfg%wire_dx = 1.0_dp
    cfg%wire_dy = 1.0_dp
    cfg%FDorder = 2

    cfg%grid%ndim = 2
    cfg%grid%nx   = 10
    cfg%grid%ny   = 10
    cfg%grid%dx   = 1.0_dp
    cfg%grid%dy   = 1.0_dp

    cfg%wire_geom%shape = 'rectangle'

    ! Two concentric regions
    cfg%numRegions = 2
    allocate(cfg%regions(2))
    cfg%regions(1)%material = 'InAs'
    cfg%regions(1)%inner = 0.0_dp
    cfg%regions(1)%outer = 4.0_dp    ! core: r <= 4
    cfg%regions(2)%material = 'GaSb'
    cfg%regions(2)%inner = 4.0_dp
    cfg%regions(2)%outer = 100.0_dp  ! shell: r > 4

    call init_wire_from_config(cfg)

    ! Grid center is at (5, 5)
    cx = 5.0_dp
    cy = 5.0_dp

    ! Verify material_id is based on 2D distance from center
    do iy = 1, 10
      do ix = 1, 10
        ij = (iy - 1) * 10 + ix
        dist = sqrt(((ix - 0.5_dp) * 1.0_dp - cx)**2 + &
                    ((iy - 0.5_dp) * 1.0_dp - cy)**2)
        if (dist <= 4.0_dp) then
#line 430 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1, cfg%grid%material_id(ij), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 430) )
  if (anyExceptions()) return
#line 431 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
        else
#line 432 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(2, cfg%grid%material_id(ij), &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 432) )
  if (anyExceptions()) return
#line 433 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
        end if
      end do
    end do

    ! All cells should have volume=1 (rectangle)
    do ij = 1, 100
#line 439 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
  call assertEqual(1.0_dp, cfg%grid%cell_volume(ij), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_geometry.pf', &
 & 439) )
  if (anyExceptions()) return
#line 440 "/data/8bandkp-fdm/tests/unit/test_geometry.pf"
    end do

    call grid_free(cfg%grid)
    deallocate(cfg%regions)
  end subroutine test_wire_rect_two_regions

end module test_geometry

module Wraptest_geometry
   use FUnit
   use test_geometry
   implicit none
   private

contains


end module Wraptest_geometry

function test_geometry_suite() result(suite)
   use FUnit
   use test_geometry
   use Wraptest_geometry
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   suite = TestSuite('test_geometry_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_rect_all_active', &
      test_rect_all_active))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_rect_material_id', &
      test_rect_material_id))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_circle_total_area', &
      test_circle_total_area))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_circle_outside_zero', &
      test_circle_outside_zero))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_circle_center_inside', &
      test_circle_center_inside))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_hexagon_total_area', &
      test_hexagon_total_area))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_polygon_triangle_area', &
      test_polygon_triangle_area))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_grid_coords_consistency', &
      test_grid_coords_consistency))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_grid_coord_arrays', &
      test_grid_coord_arrays))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_rect_ghost_map', &
      test_rect_ghost_map))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_grid_free_resets', &
      test_grid_free_resets))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_wire_rect_single_material', &
      test_wire_rect_single_material))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_wire_rect_two_regions', &
      test_wire_rect_two_regions))
   call suite%addTest(t)


end function test_geometry_suite

