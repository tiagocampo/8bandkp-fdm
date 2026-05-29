module test_defs
  use funit
  use definitions
  implicit none

contains

  !@test
  subroutine test_kronij_same()
#line 10 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1, kronij(3, 3), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 10) )
  if (anyExceptions()) return
#line 11 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 11 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1, kronij(1, 1), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 11) )
  if (anyExceptions()) return
#line 12 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 12 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1, kronij(0, 0), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 12) )
  if (anyExceptions()) return
#line 13 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_kronij_same

  !@test
  subroutine test_kronij_different()
#line 17 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(0, kronij(1, 2), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 17) )
  if (anyExceptions()) return
#line 18 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 18 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(0, kronij(3, 1), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 18) )
  if (anyExceptions()) return
#line 19 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 19 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(0, kronij(-1, 1), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 19) )
  if (anyExceptions()) return
#line 20 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_kronij_different

  !@test
  subroutine test_pi_constant()
    real(kind=dp) :: expected_pi
    expected_pi = acos(-1.0_dp)
#line 26 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(expected_pi, pi_dp, tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 26) )
  if (anyExceptions()) return
#line 27 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_pi_constant

  !@test
  subroutine test_imaginary_unit()
    complex(kind=dp) :: result
    result = IU * IU
#line 33 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(-1.0_dp, real(result, kind=dp), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 33) )
  if (anyExceptions()) return
#line 34 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 34 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(0.0_dp, aimag(result), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 34) )
  if (anyExceptions()) return
#line 35 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_imaginary_unit

  !@test
  subroutine test_unit_complex()
#line 39 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(0.0_dp, real(ZERO, kind=dp), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 39) )
  if (anyExceptions()) return
#line 40 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 40 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(0.0_dp, aimag(ZERO), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 40) )
  if (anyExceptions()) return
#line 41 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 41 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1.0_dp, real(UM, kind=dp), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 41) )
  if (anyExceptions()) return
#line 42 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 42 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(0.0_dp, aimag(UM), tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 42) )
  if (anyExceptions()) return
#line 43 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_unit_complex

  !@test
  subroutine test_sqrt_constants()
#line 47 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(sqrt(3.0_dp), SQR3, tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 47) )
  if (anyExceptions()) return
#line 48 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 48 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(sqrt(2.0_dp), SQR2, tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 48) )
  if (anyExceptions()) return
#line 49 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 49 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1.0_dp/SQR3, RQS3, tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 49) )
  if (anyExceptions()) return
#line 50 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 50 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1.0_dp/SQR2, RQS2, tolerance=1.0e-14_dp, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 50) )
  if (anyExceptions()) return
#line 51 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_sqrt_constants

  !@test
  subroutine test_conf_direction_bulk()
#line 55 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual('n', conf_direction('bulk'), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 55) )
  if (anyExceptions()) return
#line 56 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_conf_direction_bulk

  !@test
  subroutine test_conf_direction_qw()
#line 60 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual('z', conf_direction('qw'), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 60) )
  if (anyExceptions()) return
#line 61 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_conf_direction_qw

  !@test
  subroutine test_conf_direction_wire()
#line 65 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual('z', conf_direction('wire'), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 65) )
  if (anyExceptions()) return
#line 66 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_conf_direction_wire

  !@test
  subroutine test_conf_direction_landau()
#line 70 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual('x', conf_direction('landau'), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 70) )
  if (anyExceptions()) return
#line 71 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_conf_direction_landau

  ! ------------------------------------------------------------------
  ! num_layers semantics: num_layers = material layers, not regions.
  ! Wire mode stores its region count in cfg%wire%num_regions.
  ! ------------------------------------------------------------------

  !@test
  subroutine test_wire_num_layers_is_one_single_region()
    ! Wire with 1 region: num_layers must be 1, not the region count.
    type(simulation_config) :: cfg

    cfg%confinement = 'wire'
    cfg%wire%num_regions = 1
    cfg%num_layers = 1

#line 87 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1, cfg%num_layers, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 87) )
  if (anyExceptions()) return
#line 88 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 88 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1, cfg%wire%num_regions, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 88) )
  if (anyExceptions()) return
#line 89 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_wire_num_layers_is_one_single_region

  !@test
  subroutine test_wire_num_layers_is_one_multi_region()
    ! Wire with 3 regions: num_layers must be 1, regions in wire%num_regions.
    type(simulation_config) :: cfg

    cfg%confinement = 'wire'
    cfg%wire%num_regions = 3
    cfg%num_layers = 1

#line 100 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1, cfg%num_layers, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 100) )
  if (anyExceptions()) return
#line 101 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 101 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(3, cfg%wire%num_regions, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 101) )
  if (anyExceptions()) return
#line 102 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_wire_num_layers_is_one_multi_region

  !@test
  subroutine test_wire_material_array_sized_by_regions()
    ! Wire material arrays are sized by region count, not num_layers.
    ! After parser, size(material_names) == wire%num_regions, num_layers == 1.
    type(simulation_config) :: cfg

    cfg%confinement = 'wire'
    cfg%wire%num_regions = 2
    cfg%num_layers = 1

    allocate(cfg%material_names(cfg%wire%num_regions))
    cfg%material_names(1) = 'InAs'
    cfg%material_names(2) = 'GaAs'

#line 118 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1, cfg%num_layers, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 118) )
  if (anyExceptions()) return
#line 119 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 119 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(2, cfg%wire%num_regions, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 119) )
  if (anyExceptions()) return
#line 120 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 120 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(2, size(cfg%material_names), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 120) )
  if (anyExceptions()) return
#line 121 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_wire_material_array_sized_by_regions

  !@test
  subroutine test_qw_num_layers_equals_nmat()
    ! QW with 3 materials: num_layers = 3 (material layer count).
    type(simulation_config) :: cfg

    cfg%confinement = 'qw'
    cfg%num_layers = 3

#line 131 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(3, cfg%num_layers, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 131) )
  if (anyExceptions()) return
#line 132 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_qw_num_layers_equals_nmat

  !@test
  subroutine test_bulk_num_layers_is_one()
    ! Bulk: single material, num_layers = 1.
    type(simulation_config) :: cfg

    cfg%confinement = 'bulk'
    cfg%num_layers = 1

#line 142 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1, cfg%num_layers, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 142) )
  if (anyExceptions()) return
#line 143 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_bulk_num_layers_is_one

  !@test
  subroutine test_landau_num_layers_is_one()
    ! Landau: single material, num_layers = 1.
    type(simulation_config) :: cfg

    cfg%confinement = 'landau'
    cfg%num_layers = 1

#line 153 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual(1, cfg%num_layers, &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_landau_num_layers_is_one

  ! ------------------------------------------------------------------
  ! conf_direction with unrecognized confinement string returns 'n'.
  ! ------------------------------------------------------------------

  !@test
  subroutine test_conf_direction_unknown()
#line 162 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual('n', conf_direction('unknown'), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 162) )
  if (anyExceptions()) return
#line 163 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 163 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual('n', conf_direction(''), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 163) )
  if (anyExceptions()) return
#line 164 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
#line 164 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertEqual('n', conf_direction('QW'), &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 164) )
  if (anyExceptions()) return
#line 165 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_conf_direction_unknown

  ! ------------------------------------------------------------------
  ! validate_semantic: wire config with num_layers=1 and multiple regions
  ! should pass validation (material arrays sized by region count, not
  ! num_layers).  This tests the acceptance criterion for Issue 03.
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_wire_multi_region()
    type(simulation_config) :: cfg

    ! Confinement and geometry
    cfg%confinement = 'wire'
    cfg%FDorder     = 2

    ! Wire geometry: 3 regions, grid large enough for FDorder=2
    cfg%wire%nx           = 5
    cfg%wire%ny           = 5
    cfg%wire%num_regions  = 3

    ! Grid must be consistent with wire dimensions
    cfg%grid%ndim = 2
    cfg%grid%nx   = 5
    cfg%grid%ny   = 5

    ! Bands
    cfg%bands%num_cb = 2
    cfg%bands%num_vb = 6
    cfg%evnum        = 8

    ! Material layers: num_layers=1, but arrays sized by region count (3)
    cfg%num_layers = 1
    allocate(cfg%material_names(3))
    cfg%material_names(1) = 'InAs'
    cfg%material_names(2) = 'GaAs'
    cfg%material_names(3) = 'InAs'
    allocate(cfg%params(3))
    cfg%params(:)%Eg = 0.0_dp

    ! Call validate — should not error stop
    call cfg%validate()

    ! If we reach here, validation passed
#line 209 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertTrue(.true., &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 209) )
  if (anyExceptions()) return
#line 210 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_validate_wire_multi_region

  ! ------------------------------------------------------------------
  ! V1: bulk accepts evnum=8 (8x8 Hamiltonian limit)
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_bulk_evnum_ok()
    type(simulation_config) :: cfg

    cfg%confinement = 'bulk'
    cfg%FDorder     = 2
    cfg%fd_step     = 1
    cfg%grid%ndim   = 0
    cfg%grid%nx     = 1
    cfg%grid%ny     = 1

    cfg%bands%num_cb = 2
    cfg%bands%num_vb = 6
    cfg%evnum        = 8

    cfg%num_layers = 1
    cfg%wave_vector%mode = 'k0'
    allocate(cfg%material_names(1))
    cfg%material_names(1) = 'GaAs'
    allocate(cfg%params(1))
    cfg%params(:)%Eg = 1.519_dp

    call cfg%validate()
#line 239 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertTrue(.true., &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 239) )
  if (anyExceptions()) return
#line 240 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_validate_bulk_evnum_ok

  ! ------------------------------------------------------------------
  ! V2: QW accepts valid num_cb (num_cb < NUM_CB_STATES * npoints)
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_qw_num_cb_ok()
    type(simulation_config) :: cfg

    cfg%confinement = 'qw'
    cfg%FDorder     = 2
    cfg%fd_step     = 50
    cfg%grid%ndim   = 1
    cfg%grid%nx     = 1
    cfg%grid%ny     = 50

    cfg%bands%num_cb = 4
    cfg%bands%num_vb = 12
    cfg%evnum        = 16

    cfg%num_layers = 1
    cfg%wave_vector%mode = 'k0'
    allocate(cfg%material_names(1))
    cfg%material_names(1) = 'GaAs'
    allocate(cfg%params(1))
    cfg%params(:)%Eg = 1.519_dp

    call cfg%validate()
#line 269 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertTrue(.true., &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 269) )
  if (anyExceptions()) return
#line 270 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_validate_qw_num_cb_ok

  ! ------------------------------------------------------------------
  ! V3: QW accepts valid num_vb (num_vb < NUM_VB_STATES * npoints)
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_qw_num_vb_ok()
    type(simulation_config) :: cfg

    cfg%confinement = 'qw'
    cfg%FDorder     = 2
    cfg%fd_step     = 50
    cfg%grid%ndim   = 1
    cfg%grid%nx     = 1
    cfg%grid%ny     = 50

    cfg%bands%num_cb = 2
    cfg%bands%num_vb = 12
    cfg%evnum        = 14

    cfg%num_layers = 1
    cfg%wave_vector%mode = 'k0'
    allocate(cfg%material_names(1))
    cfg%material_names(1) = 'GaAs'
    allocate(cfg%params(1))
    cfg%params(:)%Eg = 1.519_dp

    call cfg%validate()
#line 299 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertTrue(.true., &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 299) )
  if (anyExceptions()) return
#line 300 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_validate_qw_num_vb_ok

  ! ------------------------------------------------------------------
  ! V4: wave_vector mode='k0' is recognized
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_wave_vector_k0()
    type(simulation_config) :: cfg

    cfg%confinement = 'bulk'
    cfg%FDorder     = 2
    cfg%fd_step     = 1
    cfg%grid%ndim   = 0
    cfg%grid%nx     = 1
    cfg%grid%ny     = 1

    cfg%bands%num_cb = 2
    cfg%bands%num_vb = 6
    cfg%evnum        = 8

    cfg%num_layers = 1
    cfg%wave_vector%mode = 'k0'
    allocate(cfg%material_names(1))
    cfg%material_names(1) = 'GaAs'
    allocate(cfg%params(1))
    cfg%params(:)%Eg = 1.519_dp

    call cfg%validate()
#line 329 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertTrue(.true., &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 329) )
  if (anyExceptions()) return
#line 330 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_validate_wave_vector_k0

  ! ------------------------------------------------------------------
  ! V5: B-sweep with positive step accepted
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_bsweep_positive()
    type(simulation_config) :: cfg

    cfg%confinement = 'bulk'
    cfg%FDorder     = 2
    cfg%fd_step     = 1
    cfg%grid%ndim   = 0
    cfg%grid%nx     = 1
    cfg%grid%ny     = 1

    cfg%bands%num_cb = 2
    cfg%bands%num_vb = 6
    cfg%evnum        = 8

    cfg%num_layers = 1
    cfg%wave_vector%mode = 'k0'
    allocate(cfg%material_names(1))
    cfg%material_names(1) = 'GaAs'
    allocate(cfg%params(1))
    cfg%params(:)%Eg = 1.519_dp

    ! Positive B-sweep step
    cfg%bdg%B_sweep = [0.0_dp, 1.0_dp, 0.1_dp]

    call cfg%validate()
#line 362 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertTrue(.true., &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 362) )
  if (anyExceptions()) return
#line 363 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_validate_bsweep_positive

  ! ------------------------------------------------------------------
  ! V6: z_min < z_max for all material layers accepted
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_zmin_lt_zmax()
    type(simulation_config) :: cfg

    cfg%confinement = 'qw'
    cfg%FDorder     = 2
    cfg%fd_step     = 50
    cfg%grid%ndim   = 1
    cfg%grid%nx     = 1
    cfg%grid%ny     = 50

    cfg%bands%num_cb = 2
    cfg%bands%num_vb = 6
    cfg%evnum        = 8

    cfg%num_layers = 1
    cfg%wave_vector%mode = 'k0'
    allocate(cfg%material_names(1))
    cfg%material_names(1) = 'GaAs'
    allocate(cfg%params(1))
    cfg%params(:)%Eg = 1.519_dp

    ! Valid material layer boundaries
    allocate(cfg%z_min(1))
    allocate(cfg%z_max(1))
    cfg%z_min(1) = 0.0_dp
    cfg%z_max(1) = 10.0_dp

    call cfg%validate()
#line 398 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertTrue(.true., &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 398) )
  if (anyExceptions()) return
#line 399 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_validate_zmin_lt_zmax

  ! ------------------------------------------------------------------
  ! V7: SC with QW confinement accepted
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_sc_qw_ok()
    type(simulation_config) :: cfg

    cfg%confinement = 'qw'
    cfg%FDorder     = 2
    cfg%fd_step     = 50
    cfg%grid%ndim   = 1
    cfg%grid%nx     = 1
    cfg%grid%ny     = 50

    cfg%bands%num_cb = 2
    cfg%bands%num_vb = 6
    cfg%evnum        = 8

    cfg%num_layers = 1
    cfg%wave_vector%mode = 'k0'
    allocate(cfg%material_names(1))
    cfg%material_names(1) = 'GaAs'
    allocate(cfg%params(1))
    cfg%params(:)%Eg = 1.519_dp

    ! SC enabled with QW confinement — valid
    cfg%sc%enabled = 1

    call cfg%validate()
#line 431 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertTrue(.true., &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 431) )
  if (anyExceptions()) return
#line 432 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_validate_sc_qw_ok

  ! ------------------------------------------------------------------
  ! V8: valid electric field with z(1) /= 0 accepted
  ! ------------------------------------------------------------------

  !@test
  subroutine test_validate_efield_z_nonzero()
    type(simulation_config) :: cfg

    cfg%confinement = 'qw'
    cfg%FDorder     = 2
    cfg%fd_step     = 50
    cfg%grid%ndim   = 1
    cfg%grid%nx     = 1
    cfg%grid%ny     = 50

    cfg%bands%num_cb = 2
    cfg%bands%num_vb = 6
    cfg%evnum        = 8

    cfg%num_layers = 1
    cfg%wave_vector%mode = 'k0'
    allocate(cfg%material_names(1))
    cfg%material_names(1) = 'GaAs'
    allocate(cfg%params(1))
    cfg%params(:)%Eg = 1.519_dp

    ! Electric field enabled with z(1) /= 0
    cfg%external_field%enabled = .true.
    cfg%external_field%type = 'EF'
    allocate(cfg%z(2))
    cfg%z(1) = 0.5_dp
    cfg%z(2) = 1.0_dp

    call cfg%validate()
#line 468 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  call assertTrue(.true., &
 & location=SourceLocation( &
 & 'test_defs.pf', &
 & 468) )
  if (anyExceptions()) return
#line 469 "/data/8bandkp-fdm/tests/unit/test_defs.pf"
  end subroutine test_validate_efield_z_nonzero

end module test_defs

module Wraptest_defs
   use FUnit
   use test_defs
   implicit none
   private

contains


end module Wraptest_defs

function test_defs_suite() result(suite)
   use FUnit
   use test_defs
   use Wraptest_defs
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   suite = TestSuite('test_defs_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_kronij_same', &
      test_kronij_same))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_kronij_different', &
      test_kronij_different))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_pi_constant', &
      test_pi_constant))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_imaginary_unit', &
      test_imaginary_unit))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_unit_complex', &
      test_unit_complex))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_sqrt_constants', &
      test_sqrt_constants))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_conf_direction_bulk', &
      test_conf_direction_bulk))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_conf_direction_qw', &
      test_conf_direction_qw))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_conf_direction_wire', &
      test_conf_direction_wire))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_conf_direction_landau', &
      test_conf_direction_landau))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_wire_num_layers_is_one_single_region', &
      test_wire_num_layers_is_one_single_region))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_wire_num_layers_is_one_multi_region', &
      test_wire_num_layers_is_one_multi_region))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_wire_material_array_sized_by_regions', &
      test_wire_material_array_sized_by_regions))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_qw_num_layers_equals_nmat', &
      test_qw_num_layers_equals_nmat))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_bulk_num_layers_is_one', &
      test_bulk_num_layers_is_one))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_landau_num_layers_is_one', &
      test_landau_num_layers_is_one))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_conf_direction_unknown', &
      test_conf_direction_unknown))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_wire_multi_region', &
      test_validate_wire_multi_region))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_bulk_evnum_ok', &
      test_validate_bulk_evnum_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_qw_num_cb_ok', &
      test_validate_qw_num_cb_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_qw_num_vb_ok', &
      test_validate_qw_num_vb_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_wave_vector_k0', &
      test_validate_wave_vector_k0))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_bsweep_positive', &
      test_validate_bsweep_positive))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_zmin_lt_zmax', &
      test_validate_zmin_lt_zmax))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_sc_qw_ok', &
      test_validate_sc_qw_ok))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_validate_efield_z_nonzero', &
      test_validate_efield_z_nonzero))
   call suite%addTest(t)


end function test_defs_suite

