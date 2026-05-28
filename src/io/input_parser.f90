module input_parser

  use definitions
  use parameters
  use geometry
  use outputFunctions
  use tomlf, only : get_value, set_value, toml_path, toml_parse, toml_load, &
    & toml_loads, toml_context, toml_parser_config, toml_level, &
    & toml_error, toml_stat, toml_serializer, toml_serialize, toml_dump, &
    & toml_dumps, toml_terminal, toml_table, toml_array, toml_keyval, &
    & toml_key, toml_value, is_array_of_tables, new_table, add_table, &
    & add_array, add_keyval, tomlf_len => len

  implicit none

  private
  public :: read_config

contains

  ! ==================================================================
  ! Top-level parser: reads input.toml and populates simulation_config.
  ! For this slice, only bulk mode is fully functional.
  ! QW/wire/landau modes will throw "not yet implemented" errors.
  ! ==================================================================
  subroutine read_config(cfg)
    type(simulation_config), intent(out) :: cfg

    type(toml_table), allocatable :: table
    type(toml_table), pointer     :: child => null()
    type(toml_array), pointer     :: arr => null()
    type(toml_error), allocatable :: parse_error
    character(len=:), allocatable :: str_val
    integer :: stat, i

    ! Load TOML file
    call toml_load(table, 'input.toml', error=parse_error)
    if (allocated(parse_error)) then
      print *, 'Error: Failed to parse input.toml'
      print *, '  ', trim(parse_error%message)
      stop 1
    end if

    ! ---- Required top-level fields ----
    call require_string(table, 'confinement', str_val, 'top-level')
    cfg%confinement = trim(str_val)

    call get_value(table, 'FDorder', cfg%FDorder, stat=stat)
    if (stat /= 0) then
      cfg%FDorder = 2  ! default
    end if

    call get_value(table, 'fd_step', cfg%fd_step, stat=stat)
    if (stat /= 0) then
      cfg%fd_step = 1  ! default
    end if

    ! ---- [wave_vector] ----
    call get_value(table, 'wave_vector', child, stat=stat)
    if (stat == 0 .and. associated(child)) then
      call require_string(child, 'mode', str_val, 'wave_vector')
      cfg%wave_vector%mode = trim(str_val)
      call get_value(child, 'max', cfg%wave_vector%max, stat=stat)
      if (stat /= 0) cfg%wave_vector%max = 0.0_dp
      call get_value(child, 'step', cfg%wave_vector%step, stat=stat)
      if (stat /= 0) cfg%wave_vector%step = 0.01_dp
      call get_value(child, 'nsteps', cfg%wave_vector%nsteps, stat=stat)
      if (stat /= 0) cfg%wave_vector%nsteps = 100
    else
      ! Defaults if section missing
      cfg%wave_vector%mode = 'k0'
      cfg%wave_vector%max = 0.0_dp
      cfg%wave_vector%step = 0.01_dp
      cfg%wave_vector%nsteps = 100
    end if

    ! ---- [bands] ----
    call get_value(table, 'bands', child, stat=stat)
    if (stat == 0 .and. associated(child)) then
      call get_value(child, 'num_cb', cfg%bands%num_cb, stat=stat)
      if (stat /= 0) cfg%bands%num_cb = 2
      call get_value(child, 'num_vb', cfg%bands%num_vb, stat=stat)
      if (stat /= 0) cfg%bands%num_vb = 6
    else
      cfg%bands%num_cb = 2
      cfg%bands%num_vb = 6
    end if
    cfg%evnum = cfg%bands%num_cb + cfg%bands%num_vb

    ! ---- Mode-specific geometry ----
    select case (trim(cfg%confinement))

    case ('bulk')
      call parse_materials_bulk(table, cfg)
      cfg%conf_dir = 'n'
      cfg%ngrid = 1
      cfg%totalSize = 0.0_dp
      cfg%delta = 0.0_dp
      cfg%dz = 0.0_dp

    case ('qw')
      call parse_materials_qw(table, cfg)
      cfg%conf_dir = 'z'

    case ('wire')
      call parse_wire(table, cfg)
      cfg%conf_dir = 'z'

    case ('landau')
      call parse_landau(table, cfg)
      cfg%conf_dir = 'x'

    case default
      print *, 'Error: Unknown confinement: ', trim(cfg%confinement)
      print *, '  Supported: bulk, qw, wire, landau'
      stop 1

    end select

    ! ---- Optional physics sections ----
    call parse_external_field(table, cfg)
    call parse_b_field(table, cfg)
    call parse_sc(table, cfg)
    call parse_topology(table, cfg)
    call parse_bdg(table, cfg)
    call parse_optics(table, cfg)
    call parse_exciton(table, cfg)
    call parse_scattering(table, cfg)
    call parse_feast(table, cfg)
    call parse_strain(table, cfg)
    call parse_gfactor(table, cfg)

    ! ---- Material parameter database ----
    call paramDatabase(cfg%materialN, cfg%num_layers, cfg%params)

    ! ---- Bulk strain substrate (optional) ----
    call get_value(table, 'strain_substrate', cfg%params(1)%strainSubstrate, &
      0.0_dp, stat=stat)

    ! ---- Initialize spatial grid ----
    call init_grid_from_config(cfg)

    ! Wire geometry initialization (2D cut-cells)
    if (cfg%confinement == 'wire') then
      call init_wire_from_config(cfg)
    end if

    ! ---- Final validation ----
    call cfg%validate()

  end subroutine read_config

  ! ==================================================================
  ! Bulk material parsing: single [[material]] entry
  ! ==================================================================
  subroutine parse_materials_bulk(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg

    type(toml_array), pointer :: materials => null()
    type(toml_table), pointer :: mat => null()
    character(len=:), allocatable :: name_val
    integer :: stat, nmat

    call get_value(table, 'material', materials, stat=stat)
    if (.not. associated(materials)) then
      print *, 'Error: [[material]] section required for bulk mode'
      stop 1
    end if

    nmat = tomlf_len(materials)
    if (nmat < 1) then
      print *, 'Error: bulk mode requires at least 1 [[material]] entry'
      stop 1
    end if

    cfg%num_layers = 1
    allocate(cfg%materialN(1))
    allocate(cfg%material_names(1))
    allocate(cfg%params(1))
    allocate(cfg%startPos(1))
    allocate(cfg%endPos(1))
    allocate(cfg%z_min(1))
    allocate(cfg%z_max(1))
    cfg%startPos(1) = 0.0_dp
    cfg%endPos(1) = 1.0_dp
    cfg%z_min(1) = 0.0_dp
    cfg%z_max(1) = 1.0_dp

    call get_value(materials, 1, mat, stat=stat)
    if (.not. associated(mat)) then
      print *, 'Error: Failed to read [[material]] entry 1'
      stop 1
    end if
    call require_string(mat, 'name', name_val, 'material')
    cfg%materialN(1) = trim(name_val)
    cfg%material_names(1) = trim(name_val)
  end subroutine parse_materials_bulk

  ! ==================================================================
  ! QW material parsing: [[material]] with z_min/z_max + grid computation
  ! ==================================================================
  subroutine parse_materials_qw(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg

    type(toml_array), pointer :: materials => null()
    type(toml_table), pointer :: mat => null()
    character(len=:), allocatable :: name_val
    integer :: stat, nmat, i

    call get_value(table, 'material', materials, stat=stat)
    if (.not. associated(materials)) then
      print *, 'Error: [[material]] section required for QW mode'
      stop 1
    end if

    nmat = tomlf_len(materials)
    if (nmat < 1) then
      print *, 'Error: QW mode requires at least 1 [[material]] entry'
      stop 1
    end if

    cfg%num_layers = nmat
    allocate(cfg%materialN(nmat))
    allocate(cfg%material_names(nmat))
    allocate(cfg%params(nmat))
    allocate(cfg%startPos(nmat))
    allocate(cfg%endPos(nmat))
    allocate(cfg%z_min(nmat))
    allocate(cfg%z_max(nmat))
    allocate(cfg%int_start_pos(nmat))
    allocate(cfg%int_end_pos(nmat))
    allocate(cfg%intStartPos(nmat))
    allocate(cfg%intEndPos(nmat))

    do i = 1, nmat
      call get_value(materials, i, mat, stat=stat)
      if (.not. associated(mat)) then
        print *, 'Error: Failed to read [[material]] entry', i
        stop 1
      end if
      call require_string(mat, 'name', name_val, 'material')
      cfg%materialN(i) = trim(name_val)
      cfg%material_names(i) = trim(name_val)
      call get_value(mat, 'z_min', cfg%z_min(i), 0.0_dp, stat=stat)
      call get_value(mat, 'z_max', cfg%z_max(i), 0.0_dp, stat=stat)
      cfg%startPos(i) = cfg%z_min(i)
      cfg%endPos(i) = cfg%z_max(i)
    end do

    ! Validate fd_step minimum for QW mode
    if (cfg%fd_step < 3) then
      print *, 'Error: QW mode requires fd_step >= 3, got ', cfg%fd_step
      stop 1
    end if

    ! Validate fd_step >= FDorder + 1 for QW mode
    if (cfg%fd_step < cfg%FDorder + 1) then
      print *, 'Error: QW fd_step must be >= FDorder + 1, got fd_step=', cfg%fd_step, &
        ' FDorder=', cfg%FDorder
      stop 1
    end if

    ! Grid computation
    cfg%ngrid = cfg%fd_step
    cfg%totalSize = cfg%endPos(1) - cfg%startPos(1)
    cfg%delta = cfg%totalSize / real(cfg%ngrid - 1, kind=dp)
    cfg%z = [ (cfg%startPos(1) + (i-1)*cfg%delta, i=1, cfg%ngrid) ]
    cfg%dz = cfg%delta

    do i = 1, cfg%num_layers
      cfg%int_start_pos(i) = int(abs((cfg%startPos(1) - cfg%startPos(i)) / (cfg%totalSize/(cfg%ngrid-1))))
      cfg%int_end_pos(i) = int(cfg%int_start_pos(1) + ((cfg%endPos(1) + cfg%endPos(i)) / cfg%delta))
    end do
    cfg%int_start_pos = cfg%int_start_pos + 1
    cfg%int_end_pos = cfg%int_end_pos + 1
    cfg%intStartPos = cfg%int_start_pos
    cfg%intEndPos = cfg%int_end_pos
  end subroutine parse_materials_qw

  ! ==================================================================
  ! Wire mode parsing
  ! ==================================================================
  subroutine parse_wire(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg

    type(toml_table), pointer :: wire_tbl => null()
    type(toml_table), pointer :: geom_tbl => null()
    type(toml_array), pointer :: regions_arr => null()
    type(toml_table), pointer :: reg => null()
    type(toml_array), pointer :: verts_arr => null()
    character(len=:), allocatable :: str_val, mat_name
    integer :: stat, i, nverts, nreg

    call get_value(table, 'wire', wire_tbl, stat=stat)
    if (.not. associated(wire_tbl)) then
      print *, 'Error: [wire] section required for wire mode'
      stop 1
    end if

    call require_int(wire_tbl, 'nx', cfg%wire%nx, 'wire')
    call require_int(wire_tbl, 'ny', cfg%wire%ny, 'wire')
    call require_real(wire_tbl, 'dx', cfg%wire%dx, 'wire')
    call require_real(wire_tbl, 'dy', cfg%wire%dy, 'wire')

    ! Wire geometry
    call get_value(wire_tbl, 'geometry', geom_tbl, stat=stat)
    if (associated(geom_tbl)) then
      call require_string(geom_tbl, 'shape', str_val, 'wire.geometry')
      cfg%wire%geom%shape = trim(str_val)

      select case (trim(cfg%wire%geom%shape))
      case ('circle', 'hexagon')
        call require_real(geom_tbl, 'radius', cfg%wire%geom%radius, 'wire.geometry')
      case ('rectangle')
        call require_real(geom_tbl, 'width', cfg%wire%geom%width, 'wire.geometry')
        call require_real(geom_tbl, 'height', cfg%wire%geom%height, 'wire.geometry')
      case ('polygon')
        call get_value(geom_tbl, 'vertices', verts_arr, stat=stat)
        if (.not. associated(verts_arr)) then
          print *, 'Error: polygon shape requires vertices array'
          stop 1
        end if
        nverts = tomlf_len(verts_arr)
        cfg%wire%geom%nverts = nverts
        allocate(cfg%wire%geom%verts(2, nverts))
        do i = 1, nverts
          ! Parse each vertex as a 2-element array
          block
            type(toml_array), pointer :: vert => null()
            call get_value(verts_arr, i, vert, stat=stat)
            if (associated(vert)) then
              call get_value(vert, 1, cfg%wire%geom%verts(1, i), stat=stat)
              call get_value(vert, 2, cfg%wire%geom%verts(2, i), stat=stat)
            end if
          end block
        end do
      case default
        print *, 'Error: Unknown wire shape: ', trim(cfg%wire%geom%shape)
        stop 1
      end select
    end if

    ! Regions
    call get_value(table, 'region', regions_arr, stat=stat)
    if (associated(regions_arr)) then
      nreg = tomlf_len(regions_arr)
      cfg%wire%num_regions = nreg
      cfg%num_layers = nreg
      allocate(cfg%wire%regions(nreg))
      allocate(cfg%materialN(nreg))
      allocate(cfg%material_names(nreg))
      allocate(cfg%params(nreg))
      allocate(cfg%startPos(1))
      allocate(cfg%endPos(1))
      allocate(cfg%z_min(nreg))
      allocate(cfg%z_max(nreg))
      do i = 1, nreg
        call get_value(regions_arr, i, reg, stat=stat)
        if (.not. associated(reg)) then
          print *, 'Error: Failed to read [[region]] entry', i
          stop 1
        end if
        call require_string(reg, 'material', mat_name, 'region')
        cfg%wire%regions(i)%material = trim(mat_name)
        call get_value(reg, 'inner', cfg%wire%regions(i)%inner, 0.0_dp, stat=stat)
        call get_value(reg, 'outer', cfg%wire%regions(i)%outer, 0.0_dp, stat=stat)
        cfg%materialN(i) = trim(mat_name)
        cfg%material_names(i) = trim(mat_name)
        cfg%z_min(i) = cfg%wire%regions(i)%inner
        cfg%z_max(i) = cfg%wire%regions(i)%outer
      end do
      cfg%startPos(1) = 0.0_dp
      cfg%endPos(1) = real(cfg%wire%ny - 1, kind=dp) * cfg%wire%dy
    end if

    ! Validation
    if (cfg%wire%nx < 3) then
      print *, 'Error: wire nx must be >= 3'
      stop 1
    end if
    if (cfg%wire%ny < 3) then
      print *, 'Error: wire ny must be >= 3'
      stop 1
    end if
    if (cfg%wire%nx < cfg%FDorder + 1) then
      print *, 'Error: wire nx must be >= FDorder + 1, got nx=', cfg%wire%nx, &
        ' FDorder=', cfg%FDorder
      stop 1
    end if
    if (cfg%wire%ny < cfg%FDorder + 1) then
      print *, 'Error: wire ny must be >= FDorder + 1, got ny=', cfg%wire%ny, &
        ' FDorder=', cfg%FDorder
      stop 1
    end if

    ! Computed grid size
    cfg%ngrid = cfg%wire%ny
    cfg%dz = cfg%wire%dy
  end subroutine parse_wire

  ! ==================================================================
  ! Landau mode parsing
  ! ==================================================================
  subroutine parse_landau(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg

    type(toml_table), pointer :: landau_tbl => null()
    type(toml_array), pointer :: materials => null()
    type(toml_table), pointer :: mat => null()
    character(len=:), allocatable :: name_val, sweep_val
    integer :: stat, i

    call get_value(table, 'landau', landau_tbl, stat=stat)
    if (.not. associated(landau_tbl)) then
      print *, 'Error: [landau] section required for landau mode'
      stop 1
    end if

    call require_int(landau_tbl, 'nx', cfg%landau%nx, 'landau')
    call require_real(landau_tbl, 'width', cfg%landau%width, 'landau')
    call get_value(landau_tbl, 'sweep', sweep_val, 'ky', stat=stat)
    if (allocated(sweep_val)) then
      cfg%landau%sweep = trim(sweep_val)
    else
      cfg%landau%sweep = 'ky'
    end if

    ! Material
    call get_value(table, 'material', materials, stat=stat)
    if (.not. associated(materials)) then
      print *, 'Error: [[material]] section required for landau mode'
      stop 1
    end if

    cfg%num_layers = 1
    allocate(cfg%materialN(1))
    allocate(cfg%material_names(1))
    allocate(cfg%params(1))
    allocate(cfg%startPos(1))
    allocate(cfg%endPos(1))
    allocate(cfg%z_min(1))
    allocate(cfg%z_max(1))
    allocate(cfg%int_start_pos(1))
    allocate(cfg%int_end_pos(1))
    allocate(cfg%intStartPos(1))
    allocate(cfg%intEndPos(1))

    call get_value(materials, 1, mat, stat=stat)
    if (.not. associated(mat)) then
      print *, 'Error: Failed to read [[material]] entry for landau'
      stop 1
    end if
    call require_string(mat, 'name', name_val, 'material')
    cfg%materialN(1) = trim(name_val)
    cfg%material_names(1) = trim(name_val)

    ! Validation
    if (cfg%landau%nx < 3) then
      print *, 'Error: landau nx must be >= 3'
      stop 1
    end if
    if (cfg%landau%width <= 0.0_dp) then
      print *, 'Error: landau width must be > 0'
      stop 1
    end if
    if (cfg%landau%nx < cfg%FDorder + 1) then
      print *, 'Error: landau nx must be >= FDorder + 1, got nx=', cfg%landau%nx, &
        ' FDorder=', cfg%FDorder
      stop 1
    end if

    ! Computed fields
    cfg%ngrid = cfg%landau%nx
    cfg%dz = cfg%landau%width / real(cfg%landau%nx - 1, kind=dp)
    cfg%totalSize = cfg%landau%width
    cfg%z = [ ((i - 1) * cfg%dz - 0.5_dp * cfg%landau%width, i=1, cfg%landau%nx) ]

    cfg%int_start_pos(1) = 1
    cfg%int_end_pos(1) = cfg%landau%nx
    cfg%intStartPos(1) = 1
    cfg%intEndPos(1) = cfg%landau%nx
    cfg%startPos(1) = 0.0_dp
    cfg%endPos(1) = cfg%landau%width
    cfg%z_min(1) = 0.0_dp
    cfg%z_max(1) = cfg%landau%width
  end subroutine parse_landau

  ! ==================================================================
  ! [external_field] section
  ! ==================================================================
  subroutine parse_external_field(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg

    type(toml_table), pointer :: ef_tbl => null()
    character(len=:), allocatable :: type_val
    integer :: stat

    call get_value(table, 'external_field', ef_tbl, stat=stat)
    if (.not. associated(ef_tbl)) return

    cfg%external_field%enabled = .true.
    call get_value(ef_tbl, 'type', type_val, 'EF', stat=stat)
    if (allocated(type_val)) then
      cfg%external_field%type = trim(type_val)
    end if
    call get_value(ef_tbl, 'value', cfg%external_field%value, 0.0_dp, stat=stat)
  end subroutine parse_external_field

  ! ==================================================================
  ! [b_field] section
  ! ==================================================================
  subroutine parse_b_field(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg

    type(toml_table), pointer :: bf_tbl => null()
    type(toml_array), pointer :: comp_arr => null()
    integer :: stat

    call get_value(table, 'b_field', bf_tbl, stat=stat)
    if (.not. associated(bf_tbl)) return

    ! Parse components array
    call get_value(bf_tbl, 'components', comp_arr, stat=stat)
    if (associated(comp_arr)) then
      call get_value(comp_arr, 1, cfg%b_field%components(1), stat=stat)
      if (stat /= 0) cfg%b_field%components(1) = 0.0_dp
      call get_value(comp_arr, 2, cfg%b_field%components(2), stat=stat)
      if (stat /= 0) cfg%b_field%components(2) = 0.0_dp
      call get_value(comp_arr, 3, cfg%b_field%components(3), stat=stat)
      if (stat /= 0) cfg%b_field%components(3) = 0.0_dp
    end if
    call get_value(bf_tbl, 'g_factor', cfg%b_field%g_factor, 2.0_dp, stat=stat)

    ! Copy to bdg%B_vec for legacy compat
    cfg%bdg%B_vec = cfg%b_field%components
    cfg%bdg%g_factor = cfg%b_field%g_factor
    if (any(cfg%b_field%components /= 0.0_dp) .and. trim(cfg%confinement) /= 'landau') then
      cfg%bdg%enabled = .true.
    end if
  end subroutine parse_b_field

  ! ==================================================================
  ! [sc] section
  ! ==================================================================
  subroutine parse_sc(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg

    type(toml_table), pointer :: sc_tbl => null()
    type(toml_array), pointer :: doping_arr => null()
    type(toml_table), pointer :: dop => null()
    character(len=:), allocatable :: fermi_mode_str, bc_str
    integer :: stat, i, ndop

    call get_value(table, 'sc', sc_tbl, stat=stat)
    if (.not. associated(sc_tbl)) return

    cfg%sc%enabled = 1
    call get_value(sc_tbl, 'max_iterations', cfg%sc%max_iterations, 100, stat=stat)
    call get_value(sc_tbl, 'tolerance', cfg%sc%tolerance, 1.0e-6_dp, stat=stat)
    call get_value(sc_tbl, 'mixing_alpha', cfg%sc%mixing_alpha, 0.3_dp, stat=stat)
    call get_value(sc_tbl, 'diis_history', cfg%sc%diis_history, 7, stat=stat)
    call get_value(sc_tbl, 'temperature', cfg%sc%temperature, 300.0_dp, stat=stat)
    call get_value(sc_tbl, 'fermi_level', cfg%sc%fermi_level, 0.0_dp, stat=stat)
    call get_value(sc_tbl, 'num_kpar', cfg%sc%num_kpar, 201, stat=stat)
    call get_value(sc_tbl, 'kpar_max', cfg%sc%kpar_max, 0.0_dp, stat=stat)

    call get_value(sc_tbl, 'fermi_mode', fermi_mode_str, 'charge_neutrality', stat=stat)
    if (allocated(fermi_mode_str)) then
      select case (trim(fermi_mode_str))
      case ('charge_neutrality')
        cfg%sc%fermi_mode = 0
      case ('fixed')
        cfg%sc%fermi_mode = 1
      case default
        cfg%sc%fermi_mode = 0
      end select
    end if

    call get_value(sc_tbl, 'bc_type', bc_str, 'DD', stat=stat)
    if (allocated(bc_str)) then
      cfg%sc%bc_type = trim(bc_str)
    end if

    call get_value(sc_tbl, 'bc_left', cfg%sc%bc_left, 0.0_dp, stat=stat)
    call get_value(sc_tbl, 'bc_right', cfg%sc%bc_right, 0.0_dp, stat=stat)

    ! Doping
    call get_value(table, 'doping', doping_arr, stat=stat)
    if (associated(doping_arr)) then
      ndop = tomlf_len(doping_arr)
      ! Allocate by num_layers (not ndop) to match layer indexing
      allocate(cfg%doping(cfg%num_layers))
      do i = 1, min(ndop, cfg%num_layers)
        call get_value(doping_arr, i, dop, stat=stat)
        if (.not. associated(dop)) cycle
        call get_value(dop, 'ND', cfg%doping(i)%ND, 0.0_dp, stat=stat)
        call get_value(dop, 'NA', cfg%doping(i)%NA, 0.0_dp, stat=stat)
        call get_value(dop, 'NS', cfg%doping(i)%NS, 0.0_dp, stat=stat)
        call get_value(dop, 'fwhm', cfg%doping(i)%delta_fwhm, 10.0_dp, stat=stat)
        call get_value(dop, 'pos', cfg%doping(i)%delta_pos, 0.0_dp, stat=stat)
        block
          character(len=:), allocatable :: dtype_str
          call get_value(dop, 'type', dtype_str, 'uniform', stat=stat)
          if (allocated(dtype_str)) then
            cfg%doping(i)%dtype = trim(dtype_str)
          end if
        end block
      end do
    end if
  end subroutine parse_sc

  ! ==================================================================
  ! [topology] section
  ! ==================================================================
  subroutine parse_topology(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg

    type(toml_table), pointer :: topo_tbl => null()
    type(toml_table), pointer :: sub_tbl => null()
    type(toml_array), pointer :: range_arr => null()
    character(len=:), allocatable :: mode_val
    integer :: stat

    call get_value(table, 'topology', topo_tbl, stat=stat)
    if (.not. associated(topo_tbl)) return

    cfg%topo%enabled = .true.
    call get_value(topo_tbl, 'mode', mode_val, 'qhe', stat=stat)
    if (allocated(mode_val)) then
      cfg%topo%mode = trim(mode_val)
    end if
    call get_value(topo_tbl, 'compute_chern', cfg%topo%compute_chern, .false., stat=stat)
    call get_value(topo_tbl, 'compute_hall', cfg%topo%compute_hall, .false., stat=stat)
    call get_value(topo_tbl, 'qwz_u', cfg%topo%qwz_u, 0.0_dp, stat=stat)
    call get_value(topo_tbl, 'compute_z2', cfg%topo%compute_z2, .false., stat=stat)
    call get_value(topo_tbl, 'extract_edge_states', cfg%topo%extract_edge_states, .false., stat=stat)
    call get_value(topo_tbl, 'edge_E_window', cfg%topo%edge_E_window, 0.01_dp, stat=stat)
    call get_value(topo_tbl, 'compute_ldos', cfg%topo%compute_ldos, .false., stat=stat)
    call get_value(topo_tbl, 'ldos_eta', cfg%topo%ldos_eta, 0.001_dp, stat=stat)
    call get_value(topo_tbl, 'ldos_num_E', cfg%topo%ldos_num_E, 200, stat=stat)

    ! LDOS E range
    call get_value(topo_tbl, 'ldos_E_range', range_arr, stat=stat)
    if (associated(range_arr)) then
      call get_value(range_arr, 1, cfg%topo%ldos_E_range(1), stat=stat)
      if (stat /= 0) cfg%topo%ldos_E_range(1) = -0.1_dp
      call get_value(range_arr, 2, cfg%topo%ldos_E_range(2), stat=stat)
      if (stat /= 0) cfg%topo%ldos_E_range(2) = 0.1_dp
    end if

    ! Gap sweep
    call get_value(topo_tbl, 'compute_gap_sweep', cfg%topo%compute_gap_sweep, .false., stat=stat)
    call get_value(topo_tbl, 'gap_sweep_B_min', cfg%topo%gap_sweep_B_min, 0.0_dp, stat=stat)
    call get_value(topo_tbl, 'gap_sweep_B_max', cfg%topo%gap_sweep_B_max, 1.0_dp, stat=stat)
    call get_value(topo_tbl, 'gap_sweep_nB', cfg%topo%gap_sweep_nB, 20, stat=stat)
    call get_value(topo_tbl, 'gap_sweep_mu_min', cfg%topo%gap_sweep_mu_min, 0.0_dp, stat=stat)
    call get_value(topo_tbl, 'gap_sweep_mu_max', cfg%topo%gap_sweep_mu_max, 0.01_dp, stat=stat)
    call get_value(topo_tbl, 'gap_sweep_nMu', cfg%topo%gap_sweep_nMu, 20, stat=stat)
    block
      character(len=:), allocatable :: sweep_model_val
      call get_value(topo_tbl, 'sweep_model', sweep_model_val, 'bhz_analytic', stat=stat)
      if (allocated(sweep_model_val)) then
        cfg%topo%sweep_model = trim(sweep_model_val)
      end if
    end block

    ! Conductance
    call get_value(topo_tbl, 'compute_conductance', cfg%topo%compute_conductance, .false., stat=stat)
    block
      character(len=:), allocatable :: cond_method_val
      call get_value(topo_tbl, 'conductance_method', cond_method_val, 'kubo_chern', stat=stat)
      if (allocated(cond_method_val)) then
        cfg%topo%conductance_method = trim(cond_method_val)
      end if
    end block
    call get_value(topo_tbl, 'berry_nk', cfg%topo%berry_nk, 50, stat=stat)
    call get_value(topo_tbl, 'landauer_energy', cfg%topo%landauer_energy, 0.0_dp, stat=stat)

    ! Spectral function
    call get_value(topo_tbl, 'compute_spectral', cfg%topo%compute_spectral, .false., stat=stat)
    call get_value(topo_tbl, 'spectral_k_min', cfg%topo%spectral_k_min, -0.1_dp, stat=stat)
    call get_value(topo_tbl, 'spectral_k_max', cfg%topo%spectral_k_max, 0.1_dp, stat=stat)
    call get_value(topo_tbl, 'spectral_nk', cfg%topo%spectral_nk, 100, stat=stat)
    call get_value(topo_tbl, 'spectral_E_min', cfg%topo%spectral_E_min, -0.05_dp, stat=stat)
    call get_value(topo_tbl, 'spectral_E_max', cfg%topo%spectral_E_max, 0.05_dp, stat=stat)
    call get_value(topo_tbl, 'spectral_nE', cfg%topo%spectral_nE, 200, stat=stat)
    call get_value(topo_tbl, 'spectral_eta', cfg%topo%spectral_eta, 0.001_dp, stat=stat)
  end subroutine parse_topology

  ! ==================================================================
  ! [bdg] section
  ! ==================================================================
  subroutine parse_bdg(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg

    type(toml_table), pointer :: bdg_tbl => null()
    type(toml_array), pointer :: b_sweep_arr => null()
    type(toml_array), pointer :: bvec_arr => null()
    character(len=:), allocatable :: gauge_val
    integer :: stat, i

    call get_value(table, 'bdg', bdg_tbl, stat=stat)
    if (.not. associated(bdg_tbl)) return

    call get_value(bdg_tbl, 'mu', cfg%bdg%mu, 0.0_dp, stat=stat)
    call get_value(bdg_tbl, 'delta_0', cfg%bdg%delta_0, 0.0_dp, stat=stat)
    call get_value(bdg_tbl, 'g_factor', cfg%bdg%g_factor, 2.0_dp, stat=stat)
    call get_value(bdg_tbl, 'gauge', gauge_val, 'landau_x', stat=stat)
    if (allocated(gauge_val)) then
      cfg%bdg%gauge = trim(gauge_val)
    end if
    call get_value(bdg_tbl, 'kz', cfg%bdg%kz, 0.0_dp, stat=stat)

    ! Read B_vec from [bdg] if present (overrides b_field copy)
    call get_value(bdg_tbl, 'B_vec', bvec_arr, stat=stat)
    if (associated(bvec_arr)) then
      do i = 1, min(3, tomlf_len(bvec_arr))
        call get_value(bvec_arr, i, cfg%bdg%B_vec(i), stat=stat)
      end do
    end if

    ! Only enable BdG if pairing gap is nonzero
    if (cfg%bdg%delta_0 > 0.0_dp) then
      cfg%bdg%enabled = .true.
    end if

    ! B sweep
    call get_value(bdg_tbl, 'B_sweep', b_sweep_arr, stat=stat)
    if (associated(b_sweep_arr)) then
      call get_value(b_sweep_arr, 1, cfg%bdg%B_sweep(1), stat=stat)
      if (stat /= 0) cfg%bdg%B_sweep(1) = 0.0_dp
      call get_value(b_sweep_arr, 2, cfg%bdg%B_sweep(2), stat=stat)
      if (stat /= 0) cfg%bdg%B_sweep(2) = 0.0_dp
      call get_value(b_sweep_arr, 3, cfg%bdg%B_sweep(3), stat=stat)
      if (stat /= 0) cfg%bdg%B_sweep(3) = 0.0_dp
    end if
  end subroutine parse_bdg

  ! ==================================================================
  ! [optics] section
  ! ==================================================================
  subroutine parse_optics(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg

    type(toml_table), pointer :: opt_tbl => null()
    type(toml_array), pointer :: e_range_arr => null()
    integer :: stat

    call get_value(table, 'optics', opt_tbl, stat=stat)
    if (.not. associated(opt_tbl)) return

    cfg%optics%enabled = .true.
    call get_value(opt_tbl, 'linewidth_lorentzian', cfg%optics%linewidth_lorentzian, 0.030_dp, stat=stat)
    call get_value(opt_tbl, 'linewidth_gaussian', cfg%optics%linewidth_gaussian, 0.005_dp, stat=stat)
    call get_value(opt_tbl, 'refractive_index', cfg%optics%refractive_index, 3.3_dp, stat=stat)
    call get_value(opt_tbl, 'temperature', cfg%optics%temperature, 300.0_dp, stat=stat)
    call get_value(opt_tbl, 'num_energy_points', cfg%optics%num_energy_points, 200, stat=stat)
    call get_value(opt_tbl, 'E_min', cfg%optics%E_min, 0.5_dp, stat=stat)
    call get_value(opt_tbl, 'E_max', cfg%optics%E_max, 2.0_dp, stat=stat)
    call get_value(opt_tbl, 'carrier_density', cfg%optics%carrier_density, 0.0_dp, stat=stat)
    call get_value(opt_tbl, 'gain_enabled', cfg%optics%gain_enabled, .false., stat=stat)
    call get_value(opt_tbl, 'gain_carrier_density', cfg%optics%gain_carrier_density, 3.0e12_dp, stat=stat)
    call get_value(opt_tbl, 'ISBT', cfg%optics%isbt_enabled, .false., stat=stat)
    call get_value(opt_tbl, 'spontaneous', cfg%optics%spontaneous_enabled, .false., stat=stat)
    call get_value(opt_tbl, 'spin_resolved', cfg%optics%spin_resolved, .false., stat=stat)
  end subroutine parse_optics

  ! ==================================================================
  ! [exciton] section
  ! ==================================================================
  subroutine parse_exciton(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg

    type(toml_table), pointer :: exc_tbl => null()
    character(len=:), allocatable :: method_val
    integer :: stat

    call get_value(table, 'exciton', exc_tbl, stat=stat)
    if (.not. associated(exc_tbl)) return

    cfg%exciton%enabled = .true.
    call get_value(exc_tbl, 'method', method_val, 'variational', stat=stat)
    if (allocated(method_val)) then
      cfg%exciton%method = trim(method_val)
    end if
  end subroutine parse_exciton

  ! ==================================================================
  ! [scattering] section
  ! ==================================================================
  subroutine parse_scattering(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg

    type(toml_table), pointer :: scat_tbl => null()
    integer :: stat

    call get_value(table, 'scattering', scat_tbl, stat=stat)
    if (.not. associated(scat_tbl)) return

    cfg%scattering%enabled = .true.
    call get_value(scat_tbl, 'phonon_energy', cfg%scattering%phonon_energy, 0.036_dp, stat=stat)
    call get_value(scat_tbl, 'eps_inf', cfg%scattering%eps_inf, 10.9_dp, stat=stat)
    call get_value(scat_tbl, 'eps_0', cfg%scattering%eps_0, 12.9_dp, stat=stat)
  end subroutine parse_scattering

  ! ==================================================================
  ! [feast] section
  ! ==================================================================
  subroutine parse_feast(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg

    type(toml_table), pointer :: feast_tbl => null()
    integer :: stat

    call get_value(table, 'feast', feast_tbl, stat=stat)
    if (.not. associated(feast_tbl)) return

    call get_value(feast_tbl, 'emin', cfg%feast%emin, 0.0_dp, stat=stat)
    call get_value(feast_tbl, 'emax', cfg%feast%emax, 0.0_dp, stat=stat)
    call get_value(feast_tbl, 'm0', cfg%feast%m0, 0, stat=stat)
  end subroutine parse_feast

  ! ==================================================================
  ! [strain] section
  ! ==================================================================
  subroutine parse_strain(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg

    type(toml_table), pointer :: strain_tbl => null()
    character(len=:), allocatable :: ref_val, solver_val
    integer :: stat

    call get_value(table, 'strain', strain_tbl, stat=stat)
    if (.not. associated(strain_tbl)) return

    cfg%strain%enabled = .true.
    call get_value(strain_tbl, 'reference', ref_val, 'substrate', stat=stat)
    if (allocated(ref_val)) then
      cfg%strain%reference = trim(ref_val)
    end if
    call get_value(strain_tbl, 'solver', solver_val, 'pardiso', stat=stat)
    if (allocated(solver_val)) then
      cfg%strain%solver = trim(solver_val)
    end if
    call get_value(strain_tbl, 'piezoelectric', cfg%strain%piezoelectric, .false., stat=stat)
  end subroutine parse_strain

  ! ==================================================================
  ! G-factor fields (top-level)
  ! ==================================================================
  subroutine parse_gfactor(table, cfg)
    type(toml_table), intent(inout) :: table
    type(simulation_config), intent(inout) :: cfg
    integer :: stat

    call get_value(table, 'which_band', cfg%which_band, 0, stat=stat)
    call get_value(table, 'band_idx', cfg%band_idx, 1, stat=stat)
  end subroutine parse_gfactor

  ! ==================================================================
  ! Helper: require a string value (fail-fast)
  ! ==================================================================
  subroutine require_string(table, key, val, section)
    type(toml_table), intent(inout) :: table
    character(len=*), intent(in) :: key
    character(len=:), allocatable, intent(out) :: val
    character(len=*), intent(in) :: section
    integer :: stat

    call get_value(table, key, val, stat=stat)
    if (stat /= 0 .or. .not. allocated(val)) then
      print *, 'Error: [', trim(section), '] requires key ''', trim(key), ''''
      stop 1
    end if
  end subroutine require_string

  ! ==================================================================
  ! Helper: require an integer value (fail-fast)
  ! ==================================================================
  subroutine require_int(table, key, val, section)
    type(toml_table), intent(inout) :: table
    character(len=*), intent(in) :: key
    integer, intent(out) :: val
    character(len=*), intent(in) :: section
    integer :: stat

    call get_value(table, key, val, stat=stat)
    if (stat /= 0) then
      print *, 'Error: [', trim(section), '] requires key ''', trim(key), ''''
      stop 1
    end if
  end subroutine require_int

  ! ==================================================================
  ! Helper: require a real value (fail-fast)
  ! ==================================================================
  subroutine require_real(table, key, val, section)
    type(toml_table), intent(inout) :: table
    character(len=*), intent(in) :: key
    real(kind=dp), intent(out) :: val
    character(len=*), intent(in) :: section
    integer :: stat

    call get_value(table, key, val, stat=stat)
    if (stat /= 0) then
      print *, 'Error: [', trim(section), '] requires key ''', trim(key), ''''
      stop 1
    end if
  end subroutine require_real

end module input_parser
