module input_parser

  use definitions
  use parameters
  use geometry
  use hamiltonianConstructor
  use confinement_init, only: confinementInitialization
  use outputFunctions

  implicit none

  private
  public :: read_and_setup

contains

  pure function to_lower_ascii(text) result(lowered)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lowered
    integer :: i, code

    lowered = text
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lowered(i:i) = achar(code + 32)
      end if
    end do
  end function to_lower_ascii

  subroutine read_next_data_line(data_unit, line, status)
    integer(kind=4), intent(in) :: data_unit
    character(len=*), intent(out) :: line
    integer, intent(out) :: status

    do
      read(data_unit, '(A)', iostat=status) line
      if (status /= 0) return
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '!') cycle
      return
    end do
  end subroutine read_next_data_line

  subroutine read_optional_logical_flag(data_unit, expected_label, value, found, label)
    integer(kind=4), intent(in) :: data_unit
    character(len=*), intent(in) :: expected_label
    logical, intent(inout) :: value
    logical, intent(out) :: found
    character(len=*), intent(out) :: label

    integer :: status
    character(len=512) :: line
    character(len=255) :: label_read
    logical :: value_read

    found = .false.
    call read_next_data_line(data_unit, line, status)
    if (status /= 0) return

    read(line, *, iostat=status) label_read
    if (status /= 0) then
      backspace(data_unit)
      return
    end if

    if (trim(to_lower_ascii(label_read)) /= trim(to_lower_ascii(expected_label))) then
      backspace(data_unit)
      return
    end if

    read(line, *, iostat=status) label_read, value_read
    if (status /= 0) then
      backspace(data_unit)
      return
    end if

    label = label_read
    value = value_read
    found = .true.
  end subroutine read_optional_logical_flag

  subroutine read_optional_real_flag(data_unit, expected_label, value, found, label)
    integer(kind=4), intent(in) :: data_unit
    character(len=*), intent(in) :: expected_label
    real(kind=dp), intent(inout) :: value
    logical, intent(out) :: found
    character(len=*), intent(out) :: label

    integer :: status
    character(len=512) :: line
    character(len=255) :: label_read
    real(kind=dp) :: value_read

    found = .false.
    call read_next_data_line(data_unit, line, status)
    if (status /= 0) return

    read(line, *, iostat=status) label_read
    if (status /= 0) then
      backspace(data_unit)
      return
    end if

    if (trim(to_lower_ascii(label_read)) /= trim(to_lower_ascii(expected_label))) then
      backspace(data_unit)
      return
    end if

    read(line, *, iostat=status) label_read, value_read
    if (status /= 0) then
      backspace(data_unit)
      return
    end if

    label = label_read
    value = value_read
    found = .true.
  end subroutine read_optional_real_flag

  subroutine read_and_setup(cfg, profile, kpterms)
    type(simulation_config), intent(out) :: cfg
    real(kind=dp), allocatable, intent(out) :: profile(:,:)
    real(kind=dp), allocatable, intent(out) :: kpterms(:,:,:)

    ! Local variables for file handling
    integer(kind=4) :: data_unit
    integer :: status
    character(len=255) :: data_filename, label
    logical :: found_optional

    ! Local iteration
    integer :: i

    ! reading input/configuration file
    call get_unit(data_unit)
    data_filename = 'input.cfg'
    open(data_unit, file=data_filename, status="old", iostat=status)
    if (status /= 0) then
      print *, 'Error: Cannot open input file: ', trim(data_filename)
      stop 1
    end if

    read(data_unit, *, iostat=status) label, cfg%waveVector
    if (status /= 0) then
      print *, 'Error: Failed to read waveVector from input.cfg'
      stop 1
    end if
    print *, trim(label), cfg%waveVector
    read(data_unit, *, iostat=status) label, cfg%waveVectorMax
    if (status /= 0) then
      print *, 'Error: Failed to read waveVectorMax from input.cfg'
      stop 1
    end if
    print *, trim(label), cfg%waveVectorMax
    read(data_unit, *, iostat=status) label, cfg%waveVectorStep
    if (status /= 0) then
      print *, 'Error: Failed to read waveVectorStep from input.cfg'
      stop 1
    end if
    print *, trim(label), cfg%waveVectorStep
    read(data_unit, *, iostat=status) label, cfg%confinement
    if (status /= 0) then
      print *, 'Error: Failed to read confinement from input.cfg'
      stop 1
    end if
    print *, trim(label), cfg%confinement
    read(data_unit, *, iostat=status) label, cfg%fdStep
    if (status /= 0) then
      print *, 'Error: Failed to read fdStep from input.cfg'
      stop 1
    end if
    print *, trim(label), cfg%fdStep
    read(data_unit, *, iostat=status) label, cfg%FDorder
    if (status /= 0) then
      print *, 'Error: Failed to read FDorder from input.cfg'
      stop 1
    end if
    print *, trim(label), cfg%FDorder
    read(data_unit, *, iostat=status) label, cfg%numLayers
    if (status /= 0) then
      print *, 'Error: Failed to read numLayers from input.cfg'
      stop 1
    end if
    print *, trim(label), cfg%numLayers

    ! Validate confinement value
    if (cfg%confinement < 0 .or. cfg%confinement > 2) then
      print *, 'Error: confinement must be 0, 1, or 2, got:', cfg%confinement
      stop 1
    end if

    ! QW fdStep guard: ensure enough grid points for finite differences
    if (cfg%confinement == 1 .and. cfg%fdStep < 3) then
      print *, 'Error: QW mode requires fdStep >= 3, got:', cfg%fdStep
      stop 1
    end if

    ! Bulk fdStep guard: force fdStep=1 for bulk mode
    if (cfg%confinement == 0 .and. cfg%fdStep /= 1) then
      print *, 'Warning: bulk mode requires fdStep=1. Forcing fdStep=1.'
      cfg%fdStep = 1
    end if

    ! Validate FDorder
    if (cfg%FDorder /= 2 .and. cfg%FDorder /= 4 .and. cfg%FDorder /= 6 &
        .and. cfg%FDorder /= 8 .and. cfg%FDorder /= 10) then
      print *, 'Error: FDorder must be 2, 4, 6, 8, or 10, got:', cfg%FDorder
      stop 1
    end if

    ! Ensure fdStep is large enough for the requested FD order
    if (cfg%confinement == 1 .and. cfg%fdStep < cfg%FDorder + 1) then
      print *, 'Error: fdStep must be >= FDorder + 1 for QW simulations'
      print *, '  fdStep=', cfg%fdStep, ' FDorder=', cfg%FDorder
      print *, '  (FDorder requires at least ', cfg%FDorder+1, ' grid points)'
      stop 1
    end if

    ! Validate input combinations
    if (cfg%confinement == 0 .and. cfg%numLayers /= 1) then
      print *, 'Error: bulk mode (confinement=0) requires nlayers=1'
      stop 1
    end if

    ! --- Mode-specific material/geometry parsing ---

    allocate(cfg%params(cfg%numLayers))
    allocate(cfg%materialN(cfg%numLayers))

    if (cfg%confinement == 0) then
      ! ---- Bulk mode ----
      read(data_unit, *, iostat=status) label, cfg%materialN(1)
      if (status /= 0) then
        print *, 'Error: Failed to read material name from input.cfg'
        stop 1
      end if
      print *, trim(label), cfg%materialN(1)
      cfg%materialN(1) = trim(cfg%materialN(1))
      cfg%confDir = 'n'
      ! Allocate startPos/endPos with dummy values for bulk mode
      allocate(cfg%startPos(1))
      allocate(cfg%endPos(1))
      cfg%startPos(1) = 0.0_dp
      cfg%endPos(1) = 1.0_dp

    else if (cfg%confinement == 1) then
      ! ---- QW mode (1D confinement along z) ----
      allocate(cfg%startPos(cfg%numLayers))
      allocate(cfg%endPos(cfg%numLayers))
      allocate(cfg%intStartPos(cfg%numLayers))
      allocate(cfg%intEndPos(cfg%numLayers))

      do i = 1, cfg%numLayers, 1
        read(data_unit, *, iostat=status) label, cfg%materialN(i), cfg%startPos(i), cfg%endPos(i)
        if (status /= 0) then
          print *, 'Error: Failed to read layer', i, 'material/positions from input.cfg'
          stop 1
        end if
        print *, trim(label), trim(cfg%materialN(i)), cfg%startPos(i), cfg%endPos(i)
      end do

      cfg%totalSize = cfg%endPos(1) - cfg%startPos(1)
      cfg%delta = cfg%totalSize / real(cfg%fdStep - 1, kind=dp)

      cfg%z = [ (cfg%startPos(1) + (i-1)*cfg%delta, i=1, cfg%fdStep) ]

      cfg%dz = cfg%delta
      cfg%confDir = 'z'

      do i = 1, cfg%numLayers, 1
        cfg%intStartPos(i) = int(abs( (cfg%startPos(1)-cfg%startPos(i)) / (cfg%totalSize/(cfg%fdStep-1))) )
        cfg%intEndPos(i) = int( cfg%intStartPos(1) + ((cfg%endPos(1)+cfg%endPos(i))/cfg%delta) )
      end do
      cfg%intStartPos = cfg%intStartPos + 1
      cfg%intEndPos = cfg%intEndPos + 1

    else if (cfg%confinement == 2) then
      ! ---- Wire mode (2D confinement in x-y plane) ----
      cfg%confDir = 'z'  ! wire axis along z, confinement in x and y



      ! Wire grid parameters
      read(data_unit, *, iostat=status) label, cfg%wire_nx
      if (status /= 0) then
        print *, 'Error: Failed to read wire_nx from input.cfg'
        stop 1
      end if
      print *, trim(label), cfg%wire_nx
      read(data_unit, *, iostat=status) label, cfg%wire_ny
      if (status /= 0) then
        print *, 'Error: Failed to read wire_ny from input.cfg'
        stop 1
      end if
      print *, trim(label), cfg%wire_ny

      read(data_unit, *, iostat=status) label, cfg%wire_dx
      if (status /= 0) then
        print *, 'Error: Failed to read wire_dx from input.cfg'
        stop 1
      end if
      print *, trim(label), cfg%wire_dx

      read(data_unit, *, iostat=status) label, cfg%wire_dy
      if (status /= 0) then
        print *, 'Error: Failed to read wire_dy from input.cfg'
        stop 1
      end if
      print *, trim(label), cfg%wire_dy

      ! Wire shape
      read(data_unit, *, iostat=status) label, cfg%wire_geom%shape
      if (status /= 0) then
        print *, 'Error: Failed to read wire_shape from input.cfg'
        stop 1
      end if
      print *, trim(label), trim(cfg%wire_geom%shape)

      ! Shape-dependent dimensions
      select case (trim(cfg%wire_geom%shape))
      case ('circle', 'hexagon')
        read(data_unit, *, iostat=status) label, cfg%wire_geom%radius
        if (status /= 0) then
          print *, 'Error: Failed to read wire_radius from input.cfg'
          stop 1
        end if
        print *, trim(label), cfg%wire_geom%radius

      case ('rectangle')
        read(data_unit, *, iostat=status) label, cfg%wire_geom%width
        if (status /= 0) then
          print *, 'Error: Failed to read wire_width from input.cfg'
          stop 1
        end if
        print *, trim(label), cfg%wire_geom%width

        read(data_unit, *, iostat=status) label, cfg%wire_geom%height
        if (status /= 0) then
          print *, 'Error: Failed to read wire_height from input.cfg'
          stop 1
        end if
        print *, trim(label), cfg%wire_geom%height

      case ('polygon')
        read(data_unit, *, iostat=status) label, cfg%wire_geom%nverts
        if (status /= 0) then
          print *, 'Error: Failed to read wire_polygon from input.cfg'
          stop 1
        end if
        print *, trim(label), cfg%wire_geom%nverts

        allocate(cfg%wire_geom%verts(2, cfg%wire_geom%nverts))
        do i = 1, cfg%wire_geom%nverts
          read(data_unit, *, iostat=status) label, &
            & cfg%wire_geom%verts(1, i), cfg%wire_geom%verts(2, i)
          if (status /= 0) then
            print *, 'Error: Failed to read wire_vertex', i, 'from input.cfg'
            stop 1
          end if
          print *, trim(label), cfg%wire_geom%verts(1, i), cfg%wire_geom%verts(2, i)
        end do

      case default
        print *, 'Error: Unknown wire shape: ', trim(cfg%wire_geom%shape)
        print *, '  Supported shapes: rectangle, circle, hexagon, polygon'
        stop 1
      end select

      ! Number of material regions
      read(data_unit, *, iostat=status) label, cfg%numRegions
      if (status /= 0) then
        print *, 'Error: Failed to read numRegions from input.cfg'
        stop 1
      end if
      print *, trim(label), cfg%numRegions

      if (cfg%numRegions < 1) then
        print *, 'Error: numRegions must be >= 1 for wire mode'
        stop 1
      end if

      allocate(cfg%regions(cfg%numRegions))

      ! Read region specifications: region = materialName inner outer
      do i = 1, cfg%numRegions
        read(data_unit, *, iostat=status) label, cfg%regions(i)%material, &
          & cfg%regions(i)%inner, cfg%regions(i)%outer
        if (status /= 0) then
          print *, 'Error: Failed to read region', i, 'from input.cfg'
          stop 1
        end if
        print *, trim(label), trim(cfg%regions(i)%material), &
          & cfg%regions(i)%inner, cfg%regions(i)%outer
        cfg%regions(i)%material = trim(cfg%regions(i)%material)
      end do

      ! Build materialN and params from regions for backward compat
      ! (the paramDatabase call expects materialN array)
      if (cfg%numLayers /= cfg%numRegions) then
        ! Reallocate to match regions if numLayers was set differently
        deallocate(cfg%params)
        deallocate(cfg%materialN)
        cfg%numLayers = cfg%numRegions
        allocate(cfg%params(cfg%numRegions))
        allocate(cfg%materialN(cfg%numRegions))
      end if
      do i = 1, cfg%numRegions
        cfg%materialN(i) = cfg%regions(i)%material
      end do

      ! Wire mode grid validation
      if (cfg%wire_nx < 3) then
        print *, 'Error: wire_nx must be >= 3, got:', cfg%wire_nx
        stop 1
      end if
      if (cfg%wire_ny < 3) then
        print *, 'Error: wire_ny must be >= 3, got:', cfg%wire_ny
        stop 1
      end if
      if (cfg%wire_dx <= 0.0_dp) then
        print *, 'Error: wire_dx must be > 0, got:', cfg%wire_dx
        stop 1
      end if
      if (cfg%wire_dy <= 0.0_dp) then
        print *, 'Error: wire_dy must be > 0, got:', cfg%wire_dy
        stop 1
      end if
      ! FD order validation for wire grid
      if (cfg%wire_nx < cfg%FDorder + 1) then
        print *, 'Error: wire_nx must be >= FDorder + 1'
        print *, '  wire_nx=', cfg%wire_nx, ' FDorder=', cfg%FDorder
        stop 1
      end if
      if (cfg%wire_ny < cfg%FDorder + 1) then
        print *, 'Error: wire_ny must be >= FDorder + 1'
        print *, '  wire_ny=', cfg%wire_ny, ' FDorder=', cfg%FDorder
        stop 1
      end if

      ! Set fdStep to wire_ny for backward compat (used by array sizing in callers)
      cfg%fdStep = cfg%wire_ny
      cfg%dz = cfg%wire_dy

      ! Allocate dummy startPos/endPos for routines that expect them
      allocate(cfg%startPos(1))
      allocate(cfg%endPos(1))
      cfg%startPos(1) = 0.0_dp
      cfg%endPos(1) = real(cfg%wire_ny - 1, kind=dp) * cfg%wire_dy

    end if

    ! --- Common fields: numcb, numvb, external field, gfactor, SC ---

    read(data_unit, *, iostat=status) label, cfg%numcb
    if (status /= 0) then
      print *, 'Error: Failed to read numcb from input.cfg'
      stop 1
    end if
    print *, trim(label), cfg%numcb
    read(data_unit, *, iostat=status) label, cfg%numvb
    if (status /= 0) then
      print *, 'Error: Failed to read numvb from input.cfg'
      stop 1
    end if
    print *, trim(label), cfg%numvb
    cfg%evnum = cfg%numcb + cfg%numvb

    read(data_unit, *, iostat=status) label, cfg%ExternalField, cfg%EFtype
    if (status /= 0) then
      print *, 'Error: Failed to read ExternalField from input.cfg'
      stop 1
    end if
    print *, trim(label), cfg%ExternalField, cfg%EFtype
    if (cfg%EFtype == "EF") then
      read(data_unit, *, iostat=status) label, cfg%Evalue
      if (status /= 0) then
        print *, 'Error: Failed to read Evalue from input.cfg'
        stop 1
      end if
      print *, trim(label), cfg%Evalue
    else
      stop "Type of external field not implemented"
    end if

      ! Try reading gfactor params; use defaults if missing (backward compatible)
      read(data_unit, *, iostat=status) label, cfg%whichBand
      if (status == 0) then
        print *, trim(label), cfg%whichBand
        read(data_unit, *, iostat=status) label, cfg%bandIdx
        if (status == 0) print *, trim(label), cfg%bandIdx
      else
        cfg%whichBand = 0
        cfg%bandIdx= 1
      end if

    ! --- SC parameters (backward compatible: uses defaults if missing) ---
    read(data_unit, *, iostat=status) label, cfg%sc%enabled
    sc_block: do
      if (status == 0 .and. cfg%sc%enabled == 1) then
        print *, trim(label), cfg%sc%enabled

        read(data_unit, *, iostat=status) label, cfg%sc%max_iterations
        if (status /= 0) then; status = 0; exit sc_block; end if
        print *, trim(label), cfg%sc%max_iterations

        read(data_unit, *, iostat=status) label, cfg%sc%tolerance
        if (status /= 0) then; status = 0; exit sc_block; end if
        print *, trim(label), cfg%sc%tolerance

        read(data_unit, *, iostat=status) label, cfg%sc%mixing_alpha
        if (status /= 0) then; status = 0; exit sc_block; end if
        print *, trim(label), cfg%sc%mixing_alpha

        read(data_unit, *, iostat=status) label, cfg%sc%diis_history
        if (status /= 0) then; status = 0; exit sc_block; end if
        print *, trim(label), cfg%sc%diis_history

        read(data_unit, *, iostat=status) label, cfg%sc%temperature
        if (status /= 0) then; status = 0; exit sc_block; end if
        print *, trim(label), cfg%sc%temperature

        read(data_unit, *, iostat=status) label, cfg%sc%fermi_mode
        if (status /= 0) then; status = 0; exit sc_block; end if
        print *, trim(label), cfg%sc%fermi_mode

        read(data_unit, *, iostat=status) label, cfg%sc%fermi_level
        if (status /= 0) then; status = 0; exit sc_block; end if
        print *, trim(label), cfg%sc%fermi_level

        read(data_unit, *, iostat=status) label, cfg%sc%num_kpar
        if (status /= 0) then; status = 0; exit sc_block; end if
        print *, trim(label), cfg%sc%num_kpar

        read(data_unit, *, iostat=status) label, cfg%sc%kpar_max
        if (status /= 0) then; status = 0; exit sc_block; end if
        print *, trim(label), cfg%sc%kpar_max

        read(data_unit, *, iostat=status) label, cfg%sc%bc_type
        if (status /= 0) then; status = 0; exit sc_block; end if
        print *, trim(label), cfg%sc%bc_type

        read(data_unit, *, iostat=status) label, cfg%sc%bc_left
        if (status /= 0) then; status = 0; exit sc_block; end if
        print *, trim(label), cfg%sc%bc_left

        read(data_unit, *, iostat=status) label, cfg%sc%bc_right
        if (status /= 0) then; status = 0; exit sc_block; end if
        print *, trim(label), cfg%sc%bc_right

        ! Read doping per layer (ND NA for each layer)
        allocate(cfg%doping(cfg%numLayers))
        do i = 1, cfg%numLayers
          read(data_unit, *, iostat=status) label, cfg%doping(i)%ND, cfg%doping(i)%NA
          if (status /= 0) then
            cfg%doping(i)%ND = 0.0_dp
            cfg%doping(i)%NA = 0.0_dp
            status = 0
            exit
          end if
          print *, trim(label), cfg%doping(i)%ND, cfg%doping(i)%NA
        end do

      else if (status == 0) then
        ! SC=0, skip remaining SC lines
        print *, trim(label), cfg%sc%enabled
      else
        ! Could not read SC flag -- use defaults (SC off)
        status = 0
      end if
      exit sc_block
    end do sc_block

    ! --- Optical spectra parameters (backward compatible: uses defaults if missing) ---
    call read_optional_logical_flag(data_unit, 'optics:', cfg%optics%enabled, found_optional, label)
    optics_block: do
      if (found_optional .and. cfg%optics%enabled) then
        print *, trim(label), cfg%optics%enabled

        read(data_unit, *, iostat=status) label, cfg%optics%linewidth_lorentzian
        if (status /= 0) then; status = 0; exit optics_block; end if
        print *, trim(label), cfg%optics%linewidth_lorentzian

        read(data_unit, *, iostat=status) label, cfg%optics%linewidth_gaussian
        if (status /= 0) then; status = 0; exit optics_block; end if
        print *, trim(label), cfg%optics%linewidth_gaussian

        read(data_unit, *, iostat=status) label, cfg%optics%refractive_index
        if (status /= 0) then; status = 0; exit optics_block; end if
        print *, trim(label), cfg%optics%refractive_index

        read(data_unit, *, iostat=status) label, cfg%optics%E_min
        if (status /= 0) then; status = 0; exit optics_block; end if
        print *, trim(label), cfg%optics%E_min

        read(data_unit, *, iostat=status) label, cfg%optics%E_max
        if (status /= 0) then; status = 0; exit optics_block; end if
        print *, trim(label), cfg%optics%E_max

        read(data_unit, *, iostat=status) label, cfg%optics%num_energy_points
        if (status /= 0) then; status = 0; exit optics_block; end if
        print *, trim(label), cfg%optics%num_energy_points

        read(data_unit, *, iostat=status) label, cfg%optics%temperature
        if (status /= 0) then; status = 0; exit optics_block; end if
        print *, trim(label), cfg%optics%temperature

        read(data_unit, *, iostat=status) label, cfg%optics%carrier_density
        if (status /= 0) then; status = 0; exit optics_block; end if
        print *, trim(label), cfg%optics%carrier_density

        ! Gain parameters
        read(data_unit, *, iostat=status) label, cfg%optics%gain_enabled
        if (status /= 0) then; status = 0; exit optics_block; end if
        print *, trim(label), cfg%optics%gain_enabled

        read(data_unit, *, iostat=status) label, cfg%optics%gain_carrier_density
        if (status /= 0) then; status = 0; exit optics_block; end if
        print *, trim(label), cfg%optics%gain_carrier_density

        ! ISBT
        read(data_unit, *, iostat=status) label, cfg%optics%isbt_enabled
        if (status /= 0) then; status = 0; exit optics_block; end if
        print *, trim(label), cfg%optics%isbt_enabled

        ! Spontaneous emission (optional — backward compatible with configs
        ! that omit these lines).  Peek at the label; if it doesn't match
        ! the expected field name, push the line back with BACKSPACE.
        read(data_unit, *, iostat=status) label
        if (status /= 0) then; status = 0; exit optics_block; end if
        if (trim(adjustl(label)) == 'SpontaneousEnabled:') then
          backspace(data_unit)
          read(data_unit, *, iostat=status) label, cfg%optics%spontaneous_enabled
          if (status /= 0) then; status = 0; exit optics_block; end if
          print *, trim(label), cfg%optics%spontaneous_enabled

          ! Spin resolution (also optional)
          read(data_unit, *, iostat=status) label
          if (status /= 0) then; status = 0; exit optics_block; end if
          if (trim(adjustl(label)) == 'SpinResolved:') then
            backspace(data_unit)
            read(data_unit, *, iostat=status) label, cfg%optics%spin_resolved
            if (status /= 0) then; status = 0; exit optics_block; end if
            print *, trim(label), cfg%optics%spin_resolved
          else
            backspace(data_unit)
          end if
        else
          backspace(data_unit)
        end if

      else if (found_optional) then
        ! Optics=F, skip remaining optics lines
        print *, trim(label), cfg%optics%enabled
      end if
      exit optics_block
    end do optics_block

    ! --- Exciton parameters (backward compatible: uses defaults if missing) ---
    call read_optional_logical_flag(data_unit, 'exciton:', cfg%exciton%enabled, found_optional, label)
    exciton_block: do
      if (found_optional .and. cfg%exciton%enabled) then
        print *, trim(label), cfg%exciton%enabled

        read(data_unit, *, iostat=status) label, cfg%exciton%method
        if (status /= 0) then; status = 0; exit exciton_block; end if
        print *, trim(label), trim(cfg%exciton%method)

      else if (found_optional) then
        print *, trim(label), cfg%exciton%enabled
      end if
      exit exciton_block
    end do exciton_block

    ! --- Scattering parameters (backward compatible: uses defaults if missing) ---
    call read_optional_logical_flag(data_unit, 'scattering:', cfg%scattering%enabled, found_optional, label)
    scattering_block: do
      if (found_optional .and. cfg%scattering%enabled) then
        print *, trim(label), cfg%scattering%enabled

        read(data_unit, *, iostat=status) label, cfg%scattering%phonon_energy
        if (status /= 0) then; status = 0; exit scattering_block; end if
        print *, trim(label), cfg%scattering%phonon_energy

        read(data_unit, *, iostat=status) label, cfg%scattering%eps_inf
        if (status /= 0) then; status = 0; exit scattering_block; end if
        print *, trim(label), cfg%scattering%eps_inf

        read(data_unit, *, iostat=status) label, cfg%scattering%eps_0
        if (status /= 0) then; status = 0; exit scattering_block; end if
        print *, trim(label), cfg%scattering%eps_0

      else if (found_optional) then
        print *, trim(label), cfg%scattering%enabled
      end if
      exit scattering_block
    end do scattering_block

    ! --- FEAST eigensolver tuning (backward compatible: uses defaults if missing) ---
    call read_optional_real_flag(data_unit, 'feast_emin:', cfg%feast_emin, found_optional, label)
    if (found_optional) then
      print *, trim(label), cfg%feast_emin
      read(data_unit, *, iostat=status) label, cfg%feast_emax
      if (status == 0) then
        print *, trim(label), cfg%feast_emax
        read(data_unit, *, iostat=status) label, cfg%feast_m0
        if (status == 0) then
          print *, trim(label), cfg%feast_m0
        else
          cfg%feast_m0 = 0
          status = 0
        end if
      else
        cfg%feast_emax = 0.0_dp
        cfg%feast_m0 = 0
        status = 0
      end if
    else
      ! feast_emin not found -- use auto window
      cfg%feast_emin = 0.0_dp
      cfg%feast_emax = 0.0_dp
      cfg%feast_m0 = 0
    end if

    ! --- Strain parameters (backward compatible: uses defaults if missing) ---
    call read_optional_logical_flag(data_unit, 'strain:', cfg%strain%enabled, found_optional, label)
    strain_block: do
      if (found_optional .and. cfg%strain%enabled) then
        print *, trim(label), cfg%strain%enabled

        read(data_unit, *, iostat=status) label, cfg%strain%reference
        if (status /= 0) then; status = 0; exit strain_block; end if
        print *, trim(label), trim(cfg%strain%reference)

        read(data_unit, *, iostat=status) label, cfg%strain%solver
        if (status /= 0) then; status = 0; exit strain_block; end if
        print *, trim(label), trim(cfg%strain%solver)

        read(data_unit, *, iostat=status) label, cfg%strain%piezoelectric
        if (status /= 0) then; status = 0; exit strain_block; end if
        print *, trim(label), cfg%strain%piezoelectric

      else if (found_optional) then
        ! strain=.false., skip remaining strain lines
        print *, trim(label), cfg%strain%enabled
      end if
      exit strain_block
    end do strain_block

    print *, ' '

    !----------------------------------------------------------------------------
    ! Material parameter setup
    call paramDatabase(cfg%materialN, cfg%numLayers, cfg%params)

    ! --- Bulk strain substrate (backward compatible: 0 = no strain) ---
    ! Must come AFTER paramDatabase so the value is not overwritten.
    read(data_unit, *, iostat=status) label, cfg%params(1)%strainSubstrate
    if (status /= 0) then
      ! Not found -- default to 0 (no bulk strain), reset status
      cfg%params(1)%strainSubstrate = 0.0_dp
      status = 0
    else
      print *, trim(label), cfg%params(1)%strainSubstrate
    end if

    ! Initialize the unified spatial grid from config fields
    call init_grid_from_config(cfg)

    ! For wire mode, initialize geometry (cut-cells, ghost map, material_id)
    if (cfg%confinement == 2) then
      call init_wire_from_config(cfg)
    end if

    ! Confinement initialization for QW mode (not wire)
    if (cfg%confDir == 'z' .and. cfg%confinement == 1) then
      allocate(kpterms(grid_ngrid(cfg%grid), grid_ngrid(cfg%grid), 10))
      kpterms = 0.0_dp
      call confinementInitialization(cfg, profile, kpterms)
      ! Guard: electric field requires z(1) /= 0
      if (cfg%ExternalField == 1 .and. cfg%EFtype == "EF") then
        if (abs(cfg%z(1)) < tolerance) then
          print *, 'Error: Electric field requires z(1) /= 0.'
          print *, '  Adjust startPos/endPos so grid does not start at z=0.'
          stop 1
        end if
      end if
      if (cfg%ExternalField == 1 .and. cfg%EFtype == "EF") then
        call externalFieldSetup_electricField(profile, cfg%Evalue, cfg%totalSize, cfg%z)
      end if
    end if

    ! Close input file
    close(data_unit)

    ! Final validation pass: catch any inconsistent state
    call cfg%validate()

  end subroutine read_and_setup

end module input_parser
