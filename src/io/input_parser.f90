module input_parser

  use definitions
  use parameters
  use hamiltonianConstructor
  use outputFunctions

  implicit none

  private
  public :: read_and_setup

contains

  subroutine read_and_setup(cfg, profile, kpterms)
    type(simulation_config), intent(out) :: cfg
    real(kind=dp), allocatable, intent(out) :: profile(:,:)
    real(kind=dp), allocatable, intent(out) :: kpterms(:,:,:)

    ! Local variables for file handling
    integer(kind=4) :: data_unit
    integer :: status
    character(len=255) :: data_filename, label

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

    allocate(cfg%params(cfg%numLayers))
    allocate(cfg%materialN(cfg%numLayers))

    ! Validate input combinations
    if (cfg%confinement == 0 .and. cfg%numLayers /= 1) then
      print *, 'Error: bulk mode (confinement=0) requires nlayers=1'
      stop 1
    end if
    if (cfg%confinement < 0 .or. cfg%confinement > 1) then
      print *, 'Error: confinement must be 0 or 1, got:', cfg%confinement
      stop 1
    end if

    if (cfg%confinement == 0 .and. cfg%numLayers == 1) then
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
    end if

    if (cfg%confinement == 1) then

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

      cfg%confDir = 'z'

      do i = 1, cfg%numLayers, 1
        cfg%intStartPos(i) = int(abs( (cfg%startPos(1)-cfg%startPos(i)) / (cfg%totalSize/(cfg%fdStep-1))) )
        cfg%intEndPos(i) = int( cfg%intStartPos(1) + ((cfg%endPos(1)+cfg%endPos(i))/cfg%delta) )
      end do
      cfg%intStartPos = cfg%intStartPos + 1
      cfg%intEndPos = cfg%intEndPos + 1

    end if

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
    if (status == 0 .and. cfg%sc%enabled == 1) then
      print *, trim(label), cfg%sc%enabled

      read(data_unit, *, iostat=status) label, cfg%sc%max_iterations
      if (status /= 0) goto 100
      print *, trim(label), cfg%sc%max_iterations

      read(data_unit, *, iostat=status) label, cfg%sc%tolerance
      if (status /= 0) goto 100
      print *, trim(label), cfg%sc%tolerance

      read(data_unit, *, iostat=status) label, cfg%sc%mixing_alpha
      if (status /= 0) goto 100
      print *, trim(label), cfg%sc%mixing_alpha

      read(data_unit, *, iostat=status) label, cfg%sc%diis_history
      if (status /= 0) goto 100
      print *, trim(label), cfg%sc%diis_history

      read(data_unit, *, iostat=status) label, cfg%sc%temperature
      if (status /= 0) goto 100
      print *, trim(label), cfg%sc%temperature

      read(data_unit, *, iostat=status) label, cfg%sc%fermi_mode
      if (status /= 0) goto 100
      print *, trim(label), cfg%sc%fermi_mode

      read(data_unit, *, iostat=status) label, cfg%sc%fermi_level
      if (status /= 0) goto 100
      print *, trim(label), cfg%sc%fermi_level

      read(data_unit, *, iostat=status) label, cfg%sc%num_kpar
      if (status /= 0) goto 100
      print *, trim(label), cfg%sc%num_kpar

      read(data_unit, *, iostat=status) label, cfg%sc%kpar_max
      if (status /= 0) goto 100
      print *, trim(label), cfg%sc%kpar_max

      read(data_unit, *, iostat=status) label, cfg%sc%bc_type
      if (status /= 0) goto 100
      print *, trim(label), cfg%sc%bc_type

      read(data_unit, *, iostat=status) label, cfg%sc%bc_left
      if (status /= 0) goto 100
      print *, trim(label), cfg%sc%bc_left

      read(data_unit, *, iostat=status) label, cfg%sc%bc_right
      if (status /= 0) goto 100
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
      ! Could not read SC flag — use defaults (SC off)
      status = 0
    end if
100 continue

    print *, ' '

    !----------------------------------------------------------------------------
    ! Material parameter setup
    call paramDatabase(cfg%materialN, cfg%numLayers, cfg%params)

    ! Confinement initialization for QW mode
    if (cfg%confDir == 'z') then
      allocate(kpterms(cfg%fdStep, cfg%fdStep, 10))
      kpterms = 0.0_dp
      call confinementInitialization(cfg%z, cfg%intStartPos, cfg%intEndPos, &
        & cfg%materialN, cfg%numLayers, cfg%params, cfg%confDir, profile, kpterms, &
        & cfg%FDorder)
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
      ! Compute dz for use by callers
      cfg%dz = cfg%z(2) - cfg%z(1)
    end if

    ! Close input file
    close(data_unit)

  end subroutine read_and_setup

end module input_parser
