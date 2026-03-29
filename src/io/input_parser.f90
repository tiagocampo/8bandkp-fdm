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

    read(data_unit, *) label, cfg%waveVector
    print *, trim(label), cfg%waveVector
    read(data_unit, *) label, cfg%waveVectorMax
    print *, trim(label), cfg%waveVectorMax
    read(data_unit, *) label, cfg%waveVectorStep
    print *, trim(label), cfg%waveVectorStep
    read(data_unit, *) label, cfg%confinement
    print *, trim(label), cfg%confinement
    read(data_unit, *) label, cfg%fdStep
    print *, trim(label), cfg%fdStep
    read(data_unit, *) label, cfg%numLayers
    print *, trim(label), cfg%numLayers

    ! Bulk fdStep guard: force fdStep=1 for bulk mode
    if (cfg%confinement == 0 .and. cfg%fdStep /= 1) then
      print *, 'Warning: bulk mode requires fdStep=1. Forcing fdStep=1.'
      cfg%fdStep = 1
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
      read(data_unit, *) label, cfg%materialN(1)
      print *, trim(label), cfg%materialN(1)
      cfg%materialN(1) = trim(cfg%materialN(1))
      cfg%confDir = 'n'
    end if

    if (cfg%confinement == 1) then

      allocate(cfg%startPos(cfg%numLayers))
      allocate(cfg%endPos(cfg%numLayers))
      allocate(cfg%intStartPos(cfg%numLayers))
      allocate(cfg%intEndPos(cfg%numLayers))

      do i = 1, cfg%numLayers, 1
        read(data_unit, *) label, cfg%materialN(i), cfg%startPos(i), cfg%endPos(i)
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

    read(data_unit, *) label, cfg%numcb
    print *, trim(label), cfg%numcb
    read(data_unit, *) label, cfg%numvb
    print *, trim(label), cfg%numvb
    cfg%evnum = cfg%numcb + cfg%numvb

    read(data_unit, *) label, cfg%ExternalField, cfg%EFtype
    print *, trim(label), cfg%ExternalField, cfg%EFtype
    if (cfg%EFtype == "EF") then
      read(data_unit, *) label, cfg%Evalue
      print *, trim(label), cfg%Evalue
    else
      stop "Type of external field not implemented"
    end if

    print *, ' '

    !----------------------------------------------------------------------------
    ! Material parameter setup
    call paramDatabase(cfg%materialN, cfg%numLayers, cfg%params)

    ! Confinement initialization for QW mode
    if (cfg%confDir == 'z') then
      allocate(kpterms(cfg%fdStep, cfg%fdStep, 10))
      kpterms = 0.0_dp
      call confinementInitialization(cfg%z, cfg%intStartPos, cfg%intEndPos, &
        & cfg%materialN, cfg%numLayers, cfg%params, cfg%confDir, profile, kpterms)
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
