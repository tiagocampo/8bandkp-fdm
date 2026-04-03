module definitions

  implicit none

  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: qp = selected_real_kind(33, 4931)
  integer, parameter :: iknd = selected_int_kind(8)

  real(kind=sp), parameter :: pi_sp = 3.14159265358979323846264338327950288_sp
  real(kind=dp), parameter :: pi_dp = 3.14159265358979323846264338327950288_dp
  real(kind=qp), parameter :: pi_qp = 3.14159265358979323846264338327950288_qp

  real(kind=dp), parameter :: a0 = 0.5291772_dp !bohr radius in AA
  real(kind=dp), parameter :: hbar = 6.582119514_dp*1E-16! eV s
  real(kind=dp), parameter :: const = 3.8099820104685644_dp !ev AA**2
  real(kind=dp), parameter :: c = 2.99792458_dp*1E8*1E10 ! AA / s
  real(kind=dp), parameter :: m0 = 0.5109989461_dp*1E6/c**2 ! eV / c^2
  real(kind=dp), parameter :: e = 1.602176565_dp*1E-19
  real(kind=dp), parameter :: hbar2O2m0 = (hbar)**2 / (2.*m0) !eV AA**2
  real(kind=dp), parameter :: e0 = 8.854187817_dp*1E-12*(1E-9) ! C V-1 nm-1
  real(kind=dp), parameter :: kB_eV = 8.617333262e-5_dp        ! Boltzmann constant (eV/K)
  real(kind=dp), parameter :: tolerance=1e-7
  logical, parameter :: renormalization = .False.

  ! Band structure constants (8-band basis: 2 CB + 6 VB per k-point)
  integer, parameter :: NUM_CB_STATES = 2
  integer, parameter :: NUM_VB_STATES = 6

  ! Free electron g-factor (CODATA)
  real(kind=dp), parameter :: g_free = 2.00231_dp

  complex(kind=dp), parameter :: IU = cmplx(0.0_dp, 1.0_dp, kind=dp)
  real(kind=dp), parameter :: SQR3 = dsqrt(3.0_dp)
  real(kind=dp), parameter :: SQR2 = dsqrt(2.0_dp)
  real(kind=dp), parameter :: SQR2o3 = dsqrt(2.0_dp/3.0_dp)
  real(kind=dp), parameter :: RQS3 = 1.0_dp/SQR3
  real(kind=dp), parameter :: RQS2 = 1.0_dp/SQR2
  complex(kind=dp), parameter :: ZERO = cmplx(0.0_dp, 0.0_dp, kind=dp)
  complex(kind=dp), parameter :: UM = cmplx(1.0_dp, 0.0_dp, kind=dp)


  type wavevector

    real(kind=dp) :: kx, ky, kz

  end type wavevector

  type paramStruct

    real(kind=dp) :: meff, gamma1, gamma2, gamma3, P, A, deltaSO, EP, Eg, EV, EC
    real(kind=dp) :: eps0  ! static dielectric constant (unitless)

  end type paramStruct

  type doping_spec
    real(kind=dp) :: ND = 0.0_dp    ! donor concentration (cm^-3)
    real(kind=dp) :: NA = 0.0_dp    ! acceptor concentration (cm^-3)
  end type doping_spec

  type sc_config
    integer       :: enabled = 0              ! 0=off, 1=on
    integer       :: max_iterations = 100     ! iteration cap
    real(kind=dp) :: tolerance = 1.0e-6_dp    ! convergence |Delta Phi|_inf in eV
    real(kind=dp) :: mixing_alpha = 0.3_dp    ! linear mixing parameter
    integer       :: diis_history = 7         ! DIIS history length
    real(kind=dp) :: temperature = 300.0_dp   ! Kelvin
    integer       :: fermi_mode = 0           ! 0=charge_neutrality, 1=fixed
    real(kind=dp) :: fermi_level = 0.0_dp     ! fixed mu (eV)
    integer       :: num_kpar = 201           ! k_par sampling points (odd for Simpson)
    real(kind=dp) :: kpar_max = 0.0_dp        ! max k_par (1/AA), 0=auto
    character(len=2) :: bc_type = 'DD'        ! 'DD' or 'DN'
    real(kind=dp) :: bc_left = 0.0_dp         ! left BC potential (eV)
    real(kind=dp) :: bc_right = 0.0_dp        ! right BC potential (eV)
  end type sc_config

  ! ------------------------------------------------------------------
  ! Wire geometry specification for confinement=2.
  ! shape: 'rectangle', 'circle', 'hexagon', 'polygon'
  ! For circle/hexagon: wire_radius is the defining size.
  ! For rectangle: wire_width (y), wire_height (z).
  ! For polygon: wire_polygon_nverts + wire_polygon_verts(:,:).
  ! ------------------------------------------------------------------
  type wire_geometry
    character(len=16)  :: shape = 'rectangle'
    real(kind=dp)      :: radius = 0.0_dp     ! AA (circle/hexagon)
    real(kind=dp)      :: width  = 0.0_dp     ! AA (rectangle, y-extent)
    real(kind=dp)      :: height = 0.0_dp     ! AA (rectangle, z-extent)
    integer            :: nverts = 0           ! polygon vertex count
    real(kind=dp), allocatable :: verts(:,:)  ! (2, nverts) y,z vertices
  end type wire_geometry

  ! ------------------------------------------------------------------
  ! Region specification for wire mode (replaces layer for 2D).
  ! Each region has a material name and a radial extent (inner, outer).
  ! ------------------------------------------------------------------
  type region_spec
    character(len=255) :: material = ''
    real(kind=dp)      :: inner   = 0.0_dp   ! inner radius/distance (AA)
    real(kind=dp)      :: outer   = 0.0_dp   ! outer radius/distance (AA)
  end type region_spec

  ! ------------------------------------------------------------------
  ! Unified spatial grid for bulk (ndim=0), QW (ndim=1), wire (ndim=2).
  !
  ! Bulk:   ny=1, nz=1, dy=dz=0, no coordinate arrays.
  ! QW:     ny=1, nz=fdStep, dz>0, z(:) holds 1D grid.
  ! Wire:   ny>1, nz>1, dy>0, dz>0, coords(:,:) holds 2D grid.
  !
  ! Cut-cell fields are only allocated for wire mode (ndim=2).  For QW
  ! and bulk they remain unallocated; code can test associated() or
  ! treat them as unity when unallocated.
  ! ------------------------------------------------------------------
  type spatial_grid
    integer :: ndim = 0                ! 0=bulk, 1=QW, 2=wire
    integer :: ny = 1                  ! grid points in y (wire axis)
    integer :: nz = 1                  ! grid points in z (growth axis)
    real(kind=dp) :: dy = 0.0_dp       ! grid spacing in y (AA)
    real(kind=dp) :: dz = 0.0_dp       ! grid spacing in z (AA)

    ! 1D coordinate arrays (allocated for ndim >= 1)
    real(kind=dp), allocatable :: y(:)       ! (ny)
    real(kind=dp), allocatable :: z(:)       ! (nz)

    ! Flattened 2D coordinate array (allocated for ndim == 2)
    ! coords(1,:) = y, coords(2,:) = z,  column-major: (j-1)*ny + i
    real(kind=dp), allocatable :: coords(:,:)    ! (2, ny*nz)

    ! Material index at each grid point (1-based layer index)
    integer, allocatable  :: material_id(:)      ! (ny*nz)

    ! Cut-cell immersed boundary fields (ndim == 2 only)
    real(kind=dp), allocatable :: cell_volume(:)       ! (ny*nz) fractional vol [0,1]
    real(kind=dp), allocatable :: face_fraction_y(:,:) ! (ny*nz, 2) left/right y-face
    real(kind=dp), allocatable :: face_fraction_z(:,:) ! (ny*nz, 2) bottom/top z-face
    ! Nearest active neighbor for each ghost/inactive point
    integer, allocatable  :: ghost_map(:,:)            ! (ny*nz, 4) N S W E
  end type spatial_grid

  type simulation_config
    integer :: confinement = 0
    integer :: fdStep = 1
    integer :: FDorder = 2
    integer :: numLayers = 1
    integer :: numcb = 2
    integer :: numvb = 6
    integer :: evnum = 8
    integer :: waveVectorStep = 100
    integer :: ExternalField = 0
    character(len=2) :: waveVector = 'k0'
    character(len=1) :: confDir = 'n'
    character(len=2) :: EFtype = '  '
    real(kind=dp) :: waveVectorMax = 0.0_dp
    real(kind=dp) :: Evalue = 0.0_dp
    real(kind=dp) :: totalSize = 0.0_dp
    real(kind=dp) :: delta = 0.0_dp
    real(kind=dp) :: dz = 0.0_dp
    integer :: whichBand = 0     ! 0=CB g-factor, 1=VB g-factor
    integer :: bandIdx = 1       ! which doublet (1=ground, 2=first excited, ...)
    real(kind=dp), allocatable :: startPos(:)
    real(kind=dp), allocatable :: endPos(:)
    real(kind=dp), allocatable :: z(:)
    integer, allocatable :: intStartPos(:)
    integer, allocatable :: intEndPos(:)
    character(len=255), allocatable :: materialN(:)
    type(paramStruct), allocatable :: params(:)
    type(doping_spec), allocatable :: doping(:)   ! per-layer doping
    type(sc_config)                :: sc           ! SC parameters
    type(spatial_grid)             :: grid         ! unified spatial grid

    ! ---- Wire-specific fields (confinement=2) ----
    integer            :: wire_ny = 0            ! grid points in y
    integer            :: wire_nz = 0            ! grid points in z
    real(kind=dp)      :: wire_dy = 0.0_dp      ! grid spacing y (AA)
    real(kind=dp)      :: wire_dz = 0.0_dp      ! grid spacing z (AA)
    type(wire_geometry) :: wire_geom             ! shape + dimensions
    integer            :: numRegions = 0         ! number of material regions
    type(region_spec), allocatable :: regions(:) ! region specifications
  end type simulation_config

  type group
      integer :: order    ! original order of unsorted data
      real(kind=dp) :: value       ! values to be sorted by
  end type group


  contains

  pure function kronij(i,j)
    integer, intent(in) :: i,j
    integer :: kronij
    kronij = merge(1, 0, i == j)
  end function

  ! ------------------------------------------------------------------
  ! Return total number of spatial grid points (ny * nz).
  ! Works for all ndim: bulk returns 1, QW returns nz, wire returns ny*nz.
  ! ------------------------------------------------------------------
  pure function grid_ngrid(grid) result(n)
    type(spatial_grid), intent(in) :: grid
    integer :: n
    n = grid%ny * grid%nz
  end function

  ! ------------------------------------------------------------------
  ! Initialize a spatial_grid from the fields already set in
  ! simulation_config.  This is the primary grid initialization path.
  !
  ! For bulk (confinement=0): ndim=0, ny=nz=1.
  ! For QW  (confinement=1): ndim=1, ny=1, nz=fdStep, z(:) copied.
  ! For wire(confinement=2): ndim=2, ny=wire_ny, nz=wire_nz,
  !   set grid dimensions only.  Coordinate arrays, material_id,
  !   cut-cell fields and ghost_map are populated by
  !   init_wire_from_config() in the geometry module.
  ! ------------------------------------------------------------------
  subroutine init_grid_from_config(cfg)
    type(simulation_config), intent(inout) :: cfg

    integer :: i, j, ngrid

    select case (cfg%confinement)
    case (0)
      ! Bulk: single point, no spatial extent
      cfg%grid%ndim = 0
      cfg%grid%ny   = 1
      cfg%grid%nz   = 1
      cfg%grid%dy   = 0.0_dp
      cfg%grid%dz   = 0.0_dp
      ! No coordinate arrays allocated for bulk

    case (1)
      ! QW: 1D confinement along z
      cfg%grid%ndim = 1
      cfg%grid%ny   = 1
      cfg%grid%nz   = cfg%fdStep
      cfg%grid%dy   = 0.0_dp
      cfg%grid%dz   = cfg%dz

      ! Copy z-coordinate array from legacy field
      if (allocated(cfg%z)) then
        cfg%grid%z = cfg%z
      end if

      ! Build material_id from intStartPos/intEndPos
      if (allocated(cfg%intStartPos) .and. allocated(cfg%intEndPos)) then
        allocate(cfg%grid%material_id(cfg%fdStep))
        cfg%grid%material_id = 0
        do i = 1, cfg%numLayers
          cfg%grid%material_id(cfg%intStartPos(i):cfg%intEndPos(i)) = i
        end do
      end if

      ! No cut-cell fields for QW
      ! No coords(:,:) or y(:) for QW (ny=1)

    case (2)
      ! Wire: 2D confinement in y-z plane.
      ! Set grid dimensions only; coordinate arrays, material_id,
      ! cut-cell fields, and ghost_map are populated by
      ! init_wire_from_config() in the geometry module.
      cfg%grid%ndim = 2
      cfg%grid%ny   = cfg%wire_ny
      cfg%grid%nz   = cfg%wire_nz
      cfg%grid%dy   = cfg%wire_dy
      cfg%grid%dz   = cfg%wire_dz

    end select

  end subroutine init_grid_from_config

  subroutine tick(t)
      integer, intent(OUT) :: t

      call system_clock(t)
  end subroutine tick

  ! returns time in seconds from now to time described by t
  real function tock(t)
      integer, intent(in) :: t
      integer :: now, clock_rate

      call system_clock(now,clock_rate)

      tock = real(now - t)/real(clock_rate)
  end function tock


end module
