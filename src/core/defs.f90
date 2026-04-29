module definitions

  implicit none
  private

  ! Kinds
  public :: sp, dp, qp, iknd

  ! Physical constants
  public :: pi_sp, pi_dp, pi_qp
  public :: a0, hbar, const, c, m0, e, e0
  public :: kB_eV, tolerance, renormalization, hbar2O2m0
  public :: NUM_CB_STATES, NUM_VB_STATES, g_free
  public :: SQR3, SQR2, SQR2o3, SQR3o2, RQS3, RQS2
  public :: IU, ZERO, UM

  ! Types
  public :: bp_scalar, bir_pikus_blocks
  public :: wavevector, paramStruct
  public :: doping_spec, sc_config, wire_geometry, region_spec
  public :: spatial_grid, strain_config
  public :: optical_transition, optics_config
  public :: exciton_config, scattering_config
  public :: simulation_config, config_validation_result, group

  ! Functions and subroutines
  public :: kronij, grid_ngrid, init_grid_from_config

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
  real(kind=dp), parameter :: SQR3 = sqrt(3.0_dp)
  real(kind=dp), parameter :: SQR2 = sqrt(2.0_dp)
  real(kind=dp), parameter :: SQR2o3 = sqrt(2.0_dp/3.0_dp)
  real(kind=dp), parameter :: RQS3 = 1.0_dp/SQR3
  real(kind=dp), parameter :: RQS2 = 1.0_dp/SQR2
  real(kind=dp), parameter :: SQR3o2 = sqrt(1.5_dp)   ! sqrt(3/2) = SQR3 * RQS2
  complex(kind=dp), parameter :: ZERO = cmplx(0.0_dp, 0.0_dp, kind=dp)
  complex(kind=dp), parameter :: UM = cmplx(1.0_dp, 0.0_dp, kind=dp)

  ! Scalar Bir-Pikus strain result for a single grid point.
  type :: bp_scalar
    real(kind=dp) :: delta_Ec    = 0.0_dp   ! CB: ac*Tr(eps)
    real(kind=dp) :: delta_EHH   = 0.0_dp   ! HH: -P_eps + Q_eps
    real(kind=dp) :: delta_ELH   = 0.0_dp   ! LH: -P_eps - Q_eps
    real(kind=dp) :: delta_ESO   = 0.0_dp   ! SO: -P_eps
    complex(kind=dp) :: R_eps    = cmplx(0.0_dp, 0.0_dp, kind=dp)
    complex(kind=dp) :: S_eps    = cmplx(0.0_dp, 0.0_dp, kind=dp)
    real(kind=dp) :: QT2_eps     = 0.0_dp   ! Q_eps - T_eps = 2*Q_eps
  end type bp_scalar

  ! Per-grid-point Bir-Pikus strain Hamiltonian components.
  ! Arrays have length ngrid; unallocated when unstrained.
  type :: bir_pikus_blocks
    real(kind=dp), allocatable :: delta_Ec(:)    ! CB bands 7,8
    real(kind=dp), allocatable :: delta_EHH(:)   ! HH bands 1,4
    real(kind=dp), allocatable :: delta_ELH(:)   ! LH bands 2,3
    real(kind=dp), allocatable :: delta_ESO(:)   ! SO bands 5,6
    complex(kind=dp), allocatable :: R_eps(:)    ! VB-VB coupling
    complex(kind=dp), allocatable :: S_eps(:)    ! VB-VB coupling
    real(kind=dp), allocatable :: QT2_eps(:)     ! VB-SO coupling (2*Q_eps)
  contains
    final :: bir_pikus_blocks_finalize
  end type bir_pikus_blocks

  type wavevector

    real(kind=dp) :: kx, ky, kz

  end type wavevector

  type paramStruct

    real(kind=dp) :: meff, gamma1, gamma2, gamma3, P, A, deltaSO, EP, Eg, EV, EC
    real(kind=dp) :: eps0  ! static dielectric constant (unitless)

    ! Elastic constants and strain parameters (Vurgaftman 2001, Winkler 2003)
    ! av uses positive sign convention: P_eps = -av*Tr(eps)
    ! (opposite sign to Vurgaftman Table XIII)
    real(kind=dp) :: C11   = 0.0_dp    ! elastic constant (GPa)
    real(kind=dp) :: C12   = 0.0_dp    ! elastic constant (GPa)
    real(kind=dp) :: C44   = 0.0_dp    ! elastic constant (GPa)
    real(kind=dp) :: a0    = 0.0_dp    ! lattice constant (Angstrom)
    real(kind=dp) :: ac    = 0.0_dp    ! CB hydrostatic deformation potential (eV)
    real(kind=dp) :: av    = 0.0_dp    ! VB hydrostatic deformation potential (eV)
    real(kind=dp) :: b_dp  = 0.0_dp    ! shear deformation potential, tetragonal (eV)
    real(kind=dp) :: d_dp  = 0.0_dp    ! shear deformation potential, rhombohedral (eV)
    real(kind=dp) :: strainSubstrate = 0.0_dp  ! substrate lattice constant for bulk biaxial strain (Angstrom)

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
  ! For rectangle: wire_width (x), wire_height (y).
  ! For polygon: wire_polygon_nverts + wire_polygon_verts(:,:).
  ! ------------------------------------------------------------------
  type wire_geometry
    character(len=16)  :: shape = 'rectangle'
    real(kind=dp)      :: radius = 0.0_dp     ! AA (circle/hexagon)
    real(kind=dp)      :: width  = 0.0_dp     ! AA (rectangle, x-extent)
    real(kind=dp)      :: height = 0.0_dp     ! AA (rectangle, y-extent)
    integer            :: nverts = 0           ! polygon vertex count
    real(kind=dp), allocatable :: verts(:,:)  ! (2, nverts) x,y vertices
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
  ! Bulk:   nx=1, ny=1, dx=dy=0, no coordinate arrays.
  ! QW:     nx=1, ny=fdStep, dy>0, z(:) holds 1D grid.
  ! Wire:   nx>1, ny>1, dx>0, dy>0, coords(:,:) holds 2D grid.
  !
  ! For wire mode (confinement=2): nx and ny are the grid points in
  ! the x and y confinement directions respectively. The free
  ! propagation direction is z (kz sweep).
  !
  ! Cut-cell fields are only allocated for wire mode (ndim=2).  For QW
  ! and bulk they remain unallocated; code can test associated() or
  ! treat them as unity when unallocated.
  ! ------------------------------------------------------------------
  type spatial_grid
    integer :: ndim = 0                ! 0=bulk, 1=QW, 2=wire
    integer :: nx = 1                  ! grid points in x (confinement, wire mode)
    integer :: ny = 1                  ! grid points in y (confinement for wire, z-confinement for QW)
    real(kind=dp) :: dx = 0.0_dp       ! grid spacing in x (AA)
    real(kind=dp) :: dy = 0.0_dp       ! grid spacing in y (AA)

    ! 1D coordinate arrays (allocated for ndim >= 1)
    real(kind=dp), allocatable :: x(:)       ! (nx) x-coordinates for wire
    real(kind=dp), allocatable :: z(:)       ! (ny) y-coords for wire, z-coords for QW

    ! Flattened 2D coordinate array (allocated for ndim == 2)
    ! coords(1,:) = x, coords(2,:) = y (stored in z(:) array)
    ! column-major: flat_idx = (iy-1)*nx + ix
    real(kind=dp), allocatable :: coords(:,:)    ! (2, nx*ny)

    ! Material index at each grid point (1-based layer index)
    integer, allocatable  :: material_id(:)      ! (nx*ny)

    ! Cut-cell immersed boundary fields (ndim == 2 only)
    real(kind=dp), allocatable :: cell_volume(:)       ! (nx*ny) fractional vol [0,1]
    real(kind=dp), allocatable :: face_fraction_x(:,:) ! (nx*ny, 2) left/right face in x-direction
    real(kind=dp), allocatable :: face_fraction_y(:,:) ! (nx*ny, 2) bottom/top face in y-direction
    ! Nearest active neighbor for each ghost/inactive point
    integer, allocatable  :: ghost_map(:,:)            ! (nx*ny, 4) N S W E
  end type spatial_grid

  type strain_config
    logical          :: enabled      = .false.
    character(len=20) :: reference   = 'substrate'
    character(len=20) :: solver      = 'pardiso'
    logical          :: piezoelectric = .false.
  end type strain_config

  type :: optical_transition
    integer :: cb_idx, vb_idx     ! band indices
    real(kind=dp) :: energy       ! transition energy (eV)
    real(kind=dp) :: px, py, pz   ! |<i|dH/dk_x|j>|^2, etc. (dH/dk units)
    real(kind=dp) :: oscillator_strength  ! dimensionless, = sum|p|^2 / (hbar2O2m0 * dE)
  end type

  ! ------------------------------------------------------------------
  ! Optical spectra parameters (Tier 1: absorption, gain, ISBT).
  ! Default values correspond to a typical GaAs-based QW at 300 K.
  ! ------------------------------------------------------------------
  type optics_config
    logical          :: enabled = .false.
    real(kind=dp)    :: linewidth_lorentzian = 0.030_dp   ! eV FWHM
    real(kind=dp)    :: linewidth_gaussian = 0.005_dp     ! eV FWHM
    real(kind=dp)    :: refractive_index = 3.3_dp
    real(kind=dp)    :: temperature = 300.0_dp            ! K
    integer          :: num_energy_points = 200
    real(kind=dp)    :: E_min = 0.5_dp                    ! eV
    real(kind=dp)    :: E_max = 2.0_dp                    ! eV
    real(kind=dp)    :: carrier_density = 0.0_dp          ! 2D cm^-2 (0=equilibrium)
    ! Gain
    logical          :: gain_enabled = .false.
    real(kind=dp)    :: gain_carrier_density = 3.0e12_dp  ! 2D cm^-2
    ! ISBT
    logical          :: isbt_enabled = .false.
    ! Spontaneous emission
    logical          :: spontaneous_enabled = .false.
    ! Spin resolution
    logical          :: spin_resolved = .false.
  end type optics_config

  ! ------------------------------------------------------------------
  ! Exciton solver parameters (Tier 2).
  ! ------------------------------------------------------------------
  type exciton_config
    logical          :: enabled = .false.
    character(len=20) :: method = 'variational'
  end type exciton_config

  ! ------------------------------------------------------------------
  ! Phonon scattering parameters (Tier 3).
  ! ------------------------------------------------------------------
  type scattering_config
    logical          :: enabled = .false.
    real(kind=dp)    :: phonon_energy = 0.036_dp   ! eV (LO phonon)
    real(kind=dp)    :: eps_inf = 10.9_dp          ! high-freq dielectric
    real(kind=dp)    :: eps_0 = 12.9_dp            ! static dielectric
  end type scattering_config

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
    character(len=4) :: waveVector = 'k0'
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
    type(strain_config)            :: strain       ! strain solver parameters

    ! ---- Wire-specific fields (confinement=2) ----
    integer            :: wire_nx = 0            ! grid points in x
    integer            :: wire_ny = 0            ! grid points in y
    real(kind=dp)      :: wire_dx = 0.0_dp      ! grid spacing x (AA)
    real(kind=dp)      :: wire_dy = 0.0_dp      ! grid spacing y (AA)
    type(wire_geometry) :: wire_geom             ! shape + dimensions
    integer            :: numRegions = 0         ! number of material regions
    type(region_spec), allocatable :: regions(:) ! region specifications

    ! ---- FEAST eigensolver tuning ----
    real(kind=dp)      :: feast_emin = 0.0_dp   ! manual energy window (0=auto)
    real(kind=dp)      :: feast_emax = 0.0_dp   ! manual energy window (0=auto)
    integer            :: feast_m0 = 0           ! manual subspace size (0=auto: 2*nev)

    type(bir_pikus_blocks)     :: strain_blocks    ! precomputed strain data

    ! ---- Optical spectra / gain / exciton / scattering ----
    type(optics_config)      :: optics       ! optical spectra parameters
    type(exciton_config)     :: exciton      ! exciton solver parameters
    type(scattering_config)  :: scattering   ! phonon scattering parameters
  end type simulation_config

  type :: config_validation_result
    logical           :: ok      = .true.
    character(len=256) :: message = ''
  end type config_validation_result

  type group
      integer :: order    ! original order of unsorted data
      real(kind=dp) :: value       ! values to be sorted by
  end type group


  contains

  ! ==================================================================
  ! Finalizer: automatically called when a bir_pikus_blocks goes out of scope.
  ! Deallocates all allocatable components.
  ! ==================================================================
  subroutine bir_pikus_blocks_finalize(bp)
    type(bir_pikus_blocks), intent(inout) :: bp

    if (allocated(bp%delta_Ec))  deallocate(bp%delta_Ec)
    if (allocated(bp%delta_EHH)) deallocate(bp%delta_EHH)
    if (allocated(bp%delta_ELH)) deallocate(bp%delta_ELH)
    if (allocated(bp%delta_ESO)) deallocate(bp%delta_ESO)
    if (allocated(bp%R_eps))     deallocate(bp%R_eps)
    if (allocated(bp%S_eps))     deallocate(bp%S_eps)
    if (allocated(bp%QT2_eps))   deallocate(bp%QT2_eps)
  end subroutine bir_pikus_blocks_finalize

  elemental pure function kronij(i,j)
    integer, intent(in) :: i,j
    integer :: kronij
    kronij = merge(1, 0, i == j)
  end function

  ! ------------------------------------------------------------------
  ! Return total number of spatial grid points (nx * ny).
  ! Works for all ndim: bulk returns 1, QW returns ny, wire returns nx*ny.
  ! ------------------------------------------------------------------
  pure function grid_ngrid(grid) result(n)
    type(spatial_grid), intent(in) :: grid
    integer :: n
    n = grid%nx * grid%ny
  end function

  ! ------------------------------------------------------------------
  ! Initialize a spatial_grid from the fields already set in
  ! simulation_config.  This is the primary grid initialization path.
  !
  ! For bulk (confinement=0): ndim=0, nx=ny=1.
  ! For QW  (confinement=1): ndim=1, nx=1, ny=fdStep, z(:) copied.
  ! For wire(confinement=2): ndim=2, nx=wire_nx, ny=wire_ny,
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
      cfg%grid%nx   = 1
      cfg%grid%ny   = 1
      cfg%grid%dx   = 0.0_dp
      cfg%grid%dy   = 0.0_dp
      ! No coordinate arrays allocated for bulk

    case (1)
      ! QW: 1D confinement along z
      cfg%grid%ndim = 1
      cfg%grid%nx   = 1
      cfg%grid%ny   = cfg%fdStep
      cfg%grid%dx   = 0.0_dp
      cfg%grid%dy   = cfg%dz

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
      ! No coords(:,:) or x(:) for QW (nx=1)

    case (2)
      ! Wire: 2D confinement in x-y plane.
      ! Set grid dimensions only; coordinate arrays, material_id,
      ! cut-cell fields, and ghost_map are populated by
      ! init_wire_from_config() in the geometry module.
      cfg%grid%ndim = 2
      cfg%grid%nx   = cfg%wire_nx
      cfg%grid%ny   = cfg%wire_ny
      cfg%grid%dx   = cfg%wire_dx
      cfg%grid%dy   = cfg%wire_dy

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
