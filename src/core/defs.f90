module definitions

  use, intrinsic :: iso_fortran_env, only: real32, real64, real128, int32

  implicit none
  private

  ! Kinds
  public :: sp, dp, qp, iknd

  ! Physical constants
  public :: pi_sp, pi_dp, pi_qp
  public :: a0, hbar, const, c, m0, e, e0
  public :: kB_eV, mu_B, tolerance, coord_tolerance, renormalization, hbar2O2m0
  public :: NUM_CB_STATES, NUM_VB_STATES, g_free
  public :: SQR3, SQR2, SQR2o3, SQR3o2, RQS3, RQS2
  public :: IU, ZERO, UM

  ! Types
  public :: bp_scalar, bir_pikus_blocks
  public :: wavevector, paramStruct
  public :: doping_spec, sc_config, wire_geometry, region_spec
  public :: spatial_grid, strain_config
  public :: optical_transition, optics_config
  public :: exciton_config, scattering_config, bdg_config, topology_config, topological_result
  public :: wave_vector_config, bands_config, external_field_config
  public :: b_field_config, solver_config, landau_config, wire_config
  public :: simulation_config, group

  ! Functions and subroutines
  public :: conf_direction
  public :: resolve_solver_defaults
  public :: to_upper_trim
  public :: kronij, grid_ngrid, init_grid_from_config
  public :: validate_semantic

  integer, parameter :: sp   = real32
  integer, parameter :: dp   = real64
  integer, parameter :: qp   = real128
  integer, parameter :: iknd = int32

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
  real(kind=dp), parameter :: mu_B = 5.7883818012e-5_dp       ! Bohr magneton (eV/T)
  real(kind=dp), parameter :: tolerance=1e-7_dp
  real(kind=dp), parameter :: coord_tolerance=1e-12_dp  ! spatial coordinate zero-check threshold
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
  ! No FINAL binding: gfortran creates temporaries for intent(in) derived-type
  ! arguments (copy-in/copy-out), and the finalizer would deallocate the
  ! allocatable components of the temporary, corrupting the original via
  ! shared allocation handles. Cleanup is done explicitly via
  ! bir_pikus_blocks_free() and in simulation_config_finalize.
  type :: bir_pikus_blocks
    real(kind=dp), allocatable :: delta_Ec(:)    ! CB bands 7,8
    real(kind=dp), allocatable :: delta_EHH(:)   ! HH bands 1,4
    real(kind=dp), allocatable :: delta_ELH(:)   ! LH bands 2,3
    real(kind=dp), allocatable :: delta_ESO(:)   ! SO bands 5,6
    complex(kind=dp), allocatable :: R_eps(:)    ! VB-VB coupling
    complex(kind=dp), allocatable :: S_eps(:)    ! VB-VB coupling
    real(kind=dp), allocatable :: QT2_eps(:)     ! VB-SO coupling (2*Q_eps)
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
    character(len=7) :: dtype = 'uniform'  ! 'uniform' or 'delta'
    real(kind=dp) :: NS = 0.0_dp         ! delta: 2D sheet density (10^11 cm^-2)
    real(kind=dp) :: delta_pos = 0.0_dp   ! delta: position along z (Angstrom)
    real(kind=dp) :: delta_fwhm = 10.0_dp ! delta: Gaussian FWHM (Angstrom)
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
  ! Wire geometry specification for confinement='wire'.
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
  contains
    final :: wire_geometry_finalize
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
  ! For wire mode (confinement='wire'): nx and ny are the grid points in
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
  contains
    procedure :: npoints => spatial_grid_npoints
    final :: spatial_grid_finalize
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
    ! Geometry context (set from cfg%confinement before k-sweep)
    character(len=8) :: confinement = 'bulk'
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

  ! ------------------------------------------------------------------
  ! BdG / topological superconductivity parameters.
  ! ------------------------------------------------------------------
  type :: bdg_config
    logical          :: enabled = .false.
    real(kind=dp)   :: mu = 0.0_dp       ! chemical potential (eV)
    real(kind=dp)   :: delta_0 = 0.0_dp  ! s-wave pairing gap (eV)
    real(kind=dp)   :: B_vec(3) = 0.0_dp ! magnetic field Bx, By, Bz (Tesla)
    real(kind=dp)   :: g_factor = 2.0_dp  ! Lande g-factor for Zeeman splitting
    character(len=20) :: gauge = 'landau_x'  ! gauge choice
    real(kind=dp)   :: B_sweep(3) = 0.0_dp ! B sweep: min, max, step
    real(kind=dp)   :: kz = 0.0_dp        ! along-wire momentum (1/A) for BdG dispersion
  end type

  ! ------------------------------------------------------------------
  ! Unified topological analysis parameters (QHE / QSHE / BdG modes).
  ! ------------------------------------------------------------------
  type :: topology_config
    logical          :: enabled = .false.
    character(len=20) :: mode = 'qhe'     ! qhe | qshe | bdg
    ! Chern / Berry curvature
    logical          :: compute_chern = .false.
    logical          :: compute_hall = .false.   ! output sigma_xy = C*e^2/h
    real(kind=dp)    :: qwz_u = 0.0_dp           ! QWZ mass parameter u
    ! Z2 invariant
    logical          :: compute_z2 = .false.
    character(len=20) :: z2_method = 'auto'  ! auto | gap | fukane
    ! BHZ wire parameters for QSHE mode
    real(kind=dp)    :: bhz_M = 10.0_dp      ! BHZ mass parameter (meV)
    real(kind=dp)    :: bhz_d = 58.0_dp      ! wire width (AA)
    ! Edge states
    logical          :: extract_edge_states = .false.
    real(kind=dp)    :: edge_E_window = 0.01_dp  ! energy window for detection
    ! LDOS
    logical          :: compute_ldos = .false.
    real(kind=dp)    :: ldos_eta = 0.001_dp  ! Lorentzian broadening
    real(kind=dp)    :: ldos_E_range(2) = [-0.1_dp, 0.1_dp]
    integer          :: ldos_num_E = 200
    ! Gap sweep / phase diagram
    logical          :: compute_gap_sweep = .false.
    real(kind=dp)    :: gap_sweep_B_min = 0.0_dp
    real(kind=dp)    :: gap_sweep_B_max = 1.0_dp
    integer          :: gap_sweep_nB = 20
    real(kind=dp)    :: gap_sweep_mu_min = 0.0_dp
    real(kind=dp)    :: gap_sweep_mu_max = 0.01_dp
    integer          :: gap_sweep_nMu = 20
    character(len=20) :: sweep_model = 'bhz_analytic'  ! bhz_analytic | wire_bdg | qw_fukane
    ! Conductance
    logical          :: compute_conductance = .false.
    character(len=20) :: conductance_method = 'kubo_chern'  ! kubo_chern | kubo_berry | landauer
    integer          :: berry_nk = 50
    real(kind=dp)    :: landauer_energy = 0.0_dp
    ! Spectral function
    logical          :: compute_spectral = .false.
    real(kind=dp)    :: spectral_k_min = -0.1_dp
    real(kind=dp)    :: spectral_k_max = 0.1_dp
    integer          :: spectral_nk = 100
    real(kind=dp)    :: spectral_E_min = -0.05_dp
    real(kind=dp)    :: spectral_E_max = 0.05_dp
    integer          :: spectral_nE = 200
    real(kind=dp)    :: spectral_eta = 0.001_dp
  end type

  ! ------------------------------------------------------------------
  ! Results from topological analysis.
  ! ------------------------------------------------------------------
  type :: topological_result
    integer          :: chern_number = 0
    integer          :: z2_invariant = 0
    real(kind=dp)    :: hall_conductance = 0.0_dp   ! in units of e^2/h
    real(kind=dp)    :: min_gap = 0.0_dp
    real(kind=dp)    :: edge_xi_min = 0.0_dp      ! min edge localization length
    real(kind=dp)    :: edge_xi = 0.0_dp           ! average edge localization length
    integer          :: n_majorana = 0
    integer          :: n_majorana_fit_failed = 0
    real(kind=dp), allocatable :: edge_energies(:)
    real(kind=dp), allocatable :: phase_boundary(:,:)  ! (B, mu) pairs
    real(kind=dp), allocatable :: berry_curvature(:,:) ! Omega(kx, ky) if computed
    real(kind=dp)              :: conductance_xy = 0.0_dp  ! Hall conductance (e^2/h)
    real(kind=dp)              :: conductance_zz = 0.0_dp  ! longitudinal (2e^2/h)
    real(kind=dp), allocatable :: spectral_function(:,:)   ! A(k, E) heatmap
    real(kind=dp), allocatable :: z2_map(:,:)              ! Z2 phase diagram (nMu x nB)
    real(kind=dp), allocatable :: gap_map(:,:)             ! gap phase diagram (nMu x nB)
  contains
    final :: topological_result_finalize
  end type

  ! ---- New sub-types mirroring TOML sections ----

  ! Replaces flat waveVector/waveVectorMax/waveVectorStep
  type :: wave_vector_config
    character(len=8) :: mode = 'k0'
    real(kind=dp)    :: max = 0.0_dp
    real(kind=dp)    :: step = 0.01_dp
    integer          :: nsteps = 100       ! number of k-points in sweep
  end type

  ! Replaces flat numcb/numvb
  type :: bands_config
    integer :: num_cb = 2
    integer :: num_vb = 6
  end type

  ! Replaces flat ExternalField/EFtype/Evalue
  type :: external_field_config
    logical          :: enabled = .false.   ! true if [external_field] section present
    character(len=2) :: type = 'EF'        ! "EF" for electric field
    real(kind=dp)    :: value = 0.0_dp
  end type

  ! Replaces flat b_field(3) and bdg%g_factor for bulk Zeeman
  type :: b_field_config
    real(kind=dp) :: components(3) = 0.0_dp  ! Bx, By, Bz (Tesla)
    real(kind=dp) :: g_factor = 2.0_dp
  end type

  ! Eigensolver configuration (parsed from [solver] TOML section)
  type :: solver_config
    character(len=10) :: method = 'AUTO'   ! AUTO, DENSE, FEAST
    character(len=10) :: mode   = 'AUTO'   ! AUTO, FULL, INDEX, ENERGY
    real(kind=dp)     :: emin   = 0.0_dp   ! 0 = auto
    real(kind=dp)     :: emax   = 0.0_dp   ! 0 = auto
    integer           :: m0     = 0         ! 0 = auto (FEAST subspace)
  end type

  ! Replaces scattered landau_* fields
  type :: landau_config
    integer            :: nx = 100
    real(kind=dp)      :: width = 2000.0_dp
    character(len=4)   :: sweep = 'ky'
    character(len=255) :: material = ''
  end type

  ! Replaces scattered wire_* fields
  type :: wire_config
    integer            :: nx = 0
    integer            :: ny = 0
    real(kind=dp)      :: dx = 0.0_dp
    real(kind=dp)      :: dy = 0.0_dp
    type(wire_geometry) :: geom
    type(region_spec), allocatable :: regions(:)
    integer            :: num_regions = 0
  contains
    final :: wire_config_finalize
  end type

  type simulation_config
    ! ---- User inputs (from TOML) ----
    character(len=8)   :: confinement = 'bulk'   ! "bulk"|"qw"|"wire"|"landau"
    integer            :: FDorder = 2
    integer            :: fd_step = 1             ! user-specified grid points (QW)
    type(wave_vector_config)   :: wave_vector
    type(bands_config)         :: bands
    type(external_field_config) :: external_field
    type(b_field_config)       :: b_field
    type(strain_config)        :: strain
    type(solver_config)        :: solver

    ! Mode-specific geometry
    type(wire_config)          :: wire              ! confinement = "wire"
    type(landau_config)        :: landau            ! confinement = "landau"

    ! Material layers (QW: num_layers = nmat; bulk/landau/wire: num_layers = 1)
    ! Wire regions tracked separately in cfg%wire%num_regions.
    ! material_names/params arrays are sized by material count (QW) or region
    ! count (wire), which may differ from num_layers.
    character(len=255), allocatable :: material_names(:)
    real(dp), allocatable :: z_min(:), z_max(:)
    integer               :: num_layers = 1

    ! Physics subsystems (section presence = enabled)
    type(sc_config)           :: sc
    type(doping_spec), allocatable :: doping(:)
    type(topology_config)     :: topo
    type(bdg_config)          :: bdg
    type(optics_config)       :: optics
    type(exciton_config)      :: exciton
    type(scattering_config)   :: scattering

    ! ---- Computed fields (derived from inputs) ----
    integer            :: evnum = 8            ! bands%num_cb + bands%num_vb
    real(kind=dp)      :: totalSize = 0.0_dp
    real(kind=dp)      :: delta = 0.0_dp
    real(kind=dp)      :: dz = 0.0_dp
    real(kind=dp), allocatable :: z(:)
    integer, allocatable   :: int_start_pos(:), int_end_pos(:)
    type(paramStruct), allocatable :: params(:)
    type(bir_pikus_blocks) :: strain_blocks
    type(spatial_grid)     :: grid
    real(kind=dp)          :: sc_potential_shift = 0.0_dp

    ! G-factor
    integer :: which_band = 0    ! 0=CB, 1=VB
    integer :: band_idx = 1      ! which doublet

  contains
    final :: simulation_config_finalize
    procedure :: validate => simulation_config_validate
  end type simulation_config

  type group
      integer :: order    ! original order of unsorted data
      real(kind=dp) :: value       ! values to be sorted by
  end type group


  contains

  ! ==================================================================

  ! ==================================================================
  ! Finalizer: automatically called when a wire_geometry goes out of scope.
  ! Deallocates the polygon vertex array.
  ! ==================================================================
  subroutine wire_geometry_finalize(wg)
    type(wire_geometry), intent(inout) :: wg

    if (allocated(wg%verts)) deallocate(wg%verts)
  end subroutine wire_geometry_finalize

  ! ==================================================================
  ! Finalizer: automatically called when a spatial_grid goes out of scope.
  ! Deallocates all allocatable components.
  ! ==================================================================
  subroutine spatial_grid_finalize(sg)
    type(spatial_grid), intent(inout) :: sg

    if (allocated(sg%x))                deallocate(sg%x)
    if (allocated(sg%z))                deallocate(sg%z)
    if (allocated(sg%coords))           deallocate(sg%coords)
    if (allocated(sg%material_id))      deallocate(sg%material_id)
    if (allocated(sg%cell_volume))      deallocate(sg%cell_volume)
    if (allocated(sg%face_fraction_x))  deallocate(sg%face_fraction_x)
    if (allocated(sg%face_fraction_y))  deallocate(sg%face_fraction_y)
    if (allocated(sg%ghost_map))        deallocate(sg%ghost_map)
  end subroutine spatial_grid_finalize

  ! ==================================================================
  ! Finalizer: automatically called when a wire_config goes out of scope.
  ! ==================================================================
  subroutine wire_config_finalize(wc)
    type(wire_config), intent(inout) :: wc

    if (allocated(wc%regions)) deallocate(wc%regions)
    ! wire_geom has its own finalizer
  end subroutine wire_config_finalize

  ! ==================================================================
  ! Finalizer: automatically called when a simulation_config goes out
  ! of scope.  Deallocates all allocatable components.
  ! Note: cfg%grid, cfg%wire%geom have their own finalizers which
  ! will fire automatically when this type is destroyed.
  ! ==================================================================
  subroutine simulation_config_finalize(cfg)
    type(simulation_config), intent(inout) :: cfg

    if (allocated(cfg%z_min))          deallocate(cfg%z_min)
    if (allocated(cfg%z_max))          deallocate(cfg%z_max)
    if (allocated(cfg%z))              deallocate(cfg%z)
    if (allocated(cfg%int_start_pos))  deallocate(cfg%int_start_pos)
    if (allocated(cfg%int_end_pos))    deallocate(cfg%int_end_pos)
    if (allocated(cfg%material_names)) deallocate(cfg%material_names)
    if (allocated(cfg%params))         deallocate(cfg%params)
    if (allocated(cfg%doping))         deallocate(cfg%doping)
    if (allocated(cfg%strain_blocks%delta_Ec))  deallocate(cfg%strain_blocks%delta_Ec)
    if (allocated(cfg%strain_blocks%delta_EHH)) deallocate(cfg%strain_blocks%delta_EHH)
    if (allocated(cfg%strain_blocks%delta_ELH)) deallocate(cfg%strain_blocks%delta_ELH)
    if (allocated(cfg%strain_blocks%delta_ESO)) deallocate(cfg%strain_blocks%delta_ESO)
    if (allocated(cfg%strain_blocks%R_eps))     deallocate(cfg%strain_blocks%R_eps)
    if (allocated(cfg%strain_blocks%S_eps))     deallocate(cfg%strain_blocks%S_eps)
    if (allocated(cfg%strain_blocks%QT2_eps))   deallocate(cfg%strain_blocks%QT2_eps)
  end subroutine simulation_config_finalize

  ! ==================================================================
  ! Finalizer: automatically called when a topological_result goes
  ! out of scope.  Deallocates all allocatable components.
  ! ==================================================================
  subroutine topological_result_finalize(res)
    type(topological_result), intent(inout) :: res

    if (allocated(res%edge_energies))    deallocate(res%edge_energies)
    if (allocated(res%phase_boundary))   deallocate(res%phase_boundary)
    if (allocated(res%berry_curvature))  deallocate(res%berry_curvature)
    if (allocated(res%spectral_function)) deallocate(res%spectral_function)
    if (allocated(res%z2_map)) deallocate(res%z2_map)
    if (allocated(res%gap_map)) deallocate(res%gap_map)
  end subroutine topological_result_finalize

  ! ==================================================================
  ! Type-bound validation for simulation_config.
  ! Checks invariants already assumed by downstream code.
  ! Calls error stop on failure (F2008).
  ! ==================================================================
  subroutine simulation_config_validate(self)
    class(simulation_config), intent(in) :: self

    associate(cfg => self)
      ! grid: must have at least one point
      if (cfg%grid%npoints() <= 0) then
        error stop 'validate_simulation_config: grid must have at least one point'
      end if

      ! bands: must be positive (at least 1 CB and 1 VB)
      if (cfg%bands%num_cb < 1) then
        error stop 'validate_simulation_config: num_cb must be >= 1'
      end if
      if (cfg%bands%num_vb < 1) then
        error stop 'validate_simulation_config: num_vb must be >= 1'
      end if

      ! evnum must equal num_cb + num_vb (set by parser)
      if (cfg%evnum /= cfg%bands%num_cb + cfg%bands%num_vb) then
        error stop 'validate_simulation_config: evnum must equal num_cb + num_vb'
      end if

      ! num_layers: must be positive
      if (cfg%num_layers < 1) then
        error stop 'validate_simulation_config: num_layers must be >= 1'
      end if

      ! confinement: must be one of the valid strings
      if (trim(cfg%confinement) /= 'bulk' .and. trim(cfg%confinement) /= 'qw' &
          .and. trim(cfg%confinement) /= 'wire' .and. trim(cfg%confinement) /= 'landau') then
        error stop 'validate_simulation_config: confinement must be bulk, qw, wire, or landau'
      end if
      if (trim(cfg%confinement) == 'qw') then
        if (cfg%fd_step < 3) then
          error stop 'validate_simulation_config: QW fd_step must be >= 3'
        end if
        if (cfg%fd_step < cfg%FDorder + 1) then
          error stop 'validate_simulation_config: QW fd_step must be >= FDorder + 1'
        end if
      end if
      if (trim(cfg%confinement) == 'wire') then
        if (cfg%wire%nx < cfg%FDorder + 1) then
          error stop 'validate_simulation_config: wire nx must be >= FDorder + 1'
        end if
        if (cfg%wire%ny < cfg%FDorder + 1) then
          error stop 'validate_simulation_config: wire ny must be >= FDorder + 1'
        end if
        if (cfg%wire%num_regions < 1) then
          error stop 'validate_simulation_config: wire num_regions must be >= 1'
        end if
      end if
      if (trim(cfg%confinement) == 'landau') then
        if (cfg%landau%nx < 3) then
          error stop 'validate_simulation_config: landau nx must be >= 3'
        end if
        if (cfg%landau%nx < cfg%FDorder + 1) then
          error stop 'validate_simulation_config: landau nx must be >= FDorder + 1'
        end if
        if (cfg%landau%width <= 0.0_dp) then
          error stop 'validate_simulation_config: landau width must be > 0'
        end if
      end if

      ! FDorder: must be one of the supported values
      if (cfg%FDorder /= 2 .and. cfg%FDorder /= 4 .and. cfg%FDorder /= 6 &
          .and. cfg%FDorder /= 8 .and. cfg%FDorder /= 10) then
        error stop 'validate_simulation_config: FDorder must be 2, 4, 6, 8, or 10'
      end if

      ! material_names must be allocated with at least num_layers entries
      ! (wire mode: sized by wire%num_regions, which may exceed num_layers)
      if (.not. allocated(cfg%material_names)) then
        error stop 'validate_simulation_config: material_names not allocated'
      end if
      if (size(cfg%material_names) < cfg%num_layers) then
        error stop 'validate_simulation_config: material_names too small for num_layers'
      end if

      ! params must be allocated with at least num_layers entries
      ! (wire mode: sized by wire%num_regions, which may exceed num_layers)
      if (.not. allocated(cfg%params)) then
        error stop 'validate_simulation_config: params not allocated'
      end if
      if (size(cfg%params) < cfg%num_layers) then
        error stop 'validate_simulation_config: params too small for num_layers'
      end if

      ! ---- V1: bulk is fully diagonalized (8x8 Hamiltonian), so the
      ! requested band window must equal 8. For bulk it is a display filter
      ! only (see CONTEXT.md "requested band window"), not a compute
      ! directive — sizing storage by a window < 8 would write out of bounds
      ! because the bulk solve always returns all 8 eigenpairs.
      if (trim(cfg%confinement) == 'bulk') then
        if (cfg%evnum /= 8) then
          block
            character(len=16) :: buf
            write(buf, '(I0)') cfg%evnum
            error stop 'validate_simulation_config: bulk band count (evnum=' // &
              trim(buf) // ') must equal 8 (8x8 Hamiltonian is fully diagonalized)'
          end block
        end if
      end if

      ! ---- V2: QW num_cb must not exceed NUM_CB_STATES * npoints ----
      if (trim(cfg%confinement) == 'qw') then
        if (cfg%bands%num_cb > NUM_CB_STATES * cfg%grid%npoints()) then
          block
            character(len=16) :: buf_cb, buf_max
            write(buf_cb, '(I0)') cfg%bands%num_cb
            write(buf_max, '(I0)') NUM_CB_STATES * cfg%grid%npoints()
            error stop 'validate_simulation_config: num_cb (=' // trim(buf_cb) // &
              ') exceeds maximum (' // trim(buf_max) // ') for QW confinement'
          end block
        end if
      end if

      ! ---- V3: QW num_vb must not exceed NUM_VB_STATES * npoints ----
      if (trim(cfg%confinement) == 'qw') then
        if (cfg%bands%num_vb > NUM_VB_STATES * cfg%grid%npoints()) then
          block
            character(len=16) :: buf_vb, buf_max
            write(buf_vb, '(I0)') cfg%bands%num_vb
            write(buf_max, '(I0)') NUM_VB_STATES * cfg%grid%npoints()
            error stop 'validate_simulation_config: num_vb (=' // trim(buf_vb) // &
              ') exceeds maximum (' // trim(buf_max) // ') for QW confinement'
          end block
        end if
      end if

      ! ---- V4: wave_vector mode must be recognized ----
      select case (trim(cfg%wave_vector%mode))
      case ('kx', 'ky', 'kz', 'kxky', 'kxkz', 'kykz', 'k0')
        ! valid
      case default
        error stop 'validate_simulation_config: wave_vector mode ''' // &
          trim(cfg%wave_vector%mode) // ''' not recognized (expected one of: ' // &
          'kx, ky, kz, kxky, kxkz, kykz, k0)'
      end select

      ! ---- V5: B-sweep step must be positive when configured ----
      if (any(cfg%bdg%B_sweep /= 0.0_dp)) then
        if (cfg%bdg%B_sweep(3) <= 0.0_dp) then
          block
            character(len=32) :: buf
            write(buf, '(ES12.4)') cfg%bdg%B_sweep(3)
            error stop 'validate_simulation_config: B_sweep step must be positive, got ' // &
              trim(buf)
          end block
        end if
      end if

      ! ---- V6: z_min < z_max for all material layers ----
      if (allocated(cfg%z_min) .and. allocated(cfg%z_max)) then
        block
          integer :: il
          character(len=16) :: buf_i, buf_zmin, buf_zmax
          do il = 1, cfg%num_layers
            if (cfg%z_min(il) >= cfg%z_max(il)) then
              write(buf_i, '(I0)') il
              write(buf_zmin, '(ES12.4)') cfg%z_min(il)
              write(buf_zmax, '(ES12.4)') cfg%z_max(il)
              error stop 'validate_simulation_config: material layer ' // trim(buf_i) // &
                ': z_min (' // trim(buf_zmin) // ') >= z_max (' // trim(buf_zmax) // ')'
            end if
          end do
        end block
      end if

      ! ---- V7: SC + bulk confinement is invalid ----
      if (cfg%sc%enabled == 1 .and. trim(cfg%confinement) == 'bulk') then
        error stop 'validate_simulation_config: SC loop requires confinement=''qw'', got ''bulk'''
      end if

      ! ---- V8: electric field requires z(1) /= 0 ----
      if (cfg%external_field%enabled .and. allocated(cfg%z)) then
        if (cfg%external_field%type == 'EF') then
          if (abs(cfg%z(1)) < coord_tolerance) then
            error stop 'validate_simulation_config: electric field requires z(1) /= 0'
          end if
        end if
      end if

      ! ---- I9: which_band must be 0 or 1 ----
      if (cfg%which_band < 0 .or. cfg%which_band > 1) then
        block
          character(len=16) :: buf
          write(buf, '(I0)') cfg%which_band
          error stop 'validate_simulation_config: which_band must be 0 or 1, got ' // trim(buf)
        end block
      end if

      ! ---- I11: solver%emin < solver%emax when both nonzero ----
      if (cfg%solver%emin /= 0.0_dp .and. cfg%solver%emax /= 0.0_dp) then
        if (cfg%solver%emin >= cfg%solver%emax) then
          block
            character(len=32) :: buf_emin, buf_emax
            write(buf_emin, '(ES12.4)') cfg%solver%emin
            write(buf_emax, '(ES12.4)') cfg%solver%emax
            error stop 'validate_simulation_config: solver%emin (' // trim(buf_emin) // &
              ') must be < solver%emax (' // trim(buf_emax) // ')'
          end block
        end if
      end if

      ! ---- I12: solver%m0 must be in [1, 1000] when explicitly set (> 0) ----
      ! m0 <= 0 is the "auto-detect" sentinel (eigensolver replaces with 2*nev).
      if (cfg%solver%m0 > 0) then
        if (cfg%solver%m0 > 1000) then
          block
            character(len=16) :: buf
            write(buf, '(I0)') cfg%solver%m0
            error stop 'validate_simulation_config: solver%m0 must be in [1, 1000], got ' // trim(buf)
          end block
        end if
      end if

      ! ---- I13: solver%method must be AUTO, DENSE, or FEAST ----
      select case (trim(cfg%solver%method))
      case ('AUTO', 'DENSE', 'FEAST')
        ! ok
      case default
        error stop 'validate_simulation_config: solver%method must be AUTO, DENSE, or FEAST, got "' // &
          trim(cfg%solver%method) // '"'
      end select

      ! ---- I14: solver%mode must be AUTO, FULL, INDEX, or ENERGY ----
      select case (trim(cfg%solver%mode))
      case ('AUTO', 'FULL', 'INDEX', 'ENERGY')
        ! ok
      case default
        error stop 'validate_simulation_config: solver%mode must be AUTO, FULL, INDEX, or ENERGY, got "' // &
          trim(cfg%solver%mode) // '"'
      end select

      ! ---- I14b: a partial energy window is ambiguous ----
      ! 0 is the auto sentinel for BOTH emin and emax. Setting exactly one
      ! of them is ambiguous (is the 0 a literal bound or "auto"?), and the
      ! window authority's OR-override would otherwise produce a degenerate
      ! [emin, 0] window. Require both-or-neither. (Review finding #7.)
      if ((cfg%solver%emin == 0.0_dp) .neqv. (cfg%solver%emax == 0.0_dp)) then
        error stop 'validate_simulation_config: set both solver emin and emax, ' // &
          'or neither (0 = auto). A partial energy window is ambiguous.'
      end if

      ! ---- I15: FEAST method cannot combine with INDEX mode ----
      ! FEAST has no INDEX (range il:iu) interface; only ENERGY or FULL are
      ! supported. Caught here at input validation so the user sees a clear
      ! message before any eigensolver is constructed. The eigensolver.f90
      ! guards (eigensolver_config_validate + feast_solve_sparse_dispatch)
      ! remain as defense-in-depth for any code path that bypasses validate().
      if (trim(cfg%solver%method) == 'FEAST' .and. trim(cfg%solver%mode) == 'INDEX') then
        error stop 'validate_simulation_config: FEAST solver does not support INDEX mode ' // &
          '(use mode = ENERGY or FULL)'
      end if

    end associate

  end subroutine simulation_config_validate

  ! ==================================================================
  ! Semantic validation: app-specific constraints dispatched on app_name.
  ! Called by each main*.f90 executable after read_config().
  ! Calls error stop on failure (F2008).
  ! ==================================================================
  subroutine validate_semantic(cfg, app_name)
    type(simulation_config), intent(in) :: cfg
    character(len=*), intent(in)        :: app_name

    select case (trim(app_name))

    case ('bandStructure')
      ! No additional semantic constraints for bandStructure.

    case ('gfactor')
      if (cfg%wave_vector%nsteps /= 0 .and. cfg%wave_vector%mode /= 'k0') then
        error stop 'validate_semantic: gfactor requires k0 mode (wave_vector%nsteps=0 or mode=k0)'
      end if
      ! gfactor needs the FULL spectrum for Lowdin partitioning. For QW,
      ! setup_solve_gamma_point (simulation_setup.f90) routes the Gamma
      ! solve through the dispersion-aware window authority (apply_solver_window,
      ! issue #03) so a FEAST backend runs in ENERGY mode with a window that
      ! covers the ENTIRE spectral range (Gershgorin-envelope bound), making
      ! FEAST return all 8N eigenvalues and Lowdin see the full spectrum
      ! (issue #08, ADR 0005). A DENSE backend still uses FULL mode (zheev
      ! returns everything). FEAST+INDEX remains rejected up front at the
      ! structural check I15 above (FEAST has no INDEX interface).

      ! S1: bandIdx in range for gfactor (bulk, QW, wire)
      ! CB (which_band=0): accesses cb_state(:, bandIdx:bandIdx+1) → range [1, num_cb-1]
      ! VB (which_band=1): accesses vb_state(:, bandIdx:bandIdx+1) → range [1, num_vb-1]
      if (cfg%which_band == 0) then
        if (cfg%band_idx < 1 .or. cfg%band_idx + 1 > cfg%bands%num_cb) then
          block
            character(len=16) :: buf_idx, buf_max
            write(buf_idx, '(I0)') cfg%band_idx
            write(buf_max, '(I0)') cfg%bands%num_cb - 1
            error stop 'validate_semantic: bandIdx (=' // trim(buf_idx) // &
              ') out of range [1, ' // trim(buf_max) // '] for CB gfactor'
          end block
        end if
      else
        if (cfg%band_idx < 1 .or. cfg%band_idx + 1 > cfg%bands%num_vb) then
          block
            character(len=16) :: buf_idx, buf_max
            write(buf_idx, '(I0)') cfg%band_idx
            write(buf_max, '(I0)') cfg%bands%num_vb - 1
            error stop 'validate_semantic: bandIdx (=' // trim(buf_idx) // &
              ') out of range [1, ' // trim(buf_max) // '] for VB gfactor'
          end block
        end if
      end if

    case ('opticalProperties')
      if (.not. cfg%optics%enabled) then
        error stop 'validate_semantic: opticalProperties requires [optics] block with enabled = true'
      end if
      ! I8: optics does not support landau confinement
      if (trim(cfg%confinement) == 'landau') then
        error stop 'validate_semantic: opticalProperties does not support confinement=''landau''' // &
          ' (supported: bulk, qw, wire)'
      end if

    case ('topologicalAnalysis')
      if (.not. cfg%topo%enabled) then
        error stop 'validate_semantic: topologicalAnalysis requires [topology] block in input.toml'
      end if
      if (len_trim(cfg%topo%mode) == 0) then
        error stop 'validate_semantic: topologicalAnalysis requires topology mode to be set'
      end if

      ! S10: topology mode must be a recognized enum
      if (trim(cfg%topo%mode) /= 'qhe' .and. &
          trim(cfg%topo%mode) /= 'qshe' .and. &
          trim(cfg%topo%mode) /= 'bdg' .and. &
          trim(cfg%topo%mode) /= 'spectral' .and. &
          trim(cfg%topo%mode) /= 'conductance' .and. &
          trim(cfg%topo%mode) /= 'sweep') then
        error stop 'validate_semantic: topology mode ''' // trim(cfg%topo%mode) // &
          ''' not recognized (expected: qhe, qshe, bdg, spectral, conductance, sweep)'
      end if

      ! S2: QSHE Z2 requires QW or wire confinement
      if (trim(cfg%topo%mode) == 'qshe') then
        if (trim(cfg%confinement) /= 'qw' .and. trim(cfg%confinement) /= 'wire') then
          error stop 'validate_semantic: QSHE Z2 requires QW or wire confinement, got ''' // &
            trim(cfg%confinement) // ''''
        end if
      end if

      ! S3: BdG mode requires [bdg] section enabled
      if (trim(cfg%topo%mode) == 'bdg') then
        if (.not. cfg%bdg%enabled) then
          error stop 'validate_semantic: topology mode ''bdg'' requires [bdg] section'
        end if
        ! S4: BdG confinement must be QW or wire
        if (trim(cfg%confinement) /= 'qw' .and. trim(cfg%confinement) /= 'wire') then
          error stop 'validate_semantic: BdG requires QW or wire confinement, got ''' // &
            trim(cfg%confinement) // ''''
        end if
      end if

      ! S5–S7: spectral function parameters
      ! Validate when mode='spectral' (dispatches to run_spectral) or
      ! compute_spectral=.true. (spectral computation requested for any mode).
      if (trim(cfg%topo%mode) == 'spectral' .or. cfg%topo%compute_spectral) then
        if (cfg%topo%spectral_eta <= 0.0_dp) then
          block
            character(len=32) :: buf_eta
            write(buf_eta, '(ES12.4)') cfg%topo%spectral_eta
            error stop 'validate_semantic: spectral_eta must be positive, got ' // trim(buf_eta)
          end block
        end if
        if (cfg%topo%spectral_nk < 1) then
          block
            character(len=16) :: buf_nk
            write(buf_nk, '(I0)') cfg%topo%spectral_nk
            error stop 'validate_semantic: spectral_nk must be >= 1, got ' // trim(buf_nk)
          end block
        end if
        if (cfg%topo%spectral_nE < 1) then
          block
            character(len=16) :: buf_nE
            write(buf_nE, '(I0)') cfg%topo%spectral_nE
            error stop 'validate_semantic: spectral_nE must be >= 1, got ' // trim(buf_nE)
          end block
        end if
      end if

      ! S8: sweep_model must be recognized and match confinement
      ! Only validate when mode='sweep' (the sole consumer of sweep_model).
      ! bhz_analytic is confinement-agnostic (purely analytical).
      if (trim(cfg%topo%mode) == 'sweep') then
        select case (trim(cfg%topo%sweep_model))
        case ('bhz_analytic')
          ! confinement-agnostic, no further check needed
        case ('qw_fukane')
          if (trim(cfg%confinement) /= 'qw') then
            error stop 'validate_semantic: sweep_model ''' // trim(cfg%topo%sweep_model) // &
              ''' requires confinement=''qw'', got ''' // trim(cfg%confinement) // ''''
          end if
        case ('wire_bdg')
          if (trim(cfg%confinement) /= 'wire') then
            error stop 'validate_semantic: sweep_model ''' // trim(cfg%topo%sweep_model) // &
              ''' requires confinement=''wire'', got ''' // trim(cfg%confinement) // ''''
          end if
        case default
          error stop 'validate_semantic: sweep_model ''' // trim(cfg%topo%sweep_model) // &
            ''' not recognized (expected: bhz_analytic, qw_fukane, wire_bdg)'
        end select
      end if

      ! S9: conductance_method must be a recognized enum
      ! Validate when mode='conductance' (dispatches to run_conductance) or
      ! compute_conductance=.true. (conductance computation requested).
      if (trim(cfg%topo%mode) == 'conductance' .or. cfg%topo%compute_conductance) then
        if (trim(cfg%topo%conductance_method) /= 'kubo_chern' .and. &
            trim(cfg%topo%conductance_method) /= 'kubo_berry' .and. &
            trim(cfg%topo%conductance_method) /= 'landauer') then
          error stop 'validate_semantic: conductance_method ''' // &
            trim(cfg%topo%conductance_method) // &
            ''' not recognized (expected: kubo_chern, kubo_berry, landauer)'
        end if
      end if

    case default
      error stop 'validate_semantic: unknown app_name ''' // trim(app_name) // ''''

    end select

  end subroutine validate_semantic

  ! ------------------------------------------------------------------
  ! Map confinement string to a direction character.
  !   'bulk'   -> 'n' (no confinement; bulk is fully diagonalized)
  !   'qw'     -> 'z' (1D confinement along the z growth axis)
  !   'landau' -> 'x' (1D confinement along x; Landau orbitals)
  !   'wire'   -> 'w' (2D confinement in the x-y plane; there is NO single
  !                    confinement axis, so this is a sentinel, not an axis.
  !                    Every wire path branches to wire-specific code and
  !                    never consumes this value. Do NOT add a generic
  !                    `== 'z'` branch without a prior wire guard, and do
  !                    not treat 'w' as 'n' — the wire IS confined.)
  !   other    -> 'n' (default; validator is the primary gatekeeper for
  !                     invalid confinement values)
  ! ------------------------------------------------------------------
  pure function conf_direction(conf) result(d)
    character(len=*), intent(in) :: conf
    character(len=1) :: d

    d = 'n'
    select case(trim(conf))
    case('qw')
      d = 'z'
    case('wire')
      d = 'w'
    case('landau')
      d = 'x'
    end select
  end function conf_direction

  ! ------------------------------------------------------------------
  ! resolve_solver_defaults: single source of truth for the
  ! (confinement, method, mode) -> (resolved method, resolved mode)
  ! mapping used by every eigensolver dispatch site.
  !
  ! AUTO contract (CONTEXT.md):
  !   1. METHOD resolves first. AUTO -> confinement default
  !      (bulk/qw/landau -> DENSE, wire -> FEAST); otherwise the
  !      user's method (trimmed, uppercased) is honored.
  !   2. MODE then resolves METHOD-AWARE. AUTO -> a default compatible
  !      with the resolved method (FEAST -> ENERGY; DENSE -> the
  !      confinement's native mode: bulk FULL, qw/landau INDEX,
  !      wire ENERGY). An explicit mode is honored (trimmed, uppercased).
  !
  ! INVARIANT: AUTO never yields an invalid combination; FEAST+INDEX is
  ! unreachable through AUTO. Explicit FEAST+INDEX set by the user is
  ! still returned faithfully (the caller / validate() rejects it).
  !
  ! defs.f90 is the root module and cannot 'use' the eigensolver module
  ! (which depends on it), so this returns UPPERCASE character strings
  ! ('DENSE'/'FEAST' and 'FULL'/'INDEX'/'ENERGY'), not EIGEN_MODE_*
  ! integer constants. Each dispatch site maps the mode string to its
  ! integer constant via eigen_mode_from_string() in simulation_setup.
  ! ------------------------------------------------------------------
  pure subroutine resolve_solver_defaults(confinement, method, mode, out_method, out_mode)
    character(len=*), intent(in)  :: confinement, method, mode
    character(len=*), intent(out) :: out_method, out_mode

    character(len=len(out_method)) :: rmeth
    character(len=len(out_mode))   :: rmode

    ! --- Resolve METHOD first (confinement default unless overridden) ---
    select case (to_upper_trim(method))
    case ('AUTO')
      select case (to_upper_trim(confinement))
      case ('BULK', 'QW', 'LANDAU')
        rmeth = 'DENSE'
      case ('WIRE')
        rmeth = 'FEAST'
      case default
        ! Unknown confinement: fall back to DENSE. validate() is the
        ! primary gatekeeper for invalid confinement values; this keeps
        ! the resolver total (never error-stops inside a pure routine).
        rmeth = 'DENSE'
      end select
    case ('DENSE', 'FEAST')
      rmeth = to_upper_trim(method)
    case default
      ! Unknown method: fall back to DENSE (validate() rejects upstream).
      rmeth = 'DENSE'
    end select

    ! --- Resolve MODE method-aware ---
    select case (to_upper_trim(mode))
    case ('AUTO')
      select case (rmeth)
      case ('FEAST')
        rmode = 'ENERGY'
      case ('DENSE')
        select case (to_upper_trim(confinement))
        case ('BULK')
          rmode = 'FULL'
        case ('QW', 'LANDAU')
          rmode = 'INDEX'
        case ('WIRE')
          rmode = 'ENERGY'
        case default
          rmode = 'FULL'
        end select
      end select
    case ('FULL', 'INDEX', 'ENERGY')
      rmode = to_upper_trim(mode)
    case default
      rmode = 'FULL'
    end select

    out_method = rmeth
    out_mode   = rmode
  end subroutine resolve_solver_defaults

  ! Uppercase + trim a short ASCII string (defs.f90 cannot use downstream
  ! helpers; this is the canonical local normalizer).
  pure function to_upper_trim(s) result(u)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: u
    integer :: i, c
    u = s
    do i = 1, len(s)
      c = iachar(s(i:i))
      if (c >= iachar('a') .and. c <= iachar('z')) then
        u(i:i) = achar(c - 32)
      end if
    end do
    u = trim(u)
  end function to_upper_trim

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
  ! Type-bound accessor: total number of spatial grid points (nx * ny).
  ! Semantically identical to the standalone grid_ngrid function, but
  ! called as grid%npoints() on polymorphic or direct references.
  ! ------------------------------------------------------------------
  elemental pure function spatial_grid_npoints(self) result(n)
    class(spatial_grid), intent(in) :: self
    integer :: n
    n = self%nx * self%ny
  end function spatial_grid_npoints

  ! ------------------------------------------------------------------
  ! Initialize a spatial_grid from the fields already set in
  ! simulation_config.  This is the primary grid initialization path.
  !
  ! For bulk (confinement='bulk'): ndim=0, nx=ny=1.
  ! For QW  (confinement='qw'): ndim=1, nx=1, ny=fdStep, z(:) copied.
  ! For wire(confinement='wire'): ndim=2, nx=wire_nx, ny=wire_ny,
  !   set grid dimensions only.  Coordinate arrays, material_id,
  !   cut-cell fields and ghost_map are populated by
  !   init_wire_from_config() in the geometry module.
  ! ------------------------------------------------------------------
  subroutine init_grid_from_config(cfg)
    type(simulation_config), intent(inout) :: cfg

    integer :: i, j

    select case (trim(cfg%confinement))
    case ('bulk')
      ! Bulk: single point, no spatial extent
      cfg%grid%ndim = 0
      cfg%grid%nx   = 1
      cfg%grid%ny   = 1
      cfg%grid%dx   = 0.0_dp
      cfg%grid%dy   = 0.0_dp
      ! No coordinate arrays allocated for bulk

    case ('qw')
      ! QW: 1D confinement along z
      cfg%grid%ndim = 1
      cfg%grid%nx   = 1
      cfg%grid%ny   = cfg%fd_step
      cfg%grid%dx   = 0.0_dp
      cfg%grid%dy   = cfg%dz

      ! Copy z-coordinate array from legacy field
      if (allocated(cfg%z)) then
        cfg%grid%z = cfg%z
      end if

      ! Build material_id from int_start_pos/int_end_pos
      if (allocated(cfg%int_start_pos) .and. allocated(cfg%int_end_pos)) then
        allocate(cfg%grid%material_id(cfg%fd_step))
        cfg%grid%material_id = 0
        do i = 1, cfg%num_layers
          cfg%grid%material_id(cfg%int_start_pos(i):cfg%int_end_pos(i)) = i
        end do
      end if

      ! No cut-cell fields for QW
      ! No coords(:,:) or x(:) for QW (nx=1)

    case ('wire')
      ! Wire: 2D confinement in x-y plane.
      ! Set grid dimensions only; coordinate arrays, material_id,
      ! cut-cell fields, and ghost_map are populated by
      ! init_wire_from_config() in the geometry module.
      cfg%grid%ndim = 2
      cfg%grid%nx   = cfg%wire%nx
      cfg%grid%ny   = cfg%wire%ny
      cfg%grid%dx   = cfg%wire%dx
      cfg%grid%dy   = cfg%wire%dy

    case ('landau')
      ! Landau: 1D x-discretization for orbital Landau quantization
      ! Grid centered at x=0 so that ky=0 places the oscillator at the
      ! center of the domain, maximising flatness of Landau levels vs ky.
      cfg%grid%ndim = 1
      cfg%grid%nx   = cfg%landau%nx
      cfg%grid%ny   = 1
      cfg%grid%dx   = cfg%landau%width / real(cfg%landau%nx - 1, kind=dp)
      cfg%grid%dy   = 0.0_dp

      ! Build x-coordinate array centered at origin [-W/2, +W/2]
      if (.not. allocated(cfg%grid%x)) then
        allocate(cfg%grid%x(cfg%landau%nx))
        do i = 1, cfg%landau%nx
          cfg%grid%x(i) = (i - 1) * cfg%grid%dx - 0.5_dp * cfg%landau%width
        end do
      end if

    end select

  end subroutine init_grid_from_config


end module
