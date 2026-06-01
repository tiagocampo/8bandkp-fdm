# Plan: Input Parser TOML Refactor

Date: 2026-05-27
Branch: feat/richardson-observables-convergence (will move to feat/input-parser-toml)

## Summary

Replace the monolithic 1,369-line `input_parser.f90` (36 BACKSPACE calls, fragile ordering, mode-dependent overwrites) with a clean TOML-based parser using the `toml-f` library. Full cleanup — no backward compatibility.

## Decisions (from grilling session)

| # | Decision | Choice |
|---|----------|--------|
| D1 | TOML library | `toml-f` v0.4.2, git submodule in `subprojects/toml-f` |
| D2 | Config filename | `input.toml` (hard rename from `input.cfg`) |
| D3 | CLI argument | No `--config` flag (hard-code `input.toml`) |
| D4 | Material layers | `[[material]]` array-of-tables |
| D5 | Wave vector | `[wave_vector]` sub-table with `mode`, `max`, `step` |
| D6 | Confinement | String enum: `"bulk"`, `"qw"`, `"wire"`, `"landau"` |
| D7 | Wire params | `[wire]` sub-table + `[[region]]` array-of-tables |
| D8 | Band counts | `[bands]` sub-table with `num_cb`, `num_vb` |
| D9 | Enable flags | Section presence = enabled (no `enabled` keys) |
| D10 | Key naming | `snake_case` for all TOML keys |
| D11 | Error handling | Fail-fast for required, silent defaults for optional |
| D12 | Echo output | Drop the per-field `print *` echo entirely |
| D13 | Type restructure | Full: sub-types mirror TOML sections |
| D14 | fdStep | Top-level singleton; add `ngrid` computed field |
| D15 | C3 approach | Keep `fdStep` as user input, add `ngrid` computed field |
| D16 | Cleanup scope | Full: all configs, tests, lectures, docs updated |
| D17 | Backward compat | None. Clean break. |

## TOML Schema (canonical)

```toml
# ===== Required top-level =====
confinement = "qw"              # "bulk" | "qw" | "wire" | "landau"
FDorder = 6                     # 2 | 4 | 6 | 8 | 10
fd_step = 200                   # QW grid points (top-level, wire/landau ignore)

[wave_vector]
mode = "kx"                     # "k0" | "kx" | "kxky" | "kpar" | ...
max = 0.1
step = 0.01

[bands]
num_cb = 2
num_vb = 6

# ===== QW geometry (confinement = "qw") =====
[[material]]
name = "GaAs"
z_min = 0.0
z_max = 100.0

[[material]]
name = "InAs"
z_min = 30.0
z_max = 70.0

# ===== Wire geometry (confinement = "wire") =====
[wire]
nx = 50
ny = 50
dx = 0.5
dy = 0.5

[wire.geometry]
shape = "circle"                # "circle" | "rectangle" | "hexagon" | "polygon"
radius = 25.0                   # for circle/hexagon

# OR:
# [wire.geometry]
# shape = "rectangle"
# width = 20.0
# height = 30.0

# OR:
# [wire.geometry]
# shape = "polygon"
# vertices = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]

[[region]]
material = "GaAs"
inner = 0.0
outer = 25.0

# ===== Landau geometry (confinement = "landau") =====
# [landau]
# nx = 100
# width = 2000.0
# sweep = "ky"                   # "ky" | "kz" | "B"
# material = "InAs"

# ===== Optional physics sections (presence = enabled) =====

[external_field]
type = "EF"
value = 0.0

[b_field]
components = [0.0, 0.0, 5.0]   # [Bx, By, Bz] in Tesla
g_factor = 2.0

[strain]
substrate = "GaAs"

[sc]
max_iterations = 100
tolerance = 1e-6
mixing_alpha = 0.3
diis_history = 5
temperature = 300.0
fermi_mode = "charge_neutrality"  # "charge_neutrality" | "fixed"
fermi_level = 0.0
num_kpar = 50
kpar_max = 0.5
bc_type = "DD"                    # "DD" | "DN" | "ND" | "NN"
bc_left = "D"
bc_right = "D"

[[doping]]
layer = 1
ND = 1e18
NA = 0.0

# [[doping]]
# layer = 2
# type = "delta"
# NS = 1e12
# fwhm = 5.0
# pos = 50.0

[topology]
mode = "qhe"                      # "qhe" | "qshe" | "bdg"
compute_chern = true
compute_hall = true
qwz_u = 1.0
compute_z2 = false
extract_edge_states = false
edge_E_window = 0.01
compute_ldos = true
ldos_eta = 0.01
ldos_E_range = [-1.0, 1.0]
ldos_num_E = 100

[topology.gap_sweep]
B_range = [0.0, 2.0, 20]         # [min, max, npts]
mu_range = [-0.5, 0.5, 10]
sweep_model = "bhz"               # "bhz" | "wire_bdg" | "qw_fukane"

[topology.conductance]
method = "landauer"
berry_nk = 100
landauer_energy = 0.0

[topology.spectral]
k_range = [-1.0, 1.0, 50]
E_range = [-2.0, 2.0, 200, 0.01]

[bdg]
mu = 0.0
delta_0 = 0.00025
gauge = "landau_x"
B_sweep = [0.0, 1.0, 0.01]

[optics]
linewidth_lorentzian = 5.0
linewidth_gaussian = 0.0
refractive_index = 3.5
E_range = [0.5, 1.5, 200]
temperature = 300.0
carrier_density = 0.0
gain_enabled = false
gain_carrier_density = 0.0
ISBT = false
spontaneous = false
spin_resolved = false

[exciton]
Ry = 0.005
a0 = 100.0

[scattering]
temperature = 300.0
deformation_potential = 7.0
phonon_energy = 0.03

[feast]
emin = 0.0                       # 0 = auto
emax = 0.0                       # 0 = auto
m0 = 0                           # 0 = auto (2*nev)
```

## simulation_config Type Restructure

### New sub-types

```fortran
! Replaces scattered wire_* fields
type :: wire_config
  integer            :: nx = 0
  integer            :: ny = 0
  real(dp)           :: dx = 0.0_dp
  real(dp)           :: dy = 0.0_dp
  type(wire_geometry) :: geom
  type(region_spec), allocatable :: regions(:)
  integer            :: num_regions = 0
end type

! Replaces scattered landau_* fields
type :: landau_config
  integer          :: nx = 100
  real(dp)         :: width = 2000.0_dp
  character(len=4) :: sweep = 'ky'
  character(len=255) :: material = ''  ! single material name
end type

! Replaces flat b_field/EFtype/Evalue
type :: external_field_config
  character(len=2) :: type = 'EF'   ! "EF" for electric field
  real(dp)         :: value = 0.0_dp
end type

! Replaces scattered b_field/bdg%B_vec
type :: b_field_config
  real(dp)         :: components(3) = 0.0_dp  ! Bx, By, Bz
  real(dp)         :: g_factor = 2.0_dp
end type

! Replaces flat numcb/numvb
type :: bands_config
  integer :: num_cb = 2
  integer :: num_vb = 6
end type

! Replaces flat waveVector/waveVectorMax/waveVectorStep
type :: wave_vector_config
  character(len=8) :: mode = 'k0'
  real(dp)         :: max = 0.0_dp
  real(dp)         :: step = 0.01_dp
end type

! Replaces flat feast_emin/feast_emax/feast_m0
type :: feast_config
  real(dp) :: emin = 0.0_dp
  real(dp) :: emax = 0.0_dp
  integer  :: m0 = 0
end type
```

### Updated simulation_config

```fortran
type simulation_config
  ! ---- User inputs (from TOML) ----
  character(len=8)   :: confinement = 'bulk'  ! "bulk"|"qw"|"wire"|"landau"
  integer            :: FDorder = 2
  integer            :: fd_step = 1            ! user-specified (QW grid points)
  type(wave_vector_config) :: wave_vector
  type(bands_config)       :: bands
  type(external_field_config) :: external_field
  type(b_field_config)    :: b_field
  type(strain_config)     :: strain
  type(feast_config)      :: feast

  ! Mode-specific geometry (only one is active)
  type(wire_config)       :: wire               ! confinement = "wire"
  type(landau_config)     :: landau             ! confinement = "landau"

  ! Material layers (QW/bulk) — populated from [[material]]
  character(len=255), allocatable :: material_names(:)
  real(dp), allocatable :: z_min(:), z_max(:)
  integer               :: num_layers = 1

  ! Physics subsystems (presence = enabled)
  type(sc_config)           :: sc
  type(doping_spec), allocatable :: doping(:)
  type(topology_config)     :: topo
  type(bdg_config)          :: bdg
  type(optics_config)       :: optics
  type(exciton_config)      :: exciton
  type(scattering_config)   :: scattering

  ! ---- Computed fields (derived from inputs) ----
  integer            :: ngrid = 1         ! computed grid size (was fdStep after overwrites)
  integer            :: evnum = 8         ! bands%num_cb + bands%num_vb
  character(len=1)   :: conf_dir = 'n'    ! 'n','z','x' based on confinement
  real(dp)           :: total_size = 0.0_dp
  real(dp)           :: delta = 0.0_dp
  real(dp)           :: dz = 0.0_dp
  real(dp), allocatable :: z(:)
  integer, allocatable   :: int_start_pos(:), int_end_pos(:)
  type(paramStruct), allocatable :: params(:)
  type(bir_pikus_blocks) :: strain_blocks
  type(spatial_grid)      :: grid
  real(dp)               :: sc_potential_shift = 0.0_dp

  ! G-factor (only for gfactorCalculation executable)
  integer :: which_band = 0
  integer :: band_idx = 1
contains
  final :: simulation_config_finalize
  procedure :: validate => simulation_config_validate
end type
```

## Migration mapping: old → new TOML keys

| Old label | New TOML key | Section |
|-----------|-------------|---------|
| `waveVector:` | `mode` | `[wave_vector]` |
| `waveVectorMax:` | `max` | `[wave_vector]` |
| `waveVectorStep:` | `step` | `[wave_vector]` |
| `confinement:` | `confinement` | top-level (integer→string) |
| `FDstep:` | `fd_step` | top-level |
| `FDorder:` | `FDorder` | top-level |
| `numLayers:` | _(derived from `[[material]]` count)_ | — |
| `materialN:` | `name`, `z_min`, `z_max` | `[[material]]` |
| `numcb:` | `num_cb` | `[bands]` |
| `numvb:` | `num_vb` | `[bands]` |
| `ExternalField:` + `EFtype:` | `type`, `value` | `[external_field]` |
| `Evalue:` | `value` | `[external_field]` |
| `b_field:` | `components` | `[b_field]` |
| `g_factor:` | `g_factor` | `[b_field]` |
| `SC:` | _(section presence)_ | `[sc]` |
| `SC_max_iter:` | `max_iterations` | `[sc]` |
| `SC_tolerance:` | `tolerance` | `[sc]` |
| `SC_mixing:` | `mixing_alpha` | `[sc]` |
| `SC_diis:` | `diis_history` | `[sc]` |
| `SC_temperature:` | `temperature` | `[sc]` |
| `SC_fermi_mode:` | `fermi_mode` | `[sc]` (integer→string) |
| `SC_fermi_level:` | `fermi_level` | `[sc]` |
| `SC_num_kpar:` | `num_kpar` | `[sc]` |
| `SC_kpar_max:` | `kpar_max` | `[sc]` |
| `SC_bc:` | `bc_type` | `[sc]` |
| `SC_bc_left:` | `bc_left` | `[sc]` |
| `SC_bc_right:` | `bc_right` | `[sc]` |
| `doping:` | `layer`, `ND`, `NA` | `[[doping]]` |
| `delta<N>:` | `layer`, `type`, `NS`, `fwhm`, `pos` | `[[doping]]` |
| `whichBand:` | `which_band` | top-level |
| `bandIdx:` | `band_idx` | top-level |
| `wire_nx:` | `nx` | `[wire]` |
| `wire_ny:` | `ny` | `[wire]` |
| `wire_dx:` | `dx` | `[wire]` |
| `wire_dy:` | `dy` | `[wire]` |
| `wire_shape:` | `shape` | `[wire.geometry]` |
| `wire_radius:` | `radius` | `[wire.geometry]` |
| `wire_width:` | `width` | `[wire.geometry]` |
| `wire_height:` | `height` | `[wire.geometry]` |
| `wire_polygon:` + vertices | `vertices` | `[wire.geometry]` |
| `numRegions:` | _(derived from `[[region]]` count)_ | — |
| `region:` | `material`, `inner`, `outer` | `[[region]]` |
| `landau_nx:` | `nx` | `[landau]` |
| `landau_width:` | `width` | `[landau]` |
| `landau_sweep:` | `sweep` | `[landau]` |
| `strainSubstrate:` | `substrate` | `[strain]` |
| `topology:` + `T/F` | _(section presence)_ | `[topology]` |
| `optics:` + `T/F` | _(section presence)_ | `[optics]` |
| `exciton:` + `T/F` | _(section presence)_ | `[exciton]` |
| `scattering:` + `T/F` | _(section presence)_ | `[scattering]` |
| `feast_emin:` | `emin` | `[feast]` |

## Implementation Phases

### Phase 1: Build infrastructure (toml-f integration)

**Step 1.1**: Add toml-f as git submodule
```bash
mkdir -p subprojects
git submodule add https://github.com/toml-f/toml-f subprojects/toml-f
git submodule update --init --recursive
```

**Step 1.2**: Update CMakeLists.txt to build toml-f as static library
- Add `add_subdirectory(subprojects/toml-f)` to root CMakeLists.txt
- Link `toml-f-lib` to the parser library/executables
- Verify `-std=f2018` compatibility (toml-f is pure F2008)

**Step 1.3**: Verify build works
```bash
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl
cmake --build build
```

### Phase 2: Restructure derived types in defs.f90

**Step 2.1**: Add new sub-types (`wire_config`, `landau_config`, `external_field_config`, `b_field_config`, `bands_config`, `wave_vector_config`, `feast_config`) in `defs.f90`

**Step 2.2**: Update `simulation_config` to use new sub-types
- Remove scattered `wire_*`, `landau_*`, `b_field`, `EFtype`, `Evalue`, `waveVector*`, `numcb`, `numvb`, `feast_*` flat fields
- Add computed `ngrid` field
- Add `conf_dir` (replaces `confDir`)
- Keep `params`, `material_names`, `z_min`, `z_max`, `z`, `int_start_pos`, `int_end_pos`, etc. as computed fields

**Step 2.3**: Update `simulation_config_finalize` for new sub-types (if any have allocatable components)

**Step 2.4**: Update `simulation_config_validate` for new field names

### Phase 3: Rewrite input_parser.f90

**Step 3.1**: Replace entire file with toml-f based parser
- `read_config(filename)` returns `type(simulation_config)`
- Uses `toml_parse` to load file into `toml_table`
- Helper: `require_value(table, key, val, section)` — fail-fast for required fields
- Helper: `get_value(table, key, val, default)` — silent defaults for optional fields
- Sub-parser functions per section:
  - `parse_wave_vector(table, cfg)`
  - `parse_bands(table, cfg)`
  - `parse_materials(table, cfg)` — handles bulk/QW `[[material]]`
  - `parse_wire(table, cfg)` — handles `[wire]` + `[[region]]`
  - `parse_landau(table, cfg)` — handles `[landau]`
  - `parse_external_field(table, cfg)`
  - `parse_b_field(table, cfg)`
  - `parse_sc(table, cfg)` + `parse_doping(table, cfg)`
  - `parse_topology(table, cfg)` + sub-parsers for gap_sweep/conductance/spectral
  - `parse_bdg(table, cfg)`
  - `parse_optics(table, cfg)`
  - `parse_exciton(table, cfg)`
  - `parse_scattering(table, cfg)`
  - `parse_feast(table, cfg)`
  - `parse_strain(table, cfg)`

**Step 3.2**: Post-parse computation
- Compute `ngrid` from mode: QW→fd_step, wire→wire%ny, landau→landau%nx, bulk→1
- Compute `evnum` = bands%num_cb + bands%num_vb
- Compute `conf_dir`: bulk→'n', qw→'z', wire→'z', landau→'x'
- For QW: compute `z(:)`, `delta`, `dz`, `int_start_pos`, `int_end_pos` from material layers + fd_step
- For wire: compute `startPos`/`endPos` from wire geometry
- Call `paramDatabase` to populate `params(:)`
- Call `simulation_config_validate` for cross-field validation

### Phase 4: Update all consumers

**Step 4.1**: Update `src/apps/main.f90`
- `cfg%fdStep` → `cfg%ngrid`
- `cfg%waveVector` → `cfg%wave_vector%mode`
- `cfg%confinement` → string comparison (`cfg%confinement == "qw"`)
- `cfg%numcb` → `cfg%bands%num_cb`
- `cfg%numvb` → `cfg%bands%num_vb`
- `cfg%wire_nx` → `cfg%wire%nx`
- `cfg%materialN` → `cfg%material_names`
- `cfg%startPos` → `cfg%z_min`
- `cfg%endPos` → `cfg%z_max`
- `cfg%ExternalField` → `cfg%external_field%type`
- `cfg%Evalue` → `cfg%external_field%value`
- `cfg%confDir` → `cfg%conf_dir`
- ~20 references to update

**Step 4.2**: Update `src/apps/main_gfactor.f90`
- Same field renaming pattern
- `cfg%whichBand` → `cfg%which_band`
- `cfg%bandIdx` → `cfg%band_idx`
- ~10 references

**Step 4.3**: Update `src/apps/main_optics.f90`
- Same field renaming pattern
- ~10 references

**Step 4.4**: Update `src/apps/main_topology.f90`
- Same field renaming pattern

**Step 4.5**: Update `src/core/simulation_setup.f90` (10 references)
- `cfg%confinement` → string comparison
- `cfg%fdStep` → `cfg%ngrid`
- `cfg%numvb` → `cfg%bands%num_vb`
- `cfg%ExternalField`/`cfg%EFtype` → `cfg%external_field%type`
- `cfg%Evalue` → `cfg%external_field%value`

**Step 4.6**: Update physics modules that consume `simulation_config` (22 references across 5 files)
- `src/physics/sc_loop.f90` (8 refs) — `cfg%fdStep` → `cfg%ngrid`, SC field renames
- `src/physics/confinement_init.f90` (7 refs) — field renames
- `src/physics/green_functions.f90` (4 refs) — field renames
- `src/physics/hamiltonianConstructor.f90` (2 refs) — `cfg%materialN` → `cfg%material_names`
- `src/physics/bdg_hamiltonian.f90` (1 ref) — field rename

**Step 4.7**: Update `src/io/outputFunctions.f90` (2 references)
- Field renames for output

**Step 4.6**: Update unit tests
- All pFUnit tests in `tests/unit/` that reference `simulation_config` fields
- Update constructor calls and field accesses

### Phase 5: Convert all config files to TOML

**Step 5.1**: Write a conversion script (`scripts/convert_cfg_to_toml.py`)
- Parse old `label: value` format
- Output TOML using mapping table above
- Handle mode-dependent branching (bulk/QW/wire/Landau)
- Handle optional blocks

**Step 5.2**: Convert all 88 configs in `tests/regression/configs/`
- Run converter script
- Rename `.cfg` → `.toml`
- Verify each config by hand

**Step 5.3**: Convert validation configs in `validation/`

### Phase 6: Update test infrastructure

**Step 6.1**: Update 32 test shell scripts in `tests/integration/`
- Change `cp $config input.cfg` → `cp $config input.toml`

**Step 6.2**: Update 15 lecture scripts in `scripts/`
- Update config file references
- Update any inline config generation

**Step 6.3**: Update `tests/regression/compare_output.py` if needed

### Phase 7: Update documentation

**Step 7.1**: Rewrite `docs/reference/input-reference.md` for TOML format

**Step 7.2**: Update README if it references `input.cfg`

**Step 7.3**: Update CLAUDE.md (input file section, running section)

**Step 7.4**: Update any validation docs

### Phase 8: Verification

**Step 8.1**: Build and run all tests
```bash
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl \
    -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-4.16/cmake
cmake --build build
ctest --test-dir build -j4 --output-on-failure
```

**Step 8.2**: Run validation pipeline
```bash
source validation/kdotpy_env/bin/activate
python3 validation/run_all.py
```

**Step 8.3**: Run lecture scripts to verify configs work
```bash
python3 scripts/lecture_00_quickstart.py
python3 scripts/lecture_01_bulk.py
# ... etc
```

## File change estimate

| Area | Files | Scope |
|------|-------|-------|
| `subprojects/toml-f/` | 1 (submodule) | git submodule add |
| `CMakeLists.txt` | 2 (root + src/io) | toml-f linking |
| `src/core/defs.f90` | 1 | new sub-types, restructured simulation_config |
| `src/io/input_parser.f90` | 1 | full rewrite (~500 lines, down from 1369) |
| `src/apps/main*.f90` | 4 | field renaming |
| `src/core/simulation_setup.f90` | 1 | field renaming (10 refs) |
| `src/physics/*.f90` | ~7 | field renaming |
| `src/io/outputFunctions.f90` | 1 | field renaming |
| `tests/unit/*.pf` | ~15 | field renaming |
| `tests/regression/configs/*` | 88 | convert to TOML |
| `tests/integration/*.sh` | 32 | `input.cfg` → `input.toml` |
| `scripts/lecture_*.py` | 15 | config references |
| `validation/` | ~12 | config conversion |
| `docs/` | 3 | input reference, README, CLAUDE.md |
| `scripts/convert_cfg_to_toml.py` | 1 | new conversion script |

**Total: ~175 files touched**

## Risks and mitigations

| Risk | Mitigation |
|------|-----------|
| toml-f build fails with gfortran + MKL flags | Test Phase 1 first, before any other changes |
| Consumer field rename breaks physics | Rename is mechanical; verify with existing tests |
| Config conversion errors | Conversion script + manual spot-check of all 88 configs |
| TOML array-of-tables ordering edge case | toml-f preserves order; test with polygon wire configs |
| `simulation_config` finalizer misses new allocatable components | Review each new sub-type for allocatable components |
| `confinement` integer→string breaks consumer `if` chains | Grep all `cfg%confinement == N` and convert to string comparisons |
