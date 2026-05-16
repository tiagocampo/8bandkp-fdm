# Test Coverage Completion

Date: 2026-05-07
Status: Draft

## Problem

Three backlog items from BACKLOG.md (#4, #8, #37) leave known gaps in test coverage. No physics changes needed — only golden data generation, test script writing, and `contiguous` attribute annotations.

## Scope

### T1: G-factor regression tests (backlog #4)

Five g-factor configs exist in `tests/regression/configs/` but lack golden reference data and CMake test registration. Only `gfactor_qw_cb.cfg` and `gfactor_qw_no_optics.cfg` have full regression coverage.

**Configs needing golden data and test registration:**

| Config | Material | System | Expected g-values |
|--------|----------|--------|-------------------|
| `gfactor_bulk_gaas_cb.cfg` | GaAs | Bulk CB | g_x = g_y = g_z = -0.315 |
| `gfactor_bulk_gaas_vb.cfg` | GaAs | Bulk VB | g_z = +5.653 |
| `gfactor_bulk_gaasw_cb.cfg` | GaAsW | Bulk CB | g_x = g_y = g_z = -0.322 |
| `gfactor_bulk_inasw_cb.cfg` | InAsW | Bulk CB | g_x = g_y = g_z = -14.858 |
| `gfactor_qw_vb.cfg` | InAs/GaSb/AlSb | QW VB | g_z = -20.891 |

**Steps per config:**
1. Run `gfactorCalculation` with the config, capture output
2. Extract g-factor values from output using the existing `compare_output.py` framework
3. Store golden data in `tests/regression/data/<config_name>/`
4. Add CMake `add_test` entry in `tests/CMakeLists.txt`

**Registration pattern:** Follow the existing `regression_gfactor_cb` test pattern:
- Test shell script: `tests/integration/test_gfactor.sh` (reuse existing, it takes config path as argument)
- Or create lightweight wrapper scripts per config if needed

### T2: Integration tests + regression datasets (backlog #8)

Three integration tests and two regression datasets for wire and SC paths.

**Integration tests to create:**

1. **Wire hexagon geometry** (`test_wire_hexagon.sh`)
   - Config: hexagonal cross-section wire
   - Verify: geometry setup, material mapping, eigenvalue count
   - Check that eigenvectors and eigenvalues match golden data within tolerance

2. **Wire strain** (`test_wire_strain.sh`)
   - Config: InAs/GaAs core-shell with strain
   - Verify: strain tensor components, band edge shifts, eigenvalues
   - Reuses existing `wire_inas_gaas_strain.cfg`

3. **SC wire** (`test_wire_sc.sh`)
   - Config: wire with self-consistent Poisson
   - Verify: convergence, charge density, potential profile, eigenvalues

**Regression datasets to generate:**
- Wire hexagon: `tests/regression/data/wire_hexagon/`
- Wire strain: already exists as `wire_inas_gaas_core_shell/`; verify completeness
- SC wire: `tests/regression/data/wire_sc/` (new)

### T3: Contiguous attribute gaps (backlog #37)

Three known hot-path array arguments missing the `contiguous` attribute, as documented in CLAUDE.md:

1. **`ZB8bandQW` profile/kpterms** (`hamiltonianConstructor.f90:50-51`)
   - `profile(:,:)` and `kpterms(:,:,:)` should be declared `contiguous`

2. **`dns(:,:)` in `dnscsr_z_mkl`** (`utils.f90:19`)
   - The `dns` array argument should be declared `contiguous`

3. **`psi(:)` in `spin_weights`** (`spin_projection.f90:14`)
   - The `psi` array argument should be declared `contiguous`

**Fix:** Add `contiguous` attribute to each assumed-shape array declaration. All callers already pass contiguous arrays (allocatable or section of contiguous), so this is a declaration-only change.

## Verification

1. All 66 existing tests still pass
2. New gfactor regression tests pass (5 new tests)
3. New integration tests pass (3 new tests)
4. Build succeeds with no new warnings from `contiguous` additions
5. Total test count increases to 74+
