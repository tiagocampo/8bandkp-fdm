# PRD: Eigensolver Standardization — Unified Interface, QW CSR Builder, ARPACK Removal

## Problem Statement

The eigensolver system has two parallel architectures that never meet. The polymorphic `eigensolver_base` interface (`feast_solver_t`, `dense_lapack_solver_t`) is used exclusively for wire mode. All bulk, QW, and Landau paths bypass it entirely — they call LAPACK routines (`zheev`, `zheevd`, `zheevx`) directly with 25+ duplicated call sites, 17 independent workspace query blocks, and 4× copy-pasted OpenMP parallel k-sweep patterns.

ARPACK infrastructure (~280 lines) is compiled behind `#ifdef USE_ARPACK` but unreachable from any user input — it is dead code that adds maintenance burden and confuses the `USE_MKL_FEAST` guard nesting. The LDOS computation (`compute_ldos_csr`) is trapped behind the same ARPACK guard even though it only uses PARDISO, not ARPACK.

The user-facing `[feast]` TOML section is FEAST-specific naming that controls both the solver method (via `m0 < 0` hack for dense fallback) and energy window parameters. There is no way to express "use dense LAPACK with energy-window mode" or "use FEAST on QW geometry" — the method and geometry are implicitly coupled.

The QW Hamiltonian is always built as a dense matrix, even though it is ~1% fill (~6,000 non-zeros in a 640,000-element 800×800 matrix). This prevents using sparse solvers (FEAST) on QW problems, which would benefit large grids.

## Solution

Replace the two parallel architectures with one unified eigensolver interface:

1. **Remove ARPACK entirely** — delete all ARPACK code, unguard LDOS/PARDISO from `USE_ARPACK`, remove CMake detection.
2. **Redesign the polymorphic interface** — two type-bound methods `solve_dense(H, config, result)` and `solve_sparse(H_csr, config, result)` on `eigensolver_base`, with config-driven mode dispatch (`FULL` / `INDEX` / `ENERGY`).
3. **Rename `[feast]` to `[solver]`** — new TOML section with `method` and `mode` fields, smart defaults per confinement mode when absent.
4. **Build a QW CSR Hamiltonian builder** — proper sparse construction from 1D kpterms, with a workspace cache for fast k-point rebuilds. Enables FEAST on QW.
5. **Route all 25+ call sites** through the polymorphic dispatch, eliminating duplicated workspace queries and LAPACK call patterns.
6. **Migrate all configs, scripts, and docs** — update every `[feast]` section, lecture document, and reference.
7. **Comprehensive test coverage** — all solver × mode × geometry combinations.

## User Stories

### ARPACK Removal

1. As a developer reading `eigensolver.f90`, I want no ARPACK code to exist, so that I don't have to understand `#ifdef USE_ARPACK` nesting when maintaining the solver
2. As a developer reading `linalg.f90`, I want `pardiso_c` available unconditionally (not behind `USE_ARPACK`), so that the LDOS computation doesn't depend on an ARPACK installation
3. As a developer reading `green_functions.f90`, I want `compute_ldos_csr` available without `USE_ARPACK`, so that LDOS functionality works with just MKL
4. As a developer building the project, I want no `find_package(ARPACK)` in CMake, so that the build system is simpler
5. As a developer running unit tests, I want no `test_arpack_fallback` test, so that the test suite doesn't reference a non-existent code path

### Unified Interface

6. As a physicist running a QW band structure calculation, I want to choose between dense LAPACK and sparse FEAST via a TOML parameter, so that I can benchmark which solver is faster for my grid size
7. As a physicist running a bulk calculation, I want the solver to default to dense full diagonalization, so that I get all 8 eigenvalues without specifying solver parameters
8. As a physicist running a QW k-sweep, I want the solver to default to dense index mode, so that I get the bands around the valence-conduction edge efficiently
9. As a physicist running a wire calculation, I want the solver to default to FEAST energy mode, so that sparse diagonalization works out of the box
10. As a physicist running a Landau calculation, I want the solver to default to dense index mode, so that the B-field sweep works as before
11. As a physicist, I want to set `method = "dense"` with `mode = "energy"` to use LAPACK with an energy window, so that I can compare dense window results with FEAST
12. As a physicist, I want to set `method = "feast"` with `mode = "full"` to get all eigenvalues from FEAST, so that I can use FEAST even when I don't know the energy window
13. As a physicist, I want `method = "feast"` + `mode = "index"` to be rejected with a clear error message, so that I don't accidentally request an impossible computation
14. As a developer adding a new confinement mode, I want to specify a default solver method and mode, so that the new mode works without user configuration
15. As a developer adding a new eigensolver variant (e.g., batched LAPACK, ELPA), I want to add it to `make_eigensolver` and have it automatically available to all confinement modes, so that I don't have to touch every call site
16. As a developer maintaining `simulation_setup.f90`, I want all eigensolve calls to go through `solver%solve_dense` or `solver%solve_sparse`, so that I don't have duplicated workspace query and LAPACK call patterns

### TOML Input Standardization

17. As a physicist writing an input config, I want a `[solver]` section with `method` and `mode` fields, so that the parameter names reflect what they control (not just FEAST)
18. As a physicist running a wire simulation, I want to omit the `[solver]` section entirely and get FEAST energy mode automatically, so that existing wire workflows don't change
19. As a physicist running a QW simulation, I want to omit the `[solver]` section entirely and get dense index mode automatically, so that existing QW workflows don't change
20. As a physicist with an old config containing `[feast]`, I want a clear deprecation warning pointing to `[solver]`, so that I know how to migrate
21. As a physicist using the config converter script, I want it to emit `[solver]` instead of `[feast]`, so that converted configs use the new format

### QW CSR Builder

22. As a physicist running a large QW simulation (N > 500), I want to use FEAST on the sparse QW Hamiltonian, so that I avoid the O(N²) memory of dense storage
23. As a physicist benchmarking QW solvers, I want the QW CSR Hamiltonian to produce identical eigenvalues to the dense builder, so that I can trust the sparse results
24. As a developer reading the QW CSR builder, I want it to reuse the same COO insertion layer as the wire builder, so that the block-table-driven assembly is consistent across geometries
25. As a developer running QW k-sweeps, I want a workspace cache that reuses CSR structure across k-points, so that the rebuild cost is O(NNZ) not O(NNZ log NNZ)

### Documentation

26. As a physicist reading the input reference, I want the `[solver]` section documented with all fields and valid combinations, so that I know what to write in my config
27. As a physicist reading the quantum wire lecture, I want TOML examples to use `[solver]` not `[feast]`, so that the tutorial matches the current interface
28. As a physicist reading the numerical methods lecture, I want the eigensolver section to describe both dense and FEAST paths with the unified mode system, so that I understand the available options
29. As a developer reading CLAUDE.md, I want the eigensolver architecture section to reflect the unified interface, so that the documentation is trustworthy
30. As a developer reading AGENTS.md, I want the solver module descriptions to reflect the new naming and capabilities, so that I can navigate the codebase correctly

### Testing

31. As a developer verifying the eigensolver refactor, I want a test for each valid solver × mode combination, so that I know all paths produce correct results
32. As a developer verifying QW CSR correctness, I want a unit test that builds both dense and CSR QW Hamiltonians and compares eigenvalues, so that the CSR builder is validated against the trusted dense path
33. As a developer verifying smart defaults, I want an integration test that runs each confinement mode without `[solver]` and checks correct results, so that the default dispatch is correct
34. As a developer verifying the FEAST+INDEX rejection, I want a unit test that confirms `error stop` fires, so that the validation is enforced

## Implementation Decisions

### Module: `eigensolver` (src/math/eigensolver.f90) — Core interface redesign

- **Remove `arpack_solver_t`**, all `#ifdef USE_ARPACK` blocks, `solve_arpack_dispatch`, `solve_arpack`, ARPACK case in `make_eigensolver`. Delete ~280 lines.
- **Add mode constants**: `EIGEN_MODE_FULL = 1`, `EIGEN_MODE_INDEX = 2`, `EIGEN_MODE_ENERGY = 3` as public named constants.
- **Extend `eigensolver_config`** with `mode = EIGEN_MODE_INDEX` (default), `il = 1`, `iu = 8` for INDEX mode.
- **Add `solve_dense` deferred method** to `eigensolver_base` accepting `complex(dp), contiguous, intent(in) :: H(:,:)`.
- **Rename existing `solve` to `solve_sparse`** (accepts `type(csr_matrix)`). Keep `solve` as a legacy alias calling `solve_sparse` during transition.
- **`dense_lapack_solver_t%solve_dense`**: dispatches on `config%mode`:
  - `EIGEN_MODE_FULL` → `zheev` (all eigenvalues/vectors)
  - `EIGEN_MODE_INDEX` → `zheevx range='I'` with `il=config%il, iu=config%iu`
  - `EIGEN_MODE_ENERGY` → `zheevx range='V'` with `vl=config%emin, vu=config%emax`
- **`dense_lapack_solver_t%solve_sparse`**: converts CSR→dense, then calls `solve_dense`.
- **`feast_solver_t%solve_sparse`**: dispatches on `config%mode`:
  - `EIGEN_MODE_FULL` → Gershgorin-bounded window, then `zfeast_hcsrev`
  - `EIGEN_MODE_ENERGY` → native FEAST with `[emin, emax]`
  - `EIGEN_MODE_INDEX` → `error stop` (rejected at validation)
- **`feast_solver_t%solve_dense`**: converts dense→CSR, then calls `solve_sparse`.
- **`eigensolver_config_validate(config, method)`**: rejects FEAST+INDEX, validates emin<emax for ENERGY, validates nev>0 for INDEX, validates il/iu range.
- **Remove `solve_dense_lapack`** public procedure (replaced by the type-bound `solve_dense`).
- **Keep `feast_workspace`** caching as-is — it already caches upper-triangle CSR structure and is orthogonal to the interface change.

### Module: `definitions` (src/core/defs.f90) — Config type rename

- **Replace `feast_config` with `solver_config`**:
  ```
  type :: solver_config
    character(len=10) :: method = 'AUTO'   ! AUTO, DENSE, FEAST
    character(len=10) :: mode   = 'AUTO'   ! AUTO, FULL, INDEX, ENERGY
    real(kind=dp)     :: emin   = 0.0_dp   ! 0 = auto
    real(kind=dp)     :: emax   = 0.0_dp   ! 0 = auto
    integer           :: m0     = 0         ! 0 = auto (FEAST subspace)
  end type
  ```
- **Update `simulation_config`**: replace `type(feast_config) :: feast` with `type(solver_config) :: solver`.
- **Update validation**: checks I11/I12 reference `solver%emin/emax/m0`. Add: FEAST+INDEX rejection, valid method/mode string validation.
- **Smart defaults** applied in `simulation_setup.f90`, not in the type defaults — the type defaults to 'AUTO' for both.

### Module: `input_parser` (src/io/input_parser.f90) — TOML section rename

- **Add `parse_solver`**: reads `[solver]` with fields `method`, `mode`, `emin`, `emax`, `m0`.
- **Keep `parse_feast` temporarily**: maps `[feast]` to `solver` fields with `method='FEAST'`, `mode='ENERGY'`. Emits deprecation warning to stdout.
- **`read_config`**: tries `[solver]` first; if absent, tries `[feast]` (backward compat).
- **After all configs migrated**: delete `parse_feast` and the backward-compat path.

### Module: `simulation_setup` (src/core/simulation_setup.f90) — Solver dispatch

- **Replace `cfg%feast` with `cfg%solver`** throughout.
- **Replace `m0 < 0 → DENSE` hack** with proper string dispatch:
  - `method = 'DENSE'` → `eigen_cfg%method = 'DENSE'`
  - `method = 'FEAST'` → `eigen_cfg%method = 'FEAST'`
  - `method = 'AUTO'` → smart default per confinement (bulk→DENSE, qw→DENSE, wire→FEAST, landau→DENSE)
- **Mode dispatch**: `mode = 'AUTO'` → smart default (bulk→FULL, qw→INDEX, wire→ENERGY, landau→INDEX). Explicit modes map to `EIGEN_MODE_*` constants.
- **Add `eigen_solver` allocation for bulk/QW/Landau** (currently only wire creates it). All modes get a polymorphic solver via `make_eigensolver`.
- **Add QW CSR path**: when method is not DENSE for QW, use `ZB8bandQW_csr` and `solve_sparse`.

### Module: `hamiltonian_qw` (src/physics/hamiltonian_qw.f90) — NEW

- **`qw_workspace` type**: pre-allocated CSR blocks for 12 kp-terms (blk_Q through blk_A, blk_diff, blk_temp), COO buffers, COO-to-CSR cache, `diag_pos` arrays for Q and A diagonal positions. No scatter maps needed (1D: each kp-term has exactly one kpterms sparsity pattern; R, PP, PM are purely diagonal).
- **10 kp-term assembly functions** (`build_kp_term_Q_1d` through `build_kp_term_A_1d`): each has slow path (first call, allocates) and fast path (subsequent calls, updates workspace in-place). Formulas derived from `ZB8bandQW` dense assembly (hamiltonianConstructor.f90 lines 104–121):
  - Q = `-(kp1+kp2)*k² - kp7` (tridiag + diagonal)
  - T = `-(kp1-kp2)*k² - kp8` (tridiag + diagonal)
  - S = `2*√3*k₋ * kp9` (scalar × off-diag tridiag)
  - SC = conjg(transpose(S))
  - R = `-√3*(kp2*(kx²-ky²) - 2i*kp3*kx*ky)` (purely diagonal)
  - RC = conjg(R) (purely diagonal)
  - PZ = `(-i) * kp6` (scalar × off-diag tridiag)
  - PP = `kp4*k₊/√2` (purely diagonal)
  - PM = `kp4*k₋/√2` (purely diagonal)
  - A = `kp5 + k²*kp10` (tridiag + diagonal)
- **`ZB8bandQW_csr`**: builds QW Hamiltonian in CSR format. Calls the 10 assembly functions, then reuses COO insertion layer from `hamiltonian_wire.f90` (insert_main_blocks, insert_profile_diagonal, insert_strain_coo, insert_zeeman_coo, finalize_coo_to_csr).
- **1D kpterms→CSR**: converts the existing dense `kpterms(nz,nz,10)` tridiagonal/diagonal matrices to CSR. Uses existing `csr_apply_variable_coeff` pattern or direct tridiagonal→CSR conversion.

### Module: `hamiltonian_wire` (src/physics/hamiltonian_wire.f90) — Expose COO helpers

- **Make public**: `insert_csr_block_scaled`, `insert_main_blocks`, `insert_profile_diagonal`, `insert_strain_coo`, `insert_zeeman_coo`, `finalize_coo_to_csr`, `wire_coo_cache`, `wire_coo_cache_free`. These are grid-dimensionality-agnostic and will be reused by `hamiltonian_qw`.

### Module: `linalg` (src/math/linalg.f90) — Unguard PARDISO

- **Move `pardiso_c` interface** out of `#ifdef USE_ARPACK` block. It is needed by LDOS computation independent of ARPACK.
- **Delete all ARPACK interfaces**: `znaupd`, `zneupd`, and their public declarations.

### Module: `green_functions` (src/physics/green_functions.f90) — Unguard LDOS

- **Remove `#ifdef USE_ARPACK` guards** around `pardiso_c` import and `compute_ldos_csr`. The subroutine is always compiled.
- **Update `cfg%feast`** → `cfg%solver` references.

### All four executables — Route through polymorphic dispatch

- **main.f90**: Replace inline `zheevx` in bulk/QW/Landau k-sweeps with `solver%solve_dense`. OMP parallel k-sweeps use per-thread solver instances (thread-local `make_eigensolver`).
- **main_gfactor.f90**: Replace bulk/QW `setup_solve_kpoint_serial` LAPACK calls with solver dispatch.
- **main_optics.f90**: Replace bulk/QW optics k-sweep inline `zheevx` with solver dispatch.
- **main_topology.f90**: Replace QW BdG/Z2 direct `zheev` calls with solver dispatch. Wire BdG already uses polymorphic FEAST.

### Physics library modules — Route through polymorphic dispatch

- **sc_loop.f90**: Replace QW SC loop `zheevx` with `solver%solve_dense`. Wire SC already uses polymorphic FEAST.
- **green_functions.f90**: Replace bulk/QW spectral `zheev` with `solver%solve_dense`. Wire spectral already uses polymorphic FEAST.
- **topological_analysis.f90**: Replace Z2 TRIM loop `zheev` with `solver%solve_dense`.

### TOML config migration pattern

```
# Before:
[feast]
emin = -1.5
emax = 2.0
m0 = -1          # negative = force dense

# After:
[solver]
method = "DENSE"
mode = "ENERGY"
emin = -1.5
emax = 2.0
```

For `m0 = 0` (auto FEAST): `method = "FEAST"`, `mode = "ENERGY"`.
For `m0 > 0` (explicit FEAST subspace): `method = "FEAST"`, `mode = "ENERGY"`, `m0 = <value>`.

### Valid solver × mode combinations

| | FULL | INDEX | ENERGY |
|---|---|---|---|
| **DENSE** | ✅ `zheev` | ✅ `zheevx range='I'` | ✅ `zheevx range='V'` |
| **FEAST** | ✅ Gershgorin window | ❌ `error stop` | ✅ native |

### Smart defaults when `[solver]` absent

| Confinement | method | mode |
|---|---|---|
| bulk | DENSE | FULL |
| qw | DENSE | INDEX |
| wire | FEAST | ENERGY |
| landau | DENSE | INDEX |

## Testing Decisions

### What makes a good test

- **External behavior only**: tests verify eigenvalues and eigenvector counts, not internal workspace allocation or LAPACK routine choice.
- **Reference comparison**: dense solver results serve as the reference for sparse solver results (eigenvalues must match to numerical tolerance).
- **Error path testing**: FEAST+INDEX rejection is tested via the existing `error stop` shell-script pattern (exit code check, no output parsing).

### Modules to test

- **`eigensolver` module** — all solver × mode combinations via the polymorphic interface
- **`hamiltonian_qw` module** — CSR vs dense eigenvalue agreement for QW geometry
- **`input_parser` module** — `[solver]` parsing with valid/invalid fields
- **`definitions` module** — solver_config validation (FEAST+INDEX, emin≥emax, invalid method/mode)

### Prior art

- `tests/unit/test_eigensolver.pf`: existing FEAST/dense tests on tridiagonal Laplacian CSR — extend with mode dispatch tests
- `tests/unit/test_strain_solver.pf`: pattern for testing against known analytical values
- `tests/integration/verify_strain_wire_profile.py`: pattern for inline TOML config generation and Fortran executable invocation
- `tests/regression/`: golden-output comparison for end-to-end validation

### Test matrix

| Test | Solver | Mode | Geometry | What it verifies |
|---|---|---|---|---|
| Dense FULL | DENSE | FULL | 8×8 bulk | all 8 eigenvalues returned |
| Dense INDEX | DENSE | INDEX | NxN QW | il:iu bands returned |
| Dense ENERGY | DENSE | ENERGY | NxN QW | eigenvalues in [emin,emax] |
| FEAST ENERGY | FEAST | ENERGY | CSR wire | eigenvalues match dense |
| FEAST FULL | FEAST | FULL | CSR wire | all eigenvalues via Gershgorin |
| FEAST rejects INDEX | FEAST | INDEX | — | error stop fires |
| QW CSR vs dense | FEAST | ENERGY | QW CSR | eigenvalues match dense QW |
| Smart defaults bulk | AUTO | AUTO | bulk | equals DENSE FULL |
| Smart defaults QW | AUTO | AUTO | qw | equals DENSE INDEX |
| Smart defaults wire | AUTO | AUTO | wire | equals FEAST ENERGY |
| Smart defaults landau | AUTO | AUTO | landau | equals DENSE INDEX |

## Out of Scope

- **LDOS/PARDISO refactoring**: `compute_ldos_csr` is unguarded from ARPACK but not refactored into the eigensolver interface — it is a linear solve, not an eigenvalue problem.
- **Workspace query extraction**: The 17 duplicated LAPACK workspace query sites are not extracted into a helper function in this PRD — they are subsumed by routing through `solver%solve_dense` which manages workspace internally.
- **strain_solver decomposition** (C8 from architecture review): splitting strain_solver.f90 into bir_pikus.f90 + navier_cauchy.f90 is independent.
- **Pipeline stage assertions** (C9): adding ordering enforcement to simulation_setup is independent.
- **BdG COO workspace reuse** (C10): adding reusable COO workspace to BdG sweep is independent.
- **Topology wire setup extraction** (C5): extracting wire_setup type from main_topology is independent but facilitated by this PRD (the wire eigensolver paths become cleaner).
- **ARPACK support**: permanently removed, not merely disabled.
- **GPU offload / batched LAPACK / ELPA**: future eigensolver variants can be added to `make_eigensolver` after this PRD establishes the unified interface.

## Further Notes

### Dependency graph

```
Phase 1: Delete ARPACK
 └→ Phase 2: Add mode field to eigensolver_config
     └→ Phase 3: Unified interface (solve_dense + solve_sparse)
         ├→ Phase 4: TOML [solver] section ──┐
         └→ Phase 5: QW CSR builder ─────────┘
              └→ Phase 6: Migrate all call sites
                  └→ Phase 7: Migrate configs + docs
                      └→ Phase 8: Comprehensive tests
```

Phases 4 and 5 are independent and can proceed in parallel. Each phase must produce a compilable, testable state (113+ tests pass).

### Branch strategy

Single feature branch `refactor/eigensolver-standardization`. Each phase is a logical commit (or set of commits). Squash-merge or merge-commit at PR time — developer's choice.

### Files with `[feast]` references requiring migration

- 17 regression TOML configs in `tests/regression/configs/`
- 3 integration test Python files (inline TOML string literals)
- 4 Python scripts (config converter, lecture scripts, figure generator)
- 9 lecture docs in `docs/lecture/`
- `docs/reference/input-reference.md`, `README.md`
- `CLAUDE.md`, `src/physics/AGENTS.md`
- 10 Fortran source files with `cfg%feast` or `feast_config` references
