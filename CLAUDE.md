# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Fortran 2018 code solving the **8-band zinc-blende k.p Hamiltonian** via finite differences. Built with `-std=f2018` enforcement. Computes electronic band structures for bulk semiconductors and quantum wells, plus Landau g-factors via second-order Lowdin partitioning with commutator-based velocity operators. Computes optical absorption, gain, spontaneous emission, and intersubband transitions using commutator-based velocity matrices $v_\alpha = -i [r_\alpha, H]$. Includes self-consistent Schrödinger-Poisson solver with DIIS acceleration. GPL v3.0, authored by Tiago de Campos.

## Intent Layer

**Before modifying code in a subdirectory, read its AGENTS.md first** to understand local patterns and invariants.

- **Physics engine** (`src/physics/AGENTS.md`): 17 modules — Hamiltonian construction, optics, strain, SC loop, topology, BdG. ~113k tokens. Contains dependency DAG, basis conventions, and single-source-of-truth contracts.
- **Integration tests** (`tests/integration/AGENTS.md`): 38 files — verification ladder (rungs 1–8), standard-star benchmarks (S1–S7), convergence tests (U4–U8), coverage matrix. ~108k tokens. Shared infrastructure in `star_helpers.py` and `convergence_helpers.py`.
- **Lecture scripts** (`scripts/AGENTS.md`): 19 files — lecture-companion scripts (L00–L14), figure generation, config converter. ~155k tokens. Shared infrastructure via `tests/integration/star_helpers.py`.

### Global Invariants

- Basis ordering: bands 1–4 = valence (HH, LH, LH, HH), 5–6 = split-off, 7–8 = conduction — never change
- Unit system: Angstroms (length), eV (energy), eV·s (time)
- Single source of truth: k.p block table (`hamiltonian_blocks.f90`), Bir-Pikus formulas (`strain_solver.f90`), Zeeman table (`magnetic_field.f90`)
- All config validation via `validate()` + `validate_semantic()` in `defs.f90` — no silent corrections

## Build Commands

```bash
# Configure (first time or after clean)
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl

# Build
cmake --build build                # builds all four executables into build/src/

# Or use the Make wrapper
make all            # Configure + build all executables
make run            # Build and run bandStructure
make clean          # Remove build/ directory
make clean_all      # Also remove output files
make test           # Configure with BUILD_TESTING=ON, build, run ctest
```

**Executables** at `build/src/bandStructure`, `build/src/gfactorCalculation`, `build/src/opticalProperties`, and `build/src/topologicalAnalysis`.

**Prerequisites:** gfortran (version supporting F2018), Intel MKL (sequential, LP64 interface), FFTW3, CMake >= 3.15, Ninja (optional). Compiler enforced to `-std=f2018` via `CMAKE_Fortran_FLAGS`. MKL defaults: `MKL_INTERFACE=lp64`, `MKL_THREADING=sequential` (set in root `CMakeLists.txt`).

**Gotcha — stale `.mod` files:** Old `.mod` in project root shadow fresh ones in `build/`. Run `rm -f *.mod` if you get type mismatch errors.

**Gotcha — `cp -i` alias:** The shell has `cp -i` alias. Use the Write tool (not `cp`) when creating `input.toml`.

## Testing

```bash
# Configure with tests enabled (pFUnit installed at $HOME/.local/pfunit/PFUNIT-<ver>)
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl \
    -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-4.16/cmake
cmake --build build

# Run tests
ctest --test-dir build                    # all tests (113: 35 unit + 44 regression + 18 verification + 9 standard-star (incl. 2 dual-labeled) + 6 convergence + 2 strain + 1 coverage + misc)
ctest --test-dir build -j4                # parallel: 4 test jobs concurrently (cuts ~21min to ~8min)
ctest --test-dir build -L unit            # pFUnit unit tests only
ctest --test-dir build -L regression      # regression/golden-output tests only
ctest --test-dir build -L verification    # 8-band verification ladder (4 rungs)
ctest --test-dir build -L standard-star   # standard-star physics benchmarks (S1-S7)
ctest --test-dir build -L strain-validation  # strain validation (bulk/QW/wire InAs/GaAs)
ctest --test-dir build -L coverage        # validation coverage matrix
ctest --test-dir build -V                 # verbose output

# With OpenMP threading for QW k-point sweeps:
OMP_NUM_THREADS=12 ctest --test-dir build -j4 --output-on-failure
```

pFUnit installs as **uppercase `PFUNIT`** — use `-DPFUNIT_DIR=.../PFUNIT-<ver>/cmake`. Built from source: `git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit`.

Regression tests use shell scripts + Python (`tests/regression/compare_output.py`) comparing numerical output against reference data. Validation rejection tests (`tests/integration/test_validate_rejects_bad_configs.sh`) verify `error stop` branches via shell exit-code checks. Python dependencies: PyYAML for coverage matrix tool; numpy for standard-star/verification scripts.

## Running

```bash
./build/src/bandStructure       # Reads input.toml, outputs to output/
./build/src/gfactorCalculation  # Reads input.toml, outputs to output/ (requires k=0)
./build/src/opticalProperties   # Reads input.toml with [optics] section, outputs to output/
./build/src/topologicalAnalysis # Reads input.toml with [topology] section, outputs to output/
```

**Parallelism:** OpenMP distributes the QW k-point sweep across threads (`main.f90:691-723`). MKL runs sequential internally to avoid oversubscription. Control thread count with `OMP_NUM_THREADS`:

```bash
OMP_NUM_THREADS=12 ./build/src/bandStructure   # use 12 of 24 cores
OMP_NUM_THREADS=$[$(nproc)/2] ./build/src/bandStructure  # half of available cores
```

Without `OMP_NUM_THREADS`, OpenMP defaults to all cores. For test runs, set it before `ctest`:

```bash
OMP_NUM_THREADS=12 ctest --test-dir build --output-on-failure
```

Write a config from `tests/regression/configs/` to `input.toml`. Configs use TOML format with sections organized as documented in `docs/reference/input-reference.md`. Sections are order-independent; optional physics blocks are enabled by their presence in the file (no separate enable flags).

**Lecture-test pair scripts:** 15 Python scripts in `scripts/` (`lecture_00_quickstart.py` through `lecture_14_excitons_scattering.py`) serve dual roles as pedagogical companions and integration tests. Each runs a Fortran executable, validates physics results, generates overlay plots, and prints PASS/FAIL per section. Run individually: `python3 scripts/lecture_01_bulk.py`. Shared infrastructure: `tests/integration/star_helpers.py` (physical constants, `run_exe`, parse helpers).

**Validation Coverage Matrix:** The `coverage` ctest label runs a cross-reference between declared physics observables and test annotations. Run with `ctest --test-dir build -L coverage`. The universe file `tests/integration/validation_universe.yml` declares cells with `(observable, geometry, material, tier)`. Test files carry annotations in the format:

```
# COVERAGE: observable=Eg geometry=bulk material=GaAs ref=Vurgaftman2001
# COVERAGE: observable=m*_e geometry=bulk material=InAs
```

Convention: new tests should include `# COVERAGE:` annotations so the matrix can track physics coverage. Shell scripts delegating to `verify_*.py` scripts should place annotations in the verifier, not the wrapper. Test dependency: PyYAML (`pip install pyyaml`) is required for the coverage matrix tool.

## Cross-Code Validation

```bash
# Run all 12 validation tests (6 cross-code against kdotpy + 6 single-code/analytical)
source validation/kdotpy_env/bin/activate
python3 validation/run_all.py

# Run individual tests
python3 validation/bulk/test_bulk_k0.py
python3 validation/strain/test_strain_bandedge.py
```

The `validation/` directory contains 12 validation tests. Six are cross-code tests comparing our Fortran solver against kdotpy (Python, plane-wave discretization): bulk k=0 gate, bulk dispersion, QW subbands, QW dispersion, QW convergence, and strain bandedge. Six are single-code/analytical tests (kdotpy lacks a comparable API for these physics): bulk Zeeman (analytical Roth formula), Landau bulk (analytical LL formula), wire subbands (single-code consistency), g-factor QW (analytical Roth formula), strain QW (analytical Bir-Pikus), and self-consistent QW (single-code convergence). Strain tests cover both compressive (InAs/GaAs) and tensile (GaAs/InP) regimes. Results saved to `validation/*/results/`. Shared infrastructure: `validation/shared/` (parameter mapper, Fortran/kdotpy runners, comparison utilities). Parameter mapping between Fortran and kdotpy conventions in `validation/shared/param_mapper.py`.

## Architecture

### Directory structure

```
src/
  core/       defs.f90, parameters.f90, simulation_setup.f90, utils.f90
  math/       finitedifferences.f90, linalg.f90, sparse_matrices.f90, eigensolver.f90, geometry.f90
  io/         outputFunctions.f90, input_parser.f90
  physics/    hamiltonian_blocks.f90, hamiltonianConstructor.f90, confinement_init.f90, hamiltonian_wire.f90, gfactor_functions.f90, optical_spectra.f90, spin_projection.f90, poisson.f90, charge_density.f90, sc_loop.f90, exciton.f90, strain_solver.f90, scattering.f90, magnetic_field.f90, topological_analysis.f90, bdg_hamiltonian.f90, green_functions.f90
  apps/       main.f90, main_gfactor.f90, main_optics.f90, main_topology.f90
tests/
  unit/       pFUnit .pf test files (34 tests: defs, FD, utils, parameters, Hamiltonian, Hamiltonian blocks, CSR, eigensolver, Poisson, charge density, SC loop, geometry, optical, topology, magnetic field, BdG, Landau, Z2, strain_solver, green_functions)
  support/    shared test-support library (CSR helpers, Krylov infrastructure, reference data)
  integration/  shell scripts + Python verification scripts for full-executable tests
  regression/   configs/, data/, compare_output.py
validation/   validation pipeline (12 tests: 6 cross-code against kdotpy + 6 analytical)
cmake/        FindFFTW3.cmake
build/        .o, .mod, executables (created by cmake)
scripts/      lecture_*.py executable lecture-companion scripts (L00-L14) + generate_all_figures.py
docs/solutions/  documented solutions to past problems (bugs, best practices, patterns), organized by category with YAML frontmatter (module, tags, problem_type)
```

### Four executables

| Executable | Entry Point | Purpose |
|---|---|---|
| `bandStructure` | `src/apps/main.f90` (`program kpfdm`) | Band structure vs. k-vector sweep |
| `gfactorCalculation` | `src/apps/main_gfactor.f90` (`program gfactor`) | g-factor at Gamma point (k=0) |
| `opticalProperties` | `src/apps/main_optics.f90` (`program opticalProperties`) | Optical spectra (absorption, gain, spontaneous emission, ISBT) for bulk/QW/wire |
| `topologicalAnalysis` | `src/apps/main_topology.f90` (`program topologicalAnalysis`) | Topological invariants: Chern number, Z2, Majorana modes |

### Module dependency graph

```
defs.f90                      (kinds, constants, derived types — no deps)
  <- parameters.f90           (material parameter database for 25+ semiconductors)
       <- utils.f90           (dense-to-sparse conversion, Simpson integration, file utilities)
            <- finitedifferences.f90  (FD stencils/orders 2-10, Toeplitz matrices, Vandermonde derivative and interpolation solvers)
                 <- hamiltonian_blocks.f90     (8x8 k.p block structure as data: 52-entry table, named constants KP_Q..KP_A)
                 <- magnetic_field.f90         (Zeeman table SSOT + splitting + Peierls phase as COO insertions)
                 <- hamiltonianConstructor.f90  (8x8 bulk & 8NxN QW Hamiltonian, commutator velocity matrices; reads kp/strain tables, Zeeman from magnetic_field)
                      <- hamiltonian_wire.f90     (2D wire Hamiltonian, CSR assembly, wire geometry; reads kp block table)
                      <- gfactor_functions.f90  (Lowdin partitioning, spin/momentum matrix elements)
                      <- optical_spectra.f90    (absorption, gain, spontaneous emission, ISBT accumulation)
                      <- spin_projection.f90     (Clebsch-Gordan spin-up/down decomposition)
                 <- poisson.f90           (box-integration Poisson solver, Thomas algorithm)
                 <- charge_density.f90    (n(z), p(z) from k.p eigenstates + k_∥ sampling)
                      <- sc_loop.f90           (self-consistent SP loop, linear + DIIS mixing)
  <- scattering.f90          (phonon scattering rates, used by main.f90)
  <- outputFunctions.f90      (eigenvalue/eigenfunction file I/O)
  <- input_parser.f90         (TOML-based input.toml parser → simulation_config type, uses toml-f library)
  <- simulation_setup.f90     (simulation orchestration: H-build, eigensolve, k-point sweep for all 4 confinement modes; used by all four executables)
```

### Key design concepts

- **Basis ordering** (fixed throughout): bands 1-4 = valence (HH, LH, LH, HH), bands 5-6 = split-off, bands 7-8 = conduction
- **Dual mode**: `confinement = "bulk"` → bulk (8x8), `confinement = "qw"` → quantum well (8N x 8N, N = grid%npoints())
- **QW Hamiltonian**: 8x8 block matrix, each block NxN. Block structure defined in `hamiltonian_blocks.f90` as a 52-entry table (`get_kp_block_table()`) with named constants (KP_Q, KP_R, etc.). Both dense (`hamiltonianConstructor.f90`) and COO (`hamiltonian_wire.f90`) builders read this table instead of hard-coding the block topology. Blocks: k.p terms (Q, R, S, T, P) + band offsets
- **`simulation_config`** derived type in `defs.f90` holds all parsed input parameters; `input_parser.f90` populates it from `input.toml` using sub-types mirroring TOML sections (`wave_vector_config`, `bands_config`, `external_field_config`, `b_field_config`, `wire_config`, `landau_config`, `solver_config`)
- **`conf_direction()`** function returns the confinement direction character: bulk → `'n'` (none), QW → `'z'`, wire → `'z'`, Landau → `'x'`. Replaces the old `conf_dir` stored field.
- **`grid%npoints()`** accessor returns the total number of spatial grid points: bulk → 1, QW → `fd_step`, wire → `nx * ny`, Landau → `landau.nx`. Replaces the old `ngrid` computed field.
- **`num_layers`** counts `[[material]]` entries (material layers for bulk/QW). Wire mode uses `cfg%wire%num_regions` for `[[region]]` entries instead.
- **`confinementInitialization`** precomputes material parameters at each z-point into `kpterms(ny, ny, 10)`
- **Last-layer-wins** in `confinement_init.f90`: later material layers overwrite earlier ones. Use 2-layer pattern (barrier covers full domain, well overwrites center). Avoid full-domain layers after a narrow well layer — they silently destroy the QW.
- **Sparse vs dense**: g-factor uses MKL SpBLAS for large QW; band structure uses dense LAPACK (`zheevx`)
- **Wire g-factor velocity operator**: commutator-based `build_velocity_matrices` computes $v_\alpha = -i [r_\alpha, H]$ element-wise on CSR. Generic interface dispatches to `_2d` (wire, 4 scalar CSR args) or `_1d` (QW, 3-element CSR array). For QW: z-velocity from commutator, in-plane from `ZB8bandQW(g='g')`. g='g3' still used for z-direction in wire. Old g='g1'/'g2' modes are dead code.
- **Optical properties**: standalone `opticalProperties` executable computes absorption, gain, spontaneous emission, and ISBT spectra for bulk (8x8, spherical 3D k-integration), QW (dense, 2D cylindrical), and wire (CSR/FEAST, 1D kz-sweep). All use commutator-based velocity matrices via CSR SpMV + zdotc pattern. Unified accumulation framework in `optical_spectra.f90` shares velocity matrices across all quantities; only the occupation factor differs. Spin-resolved spectra via Clebsch-Gordan projection in `spin_projection.f90`.
- **Foreman renormalization**: disabled by default (`renormalization = .False.` in defs.f90)
- **Topological analysis**: standalone `topologicalAnalysis` executable supports QHE (Chern number via QWZ model), QSHE (Z2 invariant via BHZ wire or Fu-Kane), and BdG (Majorana modes in superconducting wire). Uses `[topology]` and `[bdg]` TOML sections. Chern number via Fukui-Hatsugai-Suzuki algorithm; Z2 via gap criterion in 1D or parity method in 2D. BdG builds 16N x 16N Nambu-space Hamiltonian with s-wave pairing and Zeeman splitting. LDOS computed via complex PARDISO. Design doc: `docs/plans/archive/2026-04-27-bdg-topological-superconductivity-design.md`
- **W-variant materials** (GaAsW, InAsW, etc.): use Winkler's parameter set with InSb as EV reference. Non-W materials use Vurgaftman parameters. EC = EV + Eg convention for all materials with EV defined
- **Self-consistent SP**: iterative Schrödinger-Poisson loop wraps the eigenvalue solve. Modifies `profile` array before Hamiltonian construction (same interface as `externalFieldSetup_electricField`). Linear mixing warm-up + DIIS/Pulay acceleration. Works for bulk (Fermi level finder) and QW (full SP). Per-layer doping (`doping_spec`), charge neutrality or fixed Fermi level, box-integration Poisson with variable dielectric. Design doc: `docs/plans/archive/2026-03-29-self-consistent-sp-design.md`
- **Config validation** consolidated in `validate()` (19 structural + 8 new checks) and `validate_semantic(cfg, app_name)` (5 existing + 10 new app-specific checks) in `defs.f90`. Executables only handle runtime errors (LAPACK info codes, FEAST convergence, file I/O). All config-level checks use `error stop` with contextual messages. No silent corrections. Parser guard for `fd_step < 2` prevents division-by-zero before `validate()` runs. See ADR 0002 (`docs/adr/0002-consolidate-validation.md`) for the full checklist.

### Input file (`input.toml`)

TOML format parsed by `toml-f` library. Sections are order-independent. Key top-level fields: `confinement` (`"bulk"`, `"qw"`, `"wire"`, `"landau"`), `FDorder` (default 2), `fd_step` (grid points for QW; not used for wire/Landau, which derive the grid from `[wire]`/`[landau]` sections).

Sections: `[wave_vector]` (mode, max, nsteps), `[bands]` (num_cb, num_vb), `[[material]]` (name, z_min, z_max), `[wire]` + `[wire.geometry]` + `[[region]]` (wire mode), `[landau]` (Landau mode), `[external_field]` (type, value), `[b_field]` (components, g_factor), `[strain]` (reference), `[sc]` (self-consistent parameters), `[[doping]]` (ND, NA or delta doping), `[topology]` (topological analysis), `[bdg]` (BdG parameters), `[optics]` (optical spectra), `[exciton]`, `[scattering]`, `[solver]` (method, mode, emin, emax, m0).

Optional sections are enabled by presence -- no separate enable flags. G-factor uses top-level `which_band` and `band_idx`. See `docs/reference/input-reference.md` for the complete schema.

## Code Conventions

- **Fortran standard:** F2018 enforced via `-std=f2018` in CMake. No GNU extensions.
- Prefer generic intrinsics (`sqrt`) over legacy typed intrinsics (`dsqrt`, `dble`, etc.).
- Prefer `do` / `do concurrent` over `forall` (removed in F2008, restored in F2018).
- Prefer `execute_command_line` over non-standard `call system(...)`.
- Use `c_loc()` from `iso_c_binding` instead of non-standard `loc()`.
- All modules use `private` default with explicit `public` exports. When adding new modules, use `private` default and enumerate `public ::` exports.
- When adding new scalar `pure` functions, use `elemental pure` by default (F2008 requires both keywords; F2018+ implies pure).
  - **Known exceptions:** `grid_ngrid` (`defs.f90:978`), `to_lower_ascii` (`input_parser.f90:17`) — upgrade when touching those functions
- Declaration ordering: variables must be declared before use in array dimension expressions.
- No `goto` in new code — use named `do` loops with `exit` for early-return blocks.
- Use `error stop` (not `stop 1`) for all fatal error exits. Include a descriptive message string. `stop 1` without message is deprecated across the codebase.
- External BLAS/LAPACK/PARDISO declarations go through `linalg.f90` interfaces, not local `external ::`.
- All types with allocatable components have finalizers (delegating to `*_free` routines where they exist, inlined where cross-module delegation isn't possible). Keep explicit `*_free` routines public for manual control. Types with both manual `*_free` and finalizer use a `was_freed` flag for idempotent double-free protection: `*_free` returns immediately if already called.
- `do concurrent` used on proven-independent loops: velocity matrix construction, optics finalization, kpterms diagonal init.
- `csr_matrix` has type-bound `free()` and `clone_structure()` but components remain public for hot-path access.
- `contiguous` attribute on all assumed-shape hot-path array arguments (not on optional or allocatable).
  - **Known gaps** (will be fixed when touching those routines): `dns(:,:)` in `utils.f90:19`, `psi(:)` in `spin_projection.f90:14`
- `spatial_grid%npoints()` accessor preferred over raw `grid%nx * grid%ny` for total grid size.
- `iso_c_binding` used for all MKL C APIs: PARDISO as `pardiso_c`, FEAST wrappers, `mkl_set_num_threads_local`. PARDISO/FEAST scalars passed by reference (MKL C API passes pointers). `mkl_set_num_threads_local` uses `value` since the C function takes `int` by value.
- PARDISO has two `iso_c_binding` interfaces: `pardiso_c` (complex, used by topology/LDOS eigensolver) and `pardiso_real` (real-valued, used by Poisson solver). Both use `c_intptr_t` for the `pt` handle array.
- Polymorphic eigensolver: `make_eigensolver(config)` factory returns a polymorphic solver; dispatch via `solver%solve(...)`, `solver%solve_sparse(...)`, or `solver%solve_dense(...)`. Three mode constants: `EIGEN_MODE_FULL` (all eigenvalues), `EIGEN_MODE_INDEX` (range il:iu), `EIGEN_MODE_ENERGY` (range [emin, emax]). `eigensolver_config_validate` rejects invalid combos (e.g., FEAST+INDEX). Existing `solve_sparse_evp` remains available as a direct interface.
- `fortran-stdlib` as optional dependency: CMake uses `find_package(fortran_stdlib CONFIG QUIET)`; code compiles without it.
- `fpm.toml` exists as experimental alternative build system (documented as such in README). CMake/Ninja remains the primary build.
- `g='g3'` derivative builds must stay isolated from `wire_workspace` cache.
- `feast_workspace` reuse must remain pattern-validated.
- `zdotc` declared as function-return ABI in `linalg.f90` — MKL on this platform returns `complex(dp)` by value, not via subroutine hidden argument. Do not change to subroutine form.
- Eigensolver variants guarded by preprocessor: `#ifdef USE_MKL_FEAST` for FEAST solver. `make_eigensolver(config)` factory dispatches on string; add new variants behind their own `#ifdef` guard.
- Scalar `pure` functions upgraded to `elemental pure`: `kronij`, `fermi_dirac`, `flat_idx`, `wire_flat_idx`, `compute_bp_scalar`, `segment_circle_fraction`.
- Precision kinds in `defs.f90`: `sp` (real32), `dp` (real64), `qp` (real128), `iknd` (int32) — aliased from `iso_fortran_env`
- k-vectors and wavevector sweeps use `real(dp)` throughout
- File/function length guidelines: 300 lines/file, 50 lines/function
- Prioritize: Correctness > Maintainability > Readability > Type-Exactness > Performance > Minimal

## Knowledge Store

`docs/solutions/` contains documented solutions to past problems, organized by category (`logic-errors/`, `best-practices/`, `workflow/`, `patterns/`). Each doc has YAML frontmatter with `module`, `tags`, `problem_type`, and `component` fields for searchability. Search this directory before implementing features, debugging issues, or modifying FD operators — learnings cover bugs, convergence testing methodology, and workflow patterns.

## Git Workflow

- **Branch naming:** `feature/<description>`, `feat/<description>`, `fix/<description>`
- **Commits:** Conventional Commits (`feat:`, `fix:`, `refactor:`, `test:`, `docs:`)
- **PRs:** All tests must pass before merge. Use `gh pr create` from feature branches.

## Superpowers Skills

Always check for and follow applicable superpowers skills when working. In particular:
- **brainstorming** — for design decisions and scoping before implementation
- **writing-plans** — for creating implementation plans before coding
- **subagent-driven-development** — for executing plans task-by-task with review
- **test-driven-development** — subagents follow TDD for each task

## Boundaries

- **NEVER** modify material parameters in `parameters.f90` without verifying against published references (Vurgaftman 2001, Winkler 2003)
- **NEVER** change the basis ordering (bands 1-4 valence, 5-6 split-off, 7-8 conduction) — it is hardcoded throughout
- **NEVER** change the Bir-Pikus sign convention: `Q_eps = -(b_dp/2) * (eps_zz - 0.5*(eps_yy + eps_xx))` with the minus sign on b_dp matching the standard Chuang/Winkler convention. VB diagonal shifts are `delta_EHH = -P_eps + Q_eps`, `delta_ELH = -P_eps - Q_eps`, `delta_ESO = -P_eps` where `P_eps = -av * Tr(eps)`. Under compressive strain (Tr < 0, b < 0), HH shifts up and LH shifts down. Single source of truth: `compute_bp_scalar` in `strain_solver.f90`.
- **NEVER** change the strain or Zeeman block tables: `get_strain_table()` in `strain_solver.f90` is the single source of truth for Bir-Pikus strain block topology (band pairs). `get_zeeman_table()` in `magnetic_field.f90` is the single source of truth for Zeeman diagonal g-multipliers (per band index). All dense and COO builders consume these tables.
- **NEVER** change the k.p block table: `get_kp_block_table()` in `hamiltonian_blocks.f90` is the single source of truth for the 52-entry block topology (band pairs, k.p terms, complex prefactors). Both dense and COO builders consume this table.
- **NEVER** commit `input.toml` with personal test configs — use `tests/regression/configs/` for test configs
- **Require approval** for: changes to `defs.f90` derived types, k.p block table in `hamiltonian_blocks.f90`, Hamiltonian construction in `hamiltonianConstructor.f90`, FD stencil coefficients in `finitedifferences.f90`, Poisson solver in `poisson.f90`, SC loop convergence logic in `sc_loop.f90`

## Known Issues

- `bir_pikus_blocks` in `defs.f90:83` has no finalizer. This is intentional: the type contains only non-allocatable fixed-size arrays (no pointers or allocatable components), so no explicit finalizer is needed. Do not add one.
- The wire g-factor was fixed using the commutator-based velocity operator (`build_velocity_matrices` in `hamiltonianConstructor.f90`). The transverse perturbation construction uses $-i [r_\alpha, H]$ element-wise on the CSR Hamiltonian, which correctly captures all k.p term contributions.
