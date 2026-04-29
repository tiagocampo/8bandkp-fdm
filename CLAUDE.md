# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Fortran 2008 code solving the **8-band zinc-blende k.p Hamiltonian** via finite differences. Built with `-std=f2008` enforcement. Computes electronic band structures for bulk semiconductors and quantum wells, plus Landau g-factors via second-order Lowdin partitioning with commutator-based velocity operators. Computes optical absorption, gain, spontaneous emission, and intersubband transitions using commutator-based velocity matrices $v_\alpha = -i [r_\alpha, H]$. Includes self-consistent Schrödinger-Poisson solver with DIIS acceleration. GPL v3.0, authored by Tiago de Campos.

## Build Commands

```bash
# Configure (first time or after clean)
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl

# Build
cmake --build build                # builds all three executables into build/src/

# Or use the Make wrapper
make all            # Configure + build all executables
make run            # Build and run bandStructure
make clean          # Remove build/ directory
make clean_all      # Also remove output files
make test           # Configure with BUILD_TESTING=ON, build, run ctest
```

**Executables** at `build/src/bandStructure`, `build/src/gfactorCalculation`, and `build/src/opticalProperties`.

**Prerequisites:** gfortran, Intel MKL (sequential, LP64 interface), FFTW3, CMake >= 3.15, Ninja (optional). Compiler enforced to `-std=f2008` via `CMAKE_Fortran_FLAGS`. MKL defaults: `MKL_INTERFACE=lp64`, `MKL_THREADING=sequential` (set in root `CMakeLists.txt`).

**Gotcha — stale `.mod` files:** Old `.mod` in project root shadow fresh ones in `build/`. Run `rm -f *.mod` if you get type mismatch errors.

**Gotcha — `cp -i` alias:** The shell has `cp -i` alias. Use the Write tool (not `cp`) when creating `input.cfg`.

## Testing

```bash
# Configure with tests enabled (pFUnit installed at $HOME/.local/pfunit/PFUNIT-<ver>)
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl \
    -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-4.16/cmake
cmake --build build

# Run tests
ctest --test-dir build                    # all tests (15 unit + 24 regression)
ctest --test-dir build -L unit            # pFUnit unit tests only
ctest --test-dir build -L regression      # regression/golden-output tests only
ctest --test-dir build -V                 # verbose output
```

pFUnit installs as **uppercase `PFUNIT`** — use `-DPFUNIT_DIR=.../PFUNIT-<ver>/cmake`. Built from source: `git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit`.

Regression tests use shell scripts + Python (`tests/regression/compare_output.py`) comparing numerical output against reference data. No extra dependencies beyond Python 3.

## Running

```bash
./build/src/bandStructure       # Reads input.cfg, outputs to output/
./build/src/gfactorCalculation  # Reads input.cfg, outputs to output/ (requires k=0)
./build/src/opticalProperties   # Reads input.cfg with optics: block, outputs to output/
```

Write a config from `tests/regression/configs/` to `input.cfg`. Keep committed examples in the canonical order documented in `docs/reference/input-reference.md`; optional block entry labels are name-aware (`optics:`, `exciton:`, `scattering:`, `feast_emin:`, `strain:`), but parameters inside each block still follow the documented sequence.

## Architecture

### Directory structure

```
src/
  core/       defs.f90, parameters.f90, utils.f90
  math/       mkl_spblas.f90, mkl_sparse_handle.f90, finitedifferences.f90, linalg.f90, sparse_matrices.f90, eigensolver.f90, geometry.f90
  io/         outputFunctions.f90, input_parser.f90
  physics/    hamiltonianConstructor.f90, confinement_init.f90, hamiltonian_wire.f90, gfactor_functions.f90, optical_spectra.f90, spin_projection.f90, poisson.f90, charge_density.f90, sc_loop.f90, exciton.f90, strain_solver.f90
  apps/       main.f90, main_gfactor.f90, main_optics.f90
tests/
  unit/       pFUnit .pf test files (test_defs, test_finitedifferences, test_utils, test_parameters, test_hamiltonian, test_hamiltonian_2d, test_csr_spmv, test_eigensolver, test_poisson, test_charge_density, test_sc_loop)
  integration/  shell scripts for full-executable tests
  regression/   configs/, data/, compare_output.py
cmake/        FindFFTW3.cmake
build/        .o, .mod, executables (created by cmake)
scripts/      gnuplot plotting scripts
```

### Three executables

| Executable | Entry Point | Purpose |
|---|---|---|
| `bandStructure` | `src/apps/main.f90` (`program kpfdm`) | Band structure vs. k-vector sweep |
| `gfactorCalculation` | `src/apps/main_gfactor.f90` (`program gfactor`) | g-factor at Gamma point (k=0) |
| `opticalProperties` | `src/apps/main_optics.f90` (`program opticalProperties`) | Optical spectra (absorption, gain, spontaneous emission, ISBT) for bulk/QW/wire |

### Module dependency graph

```
defs.f90                      (kinds, constants, derived types — no deps)
  <- parameters.f90           (material parameter database for 25+ semiconductors)
  <- mkl_spblas.f90           (vendor: Intel MKL sparse BLAS interface)
       <- utils.f90           (dense-to-sparse conversion, Simpson integration)
            <- finitedifferences.f90  (FD stencils/orders 2-10, Toeplitz matrices)
                 <- hamiltonianConstructor.f90  (8x8 bulk & 8NxN QW Hamiltonian, commutator velocity matrices)
                      <- gfactor_functions.f90  (Lowdin partitioning, spin/momentum matrix elements)
                      <- optical_spectra.f90    (absorption, gain, spontaneous emission, ISBT accumulation)
                      <- spin_projection.f90     (Clebsch-Gordan spin-up/down decomposition)
                 <- poisson.f90           (box-integration Poisson solver, Thomas algorithm)
                 <- charge_density.f90    (n(z), p(z) from k.p eigenstates + k_∥ sampling)
                      <- sc_loop.f90           (self-consistent SP loop, linear + DIIS mixing)
  <- outputFunctions.f90      (eigenvalue/eigenfunction file I/O)
  <- input_parser.f90         (shared input.cfg parser → simulation_config type)
```

### Key design concepts

- **Basis ordering** (fixed throughout): bands 1-4 = valence (HH, LH, LH, SO), bands 5-6 = split-off, bands 7-8 = conduction
- **Dual mode**: `confinement=0` → bulk (8x8), `confinement=1` → quantum well (8N x 8N, N = FDstep)
- **QW Hamiltonian**: 8x8 block matrix, each block NxN. Blocks: k.p terms (Q, R, S, T, P) + band offsets
- **`simulation_config`** derived type in `defs.f90` holds all parsed input parameters; `input_parser.f90` populates it
- **`confinementInitialization`** precomputes material parameters at each z-point into `kpterms(fdStep, fdStep, 10)`
- **Sparse vs dense**: g-factor uses MKL SpBLAS for large QW; band structure uses dense LAPACK (`zheevx`)
- **Wire g-factor velocity operator**: commutator-based `build_velocity_matrices` computes $v_\alpha = -i [r_\alpha, H]$ element-wise on CSR. Generic interface dispatches to `_2d` (wire, 4 scalar CSR args) or `_1d` (QW, 3-element CSR array). For QW: z-velocity from commutator, in-plane from `ZB8bandQW(g='g')`. g='g3' still used for z-direction in wire. Old g='g1'/'g2' modes are dead code.
- **Optical properties**: standalone `opticalProperties` executable computes absorption, gain, spontaneous emission, and ISBT spectra for bulk (8x8, spherical 3D k-integration), QW (dense, 2D cylindrical), and wire (CSR/FEAST, 1D kz-sweep). All use commutator-based velocity matrices via CSR SpMV + zdotc pattern. Unified accumulation framework in `optical_spectra.f90` shares velocity matrices across all quantities; only the occupation factor differs. Spin-resolved spectra via Clebsch-Gordan projection in `spin_projection.f90`.
- **Foreman renormalization**: disabled by default (`renormalization = .False.` in defs.f90)
- **W-variant materials** (GaAsW, InAsW, etc.): use Winkler's parameter set with InSb as EV reference. Non-W materials use Vurgaftman parameters. EC = EV + Eg convention for all materials with EV defined
- **Self-consistent SP**: iterative Schrödinger-Poisson loop wraps the eigenvalue solve. Modifies `profile` array before Hamiltonian construction (same interface as `externalFieldSetup_electricField`). Linear mixing warm-up + DIIS/Pulay acceleration. Works for bulk (Fermi level finder) and QW (full SP). Per-layer doping (`doping_spec`), charge neutrality or fixed Fermi level, box-integration Poisson with variable dielectric. Design doc: `docs/plans/2026-03-29-self-consistent-sp-design.md`

### Input file (`input.cfg`)

Label-value format with comments (`!`). Key fields: `waveVector` (kx/ky/kz/k0), `waveVectorMax`, `waveVectorStep`, `confinement` (0/1), `FDstep`, `FDorder` (default 2), `numLayers`, `materialN` (name + z-range), `numcb`/`numvb`, `ExternalField` + `EFParams`.

SC additional fields: `SC` (enable flag), `SC_max_iter`, `SC_tolerance`, `SC_mixing`, `SC_diis`, `SC_temperature`, `SC_fermi_mode` (0=charge_neutrality, 1=fixed), `SC_fermi_level`, `SC_num_kpar`, `SC_kpar_max`, `SC_bc` (DD/DN), `SC_bc_left`, `SC_bc_right`. Per-layer `doping: ND NA` lines.

g-factor additional fields: `whichBand` (0=CB, 1=VB) and `bandIdx` (subband index, default 1).

Optics block fields: `optics:` block with `T/F` enable flag, `linewidth_lorentzian`, `linewidth_gaussian`, `refractive_index`, `E_min E_max num_energy_points`, `temperature`, `carrier_density`, `gain_enabled`, `gain_carrier_density`, `ISBT`, `SpontaneousEnabled`, `SpinResolved`. See `docs/reference/input-reference.md` for full format.

## Code Conventions

- **Fortran standard:** F2008 enforced via `-std=f2008` in CMake. No GNU extensions.
- Prefer generic intrinsics (`sqrt`) over legacy typed intrinsics (`dsqrt`, `dble`, etc.).
- Prefer `do` / `do concurrent` over `forall` (deprecated in F2008).
- Prefer `execute_command_line` over non-standard `call system(...)`.
- Use `c_loc()` from `iso_c_binding` instead of non-standard `loc()`.
- All modules use `private` default with explicit `public` exports. When adding new modules, use `private` default and enumerate `public ::` exports.
- When adding new scalar `pure` functions, use `elemental pure` by default (F2008 requires both keywords; F2018+ implies pure).
- Declaration ordering: variables must be declared before use in array dimension expressions.
- No `goto` in new code — use named `do` loops with `exit` for early-return blocks.
- External BLAS/LAPACK/PARDISO declarations go through `linalg.f90` interfaces, not local `external ::`.
- All types with allocatable components have finalizers (delegating to `*_free` routines where they exist, inlined where cross-module delegation isn't possible). Keep explicit `*_free` routines public for manual control.
- `do concurrent` used on proven-independent loops: velocity matrix construction, optics finalization, kpterms diagonal init.
- `csr_matrix` has type-bound `free()` and `clone_structure()` but components remain public for hot-path access.
- `contiguous` attribute on all assumed-shape hot-path array arguments (not on optional or allocatable).
- `spatial_grid%npoints()` accessor preferred over raw `grid%nx * grid%ny` for total grid size.
- `iso_c_binding` used for all MKL C APIs: PARDISO as `pardiso_c`, FEAST wrappers, `mkl_set_num_threads_local`. Scalars passed by reference (not `value`) since MKL C API passes everything as pointers.
- Polymorphic eigensolver: `make_eigensolver(config)` factory returns a polymorphic solver; dispatch via `solver%solve(...)`. Existing `solve_sparse_evp` remains available as a direct interface.
- `fortran-stdlib` as optional dependency: CMake uses `find_package(fortran_stdlib CONFIG QUIET)`; code compiles without it.
- `g='g3'` derivative builds must stay isolated from `wire_workspace` cache.
- `feast_workspace` reuse must remain pattern-validated.
- `simulation_config` is validated after parsing via the type-bound `config%validate()` method — use `error stop` for invalid configs.
- Scalar `pure` functions upgraded to `elemental pure`: `kronij`, `fermi_dirac`, `flat_idx`, `wire_flat_idx`, `compute_bp_scalar`, `segment_circle_fraction`.
- Precision kinds in `defs.f90`: `sp` (single), `dp` (double), `qp` (quad), `iknd` (int64)
- k-vectors and wavevector sweeps use `real(dp)` throughout
- File/function length guidelines from `.clinerules`: 300 lines/file, 50 lines/function
- Prioritize: Correctness > Maintainability > Readability > Type-Exactness > Performance > Minimal

## Git Workflow

- **Branch naming:** `feature/<description>`, `fix/<description>`
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
- **NEVER** change the Bir-Pikus sign convention: VB diagonal shifts all contain `-P_eps` (which equals `av * Tr(eps)` due to the sign flip in the helper), so delta_EHH = -P_eps + Q_eps, delta_ELH = -P_eps - Q_eps, delta_ESO = -P_eps. The helper `P_eps = -av * Tr(eps)` introduces a sign flip — every consumer must negate it. Single source of truth: `compute_bp_scalar` in `strain_solver.f90`.
- **NEVER** commit `input.cfg` with personal test configs — use `tests/regression/configs/` for test configs
- **Require approval** for: changes to `defs.f90` derived types, Hamiltonian construction in `hamiltonianConstructor.f90`, FD stencil coefficients in `finitedifferences.f90`, Poisson solver in `poisson.f90`, SC loop convergence logic in `sc_loop.f90`

## Known Issues

- None currently. The wire g-factor was fixed using the commutator-based velocity operator (`build_velocity_matrices` in `hamiltonianConstructor.f90`). The transverse perturbation construction uses $-i [r_\alpha, H]$ element-wise on the CSR Hamiltonian, which correctly captures all k.p term contributions.
