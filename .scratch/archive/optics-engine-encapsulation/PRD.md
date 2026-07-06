**Status**: COMPLETE (2026-07-05)

Status: ready-for-agent

# PRD: Encapsulate Optical Spectra State into `optics_engine` Type

## Problem Statement

The `optical_spectra` module stores 15 allocatable arrays and 3 scalar flags as module-level `save` variables. The init/accumulate/finalize/cleanup lifecycle is a calling convention that callers must memorize — it is part of the interface but not enforced by the compiler. The finalize routine conflates two responsibilities (applying physical prefactors and writing output files), preventing callers from intervening between them (e.g., applying excitonic corrections to the finalized arrays before writing). No two independent optical calculations can run concurrently, and tests cannot create isolated instances without risking global state pollution.

## Solution

Bundle all module-level state into an `optics_engine` derived type with public components, module procedures that take the type as their first argument, and a finalizer. Split `optics_finalize` into `optics_apply_prefactor` (modifies internal arrays) and `optics_write_output` (writes files), so callers control the pipeline: apply prefactor, optionally intervene, then write output. The type stores its `optics_config` at initialization, eliminating the need to pass `optcfg` to every routine. The existing `save` variables are removed, and the four public module variables (`E_grid`, `alpha_te`, `alpha_tm`, `nE`) become type components accessed via `oe%E_grid`, etc.

## User Stories

1. As a physics researcher, I want to compute optical absorption of a quantum well and then apply excitonic corrections before writing output, so that the final files contain the corrected spectra.
2. As a physics researcher, I want the optical spectrum pipeline to be clearly separated into distinct steps (accumulate, apply prefactor, write output), so that I understand what each step does.
3. As a developer, I want the optics state to live in a derived type rather than module-level `save` variables, so that multiple independent optical calculations can coexist.
4. As a developer, I want the optics lifecycle (init/accumulate/prefactor/write/free) to be enforced by the type's structure, so that misuse is caught at the interface rather than by convention.
5. As a developer, I want a finalizer on the optics type, so that resources are released automatically on scope exit.
6. As a developer, I want the optics type to store its configuration at initialization, so that I don't have to pass `optcfg` to every routine.
7. As a developer, I want the accumulate routines to take only physics data (eigenvalues, eigenvectors, velocity matrices, k-weight, band counts, Fermi level), so that the signatures reflect what the routine actually needs.
8. As a developer, I want to write unit tests for optics that create and destroy instances without polluting global state, so that tests are isolated and repeatable.
9. As a codebase maintainer, I want `E_grid`, `alpha_te`, `alpha_tm`, and `nE` to be type components rather than public module variables, so that the module's interface is the type, not scattered exports.
10. As a codebase maintainer, I want `optics_apply_prefactor` and `optics_write_output` to be separate routines, so that I can intervene between them (e.g., exciton corrections, custom post-processing).
11. As a codebase maintainer, I want the caller to control the pipeline (prefactor → optional intervention → write → free), so that the optics module does not dictate when output is written.
12. As a codebase maintainer, I want the `optics_engine` type to follow the same patterns as `wire_workspace` (module procedures, public components for hot-path access), so that the codebase is consistent.
13. As a codebase maintainer, I want `exciton.f90` to require no changes, so that the encapsulation is local to the optics module and its direct caller.
14. As a developer, I want `gain_reset` to take the type as an argument, so that the gain quasi-Fermi level state is properly scoped to the instance.
15. As a developer, I want `compute_intersubband_transitions` to remain unchanged (no module state dependency), so that it continues to work as a standalone utility.

## Implementation Decisions

### Type design: `optics_engine`

A derived type with public allocatable components, following the `csr_matrix` precedent (public components for hot-path access) and the `wire_workspace` precedent (module procedures with type as first argument). The type stores a copy of `optics_config` at initialization so that subsequent routines do not need it as a separate argument.

Key components:
- `optcfg`: stored configuration (linewidths, refractive index, temperature, feature flags)
- `nE`: number of energy grid points (default 0, used as initialization guard)
- `E_grid(nE)`: photon energy grid in eV
- `alpha_te(nE)`, `alpha_tm(nE)`, `alpha_isbt(nE)`: always-allocated accumulation arrays
- `alpha_gain_te(nE)`, `alpha_gain_tm(nE)`: conditionally allocated when gain enabled
- `spont_te(nE)`, `spont_tm(nE)`: conditionally allocated when spontaneous emission enabled
- `alpha_te_up(nE)`, `alpha_te_dw(nE)`, `alpha_tm_up(nE)`, `alpha_tm_dw(nE)`: conditionally allocated when spin-resolved enabled
- `mu_e`, `mu_h`: gain quasi-Fermi levels (default 0)
- `gain_fermi_computed`: guard flag (default `.false.`)

All allocatable components have default unallocated state. The type has a finalizer that delegates to `optics_free`.

### Calling convention: module procedures

All public routines take `type(optics_engine), intent(inout) :: oe` as their first argument. This follows the `wire_workspace` pattern and avoids Fortran type-bound procedure complexity. The caller replaces `cfg%optics` with `oe` in all call sites — a minimal diff.

### Config storage: stored in type

`optics_init(oe, optcfg)` stores a copy of `optcfg` in `oe%optcfg`. All subsequent routines (`accumulate`, `apply_prefactor`, `write_output`) read config from `oe%optcfg` rather than requiring it as a separate argument. The config is small (scalar physical parameters and feature flags) and read-only during the lifecycle.

### Finalize split: `apply_prefactor` + `write_output`

The current `optics_finalize` is split into two routines:

1. `optics_apply_prefactor(oe)` — applies the physical prefactor `(2*pi*e^2) / (n_r * c * eps0_AA * hbar^2 * E)` to all accumulation arrays. Converts from raw |dH/dk|^2 units to cm^-1. No file I/O.
2. `optics_write_output(oe)` — writes all spectra files (absorption TE/TM/ISBT, gain TE/TM, spontaneous TE/TM, spin-resolved). Reads from `oe%alpha_te`, etc.

This enables the caller pipeline:
```
optics_init(oe, cfg%optics)
  → accumulate loop
  → optics_apply_prefactor(oe)
  → optional: exciton corrections on oe%alpha_te, oe%alpha_tm
  → optics_write_output(oe)
  → optics_free(oe)
```

Files are written exactly once, with the correct (possibly corrected) values.

### Module location: same file

The `optics_engine` type and all routines remain in `optical_spectra.f90`. The file is already 1089 lines and cohesive. Splitting a Fortran module across files requires submodules (F2008), which the codebase does not use. The type definition adds ~20 lines; routine signature changes are mechanical.

### Public interface changes

Removed:
- `public :: E_grid, alpha_te, alpha_tm, nE` (module variables)
- `public :: optics_finalize` (replaced by two routines)
- `public :: optics_cleanup` (replaced by `optics_free`)

Added:
- `public :: optics_engine` (the type)
- `public :: optics_init` (takes `oe` + `optcfg`)
- `public :: optics_free` (replaces `optics_cleanup`)
- `public :: optics_apply_prefactor` (replaces first half of `optics_finalize`)
- `public :: optics_write_output` (replaces second half of `optics_finalize`)

Modified (gain `oe` first arg, lose `optcfg` arg):
- `optics_accumulate`
- `optics_accumulate_spontaneous`
- `compute_gain_qw`
- `compute_isbt_absorption`
- `gain_reset`

Unchanged (no module state dependency):
- `compute_intersubband_transitions`
- Private helpers: `lineshape_voigt`, `find_quasi_fermi`, `find_quasi_fermi_holes`, `z_dipole`

### Caller impact: `main_optics.f90`

The `main_optics.f90` changes are mechanical:
- Add `type(optics_engine) :: oe` declaration
- Replace `call optics_init(cfg%optics)` → `call optics_init(oe, cfg%optics)`
- Replace `call optics_accumulate(cfg%optics, ...)` → `call optics_accumulate(oe, ...)`
- Replace `call optics_finalize(cfg%optics)` → `call optics_apply_prefactor(oe)`
- Add optional exciton intervention (using `oe%alpha_te` instead of `alpha_te`)
- Add `call optics_write_output(oe)` after exciton (or directly if no exciton)
- Replace `call optics_cleanup()` → `call optics_free(oe)`

The exciton block (QW case) changes `E_grid` → `oe%E_grid`, `alpha_te` → `oe%alpha_te`, `alpha_tm` → `oe%alpha_tm`, `nE` → `oe%nE`. The `exciton.f90` module requires no changes — it takes arrays as arguments.

### No changes to `simulation_setup`

The `simulation_setup` module does not call any optics routines and does not need to own the `optics_engine`. The optics module and simulation_setup are orthogonal: `simulation_setup` provides physics infrastructure (Hamiltonian, velocity matrices), while `optics_engine` manages spectrum accumulation. The `main_optics.f90` app orchestrates both.

### Consistency with ADR 0001

ADR 0001 established the fat-type pattern for `simulation_setup`. This PRD applies the same pattern to a different domain: a flat derived type with allocatable components, `select case` dispatch is unnecessary here (optics routines are geometry-agnostic), and a finalizer for cleanup. The design is consistent with the codebase's direction.

## Testing Decisions

### What makes a good test

Tests assert on observable outcomes through the public interface: after init, arrays are allocated and zeroed; after accumulate with known inputs, arrays have expected values; after apply_prefactor, values are scaled correctly. Tests do not assert on internal state beyond what the public components expose.

### Modules to test

The `optics_engine` type in `optical_spectra.f90`. A new pFUnit test file `tests/unit/test_optical_spectra.pf`.

### Prior art

The existing `tests/unit/test_simulation_setup.pf` (12 tests) follows the same pattern: construct a config, call init, assert on type components, call free. The optics tests will follow this template:
1. Construct a minimal `optics_config` with known parameters (energy grid range, linewidths, etc.)
2. Call `optics_init(oe, optcfg)` — assert arrays allocated, values zeroed
3. Call accumulate routines with fixture data — assert accumulation
4. Call `optics_apply_prefactor(oe)` — assert prefactor applied
5. Call `optics_free(oe)` — assert arrays deallocated

### Test fixtures

Minimal fixtures needed:
- A small `optics_config` (3-5 energy points, narrow linewidth)
- Small eigenvalue/eigenvector arrays (e.g., 2 CB + 2 VB states)
- CSR velocity matrices with known values (use existing CSR helpers from `tests/support/`)

### Regression safety

The existing regression tests for optics (golden-output file comparison) verify end-to-end correctness. The new unit tests verify the type lifecycle. Both layers are needed: unit tests for the encapsulation contract, regression tests for physics correctness.

## Out of Scope

- **Pauli matrices in `gfactor_functions.f90`**: 4 `save` variables (3 constant 8x8 matrices + 1 guard flag). Lazy-init-once pattern, no accumulate/finalize lifecycle, no cleanup needed. Addressable as a trivial follow-up but not bundled here.
- **`hamiltonian_wire.f90` module-level state**: Already encapsulated in `wire_workspace` and `wire_coo_cache` derived types with finalizers. No work needed.
- **Candidate 4 (unified Hamiltonian block structure)**: Orthogonal refactoring. The `optics_engine` type does not touch Hamiltonian construction.
- **Candidate 5 (split simulation_config)**: Orthogonal refactoring. The `optics_engine` type stores a copy of `optics_config` but does not change its definition.
- **Changing `exciton.f90`**: It already takes arrays as arguments — no changes needed.
- **Moving optics into `simulation_setup`**: The two modules are orthogonal. No ownership change.

## Further Notes

### Relationship to the architecture review

This is Candidate 3 from the 2026-05-25 architecture review. Candidates 1 (QW/wire duplication collapse) and 2 (deepen app orchestration via `simulation_setup`) are complete. This PRD implements the next highest-leverage candidate.

### Lifecycle invariant

The `optics_engine` type enforces a lifecycle invariant through its interface design:
- `nE == 0` means uninitialized (all routines guard on this)
- After `optics_init`: `nE > 0`, arrays allocated and zeroed
- After accumulate loop: arrays contain raw accumulated values
- After `optics_apply_prefactor`: arrays contain physically scaled values (cm^-1)
- After `optics_free`: all arrays deallocated, `nE == 0`

The previous design enforced this lifecycle as a calling convention (init → accumulate → finalize → cleanup). The type makes it a structural property.

### Blast radius

Three files: `optical_spectra.f90` (type + routine signatures), `main_optics.f90` (call sites), and a new test file. No other module imports `optical_spectra`'s public variables (`E_grid`, `alpha_te`, etc.) — only `main_optics.f90` does. The `exciton.f90` module takes arrays as arguments and is unaffected.
