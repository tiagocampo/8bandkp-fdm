# ADR 0001: Fat Derived Type for Simulation Setup

## Status

Accepted

## Context

The 8-band k.p solver has three application programs (`bandStructure`, `gfactorCalculation`, `opticalProperties`) that each independently implement the same physics initialization pipeline: confinement initialization, strain computation, eigensolver configuration, and self-consistent Schrodinger-Poisson integration. This duplication makes certain physics combinations architecturally impossible (e.g., optics + SC + strain) and causes bug-fix divergence across apps.

We need a `simulation_setup` module that encapsulates this pipeline behind a single entry point. The key architectural decision is how to represent the setup state across three geometries (bulk, QW, wire) that have fundamentally different data representations (dense vs. sparse).

## Decision

Use a single **fat derived type** with allocatable components for both dense (bulk/QW) and sparse (wire) paths. Only one path's components are allocated based on `cfg%confinement`. Dispatch is via `select case` on the confinement integer.

```fortran
type :: simulation_setup
  integer :: confinement = -1
  ! Dense path (bulk/QW)
  real(dp), allocatable :: profile(:,:), kpterms(:,:,:)
  complex(dp), allocatable :: HT(:,:)
  ! Sparse path (wire)
  real(dp), allocatable :: profile_2d(:,:)
  type(csr_matrix), allocatable :: kpterms_2d(:)
  type(csr_matrix), allocatable :: HT_csr_ptr
  class(eigensolver_base), allocatable :: eigen_solver
  ! Velocity matrices (opt-in)
  type(csr_matrix) :: vel(3)
end type
```

Wire components are individually allocatable (not just their inner arrays) to prevent gfortran auto-finalization of unallocated `csr_matrix`/`wire_workspace` types when the setup is used for bulk/QW paths.

## Alternatives Considered

### Polymorphic subtypes with class hierarchy

```fortran
type, abstract :: setup_base
contains
  procedure(build_H_iface), deferred :: build_H
end type

type, extends(setup_base) :: bulk_setup
  ! 8x8 dense only
end type

type, extends(setup_base) :: qw_setup
  real(dp), allocatable :: profile(:,:), kpterms(:,:,:)
end type

type, extends(setup_base) :: wire_setup
  type(csr_matrix), allocatable :: kpterms_2d(:)
end type
```

**Pros**: Type-safe — each variant carries only its own fields. No accidental access to wire fields from bulk code.

**Cons**: Requires `select type` at every call site that accesses geometry-specific fields (profile, kpterms, CSR). Fortran's polymorphism adds boilerplate without real safety gain in a codebase with exactly 3 variants. Deferred procedures require every new operation (build_H, solve, build_velocity_matrices) to be added to the base type, creating a wide interface that changes with every new feature.

### Strategy pattern with composition

Each setup variant holds a strategy object for Hamiltonian construction, eigensolving, etc. The setup type delegates to strategies.

**Pros**: Open for extension — new strategies without changing the setup type.

**Cons**: Fortran lacks first-class procedures and closures. Strategies must be passed as procedure pointers or polymorphic objects, adding complexity. Over-engineering for a codebase with exactly 3 geometries and no extensibility requirement beyond the C4 Hamiltonian block table refactoring (now completed in `hamiltonian_blocks.f90`).

## Consequences

### Positive

- **Simplicity**: One type, one init routine, one free routine. No class hierarchy, no deferred bindings, no `select type`.
- **Single entry point**: `simulation_setup_init(cfg, setup)` replaces 3 separate initialization paths in each app.
- **New physics combinations**: Optics + SC + strain becomes possible because `simulation_setup_init` handles all pipelines internally.
- **Bug locality**: A bug in strain setup is fixed once in the setup module, not three times across apps.

### Negative

- **No compile-time type safety**: A bulk path can accidentally access `setup%kpterms_2d` (returns unallocated, caught at runtime).
- **Wide type interface**: The type exposes fields for all three geometries. Callers must know which fields are valid for their confinement.
- **Allocatable wrapper pattern**: Wire components like `HT_csr_ptr` are allocatable types wrapping allocatable components. Double indirection is verbose but necessary for lifetime control.

### Dispatch seam for C4 Hamiltonian block table (COMPLETED)

The `select case` dispatch in `setup_build_H` provided a clean seam for the C4 refactoring to a unified Hamiltonian block structure. This has been implemented: `hamiltonian_blocks.f90` defines a 52-entry k.p block table (`get_kp_block_table()`) consumed by both dense (`hamiltonianConstructor.f90`) and COO (`hamiltonian_wire.f90`) builders. Strain and Zeeman block tables (`get_strain_table()`, `get_zeeman_table()`) are defined in `strain_solver.f90`. The dispatch API (`setup_build_H(setup, cfg, kvec, HT_out)`) remains unchanged. No polymorphic types were needed.

### `read_config` split rationale

Splitting `read_and_setup` into `read_config` (pure parser) + `simulation_setup_init` (physics setup) separates parsing from initialization. This enables:
- Unit testing of setup without writing `input.cfg` to disk
- Multiple init calls with different parameters from the same parsed config
- Clear ownership: the parser owns file I/O, the setup module owns physics state

The `read_and_setup` wrapper is removed after all callers are migrated.
