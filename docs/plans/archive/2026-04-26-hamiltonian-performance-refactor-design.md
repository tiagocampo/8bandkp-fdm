# Hamiltonian Performance Refactor Design

Date: 2026-04-26

## Problem

The wire Hamiltonian builder `ZB8bandGeneralized` allocates and deallocates ~50+ CSR matrices per k-point during a kz-sweep (typically 100-500 points). Each `build_kp_term_*` call chains 3-5 `csr_add`/`csr_scale` operations, each allocating new CSR arrays (rowptr, colind, values). The COO triplet buffers (coo_rows, coo_cols, coo_vals) are also allocated/freed every call. Total: ~150+ heap allocations/deallocations per k-point, each involving arrays proportional to NNZ.

Additionally, `hamiltonianConstructor.f90` is 2440 lines with 34 procedures mixing bulk, QW, and wire concerns. The Bir-Pikus strain COO insertion is a 208-line manually unrolled loop. FEAST re-extracts the upper triangle on every k-point.

## Goals

1. Eliminate per-k-point allocation overhead in the wire hot path
2. Decompose hamiltonianConstructor.f90 into focused modules
3. Replace the unrolled strain COO insertion with a data-driven approach
4. Cache the FEAST upper-triangle extraction across k-points

Non-goals: optimizing the QW dense path (already fast via OpenMP+zheevx), changing the physics, or modifying the public API of `ZB8bandBulk`/`ZB8bandQW`.

## Design

### 1. Wire workspace type

Extend the existing `wire_coo_cache` into a `wire_workspace` type that pre-allocates all intermediate CSR matrices and COO buffers:

```fortran
type :: wire_workspace
  ! Pre-allocated kp-term CSR blocks
  type(csr_matrix) :: blk_Q, blk_T, blk_S, blk_SC
  type(csr_matrix) :: blk_R, blk_RC, blk_PZ, blk_PP, blk_PM, blk_A
  type(csr_matrix) :: blk_diff, blk_temp

  ! Pre-allocated COO buffers
  integer, allocatable          :: coo_rows(:), coo_cols(:)
  complex(kind=dp), allocatable :: coo_vals(:)
  integer                       :: coo_capacity = 0

  ! COO-to-CSR mapping (replaces wire_coo_cache)
  integer, allocatable          :: coo_to_csr(:)
  integer                       :: coo_nnz_in = 0
  logical                       :: initialized = .false.
end type wire_workspace
```

Public interface:
- `wire_workspace_init(ws, kpterms_2d, cfg)` -- first-call setup: builds structure
- `wire_workspace_free(ws)` -- releases all allocations

### 2. Compute-values-direct path

Each kp-term is a linear combination of kpterms_2d entries with kz-dependent coefficients. Example:

```
Q%values(k) = -(0.25 * kpterms_2d(7)%values(k) + 0.75 * kpterms_2d(8)%values(k)
             + kz2 * kpterms_2d(1)%values(k) - 2*kz2 * kpterms_2d(2)%values(k))
```

Since all kpterms share the same sparsity pattern (built by `confinementInitialization_2d`), the CSR structure (rowptr, colind) of each kp-term is identical on every k-point. Only the values change.

**First call** (`ws%initialized == .false.`): build all terms normally via `csr_add`/`csr_scale`, then stash the CSR structures in the workspace. Record which kpterms contribute to each term and with what base coefficients.

**Subsequent calls** (`ws%initialized == .true.`): for each kp-term, loop over the pre-built nonzero entries and compute the value directly from the kpterms values and kz-dependent coefficients. No csr_add, csr_scale, or csr_free calls.

The `build_kp_term_*` subroutines gain an optional `workspace` argument:
```fortran
subroutine build_kp_term_Q(kz2, kpterms_2d, blk, ws)
  real(kind=dp), intent(in) :: kz2
  type(csr_matrix), intent(in) :: kpterms_2d(:)
  type(csr_matrix), intent(inout) :: blk
  type(wire_workspace), intent(inout), optional :: ws

  if (present(ws) .and. ws%initialized) then
    ! Fast path: compute values directly into pre-allocated blk
    blk%values(:) = -(cmplx(0.25_dp, 0.0_dp, kind=dp) * kpterms_2d(7)%values &
                   + cmplx(0.75_dp, 0.0_dp, kind=dp) * kpterms_2d(8)%values &
                   + cmplx(kz2, 0.0_dp, kind=dp) * kpterms_2d(1)%values &
                   + cmplx(-2.0_dp*kz2, 0.0_dp, kind=dp) * kpterms_2d(2)%values)
  else
    ! Original path: allocate via csr_add
    ...
  end if
end subroutine
```

This requires that all contributing kpterms have identical sparsity. This is guaranteed by construction in `confinementInitialization_2d` where diagonal kpterms (1,2,3,4,5,10) and operator kpterms (6,7,8,11-16) each share the same sparsity within their group. Where a term mixes diagonal-only and operator kpterms, `csr_add` pads the diagonal terms to match the operator sparsity. The first call captures this padded structure.

### 3. Data-driven strain COO insertion

Replace the 208-line `insert_strain_coo` with a table of 32 entries defining the strain topology per grid point:

```fortran
type :: strain_entry
  integer     :: row_band, col_band  ! 0-based band offsets (0-7)
  integer     :: field_id            ! 1=EHH, 2=ELH, 3=ESO, 4=Ec, 5=S_eps, 6=R_eps, 7=P_eps_VBSO
  integer     :: cflag               ! 0=real, 1=conjg, 2=negate, 3=negate+conjg
  real(kind=dp) :: scale             ! scalar multiplier
end type

type(strain_entry), parameter :: strain_table(32) = [ &
  strain_entry(0, 0, 1, 0, 1.0_dp), &  ! HH diag: +EHH
  strain_entry(1, 1, 2, 0, 1.0_dp), &  ! LH diag: +ELH
  ... ]
```

The insertion loop becomes ~20 lines: iterate over `strain_table`, compute global row/col from band offsets, extract the bp field value, apply conjugation/negation flags.

### 4. FEAST upper-triangle cache

Add a `feast_workspace` type to `eigensolver.f90`:

```fortran
type :: feast_workspace
  integer, allocatable          :: rowptr_loc(:)   ! (N+1) upper-triangle row pointers
  integer, allocatable          :: colind_loc(:)   ! (nnz_upper) column indices
  integer                       :: nnz_upper = 0
  logical                       :: initialized = .false.
end type feast_workspace
```

In `solve_feast`: on first call, extract upper triangle and cache `rowptr_loc`/`colind_loc`. On subsequent calls, skip the counting pass and structure allocation — just copy values into `val_loc` using the cached structure. Saves one O(NNZ) scan + 2 allocations per k-point.

Also pre-allocate `E(M0)`, `X(N,M0)`, `res(M0)` in the workspace to avoid per-call allocation.

### 5. Three-module split

Decompose `hamiltonianConstructor.f90` into three modules:

| Module | File | Contents | Est. lines |
|--------|------|----------|-----------|
| `confinement_init` | `src/physics/confinement_init.f90` | `confinementInitialization_raw`, `_cfg`, `_2d`, FD helpers (`build_kpterm_block`, `applyVariableCoeffStaggered`, `build_diagonal_csr`), Dirichlet closures | ~650 |
| `hamiltonianConstructor` | `src/physics/hamiltonianConstructor.f90` | `ZB8bandQW`, `ZB8bandBulk`, `externalFieldSetup_electricField`, `add_bp_strain_dense`, profile helpers | ~500 |
| `hamiltonian_wire` | `src/physics/hamiltonian_wire.f90` | `ZB8bandGeneralized`, all `build_kp_term_*`, `insert_csr_block*`, `insert_strain_coo`, `insert_profile_diagonal`, `build_velocity_matrices_*`, `wire_workspace`, `csr_conjugate_transpose`, `negate_csr` | ~1300 |

**Import changes:**
- `confinement_init` uses: `definitions`, `parameters`, `finitedifferences`, `sparse_matrices`, `input_parser`
- `hamiltonianConstructor` uses: `definitions`, `parameters`, `finitedifferences`, `sparse_matrices`, `confinement_init`, `strain_solver`
- `hamiltonian_wire` uses: `definitions`, `parameters`, `sparse_matrices`, `input_parser`, `strain_solver`, `hamiltonianConstructor` (for build_velocity_matrices dispatch)
- App files add `use hamiltonian_wire` for wire mode

**CMakeLists.txt:** add the two new source files to the `8bandkp_common` library.

**Rationale for this split:** The wire module is self-contained — it never calls QW or bulk routines. The confinement init is shared infrastructure. The original module keeps the dense-path builders that non-wire code depends on.

## Implementation order

1. **Add `wire_workspace` + compute-values-direct** in the existing module (no split yet). Test with regression suite.
2. **Data-driven strain table** — replace `insert_strain_coo`. Test.
3. **FEAST workspace** — add `feast_workspace` to eigensolver. Test.
4. **3-module split** — refactor into separate files, update CMakeLists.txt and `use` statements. Test.

Each step is independently testable and revertible.

## Risk assessment

- **Correctness**: The compute-values-direct path must produce byte-identical results to the csr_add path. Regression tests (39 tests including wire-specific) validate this.
- **Sparsity assumption**: The fast path assumes all kpterms contributing to a given term share the same sparsity pattern. If this breaks (e.g., from mixed diagonal/operator terms), `csr_add` would produce different nnz counts. Verified by checking that the first-call path captures the merged structure.
- **Module split**: Low risk — purely structural, no logic changes. Risk is in circular dependency between modules, avoided by having `hamiltonian_wire` use `hamiltonianConstructor` (not vice versa).

## Expected impact

For a 100x100 wire grid with 200 k-points:
- **Before**: ~10,000 heap allocations (50/point × 200), each ~1 MB
- **After**: ~200 allocations (COO buffer + FEAST workspace per point), ~98% reduction
- **FEAST caching**: saves 2 × O(NNZ) scans per k-point
- **Wall-clock**: estimated 2-5× speedup for wire kz-sweeps (allocation overhead dominates for N < 50000)
