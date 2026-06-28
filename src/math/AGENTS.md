# Math Layer

Owns numerical primitives consumed by `src/physics/`: finite-difference operators,
linear algebra interfaces, eigensolvers, sparse matrix formats, and pure
mathematical algorithms (Pfaffian, etc.). Does NOT own: material parameters
(`src/core/parameters.f90`), Hamiltonian construction (`src/physics/`),
or input/output (`src/io/`).

## Module Inventory

| File | Module | Lines | Role |
|------|--------|------:|------|
| `finitedifferences.f90` | `finitedifferences` | — | FD stencils/orders 2-10, Toeplitz matrices, Vandermonde derivative and interpolation solvers |
| `linalg.f90` | `linalg` | — | Centralized LAPACK/BLAS/MKL/PARDISO interface declarations; FEAST under `#ifdef USE_MKL_FEAST` |
| `eigensolver.f90` | `eigensolver` | — | Polymorphic eigensolver dispatch (LAPACK / FEAST) — `make_eigensolver(config)` factory |
| `sparse_matrices.f90` | `sparse_matrices` | — | CSR/COO sparse matrix types + SpMV/SpBLAS helpers |
| `geometry.f90` | `geometry` | — | Spatial grid types (`spatial_grid`) + `npoints()` accessor; band-major index mapping |
| `pfaffian.f90` | `pfaffian` | 525 | Real/complex skew-symmetric Pfaffian (Laplace for n<=12, Parlett-Reid tridiagonalization for n>12) + `kitaev_majorana_number` wrapper (Issue 01) |

## Dependency DAG

```
Layer 0 (leaves):  definitions (upstream — owns kinds/constants),
                   finitedifferences, linalg, sparse_matrices, geometry
Layer 1:           eigensolver (uses linalg)
Layer 2:           pfaffian (uses definitions only)
```

No cycles. All consumers depend on these via `use` only — no direct calls
into the layer from outside `src/math/`.

## Contracts & Invariants

- **Skew-symmetric Pfaffian convention**: `Pf([[0,a],[-a,0]]) = +a` (Kitaev).
- **Input skew-symmetrization guard**: callers pass real/complex matrices;
  `real_pfaffian` / `complex_pfaffian` symmetrize as `(A − Aᵀ)/2` before
  reduction. This guards against tiny roundoff that would otherwise cause
  the recurrence to produce a wrong sign.
- **n parity**: Pfaffian is defined only for even-dimensional input; odd n
  returns 0. The dispatcher routes n ≤ 12 to Laplace, n > 12 to Parlett-Reid.
- **`kitaev_majorana_number`** signature:
  `kitaev_majorana_number(H_k_array, k_par_values, omega_struct)`.
  `H_k_array` is a 3D array of BdG matrices (one per PHS-invariant
  momentum); the wrapper evaluates Pf(skew(H·ω)) at each and returns the
  sign of the product (-1 topological, +1 trivial, 0 gap closure).
  Default `ω = τ_y ⊗ I_N` (Kronecker product, only diagonal of off-diagonal
  block is nonzero).

## Patterns

### Adding a new math primitive
1. Place in `src/math/` if it's pure math (no physics-specific assumptions).
2. Add row to this `AGENTS.md` Module Inventory.
3. Wire into `src/CMakeLists.txt` COMMON_SOURCES.
4. Write unit tests in `tests/unit/` with pFUnit (single-line asserts).
5. Add CMake entry in `tests/CMakeLists.txt`.

## Anti-patterns

- Never duplicate math algorithms in physics — extract to this layer first.
- Never add physics-specific assumptions to a math module (e.g., no material
  parameters, no band structure references).
- Never modify LAPACK/MKL interfaces outside `linalg.f90` (centralized).

## Pitfalls

- Stale `.mod` files shadowing fresh ones in `build/` — see root CLAUDE.md.
- `pfaffian.f90` is ~525 lines, exceeding the 300-line guideline; splitting
  was deferred as YAGNI (vendored pfapack math + thin wrapper).
- Parlett-Reid n > 12 was self-reported broken in initial implementation;
  Fix Round 1 (Issue 01 review) corrected the Householder norm (now uses
  full ||x|| = √(A(k+1,k)² + s)) and the similarity update (now `+τ·(v·wᵀ −
  w·vᵀ)` — the vᵀAv term vanishes for skew A).

## Related Context

- Physics engine: `src/physics/AGENTS.md`
- Material parameters: `src/core/parameters.f90`
- Test infrastructure: `tests/integration/AGENTS.md`