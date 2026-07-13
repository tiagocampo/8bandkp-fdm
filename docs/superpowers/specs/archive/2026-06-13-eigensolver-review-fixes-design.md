# Design: Eigensolver Standardization — Review Fixes

**Date:** 2026-06-13
**Branch:** `refactor/eigensolver-standardization`
**Status:** Approved (brainstorming complete)
**Predecessor PRD:** `.scratch/eigensolver-standardization/PRD.md`

## Context

An adversarial code review of the `refactor/eigensolver-standardization` branch (against its PRD) surfaced 12 findings in three buckets: performance regressions, correctness hardening, and cleanup / PRD-compliance gaps. This design specifies how to resolve all of them in a single pass. No physics changes — purely solver, parser, validation, and documentation plumbing.

### Findings addressed

| # | Bucket | Severity | File | Summary |
|---|--------|----------|------|---------|
| 1 | Perf | High | `hamiltonian_qw.f90:216` | QW CSR rebuilds from scratch every k-point; fast-path unimplemented (TODO at :313) |
| 2 | Perf | High | `eigensolver.f90:794` | Dense solver allocates all workspace per-call in the k-sweep hot loop |
| 3 | Hardening | Medium | `eigensolver.f90:803` | Missing `info` check after LAPACK workspace query |
| 4 | Cleanup | Low | `eigensolver.f90:745` | Dense solver's `solve_sparse` does needless CSR↔dense round-trip |
| 5 | Cleanup | Low | `eigensolver.f90:436` | Legacy `solve_dense_lapack` uses stale emin/emax heuristic |
| 6 | PRD compliance | Medium | `input_parser.f90:839` | `[feast]` section silently ignored — violates PRD US-20 |
| 7 | Cleanup | Low | `sc_loop.f90:21` | Dead import of removed `solve_sparse_evp` |
| 8 | Hardening | Low | `defs.f90:791` | FEAST+INDEX rejected at wrong validation layer |
| 9 | Hardening | Low | `main.f90:747` | Silent eigenvalue truncation in QW FEAST path |
| 10 | Cleanup | Low | `lecture_12_extending.py:86` | ARPACK references in lecture script |
| 11 | — | — | — | (Refuted: config converter correctly writes `mode`) |
| 12 | Cleanup | Low | `eigensolver.f90:400` | Legacy public procedure with conflicting semantics (same root as #5/#7) |

## Decisions (locked during brainstorming)

1. **Scope:** Full — fix all three buckets in one pass.
2. **`[feast]` backward-compat (#6):** Hard `error stop` with a migration message. No silent ignore, no mapping. Rationale: silent ignore produces wrong results with no signal; a warning can be missed in batch output; fail-fast matches the codebase "no silent corrections" convention and the project's pre-release status.
3. **Dead public procedures (#5/#7/#12):** Remove `solve_sparse_evp` and `solve_dense_lapack` entirely. Rationale: a "direct interface" with zero callers and divergent semantics is a liability; the PRD's "remains available" was a transition hedge now obsolete.
4. **Performance approach (#1/#2):** Full PRD-faithful fix — solver-object workspace caching for the dense solver, and the PRD-issue-04 10-assembly-function fast-path for QW CSR.

## Architecture

Two caching models, one per perf finding; cleanup is deletion + small guards.

### Data flow

Unchanged. Same inputs, same `eigensolver_result` outputs, same call sites. Only *where workspace lives* (dense) and *whether CSR is rebuilt* (QW) change. All golden regression eigenvalues remain bit-identical.

### Files touched

- `src/math/eigensolver.f90` — #2 (workspace caching), #3 (info check), #5/#7/#12 (deletions)
- `src/physics/hamiltonian_qw.f90` — #1 (QW CSR fast-path)
- `src/io/input_parser.f90` — #6 (`[feast]` hard error)
- `src/core/defs.f90` — #8 (FEAST+INDEX validation, new check I15)
- `src/apps/main.f90` — #9 (truncation warning)
- `src/physics/sc_loop.f90` — #7 (stale import removal)
- `scripts/lecture_12_extending.py` — #10 (ARPACK lines)
- `CLAUDE.md`, `AGENTS.md` — drop "solve_sparse_evp remains available"

No changes to: k.p block table, basis ordering, Bir-Pikus signs, Zeeman table, FD stencils, Poisson solver, SC convergence logic.

## Component Designs

### 1. Dense solver workspace caching (#2, #3)

**Type extension** in `eigensolver.f90`:

```fortran
type, extends(eigensolver_base) :: dense_lapack_solver_t
  private
  ! Cached workspace — sized on first call, reused thereafter.
  integer                              :: cached_n = 0
  complex(kind=dp), allocatable        :: A_buf(:,:), Z_buf(:,:), work(:)
  real(kind=dp), allocatable           :: W_buf(:), rwork(:)
  integer, allocatable                 :: iwork(:), ifail(:)
contains
  procedure :: solve_dense  => dense_solve_dense_dispatch
  procedure :: solve_sparse => dense_solve_sparse_dispatch
  final     :: dense_lapack_solver_final
end type
```

**`dense_solve_dense_dispatch` refactor:**
- On entry, if `N > cached_n`: deallocate-if-allocated, then allocate `A_buf(N,N)`, `Z_buf(N,N)`, `W_buf(N)`, `rwork(7N)`, `iwork(5N)`, `ifail(N)`; set `cached_n = N`. (These depend only on N, not mode.)
- `A_buf = H` each call (unavoidable — LAPACK destroys its input).
- Workspace query each call (cheap; one call on length-1 array). **#3:** check `info` immediately after the query; `error stop` with a diagnostic if nonzero *before* reading `work(1)`. If returned `lwork > size(work)`, reallocate `work`; otherwise reuse. `work` is grow-only — correct for mode switches since a larger `work` is always valid for a smaller request.
- Real solve writes into `A_buf`/`Z_buf`/`W_buf`; results copied into `result` as today. No per-call allocate/deallocate of the six N-dependent buffers.

**Finalizer** `dense_lapack_solver_final`: deallocates all seven arrays, each guarded by `if (allocated(...))`. Idempotent.

**Thread-safety:** the dense QW path allocates `solver_loc` per thread inside `!$omp parallel private(...)` (`main.f90:778`), so each thread's solver holds its own cached workspace — no shared mutable state. Bulk, Landau, optics, SC, and topology paths are serial or also per-thread; all dense callers benefit with zero interface changes.

**Cost reduction:** kills the ~384 MB–500 MB transient allocation per k-point per thread, plus one redundant workspace query per call. Restores pre-refactor hoisted-allocation performance, now inside the polymorphic layer.

### 2. QW CSR fast-path (#1)

**Goal:** on the 2nd+ k-point, recompute only values — no dense matrices, no `dense_to_csr_block` O(N²) scan, no COO→CSR sort.

**Key realization:** each kp-term block is tridiagonal + diagonal, so its CSR structure (`rowptr`/`colind`) is fixed across k-points — only `values` change. The existing infrastructure already supports this: `qw_workspace` clones block structures (lines 319–330), owns the COO buffers via `move_alloc` (lines 333–335), and has absorbed COO→CSR sort-cache fields (`coo_to_csr`, `coo_cache_valid`, lines 57–59) — but the fast path is unwired. The wire path already uses cached finalize (`hamiltonian_wire.f90:438`), proving the pattern.

**Per-term value-update routine:**

```fortran
! Walks blk%rowptr/colind, recomputes each non-zero's value from
! kpterms(i,j,·) and the current kx,ky. Structure is assumed fixed.
subroutine update_kp_term_values(blk, kpterms, kx, ky, term_id)
```

**Slow path (first call):** unchanged — build dense once (lines 216–237), convert to CSR, COO-insert, `finalize_coo_to_csr` **with the workspace sort cache populated**, clone block structures + COO buffers into `qw_workspace`. One-time O(N²) cost.

**Fast path (subsequent calls):**
1. `update_kp_term_values` for each of the 10 cached blocks (Q,T,S,SC,R,RC,PZ,PP,PM,A) — O(NNZ) each, no dense allocation. Recompute `blk_diff`/`blk_temp` via their cached `csr_add` structure.
2. Re-scatter into the pre-allocated COO buffers (`insert_main_blocks` + `insert_profile_diagonal`/`insert_strain_coo`/`insert_zeeman_coo`), resetting `coo_idx=0` — overwrites row/col/val in place, no allocation.
3. `finalize_coo_to_csr` **with the cached `coo_to_csr` mapping** — skips the sort and rowptr/colind build, just re-sums values into the pre-positioned `HT_csr%values`. O(NNZ).

**Correctness gate (risk control):** the formulas in `update_kp_term_values` are a direct refactor of the dense expressions at lines 223–237 — same arithmetic, different iteration order (CSR-walk vs dense double-loop). `verify_qw_sparse_solver.py` is extended to assert **fast-path eigenvalues == slow-path eigenvalues** to ≤1e-10 across multiple k-points, catching any formula drift. This two-paths-same-result risk is the single biggest concern; the test is the guard.

**Cost reduction:** slow path O(N²) per call → fast path O(NNZ) ≈ O(N). The 6.4 GB transient allocation churn (10 dense matrices × 1000 k-points) drops to near zero. Fulfills PRD US-25.

### 3. Cleanup & hardening (#3, #4, #5, #6, #7, #8, #9, #10, #12)

- **#3 info check** — covered in §1.
- **#4 CSR↔dense round-trip** — `dense_lapack_solver_t%solve_sparse` (CSR→dense→solve_dense) remains for interface completeness but is documented as a non-hot convenience path. No current caller triggers it on large matrices; no change needed beyond the doc note.
- **#5/#7/#12 dead code** — delete `solve_sparse_evp` + `solve_dense_lapack` from `eigensolver.f90` (bodies, `public ::` declarations, and the dispatch block at lines 173–196). Remove `solve_sparse_evp` from the `use` in `sc_loop.f90:21`. Update CLAUDE.md + AGENTS.md.
- **#6 `[feast]` hard error** — in `input_parser.f90` `read_config`, after parsing, check for a `[feast]` table; if present, `error stop '[feast] section removed — rename to [solver] (method/mode/emin/emax/m0). See docs/reference/input-reference.md'`. No mapping.
- **#8 FEAST+INDEX at validation layer** — move the cross-constraint from `eigensolver_config_validate` (eigensolver.f90:649) into `simulation_config_validate` (`defs.f90`), as new check **I15** immediately after I13/I14 (lines 791–807). Keep a defensive `error stop` in the eigensolver too (defense in depth), but the user-facing error fires from input validation with clear config context.
- **#9 truncation warning** — in `main.f90:748`, when `result_bs%nev_found > iuu-il+1`, print a warning naming the count discarded before clamping. Audit sibling QW FEAST paths in optics/gfactor for the same pattern.
- **#10 ARPACK docs** — delete the two ARPACK mentions in `scripts/lecture_12_extending.py` (lines 86, 422).

## Testing

All fixes verified by equivalence to existing behavior, not new physics values.

| Fix | Test | Type |
|---|---|---|
| #2 dense caching | Existing QW k-sweep regression tests (golden eigenvalues unchanged) + new unit test: cached solver eigenvalues identical to fresh-solver across repeated calls | regression + unit |
| #1 QW fast-path | Extend `verify_qw_sparse_solver.py`: fast-path (2nd+ call) eigenvalues match slow-path to ≤1e-10 across multiple k-points | integration |
| #3 info check | Validation rejection test (shell exit-code) for degenerate-matrix path | integration |
| #6 `[feast]` hard error | Extend `test_validate_rejects_bad_configs.sh`: TOML with `[feast]` exits nonzero | integration |
| #8 FEAST+INDEX in defs | Extend `test_validate_rejects_bad_configs.sh`: error originates from `defs` validation | integration |
| #9 truncation | New check in `verify_qw_sparse_solver.py`: wide window triggers warning (stdout grep) | integration |
| #5/#7/#12 dead code | Build succeeds with deletions; grep confirms zero remaining references | build + grep |
| #10 ARPACK | `lecture_12_extending.py` runs clean | script |

**Full-suite gate:** `ctest --test-dir build` must remain green (113+ tests). No golden-output reference data is regenerated — all comparisons are against existing references.

## Out of Scope

- LDOS/PARDISO refactoring (unchanged from original PRD).
- Workspace-query extraction into a shared helper beyond the dense solver (subsumed by solver-object caching).
- GPU offload / batched LAPACK / ELPA (future variants via `make_eigensolver`).
- Any change to material parameters, basis ordering, Bir-Pikus signs, Zeeman table, FD stencil coefficients, Poisson solver, or SC loop convergence logic.

## Risks

1. **QW fast-path formula drift** (highest) — two code paths compute the same kp-term values. Mitigated by the ≤1e-10 equivalence test and by deriving `update_kp_term_values` directly from lines 223–237.
2. **COO sort-cache wiring** — the QW path's absorbed cache fields must be correctly populated on the slow path and consumed on the fast path. Mitigated by reusing the wire path's proven `finalize_coo_to_csr(..., coo_cache=...)` mechanism.
3. **`work` grow-only correctness across mode switches** — a larger workspace is always valid for LAPACK routines requesting less. Verified by the multi-mode dense unit test.
