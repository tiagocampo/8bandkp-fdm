**Status**: COMPLETE (2026-07-13) — ticket 02 of `.scratch/archive/bdg-evaluator-pfaffian/`.

Type: research
Status: resolved
Claimed-by: Tiago (via wayfinder session 2026-07-13)
Blocked by: 01
Resolved: 2026-07-13

# 02 — Slim witness location

## Question

The slim projected Pfaffian witness (`wire_pfaffian_witness`) lives at `topological_analysis.f90`. Per `src/physics/AGENTS.md:46-58`, both `bdg_observables` and `topological_analysis` are L3 hubs; the dependency graph is acyclic. Confirm:

- Which module imports from which for the slim witness to be callable from `eval_bdg_point`.
- Whether `bdg_observables.f90` (L3 hub) calls `topological_analysis.f90` (also L3 hub) for the slim witness, or whether the witness is relocated into `bdg_observables.f90` (or a new sibling module) so the call direction is clean.
- Whether the import direction risks a cycle once the Kitaev wrapper (U3, also in `topological_analysis.f90`) is also called from `bdg_observables.f90`.
- Whether the slim witness's input shape (eigenvectors + a small descriptor for S1 vs S2 projection) is compatible with the API shape decided in ticket 01.

Also: read `test_wire_pfaffian_witness.pf` and confirm what its strict `@assertTrue(s1 == s2 .and. s1 /= 0)` exercises, whether the current `@todo U13` flag is a comment or a hard `@todo`, and whether the test will pass under the slim witness when it becomes the gate row.

Deliverable: a one-page report at the ticket resolution with file:line citations, the import-direction decision, and the noise-floor tolerance assessment (sub-decision of the fog item noted on the map's "Not yet specified").

## Constraints

- No code changes — research only. The actual move/relocation is a separate execution step downstream of the decisions in tickets 01 + 02 + 03.
- No new module unless the import-direction analysis proves no other option.
- Cite file:line for every claim.

## Out of scope for this ticket

- The retirement of `compute_z2_gap` (ticket 03).
- The call-site migration in `main_topology.f90` (ticket 04).
- The acceptance-gate wiring (ticket 05).

## Answer

Research report — dep-graph analysis + implementation placement + test `@todo` status.

### 1. Premise correction: `bdg_observables` is NOT an L3 hub

The ticket's premise ("Per `src/physics/AGENTS.md:46-58`, both `bdg_observables` and `topological_analysis` are L3 hubs") is **incorrect** — `bdg_observables` is missing from the dep DAG entirely. From `bdg_observables.f90:20`:

```fortran
use definitions, only: dp
```

That's the entire `use` block. `bdg_observables` is effectively an **L0 leaf** (definitions-only) that consumes nothing but `dp`. The Inventory table (`src/physics/AGENTS.md:44`) lists the module but the DAG (`src/physics/AGENTS.md:49-58`) omits it. This is a doc-drift gap worth flagging to ticket 07's plan close-out sweep, but it doesn't change the dep analysis — `bdg_observables` is currently a lean leaf, not a hub.

### 2. `topological_analysis` does NOT import `bdg_observables`

Verified directly at `src/physics/topological_analysis.f90:3-11`:

```fortran
use definitions, only: dp, pi_dp, simulation_config, spatial_grid, &
  & topological_result, wavevector
use sparse_matrices
use linalg, only: zgetrf
use hamiltonianConstructor, only: ZB8bandQW
use eigensolver, only: ...
use pfaffian, only: complex_pfaffian
use bdg_hamiltonian, only: pairing_sign, pairing_partner
```

No cycle risk: if `bdg_observables` gains a dep on `pfaffian` + `sparse_matrices`, `topological_analysis` does not pull `bdg_observables` back in. The graph stays acyclic.

### 3. Module layer for each dependency candidate

| Module | Current imports | Imported by what could be added to `bdg_observables`? |
|---|---|---|
| `sparse_matrices` (`csr_matrix`) | `definitions` only (`src/math/sparse_matrices.f90:14`) | Yes. L0 leaf. ✓ |
| `pfaffian` (`complex_pfaffian`, `kitaev_majorana_number`) | `definitions` only (`src/math/pfaffian.f90:24`) | Yes. L0 leaf. ✓ |
| `topological_analysis` (`wire_pfaffian_witness_sweep`) | L3 hub (sparse, linalg, hamiltonianConstructor, eigensolver, pfaffian, bdg_hamiltonian) | **No** — would explode `bdg_observables` from 88 lines + `dp` to a L3 hub. ✗ |

**Conclusion:** `bdg_observables` can host the seam siblings by importing only L0 leaves (`pfaffian`, `sparse_matrices`). It must NOT import `topological_analysis`, because that would invert the current "seam = thin" architecture into "seam = hub" without physics benefit.

### 4. Implementation placement

**Decision: slim witness body inline in `bdg_observables.f90`** (per ticket's "no new module" constraint).

The S2 path is self-contained: ~30 lines lifted from `topological_analysis.f90:1724-1764` of `wire_pfaffian_witness_sweep` (omega construction + per-site CSR row-extraction + `complex_pfaffian` call + `best_|Pf|` site pick). No eigensolver, no `hamiltonianConstructor`, no `bdg_hamiltonian` — just `csr_matrix` indexing and Pfaffian arithmetic.

```fortran
! Skeleton (to be implemented at execution step, not in this research ticket):
function eval_bdg_pfaffian_witness_csr(H_bdg_csr, n_full, params) result(s2_sign)
  use sparse_matrices, only: csr_matrix
  use pfaffian,        only: complex_pfaffian
  type(csr_matrix), intent(in)            :: H_bdg_csr
  integer,          intent(in)            :: n_full
  type(bdg_eval_params_t), intent(in)     :: params
  integer                                :: s2_sign  ! {-1, 0, +1}
  ! ... copy of omega_local + per-site idx loop + Pf( H_proj · omega )
end function
```

For the Kitaev sibling, the wrapper is even thinner — `kitaev_majorana_number` already encapsulates the `polar_decomposition_sign` + LAPACK `zheev` call internally (`src/math/pfaffian.f90:120-171`), so the seam wrapper is one-line:

```fortran
function eval_bdg_kitaev_majorana(H_k_array, k_par_values, params) result(majorana_number)
  use pfaffian, only: kitaev_majorana_number
  majorana_number = kitaev_majorana_number(H_k_array, k_par_values)
end function
```

`params` is currently unused in the body; it's part of the API for symmetry (cf. ticket 01's "leaves room for ticket 02's noise-floor sub-decision").

`wire_pfaffian_witness_sweep` (S2-only) and `wire_pfaffian_witness` (S1+S2, full-diag) **stay in `topological_analysis.f90`** — the unit test `test_wire_pfaffian_witness.pf` exercises the dense S1+S2 path directly, and the S1 path requires full diagonalization (`zheev`) which is out-of-scope for the seam. The seam only needs S2 (the swept witness); the dense per-site variant is for unit-test fidelity.

### 5. Why NOT re-export from `bdg_observables` (alternative rejected)

Re-exporting would mean `bdg_observables` adds `use topological_analysis, only: wire_pfaffian_witness_sweep` and `public :: wire_pfaffian_witness_sweep`. Two failure modes:

1. **Dep explosion.** Every importer of `bdg_observables` (currently `main_topology.f90` and `tests/unit/test_bdg_evaluator.pf`) would transitively import `topological_analysis`, which pulls in `hamiltonianConstructor` (798 lines), `eigensolver`, `bdg_hamiltonian`, `pfaffian`, `sparse_matrices`, `linalg`. The seam's "88-line leaf" identity evaporates. The unit test `test_bdg_evaluator.pf` currently builds in seconds because it has a 2-module dep; doubling this for a 30-line function is wrong order of magnitude.

2. **Subroutine vs function mismatch.** `wire_pfaffian_witness_sweep` is a `subroutine` with `intent(out) :: s2_sign`. The seam API (ticket 01) requires a function. Wrapping a subroutine in a function is `s2_sign = 0; call ...; return s2_sign` — error-prone boilerplate for no benefit when the body is already extractable.

### 6. API-shape compatibility with ticket 01

Confirmed compatible. Ticket 01 spec:
```fortran
function eval_bdg_pfaffian_witness_csr(H_bdg_csr, n_full, params) result(s2_sign)
```
takes `csr_matrix` (from `sparse_matrices`) + `n_full: integer` + `bdg_eval_params_t`, returns `integer ∈ {-1, 0, +1}`. The existing `wire_pfaffian_witness_sweep` body uses exactly these inputs and produces exactly this output via `intent(out)`. No interface translation. Same for the Kitaev wrapper.

### 7. Noise-floor tolerance (the fog sub-decision)

The Pfaffian sign-extraction uses `1.0e-12_dp` as a "spurious sign under near-zero Pf" floor — appears **3 times** in `topological_analysis.f90:1663, 1681, 1766`. It's a magic literal repeated across the same module.

This is the same pattern ADR 0008 §1 closed for `bdg_default_near_zero_frac = 0.001_dp` and `bdg_default_min_threshold = 1.0e-10_dp` (`bdg_observables.f90:35-36`). **Promote to a third module-parameter** in `bdg_observables.f90`:

```fortran
real(kind=dp), parameter, public :: bdg_default_pfaffian_floor = 1.0e-12_dp
```

(Existing lines 35-36 provide the precedent and the comment framing: "Module-level defaults (SSOT for the near-zero literals).")

Rationale: the FEAST noise floor is empirically ~1e-5 eV (per ADR 0008 / the existing BLOCKING-EMPIRICAL history at `.remember/recent.md`). 1e-12 eV is the same kind of "spurious-sign at machine precision" guard that needs to be on one knob, not three. If a future ticket needs to widen/narrow it, it's a one-line parameter change instead of a grep-and-replace.

The `params` arg in `eval_bdg_pfaffian_witness_csr` then carries the active threshold (defaulting to `bdg_default_pfaffian_floor` via `bdg_eval_params_with_delta`), so the seam signature stays uniform and future-configurability comes for free.

### 8. Test status — `test_wire_pfaffian_witness.pf`

Read in full; three tests:

| Test | Lines | Status now | Effect of slim-witness routing |
|---|---|---|---|
| `test_wire_pfaffian_witness_smoke` | 60-72 | PASSES (asserts `s ∈ [-1, +1]`) | Unaffected — uses `wire_pfaffian_witness` (dense), not the seam |
| `test_wire_pfaffian_witness_sign_agreement_synthetic` | 74-105 | **FAILS BY DESIGN** (per its own @comment) | Unaffected for the same reason |
| `test_wire_pfaffian_witness_nondiagonal_single_particle` | 107-130 | PASSES (asserts `s ∈ [-1, +1]`) | Unaffected |

**The `@todo` is a Fortran comment, not a pFUnit `@todo` macro.** From `test_wire_pfaffian_witness.pf:96`:

```fortran
! @todo replace synthetic fixture with real wire BdG once U13 lands
```

That's `!` — a comment-only marker. pFUnit does NOT recognize `@todo` outside its own macro grammar. The test is **active and runs**, and the strict assertion `@assertTrue(s1 == s2 .and. s1 /= 0)` (line 104) fails on every run because the synthetic fixture's BdG is diagonal-in-(c,c†) → Pf vanishes identically → s1=s2=0 → strict `s1 /= 0` fails.

**This `@todo` does NOT unblock when the slim witness becomes the gate row.** The failure mode is **at the synthetic-fixture level**, not at the seam level. The synthetic `build_synthetic_bdg_16` (`test_wire_pfaffian_witness.pf:31-58`) constructs a diagonal single-particle block + diagonal pairing block; under this structure no bands mix across the Nambu seam, so Pf=0 by construction. Re-routing through `eval_bdg_pfaffian_witness_csr` doesn't change the matrix — it only changes the calling convention.

**This is a sub-finding for ticket 06** (test coverage audit). Ticket 06 needs to decide one of:
- (a) Write a non-diagonal synthetic fixture (block-off-diagonal coupling between conduction bands + others) to unblock the strict assertion now, without waiting for U13.
- (b) Replace the test with a real wire BdG matrix (deferred to U13).
- (c) Soften the strict assertion to `s1 == s2` (allow either both-zero or both-nonzero) — rejects ticket 05's "no zero-escape hatch" discipline per Phase 24 follow-up A.

Recommendation: **option (a)** — non-diagonal synthetic fixture (sub-30 LOC). Decouples the test failure from U13's wire periodic/Bloch dependency and lands the strict assertion as part of this PR, not the next.

### 9. Summary of decisions landed

1. **Slim witness body**: inline S2 extraction in `bdg_observables.f90` (~30 LOC, lifted verbatim from `topological_analysis.f90:1724-1764`).
2. **Dependencies**: `bdg_observables` gains `use sparse_matrices` + `use pfaffian` (both L0). No `topological_analysis` import.
3. **No new module** (ticket constraint honored).
4. **No cycle risk** (verified: `topological_analysis` does not import `bdg_observables`).
5. **API compatibility**: confirmed against ticket 01 — no interface changes needed.
6. **Kitaev wrapper**: thin one-line facade in `bdg_observables.f90` over `pfaffian.kitaev_majorana_number`.
7. **Noise floor**: promote `1.0e-12_dp` literal to `bdg_default_pfaffian_floor` parameter in `bdg_observables.f90` (3rd SSOT, precedent set by `bdg_default_near_zero_frac` + `bdg_default_min_threshold`).
8. **`wire_pfaffian_witness_sweep` and `wire_pfaffian_witness`**: **stay** in `topological_analysis.f90`. The seam mirrors the S2 path; the dense S1+S2 path stays for unit-test fidelity.
9. **Test `@todo` sub-finding for ticket 06**: the strict sign-agreement test stays `@todo` regardless of seam routing; recommended path is non-diagonal synthetic fixture (option a, ~30 LOC) to unblock it now.

### 10. Doc-drift follow-up for ticket 07

`src/physics/AGENTS.md:46-58` dep DAG is missing `bdg_observables` (currently an L0 leaf, will become L1 with `pfaffian` + `sparse_matrices`). AGENTS.md inventory table line 44 lists the module but the DAG block doesn't reflect its layer. Ticket 07's plan close-out sweep should patch this as a 1-line DAG entry.