# T1 — Bloch-periodic H_BdG(k_par) builder

> Child of [.scratch/bdg-u13-bloch-pfaffian/map.md](../map.md). Build branch of
> this chart; unblocks T2. Materialized from the map's T1 entry.

Type: task
Status: completed
Owner: next-session (T1 chart handoff 2026-08-09) — landed on `feat/bdg-u13-bloch-pfaffian` (commits `302e092` + `1a14489`)
Blocked by: 06 (closed 2026-08-09 — T6)

## Resolution (2026-08-09)

New public subroutine `build_bdg_hamiltonian_1d_bloch` added to
`src/physics/bdg_hamiltonian.f90`. Private `dispatch_bdg_wire_builder`
helper routes the (B_vec, g_factor) optional pair (DRY/OCP). n_k=1
reproduces the fixed-kz builder bit-for-bit (regression guard pinned by
`test_bloch_n1_matches_fixed_kz` + `test_bloch_n1_matches_fixed_kz_with_bx`,
tolerance 1e-12 no-B / 1e-10 Bx-only). Bx-only Peierls guard mirrors
`defs.f90:963-969` SSOT as defense-in-depth at builder entry.

`H_k_array` is `allocatable, intent(out)` — caller no longer needs to
know Ntot ahead of time; Ntot derives from `cfg%grid%npoints()` (SSOT).

Test fixture dedup: `setup_wire_gas_fixture` + `teardown_wire_gas_fixture`
extract ~45 lines of @test-local duplication (DRY/SSOT). Per-@test body
44-52 lines, within the 50-line budget.

Verification: 52/52 unit (51 baseline + 1 T1 target with 4 `@test`s —
3 prior + 1 added for Bx parity), 5/5 wire-BdG regression, 1/1
lecture-13 acceptance gate.

## Review fixes (commits `1a14489`, post code-review)

Code review at `302e092` surfaced 5 hard violations + several smells.
Applied fixes:

- **S1**: `H_k_array` non-allocatable `intent(out)` → `allocatable, intent(out)`.
  Caller no longer needs to know Ntot.
- **S2 / Repeated Switches**: extracted `dispatch_bdg_wire_builder` private
  helper. Three-way nested if cascade replaced by single call.
- **S3 / DRY**: extracted fixture helpers, ~45 lines dedup.
- **S6 / Speculative Generality**: dropped unused `ws=` parameter.
- **S7 / Mysterious Name**: default-init `ws_slice = wire_workspace()`.
- **S8 / Assertion Roulette**: tightened third test to Hermiticity check.
- **P1**: added `test_bloch_n1_matches_fixed_kz_with_bx` — pins the
  Bx-only Peierls branch bit-for-bit (was only covered at the
  Hermiticity level).
- **W1**: per-slice `dense_slice` allocate + deallocate — eliminates
  leak-on-early-error robustness gap.
- **W2**: corrected "private from sparse_matrices" — `csr_to_dense_work`
  is public, accessed via `use sparse_matrices`.

**Out of scope (deferred to future effort)**: `bdg_hamiltonian.f90` module
size 703→723 lines (still above 300-line file budget). The module-level
SRP seam split is a parent-plan refactor, not T1 scope; flagged for
future effort.

T1 unblocks T2 (S1×S2 strict agreement seam).

## Question

Generalize `build_bdg_hamiltonian_1d` (`src/physics/bdg_hamiltonian.f90:142`)
from a **single** fixed-`kz` CSR into a **Bloch-periodic stack** over a 1D
`k_par` lattice, so the wire Pfaffian sweep can evaluate the Majorana number
across the Brillouin-zone slice instead of at "one point only" (the deficiency
the parent plan at `docs/plans/2026-06-14-001-...-plan.md` L53-54 names).

Concrete contract a session should nail down before coding (drawing on the
map's §State and the L2 decision from T6):

- **Stack shape.** Produce `H_k_array(:, :, n_k)` — a stack of `n_k` Hermitian
  `H_BdG(k_par)` matrices over `k_par_values(1:n_k)` — in exactly the layout
  `kitaev_majorana_number(H_k_array, k_par_values, omega_struct)`
  (`src/math/pfaffian.f90:120-171`) consumes. Keep the existing single-`kz`
  return path as the `n_k = 1` case so T3's `nk_par=1` default reproduces the
  fixed-kz builder bit-for-bit (scope-isolation regression guard).
- **Reuse, don't rebuild.** Each slice reuses `ZB8bandGeneralized`
  (`hamiltonian_wire.f90`, already called at `bdg_hamiltonian.f90:173/182` for
  `+kz` / `-k`) for the 8-band block and the BdG Nambu doubling. The new work is
  the *k-loop + stack assembly + periodic-boundary Peierls phase*, not a
  second Hamiltonian constructor. Preserve the canonical `-conjg(H0(-k))` form
  already in the builder.
- **Bx-only gate.** Carry the `defs.f90:957-968` Bx-only Peierls guard forward
  verbatim — whatever T6's L1 decision lands as the By/Bz scope boundary, the
  Bloch builder must **not** silently widen the Peierls path the fixed-kz
  builder rejects.
- **Finalizer/workspace hygiene.** Adhere to the global invariants in
  `CLAUDE.md`: `g='g3'` derivative builds stay isolated from
  `wire_workspace` cache; `feast_workspace` reuse stays pattern-validated;
  any new large stack allocator gets `private` default + explicit `public ::`
  exports + a finalizer (delegating to `*_free`) per the allocatable-component
  rule. The `-conjg(H0(-k))` Nambu form and `ZB8bandGeneralized` reuse mean
  a workspace-per-`k_par` cache, not a free-standing one — match the existing
  builder's `ws=` plumbing.

Approach is TDD (`superpowers:test-driven-development`): write a unit test that
 asserts the `n_k=1` stack reproduces the fixed-kz CSR bit-for-bit **red**,
generalize until green. `build_bdg_hamiltonian_1d` is named in `CLAUDE.md`
Boundaries? — **no**: that guard is on `bdg_hamiltonian.f90` Hamiltonian
construction for *Moore/SC-loop*, and BdG sits under `[bdg]` config. But the
file is in `src/physics/`, so read `src/physics/AGENTS.md` before modifying.

Resolved L2 (from T6) fixes `n_k`, BZ boundaries, and whether `k_par` reuses
`[wave_vector]` or a new `nk_par` knob (that's T3's job to wire, T1's job to
consume). Output: the builder + its `n_k=1` regression test, committed on a
`feat/bdg-u13-bloch-pfaffian` branch.

## Notes for the claiming session

- Destination row of this map: lift the 4-witness acceptance gate from
  approximate (per-B `min-|Pf|` proxy, U10) to **strict S1×S2 agreement** at
  μ≈0.6601 across B. T1 is the *foundation* — without the stack, S1 has
  no input.
- A `/research` subagent is **not** needed: `kitaev_majorana_number` and
  `zskpfa_reduction` are already in-tree (map §State). This is building, not
  researching.
