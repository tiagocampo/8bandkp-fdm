# T1 — Bloch-periodic H_BdG(k_par) builder

> Child of [.scratch/bdg-u13-bloch-pfaffian/map.md](../map.md). Build branch of
> this chart; unblocks T2. Materialized from the map's T1 entry.

Type: task
Status: claimed
Owner: next-session (T1 chart handoff 2026-08-09)
Blocked by: 06 (closed 2026-08-09 — T6)

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
