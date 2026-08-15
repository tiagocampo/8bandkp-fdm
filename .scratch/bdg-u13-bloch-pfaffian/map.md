# Wayfinder U13 Map — Bloch-periodic BdG → full Pfaffian (S1×S2)

> Map file (wayfinder:map). Lives at `.scratch/bdg-u13-bloch-pfaffian/map.md` per
> repo convention (local-markdown tracker, `.scratch` is git-tracked — maps get
> committed). Parent: Unit U13 of
> `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` (L53-54,
> "Explicitly deferred: U13 (periodic/Bloch BdG construction for the Majorana
> number; without it the wire Pfaffian sweep evaluates at one point only)").

## Charting session started 2026-08-14; tickets materialized 2026-08-09;
T2 + T5 closed 2026-08-09 (in same session, picked up after T2 handoff)

Source-grounded after landing U10 (PR #43 merged to `main@52743b3`, 2026-08-08).
U10 closed the slim S2 projected Pfaffian seam; the 4-witness acceptance gate's
`bcrit_pfaffian` row stayed **approximate-not-strict** (per-B `min-|Pf|` proxy)
because S1 (the k-product Majorana number) needs a Bloch-periodic H_BdG(k)
that the fixed-kz builder cannot produce. U13 supplies it, lifting the gate
from approximate to strict S1×S2 agreement.

## State (verified this session)

What already exists in-tree — this is *not* a green-field chart:

- **Pfaffian algorithms (avoid vendoring):** `src/math/pfaffian.f90` has
  `skpfa_reduction` (real, L291) + `zskpfa_reduction` (complex, L384),
  Parlett-Reid/Wimmer (arXiv:1102.3440), O(n³) for n>12, plus Laplace
  expansion `skpfa_laplace`/`zskpfa_laplace` for n≤~12. The Kitaev convention
  `Pf([[0,a],[-a,0]])=+a` is already the module's convention (header L17).
  → A "vendor SKBPFA" research ticket is a **no-op**: the algorithm is present.
- **S1 template (k-product Majorana number):**
  `kitaev_majorana_number(H_k_array, k_par_values, omega_struct)` at
  `pfaffian.f90:120-171` — takes a stack `H_k_array(:,:,n_k)` of Hermitian
  H_BdG(k) over `k_par_values`, polar-decomposes each, and returns
  sign(det U_odd)-product ∈ {-1, 0, +1}. This is the S1 surface U13 wires.
- **S2 seam (already shipped):** `eval_bdg_pfaffian_witness_csr` (
  `bdg_observables.f90:201`) delegates to `wire_pfaffian_witness_sweep`
  (`topological_analysis.f90:1650`) — single-point slim projected Pfaffian,
  s2 ∈ {-1,0,+1}, with optional `best_pf_abs` out-arg (U10 T1a). The forward
  reference at `bdg_observables.f90:195-199` names U13/Issue 05 explicitly.
- **Fixed-kz BdG builder:** `build_bdg_hamiltonian_1d` (`bdg_hamiltonian.f90:142`)
  builds H_BdG at a **single** `kz` (passes ±kz to `ZB8bandGeneralized`,
  L173/L182), returns CSR, Bx-only Peierls. Called at `main_topology.f90:480,
  718, 1346`. This is the builder whose output the wire Pfaffian sweep consumes
  at one point — exactly the "one point only" deficiency parent plan L54 names.
- **Bx-only Peierls guard:** `defs.f90:957-968` already rejects Bx=0 ∧ (By∨Bz
  ≠0) with "By Peierls path not yet implemented." The By-orbital Peierls path
  is an **open loose** (see §Loose) — U13's Bloch builder must keep the same
  Bx-only gate so it cannot silently evade it.
- **Acceptance gate (target surface):** the existing 3-witness gate reads the
  `z2` column from `output/z2_phase_diagram.dat` at `mu ≈ 0.6601` (z2==−1 first
  B-row). U10's 4-witness `bcrit_pfaffian` slot consumes the S2 `min-|Pf|`
  proxy; U13 promotes it to the S1×S2 agreement product so the 4th witness is
  strict, not approximate.

## Destination

Close U13 end-to-end by making the wire Pfaffian sweep evaluate across a 1D
Brillouin-zone slice (not one point) and emit the **full Pfaffian invariant**
S1×S2, so the 4-witness acceptance gate's `bcrit_pfaffian` row is strict-valid:

- A **Bloch-periodic builder** `build_bdg_hamiltonian_1d_bloch` (or a sweep-mode
  flag on the existing builder) that returns a stack `H_k_array(:,:,n_k)` of
  Hermitian H_BdG(k_par) over `k_par_values ∈ [−π/R, π/R]` for the wire
  direction. Reuses `ZB8bandGeneralized` at each k_par (the fixed-kz builder
  already calls it at ±kz; the Bloch variant generalizes kz → k_par sweep).
  The Peierls Bx-only gate (`defs.f90:963`) is preserved verbatim — the Bloch
  builder inherits it, no silent By/axial escape.
- A **S1×S2 agreement seam** `eval_bdg_pfaffian_witness_product_csr` that
  (a) calls `kitaev_majorana_number` on the k-sweep stack for S1, (b) calls the
  existing slim `wire_pfaffian_witness_sweep` per-point and reduces for S2, and
  (c) returns a single strict invariant: `agree == s1 == s2` with the
  {-1,+1,0} native convention. This replaces the `bcrit_pfaffian` proxy slot's
  approximate input.
- A regenerated (B, μ) phase diagram where the `z2` column is the strict
  S1×S2 product — open→close→reopen visible as z2 ∈ {+1, 0, −1} across B at
  μ ≈ 0.6601, the closure point's |Pf|-min now the *true* band-minimum (not the
  per-site max-projection proxy).
- The `lecture_13_topological.py` 4-witness row label switches from
  "approximate (U13 deferred)" to "strict S1×S2".
- The `@todo U13` markers at `test_wire_pfaffian_witness.pf:96` and the
  `bdg_observables.f90:195-199` forward reference are resolved into closure
  notes (GREEN via real Bloch S1, not via the synthetic-fixture escape hatch).

A "shipped" U13 turns the deferred strict-sign-agreement problem into a live
strict invariant the gate can check.

## Notes (standing decisions — lock before execution)

- **Scope isolation:** the Bloch builder must not disturb the **existing
  fixed-kz** path that the validated (B, μ) diagrams depend on. Bloch is an
  *additional mode*, gated on a config knob (e.g. `cfg%bdg%bloch_sweep` /
  `nk_par`) — the default `nk_par=1` at the existing `cfg%bdg%kz` reproduces
  today's fixed-point behavior bit-for-bit (regression guard). No change to
  `write_z2_phase_diagram`'s convention header beyond the agreement-product
  semantics.
- **Algorithm reuse (do not reinvent):** use `zskpfa_reduction` /
  `kitaev_majorana_number` from `pfaffian.f90` for all Pfaffian work. The
  earlier proposal to vendor an external SKBPFA routine is rejected — the
  in-tree module already implements it (Wimmer, arXiv:1102.3440). Any U13
  ticket proposing new skew-symmetric decomposition code must first justify why
  the existing `pfaffian.f90` routines are insufficient.
- **Native convention preserved:** the Bloch path writes Pfaffian-native
  `{-1,+1,0}` (∈ topological, trivial, closure) directly — no remap. Non-wire
  paths keep their existing convention (U10 standing decision, T3).
- **TDD entry:** the first execution ticket is **red-first**: write a test
  asserting strict S1×S2 agreement at μ ≈ 0.6601 across B — it must **fail**
  on current `main` with "wire Pfaffian evaluates at one point only" before any
  builder change. (Mirrors U10 T3's TDD-red entry.)
- **PR shape:** one scoped PR, `feat/bdg-u13-bloch-pfaffian`, branched from
  `main@52743b3`. Per CLAUDE.md Known Issues / parent plan L55 ("separate
  scoped PR") and the established U2/U10 pattern.

## Loose (open questions to resolve in a grill/chart sub-session before claiming a ticket)

<!-- L1–L3 resolved 2026-08-09 by T6 (closed). See Decisions so far for the
     locked answers; By/Bz Peierls logged in Out of scope. No loose questions
     remain before claiming T1. -->

_(L1–L3 resolved by T6 on 2026-08-09 — see [Decisions so far](#decisions-so-far). No loose questions remain before claiming T1.)_

## Files to Read first (when claiming a ticket)

1. **This map** (`.scratch/bdg-u13-bloch-pfaffian/map.md`).
2. `src/physics/bdg_hamiltonian.f90` §`build_bdg_hamiltonian_1d` (L142-443) —
   the fixed-kz builder the Bloch variant generalizes.
3. `src/math/pfaffian.f90` §`kitaev_majorana_number` (L120-171) — the S1
   k-product template, and §`zskpfa_reduction` (L384-474) — the complex Pfaffian.
4. `src/physics/bdg_observables.f90:195-201` — the S2 seam + the U13 forward
   reference being closed.
5. `src/core/defs.f90:957-968` — the Bx-only Peierls guard to preserve.
6. `.scratch/archive/bdg-u10-pfaffian-phase-diagram/HANDOFF.md` — the U10 close
   and the "deferred to U13" destination text; `tickets/T3_native_z2.md` for the
   native-convention standing decision.
7. `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` L53-55 —
   the parent spec sentence U13 exists to satisfy.

## Decisions so far

<!-- the index — one line per closed ticket: enough to judge relevance, then zoom the link -->

- [Resolve loose fog L1–L3 (T6)](issues/06-grill-loose-fog.md) — L1 preserve
  Bx-only Peierls (By/Bz deferred); L2 k_par = wire free-z + new `[bdg]`
  `nk_par` knob (defaults preserve U10 bit-for-bit); L3 disagreement → `z2 = 0`
  with `disagreement_reason` WARN logged; case (iii) strict-sign-split reserved
  for a follow-up `error stop` guard. T1 unblocked.
- [Bloch-periodic H_BdG(k_par) builder (T1)](issues/01-bloch-periodic-bdg-builder.md) —
  `build_bdg_hamiltonian_1d_bloch(H_k_array(:,:,n_k), n_k, k_par_values, ...)` in
  `src/physics/bdg_hamiltonian.f90:454-563`. Per-slice CSR via existing wire builder
  + `csr_to_dense_work` into the stack; fresh `wire_workspace` per slice. n_k=1
  reproduces fixed-kz CSR bit-for-bit (1e-12 tol). Bx-only Peierls guard mirrors
  `defs.f90:963-969` SSOT as defense-in-depth. 52/52 unit + 5/5 wire-BdG
  regression + 1/1 lecture-13 gate green. Commit `302e092`.
- [S1×S2 agreement seam `eval_bdg_pfaffian_witness_product_csr` (T2)](issues/02-s1-s2-agreement-seam.md) —
  public seam in `bdg_observables.f90:283-347` + private `map_s1_s2_to_z2` L3
  mapper (L361-389) covering all 5 L3 cases. Reads S1 via `kitaev_majorana_number`
  over T1's Bloch stack; S2 via `eval_bdg_pfaffian_witness_csr` on slice 1
  dense→CSR scratch. Native `{-1,+1,0}`. 53/53 unit + 5/5 wire-BdG regression
  + 1/1 lecture-13 gate green. Commit `449c7fb`.
- [TDD-red strict S1×S2 agreement test (T5)](issues/05-tdd-strict-agreement-test.md) —
  extended `tests/unit/test_bdg_pfaffian_witness_product_csr.pf` with 3 new
  `@test` subroutines for cases (ii) S1=S2=0, (iii) S1≠S2 sign-split,
  (v) S1=0 S2=±1. TDD red→green cycle executed (pre-T2 build error, HEAD green).
  Resolved T2's "4×4 n_odd=1 fixtures" forward-pointer to 16×16 stacks via
  `n_k=1` guard for (ii)/(v) + band-7/8 coupling sign for (iii). 53/53 unit
  + 3/3 wire-BdG regression + 1/1 lecture-13 gate green. Commit `f6498e6`.
- [Config-driven `nk_par` knob + dispatch + bit-for-bit regression guard (T3)](issues/03-bloch-sweep-config-knob.md) —
  `bdg_config` gains `nk_par` + `k_par_min`/`k_par_max` (optional, defaults preserve U10).
  `validate_semantic` rejects `nk_par<1` and `nk_par>1 ∧ k_par_max<=k_par_min` with
  `error stop`. `compute_wire_bdg_gap_sweep` branches on `nk_par`: `<=1` → existing
  `eval_wire_bdg_gap` (bit-exact U10); `>1` → new `eval_wire_bdg_gap_bloch` (T1 stack
  + T2 strict seam + `dense_to_csr` extraction in `sparse_matrices.f90`). Two
  regression guards: `regression_bdg_u13_nk_par_equiv` (byte-equivalence nk_par=1
  vs absent) + `regression_bdg_u13_nk_par_sweep` (nk_par=4 smoke near Gamma).
  53/53 unit + 4/4 wire-BdG regression + 1/1 lecture-13 gate green. Commits
  `e17fe11`, `5896455`, `9336dac`, `f4af568`, `3db44cc`. Side-finding: latency UB
  in `test_green_functions.pf` (unallocated `cfg%params`) fixed at `5896455`.
- **T4 review checkpoint (2026-08-15):** the full odd-sector determinant in
  `kitaev_majorana_number` was corrected from the production-incorrect leading
  2x2 shortcut to a full LU determinant, with a red-first 16x16 regression.
  After relinking `topologicalAnalysis`, canonical `nk_par=4` probes at B=0
  and 5 T both return native `z2=0`, `disagreement_reason=1` (S1/S2 sign
  split). T4 therefore remains claimed/open; the strict gate is not yet
  certified.

## Out of scope

<!-- scope, not sharpness, lands here — closed tickets whose answer is "this
     isn't on the route to the destination." Returns only if the destination
     is redrawn. -->

- **By/Bz Peierls (transverse-field generality)** — `defs.f90:957-968` keeps
  its `error stop` on `Bx=0 ∧ (By∨Bz≠0)`. Lifting the guard (extending
  `add_peierls_coo` to y/z directions, new fixtures, loosening `validate_semantic`)
  is a fresh effort — not part of U13. Confirmed via T6/L1.

## Open tickets

<!-- the map is an index; open tickets live as child files under issues/, the -->
<!-- frontier is the open, unblocked, unclaimed ones (first by number wins). -->

Five child tickets closed (T6, T1, T2, T5, T3 closed 2026-08-09). One
ticket remains on the frontier:
Untitled: T4 → [04](issues/04-strict-4witness-gate.md) — T4 · task · **FRONTIER NOW** —
  the **destination**: strict 4-witness gate + 4-witness label flip +
  `@todo U13` resolution in `lecture_13_topological.py`. Regenerated
  phase diagram from the strict S1×S2 product seam (T3's `nk_par>1`
  branch) is the input.

Build order: T4 → destination. T3 closed 2026-08-09 (5 commits, all
summary at T3's ticket body). T5 closed ahead of T3 (seam-test work
parallelizable; the chart's T5 ∥ T3 fork is now T4-alone).
