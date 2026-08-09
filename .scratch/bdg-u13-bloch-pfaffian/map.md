# Wayfinder U13 Map — Bloch-periodic BdG → full Pfaffian (S1×S2)

> Map file (wayfinder:map). Lives at `.scratch/bdg-u13-bloch-pfaffian/map.md` per
> repo convention (local-markdown tracker, `.scratch` is git-tracked — maps get
> committed). Parent: Unit U13 of
> `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` (L53-54,
> "Explicitly deferred: U13 (periodic/Bloch BdG construction for the Majorana
> number; without it the wire Pfaffian sweep evaluates at one point only)").

## Charting session started 2026-08-14

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

- **L1 — By-orbital Peierls:** `defs.f90:957-968` gates on Bx-only; the By
  Peierls path is "not yet implemented." Does U13's Bloch builder need By
  coupling for the (B, μ) phase-diagram physics the gate checks, or is the
  transverse-Bx axial-wire regime (the one already validated U2/U10) the
  intended physical scope? Deciding L1 bounds whether a Peierls-twist builder
  is in U13 or deferred past it.
- **L2 — k_par lattice / BZ slice:** is the 1D BZ the wire-axis `kz ∈ [−π/R,
  π/R]` (R = transverse wire-unit pitch), or does the Bloch slice mean
  something different for the Rashba-wire geometry than the QW? `ZB8bandGeneral
  ized` already parametrizes on `kz`; confirm the Bloch lattice vector maps to
  that same `kz` axis (not a separate `k_par`) so the builder reuses the
  existing k-dependent path without a new k-mesh subsystem.
- **L3 — S1/S2 disagreement regime:** at the gap-closure curve, can S1 and S2
  disagree by construction (S1 k-product zero, S2 single-point sign) such that
  strict equality is too strong? If yes, the agreement seam needs a documented
  closure-disagreement -> z2=0 branch (not a defect), matching U2's User Story
  5 "structurally valid gap-closure signal" wording (parent plan L97-99).

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

## Frontier tickets (TBD — enumeration refined 2026-08-14)

This is the *initial* chart; tickets are not yet claimed. Refined from the
proposal that preceded the source dive (the "vendor SKBPFA" research ticket is
collapsed to a verification since `zskpfa_reduction` is in-tree):

- **T1** — Bloch-periodic H_BdG(k_par) builder (generalize
  `build_bdg_hamiltonian_1d` at fixed kz → k_par stack; reuse
  `ZB8bandGeneralized`; preserve Bx-only gate).
- **T2** — S1×S2 agreement seam `eval_bdg_pfaffian_witness_product_csr`
  (wire `kitaev_majorana_number` over T1's stack × existing slim S2; native
  `{-1,+1,0}`; closure-disagreement branch per L3).
- **T3** — Config-driven `nk_par` / `bloch_sweep` knob + dispatch in
  `main_topology.f90`; `nk_par=1` default reproduces fixed-kz bit-for-bit
  (scope isolation regression guard).
- **T4** — Strict invariant → 4-witness gate: regenerate (B, μ) diagram with
  z2 = strict S1×S2 product; flip `lecture_13_topological.py` 4-witness label
  from "approximate" to "strict"; resolve `@todo U13` markers.
- **T5** — TDD-red test asserting strict agreement at μ≈0.6601 (fails on main
  "one point only") → green after T1-T4.
- **T6** — Grill/loose resolution L1-L3 before claiming T1 (By-Peierls scope,
  BZ-slice lattice, closure-disagreement semantics).

Order: T6 (grill) unblocks T1; T1 → T2 → T3 → T4; T5 is the red→green gate that
straddles T2. No research ticket — the algorithm is in-tree.
