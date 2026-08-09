# T6 — Resolve loose fog L1–L3 before claiming T1

> Child of [.scratch/bdg-u13-bloch-pfaffian/map.md](../map.md). Wayfinder
> ticket; gates the build branch of this chart (T1–T4). Materialized from
> the map's T6 entry by the charting session that finished the chart.

Type: grilling
Status: resolved

## Question

Three loose-end questions (L1–L3, recorded in the map's §Loose fog) hang over
T1 and must be **decided** — not researched — before anyone claims the
Bloch-periodic builder. Each is a convention choice this codebase has not yet
made; resolving them fixes the contract T1 builds against and the semantics T2's
agreement seam asserts. Grill them one at a time, in order, with the human who
speaks for the BdG physics:

1. **L1 — By-Peierls scope.** The Bx-only Peierls guard at `defs.f90:957-968`
   rejects `Bx=0 ∧ (By∨Bz≠0)` ("By Peierls path not yet implemented"). U13's
   Bloch-periodic `H_BdG(k_par)` must keep this exact gate so it cannot silently
   evade it. **Decide:** does U13 (a) preserve Bx-only and document `By/Bz` as a
   separate follow-up effort (scope boundary on this PR), or (b) widen the
   Peierls path to By/Bz as part of U13? The map's fog leans (a) — U13 is about
   *k*-periodicity, not *transverse-field* Peierls — but the human must confirm
   the scope cut so T1 knows what to reject.

2. **L2 — BZ-slice lattice.** A Bloch-periodic builder needs a 1D Brillouin-zone
   slice: a `k_par` lattice (count `n_k`, spacing, BZ boundaries) for the stack
   `kitaev_majorana_number` consumes. **Decide:** is the wire's free-`z`
   direction the `k_par` axis (so `k_par` ∈ `[-π/L_z_unit, +π/L_z_unit]` over the
   *wire supercell* period), and does U13 (a) sample a fixed uniform grid from
   the new `nk_par` config knob (T3's lever), or (b) reuse the existing
   `[wave_vector]` sweep machinery? The answer fixes T1's stack shape and T3's
   knob surface.

3. **L3 — Closure-disagreement semantics.** `kitaev_majorana_number` returns
   `{-1, 0, +1}`, and the slim S2 seam returns `s2 ∈ {-1, 0, +1}`. When S1 (the
   k-product) and S2 (the projected single-point Pfaffian) **disagree at
   closure** (one is 0, the other ±1, or signs differ at a gap closure), the
   T2 agreement seam must choose what `z2` to emit. **Decide:** the
   disagreement contract — (a) `z2 = 0` (treat closure-disagreement as
   "indeterminate/gap-touching," the honest answer at a topological
   transition), (b) `z2 = S1` (trust the k-product as the strict invariant, S2
   is only approximate), or (c) `error stop` (no silent disagreement). This
   fixes T2's branch and T4's gate predicate. Lean toward (a): it matches the
   native `{-1,+1,0}` convention U10 T3 already locked for the wire path and
   keeps the gate honest at the phase boundary the whole U13 effort cares about.

Resolve each by recording the chosen branch with a one-line **why** in the
Answer. Outcomes T1 (scope), T2 (BZ slice shape + knob), and T3 (disagreement
branch) build directly from these decisions.

## Answer

### L1 — By-Peierls scope: **(a) preserve Bx-only; By/Bz Peierls is out of scope for this effort.**

**Why:** U13's destination is k-periodicity for S1×S2 strict agreement at μ≈0.6601,
not transverse-field Peierls generality. The existing `add_peierls_coo`
(magnetic_field.f90:79) is gated on `abs(Bx) > 1e-12_dp` at both call sites in
`bdg_hamiltonian.f90:274-277` and `:421-425` — By/Bz paths are not even attempted
in-tree. Lifting that gate would require extending the Peierls line integral to
the y/z directions, new fixtures, and loosening `validate_semantic` at
`defs.f90:957-968` — none of which advances the U13 destination. Conflating
k-periodicity with transverse-field Peierls would mix two orthogonal extensions
and explode gate failure modes. The hard `error stop` at `defs.f90:957-968` is
the correct upstream guard; U13 inherits it.

**Trade-off documented:** any future fixture with `Bx=0 ∧ (By∨Bz≠0)` will hit
`defs.f90:957-968` and `error stop`. Lifting that guard is a **fresh effort**
(separate destination, separate chart) — not part of U13. Logged in
`map.md ## Out of scope`.

### L2 — BZ-slice lattice: **2A k_par = wire free-z; 2B new `[bdg]` `nk_par` knob (+ optional `k_par_min`/`k_par_max`), NOT reuse `[wave_vector]`.**

**Why 2A:** U13's destination is the lift from single-kz (U10's `wire_pfaffian_witness_sweep`) to a closed BZ path for S1's k-product. The wire's free-z is the only direction the BdG supercell is periodic in — the transverse xy are real-space-confined, not k-space. Any other choice either collapses to U10's single point or invents fictitious periodicity.

**Why 2B:** `[wave_vector]` drives the *eigensolve* band-sweep for bandStructure — semantically distinct from the *witness* k-stack `kitaev_majorana_number` consumes. Reusing `[wave_vector]` would couple two unrelated physics pipelines through one knob, violating the "optional sections enabled by presence" invariant (CLAUDE.md) and creating the silent-evasion failure mode the gate exists to prevent. A dedicated `nk_par` knob (with `k_par_min`/`k_par_max` defaulted to `[-π/L_z_unit, +π/L_z_unit]` if absent) keeps the BdG witness pipeline isolated and explicit.

**Defaults:** uniform sampling on the `nk_par` grid; `nk_par=1` reproduces the fixed-kz builder bit-for-bit (T3's regression guard). Old configs without `[bdg]` continue to work — U10's single-point witness path is the unchanged fallback.

**Trade-off documented:** adds one config-validation rule to `validate_semantic()` for the `topologicalAnalysis` app (presence + sign of `nk_par`; range bounds on `k_par_min`/`k_par_max`). Small cost; right surface (`[topology]` and `[bdg]` already coexist).

### L3 — Closure-disagreement semantics: **(a) `z2 = 0` on disagreement; `disagreement_reason` WARN logged to witness file; case (iv) reserved for a possible follow-up `error stop` guard.**

**Why (a):** Only honest answer at a phase boundary. The destination exists precisely to detect closures — collapsing `z2 = 0` to `z2 = S1` (option b) re-imposes the S1-vs-S2 asymmetry U10 tried to escape. `error stop` (option c) would crash `topologicalAnalysis` at every FEST-saturation point (case v) and gap closure (case vi), defeating the destination.

**Disagreement cases the seam must handle** (S1 = k-product invariant, S2 = projected Pfaffian):
- i (S1 = S2 = ±1): unambiguous strict trivial/topological — `z2 = S1`.
- ii (S1 = S2 = 0): both indeterminate (gap closure + FEST saturation) — `z2 = 0`, no WARN.
- iii (S1 ≠ S2 in sign): strict disagreement — physically impossible, signals sampling/build bug — `z2 = 0` + WARN `disagreement_reason = strict_sign_split`.
- iv (S1 = ±1, S2 = 0): S1 clean, S2 at FEST floor — `z2 = 0` + WARN `disagreement_reason = s2_fest_saturation`.
- v (S1 = 0, S2 = ±1): S1 detects closure, S2 returns clean sign at one kz — `z2 = 0` + WARN `disagreement_reason = s1_closure_s2_misses`.

**Case (iv) reservation:** strict sign-split (case iii) is the only physically-impossible disagreement and could justify a follow-up `error stop` guard. Deferred to a post-T1 ticket once empirical case (iii) data exists (noise vs build bug).

**Inheritance:** (a) matches U10 T3 native `{-1, +1, 0}` convention (`z2 = 0` = "no strict answer"). No new convention; no churn for gate predicate or lecture-13 acceptance script.

**Trade-off documented:** witness file gains a `disagreement_reason` column with five values (`none`, `s1_s2_sign_split`, `s2_fest_saturation`, `s1_closure_s2_misses`, plus `none` for the agreement cases). T3's writer follows the existing `wire_pfaffian_witness_sweep` pattern; T2's seam updates `bdg_observables.f90:6-13` call-site header if seam-face count changes.

---

## Resolution

All three loose-fog decisions (L1–L3) locked. T6 closes; T1 is now unblocked.

- **L1 (scope):** preserve Bx-only; By/Bz Peierls → `map.md ## Out of scope`.
- **L2 (k-grid):** k_par = wire free-z; new `[bdg]` `nk_par` knob (defaults preserve U10 bit-for-bit).
- **L3 (disagreement):** `z2 = 0` on disagreement with `disagreement_reason` WARN; case (iii) strict-sign-split deferred to follow-up `error stop` guard.
