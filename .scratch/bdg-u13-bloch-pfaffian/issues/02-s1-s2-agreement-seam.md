# T2 — S1×S2 agreement seam `eval_bdg_pfaffian_witness_product_csr`

> Child of [.scratch/bdg-u13-bloch-pfaffian/map.md](../map.md). Build branch;
> unblocks T3. Materialized from the map's T2 entry.

Type: task
Status: pending
Blocked by: 01

## Question

Wire the **strict S1×S2 agreement seam** that the destination row of this map
lifts the acceptance gate to. The existing slim seam
`eval_bdg_pfaffian_witness_csr` (`src/physics/bdg_observables.f90:201`)
delegates to `wire_pfaffian_witness_sweep` and returns the **S2** single-point
projected Pfaffian sign `s2 ∈ {-1, 0, +1}` (U10 shipped this). The new seam feeds
*S1* — the k-product Majorana number — by calling
`kitaev_majorana_number(H_k_array, k_par_values, omega_struct)`
(`src/math/pfaffian.f90:120-171`) over T1's Bloch-periodic stack, and produces
the strict product. Native `{-1, +1, 0}` (the convention U10 T3 locked for the
wire path; do not remap to BHZ `{0,1}`.

Contract a session should nail down (drawing on the map's §State and T6's L3
decision):

- **Name + signature.** New public function in `bdg_observables.f90` —
  `eval_bdg_pfaffian_witness_product_csr` — sibling to the existing
  `eval_bdg_pfaffian_witness_csr`. It takes T1's stack (or the inputs to build
  it) plus the existing CSR/params, returns the **strict** `s1 × s2` product as
  `{-1, 0, +1}`. Mirror the existing seam's `best_pf_abs` out-arg plumbing (the
  U10 T1a optional out-arg) so per-B `min-|Pf|` is still available to the gate
  for the degeneracy-detection branch U10 T6 added.
- **Closure-disagreement branch.** When S1 and S2 disagree at closure (one is 0,
  the other ±1, or signs differ at a gap closure), emit `z2 = 0` if T6's L3
  decision lands on branch (a) (the lean). Whichever branch L3 picks, the
  mapping is a single `select` here — do not bury disagreement handling in the
  builder (T1) or the gate (T4); it lives in this seam, the one place the two
  witness faces meet.
- **SSOT seam discipline (map §State / memory `project_bdg_evaluator_seam_ssot`).**
  Keep this a *two-call-site* seam (or update the in-source header at
  `bdg_observables.f90:6-13` that tallies call sites — U10's post-T6 cleanup
  already corrected the seam-call-site count: `eval_bdg_point` = 4 sites,
  `eval_bdg_pfaffian_witness_csr` = 1 site post-#42). Do not scatter
  `kitaev_majorana_number` calls across `main_topology.f90` call sites — the
  seam is the single face that consumes S1. Read
  `src/physics/AGENTS.md` before editing.
- **No vendoring.** `kitaev_majorana_number` and the Pfaffian reductions are
  in-tree (`src/math/pfaffian.f90`); use them, do not re-derive.

Output: the seam + a unit test driving a small synthetic stack through it. The
red-then-green **strict agreement** test (asserting S1×S2 agree at μ≈0.6601
across B) is T5 — T2 owns only the *wiring*, not the physics gate. Approach is
TDD (`superpowers:test-driven-development`).
