# T2 — S1×S2 agreement seam `eval_bdg_pfaffian_witness_product_csr`

> Child of [.scratch/bdg-u13-bloch-pfaffian/map.md](../map.md). Build branch;
> unblocks T3. Materialized from the map's T2 entry.

Type: task
Status: completed
Owner: this-session (2026-08-09 — proceeded per `/loop` continuation after handoff; landed seam on `feat/bdg-u13-bloch-pfaffian`)
Blocked by: 01 (closed 2026-08-09 — T1)

## Resolution (2026-08-09)

New public function `eval_bdg_pfaffian_witness_product_csr` added to
`src/physics/bdg_observables.f90`. Sibling to `eval_bdg_pfaffian_witness_csr`
and `eval_bdg_kitaev_majorana`. Takes T1's dense Bloch stack `H_k_array`,
computes S1 via `kitaev_majorana_number` (already on the seam as a thin
wrapper), extracts slice 1 → CSR scratch → `eval_bdg_pfaffian_witness_csr`
for S2, then maps through L3.

L3 mapping extracted to private `map_s1_s2_to_z2` helper — single
source of truth for the disagreement logic, not duplicated across call
sites. Disagreement_reason integer (0–3) is the writer-facing API for T3.

Verification (RED → GREEN):
- RED: build error (`Symbol 'eval_bdg_pfaffian_witness_product_csr' not
  found`) on main@52743b3; reverted seam, confirmed red, restored.
- GREEN: 53/53 unit (52 baseline + 1 new T2 target with 3 `@test`s for
  case (i) agreement S1=S2=-1 → z2=-1, case (i) S1=-1/S2=+1 → z2=-1,
  case (iv) S1=±1/S2=0 → z2=0 reason=2).
- 5/5 wire-BdG regression green.
- 1/1 lecture-13 acceptance gate green.

Trade-off documented: the dense→CSR conversion of slice 1 is private to
the seam — no public dense→CSR helper exists in `sparse_matrices`. The
conversion enumerates nonzeros and calls `csr_build_from_coo`. Future
refactor (parent plan ticket) should extract this to a public helper to
avoid duplication when other seams need it.

Out of scope (deferred): cases (iii) strict-sign-split and (v)
S1-closure-S2-misses — both require degenerate polar-decomposition
fixtures that the slim seam cannot read cleanly on a 16×16 stack (the
underlying `kitaev_majorana_number` hardcodes 2×2 det for n_odd>2 at
pfaffian.f90:154). T2 pins the determinable branches; T5's strict
agreement test owns the full 5-branch coverage using 4×4 (n_odd=1)
Kitaev fixtures with controlled det.

T2 unblocks T5 (TDD-red strict agreement test) and T3 (config knob).

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
