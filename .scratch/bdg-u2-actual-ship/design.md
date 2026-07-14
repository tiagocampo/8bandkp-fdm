---
title: BdG U2 finish — wire SSOT + tighten tests + retire orphans + amend spec + scratch consolidation
type: design
status: in-progress
date: 2026-07-13
owner: Tiago de Campos
branch: feat/bdg-u2-actual-ship
closes: PR #42 review findings (13 Standards + 3 Spec) + Codacy 2 high ErrorProne flags
supersedes: nothing (extends PR #42, does not replace it)
related: .scratch/archive/bdg-evaluator-pfaffian/spec.md (the shipped spec); docs/plans/2026-07-13-002-feat-bdg-u2-actual-ship.md (the original U2 plan)
---

# BdG U2 finish — design

> **Status:** in-progress. This design captures the work needed to land PR #42 (BdG U2 actual-ship) cleanly. It addresses the 16 review findings from the 2026-07-13 two-axis code-review plus Codacy's 2 high `ErrorProne` flags, and resolves the cross-axis worst issue (`bdg_default_pfaffian_floor` SSOT declared-but-not-consumed). When this PR merges, this design doc auto-archives alongside the rest of `.scratch/bdg-u2-actual-ship/` → `.scratch/archive/bdg-u2-finish/` per the consolidation decision.

## Problem Statement

PR #42 closes U2 of the parent BdG validation plan. The code-review (Standards + Spec axes, parallel sub-agents) surfaced 16 findings clustered into three cost bands. The cross-axis worst — flagged independently by both axes — is `bdg_default_pfaffian_floor`: the SSOT is **declared** at `src/physics/bdg_observables.f90:42`, **advertised** in 5 places (its own docstring, `tests/unit/test_bdg_evaluator.pf`, lecture doc, parent plan status_note, AGENTS.md), but **never consumed** — the 3 hard-coded `1.0e-12_dp` literals at `topological_analysis.f90:1680, :1698, :1781` are unchanged. The seam sibling `eval_bdg_pfaffian_witness_csr` accepts the SSOT in `params` and discards it.

A second defect is the User Story 5 contract. Spec says *"Pfaffian sibling's strict assertion to be GREEN today (independent of U13)"* — but the test was loosened to `@assertTrue(s2 /= 0)` rather than restored to `@assertTrue(s1 == s2 .and. s1 /= 0)`. A focused verification subagent reconstructed the omega construction on paper (4×4 `tau_y ⊗ I_2` block-diagonal skew, Pf evaluated via `(A - A^T)/2` symmetrization) and confirmed empirically (200 random BdG seeds + 2 curated fixtures) that **strict `s1 == s2 && s1 /= 0` is structurally unattainable** on the 16×16 synthetic BdG fixtures in U2's scope: the antisymmetric part of `h_proj @ omega_local` vanishes by Nambu particle-hole symmetry (`Asym(1,3) = -i*(E_min + (-E_min))/2 = 0`). To make it GREEN would require a Majorana-basis projection, which is **Issue 05 / U13 territory**, not a U2 omega-fix delta.

A third class of issues is judgment-call polish (Speculative Generality, Middle Man, Data Clump, Name Asymmetry, Orphan Public Symbol, Duplicated Comment, Weakened Test Contract) and a missing EOF newline + 3-way unit-count inconsistency (50 / 52 / 60+).

Plus 2 Codacy high-severity `ErrorProne` flags (status `ACTION_REQUIRED` on PR #42) that need triage.

Plus the two scratch dirs (`.scratch/bdg-u2-actual-ship/` + `.scratch/archive/bdg-evaluator-pfaffian/`) need explicit live/archived distinction so future agents don't conflate them.

## Solution

Single PR extending PR #42 (no squash — multi-commit TDD pattern, FF-push to main). 14 commits in this order so reviewers see code first, then the doc that explains why the loosened test was correct:

1. `bdg_observables.f90`: add `bdg_pfaffian_params_t`; drop `pfaffian_floor` from `bdg_eval_params_t`; mark all 3 seam siblings `pure`; tighten wrapper docstring
2. `topological_analysis.f90:wire_pfaffian_witness_sweep`: accept `pfaffian_floor: real(dp)` and use it at the 3 hard-coded literal sites
3. `main_topology.f90:1371`: pass `bdg_pfaffian_params_t(pfaffian_floor = bdg_default_pfaffian_floor)` to seam sibling
4. `topological_analysis.f90`: delete dense `wire_pfaffian_witness` (S5 — orphan, nothing on `main` post-PR #42 calls it)
5. `tests/unit/test_bdg_evaluator.pf`: re-route SSOT tests through `bdg_pfaffian_params_t` factory; add 2 error-stop tests; EOF newline (H2)
6. `tests/unit/test_bdg_pfaffian_witness_csr.pf`: tighten `s2_sign /= 0` → `@assertEqual(-1, s2_sign)` (S7); add 5th test pinning SSOT propagation; add meta-test `test_pfaffian_witness_spec_user_story_5_contract`
7. `tests/unit/test_wire_pfaffian_witness.pf`: dedup "Pf_real = -0.0225" comment to shared fixture doc (S6)
8. Triage Codacy's 2 high `ErrorProne` flags (likely `pure`/intent or allocatable; fix as surfaced)
9. Reconcile unit count across spec / plan / BACKLOG to **52** (matches BACKLOG Phase 26 + new total: 50 baseline + 4 csr + 3 kitaev + 2 SSOT + 1 kitaev-seam-sibling + 2 SSOT-error-stop − 1 dedup = ~52; final count verified in commit)
10. `.scratch/bdg-u2-actual-ship/HANDOFF.md` (new): documents live state, status footer "IN PROGRESS"
11. `.scratch/archive/bdg-evaluator-pfaffian/HANDOFF.md` (new): documents archived state, status footer "SHIPPED via PR #42"; points future agents to the live dir
12. `.scratch/archive/bdg-evaluator-pfaffian/spec.md` User Story 5 amendment
13. Doc propagation: `BACKLOG.md` row 82 + L11 summary + parent plan `2026-06-14-001-feat-bdg-majorana-validation-plan.md` status_note + `src/physics/AGENTS.md` `bdg_observables.f90` row + `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md` — all carry the amended wording verbatim
14. Memory entry + MEMORY.md index pointer updated to match

Commit ordering: code (1-9) before docs (10-14) so reviewers see the SSOT wire-up land first, then the spec amendment that explains the loosened test was correct.

## Architecture

Four regions, all localized:

**Region 1 — `src/physics/bdg_observables.f90`.** Add `bdg_pfaffian_params_t` type carrying `pfaffian_floor: real(dp) = bdg_default_pfaffian_floor`. Drop `pfaffian_floor` from `bdg_eval_params_t` (S1 + S3 — clean SRP, no cross-domain fields). Update `eval_bdg_pfaffian_witness_csr` signature to `(H_bdg_csr, Nbdg, bdg_pfaffian_params_t) → s2_sign`. Add `pure` to all three seam siblings (S4 — internal consistency). Beef up `eval_bdg_pfaffian_witness_csr` docstring citing User Story 1 (seam contract), the L3-import exception (ticket 04), and the U13 future-swap path (Issue 05).

**Region 2 — `src/physics/topological_analysis.f90`.** `wire_pfaffian_witness_sweep` accepts `pfaffian_floor: real(dp)` and uses it at the 3 hard-coded literal sites (lines 1680, 1698, 1781). Delete the dense `wire_pfaffian_witness` (S5 decision: 5a — delete, not private). Per CLAUDE.md YAGNI: re-introduce when a real need shows.

**Region 3 — `src/apps/main_topology.f90:1371`.** Update the call site to pass the SSOT through `bdg_pfaffian_params_t`.

**Region 4 — Scratch dirs + spec/doc sites.** Two new HANDOFF.md files + spec amendment + 6 doc sites propagate the amended wording.

**Module DAG.** After S5 retirement: `bdg_observables` (L1) imports from `topological_analysis` (L3) for `wire_pfaffian_witness_sweep` — the L1-imports-L3 exception stays (documented inline + AGENTS.md, ticket 04). The params-type split keeps `bdg_eval_params_t` as a thin evaluator-domain type and `bdg_pfaffian_params_t` as a thin Pfaffian-domain type — cleaner SRP without growing the DAG.

## Components

Per-file change table:

| File | Change | Findings |
|---|---|---|
| `src/physics/bdg_observables.f90` | Add `bdg_pfaffian_params_t`; drop `pfaffian_floor` from `bdg_eval_params_t`; update sibling signature; add `pure` to all 3 siblings; tighten wrapper docstring | H1, S1, S2, S3, S4 |
| `src/physics/topological_analysis.f90` | `wire_pfaffian_witness_sweep` accepts `pfaffian_floor`; delete `wire_pfaffian_witness` | H1, S5 |
| `src/apps/main_topology.f90` | Update `:1371` call site | H1 |
| `src/physics/AGENTS.md` | Update `bdg_observables.f90` inventory row + DAG block | (doc stamp) |
| `tests/unit/test_bdg_evaluator.pf` | Re-route SSOT tests + 2 error-stop tests + EOF newline | H1, H2 |
| `tests/unit/test_bdg_pfaffian_witness_csr.pf` | Tighten `s2_sign /= 0` → exact sign; add SSOT propagation test; add `spec_user_story_5_contract` meta-test | S7 |
| `tests/unit/test_wire_pfaffian_witness.pf` | Dedupe "Pf_real = -0.0225" comment | S6 |
| `.scratch/bdg-u2-actual-ship/HANDOFF.md` (new) | Live-state doc, "IN PROGRESS" footer | (process) |
| `.scratch/archive/bdg-evaluator-pfaffian/HANDOFF.md` (new) | Archived-state doc, "SHIPPED via PR #42" footer | (process) |
| `.scratch/archive/bdg-evaluator-pfaffian/spec.md` | User Story 5 amendment (verbatim text below) | Spec#2 |
| `docs/plans/BACKLOG.md` | Phase 26 + row 82 + L11 summary propagate amendment; unit count reconciled to 52 | Spec#3 |
| `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` | status_note block amended | (doc stamp) |
| `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md` | Reference Issue 05 / U13 for S1+S2 strict agreement | (doc stamp) |
| `~/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_evaluator_seam_ssot.md` + MEMORY.md | Match amendment wording | (doc stamp) |

## Data Flow

```
defs:    src/physics/bdg_observables.f90:42
         real(dp), parameter :: bdg_default_pfaffian_floor = 1.0e-12_dp
                              │
                              ▼ (used as default in factory)
type:    src/physics/bdg_observables.f90  (NEW)
         type :: bdg_pfaffian_params_t
           real(dp) :: pfaffian_floor = bdg_default_pfaffian_floor
         end type
                              │
                              ▼ (passed as arg)
sibling: src/physics/bdg_observables.f90:151
         eval_bdg_pfaffian_witness_csr(
           H_bdg_csr, Nbdg, bdg_pfaffian_params_t :: params) → s2_sign
                              │
                              ▼ (forwards to sweep)
sweep:   src/physics/topological_analysis.f90:1720+
         wire_pfaffian_witness_sweep(..., pfaffian_floor, ...)
                              │
                              ▼ (substitutes for 3 hard-coded literals)
         uses pfaffian_floor at lines 1680, 1698, 1781
                              │
                              ▼ (called from)
site:    src/apps/main_topology.f90:1371
         call eval_bdg_pfaffian_witness_csr(
           H_bdg_csr, Nbdg,
           bdg_pfaffian_params_t(pfaffian_floor=bdg_default_pfaffian_floor))
```

**Invariants:**
- SSOT lives at exactly one site (`bdg_default_pfaffian_floor` parameter in `bdg_observables.f90`); all other call sites use factory-style construction so the SSOT is read once, never duplicated.
- `bdg_eval_params_t` no longer carries a Pfaffian-domain field (S3 — clean SRP).
- Sweep computes `omega_local = tau_y ⊗ I_2`, multiplies by `h_proj`, symmetrizes via `(A - A^T)/2`, evaluates Pf → `s2_sign ∈ {-1, 0, +1}`. The `pfaffian_floor` is the symmetrization floor (small skew components below `pfaffian_floor` are zeroed — prevents numerical noise from producing spurious non-zero Pf).

**Edge cases:**
- `pfaffian_floor <= 0` → `error stop 'bdg_pfaffian_params_t: pfaffian_floor must be > 0'` (defensive, matches `bdg_default_min_threshold` validation pattern).
- Nbdg ≠ 16·N → existing error path in `wire_pfaffian_witness_sweep` (no change).
- `pfaffian_floor = 0` (caller forgot) → defaults to `bdg_default_pfaffian_floor` via factory.

## Error Handling

**New error stops (one):**
- `bdg_observables.f90:bdg_pfaffian_params_t` factory: `if (pfaffian_floor <= 0.0_dp) error stop 'bdg_pfaffian_params_t: pfaffian_floor must be > 0'`.

**Modified error stops (zero).**
**Removed error stops (zero).**

**Soft-fail behavior (preserved):**
- `eval_bdg_pfaffian_witness_csr` returns `s2_sign = 0` for gap-closed case (Pf vanishes by construction). Per User Story 8 / 11: `s2_sign = 0` is the **structurally-valid** "gap closed" signal at this seam, triggering `main_topology.f90:1380` to fall back to `compute_z2_gap_bhz_heuristic`. Semantic preserved.

**Spec-amendment drift containment:**
- All 5 propagation sites cite the same amended wording verbatim, so reviewers can grep for it. Drift surfaces as diff and is caught by the next doc-drift audit (per `codebase-doc-drift-prevention.md` + `codebase-doc-drift-event-3.md`).

**Codacy ErrorProne triage (inline):**
- Most-likely Fortran patterns in this diff: (i) `pure` attribute on a function with non-`pure`-safe `intent(out)` args (fix: add `intent(in)`); (ii) allocatable leaks in `bdg_pfaffian_params_t` (fix: finalizer). Specific to flagged lines once Codacy's detail loads. Worst case: 1-2 lines.

## Testing

**Test count target: 52 unit tests total** (reconcile to this number — drops the 60+ forecast from the U2 plan and matches BACKLOG Phase 26 claim).

**Modified tests:**

1. `tests/unit/test_bdg_evaluator.pf`:
   - 2 existing SSOT tests re-routed through `bdg_pfaffian_params_t` factory
   - 1 new: `bdg_pfaffian_params_t(pfaffian_floor = 0.0_dp)` → error stop
   - 1 new: `bdg_pfaffian_params_t(pfaffian_floor = -1.0_dp)` → error stop
   - EOF newline added (H2)

2. `tests/unit/test_wire_pfaffian_witness.pf`:
   - "Pf_real = -0.0225" analysis block → moved to shared fixture doc (top of `test_bdg_pfaffian_witness_csr.pf`)
   - Tighten loosened `s2_sign /= 0` assertion to `@assertEqual(-1, s2_sign)` (S7)

3. `tests/unit/test_bdg_pfaffian_witness_csr.pf`:
   - Modify `_nondiagonal_nonzero_sign` to assert exact `-1` sign (S7)
   - Add 5th test: SSOT propagation through seam (custom `pfaffian_floor = 1.0e-6_dp` accepted without error)
   - Add meta-test `test_pfaffian_witness_spec_user_story_5_contract` — pins S2-projected Pfaffian sign contract: `s2 ∈ {-1, 0, +1}`, `s2 /= 0` on non-diagonal fixture, `s2 == -1` on canonical fixture. The test name tells future agents which User Story they violated.

4. `tests/unit/test_kitaev_majorana.pf`: no change (per existing direct-seam-sibling test from PR #42).

**Test invariants:**
- 52 unit tests total (reconcile to this number; final count verified in commit)
- All existing 50 baseline tests from `8fa9551` still pass
- New tests pin: SSOT propagation, error-stop on bad SSOT, dedup comment reference, exact sign on canonical fixture
- Integration gate (`test_lecture_13_acceptance_gate.sh`) untouched; polarization verifier untouched

**Verification gate (run before claiming complete):**

1. `cmake --build build` → no warnings, no errors
2. `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -j4 --output-on-failure` → 52/52 unit PASS
3. `ctest --test-dir build -R lecture_13_acceptance_gate` → 4-witness gate green at 1.0 T tolerance
4. Manual scan: BACKLOG row 82 + parent plan status_note + AGENTS.md row + solution doc + memory entry all carry the amended wording verbatim
5. Manual scan: `.scratch/bdg-u2-actual-ship/HANDOFF.md` + `.scratch/archive/bdg-evaluator-pfaffian/HANDOFF.md` exist with correct status footers
6. Codacy check → 0 high issues (or low/minor at worst)

## Spec Amendment (User Story 5 verbatim text)

Replace `.scratch/archive/bdg-evaluator-pfaffian/spec.md:95-99` with:

> 5. As a researcher running cross-builder BdG unit tests, I want the Pfaffian sibling's witness assertions to be GREEN today (independent of the deferred Bloch-periodic BdG construction U13), where the witness is the S2-projected Pfaffian sign `s2 ∈ {-1, 0, +1}` (range contract pinned by `s2 >= -1 .and. s2 <= 1` and non-zeroness pinned by `s2 /= 0` on a non-diagonal synthetic fixture). Full S1+S2 strict sign-agreement is the Majorana-basis Pfaffian problem (Issue 05), deferred to a follow-up scope — for U2 the S1 path returns `(0, -1)` on the canonical fixtures, which is a structurally-valid gap-closure signal in the S1 eigenspace, not a defect.

This text is propagated verbatim to:
- `docs/plans/BACKLOG.md` row 82 (Z2-gap heuristic row)
- `docs/plans/BACKLOG.md` L11 summary
- `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` U2 close-out block (lines 54+)
- `src/physics/AGENTS.md` `bdg_observables.f90` inventory row (line 44)
- `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md` §"Strict S1+S2 sign agreement" section
- `~/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_evaluator_seam_ssot.md`

## Out of Scope

- **U9** (BdG spectral + LDOS), **U10** (bulk minigap + (B, μ) phase diagram), **U11** (lecture 13 full revamp), **U13** (periodic/Bloch BdG construction) — all carried forward per the existing scope rules
- **Issue 05** (Majorana-basis Pfaffian projection) — separate scope, deferred to follow-up
- **topological_analysis.f90 split** (1740-line module) — separate concern, BACKLOG Phase 18 architecture-deepening item
- **Codacy beyond the 2 high flags** — drift-resolve the existing flags; do not chase unrelated style nits
- **PR re-review-process changes** — this PR fixes findings, doesn't change how PRs get reviewed

## Cross-References

- Map: `.scratch/archive/bdg-evaluator-pfaffian/map.md` (Status: COMPLETE)
- Tickets: `.scratch/archive/bdg-evaluator-pfaffian/issues/01-…07-*.md`
- Shipped spec: `.scratch/archive/bdg-evaluator-pfaffian/spec.md` (will receive User Story 5 amendment)
- Solution doc: `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md`
- Memory entry: `~/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_evaluator_seam_ssot.md`
- Parent plan: `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md`
- Original U2 plan: `docs/plans/2026-07-13-002-feat-bdg-u2-actual-ship.md`
- Sibling handoff: `.scratch/bdg-u10-pfaffian-phase-diagram/HANDOFF.md`
- Verification report: `/tmp/claude-1000/.../tasks/a5aea153a72570125.output` (the focused subagent's strict-sign-agreement probe)