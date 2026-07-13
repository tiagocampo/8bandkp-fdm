# PR #41 — Blocker-Fix Addendum Design Spec

**Status**: IMPLEMENTED (PR #41 shipped 2026-07-06, all addendum decisions executed via Phase A+B+C + C.8 = 17 commits on `feat/bdg-p1-stabilization`)
**Date**: 2026-07-05
**Branch**: same as PR #41 head (`feat/bdg-validation-pass2`)
**Source PR**: #41 (BdG/Majorana ground-up validation pass, 26 commits, +11k/-22.7k)
**Source review**: 6-dimension adversarial review (plan-adherence, physics, test-coverage, architecture, documentation, anti-ponytail), 4 P0 / ~10 P1 / ~10 P2 / ~10 P3 findings
**ce-doc-review round 1**: 4 reviewers (coherence, feasibility, adversarial, scope-guardian); 6 P0 / ~12 P1 / ~6 P2 / ~3 P3 findings on the spec itself; auto-resolved per user direction.
**Companion docs**:
- `docs/superpowers/specs/2026-07-01-bdg-p1-fix-design.md` — the in-flight P1 spec (Phase 0–4); receives Section 6.1 line-23 correction
- `docs/superpowers/plans/2026-07-01-bdg-p1-fix.md` — implementation plan
- `docs/adr/0008-bdg-p1-invariants.md` — receives Section 6.3 §3 amendment

**This spec**: an ADDENDUM to the in-flight P1 spec. It addresses gaps the 6-dimension adversarial review surfaced that the existing P1 spec did not cover. Sections already handled by the existing spec are referenced but not duplicated.

## 1. Locked decisions (from brainstorming 2026-07-05, refined by ce-doc-review round 1)

| # | Decision | Choice | Rationale |
|---|---|---|---|
| D1 | Q1 (Peierls re-decision) | **KEEP symmetric `add_peierls_coo(-B)` on hole block** | Source `bdg_hamiltonian.f90:406-419` keeps the call with in-code rationale "the symmetric add_peierls_coo(-B) IS required to make test_bdg_phs pass; the double-count claim in ADR 0008 §3 was wrong". This is the AUTHORITATIVE source — it has the empirical PHS oracle verification. ADR 0008 §3 must be amended (Section 6.3). |
| D2 | PR base | **PR #41 head (`feat/bdg-validation-pass2`)** | User directive ("on the same open pr"); the local `feat/bdg-p1-stabilization` branch work rolls up into the PR via rebase/cherry-pick |
| D3 | Scope | **A + B blockers + D hygiene** (Chunks A, B, D from the decomposition in chat) | Per `/ponytail ultra`: minimum needed to merge PR #41. Chunk E (anti-ponytail trim) explicitly deferred. |
| D4 | Cross-builder identity test (PRD R4) | **Included in scope, with shared-H₀ invariant** | PRD R4 + ADR 0007 + existing P1 spec all defer the same item — that's a signal to write it. Per ce-doc-review scope-guardian: must enforce shared-H₀ (not just "synthetic 8×8") to avoid re-encountering the deleted-test failure mode per `tests/unit/test_bdg_phs.pf:678-694`. |
| D5 | Acceptance gate witness count | **3-witness, not 4-witness** | Per ce-doc-review adversarial: `output/wire_slim_pfaffian_witness.dat` has no Fortran producer (`lecture_13_topological.py:212` reads it; no source has ever written it; git history confirms). `bcrit_pfaffian = bcrit_curve` at line 224 always fires. Until U13 lands the slim-Pfaffian sweep writer, the 4th witness is unreachable — present as 3-witness honestly, not as 4-witness theater. |
| D6 | Existing P1 spec corrections | **Two line-items corrected inline in this addendum** | D6.1: spec-line-23 doc-drift (Peierls) per existing-spec §1 line 23. D6.2: ADR 0008 §3 reversal (Section 6.3). |
| D7 | Synthetic fallbacks (env-gated) | **DROP; fail-loud only** | Per ce-doc-review (3 reviewers: feasibility, scope, adversarial): env-gated synthetic re-introduces the false-PASS mode C-2/C-3 purports to fix. Real path is the default; real-path errors → exit non-zero. No silent fallback. |
| D8 | Config fixture canonical name | **`tests/regression/configs/wire_inas_gaas_bdg_topological.toml`** (no `_2d` suffix) | Per ce-doc-review feasibility: matches what `test_lecture_13_acceptance_gate.sh:128` and the dense-QW mirror's comment reference. The `_2d.toml` was the addendum's error. |

## 2. Approach: TDD-doubled on P0 path, single-commit for P1/P2

For the 3 P0 blockers (Section 3): failing test → fix → regen pattern. The acceptance gate is the regression net. For P1/P2 (Sections 4–5): single-commit refactors since the value is the doc/test surface, not the function under test.

**Layer 1 — P0 blockers (single merge-blocker path)**:
- TDD-red for slim-Pfaffian row bug (test asserts S2 projects onto conduction bands)
- TDD-green: correct row indices + integrate with the production `wire_pfaffian_witness_sweep` path
- TDD-red for `verify_majorana_polarization.py` (real eigensolve output, not synthetic)
- TDD-green: rebuild verifier against real output + add wire-polarization emitter to `main_topology.f90`

**Layer 2 — P1 + Doc fixes**:
- Cross-builder identity test (R4 from PRD) with shared-H₀ invariant
- Vacuous M=±1 test rebuild (Lutchyn-Oreg per ADR 0008 §1, covered by existing P1 spec Phase 1.4 — verification only here)
- Vacuous S1/S2 agreement test rebuild (real witness output)
- 9 issue file Status footers (per memory `codebase-doc-drift-prevention.md`)
- Lecture disclosures (cross-builder deferral, Pfaffian-iσ_y vs real ±1, slim Pfaffian plot caption)
- ADR 0008 §3 amendment (Section 6.3)

## 3. P0 Blockers (3 fixes — merge-blockers)

### 3.1 Slim Pfaffian row-index bug (C-1)

**Bug**: `src/physics/topological_analysis.f90:1722` and `:1827` use `idx = [7, 8, n_sp+7, n_sp+8]` which assumes row `r` corresponds to band `r`. Wire H₀ uses band-major layout (`hamiltonian_wire.f90:1208-1209`: `row = (band-1)*N + site`). For any wire with N>1, this reads rows 7,8 (within band 1, HH↑), not bands 7-8 (CB).

**Impact**: `main_topology::eval_wire_bdg_gap` (lines 1433-1436) uses this witness to set the Z2 invariant. Z2 colormap is unphysical for any real wire.

**Fix**: Replace `idx` with band-major equivalents across multiple sites (multi-site scan, NOT single-site per ce-doc-review adversarial). For each site `s = 1..N`:
```
idx_s = [6*N + s, 7*N + s, n_sp + 6*N + s, n_sp + 7*N + s]
```
The implementation scans all sites, computes Pf(s) for each 4×4 subblock, and returns the first non-zero Pf (or the site with maximum |Pf|, with a documented tiebreak rule). A single-site projection at `s = ⌊N/2⌋` would lose spatial information AND yield Pf≈0 if the center site sits at a barrier.

**Test (TDD-red)**:
- Build a real 3×3 wire (N=3 → 24×24 BdG) with μ in the topological gap
- Call `wire_pfaffian_witness_sweep` at B=0 and at B=B_crit
- Assert S2 projects onto conduction-band rows (not band-1 rows); assert the returned Pf is non-zero
- Assert the scan-site selector picks a non-trivial site (not just first-pass)

**Files**:
- `src/physics/topological_analysis.f90:1722, :1827` (fix — multi-site scan)
- `tests/integration/test_slim_pfaffian_witness_projection.py` (new, ~80 LOC including the 3×3 wire fixture builder)
- `tests/integration/test_slim_pfaffian_witness_projection.sh` (ctest wrapper, ~10 LOC)

### 3.2 Verifier theater (C-2) — FIXED with file producer

**Bug**: `tests/integration/verify_majorana_polarization.py:28-92` generates synthetic numpy P_M profiles, prints PASS unconditionally, never calls `topologicalAnalysis`. Lecture §13.7.0 embeds the plot without disclosure.

**Per ce-doc-review feasibility**: `output/majorana_polarization.dat` does NOT exist — actual emitter is `output/majorana_profile.dat` per `outputFunctions.f90:602`. The wire-polarization emitter does not exist in `main_topology.f90` (only QW at line 587/870 emits `majorana_profile`).

**Fix**:
1. **Add wire-polarization emitter to `main_topology.f90`** in the wire branch (analogous to dense-QW at line 870): `call majorana_polarization(...)` and `call write_majorana_polarization(...)` (new helper in `outputFunctions.f90`). New file emitted: `output/majorana_polarization.dat` with structured columns `(index, P_M, tau_z, half_wire_integral)`.
2. **Rewrite `tests/integration/verify_majorana_polarization.py`**: (a) call `topologicalAnalysis` on `tests/regression/configs/wire_inas_gaas_bdg_topological.toml` (canonical per D8), (b) parse `output/majorana_polarization.dat`, (c) assert `half_wire_mzm > 4 × half_wire_accidental` per `test_majorana_polarization.pf:317` semantics, (d) save real plot to `docs/lecture/figures/lecture_13_majorana_polarization_wire.png` with explicit caption "(real eigensolve output, Issue 04/U7)" in lecture markdown.
3. **NO env-gated synthetic fallback** (per D7). If real path fails, verifier exits non-zero with explicit error message naming the missing output file. The lecture figure generation becomes a fail-loud prerequisite.

**Files**:
- `src/apps/main_topology.f90` (add wire-polarization emitter — Phase A.3a)
- `src/io/outputFunctions.f90` (add `write_majorana_polarization` helper — Phase A.3b)
- `tests/integration/verify_majorana_polarization.py` (rewrite, ~150 LOC)
- `docs/lecture/13-topological-superconductivity.md` (figure caption update)

### 3.3 S1/S2 synthetic plot (C-3) — DEFERRED to U13

**Bug**: `tests/integration/verify_wire_bdg_topological.py:159-160` sets `s1_sign = s2_sign = np.where(B < B_crit, +1, -1)` from the 1D curve. Agreement is by construction; no real Pfaffian data feeds the plot.

**Per ce-doc-review adversarial**: `output/wire_slim_pfaffian_s1.dat` and `s2.dat` are referenced by the spec but no Fortran source emits either. S1 requires full diagonalization (per `topological_analysis.f90:1702-1703`: "S1 needs full diagonalization, deferred to U13"). The `decc6ba refactor(test): rename verify_majorana_polarization.py -> plot_*.py` commit's rationale is correct (matches the docstring role); the addendum's earlier "restore the original verifier role" claim was wrong.

**Fix**: This is **deferred to U13** (full wire Pfaffian sweep). For this PR:
1. **Add an explicit inline comment** at `verify_wire_bdg_topological.py:159-160` declaring the S1/S2 plot as "synthetic illustrative, real witness gated by U13".
2. **Update lecture §13.7.4 caption** to "(synthetic illustrative — real Pfaffian witness output requires U13)" per Section 5.2.
3. **Update lecture_13_topological.py:212** to gracefully skip the wire_pfaffian witness (line 224 fallback already exists — keep it, label it "(deferred to U13)").
4. **Reconcile D5 with §3.3**: the acceptance gate is 3-witness (wire_curve, wire_2d, qw_dense); the 4th witness (wire_pfaffian) is reserved for U13. Spec §9 acceptance criteria reflect this.

**Files**:
- `tests/integration/verify_wire_bdg_topological.py:159-185` (comment + caption only; no logic change)
- `scripts/lecture_13_topological.py:224` (re-label fallback to "deferred to U13")
- `docs/lecture/13-topological-superconductivity.md` §13.7.4 (caption update)

### 3.4 Acceptance gate tightening (C-4 + H-1) — VERIFY, don't add

**Per ce-doc-review feasibility**: The absolute window guard `0.5 ≤ BCRIT_MIN, BCRIT_MAX ≤ 6.0` is **already present** at `tests/integration/test_lecture_13_acceptance_gate.sh:102-107`. The 2.0 T tolerance at lines 87-99 is also in place. Section 3.4 work is **verification only**, no new code.

**Fix**:
1. **Verify commit `19c1da6 test: add absolute window guard + regex anchor to lecture 13 acceptance gate`** is reachable from PR #41 head.
2. **Verify commit `3ead499 test: add 'acceptance-gate' ctest label per spec §7.5`** is reachable.
3. **Tighten tolerance to 1.0 T for the 3-witness configuration** (per ce-doc-review adversarial P0): with wire_pfaffian now explicitly deferred (D5), the actual distinguishable range is ~0.5 T, so 1.0 T is defensible. Update `test_lecture_13_acceptance_gate.sh:87` from `TOL_BCRIT_RANGE="2.0"` to `TOL_BCRIT_RANGE="1.0"`.
4. **Apply D5**: update `lecture_13_topological.py` reconciliation table from 4-witness to 3-witness label, with `wire_pfaffian` row marked "(deferred to U13)".

**Files**:
- `tests/integration/test_lecture_13_acceptance_gate.sh:87` (tolerance: 2.0 → 1.0)
- `scripts/lecture_13_topological.py:212-225` (3-witness label)

**Net effect of 3.1–3.4**: The slim Pfaffian bug is fixed at the source (multi-site scan); the wire-polarization verifier runs against real output; the S1/S2 plot is honestly labeled as illustrative pending U13; the acceptance gate becomes 3-witness and tightens from 2.0 T to 1.0 T.

## 4. P1 Test Soundness

### 4.1 Cross-builder identity test (PRD R4) — with shared-H₀ invariant

**Bug**: Issue 03 AC #2 / PRD R4 ("both builders produce byte-identical hole blocks for the same H₀") was deferred by ADR 0007 line 120 and the existing P1 spec also defers it. The original test was deleted (per `tests/unit/test_bdg_phs.pf:678-694`) because it compared structurally-different H₀s (wire builder's `ZB8bandGeneralized` vs QW builder's `ZB8bandQW`).

**Per ce-doc-review feasibility + scope**: The original spec recipe ("Build a synthetic 8×8 H₀ at k=0, pass through both builders with N=1") is **infeasible** — the QW builder's `validate()` rejects `fd_step < 2`. Also, the test would not exercise the canonical form at generic k where Peierls matters.

**Fix**: Construct a 2-site wire fixture (N=2, fd_step=2 passes validate) and a QW fixture with the SAME upstream H₀, then compare hole-block entries after the canonical-form transform.

**Recipe** (per `test_bdg_phs.pf:691-694` follow-up option (a)):
1. Construct a single 8×8 H₀(k=0) using the canonical 8-band Hamiltonian at k=0 with no B-field, no Peierls (so k=0 + B=0 collapses to the regime where wire-form `-H₀ᵀ(+k)` and canonical `-conjg(H₀(-k))` are demonstrably equivalent).
2. Build a wire fixture where the 2D-confined H₀ reduces to that same 8×8 (N=2 spatial points; 2×8 = 16 single-particle, then BdG = 32×32; the hole block at site 1 must match the QW's hole block exactly).
3. Build a QW fixture with kpterms/profile chosen to reproduce that same 8×8 at fd_step=2 (which passes `validate()`).
4. Assert byte-identical hole blocks after the canonical-form transform.
5. **Also test at generic k=π/(2a), Bx≠0**: this is the regime where Peierls phase matters; the test passes only if both builders route through `build_bdg_hole_block` correctly.

**Test file**: `tests/unit/test_cross_builder_hole_block_identity.pf` (~80 LOC including 2-fixture builders)

### 4.2 Vacuous M=±1 tests (test_kitaev_majorana.pf) — verification only

**Bug**: `tests/unit/test_kitaev_majorana.pf:204-211` only asserts `M ∈ {-1, 0, +1}`. Per existing P1 spec line 13 + ADR 0008 §1, Kitaev Majorana number reformulated to Lutchyn-Oreg sign-of-det.

**Per ce-doc-review scope-guardian**: This is verification, not new work — the Lutchyn-Oreg implementation lives in existing P1 spec Phase 1.4.

**Fix**: Verify that PR #41 carries the Lutchyn-Oreg implementation (commits covering Phase 1.2–1.4). If absent, add a single Phase B.3 commit pair (TDD-red + green) covering the test rebuild. Otherwise, no commit.

### 4.3 Vacuous S1/S2 agreement test (test_wire_pfaffian_witness.pf)

**Bug**: `tests/unit/test_wire_pfaffian_witness.pf:88-94` AC test allows `if (s1==0 or s2==0) then both must be 0` escape hatch. For diagonal-in-(c,c†) BdG, both always return 0; the agreement branch is unreachable.

**Fix**: After C-1 (Section 3.1) lands, replace this test with one that asserts `s1 == s2 != 0` for a real fixture — no zero-escape hatch. Use the same 3×3 wire BdG fixture from Section 3.1.

**Style constraint (per memory `feedback_pfunit_macro_single_line.md`)**: `@assertTrue` / `@assertEqual` MUST be on a single source line (no `&` continuation). If the assertion expression is too long, factor intermediate booleans into local variables (e.g., `logical :: non_zero_agree = (s1 == s2) .and. (s1 /= 0)`) and assert the named boolean. Reference `test_majorana_polarization.pf:220` for the existing single-line pattern.

**Test file**: `tests/unit/test_wire_pfaffian_witness.pf:88-118` (rewrite, ~30 LOC)

## 5. Documentation fixes (minimal, ~30 lines of markdown added)

### 5.1 Status footers on 9 issue files + 15 archived PRDs (M-2 from doc review)

Per memory `codebase-doc-drift-prevention.md`, every PR touching behavior must update spec + plan status footers. Per ce-doc-review scope-guardian: use **uniform vocabulary** — pick `COMPLETE` per `.scratch/AGENTS.md` precedent (the existing PRD.md line 1 uses `**Status**: COMPLETE (2026-06-28)`).

**Fix**: Add `**Status**: COMPLETE (2026-07-05)` as the first line of each:

Issue files:
- `.scratch/archive/bdg-majorana-validation/issues/00-pure-bdg-evaluator.md`
- `.scratch/archive/bdg-majorana-validation/issues/01-kitaev-pfaffian-harness.md`
- `.scratch/archive/bdg-majorana-validation/issues/02-phs-oracle.md`
- `.scratch/archive/bdg-majorana-validation/issues/03-hole-block-unified.md`
- `.scratch/archive/bdg-majorana-validation/issues/04-majorana-polarization.md`
- `.scratch/archive/bdg-majorana-validation/issues/05-dense-qw-rung.md`
- `.scratch/archive/bdg-majorana-validation/issues/06-bdg-ldos-spectral.md`
- `.scratch/archive/bdg-majorana-validation/issues/07-wire-phase-diagram.md`
- `.scratch/archive/bdg-majorana-validation/issues/08-lecture-13-gate.md`

PRD files (15 total per `ls .scratch/archive/*/PRD.md | wc -l`):
- All `.scratch/archive/*/PRD.md` files lacking a Status footer.

This is 24 single-line edits; run as a single shell command with `for f in ...; do sed -i '1i **Status**: COMPLETE (2026-07-05)' "$f"; done` after a pre-check that the line isn't already present.

### 5.2 Lecture disclosures (H-2, H-11, M-14 combined)

Per the doc consistency review, three lecture claims either over-claim or hide deferred requirements. Per ce-doc-review feasibility P2: the pairing-matrix disclosure is physics-misleading.

**§13.5.1 cross-builder disclosure**:
- Current claim: "guarantees the canonical form on every call path"
- Updated claim: "the shared wrapper `build_bdg_hole_block` at `bdg_hamiltonian.f90:94-103` enforces the canonical form on every call path. Cross-builder identity (R4 from PRD) is asserted by code structure AND verified by `tests/unit/test_cross_builder_hole_block_identity.pf` (added in this PR)."

**§13.5.2 pairing-matrix iσ_y clarification** (per feasibility P2):
- Current claim: `Δ = δ₀(iσ_y ⊗ I_4)`
- Issue: `bdg_hamiltonian.f90:361` stores `cmplx(delta_0 * pairing_sign(row), 0.0_dp, kind=dp)` — purely real ±1, no imaginary factor
- Updated claim: "Pairing-sector sign: each band is +1 (spin-up) or -1 (spin-down). Kramers pairs (1,4), (2,3), (5,6), (7,8) take signs (+1,-1) respectively, matching the iσ_y⊗I₄ structure. The imaginary factor `i` is absorbed into the complex-conjugation convention; the codepath stores `cmplx(delta_0 * pairing_sign, 0.0_dp)` (purely real) per `bdg_hamiltonian.f90:361` because the upper-triangular block (electron→hole) and lower-triangular (hole→electron, adjoint) encode the imaginary factor across the diagonal (see lines 539-553)."

**§13.7.4 slim-Pfaffian plot caption**:
- After Section 3.3 lands, update the figure caption to "(real witness output via `wire_pfaffian_witness_sweep` is deferred to U13; the 3-witness acceptance gate (wire_curve, wire_2d, qw_dense) is the current regression net per spec §7.5)"

### 5.3 Empty golden data — ALREADY DONE by `7312e97`

**Per ce-doc-review scope-guardian**: Local commit `7312e97 test(reg): wire_bdg_topological_2d skip z2 column + non-flat assertion` already handles the z2-column replacement. The transitions.dat file (mentioned in existing P1 spec Phase 2.1 consequence table at line 374, NOT line 372 as the earlier spec text said per adversarial's P2 finding) is referenced separately.

**Fix**: Verify on PR #41 that commit `7312e97` lands and removes the z2 column. If `tests/regression/data/wire_inas_gaas_bdg_topological_2d_transitions.dat` is still referenced by any test config, grep and decide removal vs population. **No new commit in this addendum**.

## 6. Correcting the existing P1 spec + ADR (D6 doc-drift)

### 6.1 Line-23 doc-drift in `2026-07-01-bdg-p1-fix-design.md`

The existing spec's locked-decision table at line 23 says:
> "Symmetric Peierls philosophy | A — remove from wire, rely on `conjg()`"

But the implementation summary at line 13 says:
> "ADR 0008 §3 'double-count Peierls' claim was INVESTIGATED and disproved — symmetric `add_peierls_coo(-B)` is REQUIRED for class-D PHS, restored with documentation."

These are contradictory. The implementation summary wins (it cites PHS oracle verification; the decision table does not). Source code at `bdg_hamiltonian.f90:406-419` confirms.

**Fix**: Edit line 23 of `docs/superpowers/specs/2026-07-01-bdg-p1-fix-design.md` to read:
> "Symmetric Peierls philosophy | **A — keep symmetric (re-investigated)** — symmetric `add_peierls_coo(-B)` on hole block is REQUIRED for class-D PHS at generic k with Bx≠0. The `-conjg()` transform alone is insufficient. See ADR 0008 §3 amendment (Section 6.3 of `2026-07-05-pr41-completion-design.md`)."

### 6.2 ADR 0007 Layer D reference update

Per ce-doc-review adversarial P2: ADR 0007 §Decision Layer D describes canonical form "H_hole = -conjg(H0(-k))" but does not mention symmetric Peierls. Add a footnote citing the symmetric-Peierls call in `bdg_hamiltonian.f90:420-428`.

**Fix**: Add a footnote to ADR 0007 §Decision Layer D: "Peierls phase is applied symmetrically: electron block via `ZB8bandGeneralized` at `+kz` with `+B`, hole block via explicit `add_peierls_coo(-B)` with row filter 8N+1..16N. Canonical form transform `-conjg()` does not include Peierls. See `bdg_hamiltonian.f90:406-419` for the empirical investigation."

### 6.3 ADR 0008 §3 reversal (D6.2)

**Per ce-doc-review adversarial P0 (highest-priority fix)**: ADR 0008 §3 currently mandates "Remove this call" for the symmetric Peierls. Source code AND existing P1 spec line 13 both say keep. ADR is the load-bearing reference; this is a doc-drift event per memory `codebase-doc-drift-prevention.md`.

**Fix**: Edit `docs/adr/0008-bdg-p1-invariants.md` §3 to read:
> "**Reversed 2026-07-05**: Symmetric `add_peierls_coo(-B)` on hole block is KEPT (not removed). The `add_peierls_coo` function applies `exp(-iφ)` where `φ = e*Bx*(y_i-y_j)/hbar`. The `-conjg()` transform in `build_bdg_hole_block` does NOT include Peierls. Empirical PHS oracle verification (rel_resid 1.25e-1 → 0) on removing the call demonstrated the call is REQUIRED for class-D PHS at generic k with Bx≠0. The original 'double-count' rationale was disproved 2026-07-01; the investigation note lives at `bdg_hamiltonian.f90:406-419`."

### 6.4 Other doc-drift items from the 6-dimension review

These were not in my Section 4/5 because the existing P1 spec or local branch already addresses them:
- AGENTS.md `add_peierls_coo exp(-iφ)` note — commit `284aa93` ✓
- dense-QW block (2,1) structural+sign — commit `1984ca5` ✓
- Minigap-pattern regression replacing z2 golden — commit `7312e97` ✓
- Half-wire integral `max(left, right)` — commit `cb84f67` ✓
- `bcrit_pfaffian` real-witness wire — commit `65afb53` ✓ (NOTE: §3.3 reconciliation above; the witness is real but deferred to U13 for sweep coverage)
- Acceptance gate absolute window guard — commit `19c1da6` ✓
- Acceptance-gate ctest label — commit `3ead499` ✓
- Spec marked IMPLEMENTED + plan marked COMPLETE — commits `dfbdfa`, `df9954c` ✓

Verify each commit lands on PR #41's branch (i.e., is reachable from `origin/feat/bdg-validation-pass2` or will be brought in via rebase/cherry-pick).

## 7. Explicit deferrals (YAGNI per /ponytail ultra)

- **Chunk E (anti-ponytail trim, ~600 LOC)**: defer to a separate scoped PR after this lands. Items deferred: pfaffian.f90 Kitaev extraction, topological_analysis.f90 split, dead exports (`bdg_zero_energy_gap`, `csr_spectral_lorentzian_sum` public-stub, dead multi-k broadcast in `compute_spectral_function_bdg_wire`), `wire_pfaffian_witness` parallel-to-`_sweep`, `green_functions.f90:295-298` backstop, `bdg_hamiltonian.f90:222-227` manual CSR-to-dense. The codebase is functional; trim is cosmetic.
- **Wire Pfaffian B-sweep with periodic/Bloch BdG (U13)**: per existing P1 spec §Out-of-scope. Until U13 lands, wire_pfaffian witness is unavailable for sweep coverage; the 4th witness row in the reconciliation table is deferred (D5).
- **15 PRD Status footers beyond the 9 issue files** (per ce-doc-review scope-guardian: actual count is 15, not 11 as the earlier spec text said): INCLUDED in scope per Section 5.1. No deferral.
- **`defs.f90` enum extension (`bdq_spectral`)**: per CLAUDE.md Boundaries, this requires explicit re-confirmation in the PR description. Note the borderline status; if not re-confirmed, split into a one-line scoped PR. (No code change needed; just PR description update.)
- **`topological_analysis.f90` at 1832 lines (6× guideline)**: defer; file-size is a guideline, not a blocker.
- **`bdg_observables.f90` split-brain between factory and inline construction**: defer; works today; can be unified later without behavior change.
- **`asw_envelope`/`asw_single` (R6 P1 speculative)**: per commit `139b28b refactor(math): promote asw_evals/asw_single/asw_envelope to module-private`, already addressed locally. Verify reaches PR #41.

## 8. Implementation strategy

### 8.1 PR base + commit plan

PR #41 base is `origin/feat/bdg-validation-pass2`. The local `feat/bdg-p1-stabilization` carries ~34 commits. Strategy:

1. **Verify each local commit reaches PR #41 head**. Items that may not have made it: Sections 3.1 (C-1 row bug), 3.2 (verifier rewrite + wire-polarization emitter), 4.1 (cross-builder identity test). These are NEW in this addendum.
2. **Rebase strategy**: either fast-forward `feat/bdg-validation-pass2` to include all local commits (if the local branch is fast-forward-able), or rebase local-branch commits onto PR #41 head.
3. **Add new commits on top** for Sections 3.1, 3.2 (3a + 3b), 4.1, 4.3, 5.1, 5.2, 6.1, 6.2, 6.3.

### 8.2 Commit phasing on PR #41

**Phase A — P0 blockers (merge-blocker path):**
- A.1 [TDD-red] `tests/integration/test_slim_pfaffian_witness_projection.{py,sh}` — fails on real wire
- A.2 [TDD-green] `src/physics/topological_analysis.f90:1722, :1827` — multi-site scan band-major idx fix
- A.3a [TDD-red+green] `src/apps/main_topology.f90` + `src/io/outputFunctions.f90` — wire-polarization emitter + write helper, emit `output/majorana_polarization.dat`
- A.3b [rewrite] `tests/integration/verify_majorana_polarization.py` — call real `topologicalAnalysis`, parse real output, fail-loud on error
- A.4 [comment only] `tests/integration/verify_wire_bdg_topological.py:159-185` — inline "synthetic illustrative, U13 deferred" comment
- A.5 [verify] `tests/integration/test_lecture_13_acceptance_gate.sh` — confirm `19c1da6` + `3ead499` land; tighten tolerance 2.0 → 1.0
- A.6 [comment+label] `scripts/lecture_13_topological.py` — re-label fallback as "deferred to U13"; update reconciliation table from 4-witness to 3-witness

**Phase B — P1 test soundness:**
- B.1 [TDD-red+green] `tests/unit/test_cross_builder_hole_block_identity.pf` — shared-H₀ invariant + 2-fixture (k=0 + generic k with Bx≠0)
- B.2 [TDD-red+green] `tests/unit/test_wire_pfaffian_witness.pf:88-118` — drop zero-escape hatch (single-line `@assertTrue`)
- (B.3 Lutchyn-Oreg covered by existing P1 spec Phase 1.4 — verify only, no new commit unless missing)

**Phase C — Doc fixes:**
- C.1 24 .scratch archive files (9 issue + 15 PRD): add `**Status**: COMPLETE (2026-07-05)` footer
- C.2 lecture §13.5.1 cross-builder disclosure
- C.3 lecture §13.5.2 pairing iσ_y vs real ±1 clarification (per feasibility P2)
- C.4 lecture §13.7.4 slim-Pfaffian plot caption update
- C.5 existing P1 spec line-23 doc-drift fix (Section 6.1)
- C.6 ADR 0007 Layer D footnote (Section 6.2)
- C.7 ADR 0008 §3 reversal (Section 6.3) — load-bearing

**Phase D — Hygiene:** (none — §5.3 already covered by 7312e97)

**Total new commits**: 13 (A: 6, B: 2, C: 7, D: 0 — was previously stated as ~12; corrected per ce-doc-review coherence P3)
**Total LOC**: ~370 (Python ~230 [A.1 ~80, A.3b ~150], Fortran ~90 [A.2 ~5 + A.3a ~25 + B.1 ~30 + B.2 ~30], shell/markdown ~50)

### 8.3 Risk gates

- **(R1) PHS oracle stays 4/4 GREEN** — symmetric Peierls is untouched in this addendum
- **(R2) Acceptance gate stays 1.0 T tolerance (tightened from 2.0 T) + absolute window** — additive to current state, doesn't loosen
- **(R3) Cross-builder identity test passes for both k=0 and generic k=π/(2a) with Bx≠0** — if k=0 passes but generic k fails, Peierls phase handling is broken in one builder; if both pass, canonical-form convention is verified
- **(R4) Real-wire P_M verifier renders a real plot OR fails loudly** — no silent fallback (per D7); CI failure of real path → exit non-zero, lecture figure generation requires fix

## 9. Acceptance criteria (post-merge)

- `ctest --test-dir build -L unit` → **≥50 GREEN** (was 49; this adds ≥1 from `test_cross_builder_hole_block_identity.pf` per §4.1, removes 0; B.2's rewrite may add cases) — was previously stated as 49/49; corrected per ce-doc-review P0 (3 reviewers: feasibility, coherence, scope)
- `ctest --test-dir build -L 'regression|acceptance-gate'` → all GREEN
- `bash tests/integration/test_lecture_13_acceptance_gate.sh` → **3-witness passes WITHIN 1.0 T tolerance**, absolute window guard active
- `grep -c 'bdq_spectral' src/core/defs.f90` — already in diff, no change
- `tests/integration/verify_majorana_polarization.py` exits non-zero when real path fails (no synthetic fallback; per D7)
- Lecture §13.5.1, §13.5.2, §13.7.4 captions updated per Section 5.2
- 24 archive files (9 issue + 15 PRD) have `**Status**: COMPLETE` footers per Section 5.1
- Existing P1 spec line-23 doc-drift corrected per Section 6.1
- ADR 0008 §3 reversed per Section 6.3
- ADR 0007 Layer D footnote added per Section 6.2

## 10. Risks & out-of-scope

- **`wire_pfaffian_witness_sweep` U13 deferral** (PRD §Out-of-scope) remains: full wire Pfaffian B-sweep with periodic/Bloch BdG is not addressed. The 4th witness row in the reconciliation table is reserved for U13; acceptance gate is 3-witness until then (per D5).
- **`bdq_spectral` dense-QW support** (PRD §Out-of-scope) remains: only wire confinement supported; dense-QW deferred.
- **M=±1 on spinless chain structurally impossible** (existing P1 spec §1 line 13): production-grade M=±1 now lives on dense-QW rung, not spinless.
- **`csr_spectral_lorentzian_sum` shared-helper extraction** (existing P1 spec §Out-of-scope): deferred; internal-only in `spectral_bdg_wire.f90` is acceptable per Chunk E deferral.
- **Wire-polarization emitter** (`output/majorana_polarization.dat`) is NEW work introduced by this addendum (per ce-doc-review feasibility P0 finding §3.2): no producer exists today. Phase A.3a introduces it; Phase A.3b's verifier depends on it.

## 11. Cross-references

- **ADR 0007**: `docs/adr/0007-bdg-hole-block-canonical-convention.md` — canonical hole-block. Layer D gets a Peierls-symmetry footnote (Section 6.2).
- **ADR 0008**: `docs/adr/0008-bdg-p1-invariants.md` — Lutchyn-Oreg, pairing_sign PUBLIC, symmetric Peierls KEPT (after Section 6.3 reversal), minigap-pattern regression.
- **Existing P1 spec**: `docs/superpowers/specs/2026-07-01-bdg-p1-fix-design.md` — Phase 0–4 plans. Line 23 corrected per Section 6.1.
- **Implementation plan**: `docs/superpowers/plans/2026-07-01-bdg-p1-fix.md`.
- **PRD**: `.scratch/archive/bdg-majorana-validation/PRD.md` — Issues 00–11 + R4 cross-builder identity.
- **Memory**:
  - `codebase-doc-drift-prevention.md` — drives Status footer + spec refresh rules.
  - `feedback_no_coauthor_trailer.md` — no `Co-Authored-By:` trailer in any commit.
  - `feedback_pfunit_macro_single_line.md` — pFUnit `@assertEqual`/`@assertTrue` on single source line.

## 12. Self-review (per brainstorming skill, post-ce-doc-review)

- ✅ No TBD / TODO placeholders in body
- ✅ Internal consistency: D1 (Peierls kept) matches §6.1 (line-23 fix) AND §6.3 (ADR 0008 §3 reversal) AND source code at `bdg_hamiltonian.f90:406-419`. D5 (3-witness) matches §3.3 (deferred to U13) AND §3.4 (gate tolerance 1.0 T) AND §9 acceptance criteria.
- ✅ Scope: focused on P0 blockers + P1 soundness + Doc fixes + ADR amendment; Chunk E explicitly deferred.
- ✅ Ambiguity: each fix has named file:line + test file + LOC estimate. Acceptance criteria numeric.
- ✅ Resolved ce-doc-review P0 findings: §3.2 file path corrected (wire-polarization emitter added); §3.4 no-op work replaced with verification; 4-witness gate theater acknowledged (D5); ADR 0008 §3 amendment added (§6.3); §9 acceptance count corrected; cross-builder test recipe rewritten with shared-H₀ invariant.
- ✅ Resolved ce-doc-review P1 findings: env-gated synthetic fallbacks dropped (D7); config filename drift fixed (D8); pFUnit single-line rule explicit (§4.3); §3.3 decc6ba contradiction dropped (§3.3 keeps the rename); §5.1 IMPLEMENTED/COMPLETE unified to COMPLETE; §4.2 demoted to verification; §5.3 dropped (covered by 7312e97); §5.2 pairing-sign physics corrected; LOC breakdown accurate.
- ✅ Resolved ce-doc-review P3 findings: 12 → 13 commits; 11 → 15 PRDs.
- ✅ Three-way contradiction (ADR 0008 §3 + existing P1 spec line 23 + source code) reduced to two-way (existing P1 spec + source agree; ADR 0008 §3 amendment pending Section 6.3).