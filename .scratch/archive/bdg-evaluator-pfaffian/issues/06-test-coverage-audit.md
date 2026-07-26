**Status**: COMPLETE (2026-07-13) — ticket 06 of `.scratch/archive/bdg-evaluator-pfaffian/`.

Type: task
Status: resolved
Resolved: 2026-07-13
Blocked by: 01, 02

# 06 — Test coverage audit

## Question

For the Pfaffian-plugged evaluator, document the existing test coverage and decide what to add.

Audit:

- Does `tests/unit/test_bdg_evaluator.pf` exist? The parent plan §U2 says it was supposed to be new. The `src/physics/AGENTS.md` inventory doesn't list it. Search `tests/unit/` to confirm.
- What does `test_bdg_hamiltonian.pf` already cover for `eval_bdg_point`? (Likely nothing — the seam is in a different module.)
- What does `test_wire_pfaffian_witness.pf` cover today? Per Phase 24 follow-up ticket A, it has a strict `@assertTrue(s1 == s2 .and. s1 /= 0)` marked `@todo U13`. Map's "Not yet specified" fog: whether (d+ii) unblocks it.
- What does `test_pfaffian.pf` / `test_kitaev_majorana.pf` cover for U3's wrapper? (U3 is shipped.)
- What does `test_topological_analysis.pf` cover for `wire_pfaffian_witness`?
- What does the integration side cover? `tests/integration/test_lecture_13_acceptance_gate.sh` and `verify_majorana_polarization.py` interact with the witness.

Decide:

- Whether `test_bdg_evaluator.pf` needs to be created (likely yes — the seam's own tests live where the seam lives).
- What new tests are needed for the Pfaffian plug-in (per the API shape settled in ticket 01): at minimum, a synthetic-spectrum test that the slim witness returns the expected sign on a known-band matrix.
- Whether `test_wire_pfaffian_witness.pf`'s `@todo U13` becomes "now passing under slim witness" or stays `@todo` for the Bloch-Pfaffian future.
- Whether `tests/unit/test_bdg_hamiltonian.pf`'s existing convention-pinning tests need updates for the new seam shape.

Deliverable: a coverage table at ticket resolution (existing test → what it covers → keep / update / delete), plus a one-paragraph list of new tests to add.

## Constraints

- Read-only audit; no code changes.
- pFUnit `@assertEqual`/`@assertTrue` must be single-line (CLAUDE.md + memory feedback `feedback_pfunit_macro_single_line`).
- Cite file:line for every claim.

## Out of scope for this ticket

- Writing the new tests (downstream execution after this map's tickets close).
- The acceptance-gate wiring (ticket 05).
- The plan + backlog close-out (ticket 07).

## Answer

### Summary

Existing coverage is **adequate for what is currently shipped**, but **gaps remain for the two seam siblings locked by ticket 01**: `eval_bdg_pfaffian_witness_csr` (the CSR entry the gate row consumes) and `eval_bdg_kitaev_majorana` (the QW+Kitaev rung entry). `test_bdg_evaluator.pf` exists and covers `eval_bdg_point` thoroughly; the strict sign-agreement test in `test_wire_pfaffian_witness.pf` is a fixture-level failure (per ticket 02 §9's research, confirmed here) and does not unblock via seam routing. Five new unit tests land now; one existing test gets a non-diagonal fixture to unblock the strict sign-agreement assertion independent of U13.

### Premise corrections to the ticket

1. **`test_bdg_evaluator.pf` exists** at `tests/unit/test_bdg_evaluator.pf:1` (5.2k bytes, 8 tests). It is registered in `tests/CMakeLists.txt:182-186`. The AGENTS.md inventory omission flagged by ticket 04's doc-drift sub-finding still stands (this file is missing from the listed pFUnit files but is wired).

2. **`test_topological_analysis.pf` does not exist** (no file matches in `tests/unit/`). The `wire_pfaffian_witness` dense-path tests live in `test_wire_pfaffian_witness.pf` (registered at `tests/CMakeLists.txt:325-329`).

3. **Two Kitaev test files coexist**: `test_kitaev_majorana.pf` (11k, 10 tests, harness-builds-fixture pattern, diagonal-in-(c,c†) — wrapper M=±1 NOT distinguished, only M=0 gap-closure testable) and `test_kitaev_strict.pf` (1.5k, 2 tests, uses the separate `kitaev_bdg_fixture_2band` from `pfaffian.f90:483` — M=±1 ARE distinguished). Both are registered (`tests/CMakeLists.txt:289-298`). The strict tests call `kitaev_majorana_number` directly from `pfaffian` module, not the future seam sibling.

4. **`test_bdg_hamiltonian.pf` has zero references to `eval_bdg_point`** (verified by grep). The seam is in `bdg_observables.f90`, not `bdg_hamiltonian.f90`; the 41k test file covers the Hamiltonian builder + PHS oracle + spectrum, not the per-point evaluator. No convention-pinning update is needed.

5. **Integration side, complete inventory**:
   - `tests/integration/test_lecture_13_acceptance_gate.sh` — 4-witness agreement (currently 3-witness when Pfaffian row is `"(deferred to U13)"`; per ticket 05 becomes 4-witness always).
   - `tests/integration/verify_majorana_polarization.py` — `@todo U13` SKIP for BLOCKING-EMPIRICAL on polarization file. Ticket 05's resolution tightens the SKIP precondition (parse colormap, exit non-zero if z2==-1 row missing).
   - `tests/integration/test_slim_pfaffian_witness_projection.py` — parses `output/wire_slim_pfaffian_witness.dat` (file does not exist; test exits 0 with deferred label). This is the integration test for the CSR sibling's downstream consumer.

### Coverage table

| Existing test | File:line | What it covers | Verdict | Rationale |
|---|---|---|---|---|
| `test_bdg_evaluator.pf` (8 tests) | `tests/unit/test_bdg_evaluator.pf:11-163` | `eval_bdg_point`: minigap, near-zero count, invariant flag, threshold override, empty spectrum, purity; `q_zero_tol` floor + abs | **keep + extend** | Comprehensive coverage of the seam's pure-function contract. The 8 tests pin the seam signature, edge cases, and the `bdg_default_near_zero_frac` / `bdg_default_min_threshold` SSOT behavior. Add 1-2 tests for `bdg_default_pfaffian_floor` (the SSOT ticket 02 surfaces) when it ships. |
| `test_wire_pfaffian_witness.pf:60-72` smoke | `tests/unit/test_wire_pfaffian_witness.pf:61-72` | `wire_pfaffian_witness` returns valid sign in {-1,0,+1} on synthetic 16x16 | **keep** | Smoke guard; pins no-crash + range contract. |
| `test_wire_pfaffian_witness.pf:74-105` strict sign-agreement | `tests/unit/test_wire_pfaffian_witness.pf:75-105` | Spec §4.3 strict: `s1 == s2 .and. s1 /= 0` (no zero-escape) | **update + unblock now** | The Fortran-comment `@todo` at line 96 does NOT unblock via seam routing — fixture-level failure (`build_synthetic_bdg_16` is diagonal-in-(c,c†) → Pf=0 by construction). Per ticket 02 §9, option (a): non-diagonal synthetic fixture (~30 LOC, adds conduction-band off-diagonal coupling that crosses the Nambu seam) lands now. The strict assertion becomes GREEN independent of U13. |
| `test_wire_pfaffian_witness.pf:107-130` non-diagonal single-particle | `tests/unit/test_wire_pfaffian_witness.pf:108-130` | `wire_pfaffian_witness` range check on matrix with off-diagonal conduction-band coupling | **keep** | Already exercises the off-diagonal path. The strict-assertion fix can reuse this fixture pattern. |
| `test_pfaffian.pf` (10 tests) | `tests/unit/test_pfaffian.pf:17-222` | Real/complex Pfaffian (2x2, 4x4, random skew), `Pf² = det`, Parlett-Reid path n=14,16, skew guard correction | **keep** | The underlying math. Adequate. No BdG-specific path here. |
| `test_kitaev_majorana.pf` (10 tests) | `tests/unit/test_kitaev_majorana.pf:69-293` | Harness shape, Hermiticity, gap closure at \|μ\|=2t, M=0 detection, wrapper validity in {-1,0,+1}, M=±1 cannot be distinguished on diagonal-in-(c,c†) harness | **keep + extend** | The harness caveat (lines 184-201) is well-documented. The wrapper validity tests are correct for what they exercise. **Add 2-3 direct calls to the seam sibling `eval_bdg_kitaev_majorana`** (when ticket 01 ships) — same range checks, plus the gap-closure case via the seam entry. |
| `test_kitaev_strict.pf` (2 tests) | `tests/unit/test_kitaev_strict.pf:11-42` | M=-1 (topological, \|μ\|<2t) and M=+1 (trivial, \|μ\|>2t) via `kitaev_bdg_fixture_2band` from `pfaffian.f90:483` | **keep** | The only tests that distinguish M=±1 — uses the 2-band fixture (different from the diagonal-in-(c,c†) harness in `test_kitaev_majorana.pf`). These must be replicated through the seam sibling when ticket 01 ships, otherwise the seam sibling's M=±1 discrimination is untested. |
| `test_majorana_polarization.pf` (5+ tests) | `tests/unit/test_majorana_polarization.pf:34-?` | Sticlet P_M, half-wire integral discriminator, MZM profile peak at wire ends | **keep** | Polarization path (U6) is orthogonal to the Pfaffian plug-in. No update. |
| `test_bdg_hamiltonian.pf` (~30 tests) | `tests/unit/test_bdg_hamiltonian.pf:24-1038+` | BdG builder, Hermiticity, PHS, k=± symmetry, full vs eval-wire spectrum, QW vs wire hole-block | **keep** | Hamiltonian builder + PHS oracle. Zero `eval_bdg_point` references. No update needed. |
| `test_lecture_13_acceptance_gate.sh` | `tests/integration/test_lecture_13_acceptance_gate.sh:1-172` | 4-witness B_crit agreement (Pfaffian row deferred); false-PASS line removal; Issue 07 regression green | **update (ticket 05 downstream edits)** | Ticket 05 already locks the gate-row spec (colormap-extracted z2==-1 row at mu≈0.6601, tol 1.0 T). No additional test change beyond ticket 05's scope. |
| `verify_majorana_polarization.py` | `tests/integration/verify_majorana_polarization.py:1-149` | Polarization on real wire eigensolve | **update (ticket 05 downstream edits)** | Ticket 05 tightens SKIP precondition (parse colormap, exit non-zero if z2==-1 row missing). The `@todo U13` text on the polarization gap stays. |
| `test_slim_pfaffian_witness_projection.py` | `tests/integration/test_slim_pfaffian_witness_projection.py:1-?` | Slim Pfaffian on real wire B-sweep (currently deferred — `output/wire_slim_pfaffian_witness.dat` not emitted) | **keep + activate on U13** | The integration test for the gate row's CSR consumer. Deferred label is honest; activated when the .dat emitter ships. |

### New tests to add (this PR)

The two seam siblings locked by ticket 01 (`eval_bdg_pfaffian_witness_csr` and `eval_bdg_kitaev_majorana`) are the **only** new functions whose behavior is currently untested at the unit level. Five new tests across three files, plus one existing-test update:

1. **NEW: `tests/unit/test_bdg_pfaffian_witness_csr.pf` (4 tests)** — pins the CSR sibling the gate row consumes. Same harness conventions as `test_wire_pfaffian_witness.pf` (Fortran-comment `@todo` style, single-line `@assertTrue` per CLAUDE.md, `use bdg_observables, only: eval_bdg_pfaffian_witness_csr`). Tests:
   - **smoke**: returns valid sign in {-1,0,+1} on a constructed CSR with off-diagonal coupling (reuses the line-120 pattern).
   - **topological_sign**: pre-known Pf sign on a hand-built CSR with explicit omega = τ_y ⊗ I_2 multiplication (Pf = +a verified).
   - **gap_closure_returns_zero**: matrix with min\|E\| = 0 returns s2_sign = 0 (Pf vanishes).
   - **floor_activation**: Pf below `bdg_default_pfaffian_floor` returns 0 even when raw Pf ≠ 0.

2. **NEW: `tests/unit/test_bdg_kitaev_majorana.pf` (3 tests)** — pins the seam sibling for the QW+Kitaev rung. Distinct from `test_kitaev_majorana.pf` because it tests the **seam entry point**, not the underlying `pfaffian.f90` wrapper. Uses the same harness conventions. Tests:
   - **wrapper_returns_valid_range**: invokes seam sibling on `build_kitaev_bdg`-style matrices, checks result ∈ {-1,0,+1}.
   - **gap_closure_returns_zero**: same gap-closure semantics as `test_kitaev_majorana.pf:139-158`, but via the seam sibling.
   - **defensive_input**: empty/single-momentum input returns 0 (mirrors `test_kitaev_majorana.pf:161-181`).

3. **UPDATE: `tests/unit/test_wire_pfaffian_witness.pf:31-58` synthetic fixture** — per ticket 02 §9 option (a), ~30 LOC. Add off-diagonal coupling between conduction bands (7-8) and valence (1-4) that crosses the Nambu seam, so Pf ≠ 0 by construction. Strict sign-agreement at `:104` becomes GREEN. The `@todo` Fortran comment at `:96` is removed (replaced with "see synthetic-fixture convention below").

4. **UPDATE: `tests/unit/test_bdg_evaluator.pf` (2 tests)** — pin the new SSOT `bdg_default_pfaffian_floor` from `bdg_observables.f90` (ticket 02 §3). Tests:
   - **floor_value_pinned**: `bdg_default_pfaffian_floor == 1.0e-12_dp`.
   - **floor_consumed_by_seam_siblings**: param override via `bdg_eval_params_t`-equivalent path.

5. **UPDATE: `tests/unit/test_kitaev_majorana.pf` (1 test)** — add a single test calling the seam sibling `eval_bdg_kitaev_majorana` directly (not the underlying `kitaev_majorana_number`). Pins the seam entry point without duplicating the 10 existing tests. Test:
   - **seam_sibling_returns_valid_range**: invokes `eval_bdg_kitaev_majorana(H_k_array, k_par)` on a topological-phase harness matrix, checks ∈ {-1,0,+1}.

### Decision: `test_topological_analysis.pf`

**Do NOT create** a `test_topological_analysis.pf` umbrella file. The `wire_pfaffian_witness` dense-path coverage lives in `test_wire_pfaffian_witness.pf` (registered at `tests/CMakeLists.txt:325-329`); the CSR sibling coverage will live in the new `test_bdg_pfaffian_witness_csr.pf`. A `test_topological_analysis.pf` umbrella would mix the slim witness with the Chern/Z2/BHZ paths and dilute coverage scope — net a layering violation.

### Decision: `@todo U13` lifecycle

The Fortran-comment `@todo replace synthetic fixture with real wire BdG once U13 lands` at `test_wire_pfaffian_witness.pf:96` is **removed** by update (3) above (option-a fixture makes the strict assertion GREEN now). The pFUnit `@todo` macro does not exist outside its own grammar; the comment is documentation only.

The `@todo U13` lines in `verify_majorana_polarization.py:5, :102, :108` (the BLOCKING-EMPIRICAL polarization gap) **stay** — they are unrelated to the Pfaffian plug-in (different observable: Sticlet polarization, not invariant). Ticket 05's downstream edits tighten the SKIP precondition (colormap z2==-1 row required) but do not delete the `@todo U13`.

### Net downstream

- **4 files touched** in `tests/unit/`: 2 new (`test_bdg_pfaffian_witness_csr.pf`, `test_bdg_kitaev_majorana.pf`), 2 updated (`test_wire_pfaffian_witness.pf` fixture + `@todo` removal, `test_bdg_evaluator.pf` floor SSOT).
- **1 file touched** in `tests/unit/` for the seam-sibling-pinning minor addition (`test_kitaev_majorana.pf`).
- **2 new tests wired** in `tests/CMakeLists.txt` (~8 LOC per `add_pfunit_ctest` block, following the pattern at `:182-186` for `test_bdg_evaluator`).
- **0 integration test changes** beyond ticket 05's downstream edits (already in scope).
- **0 Fortran source changes** — this ticket is test-side only.

### Unblocks

- **Ticket 07** (plan close-out) — the test coverage map is finalized; the deck-sweep doc update can cite this ticket's coverage table.
- **Post-PR** — the gate-row's underlying CSR sibling has a dedicated unit test; future regressions in `eval_bdg_pfaffian_witness_csr` surface at the unit-test rung, not the 4-witness agreement rung.