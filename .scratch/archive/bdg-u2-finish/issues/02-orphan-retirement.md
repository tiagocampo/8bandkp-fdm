**Status**: closed + executed — ticket 02 of `.scratch/bdg-u2-actual-ship/`.

Type: task (AFK)
Status: closed (2026-07-13) + executed (2026-07-25)
Claimed-by: Claude (session 2026-07-13)
Resolved-by: Claude (session 2026-07-13) via wayfinder
Executed-by: Claude (session 2026-07-25) via wayfinder work-through
Blocked by: 01 ✓ (closed 2026-07-14, executed 2026-07-25)
Related: design.md §"Solution" item 4; §"Architecture" Region 2; §"Components" row 2 (S5)

# 02 — Orphan retirement: delete dense `wire_pfaffian_witness`

## Question

After PR #42 ships (which migrated the production wire-rung site at `main_topology.f90:1371` to the seam sibling `eval_bdg_pfaffian_witness_csr`), is the dense `wire_pfaffian_witness` subroutine at `topological_analysis.f90:1641-1707` still called by anything on the post-PR-#42 `main` branch — and if not, delete it?

## Background

PR #42 (commits `4801dd5..fe66980`) replaced the production call site `wire_pfaffian_witness_sweep` → `eval_bdg_pfaffian_witness_csr` (commit `06886b8`). The dense `wire_pfaffian_witness` subroutine (which had a different signature `(H_bdg, n_full, s1_sign, s2_sign)` taking a 2-D array rather than CSR) was retained because `tests/unit/test_wire_pfaffian_witness.pf` exercises the dense path directly (per ticket 06 of the previous map: "dense-path tests live in `test_wire_pfaffian_witness.pf`" — no `test_topological_analysis.pf` umbrella was created).

The review-finding label **S5** ("Orphan Public Symbol") flags `wire_pfaffian_witness` as no-longer-called by production code, only by tests. The ticket 01 answer of the previous map locked "delete, not private" as the resolution (option 5a per design.md §"Architecture" Region 2). Per CLAUDE.md YAGNI: re-introduce when a real need shows.

## Sub-decisions LOCKED (resolution 2026-07-13)

**Sub-decision 1 — Delete dense `wire_pfaffian_witness` (1641-1707).**
- Rationale: 0 production callers post-PR #42 (audit grep: only definition site + `public ::` export at `:40` remain in `src/`). Test-only consumers + YAGNI = delete per design.md §"Architecture" Region 2 option 5a.

**Sub-decision 2 — DELETE both helpers `s1_project` and `s2_project` (NOT keep).**
- The original ticket's sub-decision 2 ("`s1_project` is used by `wire_pfaffian_witness_sweep`"; "`s2_construct` is used by `wire_pfaffian_witness_sweep`") was empirically false. Audit grep:
  - `s1_project` called only at `topological_analysis.f90:1673` (inside the dense subroutine).
  - `s2_project` called only at `topological_analysis.f90:1691` (inside the dense subroutine).
  - The sweep at `:1720-1805` performs inline CSR extraction (`H_bdg_csr%rowptr(idx(i))`, `H_bdg_csr%colind(k)`, `H_bdg_csr%values(k)`); it does NOT call either helper.
  - Also: original ticket names `s2_construct`; actual subroutine is `s2_project` at `:1875`. Stale ticket info.
- Decision: DELETE both helpers. They are private (not in `public ::` exports) and only consumers die with the dense subroutine.

**Sub-decision 3 — DELETE all 3 dense-path test cases (NOT migrate).**
- The original ticket's sub-decision 3 premise ("if 0 hits on dense subroutine, test file already only exercises sweep") was empirically false. Audit grep on `tests/unit/test_wire_pfaffian_witness.pf`:
  - L21: `use topological_analysis, only: wire_pfaffian_witness` (import)
  - L115, L144, L172: 3 call sites to the dense subroutine (`call wire_pfaffian_witness(H_bdg, 16, s1, s2)`)
- Options evaluated:
  - **(A) DELETE** (chosen): tests 1 (`smoke`) and 3 (`nondiagonal_single_particle`) assert only `s1 ∈ [-1, +1] .and. s2 ∈ [-1, +1]` range checks — already covered by `tests/unit/test_bdg_evaluator.pf` via the parametrized floor. Test 2 (`sign_agreement_synthetic`) asserts `s2 /= 0` — covered by `test_bdg_pfaffian_params_default_uses_ssot` (ticket 01 verification). YAGNI-aligned; no physics value lost.
  - (B) MIGRATE: requires dense-fixture-to-CSR conversion + seam sibling returns `s2_sign` only (not `s1, s2` pair) → tests 1 + 3 `s1` range checks become no-ops. ~30-40 LOC test infrastructure overhead, zero added coverage.
- Decision: DELETE the test file entirely.

**Sub-decision 4 — Drop `public :: wire_pfaffian_witness` export at `topological_analysis.f90:40`.**
- Rationale: orphan export after sub-decision 1. The ticket's "Files to touch" omits this; corrected here.

**Sub-decision 5 — Update 2 stale doc/comment references.**
- `src/physics/AGENTS.md:39` — inventory row mentions `wire_pfaffian_witness`; remove the mention, preserve other Issue-07/U10 content.
- `tests/CMakeLists.txt:334-341` — comment header + `add_pfunit_ctest(test_wire_pfaffian_witness ...)` block; remove both.
- The ticket's "Files to touch" omits both; corrected here.

**Sub-decision 6 — NO private stub. NO dense wrapper sibling. NO backport of the helpers.**
- Per ticket's constraints section + YAGNI. "Delete, not private."

## Out of scope for this ticket

- The 3 hard-coded `1.0e-12_dp` literals — ticket 01 owns `:1781` (wire-up); `:1680` + `:1698` retire with sub-decision 1's deletion of the dense subroutine.
- Test changes (SSOT error-stop, s2_sign tightening, dedup, meta-test) — ticket 03.
- Module split (`topological_analysis.f90` ≈ 1740 lines) — separate concern; BACKLOG Phase 18.

## Files to touch

| File | Change |
|---|---|
| `src/physics/topological_analysis.f90` | Delete `wire_pfaffian_witness` (1630-1707 incl. banner), `s1_project` (1807-1864), `s2_project` (1866-1921); drop `public :: wire_pfaffian_witness` at `:40`. |
| `src/physics/AGENTS.md` | Edit L39 — remove `wire_pfaffian_witness` mention; adjust line-count budget (~1740 → ~1640). |
| `tests/CMakeLists.txt` | Remove comment header (`:334-336`) + `add_pfunit_ctest(test_wire_pfaffian_witness ...)` block (`:337-341`). |
| `tests/unit/test_wire_pfaffian_witness.pf` | DELETE the file entirely. |

## Constraints

- TDD per CLAUDE.md: after deletion, `ctest -L unit` MUST pass with the running total = `52 - 3 = 49` tests (per ticket 05's adjusted baseline; ticket 01's happy-path placeholder adds 1 → ticket 05 final = 50).
- Per CLAUDE.md YAGNI: this is a *delete*, not a *private*. Don't replace the public surface with a private stub.
- Per ADR 0001: no polymorphic types. The deletion is straight line removal.
- Doc-drift per `codebase-doc-drift-prevention.md`: the AGENTS.md edit carries the verbatim wording so reviewers can grep for it (ticket 07 carries the doc-propagation sites).

## Verification

- `cmake --build build` clean (no new warnings).
- `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -j4 -L unit` → 49/49 PASS (52 baseline − 3 deleted dense tests; ticket 01's happy-path placeholder adds 1 → ticket 05 will re-lock to the running total).
- `grep -nrE '\bwire_pfaffian_witness\b' src/ tests/` → 0 hits (sweep variant `wire_pfaffian_witness_sweep` survives and the regex matches it; verify separately with `grep -nrE '\bwire_pfaffian_witness_sweep\b' src/ tests/` ≥ 1 hit at definition).
- `grep -nE '1\.0e-12_dp' src/physics/topological_analysis.f90` → 0 hits (3 pre-ticket-01 → 1 post-ticket-01 parameterized `:1781` → 0 post-ticket-02 after dense subroutine deletion).
- `wc -l src/physics/topological_analysis.f90` → ~1643 (1923 baseline − ~280 deleted: 67 dense body + 56 `s1_project` + 47 `s2_project` + ~110 banners/blank-line separators; ±5 line tolerance).
- `grep -nE 'eval_bdg_pfaffian_witness_csr' src/apps/main_topology.f90` → ≥1 hit at `:1374-1375` (production site live).

## Interaction

- **Ticket 01** is blocked-by: this ticket's deletion happens AFTER ticket 01's SSOT wire-up is in place. Reason: if ticket 02 lands first and `wire_pfaffian_witness_sweep` still has the 3 hard-coded literals (because ticket 01 hasn't parameterized them yet), the SSOT promise is still broken. Strict ordering: 01 → 02 → 03 → 05.
- **Ticket 03** writes tests that pin the seam sibling's `s2_sign` contract. After ticket 02, those tests only need to exercise the sweep variant.

## Execution (2026-07-25, wayfinder work-through)

Decision-locked 2026-07-13; **executed 2026-07-25**. Code landed in working tree on `feat/bdg-u2-actual-ship` (uncommitted):

- `src/physics/topological_analysis.f90`: deleted dense `wire_pfaffian_witness` (was `:1641-1707`), `s1_project` (was `:1809-1864`), `s2_project` (was `:1866-1921`); dropped `public :: wire_pfaffian_witness` at `:40`. File 1924 → 1738 lines (deleted ~186 incl. 2 banner-comment blocks + the dense body + both helpers). The `wire_pfaffian_witness_sweep` overload survived and is now the only wire Pfaffian symbol.
- `tests/unit/test_wire_pfaffian_witness.pf`: `git rm`'d entirely.
- `tests/CMakeLists.txt`: removed the `# Issue 07 (U10)` comment header + `add_pfunit_ctest(test_wire_pfaffian_witness ...)` block (was `:334-341`).
- `src/physics/AGENTS.md:39`: replaced `wire_pfaffian_witness` mention with `wire_pfaffian_witness_sweep`, added a retirement note, updated the line-count budget `~1740 → ~1640`.

**Verified:**
- `cmake -G Ninja -B build ...` (reconfigured — CMakeLists changed) clean; `cmake --build build` clean (289/289 targets, no warnings).
- `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -j4 -L unit` → **51/51 PASS** (was 52 baseline − 1 ctest target `test_wire_pfaffian_witness`, which carried 3 `@test` subroutines as 1 ctest — ctest counts targets not `@test`s, so the drop is 1, not 3).
- `grep -nrE '\bwire_pfaffian_witness\b' src/ tests/` → 0 bare hits (matches exclude the `_sweep` variant); `wire_pfaffian_witness_sweep` survives at definition `:1650` + public export `:40` + delegation call in `bdg_observables.f90:204`.
- `grep -nE '1\.0e-12_dp' src/physics/topological_analysis.f90` → 0 hits (all 3 pre-PR literals retired: `:1781` parameterized by ticket 01, `:1680`/`:1698` deleted with the dense subroutine here).
- `eval_bdg_pfaffian_witness_csr` production call site live at `main_topology.f90:1379`.

**Count-arithmetic correction for ticket 05:** the ticket body's "52 − 3 = 49" prediction was wrong — it conflated pFUnit `@test` subroutines with ctest targets. The actual drop is **52 → 51** (one `add_pfunit_ctest` block = one ctest target, regardless of how many `@test` subroutines it bundled). Ticket 05 must re-lock the running unit total to **51** (pre-ticket-03 hardening), not 49/50.

**Out-of-scope carry-forward:** the 2 remaining hard-coded `1.0e-12_dp` literals this ticket owned are retired (verified 0 hits); no further literal retirement is owed.