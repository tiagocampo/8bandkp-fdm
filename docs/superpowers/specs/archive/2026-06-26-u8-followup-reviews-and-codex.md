---
title: "U8 Follow-up — review P1s + Codex stale-file P2 + PR27 I2 + doc/archive hygiene"
date: 2026-06-26
status: design (brainstorm-proposed; awaiting approval)
related:
  - docs/superpowers/specs/2026-06-21-u8-bdg-window-routing-design.md
  - docs/superpowers/plans/2026-06-21-u8-bdg-window-routing.md
  - docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md
  - docs/brainstorms/2026-06-14-bdg-majorana-validation-requirements.md
  - .scratch/pr27-review-fixes/REVIEW.md
  - https://github.com/tiagocampo/8bandkp-fdm/pull/40
---

# U8 Follow-up — review P1s + Codex stale-file P2 + PR27 I2 + doc/archive hygiene

## 1. Problem & reframe

The U8 fix landed on `feat/bdg-u8-window-routing` (PR40) but the deliverable is
partially complete: a comprehensive adversarial review surfaced 4 P1s + 8 P2s +
10 P3s, the Codex bot posted an inline P2 on the PR diff itself (stale
`output/bdg_eigenvalues.dat` on sentinel return), and the PR27 audit's sole
open item (I2 silent `b_field%components` fallback) is still in U8's attack
surface. The PR40 title is "fix(bdg): wire BdG window routing + validate
guard (U8)" — the follow-up commits stay on this branch and PR.

PR40 also has a Codacy "Not up to standards" bot comment (2 critical + 1
medium + 1 minor Security) but the underlying issue descriptions are not
visible in the GitHub API. Plan: ship the in-scope fixes, then re-trigger
Codacy via `@codacy review` to surface specifics; address real findings as
a C7 commit if needed.

## 2. Scope (in this PR)

**Fix commits (C1–C6, in order):**

| # | Severity | Files | Description |
|---|---|---|---|
| C1 | P1 | `src/io/input_parser.f90`, `tests/integration/test_validate_rejects_bad_configs.sh` | Replace 3 silent `if (stat /= 0) ... = 0.0_dp` fallbacks on `b_field%components(1/2/3)` (L451-456) with `check_optional_stat`. Add V13 rejection test (non-numeric `components` entry). Closes PR27 I2. |
| C2 | P1 | `src/core/defs.f90`, `tests/integration/test_topology_validate_rejects.sh` (new T4) | New `validate_semantic` check next to existing BdG window guard: reject BdG config with `B_vec(1)=0 .and. B_vec(2)=0 .and. abs(B_vec(3)) > 1e-12` (axial B only — Peierls returns silently). Error message: "BdG requires transverse B (Bx or By nonzero) for Peierls orbital coupling". |
| C3 | P1+P2 | `src/apps/main_topology.f90` | Convert sentinel branch (L525-530) from `print * + result%min_gap = -1.0_dp` to `error stop` with descriptive message (CLAUDE.md convention). Before the error stop, **delete any pre-existing `output/bdg_eigenvalues.dat`** so downstream scripts don't read stale spectra after a μ-in-gap run (Codex P2). |
| C4 | P2 | `src/core/defs.f90` | Add `real(kind=dp), parameter :: BDG_WINDOW_BOUND = 1.0_dp  ! eV — rejection ceiling for BdG solver windows` at module scope (L55 area). Replace inlined `1.0_dp` at L956 with the parameter reference. Update error message to include the symbol name. |
| C5 | P3 | `src/apps/main_topology.f90` | Drop `auto_compute_energy_window` from the `use eigensolver, only:` list (L14) — no longer called after C3. Keep exported from `eigensolver.f90` (still consumed internally). |
| C6 | P2 | `tests/integration/verify_wire_bdg_topological.py`, `tests/integration/validation_universe.yml` | Extend T2 to assert `min_gap == -1.0` within 1e-6 tolerance (currently only checks warning text — a regression that emits warning + returns `0.0` would pass). Add `minigap` cell to `validation_universe.yml`: `observable: minigap`, `geometry: wire`, `material: InAs`, `tier: required`, `reference: topological_transition`. |

**Follow-on commits (C7–C10) in this PR — folded in to prevent doc drift:**

| # | Severity | Files | Description |
|---|---|---|---|
| C7 | P1 | `src/apps/main_topology.f90:1242,1263-1264,1269-1280,1288-1291` | Route `eval_wire_bdg_gap_app` through `apply_solver_window` with the same user-override semantics as `run_bdg_wire` (review P1.1). Replace `B_vec=[0,0,B_val]` hardcode at L1242 with `B_vec=[B_val,0,0]` so the gap-sweep sits in the Peierls-active regime (transverse B — fix the design's second root cause). Replace `error stop` at L1272 + L1280 with the `-1` sentinel + warning (matching `run_bdg_wire`'s C3 path). Move cleanup (L1288-1291) into a single subroutine called from both success and sentinel paths so the memory-leak path is fixed. |
| C8 | P2 | `tests/unit/test_bdg_hamiltonian.pf` | Add PHS cross-check: build H_BdG with `B_vec=[1.0, 0, 0]` for both wire and QW builders, assert `\|C·H·C^{-1} + H\|/\|H\| < 1e-10` (review P2.2 / brainstorm R3). Closes the R3 verification gap that the parent validation plan's U4/U5 require. |
| C9 | P2/P3 | `docs/superpowers/specs/2026-06-21-u8-bdg-window-routing-design.md`, `docs/superpowers/plans/2026-06-21-u8-bdg-window-routing.md` | Resolve plan/design self-contradictions (review P2.3, P2.4, P3.9): (a) design §4.1 sentence committing to `eval_wire_bdg_gap_app` routing — replaced with status footer noting C7 landed it; (b) design §6 deferred-to-U10 sentence for `eval_wire_bdg_gap_app` — removed (no longer deferred); (c) design §4.5 / plan Task 5 "correct validation-plan U8 text" — removed as no-op (already corrected in earlier commit); (d) plan Tasks 3/4/5 — add status footer listing the commits that landed them. |
| C10 | P3 | `docs/UBIQUITOUS_LANGUAGE.md`, `src/apps/main_topology.f90:1224` | Add 3 new entries to `UBIQUITOUS_LANGUAGE.md` under a new "Numerics" section: **Gershgorin bound**, **sentinel gap value** (`-1.0_dp` for "no BdG modes found"), **FD-Nyquist tail** (cross-ref `docs/solutions/best-practices/2026-06-21-fd-nyquist-spurious-tail.md`). Rename `eval_wire_bdg_gap_app` → `eval_wire_bdg_gap` (drops meaningless `_app` suffix). Add docstring in `verify_wire_bdg_topological.py` noting the OMP thread cap relies on `OMP_NUM_THREADS` env override per CLAUDE.md ctest gotcha. |

**Out of scope (honest follow-ups — no doc-claim, no drift risk):**

- `compute_spectral_function_wire` auto-window bypass (review P2.8) — U9 per parent validation plan; no doc claim about its routing today
- `run_bdg_wire` 200-line function-size refactor (review P2.9) — pure refactor; no doc claim about its size

Both file as `U9` and `tech-debt:run_bdg_wire-size` tickets respectively at PR40 merge.

## 3. Sequencing rationale

C1 first because it's the most isolated (parser-only, no physics logic touched)
and closes the longest-standing open item (PR27 I2 from 2026-06-22 audit).
C2 next because it adds a sibling validation guard next to the existing BdG
window guard, so the diff is contextually grouped. C3 changes sentinel
behavior, which C6 then locks in via test. C4 is a trivial refactor that
touches the same file as C2 — could merge into C2; kept separate per CLAUDE.md
convention "one concern per commit". C5 is YAGNI cleanup tied to C3 (the
import is unused after fallback removal). C6 last among the original C1–C6
because it depends on C3's new sentinel semantics.

C7 follows C3 because it adopts the same sentinel pattern across both BdG
wire sites. C8 is independent — pure test addition — can run anywhere; placed
after C7 so all behavioral changes precede all test additions. C9 doc cleanup
runs after C7 because the §4.1/§6 contradiction resolves only once C7 lands.
C10 final: terminology and naming touch the same sites C7 modified, so doc-
and naming-updates go last to reflect the final state.

## 4. Doc & archive moves (after C1–C6 land and ctest is green)

**`.scratch/pr27-review-fixes/`** — entire directory moves to
`.scratch/archive/pr27-review-fixes/`. After the move:
- `REVIEW.md` updated to `Status: COMPLETE` (was INCOMPLETE)
- `prd.md` and `issues/01..04` retained as historical record
- The `mv` runs as a separate commit (no other content changes mixed in)

**`docs/superpowers/specs/2026-06-21-u8-bdg-window-routing-design.md`** — keep
in place; the spec is still the authoritative reference for the wire-BdG
window routing decision. Add a status footer:
```
Status (2026-06-26): Implementation complete via commits a4ade9d, 4c445c7,
1567d90, f2840c4 (PR40) + 6 follow-up commits (PR40 C1–C6).
```

**`docs/superpowers/plans/2026-06-21-u8-bdg-window-routing.md`** — add status
footer noting which tasks landed in which commits (Tasks 3/4/5 → `a4ade9d..c6dc762`;
Tasks 1/2/6 → PR40 follow-up commits).

**`docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md`** — U8
section gets `Status: DONE (PR40)` footer. Parent plan still has U1–U12.

**`docs/brainstorms/2026-06-14-bdg-majorana-validation-requirements.md`** —
add reprioritization annotation to R3 (per review P3.7): "Status: reprioritized
by U1 (2026-06-15) — real but tertiary; U4/U5 still required for correctness".

**`docs/UBIQUITOUS_LANGUAGE.md`** — add three new entries under a new
"Numerics" section: **Gershgorin bound**, **sentinel gap value** (`-1.0_dp`
sentinel for "no BdG modes found"), **FD-Nyquist tail** (cross-reference
`docs/solutions/best-practices/2026-06-21-fd-nyquist-spurious-tail.md`).

**`docs/CONTEXT.md`** — no change (no new ambiguity).

**`CLAUDE.md`** — no change (no new boundary; existing approval gate on
`defs.f90` derived types applies to C2/C4 but neither adds a derived type).

**`src/physics/AGENTS.md`, `tests/integration/AGENTS.md`** — no change (no
new module/concern added).

## 5. Memory entry

One new project memory under `~/.claude/projects/-data-8bandkp-fdm/memory/`:

- **Name:** `codebase-doc-drift-prevention`
- **Description:** "Every PR touching behavior must update the design spec
  status, plan status, and AGENTS.md before merge; archive closed PRDs to
  `.scratch/archive/`."
- **Why:** Past drift observed — U8 design §4.5 / plan Task 5 said the
  validation-plan U8 text needed correcting, but the validation-plan U8
  section had ALREADY been corrected by an earlier commit, making the work
  a no-op. Plan Tasks 3–5 listed as RED→GREEN when commits `a4ade9d..c6dc762`
  had already landed them. Drift produces confused future agents.
- **How to apply:**
  1. Before opening a PR that touches existing behavior: grep the related
     spec/plan for stale claims about what the PR will do; resolve
     contradictions inline (mark "aspirational" or "no-op because
     already done").
  2. After committing a PR: append a `Status (YYYY-MM-DD): committed via
     <sha list>` footer to the corresponding spec AND plan.
  3. After closing a PRD in `.scratch/pr<N>-review-fixes/`: `mv` the directory
     to `.scratch/archive/pr<N>-review-fixes/` and update any `REVIEW.md`
     to `Status: COMPLETE`.

The memory entry also adds a one-line pointer in `MEMORY.md` (the index).

## 6. Risks & approval gates

- **`defs.f90` approval gate** (C2, C4) — CLAUDE.md flags as approval-gated
  for derived-type changes. C2 adds a check function; C4 adds a `parameter`.
  Neither modifies a derived type. Same pattern as PR40's existing validate
  guard commit `4c445c7`. Pre-flagged in the PR body per PR40's "Approval
  note" precedent.
- **Codacy Security findings** — specifics unknown. Mitigation: ship C1–C6,
  re-trigger via `@codacy review`, address as C7 if real.
- **ctest cost** — 126 current + 2 new (V13, T4) = 128/128 target.
  `OMP_NUM_THREADS=$(( $(nproc)/4 ))` per CLAUDE.md gotcha.
- **Behavioral visibility** — C3 changes sentinel from silent `-1` to
  `error stop`. Users who currently see a `-1` return on misconfigured BdG
  runs will see an exit instead. Documented in commit message; intentional
  per CLAUDE.md "no silent corrections" principle.

## 7. Validation

After all commits + archive + memory:

1. `ctest --test-dir build -j4` → target 129/129 (126 current + V13, T4,
   and the new PHS test from C8).
2. `git log feat/bdg-u8-window-routing --not main` → 20 commits (10
   existing + 10 new C1–C10).
3. PR40 inline threads: Codex P2 (stale file) resolved via C3; Codacy
   specifics surfaced via `@codacy review`.
4. `.scratch/pr27-review-fixes/` no longer exists; `.scratch/archive/`
   contains it.
5. `MEMORY.md` and the new `codebase-doc-drift-prevention` file exist
   in the memory dir.
6. Design §4.1/§6 contradiction resolved; plan Tasks 3/4/5 status footer
   present; brainstorm R3 annotated.
7. `UBIQUITOUS_LANGUAGE.md` has the 3 new Numerics entries.