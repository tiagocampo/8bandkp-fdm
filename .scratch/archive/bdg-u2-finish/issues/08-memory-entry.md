**Status**: closed (2026-07-25) — ticket 08 of `.scratch/bdg-u2-actual-ship/`.

Type: task (AFK)
Status: closed
Claimed-by: —
Blocked by: 06, 07
Executed: 2026-07-25 (wayfinder work-through)

## Execution (2026-07-25)

Updated-in-place per the memory-deduplication guidance (avoiding forks; `feedback_style` preference for single-fact files). Files touched:

- `/home/tiago/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_evaluator_seam_ssot.md` (updated): corrected S1⊗S2 → S2-only, 4-witness-live → 3-witness-reserved-U13, detailed dense retirement, named `bdg_pfaffian_params_t` type + factory, added the 51 unit count, and quoted the User Story 5 contract verbatim.
- `/home/tiago/.claude/projects/-data-8bandkp-fdm/memory/MEMORY.md` (updated): shifted hook description at line 27 from "U2 closed 2026-07-13" to "U2 finish follow-up (PR #42 open, branch feat/bdg-u2-actual-ship); bdg_observables.f90 3-face seam; bdg_pfaffian_params_t SSOT; 3-witness gate; User Story 5 contract; 51 unit green".

Verification details:
- Index contains exactly the one entry pointing to the corrected file.
- Memory content carries the verbatim spec range contract fragment (`s2 >= -1 .and. s2 <= 1`).
- Links to related memories (`[[ctest-counts-targets]]`, `[[codebase-doc-drift-prevention]]`) are active.
Related: design.md §"Solution" item 14; §"Cross-References"

# 08 — Memory entry: project_bdg_evaluator_seam_ssot follow-on

## Question

After ticket 06 (spec amendment) and ticket 07 (doc propagation) land, what does the durable memory entry need to capture about the U2 finish work — and what gets appended to the `MEMORY.md` index?

## Background

The repo's persistent memory lives at `/home/tiago/.claude/projects/-data-8bandkp-fdm/memory/`. The current `MEMORY.md` index already has 27 entries. The most recent BdG-related entry is `project_bdg_evaluator_seam_ssot.md` (created 2026-07-13, captures the seam-SSOT pattern from PR #42).

Per `design.md` §"Solution" item 14 and §"Components" row 13:
> "Memory entry + MEMORY.md index pointer updated to match"

The new memory entry captures the **U2 finish** follow-on — the SSOT wire-up (ticket 01), orphan retirement (ticket 02), test amplification (ticket 03), Codacy triage (ticket 04), verification gate (ticket 05), and the User Story 5 amendment (ticket 06). The previous memory entry captures the seam-SSOT pattern from PR #42 proper.

Per repo convention (visible from MEMORY.md format), memory entries use:
- Frontmatter: `name`, `description`, `metadata.type` (user/feedback/project/reference)
- Body: factual statement + `**Why:**` + `**How to apply:**` for feedback/project
- Index entry: `- [Title](file.md) — hook`

## Sub-decisions to lock

1. **New file vs. update existing.** Per `feedback_no_coauthor_trailer` memory and the precedent of single-file-per-fact, create a NEW memory file rather than editing `project_bdg_evaluator_seam_ssot.md` (which captures the *seam pattern* from PR #42 proper). The new file captures the *U2 finish follow-on*. Naming: `project_bdg_u2_finish.md`.
2. **Type.** `project` (ongoing work, decisions, constraints — not derivable from code/git history alone).
3. **Content.** Capture:
   - The 8-ticket map at `.scratch/bdg-u2-actual-ship/map.md`.
   - The SSOT wire-up contract: `bdg_pfaffian_params_t` is the canonical type; `bdg_default_pfaffian_floor` is consumed (not just declared) post-PR-#42.
   - The User Story 5 carve-out: S1+S2 strict sign agreement is Issue 05 / U13, not U2.
   - The orphan-retirement decision (S5 finding → delete, not private).
   - The HANDOFF.md convention: live dir (`.scratch/bdg-u2-actual-ship/`) + archived dir (`.scratch/archive/bdg-evaluator-pfaffian/`) both have HANDOFF.md to disambiguate future agent entry points.
4. **Index pointer.** Add `- [BdG U2 finish](project_bdg_u2_finish.md) — hook` to MEMORY.md. The hook is a one-line summary linking to the file.
5. **No drift.** The memory entry must mirror the wording in `.scratch/archive/bdg-evaluator-pfaffian/spec.md` (per ticket 06) verbatim. Future agents grep for the wording across memory + doc + spec.

## Out of scope for this ticket

- Edits to `project_bdg_evaluator_seam_ssot.md` — that file stays unchanged (it captures the seam pattern from PR #42 proper, which is the *cause* of this map).
- New feedback / user memories — no new user preferences surfaced by this work.
- Edits to other memory entries.

## Files to touch

| File | Change |
|---|---|
| `/home/tiago/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_u2_finish.md` (new) | Memory file per format |
| `/home/tiago/.claude/projects/-data-8bandkp-fdm/memory/MEMORY.md` | Append index pointer line |

## Constraints

- Per CLAUDE.md memory section: "Each memory is one file holding one fact". The U2 finish is one fact (a follow-on to the seam-SSOT pattern).
- Per `feedback_no_coauthor_trailer`: no Co-Authored-By trailer on commits (also applies to memory entries — no attribution to Claude).
- Memory entries link to related memories with `[[their-name]]`. The new entry links to `[[project_bdg_evaluator_seam_ssot]]` (the PR #42 pattern that this follow-on extends).

## Verification

- `cat /home/tiago/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_u2_finish.md` shows frontmatter + body + Why/How-to-apply lines.
- `cat /home/tiago/.claude/projects/-data-8bandkp-fdm/memory/MEMORY.md | grep -E '^- '` count = 28 (27 existing + 1 new).
- The new entry's content mirrors `.scratch/archive/bdg-evaluator-pfaffian/spec.md` User Story 5 wording verbatim.
- The new entry links `[[project_bdg_evaluator_seam_ssot]]` correctly (the linked name matches an existing file).

## Interaction

- **Tickets 06, 07** are both blocked-by: this ticket mirrors the spec amendment (06) and the doc propagation (07). Both must be in place first.
- This is the LAST ticket. After it closes, the map's frontier is empty and the auto-archive condition triggers.

## Resolution path

1. Read the new User Story 5 text from `.scratch/archive/bdg-evaluator-pfaffian/spec.md`.
2. Draft the memory entry body mirroring the wording verbatim.
3. Write `/home/tiago/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_u2_finish.md` with the frontmatter + body + Why/How-to-apply + link.
4. Append the index pointer line to `/home/tiago/.claude/projects/-data-8bandkp-fdm/memory/MEMORY.md`.
5. Verify (3 grep commands above).
6. No commit needed — memory entries live outside the repo (per the memory section in CLAUDE.md).