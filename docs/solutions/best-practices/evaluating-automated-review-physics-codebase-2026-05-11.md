---
title: "Evaluating automated code review on physics codebases: 16 of 18 findings were false positives"
date: "2026-05-11"
category: best-practices
module: automated_code_review
problem_type: best_practice
component: tooling
severity: medium
applies_when:
  - "Running automated code review (Codex or similar LLM reviewer) on Fortran physics or scientific codebases"
  - "Triaging large batches of review findings to separate real bugs from false positives"
  - "Reviewing flagged Hamiltonian construction, basis ordering, or k.p parameter code"
tags: [codex-review, false-positive, feast-eigensolver, peierls-phase, majorana-detection, band-major-indexing, bhz-model]
---

# Evaluating automated code review on physics codebases: 16 of 18 findings were false positives

## Context

After fixing 5 real bugs found by a first-round automated Codex review (documented in `topological-magnetic-index-logic-errors-2026-05-08.md`), a second Codex review of the same modules produced 18 findings (10 P1, 8 P2). Manual investigation revealed that 16 of 18 were false positives — the code was correct but the reviewer misread domain-specific conventions. Only 2 findings were real issues (FEAST sentinel value, Peierls phase limitation).

The contrast between rounds is instructive: the first review of buggy code found 5/5 real bugs. The second review of fixed code found 2/18 real bugs. Automated reviewers are effective at catching initial bugs but produce high false-positive rates on already-correct physics code.

## Guidance

### Systematic verification of each finding

For each automated review finding, trace the actual code at the cited line numbers:

1. **Read the current code**, not the finding's description. Reviewers may describe code from an earlier commit or misread branch logic.
2. **Verify indexing claims** against the project's convention. This codebase uses band-major `(band-1)*N + site` throughout. Reviewers consistently confused this with site-major `band + (site-1)*8`.
3. **Check for existing guards**. Many "double counting" findings miss explicit delta-computation patterns like `Vz_delta = Vz_opt - Vz_cfg`.
4. **Verify against tests**. If test assertions explicitly verify a coupling pattern (e.g., BHZ hopping connects specific components), that pattern is intentional, not a bug.

### Common false-positive patterns on physics code

| Claimed Bug | Actual Code | Why It's Wrong |
|---|---|---|
| "BdG pairing adds all 64 entries" | `col = 9 - row` (antidiagonal, 8/site) | Misread inner loop as nested over all bands |
| "Band-major indexing is wrong" | `(band-1)*N + site` used consistently | Assumed site-major without checking convention |
| "Zeeman added twice in BdG wire" | `Vz_delta = Vz_opt - Vz_cfg` guard | Missed the explicit difference computation |
| "Landau HT not cleared" | `HT = ZERO` at line 486 | Didn't read far enough into the subroutine |
| "Fu-Kane evaluated once" | Inside double loop at lines 1099-1117 | Misread loop nesting structure |
| "Input parser bdg: consumes EOF" | `case default` backspace pattern | Missed the early-exit mechanism |
| "BHZ B/D hops wrong component" | `mod(row,4)+1` is intentional model variant | Assumed standard BHZ; this code uses a cross-block variant |

### When automated findings ARE trustworthy

The two real findings shared characteristics that distinguish them from false positives:

1. **FEAST min_gap=0 on solver failure** (`main_topology.f90:539`): The code set a physically meaningful value (0.0) for a failure case, creating ambiguity. This is a real design flaw — zero is a valid BdG gap at topological transitions. Fixed with sentinel `-1.0`.

2. **Peierls phase on all off-diagonals** (`magnetic_field.f90:91-124`): The phase was applied indiscriminately. In gauge A=(0,0,Bx*y), transverse hops have zero z-displacement and should not receive a phase. Documented as known limitation for multi-column wires.

## Why This Matters

Without systematic verification, teams can waste hours investigating false positives or — worse — "fix" correct code based on reviewer misunderstandings. The 16 false positives each required reading the actual source code to confirm correctness. The time investment was worthwhile because it prevented introducing bugs into already-correct physics code, but the effort could have been reduced by understanding common false-positive patterns upfront.

The FEAST sentinel bug illustrates a subtle class of issues: code that produces a valid-looking but semantically wrong result. The gap `min_gap = 0.0` looks correct (zero gap is physically meaningful) but conflates solver failure with a real physical prediction. This class of bug is hard for automated reviewers to catch precisely because the output looks reasonable.

## When to Apply

- After running any LLM-based code review tool on scientific or physics code
- When review findings reference basis ordering, indexing conventions, or Hamiltonian construction
- When a finding claims a bug in code that has passing test coverage for the exact behavior in question
- Before modifying physics-critical code (Hamiltonian construction, FD stencils, eigensolvers) based on automated review

## Examples

### Verifying a BdG pairing claim

Codex claimed: "nested col=1,8 loop inserts pairing for every electron-hole combination."

Actual code at `bdg_hamiltonian.f90:220-231`:
```fortran
do i = 1, N
  do row = 1, 8
    col = 9 - row  ! antidiagonal: pairs band 1↔8, 2↔7, etc.
    global_row = 8 * N + (row - 1) * N + i
    global_col = (col - 1) * N + i
```

The inner loop runs `row = 1..8` with `col = 9 - row` (fixed per row), producing 8 entries per site, not 64. The "col = 1, 8" in the finding was a misread — there is no nested column loop.

### Verifying a double-counting claim

Codex claimed: "Zeeman added twice — ZB8bandGeneralized adds it, then add_zeeman_coo adds it again."

Actual code at `bdg_hamiltonian.f90:172-192`:
```fortran
Vz_delta = Vz_opt - Vz_cfg  ! difference between desired and already-applied
if (maxval(abs(Vz_delta)) > 1.0e-14_dp) then
  ! Only add the correction, not the full Zeeman
```

When the same B-field is used in both places, `Vz_delta = 0` and no additional Zeeman is added.

## Related

- `docs/solutions/logic-errors/topological-magnetic-index-logic-errors-2026-05-08.md` — first-round Codex review that found 5 real bugs in these same modules (prior to the fixes that made the second round mostly false positives)
- `docs/solutions/test-failures/lecture-test-pair-review-findings-2026-05-10.md` — L13 FEAST auto-window fallback for BdG mode (related FEAST eigensolver context)
