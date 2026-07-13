# Documentation: CLAUDE.md Update + ADR 0003

**Type:** AFK
**Blocked by:** #02 (C1 must land first — docs must reflect the final code state)

## What to Build

Update project documentation to reflect the architectural changes from C6 and C1. Commit as `docs: update CLAUDE.md and add ADR 0003`.

### ADR 0003

Verify that `docs/adr/0003-landau-b-sweep-stays-in-main.md` exists with the correct content documenting the decision to leave the Landau B-sweep inline in `main.f90` rather than absorbing it into `simulation_setup`. The ADR should reference the three considered options (A: leave in main, B: add method to setup, C: polymorphic sweep) and explain why A was chosen. If the file already exists from a prior session, verify its content is accurate against the final implementation.

### CLAUDE.md Updates

Update the following sections:

1. **Line 176 — module dependency graph:** Change the `simulation_setup.f90` description from:
   ```
   simulation_setup.f90     (simulation orchestration: H-build, eigensolve, k-point sweep; used by all four executables)
   ```
   to note that it handles all 4 confinement modes:
   ```
   simulation_setup.f90     (simulation orchestration: H-build, eigensolve, k-point sweep for all 4 confinement modes; used by all four executables)
   ```

2. **Verify consistency:** Scan CLAUDE.md for any references to `setup_alloc_sweep` (should be none — it was never documented), any references to Landau init being inline in main.f90 (update if found), and any references to the old 3-mode simulation_setup pattern.

3. **Verify ADR index:** Check that `docs/adr/` numbering is sequential (0001, 0002, 0003).

## Acceptance Criteria

- [ ] `docs/adr/0003-landau-b-sweep-stays-in-main.md` exists with correct content
- [ ] CLAUDE.md line 176 updated to mention 4 confinement modes
- [ ] No stale references to deleted code (`setup_alloc_sweep`, `add_zeeman_coo`) in any docs
- [ ] ADR numbering is sequential
- [ ] `cmake --build build` succeeds (docs change should not affect build, but verify)
- [ ] Single commit with message `docs: update CLAUDE.md and add ADR 0003`

## Blocked by

- #02 (C1 must land first)

## User Stories Covered

- #15: ADR 0003 documents the B-sweep exclusion decision
