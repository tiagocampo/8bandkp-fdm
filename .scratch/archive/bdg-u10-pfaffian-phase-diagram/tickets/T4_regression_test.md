# T4 — Non-flat colormap regression test (shell + ctest)

## Question (AFK, task-ticket)

What is the smallest, robust regression test that locks the "non-flat" property of `output/z2_phase_diagram.dat` into CI, and exercises the slim Pfaffian wireup end-to-end?

Per the closer integration-test pattern (`test_topology_sweep.sh` at `tests/integration/test_topology_sweep.sh`):
- New shell script `tests/integration/test_wire_pfaffian_phase_diagram_noflat.sh`.
- Copy the (new post-T2) wire config to a temp workdir, run `topologicalAnalysis`, assert:
  1. `output/z2_phase_diagram.dat` exists and is non-empty.
  2. `awk` counts distinct `parts[3]` (z2 column) values → count ≥ 2.
  3. Existential row check: at least one row with `z2 == -1` (per Q3 native) at `mu ≈ μ_c ± δ_mu`.
  4. `output/wire_slim_pfaffian_witness.dat` (post-T1) exists and has ≥ 1 row.
- Add `add_test(...)` to `tests/integration/CMakeLists.txt` (or wherever the existing `test_topology_sweep` is registered) with `COVERAGE: observable=phase_diagram_nonflat geometry=wire material=InAs`.

**Why shell, not pFUnit?** Per `feedback_pfunit_error_stop_uncatchable`: the wire BdG path has many `error stop` exits; a pFUnit @test that hits one kills the whole ctest process. Shell + exit-code is the safe surface — same pattern as `test_wire_bdg_topological_2d.sh` and `test_topology_sweep.sh`.

## Acceptance

- [ ] `tests/integration/test_wire_pfaffian_phase_diagram_noflat.sh` exists, mode 755.
- [ ] Asserts outputs non-flat (≥ 2 distinct z2 values), tests git status clean.
- [ ] Asserts per-B `wire_slim_pfaffian_witness.dat` exists.
- [ ] ctest registration (`-L` label) wired; runs `ctest -L wire-pfaffian-phase-diagram` PASS.
- [ ] `# COVERAGE: observable=phase_diagram_nonflat geometry=wire material=InAs` annotation present.
- [ ] Status footer updates: `BACKLOG.md` Phase 26 row, parent plan `status_note`, `REVIEW.md` row 79.

## Blocks / blocked-by

- **Blocks**: T5 (gate precondition update needs the new colormap schema to assert on).
- **Blocked by**: T1 (per-B proxy file must exist), T2 (window must span both phases), T3 (native schema).

## Type
task (AFK).

## Branch
`feat/bdg-u10-pfaffian-phase-diagram`.

## Claimed by
Tiago de Campos (2026-07-27 — fresh session, continuation from T3 handoff)

## Strategy (grilled 2026-07-27)
Path B selected: probe FEST-vs-μ envelope first by widening μ-window in a probe fixture to `nMu=11, μ∈[0.650, 0.670] eV` and observing where FEST fails (U8 fail-fast). The actual sibling fixture's μ-window tightens to the widest safe envelope the probe reveals. Fallback to Path A (widen B-grid to nB=21, keep μ tight) if probe reveals μ envelope is too narrow for non-flat colormap.

## Probe fixture (in scratch)
`tests/regression/configs/wire_inas_gaas_bdg_topological_probe.toml` (temporary, deleted after envelope recorded) — siblings `_phase.toml` and `_2d.toml` with widened μ-grid only; everything else preserved.

## Probe findings (2026-07-27)
- Probe 1 (μ∈[0.650,0.670], nMu=11): FEST fails — exit 1 (eval_wire_bdg_gap line 1380 fail-fast).
- Probe 2 (μ∈[0.655,0.665], nMu=11): FEST fails somewhere mid-grid (verified running rebuild).
- Probe 3 (μ∈[0.657,0.663], nMu=11, nB=6): SUCCESS — 66 cells, all z2=-1, "no phase transitions detected".
- Probe 4 (μ∈[0.657,0.663], nMu=11, nB=21): SUCCESS — 231 cells, all z2=-1 across all B including 2.75T near B_crit≈2.85T.
- Probe 5 (probe 4 + Bx=1 instead of B_vec=[0,0,0]): SUCCESS — 231 cells, all z2=-1.

**Root cause:** the slim projected Pfaffian (per `eval_bdg_pfaffian_witness_csr`, S2 strategy on the projected band-7,8 subspace) emits a clear `s2_sign=-1` (topological) at every (B, μ) point inside the FEST envelope. The `s2_sign=0` (closure) regime requires either the full Bloch-Pfaffian with periodic Peierls-twist (U13 BLOCKING-EMPIRICAL) or an FEST window that captures the gap-closure micro-region (touches ADR 0005 window authority — explicitly out of scope per handoff). Therefore **T4's TDD-green assertion `count_distinct(parts[3]) >= 2` is not satisfiable inside the slim Pfaffian's signal regime**. Single-value colormap is the physical correct output for the current Pfaffian implementation, NOT a bug.

## Resolution (per second grilling — 2026-07-27)
**T4 closed BLOCKED — do NOT manufacture non-flat colormap by forcing the slim Pfaffian to emit `s2_sign=0` (closure) or to flip between `-1` and `+1` inside the FEST envelope.** The slim projected Pfaffian (S2 strategy, post-`eval_bdg_pfaffian_witness_csr` seam sibling) correctly classifies every (B, μ) point inside the FEST envelope as `s2_sign=-1` (topological, clear sign); the `s2_sign=0` (closure) and `s2_sign=+1` (trivial) regimes only emerge under conditions that are physically out of reach of the current implementation:
- **Closure regime** requires the full Bloch-Pfaffian with periodic Peierls-twist (U13, BLOCKING-EMPIRICAL per parent plan §U13 and CLAUDE.md Known Issues).
- **Trivial regime** requires either a wider FEST window that captures the gap-closure micro-region (touches ADR 0005 window authority — explicitly out of scope per handoff) or a geometry/material where the band inversion never occurs.

Forcing either regime to make the colormap non-flat would **bug the physics that's working**. User direction (2026-07-27, second grilling): "don't bug code that is breaking the physics" — the slim Pfaffian's single-value output IS the correct physical answer inside its regime.

**Action taken (none on source):**
- ✅ Probe fixtures `tests/regression/configs/wire_inas_gaas_bdg_topological_probe{2,3,4,5}.toml` deleted (only the working sibling `_phase.toml` retained). `wire_inas_gaas_bdg_topological_probe.toml` kept as documentation of the FEST envelope (μ∈[0.657, 0.663] is the widest safe envelope; outside of it FEST fails U8 fail-fast).
- ✅ `input.toml` restored to canonical stub (`confinement="bulk"`, top-level).
- ✅ This ticket closed BLOCKED with the pointer to U13 below.
- ✅ T5 (gate precondition update) is unaffected — gate reads `z2 == -1` natively per T3, no flat-colormap assumption in T5's logic.

**Routing:**
- Non-flat z2 colormap assertion: **deferred to U13** (full Bloch-Pfaffian + periodic Peierls-twist).
- BACKLOG.md Phase 26 row: add line "T4 closed BLOCKED 2026-07-27 — non-flat colormap assertion deferred to U13 per slim-Pfaffian signal-regime constraint (probe4 + probe5 evidence: 231/231 cells s2_sign=-1 inside FEST envelope). Do NOT widen FEST window (ADR 0005 out-of-scope) or manufacture closure regime (would bug working physics)."
- Parent plan §U10 status footer: append "T4 closed BLOCKED 2026-07-27 → deferred to U13; colormap non-flat assertion is physically unsatisfiable inside the slim Pfaffian's signal regime".
- REVIEW.md row 79: mirror the same line.

## Out of scope (locked by this resolution)
- Any code change to `eval_wire_bdg_gap`, `eval_bdg_pfaffian_witness_csr`, the FEST-window selection, the slim Pfaffian's S2 strategy, or the z2-column emitter. All current code in this scope is correct and unchanged.

## Related
- T5 (gate precondition update) — unblocked, claim next.
- PR #43 review — parallel track, covers T1b + T2 + T3 (T4 not in PR scope per this close-out).
- U13 (full Bloch-Pfaffian + periodic Peierls) — destination of the deferred non-flat colormap assertion.
