# T6 — Per-B min-|Pf| extraction Lutchyn-Oreg criterion (formalize the gate)

## Question (HITL, grill-ticket)

The 4th gate witness `bcrit_pfaffian` is currently "deferred to U13" in `lecture_13_topological.py:246`. T1 produces an open-chain projected Pfaffian `|Pf|_min` per B (cheap), but the formal Lutchyn-Oreg criterion for the topological phase boundary is:

  B_crit² = B_Z² + (μ - μ_c)² / (g μ_B / 2)²

where `μ_c` is the Rashba chemical potential at the transition and `g` is the effective g-factor. The "min-|Pf|" per B in the wire BdG sweep is a *proxy* for B_crit — it captures the B where the open-chain projected Pfaffian magnitude is smallest (near the closure), but the full Bloch-Pfaffian at PHS-invariant momentum (U13) is the rigorous invariant.

Per Q4 acceptance, T1 emits the per-B proxy. T6 formalizes the gate's consumption:
- Update `lecture_13_topological.py:222-256` to read `output/wire_slim_pfaffian_witness.dat` and compute `bcrit_pfaffian = argmin(|Pf|_min)`.
- Cross-check: `bcrit_pfaffian` should equal `bcrit_curve` (~2.85 T) within tolerance (e.g. ±0.5 T) once the slim Pfaffian drives z2 across the closure.
- Tolerance derivation: the open-chain projected Pfaffian vs. the strict Bloch-Pfaffian at k=0 (or PHS-invariant k) can differ by O(1/g) where g is the gap; for the canonical wire fixture, ε-scale expected.
- Update the gate's status string to "approximation (open-chain projected Pfaffian; full Bloch-Pfaffian strict invariant deferred U13)".

## Acceptance

- [ ] `bcrit_pfaffian` numeric value emitted in `bcrit_pfaffian` row of the L13 status table.
- [ ] Status string annotation updated.
- [ ] Cross-rung consistency check: `bcrit_pfaffian ≈ bcrit_curve ≈ bcrit_2d` within ±0.5 T (or derived tolerance — grill session).
- [ ] U13 strict invariant remains the formal source of truth — U10 proxy is a wiring milestone, not a replacement.

## Blocks / blocked-by

- **Blocks**: U11 (lecture revamp with cross-rung consistency table).
- **Blocked by**: T1 (proxy file), T4 (regression test ensures the file is non-empty).

## Type
grilling (HITL).

## Branch
`feat/bdg-u10-pfaffian-phase-diagram`.

## Claimed by
Tiago de Campos (2026-07-27 — fresh session, post-T5 closure handoff, T6 is only open frontier ticket per map.md:72).

## Empirical baseline (2026-07-27, sibling `_phase.toml` fixture, nB=6 B∈[0,5]T)

```
B=0.00000000E+00 |Pf|=4.00000000E-08
B=1.00000000E+00 |Pf|=4.00000000E-08
B=2.00000000E+00 |Pf|=4.00000000E-08
B=3.00000000E+00 |Pf|=4.00000000E-08
B=4.00000000E+00 |Pf|=4.00000000E-08
B=5.00000000E+00 |Pf|=4.00000000E-08
```

|Pf|_min is **uniformly at the floor (4.0E-08)** across the entire B-grid. The slim Pfaffian is saturated at the `bdg_default_pfaffian_floor` SSOT before the writer sees it. **argmin extraction at lecture_13_topological.py:229 would return B=0T (first occurrence of the floor) — not a physical phase boundary.** The corresponding gap from `z2_phase_diagram.dat` varies smoothly (3E-3 at B=0 → 1.4E-4 at B=3 → 3.6E-3 at B=5) and DOES track the closure near B≈2.85T, but that data path is `bcrit_curve`/`bcrit_2d`, not `bcrit_pfaffian`.

This is the empirical ground truth the grill must address. T6 ticket body assumed ±0.5 T tolerance was meaningful; the baseline shows the Pfaffian proxy has no discriminating signal at all in the FEST envelope.

## Resolution (in progress — grilled 2026-07-27)

<!-- being filled in by the grilling session; do not commit until Q5 (final review) is closed -->

### Q1 — Core finding: argmin of |Pf|_min is non-discriminating (saturated at floor)

Three resolutions offered:
- (a) drop `bcrit_pfaffian` entirely, keep status string "approximation", no numeric value
- (b) like (a) but add at-floor detection warning
- (c) re-derive `bcrit_pfaffian` from gap data in `z2_phase_diagram.dat`

User direction: try (c) first, fall back to (a) if (c) doesn't work.

**Q1 attempted: (c) — argument rejected before writing code.**

Mechanically: `bcrit_pfaffian = argmin(per_B_min_gap)` reading the same `z2_phase_diagram.dat` file. On the baseline fixture this returns **B=3.0T** — identical (to ΔB grid resolution) to `bcrit_2d` already emitted at line 178-196 of lecture_13_topological.py via the same gap-min strategy.

`bcrit_curve` is sourced from a different file (`verify_wire_bdg_topological.py` 1D curve at canonical `wire_inas_gaas_bdg_topological.toml`, fixed μ=0.6601 eV) and returns **2.85T** — different from `bcrit_2d` because of grid resolution, but `bcrit_2d == (c) bcrit_pfaffian` is structurally identical by construction.

**Failure mode:** (c) makes the 4th witness numerically equal the 2nd witness. The L13 status table reads `bcrit_pfaffian = bcrit_2d = 3.0T` and the reader is misled into thinking `bcrit_pfaffian` is an independent Pfaffian-derived witness when it's the same gap-min reading of the same file. TOLERANCE_BCRIT_RANGE check passes trivially because the 4th witness adds no new information. U10 destination line promised "4th gate witness approximate-not-deferred" — duplicating `bcrit_2d` doesn't satisfy that.

**Decision: fall back to (a) per user direction.** `bcrit_pfaffian` stays at `None` with the T5 approximation status string; the 3-witness gate (wire_curve + wire_2d + qw_dense) remains the actual acceptance surface; the slim Pfaffian witness file stays in the output as informational.

### Q2 — Implementing (a): floor-saturated branch

User choice: **(a.1)** — print a separate WARN for the at-floor case, keep the canonical T5 approximation status string at the print site (consistent with file-absent and file-empty cases).

**Decision: add a third "file-present-and-matched-but-all-at-floor" branch at `lecture_13_topological.py` lines 222-250.**

- Detect: all `pf_mags.values()` equal `bdg_default_pfaffian_floor` SSOT within tolerance (or just bitwise equal to the floor — empirically `4.0E-08`).
- Set `bcrit_pfaffian = None`.
- Print `WARN: |Pf|_min saturated at floor <value>; bcrit_pfaffian marked 'approximation (open-chain projected; full Bloch-Pfaffian deferred U13)'`.
- The "approximation" status string at line 259-260 stays canonical for all three failure modes (file absent, file empty, file all-at-floor).

### Q3 — Implementation details (relative-variance detection, hardcoded constant, no gate test change)

User direction: "do what is best based on the codebase and the problem at hand."

Decisions:

- **3.1 Detection = relative variance, not floor-equality.** Hardcode-free detection: `pmax = max(pf_mags.values()); pmin = min(pf_mags.values()); saturated = pmax > 0 and (pmax - pmin) / pmax < _PFAFFIAN_DEGENERACY_TOL`. This catches all saturation regimes (current FEST precision floor `4.0E-08`, hypothetical `bdg_default_pfaffian_floor = 1.0e-12_dp` SSOT, any future clamping) without hardcoding a specific floor value. Cross-check: `bdg_default_pfaffian_floor` is the Fortran SSOT at `src/physics/bdg_observables.f90:47` and lives in real(kind=dp). No Python-side equivalent exists. Hardcoding a numeric floor in Python would be fragile to SSOT changes; relative variance is dimensionless and survives any future SSOT update.

- **3.2 Module-level constant `_PFAFFIAN_DEGENERACY_TOL = 1e-6`** with comment `# relative-variance threshold; below this the slim Pfaffian per-B profile is numerically degenerate and cannot distinguish B_crit from noise (current FEST window saturates at ~4.0E-08, well above the bdg_default_pfaffian_floor SSOT 1.0e-12_dp)`. Tolerance `1e-6` is loose enough to catch float-rounding noise but tight enough to reject any physically meaningful per-B variation.

- **3.3 No acceptance-gate test change.** Existing `test_lecture_13_acceptance_gate.sh:86` regex `^[0-9.+-]+$` correctly handles the new "at-floor → approximation label" path (non-numeric → PFAFFIAN_DEFERRED=1 → 3-witness fallback). Existing `test_slim_pfaffian_witness_projection.py` parses signs ≠ 0 from the file directly — the file is unaffected (still 6 non-zero rows at the precision floor), so the test still passes. **The lecture 13 script's `BCRIT wire_pfaffian` line emission (lines 252-260) already covers all three failure modes (numeric → numeric; None → approximation label); the at-floor branch fits the existing `else: bcrit_pfaffian = None` slot.**

- **3.4 Update T5 doc-strings** at lecture_13_topological.py:222 (`# Per ADR 0008 §4 + spec §5.3: invoke topologicalAnalysis with slim Pfaffian mode over a B-grid, parse output, emit B_crit = argmin |Pf(kz, B)|.`) to clarify that argmin is meaningful only when the per-B profile is non-degenerate. Add: `# If the profile is numerically degenerate (relative variance < 1e-6, FEST precision floor in current regime), fall back to approximation — argmin would otherwise return the first occurrence of the floor, not a phase boundary.`

## Resolution (2026-07-27, completed)

T6 closed. Single source-code commit on `feat/bdg-u10-pfaffian-phase-diagram` (Pick-D + Pick-A — T6 has no Pick-D bookkeeping, the prior PR #43 review is unaffected and the standing decisions file from T5 carries forward):

**Commit (forthcoming):** `feat(bdg): U10 T6 — Pfaffian proxy at-floor detection (lecture 13 fallback)`

### Source change

`scripts/lecture_13_topological.py`:
1. Added module-level constant `_PFAFFIAN_DEGENERACY_TOL = 1e-6` after imports (lines 43-52), with comment explaining the relative-variance threshold rationale and SSOT cross-reference (`bdg_default_pfaffian_floor = 1.0e-12_dp` in `src/physics/bdg_observables.f90:47`).
2. Inserted "file-present-but-all-at-floor" branch in `section_wire_rung()` between the existing "pf_mags populated → argmin" branch (line 229) and the "file-present-but-empty" branch (line 236). Detection: `pmax = max(pf_mags.values()); pmin = min(pf_mags.values()); saturated = pmax > 0 and (pmax - pmin) / pmax < _PFAFFIAN_DEGENERACY_TOL`. On saturation: `bcrit_pfaffian = None` + WARN `WARN: |Pf|_min numerically degenerate (max={pmax:.3e}, min={pmin:.3e}, rel_var={...:.1e}); bcrit_pfaffian marked 'approximation (open-chain projected; full Bloch-Pfaffian deferred U13)'`.
3. Updated doc comment at lecture_13_topological.py:220-225 to clarify that argmin is meaningful only when the per-B profile is non-degenerate; fall back to approximation if numerically degenerate.

No other files touched (no Fortran change, no test-file change, no fixture change).

### Acceptance items (all closed)

- [x] `bcrit_pfaffian` numeric value emitted in `BCRIT wire_pfaffian <B>` of the L13 status table — yes, when the witness file per-B |Pf|_min profile is non-degenerate (Test 3 confirmed with synthetic varying file: `BCRIT wire_pfaffian 3.000`).
- [x] Status string annotation updated — approximation label unchanged from T5; new degenerate-detect WARN adds diagnostic context (preserves canonical T5 phrase for all three failure modes per Q2 a.1).
- [x] Cross-rung consistency check — `[Omitted per design resolution: T6 ticket body assumed `bcrit_pfaffian` would track `bcrit_curve ≈ bcrit_2d` within ±0.5 T; resolution (c)→(a) proved the proxy is non-discriminating in the current FEST envelope, so the cross-check tolerance is moot while the proxy is degenerate. The check is structurally available — `bcrit_pfaffian` can be compared to `bcrit_curve` when the slim Pfaffian is no longer saturated — but the comparison is now gated on the proxy being non-degenerate, not on a numeric tolerance against a saturated baseline.]`
- [x] U13 strict invariant remains the formal source of truth — T6 is an approximation gate-wiring milestone, not a replacement for U13 (per ticket body acceptance §4).
- [x] Empirical baseline (`wire_inas_gaas_bdg_topological_phase.toml`, nB=6 B∈[0,5]T): all 6 per-B |Pf|_min = 4.0E-08 (FEST precision floor) → relative variance = 0.0e+00, at-floor branch fires correctly with new WARN.
- [x] Three failure-mode branches verified end-to-end: file absent (T5 behavior preserved) → approximation WARN; file saturated → new degenerate WARN; file varying → numeric `BCRIT wire_pfaffian <B>`.
- [x] **P2 follow-up (2026-08-08):** added a fourth failure mode — file all-zero (closure regime, every `s2_sign=0` so `best_pf_abs=0`; `pmax == 0` short-circuits the relative-variance guard). Fix lands `if pmax == 0:` branch in `section_wire_rung()` with the same WARN shape; the test harness `tests/integration/test_pfaffian_degeneracy_detection.py` pins all four regimes. Failure-mode list is now: file absent / file empty-match / file saturated / file all-zero / file varying.

### Verification

- 51/51 unit green (41 s)
- 3/3 wire_bdg regression green: `regression_wire_bdg_topological` (111 s) + `regression_wire_bdg_topological_2d` (312 s) + `regression_wire_bdg_topological_phase` (560 s)
- 2/2 slim_pfaffian regression green: `regression_wire_slim_pfaffian_witness` (309 s) + `regression_slim_pfaffian_witness_projection` (0.1 s)
- 1/1 lecture 13 acceptance gate green (`test_lecture_13_acceptance_gate.sh` ALL PASS with 3-witness B_crit agreement 0.8 T ≤ 1.0 T tolerance)
- `python3 -c "import lecture_13_topological"` clean

### Doc-drift stamps (per `codebase-doc-drift-prevention`)

- `docs/plans/BACKLOG.md` Phase 26 row → amended: "T6 (2026-07-27) closed: slim Pfaffian per-B `|Pf|_min` saturation detection... | U2 MERGED, U10 in progress, T4 BLOCKED→U13, T6 closed"
- `docs/plans/REVIEW.md` row 79 → amended: "T5 closed 2026-07-27 (full cleanup). T6 closed 2026-07-27 (per-B `|Pf|_min` saturation detection, grilling resolution (c)→(a))"
- `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` §U10 footer → replaced T5-active + T6-unblocked blocks with T5-closed + T6-closed blocks (per-file wording)

### Out of scope (locked by this resolution)

- BHZ unit-test churn (test_edge_states.pf:80, test_z2_invariant.pf:74, test_phase_diagram.pf:178 pinning BHZ to {0,1}) — out-of-scope per T3 standing decision; T6 deliberately did not touch BHZ.
- Cross-rung consistency check between `bcrit_pfaffian` and `bcrit_curve` / `bcrit_2d` (±0.5 T tolerance from original ticket body) — moot while the slim Pfaffian proxy is degenerate; structurally available when the proxy becomes non-degenerate (which requires U13 full Bloch-Pfaffian with PHS-invariant k).
- Widening FEST window to capture the gap-closure micro-region (would manufacture the trivial regime `s2_sign=+1` for `bcrit_2d` cross-check, but touches ADR 0005 window authority — explicitly out of scope).
- Lecture 13 status table 4-witness label sweep beyond line 327 — U11 dependency, not T6.
- `bcrit_2d` deduplication (would require changing the 2D colormap extraction to source from `wire_slim_pfaffian_witness.dat` instead of `z2_phase_diagram.dat` — neither advisable nor needed; current behavior is correct: each witness reads its native data source).

### Frontier after this resolution

U10 fully shipped on `feat/bdg-u10-pfaffian-phase-diagram`:
- T1a closed (seam extended with optional `best_pf_abs` out-arg)
- T1b closed (per-B min-|Pf| proxy producer)
- T2 closed (config-driven writer paths + sibling phase fixture)
- T3 closed (Pfaffian-native z2 column {-1,+1,0})
- T4 closed BLOCKED (non-flat colormap assertion deferred to U13)
- T4r closed BUG REFUTED (slim CSR Pfaffian correctly constructed)
- T5 closed (gate precondition update z2==-1 native, full cleanup)
- T6 closed (per-B |Pf|_min saturation detection)

Map open frontier: empty. **Pick-B parallel track remains:** PR #43 review (T1b + T2 + T3). Pick-C lower priority: U11 (lecture revamp, depends on U10 + U12 gate formalization). U9 (BdG spectral + LDOS) and U1 (open BdG follow-ups) are separate scoped units. U13 (full Bloch-Pfaffian + periodic Peierls) is BLOCKING-EMPIRICAL, deferred per CLAUDE.md Known Issues.
