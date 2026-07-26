**Status**: COMPLETE (2026-07-13) — ticket 05 of `.scratch/archive/bdg-evaluator-pfaffian/`.

Type: grilling
Status: closed (2026-07-13)
Blocked by: 01, 02, 04
Resolved-by: grilling 5 questions, all answered by user

# 05 — Gate slim-witness wiring

## Question

The 3-witness acceptance gate at `tests/integration/test_lecture_13_acceptance_gate.sh` has a wire Pfaffian row currently `@todo U13` per Phase 24 follow-up ticket A. With (d+ii), that row becomes the slim projected Pfaffian witness as the documented witness for the wire rung.

Decide:

- **The witness criterion.** What threshold + tolerance + observable does the slim witness row assert? (E.g., "slim Pfaffian sign matches expected sign at the test B-field within `1.0 T` of the analytical `B_crit`" — or whatever shape the existing two gate rows use, since the gate already has a `1.0 T` `TOLERANCE_BCRIT_RANGE` per Phase 24 ticket A.)
- **The witness value source.** Is the slim witness row reading from the gate's existing `bcrit_pfaffian` parser, or from a new artifact produced by the slim-witness call site in `main_topology.f90`?
- **The `@todo U13` swap.** When (d+ii) lands, does the slim witness row replace the `@todo` outright, or stay `@todo` for the future Bloch-Pfaffian upgrade? (Map's "Not yet specified" fog lists this as a sub-decision.)
- **The lecture 13 disclosure.** What does the lecture markdown (`docs/lecture/13-topological-superconductivity.md`) and the companion script (`scripts/lecture_13_topological.py`) say about the slim witness? The 4→3 witness label sweep (Phase 24 ticket A) becomes 4-witness again — slim witness as the 4th row — and the disclosure text needs to reflect that.
- **The `verify_majorana_polarization.py` BLOCKING-EMPIRICAL SKIP**. Out of scope per the map's Out-of-scope, but the gate's wire Pfaffian row interacts with the polarization verifier's SKIP path — confirm no breakage.

Deliverable: the gate-row spec (criterion + value source + tolerance) and the lecture/script disclosure text. Both land in the resolution.

## Constraints

- The gate is a shell script — bash, not Python. The criterion must be expressible as a shell exit-code check.
- The existing two rows (likely bcrit_curve + bcrit_2d) are the precedent — match their style.
- Lecture disclosure must match the actual witness (no aspirational text). Phase 24 follow-ups established a discipline: the lecture status table is regenerated from script output, not hand-edited.
- The slim witness's noise-floor tolerance is owned by ticket 02; ticket 05's criterion depends on 02's output.

## Out of scope for this ticket

- The slim witness's noise-floor tolerance itself (ticket 02).
- The migration of the call sites that produce the witness value (ticket 04).
- The plan + backlog close-out (ticket 07).

## Answer

### Gate-row spec (the deliverable)

- **Criterion.** `bcrit_pfaffian = B at mu ≈ 0.6601 (within ±0.0001) where the colormap's `z2` column first reads -1`. The `z2 == -1` row is "slim Pf returned 0, gap-closure fallback fired" per ticket 04 site 6 (`compute_z2_gap_bhz_heuristic` at `main_topology.f90:1380`). This is **the topological invariant itself**, not a gap proxy.
- **Value source.** `output/z2_phase_diagram.dat` — produced by `compute_wire_bdg_gap_sweep` and already on disk. **No new Fortran I/O.** Parser location: inline in `scripts/lecture_13_topological.py:section_wire_rung()`, ~15 LOC new (mu filter at ±0.0001, z2 column read, first z2==-1 row lookup). Dead references to `output/wire_slim_pfaffian_witness.dat` and the silent alias `bcrit_pfaffian = bcrit_curve` are removed.
- **Tolerance.** `TOL_BCRIT_RANGE = 1.0 T` (gate + lecture script) — unchanged. Covers all four witnesses (canonical InAs wire B_crit ≈ 2.5–2.8 T).
- **Gate label.** 4-witness always; `WITNESS_LABEL = "4-witness"`; `PFAFFIAN_DEFERRED` branch stripped.

### Five grilling decisions

1. **`@todo U13` swap — (b) live derivation from colormap.** Re-derive `bcrit_pfaffian` from `output/z2_phase_diagram.dat` z2 column instead of waiting for U13 to emit a new witness file. Zero new Fortran I/O; colormap already encodes the witness from the migrated `eval_bdg_pfaffian_witness_csr` call at `main_topology.f90:1371`. Trade vs. (a) live writer: rejects untested I/O at a stabilization moment. Trade vs. (c) defer: closes U2 visually as 4-witness instead of 3.5.

2. **Criterion shape — (i) closure-row B.** `bcrit_pfaffian = arg_B at mu≈0.6601 where z2 == -1` (first such B). Disambiguates from `bcrit_2d` (which is `argmin gap`, a different observable); preserves the 4-witness design's distinct measurements.

3. **Lecture disclosure — (b) textual flip + new §13.7.5 paragraph.** Sweep every "deferred to U13" / "reserved for U13" / "synthetic illustrative" reference on the slim Pfaffian row to "live (colormap-extracted at z2==-1, mu≈0.6601)." Add §13.7.5 "Slim Pfaffian witness criterion" paragraph spelling out the criterion (`arg_B at mu≈0.6601 where z2 == -1`), the rationale (z2=-1 is the slim Pf returned 0 by construction), and a forward reference to U13 as the future Bloch-Pfaffian swap.

4. **Gate housekeeping — (a) strip dead code.** Strip `PFAFFIAN_DEFERRED` branch in the gate (`:79–87`), always compute 4-witness range. Remove `bcrit_pfaffian = bcrit_curve` alias. Update `TOL_BCRIT_RANGE` comment. Keep Plot 8 (Issue 11) as synthetic illustrative — it's the witness *concept*, U13 replaces it with a real Bloch-Pfaffian plot.

5. **Polarization SKIP — (b) tighten precondition.** Update docstring (~10 LOC) clarifying SKIP is on the polarization observable (Sticlet, n_majorana-gated), distinct from the gate's invariant row. **Plus** add a precondition in the SKIP branch: parse `output/z2_phase_diagram.dat` at mu≈0.6601 and verify z2==-1 row exists. If not, exit non-zero with "slim Pfaffian row would fail" — converts SKIP from "ignore this until U13" to "skip polarization assertion only if gate row is green."

### Files touched (downstream execution, not done here)

The decisions above specify the design; execution lands in a follow-up commit/PR (ticket 05 is a Grilling ticket — wayfinder produces decisions, not deliverables):

- `tests/integration/test_lecture_13_acceptance_gate.sh` — strip `:79–87` (`PFAFFIAN_DEFERRED` branch), update `:101` `TOL_BCRIT_RANGE` comment, hard-code `WITNESS_LABEL="4-witness"` at `:114`.
- `scripts/lecture_13_topological.py` — replace parser block at `:222–238` with colormap z2-column reader (`mu` filter ±0.0001, first z2==-1 row), remove `:230` alias, update `:289` `TOLERANCE_BCRIT_RANGE` comment.
- `scripts/lecture_13_topological.py:305, 319–323, 326–327, 350–356` — table cell for Pfaffian row flips to live numeric; `WITNESS_LABEL` summary row 4-witness.
- `docs/lecture/13-topological-superconductivity.md` — sweep lines 14, 28–33, 305, 528–533, 622–650, 836, 957, 965, 976; drop "deferred to U13" / "synthetic illustrative"; add §13.7.5 paragraph (criterion + rationale + U13 forward reference).
- `tests/integration/verify_majorana_polarization.py` — update module docstring (lines 4–28) + inline comment (lines 102–108); add precondition: in the SKIP branch, parse `output/z2_phase_diagram.dat` at mu≈0.6601, verify z2==-1 row exists, else exit non-zero.

Net: 4 files touched, 0 Fortran I/O added. The slim witness becomes a **real, colormap-derived 4th witness of the gate**, with criterion spec locked in this ticket.