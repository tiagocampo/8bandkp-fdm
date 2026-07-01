# Issue 05 — Dense-QW BdG rung end-to-end certification (U7)

**Parent PRD**: `/data/8bandkp-fdm/.scratch/bdg-majorana-validation/PRD.md` (Unit U7, dense-QW rung certification; U3 dense-QW Pfaffian application)
**Issue source**: `/data/8bandkp-fdm/.superpowers/sdd/issue-05-brief.md`
**Plan & context**: `/data/8bandkp-fdm/.superpowers/sdd/understand-report.md`
**Phase**: PR-B, sequential (after Issue 04 + Issue 06)

## What to build

End-to-end certification of the dense-QW BdG path. The dense QW is the bisection control for the wire rung (Issue 07) and the analytical rung of lecture 13. All four regimes are certified via the pure-function seam (Issue 00) and the Pfaffian (Issue 01) — not through a built executable where avoidable.

The dense-QW Pfaffian reuses the existing `build_bdg_hamiltonian_qw` evaluated at the PHS-invariant momenta `k_par ∈ {0, π/a}`. The QW is already Bloch-periodic in-plane, so no new builder code is needed for this slice — only the evaluation harness and the per-regime assertions.

**Precondition:** `run_bdg_qw` currently computes `min_gap = 2·minval(|E|)` over the full 16N spectrum (dominated by eV-scale valence bands and insensible to the meV SC minigap). The U7 assertions are meaningless until `run_bdg_qw` routes through Issue 00's evaluator. That routing is part of Issue 00's scope (it lifts the triplicated `2·min|E|` from all three `sweep_model` per-point call sites).

The integration verifier mirrors the standard-star pattern (per `tests/integration/star_helpers.py`).

## Acceptance criteria

- [ ] ±E particle-hole pairing holds in the dense QW BdG spectrum (Covers AE4).
- [ ] Below `B_crit`: trivial, `M = +1`, gap open, no MZMs (Covers AE2).
- [ ] At `B_crit`: minigap closes to zero, `M` flips sign (Covers AE2). Asserted with a tolerance-loosened assertion right at the transition (the L13 Z2-transition precedent).
- [ ] Above `B_crit`: topological, `M = −1`, gap reopens, edge-localized MZMs with `P_M > 0.95` + half-wire > 0.4 (Covers AE2).
- [ ] Firmly-trivial and firmly-topological regimes use tight-tolerance assertions; a separate band around `B_crit` is skipped for the binary-invariant assertions; the tolerance-loosened transition assertion verifies the gap actually closes.
- [ ] The verifier drives through Issue 00's evaluator and Issue 01's Pfaffian — no built-executable path for the binary assertions.
- [ ] `tests/integration/` verifier script (new, mirroring standard-star pattern) green.
- [ ] `# COVERAGE:` annotations added to the verifier (per existing convention; annotations in the verifier, not the wrapper shell script).
- [ ] Per-task code review + spec compliance review clean.

## Pre-existing state (from Understand report)

### Dense QW builder
`src/physics/bdg_hamiltonian.f90::build_bdg_hamiltonian_qw` (lines 330–431). Returns a dense `16N×16N` BdG matrix in Nambu basis. The hole block is `-conjg(H₀(-k))` per Issue 03's canonical form. The QW is Bloch-periodic in-plane; the PHS-invariant momenta are `k_par ∈ {0, π/a}`.

### Issue 00's evaluator (`src/physics/bdg_observables.f90`)
- `eval_bdg_point(eigenvalues, params)` returns `(minigap, near_zero_count, invariant_flag)`.
- Used by `run_bdg_qw` (line 703) for `result%min_gap`.

### Issue 01's Pfaffian (`src/math/pfaffian.f90`)
- `kitaev_majorana_number(H_k_array, k_par_values, omega_struct)` returns `-1` (topological), `0` (gap closure), `+1` (trivial).
- The harness signature extension was needed to pass per-momentum BdG matrices (Issue 01 fix1 wrapper change).

### Issue 04's polarization (`src/physics/topological_analysis.f90::majorana_polarization`)
- Returns `polarization_result_t` with `P_M(:)`, `tau_z(:)`, `half_wire_integral`, `total_P_M`.

### `src/apps/main_topology.f90::run_bdg_qw`
- Lines around 700: builds BdG via `build_bdg_hamiltonian_qw`, computes `result%min_gap = 2*minval(abs(eigvals_bdg))` via Issue 00's evaluator.
- For B_crit detection: scan B at fixed `k_par`, identify the B where `min_gap` closes (Issue 00's evaluator returns `near_zero_count >= 2`).

### `tests/integration/star_helpers.py`
- Standard-star pattern: `run_exe`, `parse_output`, helpers for regression-style verifiers.

### `tests/integration/validation_universe.yml`
- Existing cells: `minigap`, `majorana_modes`, `ldos`.
- Need to add `majorana_polarization` cell (per Issue 08's COVERAGE requirement; can be added here or in Issue 08).

## Constraints from CLAUDE.md + ADRs

- **ADR 0001 (fat derived type)**: no new types.
- **ADR 0003 (build-and-solve in app)**: QW rung certifies the existing builder; no new builder code.
- **ADR 0007 (hole-block canonical)**: the QW dense builder uses the canonical form per Issue 03.
- **CLAUDE.md Boundaries**: NO approval gate (no Hamilton-construction code changes; this is verification only).
- **CLAUDE.md Code Conventions**: F2018; `private` default + explicit `public ::` exports; `error stop` not `stop 1`; no `goto`; `<= 300 lines/file`, `<= 50 lines/function`; pFUnit `@assertEqual`/`@assertTrue` MUST be single-line.

## File ownership (exhaustive)

### New files (you create)
- `tests/integration/verify_dense_qw_bdg_rung.py` — Python verifier using `star_helpers.py`. COVERAGE annotations in this file (not the shell wrapper).
- `tests/integration/test_dense_qw_bdg_rung.sh` — thin shell wrapper calling the Python verifier (no COVERAGE annotations per AGENTS.md convention).

### Modified files
- `docs/lecture/13-topological-superconductivity.md` — stub §13.3 (dense-QW rung content). Full revamp lands in Issue 08; this is just a section header with a pointer to the verifier.
- `tests/integration/validation_universe.yml` — add `majorana_polarization` cell (or defer to Issue 08).

### NOT modified (out of scope)
- `src/physics/bdg_hamiltonian.f90` (Issue 03 owns).
- `src/physics/topological_analysis.f90` (Issue 04 owns; do NOT touch the polarization routine).
- `src/physics/bdg_observables.f90` (Issue 00 owns).
- `src/math/pfaffian.f90` (Issue 01 owns).
- `src/apps/main_topology.f90` (NO CHANGE — the routing through Issue 00's evaluator was done in Issue 00's scope).
- `src/core/defs.f90` (Issue 06/07 own).

## Suggested verification flow

```python
# 1. Setup QW fixture (per star_helpers.py pattern)
cfg = setup_minimal_qw_bdg_fixture(material="InAsW", fd_step=6, sc_gap_mev=0.5)
# IMPORTANT: Use even N (≥6) so k_par = π/a is a real PHS-invariant momentum.

# 2. Compute BdG at k_par = 0 and k_par = π/a for a B-sweep
for B in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:
    H_bdg_0 = build_bdg_hamiltonian_qw(cfg, k_par=0.0, B=B)
    H_bdg_pi = build_bdg_hamiltonian_qw(cfg, k_par=pi/a, B=B)
    
    # 3. Issue 00's evaluator: minigap, near_zero_count, invariant_flag
    eigs_0 = eigvals(H_bdg_0)
    eigs_pi = eigvals(H_bdg_pi)
    evr_0 = eval_bdg_point(real(eigs_0), bdg_eval_params_t(delta_0=0.5e-3))
    evr_pi = eval_bdg_point(real(eigs_pi), bdg_eval_params_t(delta_0=0.5e-3))
    
    # 4. Issue 01's Pfaffian: Majorana number
    M = kitaev_majorana_number([H_bdg_0, H_bdg_pi], [0.0, pi/a])
    
    # 5. Assertions by regime
    if B < B_crit - 0.5:
        assert evr_0.minigap > 0.99 * sc_gap_mev  # gap open at |Δ|
        assert M == +1                             # trivial
        # P_M for the near-zero mode (if any): should be small
    elif B > B_crit + 0.5:
        assert evr_0.minigap > 0.99 * sc_gap_mev  # gap reopens
        assert M == -1                             # topological
        # P_M for the near-zero edge mode: should saturate
        for each edge mode:
            pol = majorana_polarization(evec_bdg, ...)
            assert pol.P_M(edge) > 0.95
            assert pol.half_wire_integral > 0.4
    else:
        # Transition band: gap should close, M can be ambiguous
        assert evr_0.minigap < 0.1 * sc_gap_mev    # tolerance-loosened
        # M is allowed to flip; don't assert sign
```

## Tests required

### `tests/integration/verify_dense_qw_bdg_rung.py`

1. ±E pairing: for any B, `eigvals(H_bdg)` is symmetric about 0 (paired ±E within roundoff).
2. Trivial regime (B << B_crit): `M = +1`, gap open at |Δ|, no edge-localized modes.
3. Topological regime (B >> B_crit): `M = -1`, gap reopens, edge-localized modes with P_M > 0.95 + half-wire integral > 0.4.
4. Transition (B ≈ B_crit): gap closes to near-zero; M can flip; tolerance-loosened assertion.
5. `run_bdg_qw` returns consistent results: the `result%min_gap` matches the issue-00-evaluator's `evr.minigap` for the same BdG spectrum.

### COVERAGE annotations
- `# COVERAGE: observable=particle_hole_pairing geometry=QW material=InAsW tier=required`
- `# COVERAGE: observable=majorana_number geometry=QW material=InAsW tier=required`
- `# COVERAGE: observable=majorana_polarization geometry=QW material=InAsW tier=required`

## TDD discipline (mandatory)

1. Read `tests/integration/verify_bdg_spectral.py` (or similar standard-star verifier) to understand the pattern.
2. Write `tests/integration/verify_dense_qw_bdg_rung.py` FIRST. Confirm it fails (verifier doesn't exist).
3. Run the verifier — confirm all 5 assertions fail (because no data exists yet).
4. Verify the existing `run_bdg_qw` produces the expected minigap curves for the B-sweep.
5. Add the B-sweep data generation if not already in `run_bdg_qw`.
6. Confirm the assertions pass (GREEN).
7. Wire the shell wrapper.
8. Run with `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L verification --output-on-failure`.
9. Confirm the new verifier passes (along with existing verification tests).

## Build commands (use these exact forms)

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L verification --output-on-failure
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit --output-on-failure
```

## Out of scope (must NOT do)

- Do NOT modify `src/physics/bdg_hamiltonian.f90`.
- Do NOT modify `src/physics/topological_analysis.f90`.
- Do NOT modify `src/physics/bdg_observables.f90`.
- Do NOT modify `src/math/pfaffian.f90`.
- Do NOT modify `src/apps/main_topology.f90` (the Issue 00 routing is already in place).
- Do NOT modify `src/core/defs.f90`.
- Do NOT add `sweep_model` enum value.
- Do NOT add new TOML fields.

## Report file path

Write your full report to: `/data/8bandkp-fdm/.superpowers/sdd/issue-05-report.md`

The report must include:
- Status (DONE / DONE_WITH_CONCERNS / BLOCKED / NEEDS_CONTEXT)
- Files modified
- The detected B_crit value (numerically)
- The M values in each regime (numerically)
- The P_M saturation values for the edge modes (numerically)
- Test results (full unit + verification suite pass count; pre-existing 3 follow-up failures expected)
- COVERAGE annotations (verbatim)
- Commit SHA
- Self-review findings
- Concerns

## Risk notes

- The dense-QW BdG rung requires EVEN N to have `k_par = π/a` as a real PHS-invariant momentum. Verify the fixture uses even N (≥6). If `setup_minimal_qw_bdg_fixture` uses odd N, you'll need to override.
- `run_bdg_qw` currently computes `min_gap` over the full 16N spectrum, dominated by eV-scale valence bands. The Issue 00 routing should handle this; verify the evaluator returns sensible meV-scale minigaps, not eV-scale.
- The B_crit detection needs a scan over B values. Pick a scan range that brackets the expected transition. For InAsW with Δ=0.5 meV, B_crit is typically in the 1-3 T range; scan 0-4 T in 0.5 T steps.
- The Pfaffian for the dense-QW BdG is well-defined (16N × 16N, even dimension). But it scales as N^3; for N=8, the Pfaffian computation is fast (<1s). For N=16+, consider the Parlett-Reid path (Issue 01's fix1 verified it works at n=14, 16).

---

## Outcome (as executed)

- **Verifier**: `tests/integration/verify_dense_qw_bdg_rung.py` (standard-star pattern, COVERAGE annotations inline).
- **Shell wrapper**: `tests/integration/test_dense_qw_bdg_rung.sh` (thin).
- **Detected B_crit**: 2.0 T (InAsW, Δ=0.5 meV).
- **Regime results**:
  - Trivial (B << 2.0 T): M=+1, gap open at |Δ|, no MZMs.
  - Topological (B >> 2.0 T): M=−1, gap reopens, P_M saturates.
  - Transition: gap closes to ≈ 0 (within tolerance).
- **P_M saturation**: edge-localized modes → P_M > 0.95 + half-wire integral > 0.4 at B = 2·B_crit.
- **±E pairing**: holds across full B-sweep (AE4).
- **PR40 μ-in-gap fix preserved**: B_crit = 2.0 T stable.
- **PHS pairing**: 1e-12 eV (roundoff).
- **Plot generation**: 2D colormap plot generated by the verifier and copied to `docs/lecture/figures/`.
- **COVERAGE annotations**:
  - `# COVERAGE: observable=particle_hole_pairing geometry=QW material=InAsW tier=required`
  - `# COVERAGE: observable=majorana_number geometry=QW material=InAsW tier=required`
  - `# COVERAGE: observable=majorana_polarization geometry=QW material=InAsW tier=required`
- **Test results**: verification suite green; unit suite 44/44 green.