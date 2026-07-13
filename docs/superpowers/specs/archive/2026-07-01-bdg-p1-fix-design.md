# BdG/Majorana P1 stabilization — design spec

**Status**: IMPLEMENTED (27 commits on `feat/bdg-p1-stabilization`)
**Date**: 2026-07-01
**Branch**: `feat/bdg-p1-stabilization` (cut from `feat/bdg-validation-pass2`)
**Source PR**: #41 (BdG/Majorana ground-up validation pass, 26 commits, +11k/-22.7k)
**Source review**: 6-dimension adversarial review (physics, test-quality, architecture, build-hygiene, doc-drift, edge-cases), 18 P1 / 35 P2 / 19 P3 findings

**Implementation summary** (2026-07-02):
- 27 commits ahead of `feat/bdg-validation-pass2`
- All 49 unit tests + 9 BdG regression tests + acceptance gate PASS (1.75 T ≤ 2.0 T tolerance)
- Verification: `bash tests/integration/test_lecture_13_acceptance_gate.sh`
- ADR 0008 §3 "double-count Peierls" claim was INVESTIGATED and disproved — symmetric `add_peierls_coo(-B)` is REQUIRED for class-D PHS, restored with documentation. Dense-QW block (2,1) had a structural row/col swap bug beyond the sign bug.

See plan at `docs/superpowers/plans/2026-07-01-bdg-p1-fix.md` and ADR at `docs/adr/0008-bdg-p1-invariants.md`.

## 1. Locked decisions (from brainstorming)

| Decision | Choice | Rationale |
|---|---|---|
| Scope | **C — full stabilization** (P1 + P2 + P3) | Multi-week but comprehensive; matches user's chosen scope |
| Kitaev Majorana number fix | **Lutchyn-Oreg sign-of-det** (bypass Pf) | Literature-canonical, handles diagonal-in-(c,c†) BdG, makes strict tests possible |
| Symmetric Peierls philosophy | **A — keep symmetric (re-investigated)** — symmetric `add_peierls_coo(-B)` on hole block is REQUIRED for class-D PHS at generic k with Bx≠0. The `-conjg()` transform alone is insufficient. See ADR 0008 §3 amendment in `2026-07-05-pr41-completion-design.md` §6.3. |
| 2D colormap regression test | **B — replace with minigap-pattern assertion** | Robust physics assertion; avoids brittle z2=±1 numerical sign |
| Acceptance gate tolerance | **A — add absolute window guard** | Catches uniform regression that pure range check misses |
| Module split | **C — defer to follow-up** | Pre-existing 1437-line bloat is separate from P1 stabilization |
| Branching | **A — one PR with phased commits** | Matches PR #41 cadence; atomic merge; bisect-friendly |
| Half-wire integral | **A — `max(left, right)`** | Matches existing test fixture intent |
| Approach merge | **TDD-Disciplined, ADR-Designed, 4-Phase PR** (A+B+C combined) | TDD on critical path (Phase 1), ADR upfront (Phase 0), 4-phase structure (atomic merge) |

## 2. Approach: TDD-Disciplined, ADR-Designed, 4-Phase PR

**Layer 1 (from C — ADR-first design)**:
- **ADR 0008** locks the four big design decisions: Lutchyn-Oreg sign-of-det, pairing_sign DRY public API, symmetric `add_peierls_coo(-B)` on hole block KEPT (re-investigated 2026-07-05 — `conjg()` alone is insufficient for class-D PHS), replace 2D golden with minigap-pattern test.
- **This spec** captures the full design.
- **Implementation plan** at `docs/superpowers/plans/2026-07-01-bdg-p1-fix.md` enumerates ~40 commits.

**Layer 2 (from A — 4-phase structure)**:
- Phase 0 = Design artifacts (1 commit, docs only) — ADR 0008 + spec + plan
- Phase 1 = Physics correctness (13 commits, 4 TDD triples + 1 prerequisite + 2 single) — merge-blocker
- Phase 2 = Test teeth (6 commits)
- Phase 3 = Architecture cleanup (6 commits)
- Phase 4 = Hygiene + docs (8 commits, with 4.2 dependent on 3.2)

**Total**: ~34 commits including the prerequisite refactor (1.1).

**Layer 3 (from B — TDD discipline within each phase)**:
- For each P1 fix in Phase 1: failing test commit → fix commit → regen/refactor commit (3:1 ratio on critical path)
- For P2/P3 fixes: fix + test in one commit (single-commit refactors)
- All pFUnit tests added in this PR follow test-first pattern

## 3. Architecture overview

One new ADR (`0008-bdg-p1-invariants.md`), one spec (this file), one PR (`feat/bdg-p1-stabilization`), four phases. Each phase is a coherent unit that compiles + tests green at its tip. No cross-phase dependencies.

### Component map

| Phase | Module(s) touched | New module(s) | New ADR |
|---|---|---|---|
| 0 | (docs only) | — | `0008-bdg-p1-invariants.md` |
| 1 | `bdg_hamiltonian.f90`, `pfaffian.f90`, `topological_analysis.f90`, `main_topology.f90`, `magnetic_field.f90` | — | (rolled into 0008) |
| 2 | `tests/integration/*`, `tests/regression/*` | — | — |
| 3 | `eigensolver.f90`, `sparse_matrices.f90`, `outputFunctions.f90`, `bdg_observables.f90` | — | — |
| 4 | `topological_analysis.f90`, `spectral_bdg_wire.f90`, `validation_universe.yml`, `lecture.md`, PRD, AGENTS.md | — | — |

**Phase ordering note**: Phase 4.2 (pointer caching on `csr_to_dense_work`) depends on Phase 3.2 (move of `csr_to_dense_work` to `sparse_matrices.f90`). Phase 4 must execute after Phase 3 for this commit. All other Phase 4 commits are independent.

### Key architectural commitments (ADR 0008)

1. **Lutchyn-Oreg sign-of-det** replaces Kitaev's `real(Pf)` for `kitaev_majorana_number`. Bypasses Pf sign extraction entirely.
2. **`pairing_sign` is PUBLIC** in `bdg_hamiltonian.f90`. `topological_analysis.f90` imports via `use bdg_hamiltonian, only: pairing_sign, pairing_partner`.
3. **Symmetric `add_peierls_coo(-B)` on hole block is KEPT** (re-investigated 2026-07-05, reversal of the original Phase 1.10 "remove symmetric Peierls" task — see ADR 0008 §3 amendment). Both wire and dense-QW builders apply the symmetric Peierls call; the `-conjg()` transform in `build_bdg_hole_block` alone is insufficient for class-D PHS at generic k with Bx≠0.
4. **Minigap-pattern regression** replaces z2-golden test. Acceptance gate adds absolute window guard.

## 4. Components in detail

### Phase 0 — Design artifacts (1 commit, docs only)

**0.1** Write ADR 0008 + this spec + implementation plan. Commit docs only, no code change. Branched from `feat/bdg-validation-pass2`.

### Phase 1 — Physics correctness (12 commits, TDD-doubled critical path)

**1.1** [Prerequisite refactor] Verify that `add_peierls_coo` in `magnetic_field.f90:184` uses `exp(-iφ)` convention (`exp_phase = cmplx(cos(phase), -sin(phase))` with `phase = e*Bx*(y_i-y_j)/hbar`). Confirmed by code inspection 2026-07-01. Document finding in commit message so the removal in 1.10 is traceable.

**1.2** [TDD red] Add `tests/unit/test_kitaev_strict.pf` asserting `M=-1` for |μ|<2t, `M=+1` for |μ|>2t on a 2-band Lutchyn-Oreg fixture.

**1.2 [TDD red]** Add `tests/unit/test_kitaev_strict.pf` asserting `M=-1` for |μ|<2t, `M=+1` for |μ|>2t on a 2-band Lutchyn-Oreg fixture. Test fails because `kitaev_majorana_number` returns 0.

**1.3 [TDD green]** Replace `real(pf_val)` → `atan2(aimag, real)` phase extraction in `pfaffian.f90:161`. Test still fails.

**1.4 [TDD refactor]** Implement full Lutchyn-Oreg formula: `M = sign(det(U_odd))` where `U = sign(H)` from eigendecomposition of H_bdg restricted to odd-parity Nambu subspace. Add helper `polar_decomposition_sign` to `pfaffian.f90`. `kitaev_majorana_number` delegates to Lutchyn-Oreg. Test passes.

**1.5 [TDD red]** Add `tests/unit/test_pairing_sign_cross_module.pf` asserting `bdg_hamiltonian.f90::pairing_sign == topological_analysis.f90::s_sigma` byte-identical. Test fails because the latter is `private` local copy.

**1.6 [TDD green]** Add `public :: pairing_sign, pairing_partner` to `bdg_hamiltonian.f90`. Remove local copy in `topological_analysis.f90`. Add `use bdg_hamiltonian, only: pairing_sign, pairing_partner`. Test passes.

**1.7 [TDD red]** Add `tests/unit/test_dense_qw_block_21_sign.pf` asserting Hermiticity of `build_bdg_hamiltonian_qw` output: `|H(i,j) - conjg(H(j,i))| < 1e-12` for all i,j. Test fails at line 547.

**1.8 [TDD green]** Fix line 547: `pairing_sign(pairing_partner(ib))` → `pairing_sign(ib)`. Test passes. Verify PHS oracle still green.

**1.9 [TDD red]** Add `tests/unit/test_phs_qw_peierls_symmetric.pf` asserting PHS at generic k with `Bx≠0` for dense-QW builder. Test currently passes vacuously (no Peierls applied).

**1.10 [TDD green — REVERTED 2026-07-05]** Originally: "Remove the second `add_peierls_coo` call in `bdg_hamiltonian.f90:419-428` (the symmetric-Peierls fix3 over-correction). `conjg()` in `build_bdg_hole_block` is sufficient. Verify wire PHS oracle still green." — REVERTED: empirical PHS oracle verification (rel_resid 1.25e-1 → 0 on removing the call) demonstrated the symmetric call IS REQUIRED for class-D PHS at generic k with Bx≠0. The call was RESTORED; investigation note lives at `src/physics/bdg_hamiltonian.f90:406-419`. See ADR 0008 §3 amendment and §6.3 below.

**1.12** Fix `half_wire_integral = max(sum(1..N/2), sum(N/2+1..N))` in `topological_analysis.f90:939`. Existing test `test_mzm_peaks_at_wire_ends` exercises both ends; update assertion to match new convention.

**1.13** Replace `real(u·v*)` with `u·v*` (complex) at `topological_analysis.f90:919`. Add `test_majorana_polarization_soc.pf` for SOC-active fixture (synthetic BdG with imaginary coherence).

### Phase 2 — Test teeth restoration (6 commits)

**2.1** Replace z2 column in `tests/regression/data/wire_inas_gaas_bdg_topological_2d_phase.dat` with minigap-only data. Update `test_wire_bdg_topological_2d.sh` to assert minigap pattern: small near B_crit, large far from it.

**2.2** Wire `bcrit_pfaffian` in `scripts/lecture_13_topological.py:196` to real Pfaffian computation: invoke `topologicalAnalysis` with `slimp faffian` mode over a B-grid, parse output, emit `B_crit = argmin |Pf(kz, B)|`. Drop the `bcrit_curve` alias.

**2.3** Wire `bcrit_2d` in `scripts/lecture_13_topological.py:172-189` to read `output/z2_phase_diagram.dat` (live Fortran output) instead of golden.

**2.4** Add absolute window guard to `tests/integration/test_lecture_13_acceptance_gate.sh:87`: `assert 0.5 <= BCRIT_MIN and BCRIT_MAX <= 6.0` for InAs wire. Keeps range check at 2.0 T.

**2.5** Rename `tests/integration/verify_majorana_polarization.py` to `tests/integration/plot_majorana_polarization_synthetic.py`. Move from `verification` to new `plot` ctest label.

**2.6** Replace literal-substring grep in `tests/integration/test_lecture_13_acceptance_gate.sh:104` with regex anchor.

### Phase 3 — Architecture cleanup (6 commits, single-commit fixes)

**3.1** Promote `asw_evals/asw_single/asw_envelope` to module-private in `src/math/eigensolver.f90:17`. Update `green_functions.f90:283` to use `apply_solver_window` generic.

**3.2** Move `csr_to_dense_work` from `spectral_bdg_wire.f90:240-253` and `extract_block_csr` from `spectral_bdg_wire.f90:197-234` to `src/math/sparse_matrices.f90`.

**3.3** Move 6 inline `open()/write()` blocks from `main_topology.f90` to `outputFunctions.f90` as named writers: `write_majorana_profile`, `write_bdg_lowest_state_profile`, `write_topology_result`, `write_spectral_function`, `write_z2_phase_diagram`, `write_z2_transitions`.

**3.4** Add `bdg_default_near_zero_frac = 0.001_dp` module parameter and `bdg_eval_params_with_delta(delta_0)` factory to `src/physics/bdg_observables.f90`. Collapse 5 call sites.

**3.5** Update AGENTS.md rows: `bdg_hamiltonian.f90` mentions `build_bdg_hole_block` wrapper; `pfaffian.f90` listed; `lecture_13_topological.py` 5-section structure.

**3.6** Update `CLAUDE.md` Known Issues with U13 wire periodic/Bloch BdG + bdq_spectral wire-only deferred items.

### Phase 4 — Hygiene + docs (8 commits, single-commit fixes)

**4.1** Replace 4 bare `stop 1` in `topological_analysis.f90:964, 1132, 1496, 1501` with `print + error stop` (pattern in §6.3).

**4.2** Add pointer caching to `csr_to_dense_work` in `spectral_bdg_wire.f90:240-253` (the new helper introduced in this PR). Pattern matches `hamiltonian_wire.f90`.

**4.3** Add `# COVERAGE: observable=<name> geometry=<g> material=<m> tier=<t>` annotations to all 7 new pFUnit files and 5 modified files.

**4.4** Add `bdg_ldos_topological`, `bdg_spectral_function_topological`, `bdg_nambu_ldos_topological` cells to `tests/integration/validation_universe.yml`. Add the new observables to metadata list.

**4.5** Delete hand-edited reconciliation table at `docs/lecture/13-topological-superconductivity.md:38-43`. Replace with: "See embedded image above for the four B_crit values; the markdown intentionally contains no hand-edited numbers." Add grep guard to acceptance gate.

**4.6** Reconcile LDOS zero-bias peak claim at `lecture.md:480` with Issue 11 fix1 figure caption.

**4.7** Update `PRD.md`: commit count 14→26, **create** the 3 missing solutions docs from existing materials (`issue-02-report.md`, `issue-03-report.md`, the `green_functions.f90` diff in commit `328b83a`). Per MEMORY.md `codebase-doc-drift-prevention` rule, deferring is a doc-drift source — create the docs. Then `git mv .scratch/bdg-majorana-validation .scratch/archive/bdg-majorana-validation` per MEMORY.md step 3.

**4.8** `git rm` stale all-zero `docs/lecture/figures/rashba_majorana_phase_diagram.txt` + companion `.png`. Update MEMORY.md to add 3rd drift event.

## 5. Data flow

### 5.1 BdG invariant computation pipeline (after Phase 1.4)

Current broken path:
```
H_bdg(k) → assemble_skew_product → A_work (skew-symmetric, but of anti-Hermitian B)
  → complex_pfaffian → pf_val (purely imaginary)
  → real(pf_val) → 0 (silently drops sign)
  → prod = 0 → M = 0 (gap closure for everything)
```

Fixed path (Lutchyn-Oreg):
```
H_bdg(k) at k=0 and k=π → sign(H_bdg) via eigendecomposition (Hermitian → real eigvals)
  → project sign(H_bdg) to odd-parity Nambu subspace (rows/cols 1..n_sp/2)
  → det of odd-parity block
  → sign(det) ∈ {-1, +1}
  → product over k ∈ {0, π}
  → M = sign(product) ∈ {-1, +1}, or M=0 if any det is zero (gap closure)
```

**Boundary**: eigendecomposition helper `polar_decomposition_sign(H_bdg) → U_odd` lives in `src/math/pfaffian.f90`. Public API of `kitaev_majorana_number` unchanged; only internal implementation changes.

### 5.2 pairing_sign data flow (after Phase 1.6)

Before (DRY violation):
```
bdg_hamiltonian.f90 (SSOT)
  └ pairing_sign(8) = [+1,+1,-1,-1,+1,-1,+1,-1]
       ↓ exported via public module parameter

topological_analysis.f90 (local duplicate)
  └ s_sigma(8) = [+1,+1,-1,-1,+1,-1,+1,-1]   ← DRIFT RISK
       └ majorana_polarization uses s_sigma
```

After (single source):
```
bdg_hamiltonian.f90 (SSOT)
  └ pairing_sign(8), pairing_partner(8) — PUBLIC
       ├ bdg_hamiltonian internals
       ├ topological_analysis.f90 imports via `use`
       └ unit test asserts byte-identity
```

### 5.3 4-witness acceptance gate data flow (after Phase 2)

Before (3 of 4 witnesses are constants):
```
Lecture script prints BCRIT <name> <value>
  ├ wire_curve      ← live from verify_wire_bdg_topological.py (real)
  ├ wire_2d         ← reads golden.dat (CONSTANT, z2=0 trap)
  ├ wire_pfaffian   ← aliased to wire_curve (TAUTOLOGY)
  └ qw_dense        ← live from verify_dense_qw_bdg_rung.py (real)

Acceptance gate: range check only (≤ 2.0 T) → catches nothing if all regress together
```

After (4 independent live witnesses + absolute window guard):
```
Lecture script prints BCRIT <name> <value>
  ├ wire_curve      ← live (unchanged)
  ├ wire_2d         ← live from output/z2_phase_diagram.dat (real)
  ├ wire_pfaffian   ← live from topologicalAnalysis slimpfaffian mode (real)
  └ qw_dense        ← live (unchanged)

Acceptance gate:
  ├ range check (≤ 2.0 T) — catches witness disagreement
  ├ absolute window (0.5 ≤ MIN, MAX ≤ 6.0) — catches uniform regression
  └ regex anchor — catches legacy Gershgorin false-PASS values
```

### 5.4 Peierls phase path through wire vs dense-QW builders (after Phase 1.10)

Before (asymmetric, possibly double-counted):
```
wire CSR builder:
  H_e = ZB8bandGeneralized(kz, +B)            → Peierls phase exp(-iφ) on electron block
  H_h = ZB8bandGeneralized(-kz, +B)           → Peierls phase exp(-iφ) on hole block
  hole = -conjg(H_h)                          → conjg flips phase to exp(+iφ)
  + add_peierls_coo(hole, -B)                 → adds exp(+iφ) AGAIN
  = hole carries exp(+2iφ)                    → wrong by 2x

dense QW builder:
  H_e = ZB8bandQW(k, +B)
  H_h = ZB8bandQW(-k, +B)
  hole = -conjg(H_h)                          → conjg flips phase to exp(+iφ)
  (no extra Peierls call)
  = hole carries exp(+iφ)                     → correct
```

After (symmetric, both correct):
```
wire CSR builder:
  H_e = ZB8bandGeneralized(kz, +B)            → Peierls phase exp(-iφ) on electron block
  H_h = ZB8bandGeneralized(-kz, +B)           → Peierls phase exp(-iφ) on hole block
  hole = -conjg(H_h)                          → conjg flips phase to exp(+iφ)
  (no extra Peierls call)
  = hole carries exp(+iφ)                     → correct, matches dense QW

dense QW builder:
  unchanged                                    → correct
```

**Peierls phase sign convention** (per ADR 0008 §3 amendment, verified in `magnetic_field.f90:184`): `add_peierls_coo` applies `exp(-iφ)` where `φ = e*Bx*(y_i-y_j)/hbar`. For class-D PHS C·H(k,B)·C⁻¹ = -H(-k,-B), the hole block must carry Peierls(-B) which yields exp(+iφ) (sign-flipped). The `-conjg()` in `build_bdg_hole_block` alone is insufficient (REVERTED 2026-07-05): BOTH `-conjg()` AND the symmetric `add_peierls_coo(-B)` call are required to satisfy class-D PHS at generic k with Bx≠0. See ADR 0008 §3 amendment and `src/physics/bdg_hamiltonian.f90:406-419` investigation note.

## 6. Error handling

### 6.1 Sentinel vs gap closure (Phase 1.4)

The current `kitaev_majorana_number` confuses two distinct failure modes:
- **Gap closure**: det(U_odd) = 0 at some k → M = 0 (legitimate topological boundary)
- **Failed decomposition**: NaN/Inf in eigvals of H_bdg → should error-stop, not silently return 0

After Phase 1.4:
```fortran
function kitaev_majorana_number(H_k_array, k_par_values, omega_struct) result(M)
  ! ... argument validation ...
  do i = 1, n_k
    call polar_decomposition_sign(H_k_array(:,:,i), U_odd, info)
    if (info /= 0) error stop 'kitaev_majorana_number: eigendecomposition failed'
    if (any_eigenvalue_imaginary_part_nontrivial) error stop 'kitaev_majorana_number: H_bdg not Hermitian'
    det_U_odd = compute_det(U_odd)
    if (abs(det_U_odd) < tol) then
      M = 0; return  ! legitimate gap closure
    end if
    prod = prod * sign(det_U_odd)
  end do
  M = sign(prod)
end function
```

### 6.2 bdg_eval_params_t near-zero threshold consistency (Phase 3.4)

Before: `eval_bdg_point` uses `frac * δ₀`; `q_zero_tol` uses `max(1e-10, frac*δ₀)`. Inconsistent for small δ₀.

After:
- `bdg_default_near_zero_frac = 0.001_dp` module parameter
- `bdg_default_min_threshold = 1.0e-10_dp` module parameter
- Both `eval_bdg_point` and `q_zero_tol` use the same composed threshold via private helper
- `bdg_eval_params_with_delta(delta_0)` factory pre-fills both

### 6.3 Bare `stop 1` → `error stop` (Phase 4.1)

4 sites in `topological_analysis.f90`. `error stop` in F2018 accepts a literal string or a character expression; some compilers reject runtime concatenation in the stop expression itself, so we use a `print` + `error stop` pattern (the print appears on stderr before the runtime halt):

- Line 964: `print *, 'build_bhz_wire_hamiltonian: invalid parameters (N=', N, ', dz=', dz); error stop 'build_bhz_wire_hamiltonian: invalid parameters'`
- Line 1132: `print *, 'compute_phase_diagram: nB=', nB, ' nMu=', nMu, ' out of [1, 1000]'; error stop 'compute_phase_diagram: nB and nMu must be in [1, 1000]'`
- Line 1496: `print *, 'compute_z2_gap_sweep: unsupported sweep_model ', trim(cfg%topo%sweep_model); error stop 'compute_z2_gap_sweep: unsupported sweep_model (expected bhz_analytic)'`
- Line 1501: `print *, 'compute_z2_gap_sweep: evaluator failed for model ', trim(cfg%topo%sweep_model); error stop 'compute_z2_gap_sweep: evaluator failed'`

### 6.4 New failure modes introduced

| New failure mode | Detection | Response |
|---|---|---|
| Non-Hermitian H_bdg passed to `kitaev_majorana_number` | `any(imag(eigvals) > 1e-12)` | `error stop` |
| Eigenvalue imaginary part exceeds 1e-10 (numerical noise vs structural) | Unit test on real fixtures | `error stop` |
| `pairing_sign` import mismatch (e.g., new module redefines) | Compile-time check (Fortran `use` resolves) | Compile error |
| 2D golden regenerated from broken solver | Golden reproduces itself | Caught by Phase 2.1 minigap-pattern assertion (independent of golden) |
| Lecture script regex not matching new Fortran output format | `RuntimeError` in `lecture_13_topological.py` | Tighten in Phase 2 |

## 7. Testing strategy

### 7.1 TDD discipline (Phase 1 critical path)

For each P1 physics fix in Phase 1, the commit sequence is:
1. **Red**: Write a failing unit test that asserts the desired behavior. Commit as `test(red)`.
2. **Green**: Implement the fix. Existing tests must still pass; new test passes. Commit as `fix(...)` or `feat(...)`.
3. **Refactor** (if applicable): Clean up without changing behavior. Existing + new tests still pass. Commit as `refactor(...)`.

TDD commit triples:

| Commit | Test | Fix |
|---|---|---|
| 1.2 / 1.3 / 1.4 | `test_kitaev_strict.pf` (assert M=-1/+1 strict on real fixture) | Lutchyn-Oreg implementation |
| 1.5 / 1.6 | `test_pairing_sign_cross_module.pf` (byte-identity) | `public` + `use` |
| 1.7 / 1.8 | `test_dense_qw_block_21_sign.pf` (Hermiticity) | Fix line 547 |
| 1.9 / 1.10 | `test_phs_qw_peierls_symmetric.pf` (PHS with Bx≠0 generic k) | Remove symmetric Peierls [REVERTED 2026-07-05 — symmetric Peierls is REQUIRED; restored] |

**Note on numbering**: Phase 0 commit (0.1) precedes all Phase 1 commits. The prerequisite refactor (1.1) verifies the `add_peierls_coo` convention. Task 1.10 was originally "remove symmetric Peierls" but was REVERTED 2026-07-05: the symmetric `add_peierls_coo(-B)` call is REQUIRED for class-D PHS at generic k with Bx≠0 (see ADR 0008 §3 amendment). Commits 1.11 and 1.12 are single-commit refactors (the half_wire_integral and complex-coherence fixes are 1-line changes with existing test coverage).

For Phases 2-4 (less physics-critical), commits bundle fix + test in a single commit.

### 7.2 Unit test additions

**New**:
- `test_kitaev_strict.pf` (Phase 1.2) — Lutchyn-Oreg strict discriminator
- `test_pairing_sign_cross_module.pf` (Phase 1.5) — byte-identity
- `test_dense_qw_block_21_sign.pf` (Phase 1.7) — Hermiticity
- `test_phs_qw_peierls_symmetric.pf` (Phase 1.9) — PHS at generic k with Bx≠0
- `test_majorana_polarization_soc.pf` (Phase 1.12) — complex coherence
- `test_half_wire_integral_both_ends.pf` (Phase 1.11) — MZM at either end

**Modified**:
- `test_kitaev_M_topological_returns_valid.pf` — becomes `test_kitaev_wrapper_returns_valid_marker.pf`
- `test_wire_pfaffian_witness.pf` — replace s1=s2=0 vacuous acceptance with `fail`
- `test_krylov_snapshots.pf` — regenerate `init_wire_peierls_ref_n144_k6` under symmetric-Peierls-removed convention

**COVERAGE annotations** (Phase 4.3): all 7 new pFUnit files + 5 modified files.

### 7.3 Integration test changes

**Modified** (Phase 2):
- `scripts/lecture_13_topological.py:196` — drop `bcrit_pfaffian = bcrit_curve` alias
- `scripts/lecture_13_topological.py:172-189` — read live `output/z2_phase_diagram.dat` not golden
- `tests/integration/test_lecture_13_acceptance_gate.sh:87` — add absolute window guard
- `tests/integration/test_lecture_13_acceptance_gate.sh:104` — replace literal-substring grep with regex anchor
- `tests/integration/test_lecture_13_acceptance_gate.sh:41` — replace hardcoded `OMP_NUM_THREADS=4` with `nproc`-capped

**Renamed** (Phase 2.5):
- `tests/integration/verify_majorana_polarization.py` → `tests/integration/plot_majorana_polarization_synthetic.py`
- Updated `tests/CMakeLists.txt` ctest label

### 7.4 Regression test changes

**Modified** (Phase 2.1):
- `tests/regression/data/wire_inas_gaas_bdg_topological_2d_phase.dat` — drop z2 column, keep minigap only
- `tests/regression/test_wire_bdg_topological_2d.sh` — assert minigap pattern (small near B_crit, large far from it)
- `tests/regression/data/wire_inas_gaas_bdg_topological_2d_transitions.dat` — drop or replace with minigap transitions

### 7.5 CI gating

**Acceptance gate** (Phase 2.4):
- `assert 0.5 <= BCRIT_MIN` — catches uniform regression to small B_crit
- `assert BCRIT_MAX <= 6.0` — catches uniform regression to large B_crit
- `assert (BCRIT_MAX - BCRIT_MIN) <= 2.0` — keeps witness-disagreement check
- `grep -E '(0\.[12][0-9]?\s*T|Gershgorin fallback|Auto energy-window)'` — catches legacy false-PASS values

**Coverage matrix** (Phase 4.4):
- Add `bdg_ldos_topological`, `bdg_spectral_function_topological`, `bdg_nambu_ldos_topological` cells to `validation_universe.yml`
- Add the new observables to the metadata list

**Ctest labels** (Phase 2.5, 3.5):
- `verification` label no longer includes `verify_majorana_polarization.py`
- New `plot` label for `plot_*_synthetic.py` files
- New `acceptance-gate` label for `test_lecture_13_acceptance_gate.sh`

### 7.6 Performance considerations

Lutchyn-Oreg requires eigendecomposition of H_bdg at two PHS-invariant momenta per call. For 16N×16N BdG (N spatial sites), this is O((16N)³) per call.
- Wire (N≈50 → 800×800): ~5e8 flops/call, sub-second on modern CPU
- Dense-QW (N≈100 → 1600×1600): ~4e9 flops/call, a few seconds

**Mitigation**: Lutchyn-Oreg is only needed for diagnostic purposes; regression test (minigap-pattern) doesn't require it. Verifier scripts that need it compute it once per (B, μ) point, not per sweep step.

**Regression test runtime budget**: No regression test should exceed 30s. Add timing assertions to `test_wire_bdg_topological_2d.sh` and `test_lecture_13_acceptance_gate.sh`.

## 8. Out of scope (deferred to follow-up PRs)

- **Module split of `topological_analysis.f90`** (1820 lines → 4 files of ~450 each). Filed as separate follow-up. Pre-existing 1437-line bloat is not addressed in this PR.
- **U13 wire periodic/Bloch BdG** construction with Peierls-twist under Bx≠0. Required for full wire Pfaffian sweep.
- **`bdq_spectral` extension to dense-QW confinement** (currently wire-only).
- **M=±1 on spinless chain (structurally degenerate)** — already provided by Issue 05/07; the Lutchyn-Oreg fix in this PR makes strict assertions possible.
- **`csr_spectral_lorentzian_sum` shared helper extraction** (internal refactor, no behavior change).
- **Project-wide lint check for bare `stop 1`** (add as separate hygiene PR after Phase 4.1 lands).

## 9. Acceptance criteria

1. **PHS oracle green at generic k** under all four field combinations, both wire AND dense-QW. Phase 1.7, 1.8, 1.9, 1.10.
2. **Both builders produce Hermitian BdG matrices** (verified by Phase 1.7 test). Phase 1.8.
3. **Kitaev harness AE1 closed-form certification strict**: M=-1 for |μ|<2t, M=+1 for |μ|>2t on the **real Lutchyn-Oreg-compatible BdG fixture** (2-band s-wave SC with explicit hopping/onsite terms; not the synthetic diagonal-in-(c,c†) harness). Phase 1.2-1.4.
4. **Sticlet polarization with complex coherence**: P_M saturates for analytical MZM spinor; half-wire integral = max(left, right). Phase 1.12, 1.13.
5. **Dense-QW BdG rung**: ±E pairing, M=±1 strict (per AC3), Hermiticity (per AC2), edge-localized MZMs with P_M>0.95 + half-wire >0.4 at B=2·B_crit. Phase 1 (cumulative).
6. **Wire rung AE3**: min_gap traces open→close→reopen at B_crit≈2.8T; PHS oracle green at generic k with Bx≠0. Phase 1.10.
7. **2D minigap colormap regression**: pattern assertion (small near B_crit, large far from it), NOT z2=±1 numerical sign. Phase 2.1.
8. **4-witness acceptance gate**: all 4 witnesses live (not constants), range check + absolute window guard + regex anchor. Phase 2.2-2.6.
9. **Lecture 13**: hand-edited table removed; grep guard added; LDOS zero-bias claim reconciled. Phase 4.5, 4.6.
10. **All pFUnit tests have COVERAGE annotations**; validation_universe.yml has new cells. Phase 4.3, 4.4.

## 10. Risk register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Lutchyn-Oreg eigendecomposition slower than Pf | Medium | Low | Sub-second for typical sizes; verifier runtime budget verified in 7.6 |
| Phase 1.10 symmetric Peierls removal breaks wire PHS oracle [REVERTED 2026-07-05] | Medium | High | Phase 1.10 originally planned to remove the symmetric call; it was REVERTED because the call IS required for class-D PHS at generic k with Bx≠0 (rel_resid 1.25e-1 → 0 with the call, → 0.12578 without). See ADR 0008 §3 amendment. |
| Phase 2.1 minigap golden regenerated from broken solver | Low | Medium | Cross-check against `verify_wire_bdg_topological.py` live output; assertion on pattern shape (not values) |
| Phase 4.7 PRD archive breaks in-flight work | Low | Low | Confirm `.scratch/bdg-majorana-validation/` not referenced by any open issue |
| Phase 1.4 Lutchyn-Oreg API change breaks `wire_pfaffian_witness_sweep` | Medium | Medium | Wire Pfaffian witness uses the same `kitaev_majorana_number` internally; verify in Phase 1.4 TDD |
| Phase 3.4 default-near-zero-frac change breaks dense-QW regression | Low | Medium | Verify both `eval_bdg_point` and `q_zero_tol` use same threshold before commit |
| TDD on Phase 1 surfaces unexpected test infrastructure issues | Medium | Medium | Buffer 2-3 days in Phase 1 timeline for test scaffolding work |

## 11. References

- **Source review**: 6-dimension adversarial review of PR #41 (2026-07-01). Findings as P1/P2/P3.
- **ADR 0007**: `docs/adr/0007-bdg-hole-block-canonical-convention.md` — canonical hole-block convention (precedent for ADR 0008).
- **MEMORY.md**: `codebase-doc-drift-prevention` rule (3rd drift event will be added in Phase 4.8).
- **CLAUDE.md**: Engineering Principles (DRY, SOLID, YAGNI), Code Conventions (`error stop`, `iso_fortran_env`, pFUnit single-line asserts).
- **PRD**: `.scratch/bdg-majorana-validation/PRD.md` (will be archived in Phase 4.7).
- **Prior review-fix specs**: `docs/superpowers/specs/2026-06-13-eigensolver-review-fixes-design.md`, `docs/superpowers/specs/2026-06-21-pr39-review-fixes-design.md`.
- **Implementation plan**: `docs/superpowers/plans/2026-07-01-bdg-p1-fix.md` (to be created by `writing-plans` skill after this spec is approved).
- **ADR 0008**: `docs/adr/0008-bdg-p1-invariants.md` (to be created in Phase 0.1; locks Lutchyn-Oreg, pairing_sign DRY public API, symmetric Peierls KEPT [amended 2026-07-05 — `conjg()` alone is insufficient], minigap-pattern regression).

---

**SPEC-COMPLETE** — awaiting user review before transitioning to `writing-plans`.