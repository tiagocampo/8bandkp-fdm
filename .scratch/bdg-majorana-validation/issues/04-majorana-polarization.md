# Issue 04 — Majorana polarization as derivation-binding Sticlet coherence (U6)

**Parent PRD**: `/data/8bandkp-fdm/.scratch/bdg-majorana-validation/PRD.md` (Unit U6)
**Issue source**: `/data/8bandkp-fdm/.superpowers/sdd/issue-04-brief.md`
**Plan & context**: `/data/8bandkp-fdm/.superpowers/sdd/understand-report.md`
**Phase**: PR-B, fan-out with Issue 06

## What to build

A pure-function Majorana polarization that discriminates true Majorana zero modes from generic near-zero states. The polarization is the Sticlet off-diagonal electron-hole coherence `P_M(n) = 2·|Σ_σ s_σ u_{nσ} v*_{nσ}|`, with **spin-sector signs `s_σ` derived once for the code's Nambu ordering** (KTD7, ADR 0007) — not assumed from the published formula. The relative minus sign between spin sectors is the common implementation pitfall; pinning it to the codebase's Nambu layout prevents sign errors.

The polarization routine:
1. Site-resolved `P_M(n)` as a pure function of a BdG eigenvector.
2. The half-wire integral (sum over one end) as the headline scalar (saturates at ~0.5 for a true MZM).
3. The charge polarization `⟨τ_z⟩ = (|u|²−|v|²)/(|u|²+|v|²)` is computed alongside as a documented contrast — it vanishes at a true MZM and contradicts AE2's "polarization ≈ 1", so it's not the discriminator.
4. Co-located with the existing `compute_majorana_profile` (shares the `band_major_row` indexing).

Two-tier threshold for assertions:
- **Headline (well inside topological regime):** `P_M > 0.95` per site + half-wire integral > 0.4 at `B = 2·B_crit`.
- **Boundary:** assertion flips to minigap (via Issue 00's evaluator), not P_M, at `B_crit` itself.

## Acceptance criteria

- [ ] Sticlet Eq. 5 analytical MZM spinor → `P_M` saturates (Covers AE2, polarization ≈ 1).
- [ ] A pure-electron near-zero state → `P_M ≈ 0` (the discriminator works).
- [ ] Site-resolved `P_M` is peaked at the two wire ends with opposite sign for a true MZM.
- [ ] Half-wire integral ≈ 0.5 for a true MZM; ≈ 0 for an accidental near-zero mode.
- [ ] `⟨τ_z⟩` charge polarization → 0 at the MZM (documenting that it is not the discriminator).
- [ ] The spin-sector signs `s_σ` are documented in the routine header with reference to KTD7 and ADR 0007 — derived, not assumed.
- [ ] `tests/unit/test_majorana_polarization.pf` (new) green; the routine is callable from a unit test with a synthetic BdG eigenvector (no `simulation_config`).
- [ ] The polarization does **not** replace the hardcoded `0.001·δ₀` near-zero threshold in Issue 00's evaluator — the threshold stays for the near-zero count; polarization replaces the threshold as the MZM discriminator for downstream assertions.
- [ ] Per-task code review + spec compliance review clean.

## Pre-existing state (from Understand report)

### `src/physics/topological_analysis.f90` (1437 lines)
- `compute_majorana_profile` (lines 717–793): takes `evec_bdg`, `grid`, `energy_tol`, `half_n_in`, returns `xi` (edge localization length); uses `band_major_row` for indexing. Co-locate the new polarization routine near this one.
- `bdg_zero_energy_gap` (lines 891–900): `pure function` returns `minval(abs(eigenvalues))`.
- Module exports a handful of `compute_*` routines; follow the same naming convention.

### Issue 00's evaluator (`src/physics/bdg_observables.f90`)
- `eval_bdg_point(eigenvalues, params)` returns `(minigap, near_zero_count, invariant_flag)`.
- `bdg_eval_params_t` has `delta_0` and `near_zero_frac` fields.
- This is the boundary assertion (per AC #10 — Issue 00's evaluator stays for near-zero count).

## Constraints from CLAUDE.md + ADRs

- **ADR 0007 (hole-block canonical)**: Spin-sector signs `s_σ` are derived from the canonical Nambu convention (`-conjg(H₀(-k))`); both spin sectors' positions in the BdG eigenvector are determined by the Kramers pairing `pairing_partner`/`pairing_sign` from `bdg_hamiltonian.f90`.
- **ADR 0001 (fat derived type)**: Pure functions only.
- **ADR 0003 (build-and-solve in app)**: Pure polarization function in physics.
- **CLAUDE.md Boundaries**: NO approval gate (does not modify Hamilton-construction code).
- **CLAUDE.md Code Conventions**: F2018; `private` default + explicit `public ::` exports; `error stop` not `stop 1`; no `goto`; `<= 300 lines/file`, `<= 50 lines/function`; pFUnit `@assertEqual`/`@assertTrue` MUST be single-line.

## File ownership (exhaustive)

### New files (you create)
- `tests/unit/test_majorana_polarization.pf` — pFUnit unit tests for the polarization routine.

### Modified files
- `src/physics/topological_analysis.f90` — add `majorana_polarization` pure function (or co-locate as a private module function); co-locate with `compute_majorana_profile` (shares `band_major_row` indexing). Optionally add a small comment in `compute_majorana_profile` noting the polarization is the new discriminator.
- `src/physics/AGENTS.md` — one-line note about the polarization routine.
- `tests/CMakeLists.txt` — `add_pfunit_ctest(test_majorana_polarization ...)` with LABEL "unit".

### NOT modified (out of scope)
- `src/physics/bdg_hamiltonian.f90` (Issue 03 owns; do NOT touch).
- `src/physics/bdg_observables.f90` (Issue 00 owns).
- `src/apps/main_topology.f90` (Issue 05/07 own the wiring).
- `src/core/defs.f90` (Issue 06/07 own).

## Suggested API (illustrative)

```fortran
module topological_analysis
  ...
  public :: majorana_polarization, polarization_result_t

  type :: polarization_result_t
    real(dp), allocatable :: P_M(:)        ! site-resolved P_M(n)
    real(dp), allocatable :: tau_z(:)      ! site-resolved <tau_z>(n)
    real(dp) :: half_wire_integral         ! sum over one end (saturates ~0.5 for MZM)
    real(dp) :: total_P_M                  ! integral over full wire (sanity)
  end type

  ! Pure function: site-resolved P_M for a BdG eigenvector
  ! s_σ derived from KTD7 Nambu ordering: σ=1 (up) ↔ pair_partner=3, σ=2 (down) ↔ pair_partner=4
  ! (the relative sign convention is documented in the header; per ADR 0007)
  pure function majorana_polarization(evec_bdg, n_orbitals_per_site, n_sites, &
                                       kramers_pair_idx) result(pol)
    complex(dp), intent(in) :: evec_bdg(:)
    integer, intent(in) :: n_orbitals_per_site, n_sites
    integer, intent(in) :: kramers_pair_idx(:)  ! 1-based spin sector index for each orbital
    type(polarization_result_t) :: pol
  end function
end module
```

## Tests required

### `tests/unit/test_majorana_polarization.pf`
1. **Sticlet Eq. 5 analytical MZM spinor** (synthetic): `P_M` saturates to 1 at each site; half-wire integral ≈ N/2 (the "≈ 0.5 for a true MZM" headline value).
2. **Pure-electron near-zero state**: synthesize a BdG eigenvector with v=0, u nonzero; `P_M ≈ 0`.
3. **Pure-hole near-zero state**: synthesize a BdG eigenvector with u=0, v nonzero; `P_M ≈ 0`.
4. **Site-resolved peaked-at-ends**: MZM spinor → P_M peaks at the two wire ends with comparable magnitudes.
5. **Half-wire integral**: for a true MZM, `half_wire_integral > 0.4` (per AC); for an accidental near-zero mode, `half_wire_integral < 0.1`.
6. **`⟨τ_z⟩ → 0` at MZM**: charge polarization vanishes for true MZM (documented contrast).
7. **Spin-sector sign check**: verify the `s_σ` signs by constructing two spin sectors with opposite handedness and confirming P_M computes the absolute sum (per ADR 0007 derivation).

## TDD discipline (mandatory)

1. Write `tests/unit/test_majorana_polarization.pf` FIRST.
2. Confirm test fails (function doesn't exist).
3. Implement `majorana_polarization` in `topological_analysis.f90`.
4. Confirm test passes.
5. Run full unit suite to verify no regressions (38/41 expected — 3 follow-up failures from Issue 03 are pre-existing and out of scope).

## Build commands (use these exact forms)

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit --output-on-failure
```

## Out of scope (must NOT do)

- Do NOT modify `src/physics/bdg_hamiltonian.f90` (Issue 03 owns; do NOT touch).
- Do NOT modify `src/physics/bdg_observables.f90` (Issue 00 owns).
- Do NOT modify `src/apps/main_topology.f90` (wiring happens in Issue 05/07).
- Do NOT modify `src/core/defs.f90` (Issue 06/07 own).
- Do NOT add `sweep_model` enum value.
- Do NOT add new TOML fields.

## Report file path

Write your full report to: `/data/8bandkp-fdm/.superpowers/sdd/issue-04-report.md`

The report must include:
- Status (DONE / DONE_WITH_CONCERNS / BLOCKED / NEEDS_CONTEXT)
- Files modified
- TDD evidence (RED → GREEN) with exact commands and output
- Test results (full unit suite pass count; pre-existing 3 follow-up failures expected)
- Spin-sector sign derivation explicitly documented with KTD7 / ADR 0007 reference
- Commit SHA
- Self-review findings
- Concerns

## Risk notes

- The spin-sector sign `s_σ` derivation is THE most common Sticlet implementation pitfall. Per ADR 0007's Nambu convention (`-conjg(H₀(-k))`), the relative sign between the two spin-sector coherences is determined by the Kramers pairing. Document the derivation in the function header with explicit formulas; don't assume from the published Sticlet paper.
- The half-wire integral "saturates at ~0.5" applies to a continuous BdG spectrum (continuous Majorana localization length). For a discrete tight-binding lattice, the value may differ; verify with the analytical MZM spinor test (AC #1).
- The polarization does NOT replace the `0.001·δ₀` near-zero threshold in `eval_bdg_point`. The polarization is a discriminator for assertions (AE2's "polarization ≈ 1"); the threshold stays for the near-zero count (Issue 00's evaluator).

---

## Outcome (as executed)

- **Routine**: `majorana_polarization` in `src/physics/topological_analysis.f90` (co-located with `compute_majorana_profile`).
- **Spin-sector signs `s_σ`**: derived from KTD7 Nambu ordering per ADR 0007; documented in routine header with explicit KTD7 reference.
- **Headline values**: Sticlet Eq. 5 analytical MZM spinor → `P_M ≈ 1` per site; half-wire integral ≈ 0.5 for true MZM.
- **Discriminator verified**: pure-electron and pure-hole near-zero states → `P_M ≈ 0`.
- **Charge polarization `⟨τ_z⟩` → 0** at MZM (documented as contrast, not discriminator).
- **Tests**: `tests/unit/test_majorana_polarization.pf` green.
- **Fix Round 1**:
  - **Critical 1**: added distributed MZM half-wire test (tested with spatially-extended MZM, not just at boundary).
  - **Important 1**: replaced zero-fill on size mismatch with `error stop`.
  - **Size-mismatch error-stop test**: added explicit test for the error path.
  - **Critical 2**: deferred to Issue 05 (real BdG eigensolve polarization test).
- **Deferred**: Critical 2 (real BdG eigensolve polarization test on dense-QW) → Issue 05.
- **Side effects**: no impact on existing tests; co-location shares `band_major_row` indexing cleanly.