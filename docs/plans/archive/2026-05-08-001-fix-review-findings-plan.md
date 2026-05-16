---
title: "fix: Address code review findings from topological bug fixes"
type: fix
status: active
date: 2026-05-08
---

# fix: Address code review findings from topological bug fixes

## Summary

Fix 8 findings from the multi-agent code review of the topological/magnetic bug-fix commits (`46be100..f875ab1`): a BdG gap convention inconsistency, missing test coverage for the Peierls fallback branch, a Majorana test fixture that doesn't exercise the multi-column OOB bug, and several test coverage gaps. Also add missing `contiguous` attributes on hot-path COO array arguments.

---

## Requirements

- R1. BdG gap definition must be consistent across all three code paths: wire single-point, QW single-point, and wire sweep (full gap = 2 × min|E|).
- R2. Every committed code change in the Peierls, Majorana, and BHZ fixes must have test coverage that would catch a regression of the original bug.
- R3. `contiguous` attribute on hot-path assumed-shape array arguments per CLAUDE.md conventions.
- R4. All 66 existing tests must continue to pass.

---

## Scope Boundaries

- Does NOT refactor dead code or file-length issues flagged by the review (pre-existing, out of scope).
- Does NOT add QW Fu-Kane sweep unit tests (review finding #4, P1) — this requires non-trivial mock infrastructure and is deferred.
- Does NOT address advisory P3 findings (is_z2_transition over-flagging, FEAST window, BHZ unit conversion comments).

### Deferred to Follow-Up Work

- QW Fu-Kane sweep per-point unit test: requires mocking `compute_z2_fukane_qw_result` or a full QW config fixture. Deferred to a future iteration.
- `wire_bdg` sweep regression config: noted by agent-native reviewer. Deferred.

---

## Key Technical Decisions

- **Full-gap convention (2 × min|E|):** The BdG gap will use `2.0_dp * bdg_zero_energy_gap(eigenvalues)` everywhere. This matches common condensed-matter convention (Tinkham) and aligns the wire path with QW and sweep paths.
- **Column-major Peierls test fixture:** The test coordinate setup will use correct column-major flat indexing `(iy-1)*nx + ix` to derive `grid%coords`, not the coincidental `mod(i-1,ny)+1` formula.
- **Multi-column Majorana test:** Change from `nx=1` to `nx=3` so `nspatial > ny`, ensuring the test would fail if `grid%z` were used instead of `grid%coords`.

---

## Implementation Units

- U1. **Unify BdG gap convention to full gap (2 × min|E|)**

**Goal:** Make the wire single-point BdG path use the same full-gap convention as QW and sweep paths.

**Requirements:** R1

**Dependencies:** None

**Files:**
- Modify: `src/apps/main_topology.f90`
- Test: `tests/unit/test_bdg_hamiltonian.pf`

**Approach:** Update `run_bdg_wire` (line 534) from `bdg_zero_energy_gap(eigvals_bdg)` to `2.0_dp * bdg_zero_energy_gap(eigvals_bdg)`. Also update the existing test `test_bdg_gap_uses_min_abs_energy` to reflect the full-gap convention — the test should expect `0.002` (2 × 0.001) instead of `0.001`.

**Patterns to follow:** The QW path at line 671 and sweep path at line 1206 already use `2.0_dp * minval(abs(eigvals_bdg))`.

**Test scenarios:**
- Happy path: eigenvalues `[-0.8, -0.3, 0.001, 0.3, 0.8]` → gap = 0.002
- Edge case: empty array → gap = 0.0 (guard returns 0, multiplied by 2 still 0)
- Edge case: single eigenvalue `[0.5]` → gap = 1.0
- Edge case: exact zero `[-0.3, 0.0, 0.3]` → gap = 0.0 (gap closed)

**Verification:** All existing tests pass. The updated test passes with the new expected value.

---

- U2. **Add `contiguous` attributes to COO array arguments**

**Goal:** Add `contiguous` attribute to hot-path assumed-shape array arguments in `magnetic_field.f90` per CLAUDE.md conventions.

**Requirements:** R3

**Dependencies:** None

**Files:**
- Modify: `src/physics/magnetic_field.f90`

**Approach:** Add `contiguous` to `coo_vals(:)`, `coo_row(:)`, `coo_col(:)` in both `add_zeeman_coo` (around line 21) and `add_peierls_coo` (around line 61).

**Test scenarios:**
- Test expectation: none — attribute-only change, no behavioral change. Existing tests confirm correctness.

**Verification:** Build succeeds. All existing tests pass.

---

- U3. **Fix Peierls test fixture and add fallback branch test**

**Goal:** Fix the physically incorrect coordinate mapping in `test_peierls_y_axis_2d_wire` and add a test for the Peierls fallback y-index formula (integer division by nx) which currently has zero coverage.

**Requirements:** R2

**Dependencies:** None

**Files:**
- Modify: `tests/unit/test_magnetic_field.pf`

**Approach:**
1. Fix the existing `test_peierls_y_axis_2d_wire` coordinate setup to use column-major indexing: `iy = (i-1)/nx + 1`, then `grid%coords(2, i) = grid%z(iy)`.
2. Add a new test `test_peierls_fallback_y_index` that creates a 1D grid (ndim=1, no grid%coords) with nx=3, ny=5, sets up an off-diagonal COO entry between sites in different y-rows, and verifies the fallback formula `(flat_i-1)/nx + 1` produces correct y-coordinates from grid%z.

**Patterns to follow:** Column-major flat_idx: `ij = (iy-1)*nx + ix` in `geometry.f90`.

**Test scenarios:**
- Happy path (existing test, fixed fixture): 3×5 grid, site 1→site 8, verify Peierls phase matches expected value using column-major coords
- Happy path (new test): 1D grid fallback, sites at different y-positions, verify y-coordinates extracted via integer division match grid%z values
- Edge case: verify iy computation stays in [1, ny] for all flat indices in [1, nx*ny]

**Verification:** Both Peierls tests pass. Build succeeds. Full suite passes.

---

- U4. **Fix Majorana test to use multi-column grid**

**Goal:** Change the Majorana 2D wire test from nx=1 to nx=3 so the test would actually fail if grid%z were used instead of grid%coords (the original bug).

**Requirements:** R2

**Dependencies:** None

**Files:**
- Modify: `tests/unit/test_edge_states.pf`

**Approach:** Change `test_majorana_coords_2d_wire` from `nx=1, ny=20` to `nx=3, ny=7` (21 spatial points). Set `nspatial = 21`, `half_n = 8 * nspatial`. Build grid%coords using column-major ordering. Place exponential decay along y-axis. Tighten the tolerance assertion from 3.0 to 2.0.

**Patterns to follow:** Column-major `flat_idx` from `geometry.f90`. Existing `test_majorana_profile_synthetic` for 1D Majorana test patterns.

**Test scenarios:**
- Happy path: 3×7 grid with exponential decay along y, xi fitted using physical coords should match expected value
- Edge case: verify the test would fail if grid%z (length 7) were used instead of grid%coords (length 21) — out-of-bounds or wrong coordinates

**Verification:** Test passes. Build succeeds. Full suite passes.

---

- U5. **Add backward B+D hopping structural test**

**Goal:** Add explicit structural assertions for backward B+D hopping connectivity in the BHZ wire Hamiltonian, complementing the existing forward-only test.

**Requirements:** R2

**Dependencies:** None

**Files:**
- Modify: `tests/unit/test_edge_states.pf`

**Approach:** Extend `test_bhz_forward_hopping_connects_neighbors` (or add a new test) to verify that for site 2, backward B+D connects to site 1 via column `(i-2)*4 + mod(row+2,4) + 1`. Check specific orbital connections for row=1 and row=2.

**Patterns to follow:** Existing `test_bhz_forward_hopping_connects_neighbors` for forward A-term and B+D assertion pattern. `test_bhz_wire_hamiltonian_is_hermitian` implicitly covers this but structural assertions catch regressions more directly.

**Test scenarios:**
- Happy path: site 2 row 1, backward B+D connects to site 1 at column `(0)*4 + mod(3,4)+1 = 4`
- Happy path: site 2 row 2, backward B+D connects to site 1 at column `0 + mod(4,4)+1 = 1`
- Negative: verify no same-site backward B+D connection (the column should not equal `(i-1)*4 + ...`)

**Verification:** Test passes alongside existing BHZ tests. Build succeeds. Full suite passes.

---

## System-Wide Impact

- **API surface parity:** The BdG gap change affects the wire single-point output (`run_bdg_wire`) — downstream consumers comparing wire gap values against QW gap values will now see consistent units.
- **Unchanged invariants:** All 66 existing tests must pass unchanged (except the BdG gap test which gets an updated expected value).

---

## Risks & Dependencies

| Risk | Mitigation |
|------|------------|
| BdG gap convention change may confuse users comparing against old wire output | The old wire output used an incorrect convention (adjacent spacing, not min\|E\|). The new convention is consistent across all paths. |
| Majorana test fixture change (nx=1→3) may affect tolerance thresholds | Tighten tolerance conservatively (3.0→2.0) and verify the test still passes. |

---

## Sources & References

- Code review report (ce-code-review, run 20260508-080346-1e46ca2d)
- Original bug-fix plan: `docs/superpowers/plans/2026-05-07-topological-magnetic-bugfixes.md`
- Column-major flat_idx: `src/math/geometry.f90` (`flat_idx` function)
