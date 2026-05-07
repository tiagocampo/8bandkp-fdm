# Topological & Magnetic Bug Fixes

Date: 2026-05-07
Status: Draft

## Problem

Five open Codex review findings on PR #13, all in the topological analysis and magnetic field code paths. Four are confirmed bugs; one (Peierls y-axis) needs verification. Current tests pass at 66/66 because the affected code paths lack regression coverage.

## Findings

### F1 (P1): Peierls y-axis indexing for multi-column wires

**File:** `src/physics/magnetic_field.f90:90-91`
**Reporter:** Codex (follow-up to earlier fix)

The previous fix converted band-major to flat site index via `mod(row-1, ngrid)+1` and then indexes `grid%z` with that flat site. In wire mode `grid%z` has length `ny` (y-axis only) while the flat site runs over `nx*ny`. For multi-column wires (`nx > 1`) this either reads past `grid%z` or uses wrong y-coordinates.

**Fix:** Convert the flat site to its y-component: `iy = mod(flat_site - 1, ny) + 1`, then index `grid%z(iy)`. This requires passing `ny` (or the full `spatial_grid`) into `add_peierls_coo`.

### F2 (P2): BHZ wire A-term forward hopping stays on-site

**File:** `src/physics/topological_analysis.f90:770`
**Reporter:** Codex, confirmed by code analysis

The column `(i-1)*4 + 5 - row` targets the same site `i` (mapping row 1->4, 2->3, 3->2, 4->1). The forward FD hopping for `A*k_z` should couple site `i` to site `i+1`. The backward term at line 781 correctly targets `(i-2)*4 + row` (site `i-1`).

Additionally, line 792 (`B_plus_D_over_dz2` forward hopping) uses `(i-1)*4 + mod(row,4) + 1` which also stays on-site.

**Fix:** Forward A-term column should be `i*4 + 5 - row` (site `i+1` with orbital mapping). Forward B+D column should be `i*4 + mod(row,4) + 1`. After fix, the BHZ wire Z2 regression tests (`regression_topology_bhz_z2_topological`, `regression_topology_bhz_z2_trivial`) must still pass — if they break, the previous "passing" state was accidental due to the on-site coupling being numerically similar at the test's parameter scale.

### F3 (P2): Majorana fit coordinates for 2D wires

**File:** `src/physics/topological_analysis.f90:656`
**Reporter:** Codex

`xx = grid%z` copies only `ny` y-axis coordinates, but `nspatial = nx*ny` and the fit loop indexes `xx(i)` up to `nspatial`. For `nx > 1` wires this reads out of bounds or fits against wrong coordinates.

**Fix:** The `spatial_grid` type has `coords(:,:)` allocated for 2D wires, with `coords(2,:)` holding y-coordinates at each of the `nx*ny` spatial points. Use `grid%coords(2,1:nspatial)` instead of `grid%z` when `grid%ndim == 2`. For 1D wires (`nx == 1`, `ndim == 1`) the existing `grid%z` path is correct. Add a guard: if `grid%ndim == 2 .and. allocated(grid%coords)` use `coords(2,:)`, else fall back to `grid%z`.

### F4 (P2): QW Fu-Kane sweep evaluates once

**File:** `src/apps/main_topology.f90:1087-1099`
**Reporter:** Codex

`compute_qw_fukane_gap_sweep` calls `compute_z2_fukane_qw_result` once before the double loop, then fills every `(B, mu)` grid point with the same `z2_val` and `min_gap`. The sweep ranges produce a constant phase diagram.

**Fix:** Move `compute_z2_fukane_qw_result` inside the double loop. Since `cfg` is `intent(in)`, create a local copy `cfg_local = cfg_in` and modify the sweep parameters (B-field via `cfg_local%topo%B_vec`, chemical potential via the relevant field) before each call. Follow the pattern used by `compute_wire_bdg_gap_sweep` (lines 1129-1135) which correctly evaluates per grid point.

### F5 (P2): BdG wire gap computed as min-adjacent-spacing

**File:** `src/apps/main_topology.f90:540-543`
**Reporter:** Codex

The gap is the minimum spacing between any consecutive eigenvalues. For BdG with multiple states in the `+/-5*Delta` window, this reports bulk-level degeneracies near zero as the "gap" even when the true superconducting gap is open. The QW and sweep paths use the closest-to-zero definition.

**Fix:** Change to `min(|E_i|)` for `i = 1, nev_found`. The gap is the smallest absolute eigenvalue, which measures distance from zero energy. This matches the QW path convention.

## Verification

After fixes:
1. All 66 existing tests must still pass
2. New unit tests for each fix (BHZ hopping structure, Majorana coordinates, Peierls y-axis)
3. Existing BHZ Z2 regression tests must produce physically correct results (may change if previous results were wrong)
4. Manual smoke test: `topology_bhz_z2_topological.cfg` should give Z2=1, `topology_bhz_z2_trivial.cfg` should give Z2=0
