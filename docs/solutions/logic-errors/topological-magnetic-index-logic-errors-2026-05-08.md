---
title: "Five index-arithmetic and convention bugs in topological analysis and magnetic field modules"
date: 2026-05-08
category: logic-errors
module: topological_analysis_and_magnetic_field
problem_type: logic_error
component: service_object
severity: high
symptoms:
  - "Peierls phase crashes or produces wrong values for multi-column 2D wires (nx > 1)"
  - "BHZ wire Hamiltonian forward hopping couples orbitals within same site instead of neighbor"
  - "Majorana localization fit uses 1D-length coordinate array for 2D grid, causing OOB access"
  - "QW Fu-Kane Z2 phase diagram is flat (same Z2 value at all B, mu points)"
  - "BdG superconducting gap reports half-gap in wire path vs full gap in QW/sweep paths"
root_cause: logic_error
resolution_type: code_fix
related_components:
  - bdg_hamiltonian
  - green_functions
tags:
  - topological-analysis
  - bdg
  - peierls-substitution
  - majorana
  - out-of-bounds
  - fortran
  - index-arithmetic
  - bhz
---

# Five index-arithmetic and convention bugs in topological analysis and magnetic field modules

## Problem

Five logic errors in the topological analysis (`topological_analysis.f90`) and magnetic field (`magnetic_field.f90`) modules produced incorrect topological phase diagrams, silently wrong Peierls phases in 2D wires, and inconsistent BdG superconducting gaps. All five were discovered by automated code review (Codex and ce-code-review) during the BdG topological superconductivity feature development on branch `feature/bdg-topological-superconductivity`.

## Symptoms

- Peierls phase crashes or produces garbage for multi-column 2D wires (`nx > 1`): `grid%z(flat_index)` reads out of bounds when `flat_index` exceeds `ny`.
- BHZ edge states have wrong spatial structure: forward hopping terms couple orbitals within the same site instead of connecting to the neighboring site, making the Hamiltonian non-Hermitian.
- Majorana localization length is wrong in 2D wires: `grid%z` (length `ny`) is too short for `nspatial = nx*ny` sites.
- QW Fu-Kane Z2 phase diagram is flat/constant: single pre-loop evaluation broadcast to all `(B, mu)` grid points.
- BdG wire gap is half the physical value: `min|E|` instead of `2 * min|E|` (Tinkham convention).

## What Didn't Work

- **Increasing BHZ wire length to resolve edge states**: Before discovering the same-site hopping bug, increasing N from 100 to 500 did not fix incorrect Z2 invariants because the Hamiltonian topology was fundamentally wrong (intra-site rather than inter-site connections).
- **Tuning FEAST energy window for gap resolution**: The gap-closing detection failure was not a resolution issue but a definition issue (adjacent spacing vs zero-energy distance).
- **Adding more sweep points for the Fu-Kane diagram**: Finer B/mu grid resolution did not help because the single pre-loop evaluation made every point identical.
- **Peierls fallback y-index using `mod(flat-1, ny)+1`**: This formula is only correct for `nx=1`. Integer division `(flat-1)/nx + 1` is the correct column-major row extraction (session history: the plan text itself was wrong about this; the code fix used the correct formula).

## Solution

### Bug 1: Peierls y-axis OOB (`magnetic_field.f90:98-106`)

**Before** (wrong for `nx > 1`):
```fortran
iy_i = mod(flat_i - 1, grid%ny) + 1       ! only correct for nx=1
```

**After**:
```fortran
if (grid%ndim == 2 .and. allocated(grid%coords)) then
  y_i = grid%coords(2, flat_i) * 1.0e-10_dp
  y_j = grid%coords(2, flat_j) * 1.0e-10_dp
else
  iy_i = (flat_i - 1) / grid%nx + 1       ! integer division by nx
  iy_j = (flat_j - 1) / grid%nx + 1
  y_i = grid%z(iy_i) * 1.0e-10_dp
  y_j = grid%z(iy_j) * 1.0e-10_dp
end if
```

### Bug 2: BHZ forward hopping same-site (`topological_analysis.f90:769-794`)

**Before** (same-site):
```fortran
coo_col(nnz_offset) = (i-1)*4 + 5 - row   ! stays on site i
coo_col(nnz_offset) = (i-1)*4 + mod(row, 4) + 1
```

**After** (neighbor site `i+1`):
```fortran
coo_col(nnz_offset) = i*4 + 5 - row        ! connects to site i+1
coo_col(nnz_offset) = i*4 + mod(row, 4) + 1
```

Backward A-term also fixed to match for Hermiticity: `(i-2)*4 + 5 - row` with same sign rule.

### Bug 3: Majorana profile OOB (`topological_analysis.f90:656-666`)

**Before**:
```fortran
if (allocated(grid%z)) then
  xx = grid%z             ! length ny, but nspatial = nx*ny
```

**After**:
```fortran
if (grid%ndim == 2 .and. allocated(grid%coords)) then
  allocate(xx(nspatial))
  xx(1:nspatial) = grid%coords(2, 1:nspatial)   ! y-coords for all sites
else if (allocated(grid%z)) then
  xx = grid%z
```

### Bug 4: QW Fu-Kane evaluated once (`main_topology.f90:~1074-1106`)

**Before** (computation outside loop):
```fortran
call compute_z2_fukane_qw_result(cfg_in, ...)
do iB = 1, nB
  do iMu = 1, nMu
    z2_map(iMu, iB) = z2_val    ! same value everywhere
```

**After** (per-point evaluation inside loop):
```fortran
do iB = 1, nB
  B_val = cfg_in%topo%gap_sweep_B_min + real(iB - 1, kind=dp) * dB
  do iMu = 1, nMu
    mu_val = cfg_in%topo%gap_sweep_mu_min + real(iMu - 1, kind=dp) * dmu
    cfg_local = cfg_in
    cfg_local%bdg%B_vec = [0.0_dp, 0.0_dp, B_val]
    cfg_local%bdg%mu = mu_val
    call compute_z2_fukane_qw_result(cfg_local, ...)
```

### Bug 5: BdG gap convention (`main_topology.f90:545`)

**Before** (half-gap):
```fortran
result%min_gap = bdg_zero_energy_gap(eigvals_bdg)
```

**After** (full SC gap, consistent with QW and sweep paths):
```fortran
result%min_gap = 2.0_dp * bdg_zero_energy_gap(eigvals_bdg)
```

## Why This Works

**Peierls**: The column-major flat index is `flat = (iy-1)*nx + ix`. Integer division `(flat-1)/nx + 1` correctly extracts the row index `iy` for any `nx`. The 2D path uses `grid%coords(2, flat)` which provides the y-coordinate for every spatial site.

**BHZ hopping**: The BHZ wire has 4 orbitals per site. Site `i` occupies rows `(i-1)*4 + 1` through `i*4`. Forward hopping to site `i+1` must target `i*4 + offset`, not `(i-1)*4 + offset`.

**Majorana**: For 2D wires, `grid%coords(2, :)` provides y-coordinates for all `nx*ny` flat sites, while `grid%z` only has `ny` entries.

**Fu-Kane**: Moving the computation inside the loop with per-point `B_val` and `mu_val` produces the correct phase diagram with topological transitions.

**BdG gap**: The quasiparticle spectrum `E = +/- sqrt(xi^2 + Delta^2)` has minimum `|E| = Delta` at `xi=0`. The full excitation gap is `2*Delta` (one quasiparticle at +E, one at -E), matching the Tinkham convention.

## Prevention

**Test cases added** (all in `tests/unit/`):

| Bug | Test | What it catches |
|-----|------|-----------------|
| 1 | `test_peierls_y_axis_2d_wire` | 3x5 grid verifies phase uses `grid%coords(2,:)` |
| 1 | `test_peierls_fallback_y_index` | 1D fallback verifies `(flat-1)/nx+1` gives correct `iy` |
| 2 | `test_bhz_forward_hopping_connects_neighbors` | Checks H(1,8) exists, H(1,4) absent |
| 2 | `test_bhz_backward_hopping_connects_neighbors` | Checks backward A/B+D target site 1 from site 2 |
| 2 | `test_bhz_wire_hamiltonian_is_hermitian` | `max|H - H^dagger|` check (same-site hopping breaks this) |
| 3 | `test_majorana_coords_2d_wire` | 3x7 grid (nspatial=21 > ny=7), checks xi uses physical coords |
| 5 | `test_bdg_gap_full_convention` | `[-0.8, -0.3, 0.001, 0.3, 0.8]` -> gap = 0.002 |
| 5 | `test_bdg_gap_exact_zero` | `[-0.3, 0.0, 0.3]` -> gap = 0.0 |

**Strategies**:

- Use `grid%coords(dim, :)` for 2D grids instead of `grid%z`. Reserve `grid%z` for 1D-only code paths.
- For sparse matrix assembly with site-wise blocks, verify forward hopping targets `i*block_size + offset` and backward targets `(i-2)*block_size + offset`. Add structural CSR entry tests, not just eigenvalue tests.
- Before implementing parameter sweeps, write the loop skeleton first with the computation inside it. Code review should flag any computation before a loop whose result is used inside without per-iteration parameterization.
- For gap/energy quantities, add synthetic tests with known eigenvalues checked against hand-computed expected values. This catches wrong definitions and missing factors of 2.
- Build with `-fcheck=all` (gfortran) during development to catch OOB array accesses at the point of occurrence.

## Related Issues

- Original bug-fix spec: `docs/superpowers/specs/archive/2026-05-07-topological-magnetic-bugfixes-design.md`
- Bug-fix plan: `docs/superpowers/plans/archive/2026-05-07-topological-magnetic-bugfixes.md`
- Follow-up review findings plan: `docs/plans/archive/2026-05-08-001-fix-review-findings-plan.md`
- BdG topological design: `docs/plans/archive/2026-04-27-bdg-topological-superconductivity-design.md`
