# Wire g-Factor Fix: Commutator-Based Velocity Operator

Date: 2026-04-25
Branch: fix/wire-gfactor-commutator

## Problem

The wire g-factor calculation produces incorrect results that do not converge to the bulk limit for larger wires. The root cause is in the transverse perturbation construction (`g='g1'`/`g='g2'` in `ZB8bandGeneralized`), which uses the diagonal Kane P profile (`kpterms_2d(4)`) as a scalar approximation to dH/dk_x and dH/dk_y.

In the wire geometry, kx and ky are spatial operators (-i d/dx, -i d/dy), not scalar wavevectors. The scalar approximation misses contributions from Q, T, R, RC, A blocks (which give gradient operators) and is only exact in the limit of perfectly smooth envelope functions.

## Solution: Commutator Approach

Compute the velocity operator from the Heisenberg equation of motion:

```
v_alpha = [r_alpha, H] / (i hbar)
dH/dk_alpha = -i * [r_alpha, H]
```

Discretized for CSR matrices:

```
(dH/dk_alpha)_{i,j} = -i * (r_alpha_i - r_alpha_j) * H_{i,j}
```

This is an element-wise scaling of the full wire Hamiltonian by position differences. Same sparsity pattern as H. Automatically captures all contributions (PP, PM, Q, T, R, A, S, SC) without separate perturbation modes.

### Key properties

- Same-spatial-point entries: r_i - r_j = 0 -> band offsets correctly excluded
- x-neighbors: r_i - r_j = +/-dx -> FD stencil entries scaled (gives velocity)
- z-direction: all points share same z -> [z, H] = 0 -> must keep existing g3 approach

## Files to Modify

### 1. `src/physics/hamiltonianConstructor.f90`

**Add** subroutine `build_velocity_matrices(H_csr, grid, vel_x, vel_y)`:

- Input: full wire Hamiltonian CSR matrix, spatial_grid with coords
- Output: two CSR matrices vel_x, vel_y (same structure, scaled values)
- For each nonzero entry H(row, col):
  - spatial_row = mod(row - 1, Ngrid) + 1
  - spatial_col = mod(col - 1, Ngrid) + 1
  - vel_x%val(k) = cmplx(0, -1, dp) * (coords(1,spatial_row) - coords(1,spatial_col)) * H%val(k)
  - vel_y%val(k) = cmplx(0, -1, dp) * (coords(2,spatial_row) - coords(2,spatial_col)) * H%val(k)

Keep g='g3' mode unchanged. g='g1' and g='g2' become dead code (remove or leave commented).

### 2. `src/physics/gfactor_functions.f90`

**Modify** `pMatrixEleCalc_2d`:

- Change signature: accept `type(csr_matrix), intent(in) :: vel(3)` instead of profile_2d/kpterms_2d/cfg
- Replace the select-case with g1/g2/g3 by a simple `csr_spmv(vel(direction), state2, Y)`
- Remove dependency on ZB8bandGeneralized for g1/g2

**Modify** `gfactorCalculation_wire`:

- Change signature: accept `vel(3)` instead of profile_2d/kpterms_2d
- Pass vel through to compute_pele_2d -> pMatrixEleCalc_2d

**Modify** `compute_pele_2d`:

- Forward vel(3) to pMatrixEleCalc_2d

### 3. `src/apps/main_gfactor.f90`

**Modify** wire branch (confinement=2):

- After building HT_csr (the wire Hamiltonian at kz=0), call `build_velocity_matrices(HT_csr, cfg%grid, vel_x, vel_y)`
- Build vel_z using existing g3: `ZB8bandGeneralized(vel_z, 1.0_dp, profile_2d, kpterms_2d, cfg, g='g3')`
- Pass `vel = [vel_x, vel_y, vel_z]` to gfactorCalculation_wire
- Free vel matrices after computation

## What Stays Unchanged

- `sigmaElem_2d` (spin matrix elements)
- Lowdin second-order summation formula
- Tensor post-processing (-i * tensor / hbar2O2m0 - (g_free/2) * sigma)
- g-factor extraction formula
- Bulk and QW g-factor paths (already correct)
- g='g3' mode in ZB8bandGeneralized (correct for z-direction)

## Validation

1. Existing wire g-factor regression test should still pass (or be updated)
2. Bulk limit: for large wires, g-factor should converge to the bulk analytical value
3. Compare CB g-factor with Winkler/Roth analytical formulas for bulk
4. Test with InSbW wire at various diameters (55A, 100A, 200A, 500A)
