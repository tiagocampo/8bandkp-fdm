# Issue 4: TRK sum rule — QW (dense) and wire (CSR) tests

**Type**: AFK
**Blocked by**: Issue 3 (needs trk_helpers.f90 and compute_trk_sum_resolvent_dense)
**User stories**: 10, 11, 12

## What to build

Extend the TRK sum rule test to QW (dense 8N×8N) and wire (CSR 8N×8N) geometries, using the infrastructure from Issue 3.

### Part 1: CSR resolvent helper

Add to `tests/support/trk_helpers.f90`:

```fortran
subroutine compute_trk_sum_resolvent_csr(H_csr, evals, evecs, vel_csr, &
  state_idx, eps_shift, trk_sum, trk_ratio)
  type(csr_matrix), intent(in) :: H_csr
  real(dp), intent(in) :: evals(:)
  complex(dp), intent(in) :: evecs(:,:)
  type(csr_matrix), intent(in) :: vel_csr
  integer, intent(in) :: state_idx
  real(dp), intent(in) :: eps_shift
  real(dp), intent(out) :: trk_sum
  real(dp), intent(out) :: trk_ratio
```

Same algorithm as the dense version but uses CSR SpMV (`csr_spmv`) and PARDISO for the linear solve. The regularization adds ε to the diagonal element corresponding to state_idx.

For PARDISO: use `pardiso_c` from `linalg.f90` (complex, mtype=13 for complex Hermitian). Setup phase, factorization, solve, cleanup.

Also add the double commutator for real-space CSR:

```fortran
function compute_double_commutator_csr(H_csr, grid, psi0) result(dc_value)
  type(csr_matrix), intent(in) :: H_csr
  type(spatial_grid), intent(in) :: grid
  complex(dp), intent(in) :: psi0(:)
  real(dp) :: dc_value
```

Computes ½·Σ_{ij} -(x_i - x_j)² · H_ij · Re[ψ₀*(i) · ψ₀(j)] by looping over CSR entries. The position x_i is extracted from the grid based on the spatial index within the band-major ordering (sp = mod(basis_idx - 1, Ngrid) + 1).

### Part 2: QW TRK test

Add to `tests/unit/test_trk_sum_rule.pf`:

Build a 3-layer Al30Ga70As/GaAs/Al30Ga70As QW (51 grid points, same as Issue 1 setup). The z-velocity matrix comes from the commutator `build_velocity_matrices_1d`.

Test procedure:
1. Build QW H at k=0, diagonalize with zheevd
2. Build CSR version via `dnscsr_z_mkl`
3. Build z-velocity via `build_velocity_matrices_1d(H_csr, grid, vel)` — vel(3) is the z-velocity
4. Set up spatial grid: `grid%ndim=1, grid%ny=ngrid, grid%z=z`
5. Compute double commutator RHS using `compute_double_commutator_csr` with grid positions
6. Compute TRK sum LHS via resolvent for CB1 (first conduction subband)
7. Assert |LHS - RHS| / max(|RHS|, 1e-15) < 1e-10
8. Report TRK ratio

Note: For QW, only vel(3) (z-direction) is the commutator. vel(1) and vel(2) are k-derivative matrices with different algebraic structure — they are NOT tested here.

### Part 3: Wire TRK test

Add to `tests/unit/test_trk_sum_rule.pf`:

Build a uniform GaAs wire (10×10 grid, 5Å spacing). The transverse velocities come from `build_velocity_matrices_2d`.

Test procedure:
1. Build wire H via `ZB8bandGeneralized` at kz=0
2. Diagonalize: convert CSR to dense, solve with zheevd (small enough for dense solve at 10×10 → 800×800)
3. Build transverse velocities via `build_velocity_matrices_2d(H_csr, grid, vel_x, vel_y)`
4. Compute double commutator RHS for x-direction and y-direction using grid coordinates
5. Compute TRK sum LHS via CSR resolvent for CB1
6. Assert |LHS - RHS| / max(|RHS|, 1e-15) < 1e-10 for both vx and vy
7. Report TRK ratios

Note: For wire, vel(3) (z-direction) is dH/dkz, NOT the commutator. Only vel_x and vel_y are commutator-based and suitable for the TRK test.

### Grid setup pattern for wire

```fortran
type(spatial_grid) :: grid
integer, allocatable :: mat2d(:,:)
allocate(mat2d(nx, ny)); mat2d = 1
call grid_init_rect(grid, nx, ny, dx, dy, mat2d)
```

The grid's `coords` array provides (x, y) positions for computing position differences in the double commutator.

## Acceptance criteria

- [ ] CSR resolvent helper added to `trk_helpers.f90`
- [ ] Double commutator CSR helper added to `trk_helpers.f90`
- [ ] QW TRK test added to `test_trk_sum_rule.pf` — z-velocity TRK at 1e-10
- [ ] Wire TRK test added to `test_trk_sum_rule.pf` — vx, vy TRK at 1e-10
- [ ] TRK ratios reported as informational diagnostics
- [ ] All existing unit tests pass (34+)
- [ ] New tests pass

## Blocked by

Issue 3 — needs `trk_helpers.f90`, `compute_trk_sum_resolvent_dense`, and the bulk TRK pattern.
