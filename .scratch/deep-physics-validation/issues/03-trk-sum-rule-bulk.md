# Issue 3: TRK sum rule — shared helper + bulk test

**Type**: AFK
**Blocked by**: None — can start immediately
**User stories**: 9, 12, 13

## What to build

Create the shared TRK sum rule infrastructure and implement the bulk (8×8) test. This issue establishes the pattern that issues 4 (QW + wire) will follow.

### Part 1: Double commutator helper module

Create `tests/support/trk_helpers.f90` with two public subroutines:

#### `compute_double_commutator_bulk`

For the bulk 8×8 Hamiltonian in k-space, the double commutator equals d²H/dk². Compute it as:

```fortran
subroutine compute_double_commutator_bulk(d2Hdk2, params, direction)
  complex(dp), intent(out) :: d2Hdk2(8, 8)
  type(paramStruct), intent(in) :: params(1)
  integer, intent(in) :: direction  ! 1=x, 2=y, 3=z
```

This is the second k-derivative of `ZB8bandBulk(k)`. At k=0, the k·p Hamiltonian is:
- k² terms → contribute constant to d²H/dk² (these are the free-electron ℏ²/(2m₀) terms and remote band corrections)
- k terms → contribute zero at k=0 (first derivative vanishes)
- k⁰ terms → contribute zero (no k dependence)

The d²H/dk² can be computed analytically from the k·p block table: for each entry with a k-linear term (KP_PP, KP_PM, KP_PZ, KP_R, KP_RC, KP_S, KP_SC), the second derivative is zero. For k-quadratic terms (KP_Q, KP_T, KP_A), the second derivative is 2× the coefficient. The constant terms (band offsets) have zero second derivative.

Alternatively, compute numerically: evaluate H(k) at k=+δ and k=-δ, then d²H/dk² ≈ (H(+δ) - 2H(0) + H(-δ)) / δ². Use δ = 1e-4 Å⁻¹ for machine-precision accuracy.

Use the numerical approach for simplicity — it avoids deriving the analytical formula and is robust against future changes to the k.p structure.

#### `compute_trk_sum_resolvent`

Generic resolvent-based TRK computation. For dense matrices:

```fortran
subroutine compute_trk_sum_resolvent_dense(H, evals, evecs, vel, state_idx, &
  eps_shift, trk_sum, trk_ratio)
  complex(dp), intent(in) :: H(:,:)          ! Hamiltonian
  real(dp), intent(in) :: evals(:)            ! Eigenvalues
  complex(dp), intent(in) :: evecs(:,:)       ! Eigenvectors (columns)
  complex(dp), intent(in) :: vel(:,:)         ! Velocity matrix (dH/dk or commutator)
  integer, intent(in) :: state_idx            ! Initial state (1-based)
  real(dp), intent(in) :: eps_shift           ! Regularization shift
  real(dp), intent(out) :: trk_sum            ! TRK sum value
  real(dp), intent(out) :: trk_ratio          ! trk_sum / hbar2O2m0
```

Implementation:
1. b = vel · ψ₀ (matrix-vector product, where ψ₀ = evecs(:, state_idx))
2. Construct A = H - E₀·I + ε·(ψ₀·ψ₀ᴴ) (regularized system matrix)
3. Solve A·y = b using LAPACK `zgesv`
4. c = vel · y
5. trk_sum = Re[ψ₀ᴴ · c]
6. trk_ratio = trk_sum / hbar2O2m0

### Part 2: Bulk TRK test

Create `tests/unit/test_trk_sum_rule.pf` with the bulk TRK test.

For bulk at k=0, the Hamiltonian is 8×8 diagonal (GaAs: all k-linear and k-quadratic terms vanish at k=0 for a T_d symmetry material). The velocity matrices v_α = dH/dk_α at k=0 are the 8×8 Kane coupling matrices.

Test procedure:
1. Build GaAs bulk H(0) via `ZB8bandBulk` with k = (0,0,0)
2. Diagonalize with `zheevd` to get 8 eigenvalues and eigenvectors
3. Build velocity matrix vx = dH/dkx at k=0 via `ZB8bandBulk` with `g='g'` and unit kx
4. Compute double commutator RHS numerically: d²H/dkx² = (H(δ,0,0) - 2H(0,0,0) + H(-δ,0,0)) / δ²
5. Compute TRK sum LHS via resolvent for initial state = CB (band 7, highest eigenvalue)
6. Assert |LHS - RHS| / max(|RHS|, 1e-15) < 1e-10
7. Report trk_ratio as informational diagnostic

Repeat for initial state = HH (band 1, lowest valence).

Repeat for all 3 velocity directions (kx, ky, kz) — should give identical results by cubic symmetry.

Also test with InAs (different Kane parameters, larger Ep) to verify material independence.

### CMake integration

Add `trk_helpers.f90` to `tests/support/CMakeLists.txt` as a source in the pFUnit support library (alongside existing `csr_test_helpers.f90`, `krylov_helpers.f90`, etc.).

Add `test_trk_sum_rule.pf` to the pFUnit test generation in `tests/unit/CMakeLists.txt`.

## Acceptance criteria

- [ ] `tests/support/trk_helpers.f90` created with `compute_double_commutator_bulk` and `compute_trk_sum_resolvent_dense`
- [ ] `tests/unit/test_trk_sum_rule.pf` created with bulk TRK test
- [ ] Bulk TRK: LHS (resolvent) agrees with RHS (double commutator) at 1e-10 for GaAs and InAs
- [ ] All 3 velocity directions tested (kx, ky, kz)
- [ ] Initial states tested: CB (band 7) and HH (band 1)
- [ ] TRK ratio reported as informational diagnostic (not pass/fail)
- [ ] CMake updated to compile `trk_helpers.f90` and `test_trk_sum_rule.pf`
- [ ] All existing unit tests pass (34+)
- [ ] New tests pass

## Blocked by

None — can start immediately.
