# Issue 1: QW Hermiticity tests with block-level diagnostics

**Type**: AFK
**Blocked by**: None — can start immediately
**User stories**: 1, 2, 3, 7, 8

## What to build

Add three Hermiticity tests for the QW dense Hamiltonian in `tests/unit/test_hamiltonian.pf`, plus a block-level diagnostic helper. Currently, `ZB8bandQW` (the 8N×8N dense QW builder used by `bandStructure`, `gfactorCalculation`, and `opticalProperties`) has no direct H = H† test — it is only validated indirectly through LAPACK diagonalization success.

### Test 1: QW base Hermiticity (H1)

Build a 3-layer Al30Ga70As/GaAs/Al30Ga70As QW using the existing `setup_algaas_qw` pattern (or replicate the setup from `test_qw_variable_coeff_order4`). Construct the dense Hamiltonian via `ZB8bandQW`, then assert `max|H(i,j) - conjg(H(j,i))| < 1e-12` over all i,j. This is the same pattern as `test_bulk_hermitian` but for the 8N×8N dense matrix.

### Test 2: QW + strain Hermiticity (H2)

Same QW setup but enable Bir-Pikus strain. Construct a `simulation_config` with `cfg%strain%enabled = .true.` and appropriate `cfg%strain%reference` lattice constant. Pass `cfg` to `ZB8bandQW`. Assert Hermiticity at 1e-12. The strain insertion uses `apply_strain_table_dense` which adds Bir-Pikus shifts to the diagonal blocks — this test verifies those shifts preserve Hermiticity.

### Test 3: QW + Zeeman Hermiticity (H3)

Same QW setup but enable Zeeman splitting. Construct a `simulation_config` with `cfg%b_field%enabled = .true.`, `cfg%b_field%B_vec = [0.0, 0.0, 1.0]`, `cfg%b_field%g_factor = 2.0`. Pass `cfg` to `ZB8bandQW`. Assert Hermiticity at 1e-12.

### Block-level diagnostic helper

Add a subroutine `assert_hermitian_with_block_diagnostics` (either in the test file or in `tests/support/`) that:
1. Takes a dense 8N×8N matrix and grid size N
2. Checks overall Hermiticity: max|H(i,j) - conjg(H(j,i))|
3. On failure, extracts the 8×8 block structure (block(i,j) = rows [(i-1)*N+1:i*N], cols [(j-1)*N+1:j*N])
4. Reports which block pair (i,j) vs (j,i) has the largest violation, the absolute error, and the k.p terms responsible (band indices map to: 1-4 = valence HH/LH, 5-6 = split-off, 7-8 = conduction)
5. Returns the max error for the assertion

This helper is used by all three QW Hermiticity tests. The band structure is: bands 1-4 = valence (HH, LH, LH, HH), bands 5-6 = split-off, bands 7-8 = conduction.

### Setup code pattern

Follow the existing QW setup from `test_qw_variable_coeff_order4` in `test_hamiltonian.pf`:

```fortran
! Grid
integer, parameter :: ngrid = 51
real(dp), parameter :: dz_val = 2.0_dp
real(dp) :: z(ngrid)
! Layers
integer :: int_start_pos(3), int_end_pos(3)
character(len=20) :: material(3)
type(paramStruct) :: params(3)
! kpterms + profile
real(dp), allocatable :: profile(:,:), kpterms(:,:,:)
! Hamiltonian
integer :: matSize
complex(dp), allocatable :: HT(:,:)

do i = 1, ngrid
  z(i) = -50.0_dp + dble(i-1) * dz_val
end do
int_start_pos = [1, 16, 37]; int_end_pos = [15, 36, 51]
material = ["Al30Ga70As", "GaAs       ", "Al30Ga70As"]
call paramDatabase(material, 3, params)
call confinementInitialization_raw(z, int_start_pos, int_end_pos, material, &
  3, params, 'z', profile, kpterms, FDorder=2)
matSize = ngrid * 8
allocate(HT(matSize, matSize))
```

For strain: set `cfg%strain%enabled = .true.` and set a reference lattice constant.
For Zeeman: set `cfg%b_field%enabled = .true.`, `cfg%b_field%B_vec = [0,0,1]`, `cfg%b_field%g_factor = 2.0`.

## Acceptance criteria

- [ ] `test_qw_hermitian` added to `test_hamiltonian.pf` — base QW H = H† at 1e-12
- [ ] `test_qw_strain_hermitian` added to `test_hamiltonian.pf` — strained QW H = H† at 1e-12
- [ ] `test_qw_zeeman_hermitian` added to `test_hamiltonian.pf` — Zeeman QW H = H† at 1e-12
- [ ] Block-level diagnostic helper implemented and used by all three tests
- [ ] On Hermiticity failure, the helper reports: "Block (i,j) vs (j,i) error = X.XXe-XX" with band names
- [ ] All existing unit tests pass (34+)
- [ ] New tests pass

## Blocked by

None — can start immediately.
