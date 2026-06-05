# Issue 2: Wire Hermiticity tests with strain, Zeeman, and g3 mode

**Type**: AFK
**Blocked by**: None — can start immediately
**User stories**: 4, 5, 6

## What to build

Add three Hermiticity tests for the wire CSR Hamiltonian in `tests/unit/test_hamiltonian_2d.pf`. The existing wire Hermiticity tests cover only the base Hamiltonian (uniform GaAs) and a heterostructure (GaAs/InAs). No test covers strain, Zeeman, or g3 mode. These physics insertions modify the CSR structure and could break Hermiticity.

### Test 1: Wire + strain Hermiticity (H4)

Build a GaAs/InAs heterostructure wire using the existing pattern from `test_heterostructure_wire_hamiltonian_is_hermitian`. Enable strain by setting `cfg%strain%enabled = .true.` and `cfg%strain%reference = "GaAs"` (use the substrate lattice constant). The strain insertion uses `insert_strain_coo` which adds Bir-Pikus shifts via COO before finalization to CSR. Assert Hermiticity via `csr_hermitian_error(HT_csr) < 1e-10`.

Note: The strain solver (`compute_strain_wire`) solves a Navier-Cauchy PDE. For the Hermiticity test, we only need the Bir-Pikus block insertion, not the full strain field. Two options:
(a) Call the full strain pipeline and check Hermiticity of the final Hamiltonian
(b) Manually set a uniform strain tensor and check Hermiticity of the Bir-Pikus insertion only

Option (a) is more thorough but requires the full PDE solve. Option (b) is faster and isolates the Bir-Pikus algebra. Use option (a) for a realistic test.

### Test 2: Wire + Zeeman Hermiticity (H5)

Build a uniform GaAs wire. Enable Zeeman by passing `B_vec` and `g_factor` to `ZB8bandGeneralized` via the config:
```fortran
cfg%b_field%enabled = .true.
cfg%b_field%B_vec = [0.0_dp, 0.5_dp, 0.0_dp]  ! Bx = 0.5 T
cfg%b_field%g_factor = 2.0_dp
```
The Zeeman insertion adds diagonal shifts via `get_zeeman_table()`. Assert Hermiticity via `csr_hermitian_error(HT_csr) < 1e-10`.

### Test 3: Wire g3 mode Hermiticity (H6)

Build a uniform GaAs wire. Call `ZB8bandGeneralized` with `g='g3'` to build the dH/dkz matrix. This produces a CSR matrix used for the wire z-velocity (not the commutator — the k-derivative). Assert Hermiticity via `csr_hermitian_error(HT_csr) < 1e-10`.

Note: The g3 mode uses a separate workspace cache (`wire_workspace` with `g3` flag). The existing test `test_workspace_g3_does_not_corrupt_normal_cache` verifies that g3 doesn't corrupt the normal cache but does NOT check g3 matrix Hermiticity.

### Setup code pattern

Follow the existing wire setup from `test_heterostructure_wire_hamiltonian_is_hermitian` in `test_hamiltonian_2d.pf`:

```fortran
integer, parameter :: nx = 10, ny = 10
real(dp), parameter :: dx = 5.0_dp, dy = 5.0_dp
type(spatial_grid) :: grid
type(csr_matrix) :: HT_csr
type(csr_matrix), allocatable :: kpterms_2d(:)
real(dp), allocatable :: profile_2d(:,:)
integer, allocatable :: mat2d(:,:)

! Grid with material map
allocate(mat2d(nx, ny))
mat2d = 1  ! GaAs everywhere (or split for heterostructure)

call grid_init_rect(grid, nx, ny, dx, dy, mat2d)
! Materials and kpterms setup
! ... (follow existing pattern)
call ZB8bandGeneralized(HT_csr, kz_value, profile_2d, kpterms_2d, cfg)
```

Use `csr_hermitian_error(HT_csr)` from `csr_test_helpers` for the assertion.

### CTest registration

These are pFUnit unit tests — they are automatically registered by the existing CMake pFUnit infrastructure. No manual CTest registration needed.

## Acceptance criteria

- [ ] `test_wire_strain_hermitian` added to `test_hamiltonian_2d.pf` — strained wire CSR H = H† at 1e-10
- [ ] `test_wire_zeeman_hermitian` added to `test_hamiltonian_2d.pf` — Zeeman wire CSR H = H† at 1e-10
- [ ] `test_wire_g3_hermitian` added to `test_hamiltonian_2d.pf` — g3 mode CSR H = H† at 1e-10
- [ ] All existing unit tests pass (34+)
- [ ] New tests pass

## Blocked by

None — can start immediately.
