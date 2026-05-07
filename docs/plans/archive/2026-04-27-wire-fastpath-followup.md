# Wire Fast-Path Follow-Up Refactoring

Deferred from PR #11 review (feature/hamiltonian-performance-refactor).
All items verified against current code; safe to tackle in any order.

## High Priority

### 1. Deduplicate block insertion between fast/slow paths (~260 lines)

`hamiltonian_wire.f90` lines 282-478 (fast) and 570-784 (slow) contain
near-identical 8x8 block COO insertion code. The g3-mode insertion is also
duplicated (lines 255-340 and 543-633). The only difference is that fast path
uses workspace-associated arrays while slow allocates its own.

**Approach:** Extract shared insertion subroutines that take COO array
arguments (rows, cols, vals, coo_idx). Both paths call the same insertion
routine with different COO backing stores.

### 2. Pass `coo_cache` to kz-sweep callers

Neither `main.f90` nor `main_optics.f90` passes `coo_cache` alongside
`wire_ws` to `ZB8bandGeneralized`. The COO cache fast path (avoiding
O(NNZ log NNZ) sort per kz-point) is therefore never active. The sort is
likely the dominant per-kz cost for large wires.

**Approach:** Either (a) absorb `coo_cache` into `wire_workspace` so callers
only manage one object, or (b) thread `coo_cache` through the existing caller
signatures. Option (a) is cleaner — the workspace already owns the COO
arrays.

## Medium Priority

### 3. Pre-allocate `HT_csr_step` across kz-points

`main.f90` lines 393-397 does `csr_clone_structure` + `csr_free` every
kz-point. Since sparsity is kz-independent, this could be a persistent CSR
in the workspace or a module-level variable.

### 4. Cache conjugate-transpose scatter maps

`build_kp_term_SC/RC/PM` fast paths call
`csr_conjugate_transpose_to_preallocated` which does per-element binary
search for transpose positions. The transpose structure is fixed across
kz-points. Pre-compute transpose scatter maps during slow-path initialization.

### 5. Extract COO finalization into shared helper

The COO-to-CSR finalization block (lines 491-511 and 803-823) is identical
between fast and slow paths. Extract into `finalize_coo_to_csr` helper.

### 6. Extract diagonal-position finder into `sparse_matrices` utility

The same 13-line diagonal-position lookup is copy-pasted 3x in
`zb8_generalized_slow` (lines 843-888) and also needed in `eigensolver.f90:661`.
Extract to `csr_find_diag_positions(mat, diag_pos)` in `sparse_matrices.f90`.

## Low Priority

### 7. Cache `build_strain_table()` result

`insert_strain_coo` calls `build_strain_table()` every kz-point. The table
is a 32-entry compile-time constant. Store in module-level `save` variable
or in `wire_workspace`.

### 8. Move general CSR utilities out of physics modules

`csr_conjugate_transpose` and `csr_conjugate_transpose_to_preallocated`
(`hamiltonian_wire.f90:1362-1416`) and `build_diagonal_csr`
(`confinement_init.f90:660-678`) are general sparse linear algebra routines
with no physics logic. Move to `sparse_matrices.f90`.
