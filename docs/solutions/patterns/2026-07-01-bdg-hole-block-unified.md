---
module: bdg_hamiltonian
tags: [SSOT, single-source-of-truth, wrapper, before/after]
problem_type: build-divergence
component: build_bdg_hole_block
---

# BdG Hole-Block Unified Convention (Before/After)

## Problem

The wire-CSR builder and the dense-QW builder both built BdG Hamiltonians with the same class-D PHS but used DIFFERENT hole-block conventions:

- Wire CSR: `hole = -transpose(H0(+k))`
- Dense QW: `hole = -conjg(H0(-k))`

Both satisfied PHS but produced DIFFERENT spectra at generic k.

## Solution (per ADR 0007)

Extract `build_bdg_hole_block(H_h, hole_block)` as the single shared wrapper. Both builders call it. Inside the wrapper:

```fortran
hole_block = -conjg(H_h)
```

## Before/after signatures

| Builder | Before | After |
|---|---|---|
| Wire CSR | `hole = -transpose(H0)` inline | `call build_bdg_hole_block(H_h, hole_block)` |
| Dense QW | `hole = -conjg(H0(-k))` inline | `call build_bdg_hole_block(H_h, hole_block)` |

`bdg_hamiltonian.f90:60-72` exports `pairing_partner(8)` and `pairing_sign(8)` as SSOT; topological_analysis.f90 imports them via `use` (was duplicated locally — see pairing_sign DRY fix in ADR 0008).