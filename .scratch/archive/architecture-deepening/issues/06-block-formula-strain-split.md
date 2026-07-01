# 06 — Block-formula descriptor + strain split

## Parent / Source

PRD: `.scratch/architecture-deepening/PRD.md` (Thread 3, U-B + C8). ADRs 0001,
0005. ADR 0001's "COMPLETED" wording was corrected to "DATA SSOT COMPLETE" — the
data is centralized, the **interpretation is not**; this issue finishes it.

## What to build

The block-tag → block interpretation — the select-case that turns a tag into an
actual block, including the derived difference (`Q − T`) and half-sum
(`0.5·(Q + T)`) formulas — is duplicated across the dense, COO, and QW-CSR
builders. Centralize it as **one pure resolver** returning a small descriptor
(identity / difference / half-sum, the operand tags, and the prefactor) that each
builder applies generically. **No polymorphic builder types** (ADR 0001) — only a
shared resolver plus the existing builders as thin appliers.

Separately, split the strain module along its four separable concerns (Bir-Pikus
formulas + strain table; the Navier-Cauchy PDE solver; the field-id lookup
already extracted), moving the strain-table interpretation alongside the
Bir-Pikus formulas. Account for the shared SAVE-coupling variable that currently
crosses concerns.

Also add the QW-CSR Hamiltonian builder to the architecture intent layer
(dependency graph + module guide), where it is currently missing.

## Acceptance criteria

- [ ] One resolver owns the block-tag → formula mapping; the three builders apply
      its descriptor generically; the duplicated select-cases are gone.
- [ ] No polymorphic builder types introduced (ADR 0001 respected).
- [ ] Strain module split cleanly; the cross-concern SAVE coupling resolved.
- [ ] **Bit-identical eigenvalue regression across all geometries** (extend the
      verification ladder) — eigenvalues must not move.
- [ ] QW-CSR builder documented in the intent layer (dependency graph + module
      guide).

## Blocked by

- #05 (confinement-assembly unit tests) — the C7 gate; do not refactor the
  assembly without this safety net.
