# 07 — Output writers + Simpson absorption

## Parent / Source

PRD: `.scratch/architecture-deepening/PRD.md` (Thread 4, U-C2).

## What to build

Several output formats (potential profile, optical transitions, SC diagnostics,
BdG eigenvalues) are written with near-identical code across executables, and the
Simpson base-weight rule is inlined per geometry differing only in the final
Jacobian. Give each format **one writer**, and **share the Simpson base rule**
so each geometry supplies only its Jacobian.

## Acceptance criteria

- [ ] Each output format is produced by exactly one writer; the duplicated write
      blocks are gone.
- [ ] One shared Simpson base-weight rule serves every geometry that differs only
      in Jacobian.
- [ ] Output is byte-for-byte unchanged on existing configs (golden output
      regression).

## Blocked by

None — can start immediately (parallel to everything else).
