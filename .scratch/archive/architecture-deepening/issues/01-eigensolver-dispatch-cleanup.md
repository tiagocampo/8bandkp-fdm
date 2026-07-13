# 01 — Eigensolver dispatch cleanup

## Parent / Source

PRD: `.scratch/architecture-deepening/PRD.md` (Thread 1, U-A). ADR 0005.

## What to build

Two related cleanups of the eigensolver interface, end to end:

1. **Backend-aware diagnostics.** Today, non-convergence and under-delivery
   messages hardcode "FEAST" regardless of which backend ran. They must instead
   name the backend that actually ran (dense LAPACK or FEAST), read from the
   resolved config. A failure forced on each backend must produce an honest
   message.

2. **One method per matrix format.** The base type exposes three method names
   for two matrix formats, with the most-used one self-described as a "legacy
   alias." Reduce to one method per format (dense array, CSR), with the backend
   fixed at construction and never named at the call site. The format-vs-backend
   distinction (two independent axes) is already in the glossary; this makes the
   code match it.

## Acceptance criteria

- [ ] Non-convergence / under-delivery messages print the actual resolved
      backend; verifiable by forcing a failure on both backends.
- [ ] The base type exposes exactly one method per matrix format; the redundant
      alias is removed; every call site compiles and uses the canonical names.
- [ ] The existing dense-vs-sparse verification-ladder rung is extended to a
      bit-identical eigenvalue regression across all geometries.
- [ ] Golden regression outputs unchanged.

## Blocked by

None — can start immediately.
