# 10 — FEAST-enable Landau

## Parent / Source

PRD: `.scratch/architecture-deepening/PRD.md` (Thread 2, FEAST flip-on #3).
ADRs 0003, 0005.

## What to build

Landau is dense-native with no sparse path, but large-`nx` Landau pays dearly for
full dense diagonalization. Wire FEAST+ENERGY for the **Landau k-solve** using
the window authority, mirroring the **existing per-thread solver pattern the QW
parallel sweep already uses** (per-thread FEAST workspace, no oversubscription,
no cross-thread cache sharing).

The **B-sweep and fan diagram stay in the application layer** (ADR 0003); only
the k-solve backend changes.

## Acceptance criteria

- [ ] Landau k-solve runs under FEAST+ENERGY via the window authority.
- [ ] Per-thread FEAST workspace follows the existing QW parallel-sweep pattern;
      no oversubscription, no cross-thread cache sharing.
- [ ] The B-sweep and fan diagram remain application-layer (ADR 0003 respected).
- [ ] **FEAST-vs-dense regression for Landau**: eigenvalues match within
      tolerance.
- [ ] CTest parallel-run guidance (cap OMP to `nproc/N`) still passes green.

## Blocked by

- #03 (energy-window authority) — the gate.
