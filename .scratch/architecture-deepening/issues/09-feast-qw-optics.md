# 09 — FEAST-enable QW optics

## Parent / Source

PRD: `.scratch/architecture-deepening/PRD.md` (Thread 2, FEAST flip-on #2 / #8).
ADR 0005.

## What to build

QW optics today **hardcodes the dense backend and ignores `[solver].method`** (a
silent-correction violation of the project's validation rule), while wire optics
honors it — one executable, two behaviors. Route QW optics through the
solver-config derivation so it honors `[solver]`, and make FEAST work for the
optics accumulation via the window authority + ENERGY mode. This kills the silent
override and makes the sparse backend available for wide-QW optics where dense
LAPACK is slow.

## Acceptance criteria

- [ ] QW optics honors `[solver].method` (AUTO/DENSE/FEAST); the silent DENSE
      override is gone.
- [ ] **FEAST-vs-dense regression for QW optics**: spectra match within
      tolerance.
- [ ] An explicit `FEAST + INDEX` request on QW optics is still rejected up front
      (structural).
- [ ] Golden optics regression outputs unchanged for the dense default.

## Blocked by

- #02 (solver-config derivation) — QW optics crosses the seam.
- #03 (energy-window authority) — the gate.
