# 08 — FEAST-enable QW g-factor

## Parent / Source

PRD: `.scratch/architecture-deepening/PRD.md` (Thread 2, FEAST flip-on #1).
ADR 0005.

## What to build

QW g-factor currently rejects FEAST as a **temporary stopgap guard** — because
the g-factor Γ solver forces FULL mode, which routes QW+FEAST through the sparse
path and silently truncates the spectrum (FEAST is a partial-spectrum solver).
Wire g-factor keeps FEAST via a different path, so the guard is QW-only.

With the dispersion-aware window authority and ENERGY mode (the Γ solver no
longer forced to a mode that truncates under FEAST), make QW g-factor run
correctly under FEAST and **remove the temporary guard**. The "not yet supported"
framing for this combo flips to supported.

## Acceptance criteria

- [ ] QW g-factor runs under FEAST+ENERGY via the window authority; the temporary
      stopgap guard is removed.
- [ ] **FEAST-vs-dense regression for QW g-factor**: eigenvalues match within
      tolerance on the same config.
- [ ] The "not yet supported" wording for this combination is gone.
- [ ] Explicit `FEAST + INDEX` is still rejected up front (structural, permanent
      — ADR 0005).

## Blocked by

- #03 (energy-window authority) — the gate.
