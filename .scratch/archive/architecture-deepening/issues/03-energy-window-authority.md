# 03 — Energy-window authority (the #10 gate)

## Parent / Source

PRD: `.scratch/architecture-deepening/PRD.md` (Thread 1, #10 / FEAST gate).
ADR 0005.

## What to build

The FEAST energy window is currently derived at a single k-point (k=0) with a
fixed margin, so finite-k dispersion outruns it and FEAST silently under-delivers
mid-sweep. This is the root cause that makes FEAST untrustworthy in the remaining
consumers, and the **gate for the entire FEAST-parity thread**.

Provide **one window authority** that, given the user `[solver]` window (if set)
or else an automatically-derived window, returns `[emin, emax]`. The auto case is
a **dispersion-aware envelope**: the Gershgorin bound taken at the sweep's
endpoints (k=0 and k_max), unioned with a margin, yielding **one stable window
per sweep** — so k-to-k branch tracking stays valid (per-k moving windows are
rejected by design). Single-k consumers reduce to the one-k bound; spectral
consumers that already hold an eigenvalue array accept it as the bound source. A
first cut (max-k fallback) already landed in the eigensolver standardization;
this is the full envelope and the single authority.

The existing single-k Gershgorin estimator stays as the primitive beneath.

## Acceptance criteria

- [ ] One window authority is the sole source of `[emin, emax]` for every FEAST
      consumer; the user `[solver]` window is honored in every sweep path.
- [ ] The auto window is the union of the endpoint bounds (not a single k-point),
      and is **one stable window per sweep** (asserted, not per-k).
- [ ] Unit test covers user-override vs auto-envelope vs eigenvalue-array
      variants.
- [ ] A wide-QW / wide-k sweep that previously under-delivered mid-sweep now
      converges across the full range (regression).
- [ ] Golden regression outputs unchanged where the window was already adequate.

## Blocked by

- #02 (solver-config derivation) — the authority is adopted through the
  derivation seam.
