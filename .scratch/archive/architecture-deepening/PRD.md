**Status**: COMPLETE (2026-07-05)

# PRD — Architecture Deepening + FEAST Parity

Branch: `refactor/architecture-deepening` (off `main`).
Source: 2026-06-13/14 architecture-review grill session. Decisions recorded in
**ADR 0005** (FEAST parity), ADR 0001 (scope correction), ADR 0002 (validation),
ADR 0003 (Landau sweep ownership), ADR 0004 (method-aware AUTO). Domain vocabulary
from `CONTEXT.md` (*AUTO*, *requested band window*, *eigensolver dispatch —
format vs backend*).

---

## Problem Statement

The eigensolver-standardization effort unified diagonalization behind a
polymorphic interface and a method-aware `AUTO` (ADR 0004), but the sparse
(FEAST) backend only genuinely works for **some** geometries and apps. For the
others it is either silently overridden to the dense LAPACK backend or rejected
up front — so a physicist who picks the sparse solver on a large QW, a Landau
grid, or a QW-optics accumulation either gets dense LAPACK performance anyway
(silent override, violating the project's "no silent corrections" rule) or hits
an opaque abort deep in the solver.

Three forms of this gap exist today:

1. **QW optics hardcodes the dense backend** and ignores `[solver].method`, while
   wire optics honors it — one executable, two behaviors.
2. **QW g-factor rejects FEAST** as a hard guard (a stopgap for a real
   truncation hazard), and **Landau is dense-native** with no sparse path.
3. **The FEAST energy window is derived at a single k-point** (k=0, with a fixed
   margin), so finite-k dispersion outruns it and FEAST silently under-delivers
   eigenvalues mid-sweep — the root cause that makes FEAST untrustworthy in the
   remaining consumers.

Layered on top is accumulated architectural debt that makes the above hard to
fix cleanly: the solver-config assembly is hand-written per confinement, the
wire init+cleanup is copy-pasted four times (and silently skips strain in the
topology wire subroutines — a latent physics bug), the k.p block *interpretation*
is duplicated across three builders (physics-critical), and the most-used
eigensolver entry point is a self-described "legacy alias" whose error messages
always say "FEAST" regardless of the backend that ran.

## Solution

Four threads, sequenced so the physics-critical and enabling work lands first:

**Thread 1 — Solver-config & window (the FEAST-parity enabler).** Consolidate
the per-confinement eigensolver-config assembly behind one derivation, and make
the energy window **dispersion-aware** (an envelope over the sweep, not a single
k-point). The window is the gate: once FEAST converges reliably across a full
sweep, it can be trusted in every consumer.

**Thread 2 — FEAST parity (ADR 0005).** With a trustworthy window, flip FEAST on
for each remaining consumer in turn — QW g-factor, then QW optics, then Landau —
each behind a FEAST-vs-dense regression test. The only *permanent* FEAST
rejection is the structurally impossible `FEAST + INDEX`. Bulk never offers
FEAST (always 8×8, dense by nature).

**Thread 3 — Wire setup & the strain bug.** Extract the wire init+cleanup into
one type. Doing so routes the topology wire subroutines through the same setup
path as everything else — which applies strain, fixing a latent defect where
strained wire topology/BdG calculations silently ignored strain.

**Thread 4 — Block-table interpretation & output writers.** Unit-test the k.p
block assembly directly (the gate), then centralize the duplicated
block-formula interpretation. Separately, absorb the copy-pasted output writers.

The net is one solver seam, FEAST as a first-class backend everywhere it's
physical, no latent strain bug, fewer misleading errors, and a testable
foundation under the most physics-critical code.

## User Stories

### FEAST parity (the headline)
1. As a physicist running a wide QW, I want to use the sparse FEAST backend so
   the diagonalization does not scale as N² memory and N³ time.
2. As a physicist, I want FEAST to converge reliably across an entire k-sweep,
   so it does not silently drop eigenvalues partway through.
3. As a physicist, I want the `[solver]` energy window I set to be honored in
   every sweep path, not just some of them.
4. As a physicist running QW optics, I want `[solver].method = FEAST` honored
   rather than silently replaced by the dense backend.
5. As a physicist running a large-`nx` Landau grid, I want the sparse backend
   available, while the B-sweep and fan diagram remain in the application layer.
6. As a physicist, I want an explicit "not yet supported" message — never a
   silent override and never an abort deep in the eigensolver — when a
   combination is genuinely unavailable.
7. As a physicist, I want an explicit `FEAST + INDEX` request rejected up front,
   because FEAST owns an energy window, not an index range.
8. As a physicist running QW g-factor, I want FEAST as an option eventually
   (today it is a temporary stopgap rejection).

### Solver-config & window
9. As a maintainer, I want the `AUTO` resolution unit-testable per geometry
   directly, so a default-resolution regression is caught without running a full
   executable.
10. As a maintainer, I want one derivation that builds the eigensolver config, so
    the four confinement dispatch blocks cannot drift apart.
11. As a maintainer, I want the energy-window logic owned by one authority, so
    the dispersion-aware envelope is not re-implemented per consumer.
12. As a maintainer, I want a single-k consumer (g-factor at Γ, optics at fixed
    k) and a sweep consumer to use the same window authority with different
    k-contexts.
13. As a maintainer, I want one stable energy window per sweep (not a per-k
    moving window), so branch-tracking across k-points stays valid.

### Wire setup & strain
14. As a physicist running a strained wire topology or BdG calculation, I want
    strain actually applied (today the topology wire subroutines skip it
    silently).
15. As a maintainer, I want wire initialization and cleanup in one type, so a
    cleanup or strain bug is fixed once, not in four copies.
16. As a maintainer, I want sink-style callers (those that receive pre-computed
    confinement data as arguments) to accept that data rather than re-initialize.

### Block-table interpretation & tests
17. As a maintainer, I want the k.p block assembly unit-tested directly, so a
    sign error in a derived term is caught at the source — not via a shifted
    eigenvalue that still passes a loose tolerance.
18. As a maintainer, I want the block-tag → block-formula interpretation in one
    place, so adding or fixing a derived term touches one module, not three
    builders.
19. As a maintainer, I want no polymorphic builder types introduced (per ADR
    0001) while still centralizing the interpretation.

### Eigensolver interface & messages
20. As a developer debugging a non-convergence, I want the error message to name
    the backend that actually ran (FEAST vs dense LAPACK), not always say
    "FEAST."
21. As a maintainer, I want the eigensolver interface to expose one method per
    matrix format (dense array, CSR), with the backend hidden at construction —
    not a redundant "legacy alias" for one format.
22. As a maintainer, I want the format-vs-backend distinction (two independent
    axes) documented in the glossary so it stops being misdescribed.

### Output writers
23. As a maintainer, I want each output format (potential profile, optical
    transitions, SC diagnostics, BdG eigenvalues) written by one writer, so a
    format change hits one place.
24. As a maintainer, I want the Simpson base-weight rule shared across the
    geometries that differ only in their final Jacobian.

### Cross-cutting & agent-readiness
25. As a maintainer, I want every solver-touching change guarded by a
    bit-identical eigenvalue regression across geometries.
26. As an agent executing an issue, I want each user story to be an
    independently-grabbable vertical slice with a clear test.
27. As a maintainer, I want the existing verification-ladder rung that compares
    dense vs sparse eigenvalues extended to cover every newly FEAST-enabled
    consumer.
28. As a physicist, I want the band window I request (`num_cb`/`num_vb`) to be
    honored by whichever backend runs, with the backend an implementation
    detail.
29. As a maintainer, I want validation messages to distinguish a *temporary*
    "not yet supported" gap from a *permanent* "structurally impossible"
    rejection.
30. As a maintainer, I want the newly added Hamiltonian builder module
    documented in the architecture intent layer (it is currently absent from the
    dependency graph and module guide).

## Implementation Decisions

### Deep modules (simple interface, much function, rarely changes)

- **Solver-config derivation.** One pure entry point that, given the parsed
  config and the resolved geometry plus the requested band window, returns a
  fully-built, validated eigensolver config and the constructed solver. It
  *calls* the existing method-aware `AUTO` resolver (ADR 0004) rather than
  replacing it — two layers, not two functions racing to own `AUTO`. It owns the
  `nev`/`il`/`iu` bookkeeping, the subspace-size propagation, validation, and
  construction. The four confinement dispatch blocks and the apps cross this one
  seam instead of hand-assembling.

- **Energy-window authority.** One entry point that, given the user `[solver]`
  window (if set) or else an automatically-derived window, returns the
  `[emin, emax]` to use. The auto case is a **dispersion-aware envelope**: the
  Gershgorin bound taken at the sweep's endpoints and unioned with a margin,
  yielding one stable window per sweep (so k-to-k branch tracking remains
  valid). For single-k consumers it reduces to the one-k bound; for spectral
  consumers that already hold an eigenvalue array, it accepts that array as the
  bound source. The existing single-k Gershgorin estimator stays as the
  primitive beneath it. This is the gate for all of Thread 2.

- **Wire setup type.** A type owning the 2D confinement profile, the per-element
  sparse k.p terms, the wire workspace and COO cache, and the eigensolver — with
  an init and an idempotent free. Crucially, init runs the **same strain and
  electric-field setup** the QW path gets, which closes the strain-omission bug
  by construction. A variant accepts pre-computed profile/terms for sink-style
  callers.

- **Block-formula descriptor.** A pure resolver in the block-table module that
  turns a block tag into a small descriptor (identity / difference /
  half-sum, the operand block tags, and the prefactor). The three builders each
  shrink to one generic apply of the descriptor. **No polymorphic builder types**
  (respects ADR 0001) — only a shared resolver plus the existing builders.

### Non-module modifications

- The eigensolver base type drops the redundant "legacy alias" method, leaving
  one method per matrix format (dense array, CSR); the backend stays hidden at
  construction. The glossary records that format and backend are two independent
  axes.
- Error paths that today hardcode "FEAST" print the backend that actually ran.
- The QW-optics path and the three topology wire subroutines and the spectral
  function path cross the new seams instead of building configs / initializing
  confinement by hand.
- Validation re-words the temporary FEAST gaps as "not yet supported" and keeps
  the one permanent rejection (`FEAST + INDEX`) phrased as structurally
  impossible. The single-k serial solver (the g-factor Γ path) is renamed to
  reflect its true role, so it is not mistaken for a general k-sweep entry
  point.
- The newly present QW-CSR Hamiltonian builder is added to the architecture
  intent layer (dependency graph + module guide), where it is currently missing.

### Architectural decisions respected

- ADR 0001: fat setup type, `select case` dispatch, no polymorphic builder
  types. The block-formula descriptor extends the block-table work without
  introducing polymorphism.
- ADR 0002: all config-level checks via `validate` / `validate_semantic`, no
  silent corrections; the "not yet supported" messages extend this.
- ADR 0003: the Landau B-sweep and fan diagram stay in the application layer;
  only the k-solve backend changes when Landau gains FEAST.
- ADR 0004: `AUTO` stays method-aware and owned by its existing resolver; the
  new derivation calls it.
- ADR 0005: FEAST parity everywhere except bulk; the window envelope; the
  permanent-vs-temporary rejection distinction.

## Testing Decisions

**What makes a good test here.** Tests assert *external behavior* (eigenvalues,
convergence counts, whether strain shifted the spectrum, whether a config is
accepted/rejected) — never internal call structure. The bar is
**bit-identical eigenvalues** (within published tolerance) when a backend or
assembly path changes, because the physics must not move.

**Modules tested and how.**

- *Solver-config derivation* — new unit test: `AUTO` resolution per geometry ×
  method × mode, including the no-invalid-combination invariant (mirrors the
  existing table-driven eigensolver unit test that ADR 0004 added).
- *Energy-window authority* — new unit test: user-window override vs auto-envelope
  vs eigenvalue-array variant; asserts one stable window per sweep.
- *Wire setup type* — **regression test**: a strained wire topology config whose
  results shift once strain is applied (proves the bug fix). Mandatory because it
  is a latent physics defect, not a refactor.
- *Block-formula descriptor* — the **C7 gate**: a new confinement-assembly unit
  test plus direct block-tag → formula tests (difference and half-sum identities,
  Hermiticity, FD-order consistency). Must land *before* the descriptor refactor.
- *Eigensolver interface trim* — extend the existing dense-vs-sparse
  verification-ladder rung into a bit-identical eigenvalue regression across all
  geometries.
- *FEAST flip-on* — each newly-enabled consumer (QW g-factor, QW optics, Landau)
  gets its own FEAST-vs-dense regression test as its enablement gate (ADR 0005).

**Prior art in the codebase.** The table-driven eigensolver unit test
(`confinement × method × mode`), the verification ladder's dense-vs-sparse rung,
the standard-star benchmark sweeps, the golden-output regression tests, and the
validation rejection integration tests (config accepted/rejected via exit codes).

## Out of Scope

- **Bulk + FEAST** — never; bulk is 8×8, dense by nature. Rejecting it is
  permanent, not a wiring gap.
- **Per-k moving energy windows for sweeps** — rejected by design (breaks
  branch tracking); the envelope is the chosen strategy (ADR 0005).
- **Piezoelectric** — zero by symmetry for this material/orientation, excluded
  project-wide.
- **Extended cross-code validation (kdotpy)** — a separate backlog phase, blocked
  on upstream API capabilities; unrelated to this effort.
- **Reusable BdG COO workspace** and **standalone dense↔CSR converter cleanup** —
  deferred (the former behind profiling; the latter folds into the eigensolver
  work if convenient).
- **FEAST enablement beyond the flip-on order** — e.g., further consumers are
  out of scope until the window proves robust on the first three.
- **Re-litigating ADR 0001's fat-type decision** — this extends, never re-opens,
  the no-polymorphic-builders choice.

## Further Notes

- **Currency.** A first cut of the energy window (max-k auto-window fallback) and
  the QW+FEAST+g-factor stopgap guard already landed in the eigensolver
  standardization merge. The window authority above is the *full* envelope
  (union of both endpoints), of which the landed code is a partial step; the
  g-factor guard remains a temporary stopgap until QW g-factor is FEAST-enabled.
- **Hard dependency.** The block-formula descriptor (Thread 4) is **gated on the
  C7 confinement-assembly tests** — do not refactor the k.p block assembly
  without direct unit coverage first.
- **Latent bug.** The strain-omission in the topology wire subroutines is a real
  physics defect (confirmed in the grill session), not a style issue. The wire
  setup type fixes it by routing through the strain-aware path; carry a
  regression test and a documented solution entry once fixed.
- **Sequencing.** Solver-config derivation + window authority first (they enable
  everything); the message fix is unblocked and can land anytime; wire setup +
  strain fix next; C7 tests as the gate; then the block-formula descriptor
  paired with the strain-module split. The output-writer absorption runs in
  parallel throughout.
- **Intent-layer debt.** The QW-CSR Hamiltonian builder is present in the source
  but absent from the dependency graph and module guide — a documentation fix
  that should ride along with the block-table work it most affects.
