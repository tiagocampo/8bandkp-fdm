# ADR 0005: FEAST parity everywhere except bulk

## Status

Accepted

## Context

The eigensolver-standardization refactor unified every diagonalization behind a
polymorphic interface (`eigensolver_base` with `solve_dense`/`solve_sparse`) and a
`[solver]` TOML section with `method`/`mode` dispatch (ADR 0004). FEAST (sparse,
energy-window-based) currently works for: **wire** in every app, and the **QW
band-structure sweep** (via `ENERGY` mode — see ADR 0004). It is rejected or
unwired for: QW g-factor (validate rule, the "#2" rejection), QW optics
(`main_optics.f90` hardcodes `DENSE`), and Landau (DENSE-native per
`defs.f90:1030`).

These read as "FEAST is invalid here," but that framing is wrong. The only
*structurally impossible* combination is **explicit `FEAST + INDEX`** (validate
check I15, ADR 0002/0004): FEAST has no index-range (`il:iu`) interface, only an
energy-window interface. Every other "FEAST rejected" case is really **"not yet
wired"** — a temporary gap, not a design prohibition. The silent DENSE override
in QW optics (`main_optics.f90:153,322`) is, worse, a violation of the project's
"no silent corrections" rule (ADR 0002): one executable honors `[solver].method`
for wire optics and silently ignores it for QW optics.

The actual blocker for wiring FEAST into the remaining consumers is the
**energy-window derivation**. `auto_compute_energy_window`
(`eigensolver.f90:411`) computes Gershgorin bounds at a single k-point (k=0,
`main.f90:700`) plus a fixed margin (`0.1·range + 0.5 eV`) and reuses that one
window for an entire sweep. Finite-k dispersion outruns it, so FEAST finds too
few eigenvalues and reports spurious non-convergence — the "#10" defect. Until
the window is dispersion-aware, FEAST cannot be trusted in the remaining
consumers.

The QW band-sweep is the existence proof that a "requested band window"
consumer (the `il:iu` concept, `CONTEXT.md`) *can* be served by FEAST via
`ENERGY` mode with an appropriate window. g-factor, optics, and Landau have
simply not adopted that translation yet.

## Decision

**FEAST is a first-class, fully-working solver for every consumer × geometry
except bulk** (bulk is always 8×8, dense by nature — sparse iteration has no
benefit and is never offered).

1. **The only permanent FEAST rejection is explicit `FEAST + INDEX` (I15).** It
   is structural: FEAST owns an energy window, not an index range. Everything
   that resolves to `FEAST + ENERGY` (the `AUTO` resolution for FEAST, ADR 0004)
   is valid in principle.

2. **All other FEAST gaps are temporary "not yet implemented," not "invalid."**
   QW g-factor, QW optics, and Landau+FEAST flip on one at a time, each gated on
   the window (below) proving robust, each behind a regression test. Their
   validate-rejection messages must say "not yet implemented," never "invalid,"
   so a future maintainer does not mistake a wiring gap for a design prohibition.

3. **The enabling work is a dispersion-aware energy window**, owned by
   `apply_solver_window` (the window half of the `derive_solver` consolidation,
   U-C). Strategy: the **Gershgorin envelope over the consumer's visited
   k-points** — for a sweep, Gershgorin at the endpoints (k=0 and k_max) unioned
   with margin, yielding one *stable* window per sweep; for a single-k consumer
   (g-factor at Γ, optics at fixed k), the one-k bound. `auto_compute_energy_window`
   remains as the Gershgorin primitive; `apply_solver_window` orchestrates it.

4. **One stable window per sweep** is a hard constraint: the wire kz-sweep does
   branch tracking against the previous k-point (`main.f90:314`), which requires
   the *set* of eigenvalues found to be comparable point-to-point. Per-k moving
   windows are therefore rejected as the default strategy even though per-k
   Gershgorin is cheap.

## Why

- **Consistency.** A single executable must not honor `[solver].method` for one
  geometry and silently override it for another. The #8 silent override is an
  ADR-0002 violation; FEAST parity removes the inconsistency by construction.
- **Performance.** Large QW, large-`nx` Landau, and QW-optics accumulations pay
  dearly for full dense diagonalization. Sparse FEAST iteration is the reason
  the sparse path exists.
- **Honesty.** Rejections should describe reality. "Not yet implemented" is
  truth; "invalid" would be a false design claim that blocks legitimate future
  work.
- **Existence proof.** The QW band-sweep already runs FEAST+ENERGY successfully,
  so the pattern is demonstrated, not speculative.

## Consequences

- `apply_solver_window` becomes the single window authority. `auto_compute_energy_window`
  stays public as the Gershgorin primitive; the envelope logic lives above it.
- The QW g-factor (#2) and QW optics (#8) rejections are reworded "not yet
  implemented" and tracked as flip-on items, not left as silent overrides.
- **Each consumer is enabled behind a FEAST-vs-dense regression test**
  (eigenvalues bit-identical within tolerance against the dense path on the same
  config) before its rejection is removed. Order: QW g-factor → QW optics →
  Landau, following increasing integration risk.
- Landau+FEAST must respect ADR 0003 (the B-sweep and fan diagram stay in
  `main.f90`); only the k-solve backend changes.
- Bulk never offers FEAST — `[solver].method = FEAST` with `confinement = bulk`
  is rejected (validate) as pointless, not "not yet."
- **Do not re-introduce a per-k moving window** for sweeps "for accuracy." A
  future maintainer who recomputes the window every k-point — even if it looks
  more accurate — will break branch tracking. The single-envelope window and the
  branch-tracking invariant are coupled by design.

## Related

- `CONTEXT.md` — glossary entries *AUTO*, *requested band window*, and
  *eigensolver dispatch (format vs backend)*.
- ADR 0002 — config validation consolidation; the `validate()` / I15 layer that
  rejects explicit invalid combos and that the "not yet implemented" messages
  extend.
- ADR 0003 — Landau B-sweep stays in `main.f90`; Landau+FEAST changes only the
  k-solve backend, not the sweep ownership.
- ADR 0004 — `AUTO` resolves method-aware (`FEAST → ENERGY`); the resolution
  contract this ADR depends on.
