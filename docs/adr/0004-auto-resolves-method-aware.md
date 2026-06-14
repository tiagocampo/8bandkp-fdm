# ADR 0004: `AUTO` resolves method-aware

## Status

Accepted

## Context

The eigensolver-standardization refactor unified every diagonalization call site behind a polymorphic interface and a `[solver]` TOML section with `method`/`mode` dispatch. Each confinement mode (bulk, QW, wire, Landau) resolved the `AUTO` sentinel to concrete `method`/`mode` values in its own duplicated `select case` block. Two of those blocks resolved `AUTO` **method-blind** — i.e. they picked the confinement's default method, then picked the confinement's default mode *independently*, without checking whether the two were compatible.

The symptom: a physicist setting `method = "FEAST"` while leaving `mode = "AUTO"` on a QW or Landau calculation got `FEAST + INDEX`, which is an **invalid combination** — FEAST has no INDEX (range `il:iu`) interface. The program aborted deep in the eigensolver with an opaque error instead of either picking a valid mode or rejecting the input up front. This directly contradicted the design goal that "a physicist should be able to choose a solver without fighting defaults."

The four duplicated resolution blocks had also already drifted slightly inconsistent with each other, and the Landau block *unconditionally* forced `INDEX`, ignoring even an explicit user `mode` — a second smart-default defect.

## Decision

`AUTO` is **method-aware**. Resolution is centralized in a single pure subroutine, `resolve_solver_defaults(confinement, method, mode, out_method, out_mode)` in `src/core/defs.f90`, called by all four confinement dispatch blocks in `simulation_setup.f90` and by the band-structure sweep in `main.f90`. The contract (codified in the root `CONTEXT.md` glossary under *AUTO*) is:

1. **Method resolves first** — confinement default unless the user set one (bulk/Landau/QW → DENSE; wire → FEAST).
2. **Mode then resolves method-aware** — if the user set a mode, honor it; if `AUTO`, pick a default **compatible with the resolved method**: FEAST → ENERGY; DENSE → the confinement's native mode (bulk → FULL; QW/Landau → INDEX; wire → ENERGY).
3. **Invariant** — `AUTO` never returns an invalid combination. `FEAST + INDEX` is unreachable through `AUTO`.

This is a deep module: a tiny, stable interface hiding the resolution rules and the no-invalid-combo invariant. The confinement-specific `nev`/`il`/`iu` bookkeeping stays in each dispatch block; only the method/mode resolution is centralized.

A user who *explicitly* sets `method = "FEAST"` and `mode = "INDEX"` is still rejected — by `validate()` check I15 (`src/core/defs.f90`) on the raw, unresolved config — but `AUTO` itself can never produce that pair.

## Why

Method-aware resolution is the only way to honor the "choose a solver without fighting defaults" goal. A method-blind `AUTO` (resolve method, then resolve mode, independently) will always have a corner case where the two defaults are incompatible; for this eigensolver that corner is exactly `FEAST + INDEX`, which is the common benchmarking request (FEAST vs dense on a QW).

Centralizing the rule in one pure function means the four confinement modes cannot disagree about what `AUTO` means (user story: "I cannot accidentally make two confinement modes disagree"), and adding a future eigensolver variant is a one-line resolution change that takes effect for every confinement at once.

## Consequences

- `FEAST + INDEX` is unreachable via `AUTO` by construction and is asserted by a table-driven unit test (`tests/unit/test_eigensolver.pf`) covering every `confinement × method × mode` combination plus a dedicated no-invalid-combo invariant test.
- Bit-identical to prior behavior on every shipped configuration: the only resolution changes are `FEAST + AUTO` combos, which either previously aborted (QW/Landau) or were absent from the golden suite (bulk); wire `FEAST + AUTO` is unchanged (`ENERGY → ENERGY`). Landau now honors an explicit mode instead of forcing `INDEX` (unobservable today since no Landau config sets one).
- **Do not "simplify" this back to method-blind resolution.** A future maintainer who resolves the method default and the mode default independently — even if it looks equivalent for the configurations they are staring at — reintroduces the invalid-combination trap. The invariant and its unit test exist precisely to make such a regression visible.

## Related

- `CONTEXT.md` — glossary entries *AUTO* and *requested band window*.
- ADR 0002 — config validation consolidation (the `validate()` + I15 layer that rejects explicit invalid combos).
