# ADR 0006: conf_direction returns a 'w' sentinel for the wire (2D confinement)

## Status

Accepted

## Context

`conf_direction(confinement)` maps each confinement mode to a single character
naming the FDM-discretized confinement axis: bulk → `'n'`, QW → `'z'`, Landau →
`'x'`. The wire was also mapped to `'z'`, and the docblock described it as
"confinement along z."

That was physically wrong and internally contradictory:

- The wire is **2D-confined in the x–y plane** with z as the **free propagation
  axis** (grid docblock, `defs.f90`). So the wire's `'z'` names the *free* axis
  — the opposite role from QW's `'z'` (the *confinement* axis). The project's own
  docs disagreed: `defs.f90` called it "confinement along z"; `lecture 08`
  called the same `'z'` the "free propagation direction."
- The value was consumed by **no live path** — every wire path branches to
  wire-specific code before any `conf_direction` comparison (`main.f90` wire
  block stops at the top; `main_gfactor.f90` wire is the `if` branch; the wire
  uses `confinementInitialization_2d`, not the `_raw` path that takes the
  direction argument). So it was a latent footgun, not a live defect: the next
  generic `if (conf_direction == 'z')` branch added without a prior wire guard
  would silently treat the wire as a QW and flip the axis.

The deeper issue is conceptual: `conf_direction` answers "which single axis is
the system confined along?" — a 1D-confinement question that has no answer for a
2D-confined wire.

Three options were considered:

- **(A)** Keep `'z'`, fix only the docblock to admit it is a masked placeholder.
- **(B)** Return a distinct sentinel (`'w'`) the wire alone produces, so no
  generic `'z'`/`'n'` branch can catch it; update the docblock and unit test.
- **(C)** Redefine the concept (e.g. split into a 1D "confinement axis" function
  the wire is simply not a member of).

## Decision

**Option B**: `conf_direction('wire')` returns `'w'`, a sentinel meaning "wire —
2D confinement, no single confinement axis." It is documented as a sentinel, not
an axis, and the unit test (`test_defs.pf`) asserts `'w'`.

## Why

A 1D-axis function cannot honestly represent 2D confinement; any single
character returned for the wire is a fiction. Option A keeps the fiction and
leaves the footgun — a future contributor reading an accurate docblock still
sees `'z'` in the code and may trust it. Option C is the conceptually cleanest
but invasive for no behavioral gain. Option B makes the function honest *and*
fail-safe at the cost of one character and one test line: no generic axis branch
can ever accidentally consume the wire, and `'w'` self-documents that the wire is
a special case. The pre-change audit proved the wire value has no live consumer,
so the swap is behaviorally inert (all 37 unit tests green).

## Consequences

- Any future code that needs the wire's confinement geometry must branch on
  `confinement == 'wire'` explicitly; it cannot derive it from `conf_direction`.
- `'w'` must never be treated as `'n'` (the wire IS confined) or as an axis.
- The wire's free propagation axis remains z; that fact lives in the grid
  convention and `conf_direction` deliberately does not encode it.
