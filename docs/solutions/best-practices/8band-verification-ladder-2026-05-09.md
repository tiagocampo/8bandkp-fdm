---
title: "8-band verification ladder: validation hierarchy and physics insights"
date: 2026-05-09
category: best-practices
module: hamiltonianConstructor
problem_type: best_practice
component: testing_framework
severity: low
applies_when:
  - Adding new materials or modifying Hamiltonian construction
  - Debugging band structure discrepancies between model and literature
  - Extending the k.p model to new geometries or confinement dimensions
tags: [verification, effective-mass, kane-model, non-parabolicity, bulk, qw, wire, validation]
---

# 8-band verification ladder: validation hierarchy and physics insights

## Context

The 8-band zinc-blende k.p code had 24 regression tests comparing output against
golden files, but coverage was uneven and undocumented. When a test failed, the
developer had to manually trace whether the error was in parameters, Hamiltonian
construction, FD stencils, confinement, or sparse assembly. The verification
ladder creates a 4-rung hierarchy where each rung is independently confirmed,
making residual errors unambiguously attributable to the model rather than the
implementation.

## Guidance

### Rung structure

Build validation from simplest to most complex physics:

1. **Bulk k=0** (structural) — eigenvalues match diagonal parameters exactly by
   construction. Tests parameters, Hamiltonian assembly, eigenvalue solver.
2. **Bulk dispersion** (curvature) — effective masses from parabolic fitting.
   Tests off-diagonal k-coupling blocks. Self-consistency via Kane model.
3. **QW subbands** (confinement) — subband energies and degeneracies. Tests FD
   stencils, confinement initialization, band offsets.

   **Watch for overlapping material layers.** `confinement_init.f90` uses
   last-layer-wins assignment — a layer covering the full domain after a
   narrower well layer silently overwrites the well. Use the 2-layer pattern
   (barrier covers full domain, well overwrites central region) and avoid
   redundant full-domain barrier layers. See
   `docs/solutions/logic-errors/standard-star-benchmark-suite-logic-errors-2026-05-09.md`
   for a documented instance.
4. **Wire** (2D + sparse) — convergence and dense-sparse agreement. Tests CSR
   assembly, wire geometry, eigensolver paths.

### Critical physics insights

**CB effective mass is NOT the Vurgaftman tabulated value.** The 8-band model
produces inherently non-parabolic CB dispersion. GaAs gives m* ~ 0.0476 m0
vs Vurgaftman 0.067 (29% deviation). This is a known model feature, not a bug.
The correct self-consistency check is the 2-band Kane model:

```
m*_kane = Eg / (EP + Eg)
```

with ~10% tolerance to account for higher-order band-mixing beyond the 2-band
approximation. Actual results: GaAs 4.9%, InAs 1.2%, InSb 7.6%, GaSb 1.2%.

**Do not compare 8-band QW subbands against the Bastard infinite-barrier formula.**
The single-band effective mass approximation gives 56.2 meV for GaAs/AlGaAs;
the 8-band gives 9.92 meV — a 5.4x discrepancy. The VB-CB coupling
renormalizes confinement energies dramatically. Instead, compare against
benchmark configs with known 8-band results.

**Bulk mode does NOT use EV/EC.** `ZB8bandBulk` only places Eg and DeltaSO on
the diagonal; HH/LH are pinned to 0.0. GaSb non-W (which has EV=0, EC=0 in
parameters.f90) works correctly for bulk. EV/EC only matters for QW band
offsets and W-variant material alignment.

**Eigenvalue ordering is ascending.** At k=0 the 8 eigenvalues sort as:
`[-DeltaSO, -DeltaSO, 0, 0, 0, 0, Eg, Eg]` — the SO states are the lowest,
not the HH/LH states.

### Mass extraction method

Use numerical differentiation at the first nonzero k-point, not a polynomial
fit over many points:

```python
d2E_dk2 = 2.0 * (E(k1) - E(0)) / k1**2
m_star = hbar2_over_2m0 / d2E_dk2   # hbar2_over_2m0 = 3.81 eV*A^2
```

This gives the true k->0 parabolic limit without contamination from
non-parabolic curvature at larger k. Polynomial fits with a linear term
absorb non-parabolicity and pollute the quadratic coefficient.

### Wire dense-sparse comparison

Use the existing 11x11 grid configs (`wire_gaas_dense_sparse_*.cfg`) for
dense vs FEAST comparison at 1e-8 tolerance. Larger grids make dense LAPACK
impractically slow. Only compare at k=0 — FEAST is non-deterministic at
higher k-points.

## Why This Matters

Without layered validation, a failing test could mean anything: wrong
parameters, broken Hamiltonian construction, FD stencil errors, confinement
bugs, or sparse assembly issues. The ladder isolates each layer. When Rung 1
passes but Rung 2 fails, the problem is in the off-diagonal k-coupling — not
in parameters or eigenvalue sorting.

The 20-30% effective mass deviation from Vurgaftman is the single most
common source of confusion. Without understanding that this is expected
8-band physics, developers waste time chasing a non-existent bug. Documenting
it as a known model limitation turns confusion into a quantitative diagnostic.

## When to Apply

- When adding a new material to `parameters.f90` — add it to Rungs 1 and 2
- When modifying `hamiltonianConstructor.f90` — run all 4 rungs to isolate
  which physics layer is affected
- When debugging band structure discrepancies — check which rung fails first
- When extending to new geometries — add a new rung that tests only the new
  confinement dimension

## Examples

Running the full ladder:

```bash
ctest --test-dir build -L verification -j4
```

Individual rung for targeted debugging:

```bash
python3 tests/integration/verify_8band_rung1_bulk_k0.py build .
python3 tests/integration/verify_8band_rung2_dispersion.py build .
python3 tests/integration/verify_8band_rung3_qw.py build .
python3 tests/integration/verify_8band_rung4_wire.py build .
```

## Related

- `docs/reference/benchmarks.md` Section 7 — verification ladder results
- `docs/plans/archive/2026-05-08-feat-8band-verification-ladder-plan.md` — implementation plan
- `docs/brainstorms/2026-05-08-8band-verification-ladder-requirements.md` — R1-R19 requirements
- `docs/solutions/logic-errors/standard-star-benchmark-suite-logic-errors-2026-05-09.md` — standard-star benchmark bugs including layer-overlap pitfall
