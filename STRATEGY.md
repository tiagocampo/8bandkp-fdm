---
name: 8-band k.p FDM
last_updated: 2026-05-08
---

# 8-band k.p FDM Strategy

## Target problem

Semiconductor researchers need full 8-band k.p results — g-factors, optical spectra, topological invariants — but resort to approximations, patched legacy Fortran, or hacked Python notebooks because no unified codebase covers all of it. The crux: the full 8-band model is too complex for ad-hoc implementations, so researchers compromise on physics accuracy instead.

## Our approach

One Fortran codebase that eliminates the trade-off between physics completeness and computational feasibility — full 8-band k.p from bulk to confined systems, with dense and sparse solvers matched to each problem's structure.

## Who it's for

**Primary:** Academic researchers in semiconductor/condensed matter physics — PhD students, postdocs, PIs. They're hiring it to compute physically correct, end-to-end traceable results for quantum nanostructures without stitching together approximations.

## Key metrics

- **Published validation results** — number of literature benchmarks reproduced (Vurgaftman band gaps, Winkler g-factors, BHZ phase diagrams); tracked in regression tests
- **Physics coverage** — count of distinct end-to-end scenarios that work (bulk, QW, wire, g-factor, optics, SC Poisson, BdG/Majorana)
- **Reproducibility** — regression test count and pass rate; can anyone clone, build, and reproduce published results
- **External adoption** — papers published by others using or citing this code
- **Benchmark performance** — time-to-solution for standard problems vs. comparable tools

## Tracks

### Core k.p physics

Hamiltonian construction, FD stencils, material parameters, bulk/QW/wire modes.

_Why it serves the approach:_ The foundation that makes physics completeness possible.

### Spectroscopic observables

g-factors via Lowdin partitioning, optical spectra (absorption, gain, ISBT), spin-resolved outputs.

_Why it serves the approach:_ Makes results experimentally relevant and publishable.

### Topological and many-body physics

BdG, Majorana, Chern, Z2, self-consistent SP.

_Why it serves the approach:_ The frontier that no other unified k.p code covers — the differentiator.

### Performance and reproducibility

Sparse/dense solver selection, OpenMP scaling, regression tests, clean build system.

_Why it serves the approach:_ What makes it practical instead of a one-off script.

## Marketing

**One-liner:** The one to rule them all.
