# Phase 4 Figure Discrepancy Log

Run date: 2026-05-04
Script: `scripts/plotting/generate_all_figures.py --skip-build`
Result: **74 succeeded, 1 failed out of 75 registered figures**

Total PNG files generated: 77 in `docs/figures/` + 9 in `docs/lecture/figures/` = 86 PNG files
(some figure functions produce multiple PNG files, e.g. auxiliary panels)

3 stale benchmark files from 2026-04-25 not in current ALL_FIGURES dict: `benchmark_bulk_gaas.png`, `benchmark_gaas_algaas_qw.png`, `benchmark_gfactor_comparison.png`

| Figure | Issue | Root Cause | Fix | Status |
|--------|-------|------------|-----|--------|
| `qw_strained_bands` | Strained QW run segfaults (rc=-11). Unstrained panel works. Stale PNG from 2026-04-25 remains. | SIGSEGV in eigensolver after strain modification of QW Hamiltonian with `strain_ref: AlSbW` on `qw_alsbw_gasbw_inasw.cfg`. Strain calculation completes ("QW strain calculation complete") but eigensolver crashes. Likely array bounds or memory corruption when strain shifts Hamiltonian elements. Wire strain (confinement=2) works fine; only QW strain (confinement=1) fails. | Not yet fixed. Needs debugging in `hamiltonianConstructor.f90` or `strain_solver.f90` for QW strain path. | Open |
| `timing_dense_vs_sparse` | Wire (sparse) run exceeds 300s timeout. Stale PNG from 2026-04-25 remains. | `wire_gaas_rectangle.cfg` with full kz-sweep takes >300s on this hardware. The `run_executable` default timeout is 300s. | Increase timeout in `fig_timing_dense_vs_sparse` to 600s, or reduce wire grid size for benchmark. | Open |

## Physics Sanity Checks (all pass)

| Check | Expected | Computed | Status |
|-------|----------|----------|--------|
| GaAs E_g at Gamma | ~1.42 eV | 1.519 eV | OK (8-band parameter set) |
| GaAs SO splitting | ~0.34 eV | 0.341 eV | OK |
| InAs E_g at Gamma | ~0.354 eV | 0.417 eV | OK (8-band parameter set) |
| InAs SO splitting | ~0.38 eV | 0.390 eV | OK |
| InGaAs QW HH-LH split | ~50 meV | 54.0 meV | OK |
| InAs/GaAs wire gap | ~0.8 eV | 0.7749 eV | OK |
| QW g-factor (InAsW, 50A) | large negative | gx=gy=-14.9, gz=-10.3 | OK |
| Wire g-factor (GaAs, 100A) | moderate negative | gx=gy=-3.4, gz=2.1 | OK |

## All generated figures (77 PNG files in docs/figures/)

All files >10KB, no suspiciously small or zero-byte files detected.
