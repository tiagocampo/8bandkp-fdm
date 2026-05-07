# Figure Provenance Inventory

Purpose: record, for every figure under `docs/figures/`, the active generator,
the config and output dependencies, the figure type, and the current trust
level for the revamp.

Trust key:
- `validated-quantitative`: benchmarked or internally hardened enough for
  quantitative use in its current scope.
- `keep-explanatory`: acceptable only as an explanatory schematic/derived
  figure; not evidence for quantitative claims.
- `needs benchmark`: traceable figure that still needs external or
  chapter-level validation before strong claims.
- `remove/downgrade`: known false, mislabeled, overclaimed, or otherwise not
  safe to present quantitatively.

Type key:
- `computed`: driven directly from solver outputs or controlled sweeps.
- `derived`: post-processed, hybrid, cached, or manually summarized from
  solver outputs and/or chapter tables.
- `schematic`: hand-drawn explanatory figure, not solver output.

Alignment note: trust levels are keyed to
`docs/plans/2026-04-21-benchmark-matrix.md` and
`docs/plans/2026-04-21-failure-ledger.md`.

Cell style:
- `primary`: main config or output artifact.
- `secondary`: companion run or second config.
- `dynamic`: config built or modified inside the plotting script.
- `cache`: persisted sweep or benchmark cache file.
- `fallback`: alternate file used when the primary artifact is absent.
- `source`: non-solver origin such as a chapter table or literature-only block.

## Bulk And Benchmarks

| Figure | Generator | Config(s) | Output file(s) | Type | Trust |
|---|---|---|---|---|---|
| `benchmark_bulk_gaas.png` | `generate_benchmark_figures.py::benchmark_bulk_gaas_block` | `primary: tests/regression/configs/bulk_gaas_kx.cfg; cache: tests/regression/data/bulk_gaas_kx/eigenvalues.dat` | `cache: tests/regression/data/bulk_gaas_kx/eigenvalues.dat` | `derived` | `validated-quantitative` |
| `benchmark_gaas_algaas_qw.png` | `generate_benchmark_figures.py::benchmark_gaas_algaas_qw_block` | `primary: docs/benchmarks/qw_gaas_algaas.cfg; cache: docs/benchmarks/output_gaas_algaas/eigenvalues.dat` | `cache: docs/benchmarks/output_gaas_algaas/eigenvalues.dat` | `derived` | `validated-quantitative` |
| `benchmark_inasw_gasbw_broken_gap.png` | `fig_benchmark_inasw_gasbw_broken_gap` | `tests/regression/configs/qw_inasw_gasbw_broken_gap.cfg` | `output/eigenvalues.dat` | `computed` | `needs benchmark` |
| `bulk_gaas_bands.png` | `fig_bulk_gaas_bands` | `primary: tests/regression/configs/bulk_gaas_kx.cfg` | `primary: output/eigenvalues.dat` | `computed` | `validated-quantitative` |
| `bulk_gaas_bands_110.png` | `fig_bulk_gaas_bands_110` | `primary: tests/regression/configs/bulk_gaas_kxky.cfg` | `primary: output/eigenvalues.dat` | `computed` | `validated-quantitative` |
| `bulk_gaas_parts.png` | `fig_bulk_gaas_parts` | `primary: tests/regression/configs/bulk_gaas_kx.cfg` | `primary: output/parts.dat` | `computed` | `validated-quantitative` |
| `bulk_gaas_parts_vs_k.png` | `fig_bulk_gaas_parts_vs_k` | `primary: tests/regression/configs/bulk_gaas_kxky.cfg` | `primary: output/parts.dat` | `computed` | `validated-quantitative` |
| `bulk_gaas_strain_comparison.png` | `fig_bulk_gaas_strain_comparison` | `primary: tests/regression/configs/bulk_gaas_kxky.cfg; secondary: tests/regression/configs/bulk_gaas_strained.cfg` | `primary: output/eigenvalues.dat; secondary: output/eigenvalues.dat` | `computed` | `validated-quantitative` |
| `bulk_gaas_strained_bands.png` | `fig_bulk_gaas_strained_bands` | `primary: tests/regression/configs/bulk_gaas_strained.cfg` | `primary: output/eigenvalues.dat` | `computed` | `validated-quantitative` |
| `bulk_gaas_warping.png` | `fig_bulk_gaas_warping` | `primary: tests/regression/configs/bulk_gaas_kx.cfg; secondary: tests/regression/configs/bulk_gaas_kxky.cfg` | `primary: output/eigenvalues.dat; secondary: output/eigenvalues.dat` | `computed` | `validated-quantitative` |
| `bulk_inas_bands.png` | `fig_bulk_inas_bands` | `tests/regression/configs/bulk_inas_kx.cfg` | `output/eigenvalues.dat` | `computed` | `needs benchmark` |

## Quantum Wells And Wavefunctions

| Figure | Generator | Config(s) | Output file(s) | Type | Trust |
|---|---|---|---|---|---|
| `cb_parts_evolution.png` | `fig_cb_parts_evolution` | `primary: tests/regression/configs/qw_alsbw_gasbw_inasw.cfg` | `primary: output/eigenvalues.dat; secondary: output/eigenfunctions_k_*.dat` | `computed` | `needs benchmark` |
| `double_qw_anticrossing.png` | `fig_double_qw_anticrossing` | `tests/regression/configs/qw_gaas_algaas_double_qw.cfg` | `output/eigenvalues.dat` | `computed` | `needs benchmark` |
| `double_qw_potential_profile.png` | `fig_double_qw_potential_profile` | `primary: tests/regression/configs/qw_gaas_algaas_double_qw.cfg` | `primary: output/potential_profile.dat; secondary: output/eigenvalues.dat` | `computed` | `needs benchmark` |
| `double_qw_wavefunctions.png` | `fig_double_qw_wavefunctions` | `primary: tests/regression/configs/qw_gaas_algaas_double_qw.cfg` | `primary: output/eigenvalues.dat; secondary: output/eigenfunctions_k_*.dat` | `computed` | `needs benchmark` |
| `eigenvector_block_structure.png` | `fig_eigenvector_block_structure` | `source: none` | `source: none` | `schematic` | `keep-explanatory` |
| `perband_density.png` | `fig_perband_density` | `tests/regression/configs/qw_alsbw_gasbw_inasw.cfg` | `output/eigenfunctions_k_*.dat` | `computed` | `remove/downgrade` |
| `qw_alsbw_gasbw_inasw_bands.png` | `fig_qw_alsbw_gasbw_inasw_bands` | `tests/regression/configs/qw_alsbw_gasbw_inasw.cfg` | `output/eigenvalues.dat` | `computed` | `needs benchmark` |
| `qw_dispersion_broken_gap.png` | `fig_qw_dispersion_broken_gap` | `tests/regression/configs/qw_inas_gasb_broken_gap_kpar.cfg` | `output/eigenvalues.dat` | `computed` | `needs benchmark` |
| `qw_dispersion_gaas_algaas.png` | `fig_qw_dispersion_gaas_algaas` | `tests/regression/configs/qw_gaas_algaas_kpar.cfg` | `output/eigenvalues.dat` | `computed` | `needs benchmark` |
| `qw_gaas_algaas_subbands.png` | `fig_qw_gaas_algaas_subbands` | `primary: docs/benchmarks/qw_gaas_algaas.cfg` | `primary: output/eigenvalues.dat` | `computed` | `validated-quantitative` |
| `qw_parts.png` | `fig_qw_parts` | `tests/regression/configs/qw_alsbw_gasbw_inasw.cfg` | `output/parts.dat` | `computed` | `needs benchmark` |
| `qw_parts_gaas.png` | `fig_qw_parts_gaas` | `primary: docs/benchmarks/qw_gaas_algaas.cfg` | `primary: output/parts.dat` | `computed` | `validated-quantitative` |
| `qw_potential_profile.png` | `fig_qw_potential_profile` | `tests/regression/configs/qw_alsbw_gasbw_inasw.cfg` | `output/potential_profile.dat` | `computed` | `needs benchmark` |
| `qw_potential_profile_gaas.png` | `fig_qw_potential_profile_gaas` | `tests/regression/configs/qw_gaas_algaas_kpar.cfg` | `output/potential_profile.dat`; `output/eigenvalues.dat` | `computed` | `needs benchmark` |
| `qw_wavefunctions.png` | `fig_qw_wavefunctions` | `tests/regression/configs/qw_alsbw_gasbw_inasw.cfg` | `output/eigenfunctions_k_*.dat` | `computed` | `remove/downgrade` |
| `qw_wavefunctions_gaas.png` | `fig_qw_wavefunctions_gaas` | `primary: docs/benchmarks/qw_gaas_algaas.cfg` | `primary: output/eigenfunctions_k_*.dat` | `computed` | `validated-quantitative` |
| `vb_hh_lh_mixing.png` | `fig_vb_hh_lh_mixing` | `primary: tests/regression/configs/qw_alsbw_gasbw_inasw.cfg` | `primary: output/eigenvalues.dat; secondary: output/eigenfunctions_k_*.dat` | `computed` | `needs benchmark` |

## Strain

| Figure | Generator | Config(s) | Output file(s) | Type | Trust |
|---|---|---|---|---|---|
| `bir_pikus_band_shifts.png` | `fig_bir_pikus_band_shifts` | `source: none` | `source: none` | `derived` | `keep-explanatory` |
| `hh_lh_ordering.png` | `fig_hh_lh_ordering` | `source: none` | `source: none` | `derived` | `keep-explanatory` |
| `hh_lh_splitting_vs_mismatch.png` | `fig_hh_lh_splitting_vs_mismatch` | `source: Chapter 04 strain table` | `source: Chapter 04 table values` | `derived` | `needs benchmark` |
| `qw_strained_band_edges.png` | `fig_qw_strained_band_edges` | `dynamic: temporary unstrained QW config; source: manual Bir-Pikus shifts` | `primary: output/potential_profile.dat` | `derived` | `needs benchmark` |
| `qw_strained_bands.png` | `fig_qw_strained_bands` | `primary: tests/regression/configs/qw_alsbw_gasbw_inasw.cfg; dynamic: temporary strain-on copy` | `primary: output/eigenvalues.dat` | `computed` | `needs benchmark` |
| `strain_biaxial_tensor.png` | `fig_strain_biaxial_tensor` | `source: none` | `source: none` | `schematic` | `keep-explanatory` |
| `strain_lattice_mismatch.png` | `fig_strain_lattice_mismatch` | `source: none` | `source: none` | `schematic` | `keep-explanatory` |
| `strained_unit_cell.png` | `fig_strained_unit_cell` | `source: none` | `source: none` | `schematic` | `keep-explanatory` |

## g-Factor

| Figure | Generator | Config(s) | Output file(s) | Type | Trust |
|---|---|---|---|---|---|
| `benchmark_gfactor_comparison.png` | `generate_benchmark_figures.py::benchmark_gfactor_comparison_block` | `source: literature-only hard-coded values` | `source: none` | `derived` | `remove/downgrade` |
| `gfactor_components.png` | `fig_gfactor_components` | `bulk-cb: tests/regression/configs/gfactor_bulk_gaas_cb.cfg; bulk-vb: gfactor_bulk_gaas_vb.cfg; qw-cb: gfactor_qw_cb.cfg; qw-vb: gfactor_qw_vb.cfg` | `primary: output/gfactor.dat` | `computed` | `needs benchmark` |
| `gfactor_zeeman.png` | `fig_gfactor_zeeman` | `tests/regression/configs/gfactor_bulk_gaas_cb.cfg` | `output/gfactor.dat` | `derived` | `remove/downgrade` |

## Self-Consistent Potential

| Figure | Generator | Config(s) | Output file(s) | Type | Trust |
|---|---|---|---|---|---|
| `sc_charge_density.png` | `fig_sc_charge_density` | `primary: tests/regression/configs/sc_gaas_alas_qw.cfg` | `primary: output/sc_charge.dat` | `computed` | `validated-quantitative` |
| `sc_convergence.png` | `fig_sc_convergence` | `primary: tests/regression/configs/sc_gaas_alas_qw.cfg` | `primary: stdout convergence history` | `computed` | `validated-quantitative` |
| `sc_inas_alsb_convergence.png` | `fig_sc_inas_alsb_convergence` | `primary: tests/regression/configs/sc_qw_inas_alsb.cfg` | `primary: stdout convergence history; fallback: output/sc_potential_profile.dat; secondary: output/eigenvalues.dat` | `computed` | `needs benchmark` |
| `sc_inas_alsb_potential.png` | `fig_sc_inas_alsb_potential` | `primary: tests/regression/configs/sc_qw_inas_alsb.cfg` | `primary: output/sc_potential_profile.dat; fallback: output/potential_profile.dat; secondary: output/eigenvalues.dat` | `computed` | `needs benchmark` |
| `sc_potential.png` | `fig_sc_potential` | `primary: tests/regression/configs/sc_gaas_alas_qw.cfg` | `primary: output/sc_potential_profile.dat; fallback: output/potential_profile.dat` | `computed` | `validated-quantitative` |

## Optics, Excitons, Scattering, QCSE

| Figure | Generator | Config(s) | Output file(s) | Type | Trust |
|---|---|---|---|---|---|
| `absorption_excitonic_TE.png` | `fig_absorption_excitonic_TE` | `primary: tests/regression/configs/qw_gaas_algaas_optics_full.cfg` | `primary: output/absorption_TE.dat; secondary: output/exciton.dat` | `derived` | `remove/downgrade` |
| `absorption_with_exciton.png` | `fig_absorption_with_exciton` | `source: no run; uses existing TE absorption only` | `primary: output/absorption_TE.dat` | `derived` | `remove/downgrade` |
| `exciton_binding_vs_width.png` | `fig_exciton_binding_vs_width` | `dynamic: Python-built width sweep` | `primary: output/exciton.dat; cache: output/exciton_width_sweep.dat` | `computed` | `needs benchmark` |
| `exciton_bohr_vs_width.png` | `fig_exciton_bohr_vs_width` | `dynamic: Python-built width sweep` | `primary: output/exciton.dat; cache: output/exciton_width_sweep.dat` | `computed` | `needs benchmark` |
| `gain_strained_comparison.png` | `fig_gain_strained_comparison` | `primary: tests/regression/configs/qw_gaas_algaas_optics_full.cfg` | `primary: output/gain_TE.dat; secondary: output/gain_TM.dat` | `computed` | `needs benchmark` |
| `isbt_absorption.png` | `fig_isbt_absorption` | `tests/regression/configs/qw_gaas_algaas_isbt.cfg` | `output/absorption_ISBT.dat` | `computed` | `needs benchmark` |
| `isbt_dipole_moments.png` | `fig_isbt_dipole_moments` | `primary: tests/regression/configs/qw_gaas_algaas_isbt.cfg; fallback: existing output only` | `primary: output/isbt_transitions.dat` | `computed` | `needs benchmark` |
| `qcse_stark_shift.png` | `fig_qcse_stark_shift` | `primary: tests/regression/configs/sc_qcse_gaas_algaas.cfg; secondary: sc_qcse_gaas_algaas_ef.cfg` | `primary: output/potential_profile.dat; secondary: output/eigenvalues.dat` | `derived` | `remove/downgrade` |
| `qw_absorption_spectrum.png` | `fig_qw_absorption_spectrum` | `primary: tests/regression/configs/qw_gaas_algaas_optics_full.cfg` | `primary: output/absorption_TE.dat; secondary: output/absorption_TM.dat` | `computed` | `remove/downgrade` |
| `qw_absorption_strained.png` | `fig_qw_absorption_strained` | `primary: tests/regression/configs/qw_gaas_algaas_optics_full.cfg; secondary: tests/regression/configs/qw_ingaas_algaas_strained_optics.cfg` | `unstrained: output/absorption_TE_unstrained.dat + output/absorption_TM_unstrained.dat; strained: output/absorption_TE_strained.dat + output/absorption_TM_strained.dat` | `computed` | `remove/downgrade` |
| `qw_absorption_vs_width.png` | `fig_qw_absorption_vs_width` | `dynamic: manual width sweep based on optics configs` | `primary: output/absorption_TE_W*.dat; secondary: output/absorption_TM_W*.dat; fallback: output/absorption_TE.dat; fallback: output/absorption_TM.dat` | `computed` | `remove/downgrade` |
| `qw_optical_matrix_elements.png` | `fig_qw_optical_matrix_elements` | `tests/regression/configs/qw_gaas_algaas_optics.cfg` | `output/optical_transitions.dat` | `computed` | `needs benchmark` |
| `scattering_lifetime_vs_field.png` | `fig_scattering_lifetime_vs_field` | `primary: tests/regression/configs/qw_gaas_algaas_qcse_scattering.cfg; dynamic: field sweep override` | `primary: output/scattering_rates.dat; cache: output/scattering_field_sweep.dat` | `computed` | `remove/downgrade` |
| `scattering_lifetime_vs_width.png` | `fig_scattering_lifetime_vs_width` | `primary: tests/regression/configs/qw_gaas_algaas_qcse_scattering.cfg; fallback: equivalent scattering run` | `primary: output/scattering_rates.dat` | `derived` | `remove/downgrade` |

## Numerics

| Figure | Generator | Config(s) | Output file(s) | Type | Trust |
|---|---|---|---|---|---|
| `convergence_fd_order.png` | `fig_convergence_fd_order` | `primary: tests/regression/configs/qw_alsbw_gasbw_inasw.cfg; dynamic: FDorder sweep` | `primary: output/eigenvalues.dat` | `computed` | `validated-quantitative` |
| `convergence_grid_spacing.png` | `fig_convergence_grid_spacing` | `primary: tests/regression/configs/qw_alsbw_gasbw_inasw.cfg; dynamic: FDstep sweep` | `primary: output/eigenvalues.dat` | `computed` | `validated-quantitative` |
| `timing_dense_vs_sparse.png` | `fig_timing_dense_vs_sparse` | `primary: tests/regression/configs/qw_alsbw_gasbw_inasw.cfg; secondary: tests/regression/configs/wire_gaas_rectangle.cfg` | `source: wall-clock timing only` | `derived` | `remove/downgrade` |

## Quantum Wires

| Figure | Generator | Config(s) | Output file(s) | Type | Trust |
|---|---|---|---|---|---|
| `wire_density_2d.png` | `fig_wire_density_2d` | `primary: tests/regression/configs/wire_gaas_rectangle.cfg` | `primary: output/eigenvalues.dat; secondary: output/eigenfunctions_k_*.dat` | `computed` | `needs benchmark` |
| `wire_gfactor_vs_size.png` | `fig_wire_gfactor_vs_size` | `dynamic: Python-built GaAs wire sweep; secondary: tests/regression/configs/wire_gaas_rectangle.cfg` | `primary: output/gfactor.dat` | `computed` | `remove/downgrade` |
| `wire_inas_gaas_profile.png` | `fig_wire_inas_gaas_profile` | `tests/regression/configs/wire_inas_gaas_strain.cfg` | `output/potential_profile.dat` | `computed` | `needs benchmark` |
| `wire_inas_gaas_subbands.png` | `fig_wire_inas_gaas_subbands` | `primary: tests/regression/configs/wire_inas_gaas_strain.cfg; dynamic: modified kz sweep copy` | `primary: output/eigenvalues.dat` | `computed` | `needs benchmark` |
| `wire_inas_gaas_wavefunctions.png` | `fig_wire_inas_gaas_wavefunctions` | `primary: tests/regression/configs/wire_inas_gaas_strain.cfg` | `primary: output/eigenvalues.dat; secondary: output/eigenfunctions_k_*.dat` | `computed` | `needs benchmark` |
| `wire_strain_2d.png` | `fig_wire_strain_2d` | `tests/regression/configs/wire_inas_gaas_strain.cfg` | `output/strain.dat` | `computed` | `remove/downgrade` |
| `wire_strain_tensor.png` | `fig_wire_strain_tensor` | `tests/regression/configs/wire_inas_gaas_strain.cfg` | `output/strain.dat` | `computed` | `needs benchmark` |
| `wire_subbands.png` | `fig_wire_subbands` | `tests/regression/configs/wire_gaas_rectangle.cfg` | `output/eigenvalues.dat` | `computed` | `needs benchmark` |
