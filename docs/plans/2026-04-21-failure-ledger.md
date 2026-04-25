# Failure Ledger for Docs Physics Revamp

Purpose: record the known failures and high-risk items already identified before
the solver/doc revamp starts.

Operating rule:
- If an item below affects a quantitative figure or strong physics claim, the
  docs should be downgraded or the figure removed until the item is closed.

Status key:
- `open-blocking`: known false or misleading output; do not present as validated.
- `open-high-risk`: strong chance of solver or post-processing error; benchmark before keeping claims.
- `open-doc`: presentation/completeness issue that does not by itself imply bad physics.

## Open Issues

| ID | Category | Artifact/location | Affects chapters | Current problem | Impact | Required action | Status |
|---|---|---|---|---|---|---|---|
| `FL-001` | broken units | `docs/figures/gfactor_zeeman.png`; `scripts/plotting/generate_all_figures.py` | `05` | Generator uses `E_up = 0.5 * g * B` and `E_dn = -0.5 * g * B` without the Bohr magneton conversion, but `docs/lecture/05-gfactor.md` presents the figure as quantitative. | The plotted slope implies the wrong `g` magnitude and makes the figure physically false. | Fix the units or relabel/remove the figure; add a slope acceptance check. | `open-blocking` |
| `FL-002` | broken units | `docs/figures/qcse_stark_shift.png`; `docs/lecture/10-qcse.md`; `tests/regression/configs/sc_qcse_gaas_algaas_ef.cfg` | `10` | The chapter and figure label a `-70 kV/cm` case, but the config uses `EFParams: -0.007 eV/Angstrom`, which is about `-700 kV/cm`. The generator also carries a field-conversion comment that is off by `10x`. | The QCSE chapter currently overstates a mislabeled example and mixes physical regimes. | Correct the field conversion semantics, regenerate the plot, and benchmark a moderate-field case. | `open-blocking` |
| `FL-003` | bad state indexing | `docs/figures/qw_wavefunctions.png`; `docs/lecture/03-wavefunctions.md`; `scripts/plotting/generate_all_figures.py` | `03` | The chapter says the figure shows the first four conduction states and references `CB1` as state 11, but the generator loads `n_ev = 4` and therefore plots states `1-4`. | Chapter 3 state labels are not traceable to the output. | Load the claimed conduction states or relabel the figure and prose. | `open-blocking` |
| `FL-004` | bad state indexing | `docs/figures/perband_density.png`; `docs/lecture/03-wavefunctions.md`; `scripts/plotting/generate_all_figures.py` | `03` | The rendered figure appears HH-dominated while the caption says the selected state is `CB1` and CB components dominate. | The per-band composition discussion is not trustworthy. | Re-check eigenstate ordering and the selected `wf[...]` index; regenerate or downgrade. | `open-blocking` |
| `FL-005` | parser errors | `docs/figures/scattering_lifetime_vs_width.png`; `output/scattering_rates.dat`; `scripts/plotting/generate_all_figures.py` | `06` | The plot is labeled lifetime vs well width, but the generator reads a single `scattering_rates.dat` file and appears to use transition-index columns as if they were width-dependent data. | The plotted lifetime trend is not derived from the quantity the chapter claims. | Replace with a real width sweep parser or remove the figure. | `open-blocking` |
| `FL-006` | broken units; parser errors | `docs/figures/scattering_lifetime_vs_field.png`; `scripts/plotting/generate_all_figures.py` | `06`, `10` | Field conversion semantics are inconsistent, and the plotted lifetimes are `1e8-1e9 ps`, far from the `6-12 ps` range claimed later in the chapter. | Chapter 6 scattering discussion and validation table are currently not defensible, and the QCSE-linked field path is untrusted. | Fix field conversion, verify rate parsing, and rerun a controlled field sweep. | `open-blocking` |
| `FL-007` | plotting/provenance risk | `docs/figures/qw_absorption_strained.png`; `qw_absorption_vs_width.png`; legacy optics narratives | `06` | The unstrained GaAs/AlGaAs and stripped-down strained InGaAs absorption configs both preserve the expected TE-dominant edge. The remaining risk is legacy figure provenance and any claims built from older manually staged or mixed-workflow outputs. | This is no longer evidence of a broken TE/TM solver path, but old figures/prose can still overstate what was actually benchmarked. | Keep using dedicated absorption-only configs and deterministic figure generation; regenerate or retire any legacy optics figures that did not come from that path. | `open-doc` |
| `FL-008` | suspect solver behavior | Chapter 6 optical validation table and prose | `06` | The chapter claims nextnano-consistent TE/TM splitting and physically reasonable absorption scales, but the current figures do not support those claims. | Validation prose currently runs ahead of the evidence. | Rewrite Chapter 6 only after the optics path is benchmarked and re-plotted. | `open-high-risk` |
| `FL-009` | parser/figure semantics | `docs/figures/absorption_with_exciton.png`; `docs/figures/absorption_excitonic_TE.png`; `scripts/plotting/generate_all_figures.py` | `06` | One figure hardcodes `E_b = 10 meV` while the output file reports about `4.82 meV`; the other only draws a vertical resonance marker rather than a real excitonic peak. | Excitonic corrections are currently overclaimed visually. | Tie the figure to computed exciton outputs or downgrade to schematic language. | `open-blocking` |
| `FL-010` | missing embeds | `docs/lecture/07-self-consistent-sp.md` | `07` | The GaAs/AlAs SP example lists `sc_potential.png`, `sc_charge_density.png`, and `sc_convergence.png` as bullet paths instead of embedding them as figures. | The chapter is incomplete and harder to review against outputs. | Embed the figures and ensure captions reference the actual outputs. | `open-doc` |
| `FL-011` | completeness debt | `docs/lecture/12-extending-the-code.md` | `12` | The chapter is thin, has zero figures, and currently reads like implementation notes rather than a finished lecture chapter. | The lecture set lacks a complete closing chapter under the new evidence standard. | Add one worked extension example plus one architecture/provenance figure, or explicitly mark the chapter as provisional. | `open-doc` |
| `FL-012` | wire-risk items | `src/apps/main_gfactor.f90`; wire optics path discussed in `docs/lecture/06-optical-properties.md` | `06`, `08` | Broad-window wire state selection was using fixed sorted positions rather than actual band-edge adjacency, which could select deep valence states. A gap-aware selection path and dedicated regression now exist, but only for the current wire g-factor/optics driver. | Older wire optical conclusions may still reflect pre-fix runs, and the broader wire benchmark story remains incomplete. | Keep rerun-derived wire optics claims provisional, but treat the core selection bug as fixed for the current driver and regression path. | `open-doc` |
| `FL-013` | wire-risk items | `docs/lecture/08-quantum-wire.md`; `tests/regression/configs/wire_gaas_rectangle.cfg` and related wire configs | `08` | The chapter contains strong claims about wire confinement, strain, and `g`-factor behavior without a completed external reproduction suite. | The sparse wire path is the largest remaining physics-validation gap in the docs. | Build a dedicated wire benchmark suite before retaining strong claims. | `open-high-risk` |
| `FL-014` | wire-risk items | `docs/figures/wire_strain_2d.png`; `docs/lecture/04-strain.md` | `04`, `08` | The caption claims strong interface-corner concentration, but the current figure reads mostly piecewise uniform and has not been benchmarked externally. | The figure may be qualitatively right or merely over-described; either way it is not benchmarked. | Benchmark wire strain maps against symmetry/analytic expectations and soften the caption until then. | `open-high-risk` |
| `FL-015` | suspect solver behavior | `docs/figures/wire_gfactor_vs_size.png`; `docs/lecture/08-quantum-wire.md` | `05`, `08` | The qualitative anisotropy trend may be plausible, but it has not been benchmarked and the chapter presents it too strongly relative to current validation depth. | Wire `g`-factor claims could overstate what the current sparse path really supports. | Keep only as provisional until an InSb radius sweep is reproduced against literature. | `open-high-risk` |
| `FL-016` | suspect solver behavior | `src/physics/scattering.f90`; `docs/lecture/06-optical-properties.md`; `docs/figures/scattering_lifetime_vs_width.png`; `docs/figures/scattering_lifetime_vs_field.png` | `06` | The scattering module currently treats spin-degenerate Kramers partners as distinct CB subbands and builds the Froehlich form factor from band-projected densities rather than wavefunction amplitudes. In the current GaAs QW run this yields near-zero `CB1->CB2` splittings and lifetimes on the order of `10^9 ps`, far outside the cited literature scale. | The scattering figures and any quoted lifetime numbers are not physically trustworthy even before plotting/parser issues are considered. | Rework the scattering-state selection and overlap kernel, then rebuild the sweep/plot path from corrected outputs. | `open-blocking` |

## Cross-Cutting Rules Derived from the Ledger

- Do not regenerate validation tables from memory or old prose; generate them
  only from benchmarked outputs.
- Any figure driven by `generate_all_figures.py` must be tagged internally as
  `computed`, `derived`, or `schematic` before it is allowed back into the docs.
- Quantum-wire material stays provisional until sparse assembly, state
  selection, strain mapping, and `g`-factor trends all have external anchors.

## First Closure Targets

| Order | Target | Reason |
|---|---|---|
| `1` | `FL-001`, `FL-002`, `FL-003`, `FL-004` | These are clear falsehoods or label/index mismatches with straightforward ownership. |
| `2` | `FL-005`, `FL-006`, `FL-007`, `FL-009` | These block the entire optics/QCSE/scattering narrative and likely expose deeper code issues. |
| `3` | `FL-012`, `FL-013`, `FL-014`, `FL-015` | This is the dedicated quantum-wire hardening stream. |
| `4` | `FL-010`, `FL-011` | Docs completeness issues after the physics layer is stable. |
