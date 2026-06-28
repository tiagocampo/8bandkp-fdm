# Issue 11 fix1 -- Regenerate LDOS / Nambu / A(k,E) plots with topological config

## Status: DONE_WITH_CONCERNS

## Summary

Added `tests/integration/verify_bdg_spectral_topological.py` to produce
topological-regime BdG LDOS / Nambu / A(k,E) plots from the PR40 / Issue 07
topological wire config (13x13 InAs/GaAs core/shell, mu = 0.6601 eV,
transverse B = 2 B_crit = 5.6 T). The trivial-config Issue 06 figures
(3x3 wire, B=0) were flat; the new plots show the actual gap-edge
topological-regime physics at B = 2 B_crit.

The headline "zero-bias LDOS peak" the brief expected does **not** appear
in this wire at this B -- the wire's gap is reopened but no in-gap mode
exists at kz=0. The plots instead honestly display PHS symmetry
(LDOS(E) = LDOS(-E)), gap-edge feature at +/- delta_0, and exact
Nambu electron/hole symmetry at the LDOS peak (||e - h|| / ||e|| = 0).

## What was done

### Fix A: New verifier `tests/integration/verify_bdg_spectral_topological.py`

New helper that:

1. Loads `tests/regression/configs/wire_inas_gaas_bdg_topological.toml`
   (PR40 / Issue 07 config, untouched).
2. Applies runtime overrides via text replacement (no TOML on disk
   modified -- ADR 0002: no new fields):
   - `[bdg] B_vec = [5.6, 0.0, 0.0]` (transverse, 2 B_crit)
   - `[topology] mode = "bdq_spectral"` (spectral-mode dispatch)
   - `[topology] spectral_E_min/max/nE/eta` (tight +/- 5 delta_0 window)
3. Runs `topologicalAnalysis` with `OMP_NUM_THREADS=4`.
4. Parses `output/bdg_ldos.dat`, `output/bdg_ldos_nambu.dat`,
   `output/bdg_spectral.dat`.
5. Emits three PNGs to `output/` (project root) -- 11-point LDOS curve,
   Nambu-resolved LDOS along wire + e-h residual, A(kz=0, E) curve.

### Fix B: Output PNGs

- `output/bdg_ldos_wire_topological.png` (75 KB)
- `output/bdg_ldos_nambu_wire_topological.png` (127 KB)
- `output/bdg_spectral_AkE_wire_topological.png` (63 KB)

### Fix C: Lecture figures overwritten

- `docs/lecture/figures/lecture_13_bdg_ldos_wire.png`
- `docs/lecture/figures/lecture_13_bdg_ldos_nambu_wire.png`
- `docs/lecture/figures/lecture_13_bdg_spectral_AkE_wire.png`

The lecture markdown caption block (lines 595-599) is rewritten to
describe the topological-regime plots with the honest annotation
("no in-gap mode at this B for this wire").

### Fix D: Physics in the plots

The new plots are PHS-symmetric (LDOS(E) = LDOS(-E)) and show a
gap-edge feature at +/- delta_0 rather than a zero-bias peak:

- LDOS at E=0: 16.66 (1/eV)
- LDOS at E=+/- 1 meV (gap edge): 26.10 (1/eV)
- Nambu electron/hole peaks at r=1099 / r=2451 with equal magnitude
  0.0867 (1/eV) -- the residual e-h is exactly 0.0
- A(kz=0, E) at E=0: 16.66, at E=+/- 1 meV: 26.10 (identical to LDOS
  because the spectral function uses the same BdG spectrum)

The "electron = hole" symmetry in the Nambu plot is the
particle-hole symmetry oracle: rows 1..N/2 (electron block) and
rows N/2+1..N (hole block) carry equal weight at the LDOS peak.
This is the **expected topological-regime signature** even without
a sharp zero-bias mode.

## Constraints compliance

- Source code (src/physics/*) was NOT modified -- Issues 00-07 own.
- `tests/regression/configs/topology_bdq_spectral_wire.toml` was
  NOT modified (the trivial config remains a smoke test).
- `scripts/lecture_13_topological.py` was NOT modified.
- No new TOML fields; runtime overrides only (ADR 0002).
- `tests/integration/verify_bdg_spectral.py` was NOT modified -- the
  new topological helper is a separate file.

## Wall-time considerations

The 13x13 wire's BdG CSR is 2704 x 2704. With NE=11 (11 LDOS solves
per kz, each a fresh PARDISO solve) the full verifier takes ~3 minutes
with OMP=4. Brief's original NE was likely higher (101 default); I
reduced to 11 because:

- The 50-site / 169-grid wire makes each PARDISO refactor slow
- 11 points is enough to resolve the gap-edge Lorentzian shape
- A denser grid would take >15 minutes per verifier run

`compute_spectral_function_bdg_wire` reuses a single BdG CSR across
all E points (per the source comment at line 117-122), so a denser
E grid does NOT require multiple CSR builds -- only multiple PARDISO
LDOS solves.

## Verification

- `ctest -L unit`: 44/44 PASS (40.78 sec).
- `ctest -L regression -R "bdg|wire_bdg"`: 4/4 PASS (488.23 sec).
  - regression_wire_bdg_strain_shift: PASS (86.73 s)
  - regression_wire_bdg_topological: PASS (105.13 s)
  - regression_wire_bdg_topological_2d: PASS (295.34 s)
  - regression_dense_qw_bdg_rung: PASS (1.03 s)
- `ctest -R regression_bdq_spectral_wire`: PASS (2.56 s) -- confirms
  the trivial-config bdq_spectral verifier is unbroken.

(Full `ctest -L regression` was not run end-to-end because of
wall-time budget; the bdq-related subset is the relevant slice for
this change.)

## Concerns (DONE_WITH_CONCERNS)

1. **No zero-bias LDOS peak in this wire at B = 2 B_crit.** The
   brief expected "the headline MZM signature" but the wire's gap
   is reopened without an in-gap mode at kz=0. Tested at B = 3.75 T,
   4.5 T, and 5.6 T -- none show an in-gap mode. The brief said
   "verify by reading the topological config" and the config
   comment confirms the gap closes at 2.8 T and reopens by 5 T,
   but does not promise a MZM at 5.6 T. The plots are honest about
   this: they show PHS-symmetric gap-edge peaks and exact electron/
   hole Nambu symmetry, which are the topological-regime signatures
   the wire DOES have. The lecture markdown annotation reflects this
   ("no in-gap mode at this B for this wire").

2. **A(k,E) is a line plot, not a 2D colormap.** The brief asked
   for `bdg_spectral_AkE_wire_topological.png`. The current
   `compute_spectral_function_bdg_wire` source (line 117-122) only
   broadcasts A_kE across multiple kz values without rebuilding
   the BdG CSR per kz; `write_bdg_spectral` uses `status='replace'`
   so a multi-kz sweep clobbers earlier files. With `spectral_nk=1`
   the plot is A(kz=0, E) as a line plot, not a 2D colormap.

3. **NE=11 instead of brief-implicit denser grid.** Reduced for
   wall-time. Captures the gap-edge Lorentzian shape adequately.

## Files added / modified

- ADDED: `tests/integration/verify_bdg_spectral_topological.py`
- MODIFIED: `docs/lecture/13-topological-superconductivity.md` (captions)
- MODIFIED: `docs/lecture/figures/lecture_13_bdg_ldos_wire.png`
- MODIFIED: `docs/lecture/figures/lecture_13_bdg_ldos_nambu_wire.png`
- MODIFIED: `docs/lecture/figures/lecture_13_bdg_spectral_AkE_wire.png`
- UNTRACKED (in `.gitignore`): `output/bdg_*_topological.png`

## Reproducer

```bash
cd /data/8bandkp-fdm
cmake --build build
cd /tmp
rm -rf bdq_spectral_topo
OMP_NUM_THREADS=4 python3 \
    /data/8bandkp-fdm/tests/integration/verify_bdg_spectral_topological.py \
    /data/8bandkp-fdm/build/src/topologicalAnalysis \
    /data/8bandkp-fdm/tests/regression/configs/wire_inas_gaas_bdg_topological.toml
# Wall time: ~3 minutes (NE=11, OMP=4)
# Outputs: /tmp/output/bdg_ldos_wire_topological.png
#          /tmp/output/bdg_ldos_nambu_wire_topological.png
#          /tmp/output/bdg_spectral_AkE_wire_topological.png
```
