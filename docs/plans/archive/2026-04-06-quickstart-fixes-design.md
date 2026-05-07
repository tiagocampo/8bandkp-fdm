# Quickstart Guide Fixes — Design

**Date:** 2026-04-06
**File:** `docs/lecture/00-quickstart.md`
**Branch:** `feature/docs-overhaul`

## Problem

Two issues in the quickstart guide:

1. **Lines 65-66** dismiss `gfactorCalculation` with "You can safely ignore... we will
   only use `bandStructure`." This is wrong — the quickstart should introduce both
   executables and walk the reader through each one.
2. **Lines 257-260** describe a `cp -i` shell alias issue that only applies to the
   author's specific setup, not to general readers. Should be removed.

## Changes

### 1. Fix executable introduction (Section 2, lines 58-66)

**Before:**

```
When the build succeeds you will find two executables:

    build/src/bandStructure        # band-structure sweeps
    build/src/gfactorCalculation   # Landau g-factor at Gamma

You can safely ignore the `gfactorCalculation` executable for now; we will
only use `bandStructure` in this guide.
```

**After:**

```
When the build succeeds you will find two executables:

    build/src/bandStructure        # band-structure E(k) sweeps
    build/src/gfactorCalculation   # Landau g-factor at Gamma

Both read the same `input.cfg` format. We will use `bandStructure` in
Section 3 and `gfactorCalculation` in Section 4.
```

### 2. Add Section 4 — "Second Run: Bulk GaAs g-Factor"

Insert between current Section 3 (First Run) and Section 4 (Interpreting the
Output). Current sections 4-6 renumber to 5-7.

Content outline:

- New input.cfg for bulk GaAs CB g-factor (based on
  `tests/regression/configs/gfactor_bulk_gaas_cb.cfg`)
- Key parameter differences from the band-structure run:
  - `waveVector: k0` (Gamma point only)
  - `waveVectorStep: 0` (single point)
  - `whichBand: 0` (conduction band)
  - `bandIdx: 1` (first subband)
- Run command: `./build/src/gfactorCalculation`
- Expected stdout: parameter echo, spin matrices, 2x2 g-tensor per direction,
  g-factor eigenvalues for gx, gy, gz
- Expected `output/gfactor.dat`: three isotropic values ≈ -0.315
- Brief physics note: isotropic because bulk zincblende has cubic symmetry;
  the 8-band value (-0.315) excludes remote-band contributions (experiment: -0.44)

Input.cfg for the worked example:

```
waveVector: k0
waveVectorMax: 0.1
waveVectorStep: 0
confinement:  0
FDstep: 1
FDorder: 2
numLayers:  1
material1: GaAs
numcb: 2
numvb: 6
ExternalField: 0  EF
EFParams: 0.0005
whichBand: 0
bandIdx: 1
```

### 3. Remove `cp -i` entry from Common Issues (Section 7)

Delete lines 257-260:

```markdown
### `cp -i` shell alias

If your shell aliases `cp` to `cp -i` (interactive), copying `input.cfg` via a
script will hang waiting for confirmation. Use your editor or write the file
directly instead of `cp`.
```

### 4. Update section numbers and Next Steps table

- Old Section 4 → Section 5 (Interpreting the Output)
- Old Section 5 → Section 6 (Next Steps)
- Old Section 6 → Section 7 (Common Issues)
- Update the "jump straight to Section 6" link in the intro to "Section 7"
- Update the Next Steps table to include the new g-factor chapter reference
  alongside the existing Chapter 05 link
