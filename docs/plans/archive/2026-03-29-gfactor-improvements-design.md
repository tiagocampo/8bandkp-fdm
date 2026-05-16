# G-Factor Calculation Improvements Design

Date: 2026-03-29

## Problem Statement

The g-factor calculation (`gfactorCalculation`) has two limitations:
1. **CB subband truncation**: `main_gfactor.f90` hard-codes `numcb=2`, so only 2 CB eigenvectors are The CB-CB perturbation sum (remote conduction band contributions) never executes for QW systems.
2. **No hole g-factor**: Only `whichBand == 0` (CB) is Valence band g-factors are important for spin physics, are not implemented.

## Design

### 1. Fix CB Subband Truncation

**Problem:** Lines 78 and 186 of `main_gfactor.f90` set `cfg%numcb = 2`, extracting only 2 CB eigenvectors. The CB-CB perturbation loop in `gfactorCalculation` (line 447) is guarded by `if (numcb > 2)` and never fires.

**Fix in `main_gfactor.f90`:**
- Line 78: Change `cfg%numcb = 2` to `cfg%numcb = 2*cfg%fdStep`
- Line 186: Same change
- Line 85: The assertion `evnum /= N` still holds: `2*fdStep + 6*fdStep = 8*fdStep = N`

**Fix in `gfactorCalculation` (`gfactor_functions.f90`):**
- CB-CB loop (line 454): Add self-term exclusion:
  ```fortran
  ! Skip self-interaction (l==n or l==m gives zero contribution)
  if (l == n .or. l == m) cycle
  ```
- The guard `if (numcb > 2)` can be relaxed to always execute; for bulk `numcb=2` and the loop body skips both values (correct behavior).

**Bulk correctness:** For `nlayers==1`, the matrix is 8x8. With `numcb=2` the CB-CB sum has only 2 entries, both skipped. The Roth formula only involves VB intermediate states. No change needed.

### 2. Hole G-Factor (`whichBand == 1`)

**Theory:** Same Lowdin partitioning formula, but with:
- **Doublet:** Two consecutive VB states (Kramers partners) at `bandIdx` and `bandIdx+1`
- **Intermediate states:** All CB states + all VB states except the doublet
- **Energy denominators:** `denom = (E_vb(n) - E_int(l)) + (E_vb(m) - E_int(l))` — automatically handles signs

**State ordering:** `vb_state(:,1)` = highest VB energy (topmost), `vb_state(:,numvb)` = lowest. `cb_state(:,1)` = lowest CB, `cb_state(:,numcb)` = highest.

**Implementation -- add `else if (whichBand == 1)` block in `gfactorCalculation`:**

1. **Sigma matrix:** Evaluate between VB doublet states:
   ```fortran
   sigma(ii,jj,:) = sigmaElem(vb_state(:,n), vb_state(:,m), ...)
   ```
   where `n` and `m` iterate over `bandIdx` and `bandIdx+1`.

2. **CB intermediate sum** over `l=1,numcb`:
   ```fortran
   Pele1 = <vb_n | H(mod1) | cb_l>
   Pele2 = <cb_l | H(mod2) | vb_m>
   Pele3 = <vb_n | H(mod2) | cb_l>
   Pele4 = <cb_l | H(mod1) | vb_m>
   denom = (vb_value(n) - cb_value(l)) + (vb_value(m) - cb_value(l))
   ```
   Note: `cb_value(l) > vb_value(n)`, so denom is always negative (physically correct).

3. **VB intermediate sum** over `l=1,numvb`, skipping `l==n` and `l==m`:
   ```fortran
   if (l == n .or. l == m) cycle
   Pele1 = <vb_n | H(mod1) | vb_l>
   Pele2 = <vb_l | H(mod2) | vb_m>
   Pele3 = <vb_n | H(mod2) | vb_l>
   Pele4 = <vb_l | H(mod1) | vb_m>
   denom = (vb_value(n) - vb_value(l)) + (vb_value(m) - vb_value(l))
   ```

4. **Final rescaling** (lines 491-492): Same as CB -- `-i/hbar2O2m0` and `g0/2 * sigma` already applied at the end.

### 3. Input & Config Changes

**Add to `simulation_config` type (`defs.f90`):**
```fortran
integer :: whichBand = 0     ! 0=CB g-factor, 1=VB g-factor
integer :: bandIdx = 1       ! which doublet (1=ground, 2=first excited, ...)
```

**Add to `input_parser.f90`:**
After the ExternalField/EFParams block, read optional gfactor params with fallback:
```fortran
! Try reading gfactor params; use defaults if missing (backward compatible)
read(data_unit, *, iostat=status) label, cfg%whichBand
if (status == 0) then
  print *, trim(label), cfg%whichBand
  read(data_unit, *, iostat=status) label, cfg%bandIdx
  if (status == 0) print *, trim(label), cfg%bandIdx
else
  cfg%whichBand = 0
  cfg%bandIdx = 1
end if
```

**`main_gfactor.f90` changes:**
- Remove hardcoded `whichBand = 0` and `bandIdx = 1` (lines 184-185)
- Use `cfg%whichBand` and `cfg%bandIdx` instead
- Keep `numcb = 2*fdStep` override (gfactor needs all states)

**Updated `gfactor.example`:**
```
waveVector: k0
waveVectorMax: 0.1
waveVectorStep: 0
confinement:  1
FDstep: 101
numLayers:  3
material1: AlSb -250  250 0
material2: GaSb -135  135 0.2414
material3: InAs  -35   35 -0.0914
numcb: 32
numvb: 32
ExternalField: 0  EF
EFParams: 0.0005
whichBand: 0
bandIdx: 1
```

## Files to Modify

| File | Change |
|------|--------|
| `src/core/defs.f90` | Add `whichBand`, `bandIdx` to `simulation_config` |
| `src/io/input_parser.f90` | Read `whichBand`, `bandIdx` with fallback |
| `src/apps/main_gfactor.f90` | Use `cfg%whichBand`/`bandIdx`, set `numcb=2*fdStep` |
| `src/physics/gfactor_functions.f90` | Add `whichBand==1` block with VB g-factor, fix CB-CB self-term skip |
| `gfactor.example` | Add `whichBand` and `bandIdx` fields |
