# Docs Physics Revamp Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Rebuild the lecture suite as research-grade internal documentation by fixing broken physics/code paths, reproducing the maximum practical subset of relevant nextnano examples, regenerating trustworthy figures, and rewriting the chapters around validated evidence.

**Architecture:** Execute the revamp in five phases: control documents, known-false figure triage, benchmark-driven solver repair, quantum-wire hardening, and chapter rebuild plus regressionization. The critical path runs through optics/QCSE correctness and the sparse wire stack; documentation updates happen only after the underlying outputs are trustworthy.

**Tech Stack:** Fortran 90, Python 3, matplotlib, numpy, shell integration tests, pFUnit unit tests, CMake/CTest, Markdown docs.

---

## Task 1: Create the benchmark matrix and failure ledger

**Files:**
- Create: `docs/plans/2026-04-21-benchmark-matrix.md`
- Create: `docs/plans/2026-04-21-failure-ledger.md`
- Reference: `docs/lecture/01-bulk-band-structure.md`
- Reference: `docs/lecture/02-quantum-well.md`
- Reference: `docs/lecture/04-strain.md`
- Reference: `docs/lecture/06-optical-properties.md`
- Reference: `docs/lecture/07-self-consistent-sp.md`
- Reference: `docs/lecture/08-quantum-wire.md`
- Reference: `docs/lecture/10-qcse.md`

**Step 1: Write the benchmark matrix**

Document, per chapter:

- internal example/config,
- external benchmark source,
- observables,
- pass/fail tolerance,
- current status.

**Step 2: Write the failure ledger**

Record the known issues already identified:

- broken units,
- bad state indexing,
- parser errors,
- missing embeds,
- suspect solver behavior,
- wire-risk items.

**Step 3: Verify the files exist**

Run:
```bash
ls -la docs/plans/2026-04-21-benchmark-matrix.md docs/plans/2026-04-21-failure-ledger.md
```

Expected: both files listed.

**Step 4: Commit**

```bash
git add docs/plans/2026-04-21-benchmark-matrix.md docs/plans/2026-04-21-failure-ledger.md
git commit -m "docs: add benchmark matrix and revamp failure ledger"
```

---

## Task 2: Add figure provenance inventory

**Files:**
- Create: `docs/plans/2026-04-21-figure-provenance.md`
- Modify: `scripts/plotting/generate_all_figures.py`

**Step 1: Write the provenance inventory**

For each figure in `docs/figures/`, record:

- generator function name,
- config(s),
- output file(s),
- figure type (`computed`, `derived`, `schematic`),
- trust level (`validated`, `needs benchmark`, `remove/downgrade`).

**Step 2: Mark schematic outputs explicitly**

In `scripts/plotting/generate_all_figures.py`, update docstrings or comments for
all non-computed figures so they cannot be mistaken for validated output.

**Step 3: Verify no syntax errors were introduced**

Run:
```bash
python3 -m py_compile scripts/plotting/generate_all_figures.py
```

Expected: no output.

**Step 4: Commit**

```bash
git add docs/plans/2026-04-21-figure-provenance.md scripts/plotting/generate_all_figures.py
git commit -m "docs: inventory figure provenance and trust level"
```

---

## Task 3: Fix Zeeman figure correctness and g-factor presentation

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py`
- Modify: `docs/lecture/05-gfactor.md`
- Test: `tests/integration/test_gfactor.sh`

**Step 1: Write the failing acceptance check**

Document an acceptance check in the task branch: the slope in
`docs/figures/gfactor_zeeman.png` must match `mu_B * g / 2` in meV/T, not `g/2`.

**Step 2: Run the current integration test**

Run:
```bash
ctest --output-on-failure -R gfactor
```

Expected: current tests pass or expose nearby issues, but they do not yet catch
the bad plotting semantics.

**Step 3: Implement the minimal fix**

Update `fig_gfactor_zeeman` so the plotted energies use the correct Bohr
magneton conversion and label the figure as `computed` rather than `schematic`
only if the output is actually quantitative. If it remains schematic, relabel it
as schematic in the chapter.

**Step 4: Regenerate the figure**

Run:
```bash
python3 scripts/plotting/generate_all_figures.py gfactor_zeeman
```

Expected: `docs/figures/gfactor_zeeman.png` is regenerated with physically
consistent slopes.

**Step 5: Update chapter prose**

Adjust `docs/lecture/05-gfactor.md` so the text matches the new figure semantics.

**Step 6: Re-run verification**

Run:
```bash
ctest --output-on-failure -R gfactor
```

Expected: pass.

**Step 7: Commit**

```bash
git add scripts/plotting/generate_all_figures.py docs/lecture/05-gfactor.md docs/figures/gfactor_zeeman.png
git commit -m "fix: correct Zeeman figure units and chapter framing"
```

---

## Task 4: Fix QCSE field units, labeling, and Stark validation

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py`
- Modify: `docs/lecture/10-qcse.md`
- Modify: `tests/integration/test_qcse_with_ef.sh`
- Modify: `tests/integration/verify_stark_shift.py`
- Reference: `tests/regression/configs/sc_qcse_gaas_algaas_ef.cfg`
- Reference: `tests/regression/configs/qw_gaas_algaas_qcse_scattering.cfg`

**Step 1: Write the failing test/update target**

Extend the QCSE verification path so it checks the field conversion and the
reported field label against the config value.

**Step 2: Run the current QCSE checks**

Run:
```bash
ctest --output-on-failure -R qcse
```

Expected: current checks pass or expose gaps, but they do not yet catch the
`-70` vs `-700 kV/cm` mismatch.

**Step 3: Implement the minimal code and plot fix**

Correct the field conversion comments and logic in
`scripts/plotting/generate_all_figures.py`, and update the plot labels and
chapter prose so they describe the actual simulated field.

**Step 4: Regenerate the QCSE figure**

Run:
```bash
python3 scripts/plotting/generate_all_figures.py qcse_stark_shift
```

Expected: the figure title, legend, and chapter explanation agree with the
actual config.

**Step 5: Re-run QCSE verification**

Run:
```bash
ctest --output-on-failure -R qcse
```

Expected: pass with the stronger field/unit check in place.

**Step 6: Commit**

```bash
git add scripts/plotting/generate_all_figures.py docs/lecture/10-qcse.md tests/integration/test_qcse_with_ef.sh tests/integration/verify_stark_shift.py docs/figures/qcse_stark_shift.png
git commit -m "fix: align QCSE field units, labels, and verification"
```

---

## Task 5: Fix QW state indexing and per-band density semantics

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py`
- Modify: `docs/lecture/03-wavefunctions.md`
- Test: `tests/integration/test_qw_bandstructure.sh`
- Test: `tests/unit/test_optical_qw.pf`

**Step 1: Write the failing acceptance note**

The first wavefunction figure in Chapter 3 must either show `states 1-4` or
load the conduction states it claims to show. The per-band CB1 figure must
visually match its caption.

**Step 2: Run the current QW checks**

Run:
```bash
ctest --output-on-failure -R "qw_bandstructure|optical_qw"
```

Expected: existing tests do not yet validate figure-state mapping.

**Step 3: Implement the minimal fix**

Update the generator to:

- load the actual conduction states for `qw_wavefunctions.png`, or relabel the
  figure and chapter if that is the intended output,
- verify the selected state for `perband_density.png`,
- correct captions and explanatory prose in `docs/lecture/03-wavefunctions.md`.

**Step 4: Regenerate the affected figures**

Run:
```bash
python3 scripts/plotting/generate_all_figures.py qw_wavefunctions perband_density
```

Expected: both figures regenerate and the visual semantics match the chapter.

**Step 5: Add or tighten verification**

Where practical, add a small automated check that the selected eigenstate index
used by the generator matches the captioned state.

**Step 6: Re-run checks**

Run:
```bash
ctest --output-on-failure -R "qw_bandstructure|optical_qw"
```

Expected: pass.

**Step 7: Commit**

```bash
git add scripts/plotting/generate_all_figures.py docs/lecture/03-wavefunctions.md docs/figures/qw_wavefunctions.png docs/figures/perband_density.png
git commit -m "fix: align QW wavefunction figures with actual state selection"
```

---

## Task 6: Repair optics, exciton, and scattering figure trustworthiness

**Files:**
- Modify: `src/physics/optical_spectra.f90`
- Modify: `src/physics/scattering.f90`
- Modify: `src/physics/exciton.f90`
- Modify: `scripts/plotting/generate_all_figures.py`
- Modify: `docs/lecture/06-optical-properties.md`
- Test: `tests/unit/test_optical.pf`
- Test: `tests/unit/test_optical_qw.pf`

**Step 1: Write failing benchmark targets**

Define benchmark acceptance checks for:

- TE/TM contrast near the HH edge,
- exciton peak location vs computed binding energy,
- scattering lifetime extraction from `scattering_rates.dat`.

**Step 2: Run current optical tests**

Run:
```bash
ctest --output-on-failure -R "optical|optical_qw"
```

Expected: current unit tests pass or expose nearby issues, but do not yet prove
the plotted Chapter 6 behavior.

**Step 3: Diagnose whether each failure is solver-side or parser-side**

For each broken figure, decide whether the defect is in:

- `src/physics/*` output,
- output formatting,
- or `generate_all_figures.py`.

**Step 4: Implement minimal repairs**

Fix the responsible code paths, including:

- TE/TM semantics,
- exciton binding usage in plotted spectra,
- scattering-rate column interpretation and lifetime extraction.

**Step 5: Regenerate Chapter 6 figures**

Run:
```bash
python3 scripts/plotting/generate_all_figures.py qw_absorption_spectrum qw_absorption_strained qw_absorption_vs_width absorption_with_exciton absorption_excitonic_TE scattering_lifetime_vs_width scattering_lifetime_vs_field
```

Expected: the regenerated figures support the chapter claims or force those
claims to be weakened.

**Step 6: Rewrite Chapter 6 around the validated outputs**

Update `docs/lecture/06-optical-properties.md` so every retained claim is backed
by the repaired figures and benchmark evidence.

**Step 7: Re-run verification**

Run:
```bash
ctest --output-on-failure -R "optical|optical_qw|qcse"
```

Expected: pass.

**Step 8: Commit**

```bash
git add src/physics/optical_spectra.f90 src/physics/scattering.f90 src/physics/exciton.f90 scripts/plotting/generate_all_figures.py docs/lecture/06-optical-properties.md docs/figures/
git commit -m "fix: harden optics, exciton, and scattering outputs"
```

---

## Task 7: Reproduce and lock core bulk, QW, and strain benchmarks

**Files:**
- Modify: `tests/integration/test_bulk_bandstructure.sh`
- Modify: `tests/integration/test_qw_bandstructure.sh`
- Modify: `tests/unit/test_hamiltonian.pf`
- Modify: `tests/unit/test_hamiltonian_2d.pf`
- Modify: `docs/lecture/01-bulk-band-structure.md`
- Modify: `docs/lecture/02-quantum-well.md`
- Modify: `docs/lecture/04-strain.md`

**Step 1: Encode the benchmark observables**

Add checks for the specific observables that map to the benchmark matrix:

- bulk GaAs gap and SO split,
- subband ordering in the reference QWs,
- strain-induced HH/LH ordering and sign conventions.

**Step 2: Run the current checks**

Run:
```bash
ctest --output-on-failure -R "bulk_bandstructure|qw_bandstructure|strain"
```

Expected: passing or failing behavior becomes the baseline for this task.

**Step 3: Repair any benchmark misses**

Touch only the code paths needed to align the outputs with the benchmark matrix.

**Step 4: Update the benchmark-backed chapter sections**

Revise the benchmark sections in Chapters `01`, `02`, and `04` to use the
actual validated numbers and references.

**Step 5: Re-run checks**

Run:
```bash
ctest --output-on-failure -R "bulk_bandstructure|qw_bandstructure|strain"
```

Expected: pass.

**Step 6: Commit**

```bash
git add tests/integration/test_bulk_bandstructure.sh tests/integration/test_qw_bandstructure.sh tests/unit/test_hamiltonian.pf tests/unit/test_hamiltonian_2d.pf docs/lecture/01-bulk-band-structure.md docs/lecture/02-quantum-well.md docs/lecture/04-strain.md
git commit -m "test: lock core bulk QW and strain benchmarks"
```

---

## Task 8: Reproduce and harden self-consistent SP and QCSE benchmarks

**Files:**
- Modify: `tests/integration/test_sc_loop.sh`
- Modify: `tests/integration/test_sc_mod_doped.sh`
- Modify: `tests/integration/test_sc_qw_inas_alsb.sh`
- Modify: `tests/unit/test_sc_loop.pf`
- Modify: `docs/lecture/07-self-consistent-sp.md`
- Modify: `docs/lecture/10-qcse.md`

**Step 1: Add benchmark observables**

Check the key SP/QCSE observables identified in the matrix:

- convergence behavior,
- sign and shape of band bending,
- subband shifts under field,
- consistency between chapter tables and actual outputs.

**Step 2: Run the current SP/QCSE checks**

Run:
```bash
ctest --output-on-failure -R "sc_|qcse"
```

Expected: passing or failing baseline established.

**Step 3: Repair code or prose mismatches**

When a benchmark misses, fix the solver or the chapter, not just the caption.

**Step 4: Embed missing SC figures**

Replace path bullets in `docs/lecture/07-self-consistent-sp.md` with actual
figure embeds where the outputs are already trusted.

**Step 5: Re-run checks**

Run:
```bash
ctest --output-on-failure -R "sc_|qcse"
```

Expected: pass.

**Step 6: Commit**

```bash
git add tests/integration/test_sc_loop.sh tests/integration/test_sc_mod_doped.sh tests/integration/test_sc_qw_inas_alsb.sh tests/unit/test_sc_loop.pf docs/lecture/07-self-consistent-sp.md docs/lecture/10-qcse.md
git commit -m "fix: harden SP and QCSE benchmarks and docs"
```

---

## Task 9: Harden the quantum-wire solver and benchmark suite

**Files:**
- Modify: `src/math/geometry.f90`
- Modify: `src/math/sparse_matrices.f90`
- Modify: `src/math/eigensolver.f90`
- Modify: `src/physics/strain_solver.f90`
- Modify: `src/physics/gfactor_functions.f90`
- Modify: `src/io/outputFunctions.f90`
- Modify: `tests/integration/test_wire_bandstructure.sh`
- Modify: `tests/integration/test_wire_core_shell.sh`
- Modify: `tests/integration/test_wire_insb_gfactor.sh`
- Modify: `tests/unit/test_geometry.pf`
- Modify: `tests/unit/test_hamiltonian_2d.pf`
- Create or update: `tests/regression/data/`
- Modify: `docs/lecture/08-quantum-wire.md`

**Step 1: Write the wire benchmark checklist**

Define what must be true before strong wire claims survive:

- geometry/material masks are correct,
- sparse Hamiltonian assembly is stable,
- eigenpairs are consistent,
- wire strain maps correctly,
- wire g-factor and optical outputs use the correct states and axes.

**Step 2: Run the current wire checks**

Run:
```bash
ctest --output-on-failure -R "wire_"
```

Expected: current wire stability and coverage baseline is established.

**Step 3: Repair the sparse path**

Fix the code paths indicated by the benchmark failures. Focus on correctness,
not feature expansion.

**Step 4: Add or refresh regression data**

Update `tests/regression/data/` with wire outputs that represent repaired
behavior and can be compared in future runs.

**Step 5: Rewrite Chapter 8 conservatively**

Remove or downgrade any wire claims that remain under-validated after the code
repairs. Keep only what the sparse path can defend.

**Step 6: Re-run wire verification**

Run:
```bash
ctest --output-on-failure -R "wire_"
```

Expected: pass.

**Step 7: Commit**

```bash
git add src/math/geometry.f90 src/math/sparse_matrices.f90 src/math/eigensolver.f90 src/physics/strain_solver.f90 src/physics/gfactor_functions.f90 src/io/outputFunctions.f90 tests/integration/test_wire_bandstructure.sh tests/integration/test_wire_core_shell.sh tests/integration/test_wire_insb_gfactor.sh tests/unit/test_geometry.pf tests/unit/test_hamiltonian_2d.pf tests/regression/data docs/lecture/08-quantum-wire.md
git commit -m "fix: harden quantum wire solver path and benchmarks"
```

---

## Task 10: Rebuild the lecture chapters around validated evidence

**Files:**
- Modify: `docs/lecture/00-quickstart.md`
- Modify: `docs/lecture/01-bulk-band-structure.md`
- Modify: `docs/lecture/02-quantum-well.md`
- Modify: `docs/lecture/03-wavefunctions.md`
- Modify: `docs/lecture/04-strain.md`
- Modify: `docs/lecture/05-gfactor.md`
- Modify: `docs/lecture/06-optical-properties.md`
- Modify: `docs/lecture/07-self-consistent-sp.md`
- Modify: `docs/lecture/08-quantum-wire.md`
- Modify: `docs/lecture/09-numerical-methods.md`
- Modify: `docs/lecture/10-qcse.md`
- Modify: `docs/lecture/11-convergence.md`
- Modify: `docs/lecture/12-extending-the-code.md`

**Step 1: Apply the chapter template**

For each chapter, ensure the retained content is organized as:

- validated example,
- relevant config or output path,
- quantitative result,
- benchmark/reference,
- limitations.

**Step 2: Remove weak or unsupported claims**

Delete or downgrade sections that still fail the acceptance bar.

**Step 3: Upgrade thin chapters**

Bring Chapters `05`, `10`, and `12` up to the new standard with either stronger
evidence or tighter, more honest scope.

**Step 4: Verify link and image integrity**

Run:
```bash
rg -n "!\[|Figure|nextnano|benchmark|TODO|FIXME" docs/lecture
```

Expected: no stale TODO/FIXME markers in finished content; figure references
match the rebuilt chapters.

**Step 5: Commit**

```bash
git add docs/lecture
git commit -m "docs: rebuild lecture chapters around validated outputs"
```

---

## Task 11: Add final validation summary and revamp closeout report

**Files:**
- Create: `docs/plans/2026-04-21-validation-summary.md`
- Create: `docs/plans/2026-04-21-wire-validation-summary.md`

**Step 1: Write the global validation summary**

Summarize:

- which nextnano examples were reproduced,
- which observables were matched,
- which chapters are now benchmark-backed,
- what remains intentionally downgraded or removed.

**Step 2: Write the wire validation summary**

Document the final trust envelope for the sparse wire path, including any
remaining limits.

**Step 3: Verify files exist**

Run:
```bash
ls -la docs/plans/2026-04-21-validation-summary.md docs/plans/2026-04-21-wire-validation-summary.md
```

Expected: both files listed.

**Step 4: Commit**

```bash
git add docs/plans/2026-04-21-validation-summary.md docs/plans/2026-04-21-wire-validation-summary.md
git commit -m "docs: add validation closeout summaries for revamp"
```

---

## Task 12: Final verification

**Files:**
- Verify: entire worktree

**Step 1: Build from a clean state**

Run:
```bash
cmake --build build
```

Expected: successful build.

**Step 2: Run focused tests**

Run:
```bash
ctest --output-on-failure -R "bulk_bandstructure|qw_bandstructure|gfactor|qcse|sc_|wire_|optical"
```

Expected: all targeted physics tests pass.

**Step 3: Regenerate the full figure suite**

Run:
```bash
python3 scripts/plotting/generate_all_figures.py all
```

Expected: all retained figures regenerate successfully.

**Step 4: Inspect git diff**

Run:
```bash
git status --short
```

Expected: only intentional revamp changes remain.

**Step 5: Final commit**

```bash
git add -A
git commit -m "docs: complete physics-backed lecture revamp"
```
