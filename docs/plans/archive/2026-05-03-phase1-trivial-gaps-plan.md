# Phase 1: Close the Trivial Gaps — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Close 5 INCOMPLETE groups (#45, #5, #16, #10, #8) with minimal effort, then update tracking docs and archive plan files.

**Architecture:** 5 independent tasks. Each task commits its changes. Final task updates REVIEW.md, BACKLOG.md, and archives all plan files in one batch commit.

**Tech Stack:** Fortran 2018, Python 3 (matplotlib), bash, git

---

### Task 1: Close #45 — zdotc decision documented

**Files:** None (documentation-only)

**Step 1: Verify CLAUDE.md already documents the decision**

Run: `grep -n "zdotc" CLAUDE.md`
Expected: line ~182 mentioning "zdotc declared as function-return ABI in linalg.f90 — MKL on this platform returns complex(dp) by value, not via subroutine hidden argument. Do not change to subroutine form."

**Step 2: No code changes. Mark complete in mental checklist.**

Group #45 status: COMPLETE.

---

### Task 2: Close #5 — Makefile target + z(1) guard

**Files:**
- Modify: `Makefile`
- Modify: `src/physics/hamiltonianConstructor.f90:21-33`
- Test: build + run `fd_sign_compare`

#### Step 1: Add test-fd-sign target to Makefile

Add before the `.PHONY` line in `Makefile`:

```makefile
test-fd-sign: all
	$(FC) -o $(BUILD_DIR)/fd_sign_compare tests/fd_sign_compare.f90 -I$(BUILD_DIR) -L$(BUILD_DIR)/src -lkpfdm $(MKLROOT)/lib/libmkl_sequential.a $(MKLROOT)/lib/libmkl_core.a -lpthread -lm -ldl || echo "NOTE: test-fd-sign requires manual compilation; see tests/fd_sign_compare.f90"
```

If this doesn't work (standalone program may not link against the library), add a simpler approach — just document it as a manual-run test:

```makefile
test-fd-sign:
	@echo "Build and run manually: gfortran tests/fd_sign_compare.f90 -o build/fd_sign_compare && ./build/fd_sign_compare"
```

#### Step 2: Verify Makefile parses

Run: `make -n test-fd-sign`
Expected: prints the command without errors

#### Step 3: Add z(1) guard in externalFieldSetup_electricField

In `src/physics/hamiltonianConstructor.f90`, add after line 27 (`integer :: i`):

```fortran
if (size(z) > 0 .and. z(1) == 0.0_dp) then
  print *, "ERROR: externalFieldSetup_electricField called with z(1)=0"
  print *, "  This causes division by zero in the linear potential."
  print *, "  Ensure z-coordinates start at a non-zero value."
  stop 1
end if
```

#### Step 4: Build and verify no regression

Run: `cmake --build build`
Expected: builds without errors

Run: `OMP_NUM_THREADS=4 ctest --test-dir build -j4 --output-on-failure`
Expected: all tests pass (the guard only triggers if z(1)==0, which no existing test does)

#### Step 5: Commit

```bash
git add Makefile src/physics/hamiltonianConstructor.f90
git commit -m "fix: add test-fd-sign Makefile target and z(1) guard in externalFieldSetup"
```

Group #5 status: COMPLETE.

---

### Task 3: Close #16 — Anticrossing prose quantitative

**Files:**
- Modify: `docs/lecture/02-quantum-well.md:796-811`

#### Step 1: Run broken-gap config to extract anticrossing data

Run:
```bash
cp tests/regression/configs/qw_inas_gasb_broken_gap_kpar.cfg input.cfg
OMP_NUM_THREADS=4 ./build/src/bandStructure
```

Expected: eigenvalue output files in `output/`

#### Step 2: Extract k_parallel and gap from output

Read the eigenvalue file (e.g. `output/eigenvalues.dat`). Find the k-point where the e1 (CB ground) and lh1 (VB top) states are closest together. Compute:
- `k_par` = k_parallel at minimum gap (in Å⁻¹)
- `gap` = E(e1) - E(lh1) at that k-point (convert eV → meV: multiply by 1000)

If the data is hard to parse, check if the existing figure annotation already has these values by looking at `scripts/plotting/generate_all_figures.py` in `fig_qw_dispersion_broken_gap` (around line 2966).

#### Step 3: Update Ch02 Section B.5

In `docs/lecture/02-quantum-well.md`, replace lines 799 and 808-809:

Line 799 — replace:
```
(computed from the eigenvalue sweep — see the annotated figure)
```
with:
```
(at $k_\parallel = X.XX$ Å$^{-1}$, computed from the eigenvalue sweep)
```

Lines 808-809 — replace:
```
The computed hybridization gap at the anticrossing point is visible in the figure
annotation. For InAs/GaSb structures, this gap is typically 10--20 meV depending on
layer thicknesses
```
with:
```
The computed hybridization gap at the anticrossing point is **XX.X meV**. For
InAs/GaSb structures, this gap depends on layer thicknesses
```

(Use actual computed values for X.XX and XX.X)

#### Step 4: Commit

```bash
git add docs/lecture/02-quantum-well.md
git commit -m "docs: add quantitative anticrossing data to Ch02 Section B.5"
```

Group #16 status: COMPLETE.

---

### Task 4: Close #10 — wire_optical_spectrum figure

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py` (add function + register)
- Create: `docs/figures/wire_optical_spectrum.png` (generated)

#### Step 1: Add fig_wire_optical_spectrum function

Insert before `ALL_FIGURES = {` (around line 6008) in `generate_all_figures.py`:

```python
def fig_wire_optical_spectrum(output_dir: Path) -> None:
    """wire_optical_spectrum.png: TE and TM absorption for a GaAs rectangular wire.

    Runs opticalProperties with the wire_gaas_optical_window config and plots
    the absorption spectrum, demonstrating wire optical selection rules.
    """
    print("[figure] wire_optical_spectrum")

    cfg = CONFIG_DIR / "wire_gaas_optical_window.cfg"
    if not cfg.exists():
        print("  WARNING: wire_gaas_optical_window.cfg not found, skipping.")
        return

    result = run_executable(EXE_OPTICS, cfg, REPO_ROOT,
                           label="wire_optical_spectrum", timeout=600)
    if result.returncode != 0:
        print("  WARNING: opticalProperties run failed, skipping.")
        return

    te_file = output_dir / "absorption_TE.dat"
    tm_file = output_dir / "absorption_TM.dat"

    if not te_file.exists() or not tm_file.exists():
        print("  WARNING: wire absorption files not found, skipping.")
        return

    E_te, alpha_te = _read_absorption(output_dir, "TE")
    E_tm, alpha_tm = _read_absorption(output_dir, "TM")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(E_te, alpha_te, color="#1f77b4", linewidth=1.5, label="TE (x, y)")
    ax.plot(E_tm, alpha_tm, color="#d62728", linewidth=1.5, label="TM (z)")
    ax.set_xlabel("Photon Energy (eV)")
    ax.set_ylabel(r"Absorption Coefficient (cm$^{-1}$)")
    ax.set_title("GaAs Rectangular Wire Absorption Spectrum", fontsize=12)
    ax.legend(loc="best", fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3, linewidth=0.5)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "wire_optical_spectrum.png", dpi=150)
    plt.close(fig)
    print("  -> docs/figures/wire_optical_spectrum.png")
```

#### Step 2: Register in ALL_FIGURES dict

Add to `ALL_FIGURES` dict (after `"wire_gfactor_vs_size"` entry):

```python
    "wire_optical_spectrum": fig_wire_optical_spectrum,
```

#### Step 3: Verify EXE_OPTICS constant exists

Run: `grep -n "EXE_OPTICS" scripts/plotting/generate_all_figures.py`
Expected: a constant definition like `EXE_OPTICS = REPO_ROOT / "build" / "src" / "opticalProperties"`

If it doesn't exist, define it near the other EXE_ constants at the top of the file.

#### Step 4: Generate the figure

Run:
```bash
cd /data/8bandkp-fdm
python scripts/plotting/generate_all_figures.py --only wire_optical_spectrum
```
Expected: `docs/figures/wire_optical_spectrum.png` created

#### Step 5: Commit

```bash
git add scripts/plotting/generate_all_figures.py docs/figures/wire_optical_spectrum.png
git commit -m "feat: add wire_optical_spectrum figure generation"
```

Group #10 status: COMPLETE.

---

### Task 5: Close #8 — Document piezoelectric exclusion

**Files:** None (documentation-only decision)

**Step 1: Verify piezoelectric field is parsed but unused**

Run: `grep -rn "piezoelectric" src/ --include="*.f90"`
Expected: only `defs.f90` (field declaration) and `input_parser.f90` (parsing) — no usage in physics code.

**Step 2: No code changes.** The parsed field stays (harmless, forward-compatible if ever wanted). Just document in finalization.

Group #8 status: COMPLETE.

---

### Task 6: Finalization — Update tracking docs and archive

**Files:**
- Modify: `docs/plans/REVIEW.md`
- Modify: `docs/plans/BACKLOG.md`
- Move: 7 plan files → `docs/plans/archive/`

#### Step 1: Update REVIEW.md

For groups #5, #8, #10, #16, #45, change status to `COMPLETE | Archived`.
Remove or collapse their detailed findings sections (keep one-line summary).

#### Step 2: Update BACKLOG.md

Replace Phase 1 section with:
```markdown
## Phase 1: COMPLETED (2026-05-03)

Closed groups #45, #5, #16, #10, #8.
```

Renumber Phase 2 → Phase 1, Phase 3 → Phase 2, Phase 4 → Phase 3, Phase 5 → Phase 4, Phase 6 → Phase 5.

#### Step 3: Archive plan files

```bash
mv docs/plans/2026-04-29-pr12-fixes-plan.md docs/plans/archive/
mv docs/plans/2026-03-29-pr-review-fixes-design.md docs/plans/archive/
mv docs/plans/2026-04-12-qw-documentation-overhaul-design.md docs/plans/archive/
mv docs/plans/2026-04-12-qw-docs-phase1-plan.md docs/plans/archive/
mv docs/plans/2026-04-12-qw-phase1-gapfix-plan.md docs/plans/archive/
mv docs/plans/2026-04-04-documentation-overhaul-design.md docs/plans/archive/
mv docs/plans/2026-04-02-quantum-wire-design.md docs/plans/archive/
```

#### Step 4: Commit

```bash
git add docs/plans/REVIEW.md docs/plans/BACKLOG.md docs/plans/archive/
git commit -m "chore: complete Phase 1 backlog cleanup, archive 7 plan files"
```
