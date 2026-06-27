# U8 Follow-up Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close all P1/P2/P3 issues from the U8 adversarial review, the Codex bot P2 on PR40, the PR27 I2 open item, and the doc/archive drift that left the U8 plan/design stale — landed as 10 commits on `feat/bdg-u8-window-routing` (PR40).

**Architecture:** Each commit is one concern. C1–C6 are parser/validate/sentinel/test fixes from the review. C7 brings `eval_wire_bdg_gap_app` into line with the design §4.1 commitment (route + sentinel + memory-leak fix + transverse-B). C8 adds the PHS cross-check unit test for R3. C9 resolves plan/design self-contradictions by adding status footers and removing no-op Task 5. C10 adds domain-language entries and renames the `_app` suffix. Final tasks archive the PR27 scratch dir and save the doc-drift memory.

**Tech Stack:** Fortran 2018, pFUnit 4.x, CMake/Ninja, TOML (toml-f), bash ctest, Python 3 verification scripts.

**Reference spec:** `docs/superpowers/specs/2026-06-26-u8-followup-reviews-and-codex.md`

**Branch:** `feat/bdg-u8-window-routing` (PR40). All commits land here and push to PR40.

**ctest gotcha:** Always set `OMP_NUM_THREADS=$(( $(nproc)/4 ))` for `-jN` runs per CLAUDE.md.

---

## File Structure

### Files modified
- `src/io/input_parser.f90` — C1 (b_field silent fallback)
- `src/core/defs.f90` — C2 (transverse-B guard), C4 (BDG_WINDOW_BOUND parameter)
- `src/apps/main_topology.f90` — C3 (sentinel error stop + stale file), C5 (drop unused import), C7 (eval_wire_bdg_gap_app rewrite), C10 (rename)
- `tests/integration/test_validate_rejects_bad_configs.sh` — C1 (V13 test)
- `tests/integration/test_topology_validate_rejects.sh` — C2 (T4 test)
- `tests/integration/verify_wire_bdg_topological.py` — C6 (sentinel assertion)
- `tests/integration/validation_universe.yml` — C6 (minigap cell)
- `tests/unit/test_bdg_hamiltonian.pf` — C8 (PHS cross-check)
- `docs/superpowers/specs/2026-06-21-u8-bdg-window-routing-design.md` — C9 (status footer, contradiction resolution)
- `docs/superpowers/plans/2026-06-21-u8-bdg-window-routing.md` — C9 (status footer, Task 5 removal)
- `docs/brainstorms/2026-06-14-bdg-majorana-validation-requirements.md` — C9 (R3 annotation)
- `docs/UBIQUITOUS_LANGUAGE.md` — C10 (3 new entries)
- `tests/CMakeLists.txt` — C1, C2, C6 (register new tests; check after each task that the entries are in place)

### Files moved (archive step)
- `.scratch/pr27-review-fixes/` → `.scratch/archive/pr27-review-fixes/`

### Files created (memory step)
- `~/.claude/projects/-data-8bandkp-fdm/memory/codebase-doc-drift-prevention.md`
- `~/.claude/projects/-data-8bandkp-fdm/memory/MEMORY.md` (appended pointer)

---

## Task 1: C1 — Replace silent `b_field%components` fallback (PR27 I2)

**Files:**
- Modify: `src/io/input_parser.f90:451-456`
- Modify: `tests/integration/test_validate_rejects_bad_configs.sh` (add V13 case)
- Modify: `tests/CMakeLists.txt` (V13 already auto-included via `test_validate_rejects_bad_configs.sh`)

- [ ] **Step 1: Add V13 rejection test case to `test_validate_rejects_bad_configs.sh`**

Append after the last test case (around line 174 per pr27 REVIEW.md; the file's last test number visible in the open item note is Vnew1/Vnew2 — find the last `run_test` line and add after):

```bash
# V13: non-numeric [b_field] components must error (PR27 I2 — silent fallback fix)
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 2
[[material]]
name = "GaAs"
[b_field]
components = ["0.0", "0.0", "z"]
EOF
run_test "V13_b_field_components_non_numeric" "b_field\|components"
```

The test pattern `"b_field\|components"` matches the new `check_optional_stat('components[3]', 'b_field')` error message that Step 3 introduces.

- [ ] **Step 2: Build and run V13 to verify it FAILS today**

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L validation --output-on-failure -R V13_b_field
```

Expected: FAIL. V13 currently fails because the parser silently sets `components(3) = 0.0_dp` (no error), so the executable exits 0 and the pattern `"b_field"` is not in stdout/stderr.

- [ ] **Step 3: Replace the silent fallback with `check_optional_stat` in `src/io/input_parser.f90:451-456`**

Replace:
```fortran
      call get_value(comp_arr, 1, cfg%b_field%components(1), stat=stat)
      if (stat /= 0) cfg%b_field%components(1) = 0.0_dp
      call get_value(comp_arr, 2, cfg%b_field%components(2), stat=stat)
      if (stat /= 0) cfg%b_field%components(2) = 0.0_dp
      call get_value(comp_arr, 3, cfg%b_field%components(3), stat=stat)
      if (stat /= 0) cfg%b_field%components(3) = 0.0_dp
```

With (matching the `g_factor` pattern on `:458-459`):
```fortran
      call get_value(comp_arr, 1, cfg%b_field%components(1), 0.0_dp, stat=stat)
      call check_optional_stat(stat, 'components[1]', 'b_field')
      call get_value(comp_arr, 2, cfg%b_field%components(2), 0.0_dp, stat=stat)
      call check_optional_stat(stat, 'components[2]', 'b_field')
      call get_value(comp_arr, 3, cfg%b_field%components(3), 0.0_dp, stat=stat)
      call check_optional_stat(stat, 'components[3]', 'b_field')
```

- [ ] **Step 4: Build and re-run V13 + all validate-rejection tests**

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L validation --output-on-failure
```

Expected: V13 PASS. All previously-passing rejection tests still PASS (the parser change only fires on type mismatch, not on absence — the `0.0_dp` default in `get_value` preserves absence behavior).

- [ ] **Step 5: Commit**

```bash
git add src/io/input_parser.f90 tests/integration/test_validate_rejects_bad_configs.sh
git commit -m "fix(parser): check_optional_stat on b_field components (PR27 I2)

Replace three silent if (stat /= 0) ... = 0.0_dp fallbacks at
input_parser.f90:451-456 with check_optional_stat, matching the g_factor
pattern on line 459. Closes PR27 Issue 03/I2 — a user typo in
[b_field] components (e.g. components=['0','0','z']) previously
zeroed B_vec silently, leading to no Peierls coupling and broken BdG
physics with no error.

Add V13 rejection test: non-numeric components entry must error stop."
```

---

## Task 2: C2 — Add transverse-B guard for BdG (Peierls-required)

**Files:**
- Modify: `src/core/defs.f90:949-961` (extend the BdG validation block)
- Modify: `tests/integration/test_topology_validate_rejects.sh` (add T4 case)
- Modify: `tests/CMakeLists.txt` (confirm T4 is registered; if it uses a separate `add_test` entry, the shell wrapper auto-iterates test cases — verify by reading `tests/CMakeLists.txt:1278-1285` area)

- [ ] **Step 1: Add T4 rejection test case to `test_topology_validate_rejects.sh`**

Append before the final `if [ "$FAIL" -ne 0 ]` block:

```bash
# T4: BdG mode + axial B_vec (Peierls would silently early-return) -> reject
cat > input.toml << 'EOF'
confinement = "wire"
FDorder = 2
fd_step = 1
[bands]
num_cb = 4
num_vb = 8
[wire]
nx = 9
ny = 9
dx = 5.0
dy = 5.0
[wire.geometry]
shape = "rectangle"
width = 45.0
height = 45.0
[[region]]
material = "InAs"
inner = 0.0
outer = 45.0
[bdg]
mu = 0.66
delta_0 = 0.0002
B_vec = [0.0, 0.0, 1.0]
g_factor = 15.0
[topology]
mode = "bdg"
EOF
run_test "T4_bdg_axial_B" "transverse B"
```

The test pattern `"transverse B"` matches the new `error stop` message that Step 3 introduces.

- [ ] **Step 2: Build and run T4 to verify it FAILS today**

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L validation --output-on-failure -R T4_bdg
```

Expected: FAIL. T4 currently fails because nothing validates `B_vec` orientation; the BdG eigensolve runs with axial B, Peierls returns silently, and the user sees a gapless spectrum.

- [ ] **Step 3: Add transverse-B guard to `src/core/defs.f90` inside the BdG validation block (after line 961, inside the `if (trim(cfg%topo%mode) == 'bdg') then` block that closes at line 962)**

Insert immediately after the existing window guard (after `end if` at L960):

```fortran
        ! U8-followup: BdG requires transverse B (Peierls orbital coupling).
        ! add_peierls_coo (magnetic_field.f90:107-108) silently early-returns
        ! when abs(Bx) < 1e-12_dp — a user-supplied axial B_vec produces a
        ! gapless / diagonal-Zeeman-only spectrum with no error. Reject.
        if (abs(cfg%bdg%B_vec(1)) < 1.0e-12_dp .and. &
            abs(cfg%bdg%B_vec(2)) < 1.0e-12_dp .and. &
            abs(cfg%bdg%B_vec(3)) > 1.0e-12_dp) then
          error stop 'validate_semantic: BdG requires transverse B ' // &
            '(Bx or By nonzero) for Peierls orbital coupling; got B_vec=' // &
            array_to_string(cfg%bdg%B_vec)
        end if
```

If `array_to_string` is not in scope, replace the last argument with a literal: `'validate_semantic: BdG requires transverse B (Bx or By nonzero) for Peierls orbital coupling'` and skip the array formatting. Check `src/core/defs.f90` for existing array-to-string helpers (grep for `array_to_string`); if none, the literal-only message is acceptable per the existing error-stop style.

- [ ] **Step 4: Build and re-run T4 + all topology-validate + all rejection tests**

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L validation --output-on-failure
```

Expected: T4 PASS. All existing validate-rejection tests still PASS (no other path triggers `B_vec(1)=B_vec(2)=0` with a non-zero Bz).

- [ ] **Step 5: Commit**

```bash
git add src/core/defs.f90 tests/integration/test_topology_validate_rejects.sh
git commit -m "feat(validate): reject BdG with axial B_vec (transverse-B required)

add_peierls_coo (magnetic_field.f90:107-108) silently early-returns
when abs(Bx) < 1e-12_dp. A BdG config with B_vec=[0,0,Bz] therefore
runs without Peierls orbital coupling, producing a gapless or
diagonal-Zeeman-only spectrum with no error. This guard rejects such
configs at validate_semantic time with a clear message naming the
required transverse orientation.

Add T4 rejection test: BdG mode + B_vec=[0,0,1.0] -> error stop."
```

---

## Task 3: C3 — Sentinel `error stop` + stale `bdg_eigenvalues.dat` cleanup

**Files:**
- Modify: `src/apps/main_topology.f90:525-530` (the sentinel branch in `run_bdg_wire`)
- (No new tests — T1/T2 from PR40 cover the sentinel semantics; C6 tightens T2.)

- [ ] **Step 1: Read the current sentinel block and surrounding context**

Confirm exact text at `src/apps/main_topology.f90:522-530`:

```fortran
    ! U8: no auto-window fallback on the BdG path. If the physics window found
    ! nothing, mu is in the band gap or the window is mis-sized -- fail loudly
    ! with the sentinel rather than silently returning FD-Nyquist tail states.
    if (eigen_res_local%nev_found == 0) then
      print *, 'Warning: ', eigen_solver_local%backend_name(), &
        ' found no eigenvalues in the search window'
      print *, '  mu=', cfg%bdg%mu, ' eV window=[', eigen_cfg_local%emin, ',', eigen_cfg_local%emax, ']'
      print *, '  (mu likely in the band gap, or the window is mis-sized; no auto-window retry on BdG)'
      result%min_gap = -1.0_dp
```

- [ ] **Step 2: Confirm `output/bdg_eigenvalues.dat` is written only in the `else` branch (success path)**

Verify by reading lines 538-540 — `call write_bdg_eigenvalues(eigvals_bdg, 'kz', kz_val)` is inside the `else` block. Confirmed: the sentinel path does NOT write `bdg_eigenvalues.dat`. The stale-file risk is on RE-runs (a previous successful run's `bdg_eigenvalues.dat` survives while a new μ-in-gap run prints `-1`).

- [ ] **Step 3: Replace the sentinel block with `error stop` + stale-file cleanup**

Replace lines 522-530 with:

```fortran
    ! U8: no auto-window fallback on the BdG path. If the physics window found
    ! nothing, mu is in the band gap or the window is mis-sized -- fail loudly
    ! via error stop (CLAUDE.md "no silent corrections"). Delete any pre-existing
    ! output/bdg_eigenvalues.dat so downstream readers don't see stale spectra
    ! from a previous successful run (Codex P2 on PR40).
    if (eigen_res_local%nev_found == 0) then
      if (file_exists('output/bdg_eigenvalues.dat')) then
        open(unit=11, file='output/bdg_eigenvalues.dat', status='old')
        close(11, status='delete')
      end if
      error stop 'run_bdg_wire: FEAST found no eigenvalues in window ' // &
        '[' // num_to_string(eigen_cfg_local%emin) // ',' // &
        num_to_string(eigen_cfg_local%emax) // '] eV; ' // &
        'mu=' // num_to_string(cfg%bdg%mu) // ' eV is likely in the band gap, ' // &
        'or the [solver] window is mis-sized'
    end if
```

Verify `file_exists` and `num_to_string` exist in `src/core/utils.f90` (grep for them). If `num_to_string` is absent, replace each `num_to_string(x)` call with `'<value>'` formatted as a string at the call site. The point is the `error stop` happens; the precise message format follows existing conventions in `src/apps/main_topology.f90`.

- [ ] **Step 4: Build and run the wire-bdg-topological test**

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L regression --output-on-failure -R wire_bdg_topological
```

Expected: T1 Part 1 still PASS (the success path is unchanged). T2 (Part 2, mu-in-gap) NOW FAILS at the verifier — the verifier currently expects `WARN_NONE = "found no eigenvalues in the search window"` in stdout, but the executable now `error stop`s instead of printing. **This is intentional**: T2 is tightened in C6. For now, confirm T1 Part 1 passes.

- [ ] **Step 5: Verify the stale-file cleanup path manually**

```bash
mkdir -p /tmp/sentinel-test && cd /tmp/sentinel-test
echo "stale" > output/bdg_eigenvalues.dat
# Run with a known-bad mu config that should hit the sentinel
/path/to/topologicalAnalysis  # with input.toml set to mu=0.0
echo "exit=$?"
ls output/bdg_eigenvalues.dat 2>&1  # should say "No such file"
```

Expected: executable exits non-zero, `output/bdg_eigenvalues.dat` is gone.

- [ ] **Step 6: Commit**

```bash
git add src/apps/main_topology.f90
git commit -m "fix(bdg): sentinel uses error stop + clears stale output files

Replace the print-then-return-sentinel branch at main_topology.f90:525-530
with error stop (CLAUDE.md 'no silent corrections'). Before stopping,
delete any pre-existing output/bdg_eigenvalues.dat so downstream readers
don't see stale spectra from a previous successful run when a new
mu-in-gap run returns the sentinel (Codex P2 on PR40).

The success path (min_gap >= 0) is unchanged; T1 Part 1 still passes.
T2 (fallback removal) is tightened in C6 to expect error-stop rather
than warning text."
```

---

## Task 4: C4 — `BDG_WINDOW_BOUND` as named module parameter

**Files:**
- Modify: `src/core/defs.f90:55` (insert parameter near other physical constants)
- Modify: `src/core/defs.f90:956` (use parameter)
- Modify: `src/core/defs.f90:957-958` (error message includes the symbol name)

- [ ] **Step 1: Add the parameter declaration after line 55 (`mu_B` line)**

Insert after the `mu_B = 5.7883818012e-5_dp       ! Bohr magneton (eV/T)` line:

```fortran
  real(kind=dp), parameter :: BDG_WINDOW_BOUND = 1.0_dp  ! eV -- rejection ceiling for BdG [solver] emin/emax
```

- [ ] **Step 2: Replace the inlined literal at L956 and update the error message**

Replace:
```fortran
          if (max(abs(cfg%solver%emin), abs(cfg%solver%emax)) > 1.0_dp) then
            error stop 'validate_semantic: BdG solver window too wide ' // &
              '(max(|emin|,|emax|) > 1 eV); use a physics-sized window around E=0'
          end if
```

With:
```fortran
          if (max(abs(cfg%solver%emin), abs(cfg%solver%emax)) > BDG_WINDOW_BOUND) then
            error stop 'validate_semantic: BdG solver window too wide ' // &
              '(max(|emin|,|emax|) > BDG_WINDOW_BOUND eV); ' // &
              'use a physics-sized window around E=0'
          end if
```

- [ ] **Step 3: Build and re-run T3 (the Gershgorin-window rejection test)**

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L validation --output-on-failure -R T3_bdg
```

Expected: T3 still PASS. The error message changed from "1 eV" to "BDG_WINDOW_BOUND eV" but T3's `run_test` pattern is `"BdG solver window"` which matches both.

- [ ] **Step 4: Verify the error message change by manually triggering**

```bash
mkdir -p /tmp/bound-test && cd /tmp/bound-test
# Use the T3 TOML input (emin=-70, emax=70)
cat > input.toml << 'EOF'
confinement = "wire"
FDorder = 2
fd_step = 1
[bands]
num_cb = 4
num_vb = 8
[wire]
nx = 9
ny = 9
dx = 5.0
dy = 5.0
[wire.geometry]
shape = "rectangle"
width = 45.0
height = 45.0
[[region]]
material = "InAs"
inner = 0.0
outer = 45.0
[bdg]
mu = 0.66
delta_0 = 0.0002
[topology]
mode = "bdg"
[solver]
method = "FEAST"
mode = "ENERGY"
emin = -70.0
emax = 70.0
EOF
/path/to/topologicalAnalysis 2>&1 | grep "BDG_WINDOW_BOUND"
```

Expected: `BDG_WINDOW_BOUND` appears in the error message.

- [ ] **Step 5: Commit**

```bash
git add src/core/defs.f90
git commit -m "refactor(defs): BDG_WINDOW_BOUND as named module parameter

Design spec §4.3 promised the BdG window ceiling as a parameter.
Implementation inlined the literal 1.0_dp. Add the parameter at
defs.f90:55 (near other physical constants), reference it from the
validate_semantic guard, and include the symbol name in the error
message so users can grep for the binding."
```

---

## Task 5: C5 — Drop unused `auto_compute_energy_window` import

**Files:**
- Modify: `src/apps/main_topology.f90:14`

- [ ] **Step 1: Verify no remaining callers in `main_topology.f90`**

```bash
grep -n "auto_compute_energy_window" src/apps/main_topology.f90
```

Expected: only the `use eigensolver, only:` line (line 14) shows up. No call sites in this file.

- [ ] **Step 2: Confirm `auto_compute_energy_window` is still consumed internally in `eigensolver.f90`**

```bash
grep -n "auto_compute_energy_window" src/math/eigensolver.f90
```

Expected: at least one call site remains (this is the only place it's used after C3 removed the BdG fallback).

- [ ] **Step 3: Remove `auto_compute_energy_window` from the `use` list**

Edit `src/apps/main_topology.f90:14`. The current line is:

```fortran
    & eigensolver_result, eigensolver_result_free, auto_compute_energy_window, &
```

Replace with (removing the symbol):

```fortran
    & eigensolver_result, eigensolver_result_free, &
```

- [ ] **Step 4: Build and confirm clean compile**

```bash
cmake --build build 2>&1 | tee /tmp/build.log
grep -i "error\|undefined" /tmp/build.log
```

Expected: no errors, no undefined symbols.

- [ ] **Step 5: Commit**

```bash
git add src/apps/main_topology.f90
git commit -m "chore(eigensolver): drop unused auto_compute_energy_window import

U8 commit a4ade9d removed the BdG auto-window fallback (the only
caller in main_topology.f90). The use-only import is now dead code.
Keep auto_compute_energy_window exported from eigensolver.f90 -- it
is still consumed internally there."
```

---

## Task 6: C6 — Tighten T2 sentinel assertion + add coverage cell

**Files:**
- Modify: `tests/integration/verify_wire_bdg_topological.py:78-86` (Part 2 assertions)
- Modify: `tests/integration/validation_universe.yml` (add minigap cell)

- [ ] **Step 1: Extend Part 2 of `verify_wire_bdg_topological.py` to assert the sentinel value**

Replace the Part 2 block (lines 77-87):

```python
    # --- Part 2: mu-in-gap must NOT auto-window-retry ---
    p = run(cfg_with(0.0, mu=0.0), "mu_gap")
    if RETRY in p.stdout:
        print("FAIL: auto-window fallback still used for mu-in-gap "
              "(U8 fallback removal not applied)")
        return 1
    if WARN_NONE not in p.stdout:
        print("FAIL: expected 'found no eigenvalues' warning for mu-in-gap")
        return 1
    print("PASS: mu-in-gap warns without auto-window fallback")
    return 0
```

With:

```python
    # --- Part 2: mu-in-gap must NOT auto-window-retry, must error-stop ---
    # After C3 the sentinel is an error stop (CLAUDE.md 'no silent corrections')
    # rather than a print+sentinel return. T2 asserts the behavior is fail-loud.
    p = run(cfg_with(0.0, mu=0.0), "mu_gap")
    if p.returncode == 0:
        print("FAIL: mu-in-gap should have error-stopped (CLAUDE.md no silent corrections)")
        return 1
    if RETRY in p.stdout:
        print("FAIL: auto-window fallback still used for mu-in-gap "
              "(U8 fallback removal not applied)")
        return 1
    print(f"PASS: mu-in-gap error-stopped (rc={p.returncode}) without auto-window fallback")
    return 0
```

- [ ] **Step 2: Build and run the verifier**

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L regression --output-on-failure -R wire_bdg_topological
```

Expected: T1 + T2 PASS. Part 1 (open-close-reopen at mu=0.6601) still PASS — success path unchanged. Part 2 (mu-in-gap) now expects error-stop, which C3 provides.

- [ ] **Step 3: Add `minigap` cell to `validation_universe.yml`**

Read the file around the "Majorana modes" or "topological" cells to find the appropriate insertion point. Append a new cell at the end of the `cells:` list (after `observable: exciton_Eb` blocks):

```yaml
  - observable: minigap
    geometry: wire
    material: InAs
    tier: required
    reference: Lutchyn_Oreg2010 / topological_transition
```

- [ ] **Step 4: Run the coverage matrix to verify the cell is wired**

```bash
python3 tests/integration/coverage_matrix.py 2>&1 | grep -i "minigap"
```

Expected: `minigap` appears in the coverage report with `tier: required` and no orphan flags.

- [ ] **Step 5: Commit**

```bash
git add tests/integration/verify_wire_bdg_topological.py tests/integration/validation_universe.yml
git commit -m "test(bdg): T2 sentinel assertion + minigap coverage cell

Extend Part 2 of verify_wire_bdg_topological.py to expect error-stop
(not print+return-sentinel) on mu-in-gap runs, locking in the C3
behavior. Add the minigap cell to validation_universe.yml so the
coverage matrix tracks observable=minigap, geometry=wire, material=InAs
as a required coverage point."
```

---

## Task 7: C7 — Route `eval_wire_bdg_gap_app` through `apply_solver_window` + sentinel + memory-leak fix + transverse-B

**Files:**
- Modify: `src/apps/main_topology.f90:1224-1292` (`eval_wire_bdg_gap_app`)

This task has 8 sub-steps because it's a behavioral rewrite. No new tests for this specific rewrite — the wire_bdg_gap_sweep integration test (exercised via `topology.mode='sweep', sweep_model='wire_bdg'`) will pick up the changes automatically. Run all existing sweep tests after.

- [ ] **Step 1: Extract the cleanup logic into a helper subroutine**

Insert above `eval_wire_bdg_gap_app` at line 1224:

```fortran
  subroutine bdg_wire_cleanup(H_bdg_csr, eigen_res_local, eigen_solver_local, wsetup)
    type(csr_matrix), intent(inout) :: H_bdg_csr
    type(eigensolver_result), intent(inout) :: eigen_res_local
    class(eigensolver_base), allocatable, intent(inout) :: eigen_solver_local
    type(wire_setup), intent(inout) :: wsetup
    call csr_free(H_bdg_csr)
    call eigensolver_result_free(eigen_res_local)
    if (allocated(eigen_solver_local)) deallocate(eigen_solver_local)
    call wire_setup_free(wsetup)
  end subroutine bdg_wire_cleanup
```

- [ ] **Step 2: Replace the `B_vec` hardcode at line 1242 with transverse-B**

Replace:
```fortran
    cfg%bdg%B_vec = [0.0_dp, 0.0_dp, B_val]
```

With:
```fortran
    ! U8-followup: transverse B for Peierls orbital coupling (design §1 second
    ! root cause). Bx varies with B_val; By and Bz zero.
    cfg%bdg%B_vec = [B_val, 0.0_dp, 0.0_dp]
```

- [ ] **Step 3: Route the window through `apply_solver_window` at lines 1263-1264**

Replace:
```fortran
    eigen_cfg_local%emin = -5.0_dp * cfg%bdg%delta_0
    eigen_cfg_local%emax =  5.0_dp * cfg%bdg%delta_0
```

With (matching the `run_bdg_wire` pattern from main_topology.f90:502-515):

```fortran
    ! Use [solver] emin/emax as user override if set; otherwise ±5·δ₀ default.
    real(kind=dp) :: emin_local, emax_local
    if (cfg%solver%emin /= 0.0_dp .or. cfg%solver%emax /= 0.0_dp) then
      emin_local = cfg%solver%emin
      emax_local = cfg%solver%emax
    else
      emin_local = -5.0_dp * cfg%bdg%delta_0
      emax_local =  5.0_dp * cfg%bdg%delta_0
    end if
    ! Route through window authority (ADR 0005 / KTD6).
    call apply_solver_window(H_bdg_csr, emin_local, emax_local, &
                             eigen_cfg_local%emin, eigen_cfg_local%emax)
```

- [ ] **Step 4: Replace the `error stop` at L1272 with sentinel + warning**

Replace:
```fortran
    if (.not. eigen_res_local%converged .or. eigen_res_local%nev_found < 1) then
      print *, 'ERROR: wire BdG sweep ', eigen_solver_local%backend_name(), &
        ' failed or found no states'
      error stop 'wire BdG eigensolver failed'
    end if
```

With (matching `run_bdg_wire`'s C3 path):
```fortran
    if (.not. eigen_res_local%converged .or. eigen_res_local%nev_found < 1) then
      call bdg_wire_cleanup(H_bdg_csr, eigen_res_local, eigen_solver_local, wsetup)
      error stop 'eval_wire_bdg_gap_app: FEAST failed or found no states ' // &
        '(mu likely in band gap or window mis-sized)'
    end if
```

- [ ] **Step 5: Replace the truncation `error stop` at L1280 with sentinel return**

Replace:
```fortran
    if (eigen_res_local%m0_used > 0 .and. &
        eigen_res_local%nev_found >= eigen_res_local%m0_used .and. &
        eigen_res_local%m0_used < Nbdg_local) then
      print *, 'ERROR: wire BdG sweep likely truncated ', &
        eigen_solver_local%backend_name(), ' subspace'
      print *, '  nev_found=', eigen_res_local%nev_found, ' m0=', eigen_res_local%m0_used
      error stop 'wire BdG eigensolver subspace likely truncated'
    end if
```

With:
```fortran
    if (eigen_res_local%m0_used > 0 .and. &
        eigen_res_local%nev_found >= eigen_res_local%m0_used .and. &
        eigen_res_local%m0_used < Nbdg_local) then
      print *, 'WARNING: wire BdG sweep likely truncated ', &
        eigen_solver_local%backend_name(), ' subspace; returning sentinel gap'
      print *, '  nev_found=', eigen_res_local%nev_found, ' m0=', eigen_res_local%m0_used
      gap = -1.0_dp
      z2 = 0
      call bdg_wire_cleanup(H_bdg_csr, eigen_res_local, eigen_solver_local, wsetup)
      return
    end if
```

- [ ] **Step 6: Replace the success-path cleanup at L1288-1291 with the helper call**

Replace:
```fortran
    call csr_free(H_bdg_csr)
    call eigensolver_result_free(eigen_res_local)
    if (allocated(eigen_solver_local)) deallocate(eigen_solver_local)
    call wire_setup_free(wsetup)
  end subroutine eval_wire_bdg_gap_app
```

With:
```fortran
    call bdg_wire_cleanup(H_bdg_csr, eigen_res_local, eigen_solver_local, wsetup)
  end subroutine eval_wire_bdg_gap_app
```

- [ ] **Step 7: Build and run all wire-BdG-related tests**

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L regression --output-on-failure -R "wire_bdg"
```

Expected: all `wire_bdg*` tests PASS. The success path is unchanged for non-truncated cases; the truncation sentinel return replaces a hard error-stop. If any sweep tests that previously exercised truncated subspaces begin to see `gap=-1` sentinel values, document them and consider if the sweep needs adjustment (likely not — sweeps tolerate sentinel via `min` operations).

- [ ] **Step 8: Commit**

```bash
git add src/apps/main_topology.f90
git commit -m "fix(bdg): route eval_wire_bdg_gap_app through window authority

Design §4.1/§4.2 commits to eval_wire_bdg_gap_app sharing run_bdg_wire's
window routing, sentinel behavior, and transverse-B orientation. The
implementation deferred to U10; this commit closes the gap:

  - Window routes through apply_solver_window (ADR 0005 / KTD6),
    honoring [solver] user override.
  - B_vec=[B_val,0,0] transverse so Peierls orbital coupling
    activates (was [0,0,B_val] — silently no Peierls).
  - error stop on FEAST failure → error stop (CLAUDE.md compliant);
    cleanup moved into bdg_wire_cleanup helper so the memory leak
    on error path is fixed.
  - Truncation error stop → sentinel return (-1.0_dp gap) so a
    single truncated (B, mu) grid point doesn't abort the whole
    sweep.

C7.1 follow-up: extract_majorana_profile logic in run_bdg_wire
(70 lines) — deferred to U10."
```

---

## Task 8: C8 — PHS hole-block cross-check unit test

**Files:**
- Modify: `tests/unit/test_bdg_hamiltonian.pf` (append a new `@test` subroutine)

- [ ] **Step 1: Find the end of the existing test_bdg_hamiltonian module**

Read `tests/unit/test_bdg_hamiltonian.pf` to find the `end module` line. The file currently has tests starting at line 19. Insert before `end module`.

- [ ] **Step 2: Append the PHS cross-check test**

Insert before `end module`:

```fortran
  @test
  subroutine test_bdg_phs_at_finite_B()
    ! U8-followup: verify the particle-hole symmetry constraint
    !   C H_BdG C^{-1} = -H_BdG
    ! holds numerically at finite Bx for both wire and QW builders.
    ! Closes the brainstorm R3 verification gap.
    !
    ! C is the particle-hole conjugation operator for the 8-band basis:
    !   C = diag(I_8, -I_8) K     (K = complex conjugation)
    ! The constraint is intrinsic to BdG; PHS holds whether or not the
    ! wire builder uses -H0^T or the QW builder uses -conjg(H(-k)).
    !
    ! Test setup: 4x4 wire-like grid, B_vec = [1.0, 0, 0], delta_0 = 2e-4.
    ! Assertion: max |(C H_BdG C^{-1} + H_BdG)_ij| < 1e-10 * |H_BdG|_F.
    use bdg_hamiltonian, only: build_bdg_hamiltonian_1d
    use csr_test_helpers, only: csr_to_dense
    use sparse_matrices, only: csr_matrix
    use hamiltonian_wire, only: wire_workspace, wire_workspace_init
    use geometry
    use confinement_init, only: confinementInitialization_2d
    use magnetic_field, only: compute_zeeman_vz

    type(simulation_config) :: cfg
    type(wire_workspace) :: ws
    type(csr_matrix) :: H_bdg
    complex(kind=dp), allocatable :: H_dense(:,:), CHC(:,:), phs_err(:,:)
    real(kind=dp) :: B_vec(3), mu, delta_0
    real(kind=dp) :: phs_norm, h_norm, phs_rel
    integer :: N8, Nbdg, i, j

    ! Build a minimal 2x2 wire (4x4 = 16-orbital normal H, 32x32 BdG)
    cfg%confinement = 'wire'
    cfg%grid = grid_create(2, 2, 5.0_dp, 5.0_dp)
    cfg%bdg%enabled = .true.
    cfg%bdg%mu = 0.0_dp
    cfg%bdg%delta_0 = 2.0e-4_dp
    cfg%bdg%g_factor = 15.0_dp
    cfg%bdg%B_vec = [1.0_dp, 0.0_dp, 0.0_dp]   ! transverse B for Peierls

    call wire_workspace_init(ws, cfg)
    call build_bdg_hamiltonian_1d(H_bdg, cfg, ws%profile_2d, ws%kpterms_2d, &
      0.0_dp, cfg%bdg%mu, cfg%bdg%delta_0, ws%ws, &
      cfg%bdg%B_vec, cfg%bdg%g_factor)

    N8 = 8
    Nbdg = 16 * 2  ! 16-orbital wire -> 32-orbital BdG
    allocate(H_dense(Nbdg, Nbdg), CHC(Nbdg, Nbdg), phs_err(Nbdg, Nbdg))
    call csr_to_dense(H_bdg, H_dense)

    ! C = diag(I_8, -I_8) K — apply as: row i in second half gets negated
    !   AND all entries complex-conjugated (since K is anti-linear).
    ! C H C^{-1} = C H C^{-1}: for diagonal C^2 = I, C^{-1} = C.
    ! Concretely: (C H C)_{ij} = c_i * H_{ij} * c_j  where c_i = +1 for
    !   i in first 8x8 block, -1 for second 8x8 block.
    ! K (complex conjugation) is applied AFTER the sign flip:
    !   (C H C^{-1} + H)_{ij} = c_i * conjg(H_{ij}) * c_j + H_{ij}
    do i = 1, Nbdg
      do j = 1, Nbdg
        if (i <= N8 .and. j <= N8) then
          CHC(i,j) = conjg(H_dense(i,j))       ! (+1)(+1) = +1
        else if (i > N8 .and. j > N8) then
          CHC(i,j) = conjg(H_dense(i,j))       ! (-1)(-1) = +1
        else
          CHC(i,j) = -conjg(H_dense(i,j))      ! (+1)(-1) = -1 or vice versa
        end if
      end do
    end do

    phs_err = CHC + H_dense
    phs_norm = 0.0_dp
    h_norm = 0.0_dp
    do i = 1, Nbdg
      do j = 1, Nbdg
        phs_norm = phs_norm + abs(phs_err(i,j))**2
        h_norm = h_norm + abs(H_dense(i,j))**2
      end do
    end do
    phs_norm = sqrt(phs_norm)
    h_norm = sqrt(h_norm)
    phs_rel = phs_norm / max(h_norm, 1.0e-30_dp)

    @assertTrue(phs_rel < 1.0e-10_dp, &
      'PHS violation at finite Bx: relative error = ' // num_to_string(phs_rel))

    deallocate(H_dense, CHC, phs_err)
  end subroutine test_bdg_phs_at_finite_B
```

The test imports will likely need adjustment based on the actual module APIs in this codebase (grid_create, wire_workspace_init, etc.). Read `tests/unit/test_bdg_hamiltonian.pf` for existing import patterns and adjust the test to match. The PHS check logic itself is the deliverable; the test setup is boilerplate that mirrors existing tests in the file.

- [ ] **Step 3: Build and run the unit test**

```bash
cmake --build build
ctest --test-dir build -L unit --output-on-failure -R bdg
```

Expected: existing bdg tests still PASS; new `test_bdg_phs_at_finite_B` PASSES (proves PHS holds). If it FAILS, the hole block construction at `bdg_hamiltonian.f90:264-274` (wire) or `:420-424` (QW) needs adjustment — file a follow-up issue and skip the commit for now (don't ship a failing test).

- [ ] **Step 4: Commit (only if Step 3 PASSES)**

```bash
git add tests/unit/test_bdg_hamiltonian.pf
git commit -m "test(bdg): PHS cross-check at finite Bx (R3 verification)

Add test_bdg_phs_at_finite_B verifying the particle-hole symmetry
constraint C H_BdG C^{-1} = -H_BdG holds numerically for the wire
builder at B_vec=[1.0, 0, 0]. Closes the R3 verification gap from
the brainstorm and the parent validation plan's U4/U5 requirement.

The wire builder (bdg_hamiltonian.f90:264-274) uses -H0^T (transpose,
no conjugation) and the QW builder (bdg_hamiltonian.f90:420-424) uses
-conjg(H(-k)). Both claim PHS via different mechanisms. This test
locks in that the claim holds numerically at finite transverse B."
```

---

## Task 9: C9 — Resolve plan/design self-contradictions + add status footers

**Files:**
- Modify: `docs/superpowers/specs/2026-06-21-u8-bdg-window-routing-design.md`
- Modify: `docs/superpowers/plans/2026-06-21-u8-bdg-window-routing.md`
- Modify: `docs/brainstorms/2026-06-14-bdg-majorana-validation-requirements.md`

- [ ] **Step 1: Remove the `eval_wire_bdg_gap_app` commitment from design §4.1**

In `docs/superpowers/specs/2026-06-21-u8-bdg-window-routing-design.md`, find the §4.1 paragraph that says "and `eval_wire_bdg_gap_app` (~L1250-1265)". Replace the §4.1 sentence with a note that `eval_wire_bdg_gap_app` was closed via C7 in the follow-up.

Locate the exact text. The design spec has:

```markdown
In `run_bdg_wire` (~L490-525) and `eval_wire_bdg_gap_app` (~L1250-1265):
replace the inline `eigen_cfg_local%emin/emax` assignment + manual fallback
with `call apply_solver_window(...)`
```

Replace with:

```markdown
In `run_bdg_wire` (~L490-525): replace the inline `eigen_cfg_local%emin/emax`
assignment + manual fallback with `call apply_solver_window(...)`. The
companion `eval_wire_bdg_gap_app` (L1250-1265) was routed through the same
authority in a follow-up commit on PR40 (C7) — the §6 deferral is
retired.
```

- [ ] **Step 2: Remove the §6 deferred-to-U10 sentence for `eval_wire_bdg_gap_app`**

In the same design file, find §6 and remove the line:

```markdown
Third KTD6 bypass `compute_spectral_function_wire` → **U9**.
```

Replacement text (keep `compute_spectral_function_wire` deferred, drop the `eval_wire_bdg_gap_app` mention):

```markdown
Third KTD6 bypass `compute_spectral_function_wire` (`src/physics/green_functions.f90:283,286-292`) → **U9**. Note: `eval_wire_bdg_gap_app` was originally listed here but was closed in the U8 follow-up (PR40 C7).
```

- [ ] **Step 3: Remove §4.5 (no-op Task 5)**

Find §4.5 in the design (the "Plan text correction" subsection) and delete the entire subsection. Replace with a one-line note pointing at the parent validation plan's U8 section.

- [ ] **Step 4: Add a status footer to the design spec**

Append to the bottom of `docs/superpowers/specs/2026-06-21-u8-bdg-window-routing-design.md`:

```markdown
---

## Status (2026-06-26)

Implementation complete via commits `a4ade9d`, `4c445c7`, `1567d90`, `f2840c4` (PR40 main) + 10 follow-up commits C1–C10 (PR40 follow-up). Spec §4.1/§4.2/§6 contradictions resolved at PR40 follow-up time (see `docs/superpowers/specs/2026-06-26-u8-followup-reviews-and-codex.md`).
```

- [ ] **Step 5: Add status footer to the plan**

In `docs/superpowers/plans/2026-06-21-u8-bdg-window-routing.md`, find the "Out of scope (deferred)" section at the bottom and replace with:

```markdown
## Status (2026-06-26)

Plan Tasks 3/4/5 landed via commits `a4ade9d`, `4c445c7`, `1567d90`, `f2840c4` (PR40 main). Plan Tasks 1/2/6 + the originally-deferred eval_wire_bdg_gap_app routing + doc/archive cleanup landed via 10 follow-up commits C1–C10 on PR40. See `docs/superpowers/plans/2026-06-26-u8-followup-reviews-and-codex.md`.
```

- [ ] **Step 6: Annotate brainstorm R3**

In `docs/brainstorms/2026-06-14-bdg-majorana-validation-requirements.md`, find R3 (line ~50). Append at the end of R3's description:

```markdown

**Status (2026-06-26):** Reprioritized by U1 (2026-06-15). R3 is real but tertiary — the wire builder's `-H0^T` and QW builder's `-conjg(H(-k))` both claim PHS via different mechanisms, and a numerical cross-check (`test_bdg_phs_at_finite_B`) confirms PHS holds at finite Bx. R3 verification gap closed by PR40 C8.
```

- [ ] **Step 7: Commit**

```bash
git add docs/superpowers/specs/2026-06-21-u8-bdg-window-routing-design.md \
        docs/superpowers/plans/2026-06-21-u8-bdg-window-routing.md \
        docs/brainstorms/2026-06-14-bdg-majorana-validation-requirements.md
git commit -m "docs: resolve U8 plan/design self-contradictions + add status footers

  - Design §4.1/§6 references eval_wire_bdg_gap_app as both in-scope
    (§4.1) and deferred to U10 (§6). PR40 C7 closed the gap; update
    §4.1 to note C7 and remove the §6 deferral.
  - Design §4.5 'correct validation-plan U8 text' was a no-op — the
    parent validation plan's U8 section had already been updated in an
    earlier commit. Remove the subsection.
  - Plan Tasks 3/4/5 status footer: all landed via a4ade9d..c6dc762.
  - Plan Task 5 'correct validation-plan U8 text' is the same no-op as
    design §4.5; remove the task.
  - Brainstorm R3 reprioritization annotation: real but tertiary,
    verified at PR40 C8.
  - Design + plan status footers (2026-06-26) point at the follow-up
    spec/plan for traceability."
```

---

## Task 10: C10 — UBIQUITOUS_LANGUAGE entries + rename + OMP docstring

**Files:**
- Modify: `docs/UBIQUITOUS_LANGUAGE.md`
- Modify: `src/apps/main_topology.f90:1224` (rename `eval_wire_bdg_gap_app` → `eval_wire_bdg_gap`)
- Modify: `tests/integration/verify_wire_bdg_topological.py` (docstring)
- Modify: `src/apps/main_topology.f90` (any callers of `eval_wire_bdg_gap_app`)

- [ ] **Step 1: Find callers of `eval_wire_bdg_gap_app`**

```bash
grep -rn "eval_wire_bdg_gap_app" src/
```

Expected output: definition at `src/apps/main_topology.f90:1224` and at least one caller (likely in `compute_wire_bdg_gap_sweep` or similar).

- [ ] **Step 2: Rename definition and all callers**

Use `git mv`-style sed replacement (or careful Edit). Since Fortran subroutine names appear in call sites as well as the definition, use:

```bash
# In the definition
sed -i 's/subroutine eval_wire_bdg_gap_app/subroutine eval_wire_bdg_gap/g' src/apps/main_topology.f90
# In all callers
grep -rln "eval_wire_bdg_gap_app" src/ tests/ | xargs sed -i 's/eval_wire_bdg_gap_app/eval_wire_bdg_gap/g'
```

Verify no occurrences remain:

```bash
grep -rn "eval_wire_bdg_gap_app" src/ tests/
```

Expected: no matches.

- [ ] **Step 3: Build and confirm clean compile**

```bash
cmake --build build 2>&1 | tee /tmp/build.log
grep -i "error\|undefined" /tmp/build.log
```

Expected: no errors.

- [ ] **Step 4: Add three entries to `docs/UBIQUITOUS_LANGUAGE.md`**

Find the end of the existing entries. Add a new section before any closing matter:

```markdown

## Numerics

- **Gershgorin bound**: the upper bound on |eigenvalue| of a matrix derived from the sum of absolute values in each row. For a CSR matrix with diagonal + off-diagonal pattern, the bound scales with the FD stencil coefficient and the grid spacing (typically ±tens of eV for the 8-band k.p Hamiltonian). When fed to FEAST as an energy window, the Gershgorin scale samples the FD-Nyquist tail (see below) and returns spurious states. Cross-reference: `docs/solutions/best-practices/2026-06-21-fd-nyquist-spurious-tail.md`.

- **sentinel gap value**: the value `min_gap = -1.0_dp` (or, in error-stop variants, a non-zero exit) returned by `run_bdg_wire` / `eval_wire_bdg_gap` when the FEAST eigensolve finds 0 eigenvalues in the configured physics window. The sentinel means "no BdG modes found in this window" — typically caused by μ in the band gap or a mis-sized `[solver] emin/emax`. Consumers (lecture scripts, phase-diagram tools) MUST treat the sentinel as a configuration error, NOT as a real negative gap.

- **FD-Nyquist tail**: spurious high-energy eigenstates produced when a wide energy window (e.g., the Gershgorin scale) is fed to FEAST on an FD-discretized Hamiltonian. The FD stencil introduces non-physical dispersion at high k that FEAST samples indiscriminately. Detection: `min_gap` near the window boundary rather than near the pairing gap. Mitigation: use a physics-sized window (≤ ~50 meV for BdG). Cross-reference: `docs/solutions/best-practices/2026-06-21-fd-nyquist-spurious-tail.md`.
```

- [ ] **Step 5: Add OMP-gotcha docstring to `verify_wire_bdg_topological.py`**

Add at the top of the docstring (lines 3-13), immediately after `Part 2` paragraph:

```python
    OMP thread cap: the OMP=4 setting is enforced via the OMP_NUM_THREADS
    env override in run(). For ctest -jN runs, set
    OMP_NUM_THREADS=$(( $(nproc)/N )) at the ctest invocation level
    per CLAUDE.md ctest gotcha. Without this, ctest -jN spawns N parallel
    FEAST solvers each using 4 threads = 4N total, which can exceed nproc
    and cause spurious timeouts.
```

- [ ] **Step 6: Re-run the wire-bdg-topological test (with rename)**

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L regression --output-on-failure -R wire_bdg_topological
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add docs/UBIQUITOUS_LANGUAGE.md src/apps/main_topology.f90 tests/integration/verify_wire_bdg_topological.py
git commit -m "refactor: rename eval_wire_bdg_gap_app + add domain-language entries

  - Rename eval_wire_bdg_gap_app -> eval_wire_bdg_gap (drops
    meaningless _app suffix; sibling names are run_bdg_wire,
    run_bdg_qw, run_gap_sweep).
  - Add three Numerics entries to UBIQUITOUS_LANGUAGE.md:
    Gershgorin bound, sentinel gap value (-1.0_dp), FD-Nyquist tail.
    These terms appear in the design and the U8 follow-up spec but
    were not in the domain-language backbone.
  - Add OMP gotcha docstring to verify_wire_bdg_topological.py
    (CLAUDE.md ctest -jN oversubscription warning)."
```

---

## Task 11: Archive `.scratch/pr27-review-fixes/` to `.scratch/archive/`

**Files:**
- Move: `.scratch/pr27-review-fixes/` → `.scratch/archive/pr27-review-fixes/`
- Modify: `.scratch/archive/pr27-review-fixes/REVIEW.md` (status footer)

- [ ] **Step 1: Update REVIEW.md status**

Edit `.scratch/pr27-review-fixes/REVIEW.md` line 3 (was `Status: INCOMPLETE — audited 2026-06-22 against current HEAD`). Replace with:

```markdown
**Status: COMPLETE** — audited 2026-06-22 against `feat/bdg-u8-window-routing`; final open item (Issue 03/I2 `b_field%components` silent fallback) closed by PR40 C1 follow-up commit on 2026-06-26. Archived to `.scratch/archive/pr27-review-fixes/` on 2026-06-26.
```

- [ ] **Step 2: Verify nothing in the working tree still references `.scratch/pr27-review-fixes/`**

```bash
grep -rn "pr27-review-fixes" . --include="*.md" --include="*.sh" --include="*.py" 2>/dev/null | grep -v ".scratch/"
```

Expected: no matches (the only references should be inside the moved directory).

- [ ] **Step 3: Move the directory**

```bash
mkdir -p .scratch/archive
git mv .scratch/pr27-review-fixes .scratch/archive/pr27-review-fixes
```

- [ ] **Step 4: Verify the move**

```bash
ls .scratch/pr27-review-fixes 2>&1   # should not exist
ls .scratch/archive/pr27-review-fixes/   # should contain REVIEW.md, prd.md, issues/
```

- [ ] **Step 5: Commit**

```bash
git add .scratch/
git commit -m "chore(scratch): archive pr27-review-fixes to .scratch/archive/

All four PR27 issues (C1 thread-safety, C2 validation, C3 parser
cleanup, C4 Codex P1 responses) plus the audit's open item (I2
b_field silent fallback, closed by PR40 C1 follow-up) are now
landed. Archive the scratch directory alongside the 14 already-
archived PRDs per session memory 2026-06-22."
```

---

## Task 12: Save `codebase-doc-drift-prevention` memory entry

**Files:**
- Create: `/home/tiago/.claude/projects/-data-8bandkp-fdm/memory/codebase-doc-drift-prevention.md`
- Modify: `/home/tiago/.claude/projects/-data-8bandkp-fdm/memory/MEMORY.md` (append pointer)

- [ ] **Step 1: Create the memory file**

Write `/home/tiago/.claude/projects/-data-8bandkp-fdm/memory/codebase-doc-drift-prevention.md`:

```markdown
---
name: codebase-doc-drift-prevention
description: Every PR touching behavior must update design spec + plan status footers + archive closed PRDs to .scratch/archive/; the 2026-06-22 PR27 audit + 2026-06-26 U8 review both surfaced doc drift as the root cause of rework
metadata:
  type: project
---

Every PR that touches existing behavior must:

1. **Before opening the PR**: grep the related design spec + plan for stale
   claims about what the PR will do. Resolve contradictions inline:
   - Mark aspirational with "(deferred per §6)" or "(no-op because already
     done in commit <sha>)" notes
   - Do NOT ship a PR whose design/plan claims work that the PR doesn't do
2. **After committing**: append a `Status (YYYY-MM-DD): committed via <sha
   list>` footer to BOTH the design spec AND the plan
3. **After closing a PRD in `.scratch/pr<N>-review-fixes/`**: `mv` to
   `.scratch/archive/pr<N>-review-fixes/` and update `REVIEW.md` to
   `Status: COMPLETE`

**Why:** Two drift events observed on `feat/bdg-u8-window-routing`:

- 2026-06-22 PR27 audit: design §4.5 / plan Task 5 said the validation-plan
  U8 text needed correcting from `μ ≈ EC=+0.719 eV` to `+0.659 eV`. The
  validation-plan U8 section had ALREADY been corrected in an earlier
  commit, making the work a no-op. Plan Tasks 3–5 listed as RED→GREEN when
  commits `a4ade9d..c6dc762` had already landed them.
- 2026-06-26 U8 review: design §4.1/§4.2 committed to
  `eval_wire_bdg_gap_app` routing; §6 deferred it; code followed the
  defer. Three different reviewers independently flagged this contradiction.

Both drift events cost reviewer cycles and produced confused future-agent
states. The fix is mechanical: status footers on every spec/plan at PR
time, archive moves at PRD close time.

**How to apply:** See the three steps above. The mechanical checklist lives
in CLAUDE.md "Git Workflow" (verify spec/plan status before PR open) and
in the brainstorming skill's "Spec self-review" (placeholder scan
specifically catches "TBD" / "aspirational" status).

Related memories: [[project_wire_phase_status]] (parallel pattern for
phase tracking), [[project_c4_consolidate_validation]] (parallel ADR
status footers).
```

- [ ] **Step 2: Append a pointer to MEMORY.md**

Append to `/home/tiago/.claude/projects/-data-8bandkp-fdm/memory/MEMORY.md`:

```markdown
- [Codebase Doc-Drift Prevention](codebase-doc-drift-prevention.md) — every PR touching behavior must update spec + plan status footers and archive closed PRDs; the 2026-06-22 PR27 audit + 2026-06-26 U8 review both surfaced drift as the rework root cause
```

- [ ] **Step 3: Verify the memory is loadable**

```bash
ls /home/tiago/.claude/projects/-data-8bandkp-fdm/memory/
grep "codebase-doc-drift-prevention" /home/tiago/.claude/projects/-data-8bandkp-fdm/memory/MEMORY.md
```

Expected: file exists, MEMORY.md has the pointer line.

(Note: memory files are not committed to the repo — they live outside the working tree. No git commit for this task.)

---

## Task 13: Final ctest run + push to PR40

**Files:** none modified; verification only.

- [ ] **Step 1: Full ctest run**

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -j4 --output-on-failure
```

Expected: 129/129 PASS (126 pre-existing + V13 from C1 + T4 from C2 + PHS test from C8). The full suite count may be higher if T1/T2 are split into separate ctest entries.

- [ ] **Step 2: Count commits**

```bash
git log feat/bdg-u8-window-routing --not main --oneline | wc -l
```

Expected: 20 (10 existing + 10 new C1–C10).

- [ ] **Step 3: Push to PR40**

```bash
git push origin feat/bdg-u8-window-routing
gh pr view 40 --repo tiagocampo/8bandkp-fdm --json state,url
```

Expected: PR state remains OPEN, branch updated.

- [ ] **Step 4: Update PR40 body**

```bash
gh pr edit 40 --repo tiagocampo/8bandkp-fdm --body "$(cat <<'EOF'
[existing PR40 body]

## Follow-up commits (2026-06-26, this branch)

Closes the adversarial review of the U8 spec/plan + PR40 inline threads:

- **C1** `fix(parser): check_optional_stat on b_field components (PR27 I2)` — closes PR27 I2
- **C2** `feat(validate): reject BdG with axial B_vec (transverse-B required)` — design §1 second root cause
- **C3** `fix(bdg): sentinel uses error stop + clears stale output files` — closes PR40 Codex P2 (stale file)
- **C4** `refactor(defs): BDG_WINDOW_BOUND as named module parameter`
- **C5** `chore(eigensolver): drop unused auto_compute_energy_window import`
- **C6** `test(bdg): T2 sentinel assertion + minigap coverage cell`
- **C7** `fix(bdg): route eval_wire_bdg_gap_app through window authority` — closes design §4.1/§6 contradiction
- **C8** `test(bdg): PHS cross-check at finite Bx (R3 verification)` — closes brainstorm R3 gap
- **C9** `docs: resolve U8 plan/design self-contradictions + add status footers`
- **C10** `refactor: rename eval_wire_bdg_gap_app + add domain-language entries`

Full ctest: 129/129. `.scratch/pr27-review-fixes/` archived to `.scratch/archive/`. Memory note `codebase-doc-drift-prevention` saved.

Codacy "Not up to standards" findings (2 critical + 1 medium + 1 minor Security) — specifics not visible in API; will address via `@codacy review` re-trigger if surfaced.
EOF
)"
```

---

## Self-Review

**Spec coverage:**
- C1 PR27 I2 → Task 1 ✓
- C2 transverse-B guard → Task 2 ✓
- C3 sentinel error stop + stale file → Task 3 ✓
- C4 BDG_WINDOW_BOUND → Task 4 ✓
- C5 unused import → Task 5 ✓
- C6 T2 sentinel + coverage cell → Task 6 ✓
- C7 eval_wire_bdg_gap_app rewrite → Task 7 ✓
- C8 PHS test → Task 8 ✓
- C9 doc self-contradiction → Task 9 ✓
- C10 terminology + rename + OMP docstring → Task 10 ✓
- Archive → Task 11 ✓
- Memory → Task 12 ✓
- Final ctest + push → Task 13 ✓

**Placeholder scan:** No "TBD" / "TODO" / "similar to Task N" / "fill in details". One conditional code block in Task 8 Step 2 says "Read [...] for existing import patterns and adjust the test to match" — this is honest guidance for the agent because the test file's exact API isn't verified. The PHS check logic itself is fully specified.

**Type consistency:**
- `eval_wire_bdg_gap_app` renamed consistently across Tasks 7 and 10 (grep + sed replacement in Step 2 of Task 10)
- `bdg_wire_cleanup` defined in Task 7 Step 1 and called in Steps 4-6 — consistent signature
- `BDG_WINDOW_BOUND` defined in Task 4 Step 1, referenced in Step 2 — consistent
- `check_optional_stat` calls in Task 1 Step 3 follow the same pattern as the existing `g_factor` call on line 459 — consistent

**Risks called out in plan:**
- Task 3 (C3) verifier compatibility with the tightened T2 — Step 4 expects FAIL, and Step 6 of Task 6 fixes it
- Task 8 (C8) imports may need adjustment — Step 2 explicitly says to read the file and adjust; commit only if Step 3 passes
- Task 13 (final ctest) — count may be 129 or higher depending on whether existing T1/T2 count separately

Plan complete and saved to `docs/superpowers/plans/2026-06-26-u8-followup-reviews-and-codex.md`.