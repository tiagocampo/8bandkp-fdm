---
title: "U8 — BdG wire path: corrected μ + window routing + non-zero-minigap regression"
date: 2026-06-21
status: design (brainstorm-approved 2026-06-21)
scope_choice: "B (lock-in + window routing)"
related:
  - docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md  # U8 task
  - docs/solutions/best-practices/2026-06-21-fd-nyquist-spurious-tail.md
  - /tmp/handoff-bdg-majorana-2026-06-21.md
---

# U8 — Lock in the confirmed BdG fix + close the window bypasses

## 1. Problem & reframe

The BdG/Majorana "all-zero" minigap (zero for all B) on the sparse wire path is
**not** a build, solver, or hole-block bug. It has two root causes, both verified
in code and now confirmed experimentally:

1. **μ in the solver's gap.** The solver zeroes energy at EV=0 (EC=+Eg). The
   all-zero lineage config placed μ mid-gap → no Fermi surface → no near-μ BdG
   states → `min_gap=0` for all B. The 06-15 prescription "μ≈EC=0.719
   (parameters.f90 scale)" was **refuted**; μ must sit at a real *wire* subband
   edge.
2. **Gauge/SOC absent with B along the wire.** Peierls orbital coupling
   (`add_peierls_coo`, `src/physics/magnetic_field.f90`) early-returns when
   `Bx=0`; the `[bdg] gauge` field is parsed (`src/io/input_parser.f90:718`) but
   never consumed (dead knob). A symmetric single-material wire has ~0 SIA-Rashba.
   So `B_vec=[0,0,Bz]` yields diagonal Zeeman only → no topological gap. **Fix:
   transverse `B_vec=[Bx,0,0]`** so Zeeman (|B|) **and** Peierls (Bx) both
   activate.

A third, compounding issue made diagnosis painful: the **auto-window fallback**
in `run_bdg_wire` retries with `auto_compute_energy_window`, which on the BdG
matrix returns the Gershgorin ±~78 eV (the FD-Nyquist tail) — slow (~150 s) and
unphysical. This must be removed on the BdG path.

## 2. Confirmed experimental result (2026-06-21, core/shell InAs/GaAs wire)

Subband edges located directly via a DENSE-FULL normal solve reading
`output/eigenvalues.dat`: **conduction bottom = +0.6586 eV**, valence top =
−0.5862 eV (gap 1.245 eV). (The slow BdG μ-scan is superseded by this fast,
precise, convention-agnostic readout.)

B-sweep (μ=0.6601 eV → μ_eff≈1.5 meV, g_factor=15, δ₀=2·10⁻⁴, ±5 meV window,
transverse B) — **open → close → reopen**:

| Bx (T) | 0 | 1 | 2 | 2.5 | **2.8** | 3.0 | 3.5 | 4.5 | 6 |
|---|---|---|---|---|---|---|---|---|---|
| minigap (meV) | 2.91 | 1.21 | 0.49 | 0.54 | **0.019** | 0.33 | 1.20 | 2.93 | 5.53 |

Gap closes to **0.019 meV ≈ 0.1·δ₀** at Bx=2.8 T (below the pairing gap → a
genuine closing, not an avoided crossing), then reopens.

μ-shift (gold-standard Lutchyn–Oreg signature, B_crit ∝ μ_eff):

| μ_eff (meV) | μ (eV) | B_crit (T) | gap at closing (meV) |
|---|---|---|---|
| 0.5 | 0.6591 | 0.6 | 0.061 |
| 1.5 | 0.6601 | 2.8 | 0.019 |
| 3.0 | 0.6616 | 3.6 | 0.102 |

B_crit increases monotonically with μ_eff (sub-linear at the top end because
μ=0.6616 eV nears the second Kramers subband at 0.6625 eV — multi-band regime;
pinning the exact phase boundary uses the Z2/Pfaffian invariant, deferred to U10).

Scratch artifacts: `/tmp/bdg-probe/coreshell_bsweep/`, `/tmp/bdg-probe/coreshell_normal/`;
drivers `/tmp/bdg-probe/coreshell_uscan.py`, `coreshell_bsweep.py`.

## 3. Scope (approved: B — lock-in + window routing)

- **Bake in the corrected prescription** (μ at conduction edge + transverse B) via
  a config + a regression test.
- **Route** the two wire-BdG eigensolve sites through `apply_solver_window`
  (ADR 0005 / KTD6) and **remove the auto-window fallback** on the BdG path.
- **Add a validation guard** rejecting explicitly Gershgorin-scale windows for BdG.
- **Correct the plan's U8 text** (μ≈0.719 → corrected prescription).

## 4. Design

### 4.1 Window routing (KTD6) — `src/apps/main_topology.f90`
In `run_bdg_wire` (~L490-525) and `eval_wire_bdg_gap_app` (~L1250-1265): replace
the inline `eigen_cfg_local%emin/emax` assignment + manual fallback with
`call apply_solver_window(H_bdg_csr, user_emin, user_emax, emin, emax)` (the
`asw_single` variant, `src/math/eigensolver.f90:561`). Pass the physics window
(default ±50·δ₀, or `cfg%solver%emin/emax` if set) as the user override, which
`asw_single` honors verbatim (eigensolver.f90:566-570). Behavior-preserving for
the window value; makes `apply_solver_window` the single authority and removes
the manual `auto_compute_energy_window` bypass. The third KTD6 bypass
(`compute_spectral_function_wire`) is **deferred to U9**.

### 4.2 Fallback removal (i) — `src/apps/main_topology.f90`
On FEAST finding 0 eigenvalues in the tight window: **delete the
`auto_compute_energy_window` retry** (current `run_bdg_wire` L517-524). Instead
emit a warning naming μ and the window, and set `result%min_gap = -1.0_dp`
(sentinel — `run_bdg_wire` already has this branch at L528). Make
`eval_wire_bdg_gap_app` consistent: it currently `error stop`s on failure
(L1270) — change to sentinel + warn so a B-sweep survives one μ-in-gap point
without aborting. The 0-found case is almost always a μ-in-gap or wrong-window
config error; fail loudly, never silently return tail states.

### 4.3 Validation guard (iii) — `src/core/defs.f90` `validate_semantic`
New check (next free S-number, slotted next to the existing BdG checks S3/S4 at
L939-949): if `cfg%bdg%enabled` **and** the solver window is explicitly set
(`cfg%solver%emin /= 0 .or. cfg%solver%emax /= 0`), require
`max(abs(emin), abs(emax)) <= BDG_WINDOW_BOUND`; else `error stop` with a
contextual message. `BDG_WINDOW_BOUND` is a rejection ceiling (a `parameter`,
~1.0 eV) — physics BdG windows are ≤~50 meV, so 1 eV rejects only Gershgorin-scale
(±tens of eV) windows while never touching a legitimate one. Defense-in-depth
against future configs
that explicitly set a Gershgorin-scale BdG window. Reads only existing fields →
**no new TOML field** (ADR 0002 ✓). *(Editing `defs.f90` is approval-gated per
CLAUDE.md; flag at implementation time.)*

### 4.4 Config + regression test
- New config `tests/regression/configs/wire_inas_gaas_bdg_topological.toml`:
  core/shell InAs/GaAs, μ=0.6601 eV, transverse `B_vec=[Bx,0,0]`, g_factor=15,
  δ₀=2·10⁻⁴, ±5 meV window, strain-enabled (the confirmed working lineage).
- New verifier `tests/integration/verify_wire_bdg_topological.py` + shell wrapper
  `tests/integration/test_wire_bdg_topological.sh`: run `topologicalAnalysis` at
  Bx ∈ {0, 2.8, 5} T, parse `Min gap`, assert **robust shape**:
  `minigap(2.8) < 0.5·minigap(0)` **and** `minigap(5) > minigap(2.8)`
  (not exact B_crit — robust to small numerical drift). Add a `# COVERAGE:` cell
  `observable=minigap geometry=wire material=InAs` and a ctest entry (label
  `regression` + `bdg`). ~3 solves × ~50 s ≈ 2-3 min ctest (within norms).

### 4.5 Plan text correction — `docs/plans/2026-06-14-001-…md`
Update the U8 section: replace "μ ≈ EC=+0.719 eV or EV=−0.8 eV" with "μ at the
conduction subband edge in the solver's EV=0/EC=Eg convention, located via a
DENSE-FULL normal solve (`bandStructure`, read `output/eigenvalues.dat`); for the
core/shell wire the edge is +0.659 eV. Use transverse `B_vec=[Bx,0,0]` so Peierls
orbital coupling activates (B along the wire gives Zeeman only)."

## 5. Tests (TDD)

- **T1 — lock-in (passes before & after code changes):** the §4.4 regression
  test. Proves the fix and guards it against regressions.
- **T2 — fallback removal (RED→GREEN):** a μ-in-gap BdG config (e.g. μ=0.0) run
  through `topologicalAnalysis` must return the `min_gap=-1` sentinel and emit
  the warning, **not** a tail-spectrum minigap. Fails today (returns ~42 meV
  from the ±78 eV tail), passes after §4.2.
- **T3 — validate guard (RED→GREEN):** add to
  `tests/integration/test_validate_rejects_bad_configs.sh` a case: BdG mode +
  `[solver] emin=-70, emax=70` → `error stop`. Fails today, passes after §4.3.

## 6. Out of scope (deferred)

- Third KTD6 bypass `compute_spectral_function_wire` → **U9**.
- Dead `[bdg] gauge` knob cleanup + principled `g_factor` → own small task.
- Pure per-point BdG evaluator extraction → **U2**.
- (B, μ) phase diagram regeneration → **U10/U11** (needs the Pfaffian invariant).

## 7. Risks / approval gates

- **`defs.f90` touched** (§4.3 validate guard) — CLAUDE.md flags `defs.f90`
  derived-type changes as approval-gated. This adds a check to an existing
  validator (not a derived-type change), but flag explicitly at implementation.
- **ADR compliance:** 0002 (no new TOML fields — guard reads existing fields ✓),
  0005 (window authority — routing complies ✓), 0001 (no new polymorphic types ✓).
- **ctest cost:** ~2-3 min for T1 (within the range of existing standard-star
  tests; ensure OMP thread cap per the ctest gotcha).
- **Behavior change visibility:** removing the auto-fallback changes observable
  behavior for μ-in-gap configs (sentinel instead of tail minigap). This is
  intended and covered by T2; document in the commit message.
