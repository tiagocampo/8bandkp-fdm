---
title: "Wire topology/BdG/spectral paths silently dropped strain (copy-paste omitted compute_strain)"
date: "2026-06-15"
category: logic-errors
module: src/physics/wire_setup.f90
problem_type: logic_error
component: service_object
severity: critical
symptoms:
  - "Strained wire topology/BdG/spectral calculations produced spectra identical to the unstrained run (zero strain shift)"
  - "No runtime error -- strain was silently ignored on the wire topology/BdG/spectral paths"
  - "cfg%strain_blocks%delta_Ec never allocated on these paths, so the strain insertion in ZB8bandGeneralized was gated off"
  - "Only the canonical bandStructure path (simulation_setup case 'wire') applied strain correctly"
root_cause: copy_paste_omission
resolution_type: refactor_fix
tags: [wire, strain, bir-pikus, topology, bdg, spectral, copy-paste, encapsulation, architecture-deepening]
related_components:
  - src/apps/main_topology.f90
  - src/physics/green_functions.f90
  - src/core/simulation_setup.f90
  - src/physics/hamiltonian_wire.f90
  - tests/integration/verify_wire_bdg_strain_shift.py
---

# Wire topology/BdG/spectral strain-omission bug (#04)

## Problem

The 2D-wire confinement initialization (`confinementInitialization_2d` building
`profile_2d`/`kpterms_2d`) was copy-pasted across four sites:

1. `simulation_setup_init` `case('wire')` — **canonical** (applies strain).
2. `run_bdg_wire` in `main_topology.f90` — copy-pasted, **strain skipped**.
3. `eval_wire_bdg_gap_app` in `main_topology.f90` — copy-pasted, **strain skipped**.
4. `compute_spectral_function_wire` in `green_functions.f90` — copy-pasted, **strain skipped**.

The canonical path (site 1) follows `confinementInitialization_2d` with, when
`cfg%strain%enabled`, `compute_strain` + `compute_bir_pikus_blocks`, which
populates `cfg%strain_blocks`. The strain insertion in `ZB8bandGeneralized`
is gated on `allocated(cfg%strain_blocks%delta_Ec)`, so the copy-pasted paths
(2-4) silently dropped strain: `cfg%strain_blocks` was never populated and the
strain insertion was skipped. Strained wire topology/BdG/spectral calculations
were physically wrong (identical to the unstrained run).

A fifth site, `run_qshe_wire`, also copy-pasted the init but built a **BHZ**
4-band Hamiltonian (`build_bhz_wire_hamiltonian`) that never reads
`profile_2d`/`kpterms_2d` — the init there was dead boilerplate.

## Root cause

Copy-paste of the wire init boilerplate without the strain step. The strain
gating in `ZB8bandGeneralized` (`if (allocated(cfg%strain_blocks%delta_Ec))`)
made the omission silent rather than catastrophic.

## Fix

Encapsulated wire init + cleanup in one type, `wire_setup` (new module
`src/physics/wire_setup.f90`), whose `wire_setup_init` runs the **same strain
step as the canonical path** (`confinementInitialization_2d` +
`compute_strain` + `compute_bir_pikus_blocks` when `cfg%strain%enabled`).
Routed the buggy sites through it:

- `run_bdg_wire` → `wire_setup_init` + `wire_setup_free` (fixes the bug).
- `eval_wire_bdg_gap_app` → `wire_setup_init` + `wire_setup_free` (fixes the bug).
- `compute_spectral_function_wire` → `wire_setup_init` + `wire_setup_free` (fixes the bug).
- `run_qshe_wire` (BHZ) → dead `confinementInitialization_2d` + cleanup removed
  (BHZ is a 4-band model; 8-band k.p strain does not apply).

The fix is "by construction": the strain step cannot be omitted because it
lives inside `wire_setup_init`, not at each call site.

A second init variant, `wire_setup_adopt_precomputed`, wraps caller-supplied
already-strain-applied profile/terms for sink-style callers that receive that
data as arguments (no re-init, no strain step).

## Evidence (regression test)

`tests/integration/verify_wire_bdg_strain_shift.py` runs
`topologicalAnalysis` (mode=bdg) on a strained InAs/GaAs core/shell wire
twice — once with `[strain]` enabled, once with the block stripped — and
asserts the BdG spectra differ.

Before the fix: `max |E_strained - E_unstrained| = 0.000000e+00 eV` (bug).
After the fix:  `max |E_strained - E_unstrained| = 3.923149e-02 eV` (39 meV
Bir-Pikus shift, physically correct for InAs/GaAs lattice mismatch).

## Golden handling

No golden regression config exercised the strained wire-topology/BdG/spectral
path (all topology/spectral golden configs have strain disabled; the strained
wire configs run only through the canonical `bandStructure` path which already
applied strain). So no golden reference needed updating — the ideal outcome.

## What did NOT change

- Bir-Pikus sign convention, `get_strain_table` / `compute_bp_scalar` SSOTs.
- Basis ordering, k.p block table, material parameters.
- The canonical `simulation_setup` wire path (already correct).
- BHZ 4-band model (strain in the 8-band k.p sense does not apply).

## Design constraints honored (ADRs 0001 + 0005)

- Concrete type (no polymorphic builder hierarchy).
- New type in a physics module (`wire_setup.f90`), NOT in approval-gated `defs.f90`.
- F2018, `error stop` with message, `private` default + explicit `public ::`.
- Idempotent `free` via `was_freed` flag; finalizer delegates to `wire_setup_free`.
