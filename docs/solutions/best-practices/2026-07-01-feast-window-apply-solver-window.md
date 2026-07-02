---
module: eigensolver
tags: [FEAST, energy-window, KTD6, ADR-0005, encapsulation]
problem_type: manual-window-bypass
component: apply_solver_window
---

# FEAST Window Routed Through apply_solver_window (KTD6 Closure)

## Problem

ADR 0005 mandates: "one stable window per sweep via apply_solver_window". New sweep modes (e.g., bdq_spectral) needed FEAST windows for two consumers. Both could fall back to manual windows if `cfg%solver%emin/emax` weren't set.

## Solution

`compute_spectral_function_wire` calls `apply_solver_window(E_arr, cfg%solver%emin, cfg%solver%emax, eigen_cfg%emin, eigen_cfg%emax)` rather than constructing windows inline. This routes through the generic interface (`asw_envelope` / `asw_single` / `asw_evals`) which dispatches on the solver type.

## Test for KTD6

A test was added to assert that calling `compute_spectral_function_wire` with empty `cfg%solver%emin/emax` produces the SAME auto-window result as calling it with the explicit window — proving the auto-window path is the authority, not a manual override.

## When to use

Whenever adding a new eigensolver consumer in the topology path: forward `cfg%solver%emin/emax` through `apply_solver_window` rather than computing windows inline. Encapsulation prevents per-call divergence.