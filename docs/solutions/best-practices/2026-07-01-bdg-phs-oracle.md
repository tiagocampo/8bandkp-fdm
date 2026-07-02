---
module: bdg_hamiltonian
tags: [PHS, particle-hole-symmetry, oracle, regression-test]
problem_type: false-pass-prevention
component: bdg_observables
---

# BdG PHS Oracle as Teeth-First Regression Baseline

## Problem

The BdG/TopologicalSC machinery historically asserted PHS indirectly via ±E pairing tests and gap-closing heuristics. Both can pass even when the particle-hole symmetry is silently broken.

## Solution

Issue 02 added a direct numerical PHS oracle: `C H_BdG(k) C^{-1} + H_BdG(-k) ≈ 0` with `C = τ_x K`. The oracle has 4 sub-tests across all field combinations (Zeeman on/off, Peierls on/off).

Pre-Issue-03, the wire/dense-QW builders used different hole-block conventions (`-H0^T(+k)` vs `-conjg(H0(-k))`); the oracle's relative residual of ~0.12578 flagged this divergence. After Issue 03 unifies both builders on `-conjg(H0(-k))` per ADR 0007, residual drops to round-off.

## When to use

Add a PHS oracle to any non-Hermitian physical symmetry (chiral, charge-conjugation, time-reversal). The oracle's "teeth-first" property — proven to catch the bug it was added for — is the regression test design pattern of choice for symmetry invariants.