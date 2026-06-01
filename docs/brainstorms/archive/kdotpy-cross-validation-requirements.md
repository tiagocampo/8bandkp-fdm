---
date: 2026-05-17
topic: kdotpy-cross-validation
---

# kdotpy Cross-Code Validation Pipeline

## Summary

Build an automated comparison pipeline between our Fortran 8-band k.p code and kdotpy (Wuerzburg, v1.3.1, Python, plane-wave discretization). Map shared capabilities, port materials bidirectionally (III-V to kdotpy, II-VI to our code), and run layered verification from bulk (exact agreement) through confined systems (convergence studies). All discrepancies are investigated and resolved — fixes applied to whichever code is wrong, model differences documented.

---

## Problem Frame

Our 8-band k.p code has 91 tests validating against analytical formulas and literature values, but no cross-code comparison. When the 8-band CB effective mass deviates 29% from Vurgaftman's 2-band value, our tests tolerate that deviation — but also tolerate a wrong 8-band implementation. An independent 8-band code using different discretization (plane-wave vs finite differences) provides ground truth that 2-band analytical formulas cannot. The CECAM 2024 workshop identified this gap: no shared inter-code benchmark suite exists for k.p codes. kdotpy is the closest open-source comparator — same Vurgaftman parameters, different discretization, actively maintained, and now installed locally.

---

## Requirements

**Infrastructure**

- R1. kdotpy v1.3.1 installed in `validation/kdotpy_env/` Python virtual environment (DONE). All comparison scripts run from this venv.
- R2. `validation/` directory contains structured subdirectories for each physics area: `bulk/`, `qw/`, `wire/`, `landau/`, `gfactor/`, `strain/`, `berry/`, `selfconsistent/`, plus `shared/` for param_mapper.py and run_comparison.py. Each subdirectory holds input configs (both formats), output data (both codes), and comparison results.
- R3. Shared parameter mapping module (`validation/shared/param_mapper.py`) translates between our parameter format (eV, Å, Vurgaftman naming) and kdotpy's (meV, nm, kdotpy naming). Must include a **Parameter Convention Analysis** that explicitly derives the mathematical mapping for each parameter:
  - **Standard conversions**: Eg (eV→meV ×1000), delta_SO (eV→meV ×1000), lattice constant (Å→nm /10), dielectric constant (dimensionless, direct copy), elastic constants (GPa, verify kdotpy uses same).
  - **EP → P**: kdotpy's `P` is the Kane matrix element in meV·nm. Our `EP` is in eV. Formula: `P = sqrt(EP_eV * 1000 * hbarm0_kdotpy)` where `hbarm0_kdotpy` is kdotpy's hbar²/(2m0) in meV·nm². Must use kdotpy's constant value, not ours, to avoid unit mismatch.
  - **A → F**: Our `paramStruct` has no `F` (remote-band Kane parameter). Our `A = 1/m_eff` is the Vurgaftman effective mass, which INCLUDES P-coupling. kdotpy's `F` represents ONLY the remote-band contribution: `F = (1/m_eff - EP/(Eg + delta_SO)) / (2 * hbarm0) - 1/2` (derived from the 8-band relation `1/m_eff = 1 + 2F + EP*(2/(Eg+delta_SO) + 1/(Eg))`). This derivation must be verified against both codes' CB diagonal matrix elements before bulk dispersion comparison (R8). WRONG F MAPPING = WRONG BULK DISPERSION.
  - **kappa**: Our `paramStruct` has no `kappa`. Derive from Luttinger parameters: `kappa = -(gamma3 - 2*gamma2 - 1) / 3` (Winkler convention). Verify against kdotpy's usage in `blocks.py` h1z (lines 278-282) where kappa appears in VB off-diagonal blocks.
  - **ge**: Our code has no explicit `ge` (CB g-factor). kdotpy uses `ge` for Zeeman. Default `ge = 2.0` (free electron). For III-V materials with strong spin-orbit, the effective ge may differ. Start with `ge = 2.0` and adjust if bulk Zeeman comparison (R9) fails.
  - **Strain**: `strain_C1 = ac * 1000` (eV→meV), `strain_Dd = av * 1000`, `strain_Du = -1.5 * b_dp * 1000`, `strain_Duprime = -0.5 * sqrt(3) * d_dp * 1000`. Verify sign conventions against kdotpy's `hstrain` block (blocks.py line 1653).
- R4. Shared comparison runner (`validation/run_comparison.py`) executes both codes with matched configurations, parses outputs, and produces structured comparison reports (JSON + human-readable summary).

**Material porting**

- R5. kdotpy INI material files created for all our standard-star materials: GaAs, InAs, InSb, AlAs, GaSb, AlSb, InP, AlGaAs alloys (20%, 30%), InAsSb alloys. All III-V materials are absent from kdotpy — every file must be authored from scratch. Parameters sourced from our `parameters.f90` (Vurgaftman 2001) via the param_mapper (R3). INI file must include all required kdotpy fields: Ec, Ev, P (via `sqrt(EP * hbarm0)` using kdotpy's constant), F (derived from A/m_eff), gamma1/2/3, kappa, ge (default 2.0), delta_so, a (nm), strain_C1/Dd/Du/Duprime, diel_epsilon, elasticity_c11/c12/c44. Each file validated by running kdotpy bulk calculation and comparing band gap at Gamma against our code's band gap to < 0.01 meV.
- R6. New entries added to our `parameters.f90` for II-VI materials: HgTe, CdTe, HgCdTe (x=0.3), CdZnTe. Parameters sourced from kdotpy's default material database. Requires approval per CLAUDE.md boundary on parameter modifications. Each material entry must cite the primary published source (kdotpy references: Pfeuffer-Jeschke PhD thesis, Becker et al. PRB 62, 10353), not copied from kdotpy's database without independent verification against those primary sources.

**Bulk band structure**

- R7. Bulk 8-band Hamiltonian comparison at Gamma (k=0): all 8 eigenvalues must agree to < 0.01 meV for all ported materials. This is an exact matrix comparison — both codes build the same 8x8 Hamiltonian at k=0.
- R8. Bulk dispersion E(k) along [100], [110], [111] directions. Eigenvalues compared at 20+ k-points per direction. Agreement expected to < 0.1 meV for |k| < 0.1 nm^-1 (exact k.p, no discretization). Effective masses from parabolic fit compared to < 1%. This tolerance assumes parameter mapping (R3) is correct — if R7 passes but R8 fails, the root cause is the F/kappa mapping.
- R9. Bulk with magnetic field (Zeeman only): compare Zeeman-split eigenvalues at Gamma for B along x, y, z. Agreement to < 0.01 meV. Prerequisite: kappa and ge parameter mapping validated (from R3 convention analysis). If our Zeeman formalism differs from kdotpy's, document the mathematical equivalence before comparing. Conditional on R3 resolution — skip if kappa/ge cannot be reliably derived from our parameters.

**QW subbands**

- R10. QW subband energies E(k_par=0) for standard-star systems: GaAs/AlGaAs, InAs/GaSb, InAs/GaAs strained QW. Compare at multiple well widths (5, 7, 10, 15 nm). Agreement target: < 1 meV at fine grid spacing (both codes converged).
- R11. QW dispersion E(k_par) for in-plane sweep. Compare subband curvature (effective masses) to < 2%. Non-parabolicity effects included.
- R12. Convergence study: run both codes at multiple grid resolutions, extract continuum-limit values via Richardson extrapolation, verify both codes converge to the same extrapolated value. Minimum 4 grid spacings per code. Our code: FD orders 2, 4, 6, 8 at fixed grid spacing, plus FDstep sweep at fixed order. kdotpy: zres = 0.5, 0.25, 0.125, 0.0625 nm. Both sets extrapolated independently; extrapolated values compared. Grid resolution matching for R10/R11 means: both codes at their respective converged limits (as determined by this study), not at any particular resolution.

**Wire subbands**

- R13. Wire subband energies for rectangular wire geometry at k_z=0. kdotpy's 1D mode (`kdotpy-1d`) simulates a strip with hard-wall or periodic y-confinement, analogous to our wire geometry. The plane-wave vs FD discretizations have fundamentally different convergence behavior in 2D. Before tolerance comparison, run a convergence study on a single wire geometry in both codes to establish achievable agreement. Target: < 2 meV at convergence, < 5 meV at moderate resolution.
- R14. Wire dispersion E(k_z) along wire axis. Subband effective masses compared to < 5%. Same convergence caveats as R13.

**Landau levels**

- R15. Landau level fan diagram E(B) for GaAs QW with B perpendicular. Compare LL energies at 10+ B-field values (0.1 T to 10 T). Agreement target: < 1 meV at matching grid resolution.
- R16. Bulk Landau levels for GaAs. Compare LL energies at representative B values. Agreement to < 0.5 meV.

**BHZ / g-factor**

- R17. BHZ parameter extraction comparison: run kdotpy's BHZ extraction on QW systems, compare A, B, C, D, M parameters against our g-factor Lowdin partitioning results. kdotpy's BHZ uses numerical dH/dk derivatives; our g-factor uses commutator velocity matrices [r_alpha, H]. These give identical results in a complete basis but may differ by O((P/Eg)^4) in a truncated 8-band basis. Prerequisite: kappa/ge parameter mapping validated (resolves R3 convention analysis). Expected agreement: < 10% for wide QWs (> 10 nm), larger deviations acceptable for narrow QWs (< 5 nm).
- R18. QW g-factor: compare our Lowdin g-factors (g_x, g_y, g_z) against kdotpy's Zeeman splitting derivative dg/dB for standard-star QWs. Prerequisite: kappa/ge mapping validated. Agreement target: < 5%.

**Strain**

- R19. QW biaxial strain: compare strained band edge shifts (delta_Ec, delta_EHH, delta_ELH, delta_ESO) against analytical Bir-Pikus formulas. Both codes should reproduce the analytical result to < 0.1 meV. Then compare strained subband energies between codes to < 1 meV.
- R20. Strained QW subband energies: InAs/GaAs (compressive) and InAs/GaSb systems with strain enabled. Agreement to < 1 meV.

**Berry curvature and Chern numbers**

- R21. Berry curvature at k-points in QW Brillouin zone. Both codes use Kubo formula with numerical Hamiltonian derivatives. Compare Omega(k) values to < 10% (sensitive to numerical differentiation step size and eigenvalue spacing).
- R22. LL Chern number comparison: kdotpy's Chern numbers (`chernnumber_ll`) are LL-based (magnetic-field-dependent), not BZ-integrated like our Fukui-Hatsugai-Suzuki algorithm. Scope R22 to LL Chern numbers at matching B-field values (same integer result). Both codes should produce the same integer for C = 0 (trivial insulator) and C = ±1 (QHE-like) LL systems. BZ-integrated Chern comparison is not feasible — kdotpy has no BZ Chern calculator.

**Self-consistent**

- R23. Self-consistent QW: compare converged potential profiles V(z), subband energies, and charge densities n(z) for a doped GaAs/AlGaAs QW. Both codes use iterative SP loops with different mixing strategies (ours: DIIS/Pulay, kdotpy: dynamic time stepping). SC results depend on k-parallel sampling, temperature, convergence criterion, and boundary conditions. Both codes must use matched: temperature, k_parallel grid (number and max), convergence tolerance, and boundary conditions (Dirichlet/Neumann). Before comparing absolute values, run both SC solvers on a single GaAs/AlGaAs QW and verify they converge to the SAME potential profile shape. Converged quantities compared to < 5 meV (potential), < 5 meV (subband energies), < 10% (charge density).

**Discrepancy resolution**

- R24. Every discrepancy above tolerance is logged with: observable, both values, absolute and relative difference, configuration details. A discrepancy report is generated per comparison run.
- R25. Each logged discrepancy must be resolved before proceeding to the next verification layer. Resolution options: (a) fix our code, (b) fix kdotpy config, (c) document as model difference with physics justification.
- R26. Layered verification order with dependency chain:
  1. Parameter mapping validation (R3 convention analysis) — GATE for all subsequent work
  2. Bulk at k=0 (R7) — validates material parameter porting
  3. Bulk dispersion (R8) — validates F/kappa mapping
  4. Bulk Zeeman (R9) — conditional on kappa/ge mapping
  5. QW subbands (R10-R12) — after bulk passes
  6. Wire (R13-R14), Landau levels (R15-R16) — after QW passes
  7. g-factor/BHZ (R17-R18) — after QW + kappa/ge mapping validated
  8. Strain (R19-R20) — after QW passes
  9. Berry/Chern (R21-R22) — after Landau levels pass
  10. Self-consistent (R23) — after QW and strain pass

---

## Success Criteria

- Bulk eigenvalues at k=0 agree to < 0.01 meV for all ported materials (exact matrix comparison).
- Bulk dispersion at non-zero k agrees to < 0.1 meV (validates F/kappa mapping).
- QW subband energies agree to < 1 meV at converged grid resolution (validated by Richardson extrapolation).
- All discrepancy reports are resolved (fixed or documented with physics justification).
- The comparison pipeline is reproducible: a single command reruns the full comparison suite.
- New II-VI materials in our code pass the same bulk verification as existing materials.
- F/kappa parameter derivation is validated — bulk dispersion agreement confirms the convention analysis.

---

## Scope Boundaries

- Physics only our code has (no kdotpy comparator): BdG/Majorana, Z2 (Fu-Kane parity), optical absorption/gain/spontaneous emission, ISBT, exciton binding energy, scattering rates, spin-resolved spectra, LDOS via Green's functions, topological phase diagrams, spectral functions.
- Physics only kdotpy has (not in our comparison scope): 6-band mode, BIA/Dresselhaus, exchange interaction (Mn-doped), cubic Zeeman, arbitrary crystal orientation, lattice regularization, symbolic Hamiltonian, GPU solvers.
- Code quality, style, or performance comparison between codes.
- Publishing a joint benchmark paper (potential future work).
- Parameter optimization or fitting to experimental data.

---

## Key Decisions

- Bidirectional material porting (both directions): maximizes testable physics, validates parameter mapping itself.
- Strict discrepancy protocol (investigate all): avoids tolerance masks hiding real bugs, at cost of more investigation effort.
- Layered verification starting from bulk: bulk is exact (8x8 matrix, no discretization), so parameter mapping bugs are caught before confined systems introduce convergence complexity.
- kdotpy's plane-wave vs our FD: two independent discretizations converging to the same physics provides stronger validation than either code alone.

---

## Dependencies / Assumptions

- kdotpy v1.3.1 is stable and its physics is correct for the tested regimes. We validate agreement, not correctness of either code in isolation.
- Both codes implement the standard 8-band zinc-blende k.p Hamiltonian with the same physical terms (k.p coupling, SOC, strain, Zeeman), but with different basis ordering (kdotpy: CB-first, ours: VB-first) and potentially different treatments of remote-band parameters (F, kappa). The permutation doesn't affect eigenvalues.
- Deformation potential convention: kdotpy uses D_u (= -1.5*b) and D_u' (= -0.5*sqrt(3)*d). Our code uses b and d directly. The param_mapper handles this conversion. Explicit formulas in R3.
- kdotpy uses meV/nm units; our code uses eV/Å internally. Conversion factors: 1 eV = 1000 meV, 1 nm = 10 Å.
- Physical constants: both reference CODATA 2014. Values should agree to full precision within their respective unit systems, but floating-point representations may differ by a few ULP. The param_mapper must use kdotpy's constant values (not ours) when computing derived parameters for kdotpy INI files (e.g., P from EP).
- kdotpy's plane-wave discretization uses a single Fourier mode per z-layer. Our FD uses multi-point stencils. Both converge to continuum limit but at different rates — convergence studies (R12) are essential for QW/wire before tolerance comparison is meaningful.
- kdotpy's BHZ extraction and Chern numbers are LL-based (magnetic-field-dependent), while our Chern numbers use BZ integration. Only LL Chern comparison is feasible (R22).

---

## Outstanding Questions

### Resolve Before Planning

- [Affects R3, R8, R9, R17][Technical] **F/kappa parameter derivation** — Our code has no F or kappa fields. The derivation `F = (1/m_eff - EP*(2/(Eg+delta_SO) + 1/Eg)) / (2*hbarm0) - 1/2` and `kappa = -(gamma3 - 2*gamma2 - 1)/3` must be verified by comparing bulk Hamiltonian matrix elements between both codes at non-zero k. If the derivation is wrong, bulk dispersion (R8) will fail and all downstream comparisons are blocked. This is the highest-risk item.
- [Affects R1][Technical] **kdotpy smoke test** — Run `kdotpy-bulk` with a trivial configuration (HgTe) to verify the installation works end-to-end before building any comparison infrastructure.

### Deferred to Planning

- [Affects R5][Needs research] Whether kdotpy INI material files support all Vurgaftman parameters we need (F, kappa, q) or if some require code-level configuration. Check kdotpy material parameter schema completeness in `materials/materials.py`.
- [Affects R12][Technical] Richardson extrapolation implementation details: exact grid spacing range, whether to use our FD order sweep vs FDstep sweep, how to handle kdotpy's different convergence order.
- [Affects R17][Technical] BHZ vs commutator velocity correspondence: numerical dH/dk vs -i[r,H] give identical results in complete basis but differ by O((P/Eg)^4) in truncated 8-band. Need to quantify this for our test systems.
- [Affects R13][Technical] Wire grid resolution matching: plane-wave in kdotpy vs FD in our code have fundamentally different 2D convergence behavior. Pilot study needed before setting tolerance.
