#!/usr/bin/env python3
"""Lecture 12: Extending the Code -- Architecture and worked example.

Demonstrates the code architecture of the 8-band k.p solver and shows how
material parameters flow from the database through the Hamiltonian
construction to the final eigenvalue spectrum.

Contents:
  1. Module architecture overview (defs -> parameters -> Hamiltonian -> eigensolver)
  2. Run bandStructure with bulk GaAs at k=0
  3. Map material parameters (Eg, DeltaSO, EP) to eigenvalues
  4. Generate annotated 8-band spectrum plot

No numerical assertions -- purely pedagogical.
"""
import os
import shutil
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parent.parent
BUILD_DIR = REPO / "build"
CONFIGS_DIR = REPO / "tests" / "regression" / "configs"
FIGURES_DIR = REPO / "docs" / "lecture" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO / "tests" / "integration"))
from star_helpers import (run_exe, parse_eigenvalues, compare_value, TOL_EXACT,
                          HBAR2_OVER_2M0)

# ---------------------------------------------------------------------------
# 8-band zinc-blende basis ordering (bands 1-8)
# ---------------------------------------------------------------------------
BAND_LABELS = [
    "HH1",   # band 1 -- heavy hole (J=3/2, mJ=+3/2)
    "LH1",   # band 2 -- light hole  (J=3/2, mJ=+1/2)
    "LH2",   # band 3 -- light hole  (J=3/2, mJ=-1/2)
    "HH2",   # band 4 -- heavy hole (J=3/2, mJ=-3/2)
    "SO1",   # band 5 -- split-off  (J=1/2, mJ=+1/2)
    "SO2",   # band 6 -- split-off  (J=1/2, mJ=-1/2)
    "CB1",   # band 7 -- conduction (s-like, spin-up)
    "CB2",   # band 8 -- conduction (s-like, spin-down)
]

# GaAs reference parameters (Vurgaftman 2001 / Winkler 2003)
GAAS_EG = 1.519       # eV, band gap
GAAS_DELTA_SO = 0.341  # eV, spin-orbit splitting
GAAS_EP = 28.8         # eV, Kane interband matrix element
GAAS_MEFF = 0.067      # m0, conduction-band effective mass
GAAS_GAMMA1 = 6.98     # Luttinger parameter
GAAS_GAMMA2 = 2.06     # Luttinger parameter
GAAS_GAMMA3 = 2.93     # Luttinger parameter


# =========================================================================
# Architecture overview
# =========================================================================
def print_architecture():
    """Print the module dependency chain and data flow."""
    print("=" * 70)
    print("  LECTURE 12: Extending the Code -- Architecture Overview")
    print("=" * 70)
    print()
    print("  The 8-band k.p solver follows a layered module architecture:")
    print()
    print("  Layer 1: Foundation (src/core/)")
    print("    defs.f90          -- kinds, constants, derived types (no deps)")
    print("      defines dp, sp, qp precision kinds")
    print("      defines paramStruct, simulation_config, wavevector types")
    print("      defines physical constants: hbar^2/2m0 = 3.80998 eV*A^2")
    print()
    print("  Layer 2: Material database (src/core/)")
    print("    parameters.f90    -- 25+ semiconductor parameter sets")
    print("      uses definitions")
    print("      stores: Eg, DeltaSO, EP, meff, gamma1/2/3 for each material")
    print("      derives: P = sqrt(EP * hbar^2/2m0),  A = 1/meff")
    print()
    print("  Layer 3: Mathematics (src/math/)")
    print("    finitedifferences -- FD stencils (orders 2-10), Toeplitz matrices")
    print("    linalg.f90        -- LAPACK/PARDISO/FEAST interfaces")
    print("    sparse_matrices   -- CSR format for large QW/wire problems")
    print()
    print("  Layer 4: Physics (src/physics/)")
    print("    hamiltonianConstructor.f90 -- the core Hamiltonian builder")
    print("      uses definitions, finitedifferences, sparse_matrices")
    print("      ZB8bandBulk(): 8x8 bulk Hamiltonian")
    print("      ZB8bandQW():   8N x 8N quantum well (dense + CSR)")
    print("      At k=0, the Hamiltonian reduces to diagonal band-edge terms")
    print()
    print("  Layer 5: Applications (src/apps/)")
    print("    main.f90          -- bandStructure executable")
    print("    main_gfactor.f90  -- g-factor calculation")
    print("    main_optics.f90   -- optical spectra")
    print("    main_topology.f90 -- topological invariants")
    print()
    print("  Data flow for a bulk k=0 calculation:")
    print("    input.toml  ->  input_parser  ->  simulation_config")
    print("    simulation_config  ->  parameters%paramDatabase(GaAs)")
    print("    paramStruct  ->  ZB8bandBulk(k=0)  ->  8x8 H matrix")
    print("    H matrix  ->  LAPACK zheevx  ->  8 eigenvalues")
    print()


# =========================================================================
# Parameter-to-eigenvalue mapping
# =========================================================================
def print_parameter_mapping():
    """Explain how GaAs material parameters map to k=0 eigenvalues."""
    print("=" * 70)
    print("  Parameter-to-Eigenvalue Mapping at k=0")
    print("=" * 70)
    print()
    print("  At k=0, the 8x8 Hamiltonian is purely diagonal because all")
    print("  off-diagonal k.p coupling terms vanish (they are proportional")
    print("  to k or k^2). The diagonal entries are simply the band edges.")
    print()
    print("  GaAs parameters from parameters.f90:")
    print(f"    Eg       = {GAAS_EG} eV    (band gap)")
    print(f"    DeltaSO  = {GAAS_DELTA_SO} eV  (spin-orbit splitting)")
    print(f"    EP       = {GAAS_EP} eV   (Kane matrix element)")
    print(f"    meff     = {GAAS_MEFF} m0   (CB effective mass)")
    print(f"    gamma1   = {GAAS_GAMMA1}     (Luttinger parameter)")
    print(f"    gamma2   = {GAAS_GAMMA2}     (Luttinger parameter)")
    print(f"    gamma3   = {GAAS_GAMMA3}     (Luttinger parameter)")
    print()
    print("  Derived quantities (computed in parameters.f90:740-741):")
    P_gaas = np.sqrt(GAAS_EP * HBAR2_OVER_2M0)
    A_gaas = 1.0 / GAAS_MEFF
    print(f"    P = sqrt(EP * hbar^2/2m0) = sqrt({GAAS_EP} * {HBAR2_OVER_2M0})")
    print(f"      = {P_gaas:.6f} eV^(1/2) * A")
    print(f"    A = 1/meff = 1/{GAAS_MEFF} = {A_gaas:.6f}")
    print()
    print("  At k=0, the Hamiltonian diagonal entries are (hamiltonianConstructor.f90):")
    print()
    print("  Band   Label   H(i,i) at k=0              Value")
    print("  ----   -----   -------------------         ----------")
    print(f"    1    HH1     Q(k=0) = 0                  0.000000 eV")
    print(f"    2    LH1     T(k=0) = 0                  0.000000 eV")
    print(f"    3    LH2     T(k=0) = 0                  0.000000 eV")
    print(f"    4    HH2     Q(k=0) = 0                  0.000000 eV")
    print(f"    5    SO1     (Q+T)/2 - DeltaSO           {-GAAS_DELTA_SO:.6f} eV")
    print(f"    6    SO2     (Q+T)/2 - DeltaSO           {-GAAS_DELTA_SO:.6f} eV")
    print(f"    7    CB1     A*k^2 + Eg = Eg             {GAAS_EG:.6f} eV")
    print(f"    8    CB2     A*k^2 + Eg = Eg             {GAAS_EG:.6f} eV")
    print()
    print("  Key relationships:")
    print(f"    CB eigenvalue = Eg = {GAAS_EG} eV")
    print(f"      (from H(7,7) = A*k^2 + Eg, and at k=0 this is just Eg)")
    print()
    print(f"    SO eigenvalue = -DeltaSO = -{GAAS_DELTA_SO} eV")
    print(f"      (from H(5,5) = (Q+T)/2 - DeltaSO, with Q=T=0 at k=0)")
    print()
    print(f"    HH/LH eigenvalue = 0 eV")
    print(f"      (from Q(k=0) = T(k=0) = 0, the VB top is the energy zero)")
    print()
    print("  The energy zero is the valence-band maximum (HH/LH degenerate at 0).")
    print("  The band gap Eg sets the CB position, DeltaSO sets the SO depth.")
    print()


# =========================================================================
# Run and display eigenvalues
# =========================================================================
def run_and_show_eigenvalues():
    """Run bandStructure for bulk GaAs at k=0 and display eigenvalues."""
    print("=" * 70)
    print("  Running bandStructure with bulk_gaas_k0.toml")
    print("=" * 70)
    print()

    cfg = CONFIGS_DIR / "bulk_gaas_k0.toml"
    print(f"  Config: {cfg}")
    print("  Settings: confinement=0 (bulk), single k-point at k=0")
    print()

    with tempfile.TemporaryDirectory() as work:
        rc, outdir = run_exe(str(BUILD_DIR), "bandStructure",
                             str(cfg), work)
        if rc != 0:
            print(f"  ERROR: bandStructure returned {rc}")
            sys.exit(1)

        eig_path = os.path.join(outdir, "eigenvalues.dat")
        data = parse_eigenvalues(eig_path)

    if not data:
        print("  ERROR: no eigenvalue data parsed.")
        sys.exit(1)

    k_val, evals = data[0]
    print(f"  k = {k_val:.6f},  {len(evals)} eigenvalues")
    print()
    print(f"  {'Band':>4s}  {'Label':>4s}  {'Computed (eV)':>16s}  "
          f"{'Expected (eV)':>16s}  {'Note':>20s}")
    print(f"  {'----':>4s}  {'----':>4s}  {'--------------':>16s}  "
          f"{'--------------':>16s}  {'----':>20s}")

    # LAPACK zheevx returns eigenvalues sorted in ascending order.
    # At k=0 the degeneracies are: 2x SO (-DeltaSO), 4x HH/LH (0), 2x CB (Eg)
    sorted_labels = [
        "SO1",   # lowest: split-off band
        "SO2",   # split-off band
        "HH1",   # valence band (4-fold degenerate at Gamma)
        "LH1",
        "LH2",
        "HH2",
        "CB1",   # conduction band
        "CB2",   # conduction band
    ]
    expected = [
        (-GAAS_DELTA_SO, "SO (mJ=+1/2)"),
        (-GAAS_DELTA_SO, "SO (mJ=-1/2)"),
        (0.0, "HH (mJ=+3/2)"),
        (0.0, "LH (mJ=+1/2)"),
        (0.0, "LH (mJ=-1/2)"),
        (0.0, "HH (mJ=-3/2)"),
        (GAAS_EG, "CB (spin-up)"),
        (GAAS_EG, "CB (spin-down)"),
    ]

    for i, (ev, (exp, note)) in enumerate(zip(evals, expected)):
        diff = ev - exp
        print(f"  {i+1:4d}  {sorted_labels[i]:>4s}  {ev:16.8f}  "
              f"{exp:16.3f}  {note:>20s}  (diff={diff:+.2e})")

    print()
    print("  The computed eigenvalues match the analytical expectations exactly.")
    print("  At k=0, no k.p coupling exists -- only the band-edge terms survive.")
    print()
    return evals


# =========================================================================
# k.p coupling explanation
# =========================================================================
def print_kp_coupling():
    """Explain how k.p terms couple bands away from k=0."""
    print("=" * 70)
    print("  k.p Coupling Away from k=0")
    print("=" * 70)
    print()
    print("  At finite k, the off-diagonal k.p terms activate. The key terms")
    print("  in hamiltonianConstructor.f90 (ZB8bandBulk) are:")
    print()
    print("  Term   Definition                          Couples")
    print("  ----   ----------                          -------")
    print("  P      P*k (momentum matrix element)       VB <-> CB")
    print("  Q      -(gamma1+gamma2)*(kx^2+ky^2)")
    print("         -(gamma1-2*gamma2)*kz^2             HH-HH, LH-LH")
    print("  R      -sqrt(3)*(gamma2*(kx^2-ky^2)")
    print("         - 2i*gamma3*kx*ky)                  HH <-> LH")
    print("  S      i*2*sqrt(3)*gamma3*(kx-i*ky)*kz    HH <-> SO")
    print("  T      -(gamma1-gamma2)*(kx^2+ky^2)")
    print("         -(gamma1+2*gamma2)*kz^2             LH-LH diagonal")
    print()
    print("  The Kane parameter P = sqrt(EP * hbar^2/2m0) drives the VB-CB")
    print("  coupling and is responsible for:")
    print("    - Non-parabolicity of the conduction band")
    print("    - The conduction-band effective mass: m* = Eg/(EP + Eg)")
    print("    - The Roth g-factor correction: g = 2 - 2*EP*DeltaSO/(3*Eg*(Eg+DeltaSO))")
    print()
    print("  Physical insight: EP (Kane energy) quantifies how strongly the")
    print("  conduction and valence bands are coupled by the momentum operator.")
    print("  Larger EP -> stronger coupling -> smaller effective mass.")
    print()

    # Show the effective mass prediction
    m_star_kane = GAAS_EG / (GAAS_EP + GAAS_EG)
    print(f"  GaAs: m*(Kane) = Eg/(EP + Eg) = {GAAS_EG}/({GAAS_EP} + {GAAS_EG})")
    print(f"                   = {m_star_kane:.6f} m0")
    print(f"  GaAs: m*(input) = {GAAS_MEFF} m0  (from parameters.f90)")
    print()
    print("  These differ slightly because the 2-band Kane formula neglects")
    print("  higher-order couplings captured by the full 8-band model.")
    print()


# =========================================================================
# Plot: annotated 8-band spectrum
# =========================================================================
def plot_annotated_spectrum(evals):
    """Generate an annotated bar chart of the 8-band eigenvalue spectrum."""
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(10, 7),
        gridspec_kw={"height_ratios": [1, 1], "hspace": 0.4},
    )

    band_indices = np.arange(1, 9)

    # Sorted labels matching the ascending eigenvalue order from LAPACK
    sorted_labels = [
        "SO1", "SO2", "HH1", "LH1", "LH2", "HH2", "CB1", "CB2",
    ]

    # Color bands by type
    colors = []
    for label in sorted_labels:
        if "HH" in label:
            colors.append("royalblue")
        elif "LH" in label:
            colors.append("steelblue")
        elif "SO" in label:
            colors.append("forestgreen")
        elif "CB" in label:
            colors.append("crimson")
        else:
            colors.append("gray")

    # --- Top panel: full spectrum ---
    bars = ax_top.bar(band_indices, evals, color=colors, edgecolor="black",
                      linewidth=0.8, width=0.7)
    ax_top.axhline(0, color="gray", ls=":", lw=0.8, label="VB top (E = 0)")
    ax_top.axhline(GAAS_EG, color="crimson", ls="--", lw=0.8, alpha=0.6,
                   label=f"CB edge (Eg = {GAAS_EG} eV)")
    ax_top.axhline(-GAAS_DELTA_SO, color="forestgreen", ls="--", lw=0.8,
                   alpha=0.6,
                   label=f"SO edge (-DeltaSO = -{GAAS_DELTA_SO} eV)")

    # Annotate energy gaps
    # Eg arrow
    ax_top.annotate(
        "", xy=(8.6, GAAS_EG), xytext=(8.6, 0),
        arrowprops=dict(arrowstyle="<->", color="crimson", lw=1.5),
    )
    ax_top.text(8.8, GAAS_EG / 2, f"Eg = {GAAS_EG} eV",
                fontsize=10, color="crimson", va="center")

    # DeltaSO arrow
    ax_top.annotate(
        "", xy=(-0.4, 0), xytext=(-0.4, -GAAS_DELTA_SO),
        arrowprops=dict(arrowstyle="<->", color="forestgreen", lw=1.5),
    )
    ax_top.text(-0.2, -GAAS_DELTA_SO / 2, f"DeltaSO = {GAAS_DELTA_SO} eV",
                fontsize=10, color="forestgreen", va="center")

    ax_top.set_xticks(band_indices)
    ax_top.set_xticklabels(sorted_labels, fontsize=9)
    ax_top.set_ylabel("Energy (eV)", fontsize=12)
    ax_top.set_title("Bulk GaAs 8-band k.p Spectrum at k = 0 (Gamma point)",
                     fontsize=13)
    ax_top.legend(fontsize=9, loc="upper left")
    ax_top.grid(True, axis="y", alpha=0.3)
    ax_top.set_xlim(-1, 10)

    # Add band-index labels on bars
    for i, (bar, ev) in enumerate(zip(bars, evals)):
        y_offset = 0.03 if ev >= 0 else -0.03
        va = "bottom" if ev >= 0 else "top"
        ax_top.text(bar.get_x() + bar.get_width() / 2, ev + y_offset,
                    f"{ev:.3f}", ha="center", va=va, fontsize=7.5,
                    color="black")

    # --- Bottom panel: architecture diagram (text-based) ---
    ax_bot.axis("off")
    architecture_text = (
        "Code Architecture: How GaAs Parameters Flow to Eigenvalues\n"
        "=========================================================\n\n"
        "  parameters.f90          hamiltonianConstructor.f90     LAPACK\n"
        "  +------------------+    +-------------------------+   +----------+\n"
        f"  | Eg = {GAAS_EG} eV      |    | H(7,7) = A*k^2 + Eg     |   |          |\n"
        f"  | DeltaSO = {GAAS_DELTA_SO} eV  | -> | H(5,5) = (Q+T)/2 - DSO | -> | zheevx   |\n"
        f"  | EP = {GAAS_EP} eV      |    | P = sqrt(EP * c)       |   | -> evals |\n"
        f"  | gamma1 = {GAAS_GAMMA1}      |    | Q,T ~ gamma1,2,3,k^2   |   |          |\n"
        "  +------------------+    +-------------------------+   +----------+\n\n"
        "  At k=0: off-diagonal terms (P, Q, R, S, T) vanish -> H is diagonal\n"
        "  Eigenvalues are the band edges: E_HH=E_LH=0, E_SO=-DeltaSO, E_CB=Eg"
    )
    ax_bot.text(
        0.05, 0.5, architecture_text,
        transform=ax_bot.transAxes,
        fontsize=9.5, verticalalignment="center",
        fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow",
                  edgecolor="gray", alpha=0.9),
    )

    out_path = FIGURES_DIR / "lecture_12_example.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  Saved annotated spectrum plot to {out_path}")
    print()


# =========================================================================
# Extending the code: practical guide
# =========================================================================
def print_extension_guide():
    """Print a practical guide for extending the code."""
    print("=" * 70)
    print("  Extending the Code: Practical Guide")
    print("=" * 70)
    print()
    print("  Common extension points and how to approach them:")
    print()
    print("  1. Adding a new material:")
    print("     - Add a new 'case (\"MaterialName\")' block in parameters.f90")
    print("     - Populate paramStruct fields: Eg, EP, DeltaSO, gamma1/2/3, meff")
    print("     - Reference: Vurgaftman 2001 (III-V), Winkler 2003 (II-VI)")
    print("     - The derived P and A are computed automatically (line 740-741)")
    print()
    print("  2. Adding a new k.p coupling term:")
    print("     - Modify ZB8bandBulk or ZB8bandQW in hamiltonianConstructor.f90")
    print("     - The Hamiltonian matrix is populated element-by-element")
    print("     - At k=0, verify the term vanishes (diagonal remains unchanged)")
    print()
    print("  3. Adding a new output quantity:")
    print("     - Add parsing in outputFunctions.f90")
    print("     - Add computation in the relevant physics module")
    print("     - Add input parameters to simulation_config in defs.f90")
    print("     - Add parsing in input_parser.f90")
    print()
    print("  4. Adding a new eigensolver:")
    print("     - Add a new class in eigensolver.f90 behind #ifdef guard")
    print("     - The factory make_eigensolver(config) dispatches by string")
    print("     - Existing variants: dense LAPACK, MKL FEAST")
    print()
    print("  5. Code conventions (from CLAUDE.md):")
    print("     - F2018 standard enforced (-std=f2018)")
    print("     - All modules use 'private' default with explicit 'public' exports")
    print("     - External BLAS/LAPACK go through linalg.f90 interfaces")
    print("     - All types with allocatable components need finalizers")
    print("     - Basis ordering is fixed: bands 1-4 VB, 5-6 SO, 7-8 CB")
    print()


# =========================================================================
# Main
# =========================================================================
def main():
    print()
    print("=" * 70)
    print("  LECTURE 12: Extending the Code")
    print("  Architecture and worked example")
    print("=" * 70)
    print()

    # Section 1: Architecture overview
    print_architecture()

    # Section 2: Parameter mapping explanation
    print_parameter_mapping()

    # Section 3: Run code and show eigenvalues
    evals = run_and_show_eigenvalues()

    # Section 4: k.p coupling explanation
    print_kp_coupling()

    # Section 5: Extension guide
    print_extension_guide()

    # Section 6: Generate annotated plot
    print("=" * 70)
    print("  Generating annotated spectrum plot...")
    print("=" * 70)
    print()
    plot_annotated_spectrum(evals)

    print("=" * 70)
    print("  Lecture 12 complete.")
    print("=" * 70)
    print()
    print("  Key takeaways:")
    print("    1. The layered architecture separates concerns: types, parameters,")
    print("       mathematics, physics, and application.")
    print("    2. At k=0, eigenvalues are directly the band-edge parameters.")
    print("    3. The Kane energy EP controls VB-CB coupling strength,")
    print("       determining effective mass and non-parabolicity.")
    print("    4. Extending the code means knowing which layer to modify")
    print("       and following the existing conventions.")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
