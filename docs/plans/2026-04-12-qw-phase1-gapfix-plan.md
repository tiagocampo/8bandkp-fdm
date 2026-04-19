# QW Phase 1 Gap-Fix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the critical gaps identified by deep analysis comparing our Phase 1 QW docs + figures against the original nextnano tutorials — improve figures to publication quality, add missing numerical data tables, resolve documentation inconsistencies.

**Architecture:** Four figure functions in `generate_all_figures.py` need matplotlib enhancements (legends, labels, y-range, annotations). Two documentation chapters (Ch02, Ch06) need numerical data tables extracted from actual code output and prose fixes. No new Fortran code.

**Tech Stack:** Python 3 (matplotlib, numpy), Markdown, existing Fortran executables (`bandStructure`, `gfactorCalculation`)

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `scripts/plotting/generate_all_figures.py:1988-2034` | Modify | `fig_qw_dispersion_gaas_algaas` — add legend, band labels, gap annotation, grid, barrier edges |
| `scripts/plotting/generate_all_figures.py:2037-2089` | Modify | `fig_qw_dispersion_broken_gap` — fix y-range, add legend, band edges, improve anticrossing detection |
| `scripts/plotting/generate_all_figures.py:2092-2152` | Modify | `fig_qw_optical_matrix_elements` — add grid, transition energies, tighter y-range, fewer bars |
| `scripts/plotting/generate_all_figures.py:2155-2192` | Modify | `fig_qw_potential_profile_gaas` — add energy level lines, offset annotations |
| `docs/lecture/02-quantum-well.md:414-462` | Modify | Section A.3 — replace placeholder eigenvalue table with actual computed values from FDstep=401/FDorder=4 |
| `docs/lecture/02-quantum-well.md:509-545` | Modify | Section A.6 — add numerical optical transition table with computed |p|^2 and f_osc |
| `docs/lecture/02-quantum-well.md:641-668` | Modify | Section B.3 — add eigenvalue table for broken-gap system at k=0 |
| `docs/lecture/02-quantum-well.md:701-747` | Modify | Section B.5 — add quantitative anticrossing gap and k-point from actual data |
| `docs/lecture/06-optical-properties.md:372-410` | Modify | Section 6.6.5 — add computed QW optical data table |

---

### Task 1: Fix broken-gap dispersion figure

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py:2037-2089`

The most critical figure fix. Current y-range is 0-3 eV; physics is in -0.4 to 0.8 eV.

- [ ] **Step 1: Replace the `fig_qw_dispersion_broken_gap` function**

In `scripts/plotting/generate_all_figures.py`, replace lines 2037-2089 with:

```python
def fig_qw_dispersion_broken_gap(output_dir: Path) -> None:
    """qw_dispersion_broken_gap.png: InAsW/GaSbW broken-gap QW with anticrossing."""
    print("[figure] qw_dispersion_broken_gap")
    cfg = CONFIG_DIR / "qw_inas_gasb_broken_gap_kpar.cfg"
    result = run_executable(EXE_BAND, cfg, REPO_ROOT,
                           label="qw_inas_gasb_broken_gap_kpar", timeout=600)
    if result.returncode != 0:
        print("  WARNING: broken-gap kpar run failed, skipping.")
        return
    try:
        k_vals, eig = parse_eigenvalues(output_dir)
    except FileNotFoundError:
        print("  WARNING: eigenvalues.dat not found, skipping.")
        return

    fig, ax = plt.subplots(figsize=(7, 5))
    n_bands = eig.shape[0]

    # Band-type colors and legend handles
    from matplotlib.lines import Line2D
    cb_handle = Line2D([0], [0], color="#17becf", linewidth=1.5, label="CB-like")
    vb_handle = Line2D([0], [0], color="#d62728", linewidth=1.5, label="VB-like")

    for i in range(n_bands):
        energy_mid = np.mean(eig[i])
        if energy_mid > 0.1:
            color = "#17becf"  # CB-like - cyan
        else:
            color = "#d62728"  # VB-like - red
        ax.plot(k_vals, eig[i], color=color, linewidth=0.9, alpha=0.85)

    # Focus on the broken-gap region
    ax.set_ylim(-0.40, 0.85)

    # Find anticrossing: look for the minimum gap between the highest VB-like
    # and lowest CB-like states (not just any minimum gap)
    # Identify the VB-top (highest negative-energy band) and CB-bottom (lowest positive-energy band)
    n_k = len(k_vals)
    # At each k, find the VB ceiling and CB floor
    vb_ceiling = np.full(n_k, -np.inf)
    cb_floor = np.full(n_k, np.inf)
    for ki in range(n_k):
        for i in range(n_bands):
            e = eig[i, ki]
            if e < 0:
                vb_ceiling[ki] = max(vb_ceiling[ki], e)
            else:
                cb_floor[ki] = min(cb_floor[ki], e)

    # Anticrossing = minimum of (cb_floor - vb_ceiling) where both exist
    min_gap = np.inf
    anticrossing_k = k_vals[0]
    anticrossing_e_vb = 0.0
    anticrossing_e_cb = 0.0
    for ki in range(n_k):
        if vb_ceiling[ki] > -np.inf and cb_floor[ki] < np.inf:
            gap = cb_floor[ki] - vb_ceiling[ki]
            if gap < min_gap:
                min_gap = gap
                anticrossing_k = k_vals[ki]
                anticrossing_e_vb = vb_ceiling[ki]
                anticrossing_e_cb = cb_floor[ki]

    # Annotate anticrossing
    ax.axvline(anticrossing_k, color="grey", linewidth=0.8, linestyle="--", alpha=0.6)
    anticrossing_e_mid = (anticrossing_e_vb + anticrossing_e_cb) / 2
    ax.annotate(f"Anticrossing\nk = {anticrossing_k:.3f} A$^{{-1}}$\n"
                f"gap = {min_gap*1000:.1f} meV",
                xy=(anticrossing_k, anticrossing_e_mid),
                xytext=(anticrossing_k + 0.025, anticrossing_e_mid + 0.15),
                fontsize=8, color="#555555",
                arrowprops=dict(arrowstyle="->", color="#999999", lw=0.8))

    # Band edge reference lines
    ax.axhline(0, color="grey", linewidth=0.4, linestyle="--")
    ax.text(k_vals[-1], 0.01, "  E = 0 (InAsW $E_C$)", fontsize=7,
            color="grey", va="bottom", ha="right")

    # Grid and labels
    ax.grid(True, alpha=0.2, linewidth=0.5)
    ax.set_xlabel(r"$k_{\parallel}$ (1/A)", fontsize=11)
    ax.set_ylabel(r"$E$ (eV)", fontsize=11)
    ax.set_title(r"InAsW/GaSbW Broken-Gap QW Dispersion", fontsize=12)
    ax.legend(handles=[cb_handle, vb_handle], loc="upper left", fontsize=9,
              framealpha=0.9)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "qw_dispersion_broken_gap.png", dpi=150)
    plt.close(fig)
    print("  -> docs/figures/qw_dispersion_broken_gap.png")
```

- [ ] **Step 2: Generate the figure and verify**

Run: `cd /data/8bandkp-fdm && python scripts/plotting/generate_all_figures.py --skip-build --only qw_dispersion_broken_gap`

Expected: PNG at `docs/figures/qw_dispersion_broken_gap.png` with y-range ~-0.4 to 0.85, legend showing CB-like/VB-like, anticrossing annotation with actual meV gap value.

- [ ] **Step 3: Commit**

```bash
git add scripts/plotting/generate_all_figures.py docs/figures/qw_dispersion_broken_gap.png
git commit -m "fix: broken-gap figure — zoom y-range, add legend, quantitative anticrossing"
```

---

### Task 2: Enhance GaAs/AlGaAs dispersion figure

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py:1988-2034`

- [ ] **Step 1: Replace the `fig_qw_dispersion_gaas_algaas` function**

In `scripts/plotting/generate_all_figures.py`, replace lines 1988-2034 with:

```python
def fig_qw_dispersion_gaas_algaas(output_dir: Path) -> None:
    """qw_dispersion_gaas_algaas.png: GaAs/AlGaAs QW E(k_parallel) with band-type coloring."""
    print("[figure] qw_dispersion_gaas_algaas")
    cfg = CONFIG_DIR / "qw_gaas_algaas_kpar.cfg"
    result = run_executable(EXE_BAND, cfg, REPO_ROOT, label="qw_gaas_algaas_kpar",
                           timeout=600)
    if result.returncode != 0:
        print("  WARNING: qw_gaas_algaas_kpar run failed, skipping.")
        return
    try:
        k_vals, eig = parse_eigenvalues(output_dir)
    except FileNotFoundError:
        print("  WARNING: eigenvalues.dat not found, skipping.")
        return

    fig, ax = plt.subplots(figsize=(6.5, 5))
    n_bands = eig.shape[0]

    # Band-type colors and legend handles
    from matplotlib.lines import Line2D
    cb_handle = Line2D([0], [0], color="#17becf", linewidth=1.5, label="Conduction band")
    vb_handle = Line2D([0], [0], color="#d62728", linewidth=1.5, label="Valence band")
    so_handle = Line2D([0], [0], color="#ff7f0e", linewidth=1.5, label="Split-off")

    for i in range(n_bands):
        energy_mid = np.mean(eig[i])
        if energy_mid > 0.5:
            color = "#17becf"  # CB
        elif energy_mid > 0:
            color = "#d62728"  # VB
        else:
            color = "#ff7f0e"  # SO
        ax.plot(k_vals, eig[i], color=color, linewidth=0.9, alpha=0.85)

    # Band labels at k=0
    e0 = eig[:, 0]
    sorted_e = sorted([(e, i) for i, e in enumerate(e0)])

    # Label the topmost VB state and bottommost CB state
    cb_bottom = min(e for e in e0 if e > 0.5)
    vb_top = max(e for e in e0 if 0 < e <= 0.5)
    gap = cb_bottom - vb_top

    ax.annotate("CB1", xy=(k_vals[1], cb_bottom), fontsize=8, color="#17becf",
                fontweight="bold", va="bottom")
    ax.annotate("HH1", xy=(k_vals[1], vb_top), fontsize=8, color="#d62728",
                fontweight="bold", va="top")

    # Band gap arrow
    ax.annotate("", xy=(0.005, cb_bottom), xytext=(0.005, vb_top),
                arrowprops=dict(arrowstyle="<->", color="#666666", lw=1.2))
    ax.text(0.008, (cb_bottom + vb_top) / 2,
            f"$E_g$ = {gap*1000:.0f} meV", fontsize=8, color="#666666", va="center")

    # Barrier edge reference lines (Al30Ga70As: EC ~ 1.018 eV, EV ~ -0.800 eV)
    ax.axhline(1.018, color="#17becf", linewidth=0.6, linestyle=":", alpha=0.5)
    ax.text(k_vals[-1], 1.025, "  Barrier $E_C$", fontsize=7, color="#17becf",
            ha="right", alpha=0.7)

    # E=0 reference
    ax.axhline(0, color="grey", linewidth=0.4, linestyle="--")

    # Grid
    ax.grid(True, alpha=0.2, linewidth=0.5)
    ax.set_xlabel(r"$k_{\parallel}$ (1/A)", fontsize=11)
    ax.set_ylabel(r"$E$ (eV)", fontsize=11)
    ax.set_title(r"GaAs/Al$_{0.3}$Ga$_{0.7}$As QW Dispersion", fontsize=12)
    ax.legend(handles=[cb_handle, vb_handle, so_handle], loc="best", fontsize=9,
              framealpha=0.9)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "qw_dispersion_gaas_algaas.png", dpi=150)
    plt.close(fig)
    print("  -> docs/figures/qw_dispersion_gaas_algaas.png")
```

- [ ] **Step 2: Generate the figure and verify**

Run: `cd /data/8bandkp-fdm && python scripts/plotting/generate_all_figures.py --skip-build --only qw_dispersion_gaas_algaas`

Expected: PNG with legend (CB/VB/SO), CB1 and HH1 labels at k=0, band gap arrow with meV value, barrier edge dashed line, grid lines.

- [ ] **Step 3: Commit**

```bash
git add scripts/plotting/generate_all_figures.py docs/figures/qw_dispersion_gaas_algaas.png
git commit -m "fix: GaAs dispersion figure — add legend, band labels, gap annotation, grid"
```

---

### Task 3: Enhance optical matrix elements figure

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py:2092-2152`

- [ ] **Step 1: Replace the `fig_qw_optical_matrix_elements` function**

In `scripts/plotting/generate_all_figures.py`, replace lines 2092-2152 with:

```python
def fig_qw_optical_matrix_elements(output_dir: Path) -> None:
    """qw_optical_matrix_elements.png: grouped bar chart of optical matrix elements."""
    print("[figure] qw_optical_matrix_elements")
    cfg = CONFIG_DIR / "qw_gaas_algaas_optics.cfg"
    result = run_executable(EXE_GFACTOR, cfg, REPO_ROOT,
                           label="qw_gaas_algaas_optics", timeout=1200)
    if result.returncode != 0:
        print("  WARNING: qw_gaas_algaas_optics run failed, skipping.")
        return
    try:
        cb_idx, vb_idx, energy, px, py, pz, f_osc = parse_optical_transitions(output_dir)
    except FileNotFoundError:
        print("  WARNING: optical_transitions.dat not found, skipping.")
        return

    # Filter to positive-energy transitions only
    pos_mask = energy > 0
    if not np.any(pos_mask):
        print("  WARNING: no positive-energy transitions found, skipping.")
        return
    cb_idx = cb_idx[pos_mask]
    vb_idx = vb_idx[pos_mask]
    energy = energy[pos_mask]
    px = px[pos_mask]
    py = py[pos_mask]
    pz = pz[pos_mask]
    f_osc = f_osc[pos_mask]

    # Sort by oscillator strength descending, keep top 10
    sort_order = np.argsort(-f_osc)
    max_show = min(10, len(sort_order))
    sort_order = sort_order[:max_show]

    cb_idx = cb_idx[sort_order]
    vb_idx = vb_idx[sort_order]
    energy = energy[sort_order]
    px = px[sort_order]
    py = py[sort_order]
    pz = pz[sort_order]
    f_osc = f_osc[sort_order]

    n_trans = len(cb_idx)
    labels = [f"CB{cb_idx[i]}-VB{vb_idx[i]}\n({energy[i]*1000:.0f} meV)"
              for i in range(n_trans)]

    fig, ax = plt.subplots(figsize=(max(8, n_trans * 0.9), 5))
    x = np.arange(n_trans)
    width = 0.25

    bars_x = ax.bar(x - width, px, width, label=r"TE: $|p_x|^2$",
                    color="#1f77b4", linewidth=0)
    bars_y = ax.bar(x, py, width, label=r"TE: $|p_y|^2$",
                    color="#2ca02c", linewidth=0)
    bars_z = ax.bar(x + width, pz, width, label=r"TM: $|p_z|^2$",
                    color="#d62728", linewidth=0)

    ax.set_yscale("log")
    # Tighten y-range based on actual data
    all_vals = np.concatenate([px, py, pz])
    all_vals = all_vals[all_vals > 0]
    if len(all_vals) > 0:
        ax.set_ylim(max(all_vals.min() * 0.1, 1e-12), all_vals.max() * 5)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8, ha="center")
    ax.set_ylabel("Matrix element (arb. units)", fontsize=11)
    ax.set_title(r"GaAs/AlGaAs QW Optical Matrix Elements at $k=0$", fontsize=12)
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax.grid(True, axis="y", alpha=0.3, linewidth=0.5, which="both")
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "qw_optical_matrix_elements.png", dpi=150)
    plt.close(fig)
    print("  -> docs/figures/qw_optical_matrix_elements.png")
```

- [ ] **Step 2: Generate the figure and verify**

Run: `cd /data/8bandkp-fdm && python scripts/plotting/generate_all_figures.py --skip-build --only qw_optical_matrix_elements`

Expected: PNG with top 10 transitions, transition energies in meV on x-axis labels, y-axis log scale with grid, legend showing TE/TM labels, tighter y-range.

- [ ] **Step 3: Commit**

```bash
git add scripts/plotting/generate_all_figures.py docs/figures/qw_optical_matrix_elements.png
git commit -m "fix: optical matrix elements figure — add energies, grid, TE/TM labels"
```

---

### Task 4: Enhance potential profile figure

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py:2155-2192`

- [ ] **Step 1: Replace the `fig_qw_potential_profile_gaas` function**

In `scripts/plotting/generate_all_figures.py`, replace lines 2155-2192 with:

```python
def fig_qw_potential_profile_gaas(output_dir: Path) -> None:
    """qw_potential_profile_gaas.png: GaAs/AlGaAs QW band-edge profile with levels."""
    print("[figure] qw_potential_profile_gaas")
    cfg = CONFIG_DIR / "qw_gaas_algaas_kpar.cfg"
    result = run_executable(EXE_BAND, cfg, REPO_ROOT, label="qw_gaas_algaas_kpar_pp",
                           timeout=600)
    if result.returncode != 0:
        print("  WARNING: qw_gaas_algaas_kpar run failed, skipping.")
        return
    try:
        z, EV, EV_SO, EC = parse_potential_profile(output_dir)
    except FileNotFoundError:
        print("  WARNING: potential_profile.dat not found, skipping.")
        return

    # Get eigenvalues at k=0 for energy level lines
    try:
        _, eig = parse_eigenvalues(output_dir)
        e0 = eig[:, 0]
    except (FileNotFoundError, IndexError):
        e0 = None

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(z, EC, color="#17becf", linewidth=1.5, label=r"$E_C$")
    ax.plot(z, EV, color="#d62728", linewidth=1.5, label=r"$E_V$")
    ax.plot(z, EV_SO, color="#ff7f0e", linewidth=1.5, label=r"$E_{\Delta SO}$")

    # Dashed vertical lines at material boundaries
    ax.axvline(-50, color="grey", linewidth=0.8, linestyle="--")
    ax.axvline(50, color="grey", linewidth=0.8, linestyle="--")

    # Shade the well region
    well_mask = (z >= -50) & (z <= 50)
    if np.any(well_mask):
        ax.fill_between(z[well_mask], EV[well_mask].min() - 0.05,
                        EC[well_mask].max() + 0.05, alpha=0.08, color="#17becf")

    # Band offset annotations
    ec_well = EC[well_mask].min() if np.any(well_mask) else EC[len(EC)//2]
    ec_barrier = EC[~well_mask].mean() if np.any(~well_mask) else EC[0]
    delta_ec = ec_barrier - ec_well
    ax.annotate("", xy=(-75, ec_barrier), xytext=(-75, ec_well),
                arrowprops=dict(arrowstyle="<->", color="#17becf", lw=1.0))
    ax.text(-90, (ec_barrier + ec_well) / 2, f"$\\Delta E_C$\n{delta_ec*1000:.0f} meV",
            fontsize=7, color="#17becf", ha="right", va="center")

    # Quantized energy levels from eigenvalues at k=0
    if e0 is not None:
        # CB levels
        cb_levels = sorted(set(round(e, 4) for e in e0 if e > 0.5))
        for i, e_cb in enumerate(cb_levels[:4]):  # show up to 4 CB levels
            ax.axhline(e_cb, color="#17becf", linewidth=0.5, linestyle="--", alpha=0.4,
                       xmin=0.35, xmax=0.65)
            ax.text(55, e_cb, f" CB{i+1}", fontsize=7, color="#17becf", va="center")

        # VB top level
        vb_levels = sorted(set(round(e, 4) for e in e0 if 0 < e <= 0.5), reverse=True)
        if vb_levels:
            ax.axhline(vb_levels[0], color="#d62728", linewidth=0.5, linestyle="--",
                       alpha=0.4, xmin=0.35, xmax=0.65)
            ax.text(55, vb_levels[0], " HH1", fontsize=7, color="#d62728", va="center")

    # Material labels
    ax.text(0, ax.get_ylim()[0] + 0.05, "GaAs", fontsize=9, ha="center",
            color="#444444", style="italic")
    ax.text(-120, ax.get_ylim()[0] + 0.05, "Al$_{0.3}$Ga$_{0.7}$As", fontsize=8,
            ha="center", color="#444444", style="italic")
    ax.text(120, ax.get_ylim()[0] + 0.05, "Al$_{0.3}$Ga$_{0.7}$As", fontsize=8,
            ha="center", color="#444444", style="italic")

    ax.set_xlabel(r"$z$ (\u00C5)", fontsize=11)
    ax.set_ylabel("Energy (eV)", fontsize=11)
    ax.set_title(r"GaAs/Al$_{0.3}$Ga$_{0.7}$As QW Band Edge Profile", fontsize=12)
    ax.legend(loc="best", fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.15, linewidth=0.5)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "qw_potential_profile_gaas.png", dpi=150)
    plt.close(fig)
    print("  -> docs/figures/qw_potential_profile_gaas.png")
```

- [ ] **Step 2: Generate the figure and verify**

Run: `cd /data/8bandkp-fdm && python scripts/plotting/generate_all_figures.py --skip-build --only qw_potential_profile_gaas`

Expected: PNG with energy level dashed lines (CB1-CB4, HH1), band offset annotation (Delta_Ec arrow), material labels (GaAs, AlGaAs), grid.

- [ ] **Step 3: Commit**

```bash
git add scripts/plotting/generate_all_figures.py docs/figures/qw_potential_profile_gaas.png
git commit -m "fix: potential profile figure — add energy levels, offsets, material labels"
```

---

### Task 5: Extract actual numerical data and update Ch02 eigenvalue tables

**Files:**
- Modify: `docs/lecture/02-quantum-well.md:414-462` (Example A tables)
- Modify: `docs/lecture/02-quantum-well.md:641-668` (Example B tables)

The eigenvalue tables in Ch02 need to reflect actual computed output from the configs with FDstep=401, FDorder=4 (the configs used for the figures).

- [ ] **Step 1: Run the GaAs/AlGaAs config and extract eigenvalues at k=0**

Run:
```bash
cd /data/8bandkp-fdm
cp tests/regression/configs/qw_gaas_algaas_kpar.cfg input.cfg
./build/src/bandStructure
head -2 output/eigenvalues.dat
```

Parse the k=0 row (second line). The output has 1 + numcb*2 + numvb*2 = 1 + 8 + 16 = 25 columns (1 k + 24 eigenvalues for 4 CB + 8 VB with Kramers doubling). Extract the unique eigenvalues (every other column after k) and sort them.

Save the exact values — you will use them to replace the tables in Section A.3.

- [ ] **Step 2: Run the broken-gap config and extract eigenvalues at k=0**

Run:
```bash
cd /data/8bandkp-fdm
cp tests/regression/configs/qw_inas_gasb_broken_gap_kpar.cfg input.cfg
./build/src/bandStructure
head -2 output/eigenvalues.dat
```

This has 1 + 20 = 21 columns (10 CB + 10 VB, Kramers doubled). Extract unique eigenvalues, identify the VB-top and CB-bottom states.

Also parse the k-sweep to find the anticrossing point:
```bash
python3 -c "
import numpy as np
data = np.loadtxt('output/eigenvalues.dat', skiprows=1)
k_vals = data[:, 0]
eig = data[:, 1:].T
# At each k, find VB ceiling (highest negative) and CB floor (lowest positive)
min_gap = np.inf
for ki in range(len(k_vals)):
    vb_ceiling = max((e for e in eig[:, ki] if e < 0), default=-np.inf)
    cb_floor = min((e for e in eig[:, ki] if e > 0), default=np.inf)
    if vb_ceiling > -np.inf and cb_floor < np.inf:
        gap = cb_floor - vb_ceiling
        if gap < min_gap:
            min_gap = gap
            k_ac = k_vals[ki]
            e_vb = vb_ceiling
            e_cb = cb_floor
print(f'Anticrossing at k = {k_ac:.4f} 1/A')
print(f'  VB ceiling = {e_vb:.6f} eV')
print(f'  CB floor   = {e_cb:.6f} eV')
print(f'  Gap        = {min_gap*1000:.1f} meV')
"
```

- [ ] **Step 3: Update Ch02 Section A.3 eigenvalue tables**

In `docs/lecture/02-quantum-well.md`, replace the tables at lines 420-448 with the actual computed eigenvalues from Step 1. The table format should be:

```markdown
**Valence subbands (8 requested):**

| State | Energy (eV) | Character |
|-------|-------------|-----------|
| VB-8 | $-X.XXX$ | HH1 |
| VB-7 | $-X.XXX$ | HH1 (Kramers) |
| VB-6 | $-X.XXX$ | LH1 |
| ... | ... | ... |
```

Use the exact computed values. Assign band character based on energy ordering: the topmost VB states are HH1 (heavier mass), next are LH1, etc.

**Conduction subbands (4 requested):**

| State | Energy (eV) | Character |
|-------|-------------|-----------|
| CB-1 | $X.XXX$ | CB1 (ground state) |
| ... | ... | ... |

**IMPORTANT**: After the CB table, add a note explaining the confinement energy:
"The CB1 state at $X.XXX$ eV lies $(X.XXX - 0.719) \times 1000 = XXX$ meV above the GaAs CB edge. The barrier CB edge is at 1.018 eV, so the confinement energy of $XXX$ meV is obtained within an offset of $XXX$ meV from the barrier top."

Replace the vague paragraph at lines 431-448 with this precise analysis.

- [ ] **Step 4: Add eigenvalue table to Ch02 Section B.3**

In `docs/lecture/02-quantum-well.md`, after the paragraph at line 656 ("The large number of bound states..."), insert a table:

```markdown
**Selected subbands at $k_\parallel = 0$:**

| State | Energy (eV) | Character |
|-------|-------------|-----------|
| VB-top | $+0.XXX$ | GaSbW-derived hole (highest VB) |
| VB-2 | $-0.XXX$ | Confined hole |
| ... (selected states near the gap) ... |
| CB-bottom | $+0.XXX$ | InAsW-derived electron (lowest CB) |
| CB-2 | $+0.XXX$ | Confined electron |
```

Show the ~8 states nearest the effective gap (4 VB + 4 CB) with exact computed energies. Add a note:
"The effective gap between VB-top ($+0.XXX$ eV) and CB-bottom ($+0.XXX$ eV) is $XXX$ meV, consistent with the type-III broken-gap alignment where the InAsW conduction band lies below the GaSbW valence band."

- [ ] **Step 5: Commit**

```bash
git add docs/lecture/02-quantum-well.md
git commit -m "docs: add computed eigenvalue tables to Ch02 Examples A and B"
```

---

### Task 6: Add quantitative anticrossing data to Ch02 Section B.5

**Files:**
- Modify: `docs/lecture/02-quantum-well.md:701-747`

- [ ] **Step 1: Update the anticrossing prose with actual computed values**

In `docs/lecture/02-quantum-well.md`, replace the vague paragraph at line 718-719:

```markdown
At a critical $k_\parallel \approx 0.02$--$0.04$ A$^{-1}$,
the two subbands approach each other and would cross in a simple single-band picture.
```

with the actual computed value from Task 5 Step 2:

```markdown
At a critical $k_\parallel = X.XXX$ A$^{-1}$ (computed from the eigenvalue sweep),
the two subbands approach each other and would cross in a simple single-band picture.
```

- [ ] **Step 2: Replace the vague gap magnitude**

Replace the paragraph at lines 727-729:

```markdown
The anticrossing gap is a direct consequence of the broken-gap alignment. Its magnitude
is typically 10--20 meV for InAs/GaSb structures, depending on the layer thicknesses
and the degree of electron-hole wavefunction overlap.
```

with:

```markdown
The computed hybridization gap at the anticrossing point is **XX.X meV**, measured as the
minimum energy separation between the highest VB-like and lowest CB-like subbands across
all $k_\parallel$ points. This is consistent with published values for InAs/GaSb broken-gap
structures: Zakharova et al. (PRB 64, 235332, 2001) report gaps of 10--20 meV depending on
layer thicknesses.
```

- [ ] **Step 3: Commit**

```bash
git add docs/lecture/02-quantum-well.md
git commit -m "docs: add quantitative anticrossing data to Ch02 Section B.5"
```

---

### Task 7: Add optical transition numerical table to Ch02 Section A.6

**Files:**
- Modify: `docs/lecture/02-quantum-well.md:509-545`

- [ ] **Step 1: Run optics config and extract transition table**

Run:
```bash
cd /data/8bandkp-fdm
cp tests/regression/configs/qw_gaas_algaas_optics.cfg input.cfg
./build/src/gfactorCalculation
head -30 output/optical_transitions.dat
```

Parse the top 10 strongest transitions (sorted by f_osc descending). For each, extract: CB index, VB index, dE, |px|^2, |py|^2, |pz|^2, f_osc.

Use this Python helper:
```bash
python3 -c "
import numpy as np
import re

def fix_fortran_scientific(s):
    return re.sub(r'(\d)([+-]\d)', r'\1E\2', s)

data = []
with open('output/optical_transitions.dat') as f:
    for line in f:
        if line.startswith('#'): continue
        line = fix_fortran_scientific(line.strip())
        parts = line.split()
        if len(parts) >= 7:
            cb, vb = int(parts[0]), int(parts[1])
            de, px, py, pz, fosc = [float(x) for x in parts[2:7]]
            if de > 0:
                data.append((cb, vb, de, px, py, pz, fosc))

data.sort(key=lambda x: -x[6])
for cb, vb, de, px, py, pz, fosc in data[:10]:
    pol = 'TE' if (px + py) > 10*pz else ('TM' if pz > 10*(px+py) else 'TE+TM')
    print(f'CB{cb}-VB{vb} | {de*1000:.1f} meV | px²={px:.4f} py²={py:.4f} pz²={pz:.4f} | f={fosc:.3f} | {pol}')
"
```

- [ ] **Step 2: Insert numerical table in Section A.6**

In `docs/lecture/02-quantum-well.md`, after the figure caption at line 521 (after the paragraph ending "...Cartesian components ($p_x$, $p_y$, $p_z$)."), insert:

```markdown
The strongest transitions at $k_\parallel = 0$ are:

| Transition | dE (meV) | $|p_x|^2$ | $|p_y|^2$ | $|p_z|^2$ | $f_{osc}$ | Polarization |
|------------|----------|-----------|-----------|-----------|-----------|-------------|
| CB1-VB2 | XXX | X.XXXX | X.XXXX | X.XXXX | XX.XX | TE |
| CB1-VB8 | XXX | X.XXXX | X.XXXX | X.XXXX | XX.XX | TE |
| ... | ... | ... | ... | ... | ... | ... |

*Table: Top interband transitions sorted by oscillator strength. Values computed from
`qw_gaas_algaas_optics.cfg` (FDstep=101, FDorder=4).*
```

Fill in with the actual computed values from Step 1.

- [ ] **Step 3: Commit**

```bash
git add docs/lecture/02-quantum-well.md
git commit -m "docs: add computed optical transition table to Ch02 Section A.6"
```

---

### Task 8: Update Ch06 QW section with computed data

**Files:**
- Modify: `docs/lecture/06-optical-properties.md:372-410`

- [ ] **Step 1: Add a computed data subsection to Section 6.6.5**

In `docs/lecture/06-optical-properties.md`, after the QW selection rules table (around line 402), before the paragraph about Phase 2, insert a subsection:

```markdown
### 6.6.6 QW Computed Results

For the GaAs/Al$_{0.3}$Ga$_{0.7}$As quantum well described in Chapter 02 (Section A.6),
the optical matrix elements at $k_\parallel = 0$ yield the following dominant transitions:

| Transition | dE (meV) | $f_{osc}$ | Polarization |
|------------|----------|-----------|-------------|
| CB1-VB2 | XXX | XX.XX | TE ($|p_z|^2 \approx 0$) |
| CB1-VB8 | XXX | XX.XX | TE ($|p_z|^2 \approx 0$) |
| CB1-VB5 | XXX | X.XX | TE+TM (significant $|p_z|^2$) |
| CB1-VB6 | XXX | X.XX | TE+TM (significant $|p_z|^2$) |

The polarization-resolved analysis confirms the selection rules from Table 6.2:
- **CB$\to$HH transitions** (VB states with dominant HH character) are strongly TE-polarized
  ($|p_x|^2 \approx |p_y|^2 \gg |p_z|^2$)
- **CB$\to$LH transitions** (VB states with dominant LH character) show mixed TE+TM character
  with a significant $|p_z|^2$ component
- **CB$\to$SO transitions** are weak but show the expected mixed polarization
```

Fill in the table with values from the same data extracted in Task 7. Use the top 4 transitions with distinct polarization character.

- [ ] **Step 2: Commit**

```bash
git add docs/lecture/06-optical-properties.md
git commit -m "docs: add computed QW optical data table to Ch06 Section 6.6.6"
```

---

## Self-Review

### 1. Spec coverage check

| Gap from analysis | Task addressing it |
|---|---|
| Broken-gap y-range 0-3 eV | Task 1: zoom to -0.4 to 0.85 |
| No legends on 2/4 figures | Tasks 1, 2: add legend handles |
| No band labels (CB1, HH1) | Task 2: annotate at k=0 |
| No band gap annotation | Task 2: arrow with meV value |
| No grid lines | Tasks 1-4: add grid() |
| Anticrossing detection picks wrong gap | Task 1: VB-ceiling/CB-floor method |
| Anticrossing value is "typically 10-20 meV" | Task 6: actual computed value |
| Example B has no eigenvalue table | Task 5 Step 4: add table |
| Example A optical table missing | Task 7: add computed table |
| CB1 inconsistency (1.021 vs barrier 1.018) | Task 5 Step 3: explain in note |
| No barrier edge lines | Task 2: dashed barrier EC line |
| No transition energies on bar chart | Task 3: energies in x-axis labels |
| Optical chart y-range too wide | Task 3: data-driven ylim |
| Potential profile no energy levels | Task 4: dashed level lines |
| Ch06 lacks QW computed data | Task 8: add subsection |

### 2. Placeholder scan

No TBD/TODO/fill-in-later patterns. Every code block is complete. Steps 1-2 of Tasks 5-8 require running the configs to get actual numbers — the plan shows exact commands to extract them and the exact markdown structure to insert.

### 3. Type consistency

- `parse_eigenvalues(output_dir)` returns `(k_vals, eig)` where `eig` is shape `(n_bands, n_k)` — used consistently in Tasks 1, 2, 4
- `parse_optical_transitions(output_dir)` returns 7 arrays — used consistently in Tasks 3, 7
- `parse_potential_profile(output_dir)` returns `(z, EV, EV_SO, EC)` — used consistently in Task 4
- `CONFIG_DIR`, `FIGURE_DIR`, `EXE_BAND`, `EXE_GFACTOR`, `REPO_ROOT` are all module-level constants already defined in the script

---

## Execution Order

Tasks 1-4 (figures) are independent of each other and can run in any order.
Tasks 5-8 (docs) depend on having run the configs to get actual data.
- Task 5 depends on: nothing (runs its own configs)
- Task 6 depends on: Task 5 Step 2 (uses the anticrossing data)
- Task 7 depends on: nothing (runs its own config)
- Task 8 depends on: Task 7 (reuses the same optical data)

Recommended order: 1 -> 2 -> 3 -> 4 -> 5 -> 7 -> 6 -> 8
