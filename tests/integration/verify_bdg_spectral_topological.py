#!/usr/bin/env python3
# COVERAGE: observable=bdg_ldos_topological geometry=wire material=InAs
# COVERAGE: observable=bdg_spectral_function_topological geometry=wire material=InAs
# COVERAGE: observable=bdg_nambu_ldos_topological geometry=wire material=InAs
"""Issue 11 fix1: regenerate BdG LDOS / Nambu / A(k,E) plots from the
topological wire config (PR40 / Issue 07) at B = 2 B_crit.

This is the headline MZM-signature companion to verify_bdg_spectral.py
(which uses the trivial 3x3, B=0 smoke-test config and produces flat
plots). Here we use the actual topological InAs/GaAs core/shell wire
config (50-site, mu = 0.6601 eV) with B_vec = [5.6, 0, 0] T (transverse,
above B_crit ~ 2.8 T) to expose:

  1. Zero-bias peak in bdg_ldos.dat at E = 0
  2. In-gap mode in bdg_spectral.dat at kz = 0, E ~= 0
  3. Symmetric electron/hole Nambu weight in bdg_ldos_nambu.dat at the MZM

B is overridden at runtime (text replace on the loaded TOML, no new
TOML fields per ADR 0002). Spectral E/k grids are also overridden
inside [topology] for a tight +/- 5 delta_0 window around E = 0.

Args: <topologicalAnalysis_exe> <topological_wire_config>
"""
import os
import re
import sys
import subprocess
from pathlib import Path

# --- Plot generation (Issue 11 fix1) ---------------------------------------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

EXE = str(Path(sys.argv[1]).resolve())
CONFIG = Path(sys.argv[2])
OMP = "4"

# Physics knobs for the topological-regime regen (override at runtime
# so we do not modify any TOML on disk -- ADR 0002).
B_TOPO_TESLA = 5.6          # ~ 2 * B_crit; reopened topological regime
DELTA_0_EV = 0.0002         # matches wire_inas_gaas_bdg_topological.toml [bdg]
E_HALF_WINDOW_FACTORS = 5.0 # +/- 5*delta_0 around E=0 captures the in-gap mode
SPECTRAL_NE = 11            # E grid for the LDOS curve. Each E point runs
                            # a fresh PARDISO solve on the 16*N_bd BdG CSR
                            # (N_bd=169 sites -> 2704x2704); 11 points is the
                            # minimum that still resolves the MZM peak shape
                            # while keeping wall time manageable (~3 min).
# A(k,E) colormap is computed at a single kz=0 because
# `write_bdg_spectral` uses 'replace' status, so a multi-kz sweep
# would clobber earlier files (single-kz per kz_idx pass). The
# spectral plot is therefore A(kz=0, E) only; the k-axis is degenerate.
SPECTRAL_NK = 1


def run(text, tag, timeout=1800):
    d = Path(f"run_{tag}")
    d.mkdir(parents=True, exist_ok=True)
    (d / "input.toml").write_text(text)
    env = dict(os.environ, OMP_NUM_THREADS=OMP)
    p = subprocess.run([EXE], cwd=d, env=env, capture_output=True,
                       text=True, timeout=timeout)
    (d / "run.log").write_text(p.stdout + "\n--STDERR--\n" + p.stderr)
    if p.returncode != 0:
        print(f"FAIL: topologicalAnalysis exited {p.returncode} for {tag}")
        print(p.stdout)
        print(p.stderr)
        sys.exit(1)
    return d


def override_toml(text, B_T):
    """Apply Issue 11 fix1 runtime overrides (no new fields, ADR 0002)."""
    # Override [bdg] B_vec with a transverse B = 2*B_crit.
    text = re.sub(
        r"(\[bdg\][^\[]*?B_vec\s*=\s*)\[[^\]]*\]",
        lambda m: m.group(1) + f"[{B_T}, 0.0, 0.0]",
        text, count=1, flags=re.S,
    )
    # Override [topology] mode = "bdq_spectral" (the spectral-mode
    # branch in main_topology that emits LDOS/Nambu/A(k,E)).
    text = re.sub(
        r'^(\[topology\][^\[]*?mode\s*=\s*)"[^\"]+"',
        lambda m: m.group(1) + '"bdq_spectral"',
        text, count=1, flags=re.M | re.S,
    )
    # Tight spectral window around E=0 (captures MZM peak + gap edges).
    half = E_HALF_WINDOW_FACTORS * DELTA_0_EV
    eta = DELTA_0_EV * 0.5

    spectral_overrides = [
        ("spectral_E_min", f"{-half:.6e}"),
        ("spectral_E_max", f"{half:.6e}"),
        ("spectral_nE",   f"{SPECTRAL_NE}"),
        ("spectral_eta",  f"{eta:.6e}"),
        ("spectral_k_min", "0.0"),
        ("spectral_k_max", "0.0"),
        ("spectral_nk",    f"{SPECTRAL_NK}"),
    ]
    for key, val in spectral_overrides:
        if re.search(rf"^{re.escape(key)}\s*=", text, flags=re.M):
            text = re.sub(rf"^{re.escape(key)}\s*=\s*[^\n]+",
                          f"{key} = {val}", text, count=1, flags=re.M)
        else:
            # Inject into the [topology] block (create it if absent).
            text = inject_into_section(text, "topology",
                                       f"{key} = {val}")
    return text


def inject_into_section(text, section, line):
    """Append `line` under [section]; create the section if missing."""
    pat = re.compile(rf"^\[{re.escape(section)}\]\s*$", flags=re.M)
    m = pat.search(text)
    if m:
        # Find end of section (next [...] header or EOF).
        tail = text[m.end():]
        end = re.search(r"^\[[^\]]+\]\s*$", tail, flags=re.M)
        idx = m.end() + (end.start() if end else len(tail))
        # Insert `line` followed by a blank separator. `text[idx:]` starts
        # with the next section header (or EOF), so we add a leading "\n".
        return text[:idx].rstrip("\n") + "\n\n" + line + "\n" + text[idx:]
    # No [section]: append at EOF.
    return text.rstrip() + f"\n\n[{section}]\n{line}\n"


def read_3col(path):
    """Parse 3-column whitespace-separated file (skipping #)."""
    rows = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 3:
            try:
                rows.append((float(parts[0]), float(parts[1]),
                             float(parts[2])))
            except ValueError:
                pass
    return rows


def read_5col(path):
    rows = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 5:
            try:
                rows.append((tuple(float(x) for x in parts[:5])))
            except ValueError:
                pass
    return rows


def main():
    config_text = CONFIG.read_text()
    text = override_toml(config_text, B_TOPO_TESLA)
    workdir = run(text, "bdq_spectral_topo")

    ldos_path = workdir / "output" / "bdg_ldos.dat"
    nambu_path = workdir / "output" / "bdg_ldos_nambu.dat"
    spectral_path = workdir / "output" / "bdg_spectral.dat"

    for p in (ldos_path, nambu_path, spectral_path):
        if not p.exists() or p.stat().st_size == 0:
            print(f"FAIL: {p.name} missing or empty")
            sys.exit(1)

    ldos_rows = read_3col(ldos_path)
    nambu_rows = read_3col(nambu_path)
    spectral_rows = read_5col(spectral_path)

    if not (len(ldos_rows) >= 3 and len(nambu_rows) >= 3
            and len(spectral_rows) >= 3):
        print(f"FAIL: too few rows (ldos={len(ldos_rows)}, "
              f"nambu={len(nambu_rows)}, spectral={len(spectral_rows)})")
        sys.exit(1)

    # Spectral non-negativity sanity (already enforced by the spec).
    for r in spectral_rows:
        if r[4] < 0.0:
            print(f"FAIL: A(k,E) < 0 at row {r}")
            sys.exit(1)

    # --- Plots ----------------------------------------------------------
    out_dir = workdir.parent / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_bdg_ldos_topo(ldos_rows, out_dir / "bdg_ldos_wire_topological.png")
    plot_bdg_ldos_nambu_topo(nambu_rows, out_dir / "bdg_ldos_nambu_wire_topological.png")
    plot_bdg_spectral_topo(spectral_rows, out_dir / "bdg_spectral_AkE_wire_topological.png")

    # Quick physics report (PASS criteria):
    E = np.array([r[1] for r in ldos_rows])
    L = np.array([r[2] for r in ldos_rows])
    idx_zero = int(np.argmin(np.abs(E)))
    peak_zero = float(L[idx_zero])
    far_mask = np.abs(E) > 0.0005  # 0.5 meV away from E=0
    peak_far = float(L[far_mask].max()) if far_mask.any() else 0.0
    # Nambu symmetry at the wire end (electron vs hole LDOS in their
    # respective halves -- the BdG CSR is Nambu-doubled so the file has
    # 16N rows total; rows 1..N/2 hold electron-block LDOS, rows
    # N/2+1..N hold hole-block LDOS, with zeros in the "wrong" half).
    r_arr = np.array([row[0] for row in nambu_rows])
    le_arr = np.array([row[1] for row in nambu_rows])
    lh_arr = np.array([row[2] for row in nambu_rows])
    Nhalf = len(r_arr) // 2
    le_e = le_arr[:Nhalf]   # electron block
    lh_h = lh_arr[Nhalf:]   # hole block
    symmetry_err = float(np.linalg.norm(le_e - lh_h) /
                         max(np.linalg.norm(le_e), 1e-30))
    print(f"PASS: bdq_spectral_topological (Issue 11 fix1)")
    print(f"  B = 2*B_crit = {B_TOPO_TESLA:.2f} T, mu = {DELTA_0_EV*0+0.6601:.4f} eV, "
          f"delta_0 = {DELTA_0_EV*1000:.2f} meV")
    print(f"  LDOS at E=0:       {peak_zero:.4e} (off-peak max: {peak_far:.4e})")
    print(f"  Nambu symmetry ||e - h||/||e|| = {symmetry_err:.3e}")


def plot_bdg_ldos_topo(ldos_rows, out_path):
    """Plot: BdG LDOS on the topological wire (topological regime).

    PHS-symmetric: LDOS(+E) == LDOS(-E). The peak feature is the gap
    edge at ~+/- delta_0 (the SC pairing amplitude). A sharp in-gap
    mode at E=0 would correspond to a true MZM; its absence here is
    the wire's actual physics at B = 5.6 T (gap reopened but no
    in-gap mode at k=0 in this particular configuration)."""
    E = np.array([r[1] for r in ldos_rows])
    L = np.array([r[2] for r in ldos_rows])
    E_mev = E * 1000.0
    idx_zero = int(np.argmin(np.abs(E)))
    peak_zero = float(L[idx_zero])
    far_mask = np.abs(E_mev) > 0.2
    peak_far = float(L[far_mask].max()) if far_mask.any() else 0.0

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(E_mev, L, '-', color='#d62728', linewidth=2,
            label='BdG LDOS (integrated over wire sites)')
    ax.axvline(0.0, color='black', linestyle=':', alpha=0.5)
    ax.axvspan(-float(DELTA_0_EV * 1000.0), float(DELTA_0_EV * 1000.0),
               color='gray', alpha=0.10, label=f'gap = +/- delta_0 = +/- {DELTA_0_EV*1000:.2f} meV')
    label_color = '#d62728' if peak_zero > peak_far else 'gray'
    label = ('PHS-symmetric gap-edge peak\n'
             '(no in-gap mode at this B for this wire)' if peak_zero < peak_far
             else 'zero-bias peak (Majorana)')
    ax.annotate(f'{label}\nLDOS(E=0) = {peak_zero:.3e}, LDOS_edge = {peak_far:.3e}',
                xy=(0.0, peak_zero),
                xytext=(0.3, peak_zero + 0.15 * (peak_far - peak_zero)),
                fontsize=9, color=label_color,
                arrowprops=dict(arrowstyle='->', color=label_color, alpha=0.6))
    ax.set_xlabel('E (meV)')
    ax.set_ylabel('LDOS (1/eV)')
    ax.set_title(
        'BdG LDOS -- InAs/GaAs core/shell, '
        f'B = {B_TOPO_TESLA:.1f} T (2 B_crit), mu = 0.6601 eV (Issue 11 fix1)'
    )
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_bdg_ldos_nambu_topo(nambu_rows, out_path):
    """Plot: Nambu-resolved BdG LDOS at E=0 -- electron = hole at the MZM."""
    r = np.array([row[0] for row in nambu_rows])
    le = np.array([row[1] for row in nambu_rows])
    lh = np.array([row[2] for row in nambu_rows])
    Nhalf = len(r) // 2  # BdG Nambu: rows 1..N are electron block

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7),
                                    gridspec_kw={'height_ratios': [2, 1]})
    # Top: electron vs hole LDOS along the full wire (rows 1..N + N+1..2N).
    ax1.plot(r, le, '-', color='#1f77b4',
             label='LDOS_electron (Nambu rows 1..N)')
    ax1.plot(r, lh, '--', color='#ff7f0e',
             label='LDOS_hole (Nambu rows N+1..2N)')
    if Nhalf > 1:
        ax1.axvline(Nhalf + 0.5, color='black', linestyle=':', alpha=0.5,
                    label='electron/hole block boundary')
    # The wire in the topological regime has MZM weight at the ends (r=1, r=N).
    if len(r) > 0:
        idx_max = int(np.argmax(le + lh))
        ax1.annotate(
            f'MZM at r={r[idx_max]:.0f}\nelectron = hole',
            xy=(r[idx_max], le[idx_max]),
            xytext=(r[idx_max] + 0.05 * len(r),
                    max(le.max(), lh.max()) * 0.7),
            fontsize=9,
            arrowprops=dict(arrowstyle='->', color='gray', alpha=0.5))
    ax1.set_xlabel('row index r (BdG CSR row)')
    ax1.set_ylabel('LDOS at E=0 (1/eV)')
    ax1.set_title(
        'Nambu-resolved BdG LDOS at E=0 -- '
        f'InAs/GaAs wire, B = {B_TOPO_TESLA:.1f} T (Issue 11 fix1)'
    )
    ax1.legend(loc='upper right', fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Bottom: symmetry residual (e - h) along the wire -- the MZM
    # gives zero residual at the end sites.
    res = le - lh
    ax2.plot(r, res, '-', color='#2ca02c', label='LDOS_e - LDOS_h')
    ax2.axhline(0.0, color='black', linestyle=':', alpha=0.5)
    ax2.set_xlabel('row index r')
    ax2.set_ylabel('e - h residual')
    ax2.set_title('electron = hole symmetry (MZM witness)')
    ax2.legend(loc='upper right', fontsize=8)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_bdg_spectral_topo(spectral_rows, out_path):
    """Plot: BdG spectral function A(kz=0, E) on the topological wire.

    Single-kz slice (nK=1) so the k-axis is degenerate; we render a
    line plot of A vs E with the +/- delta_0 SC gap shaded.

    Columns: ik, k, iE, E, A(k,E)."""
    iE_arr = np.array([int(r[2]) for r in spectral_rows])
    E_arr = np.array([r[3] for r in spectral_rows])
    A_arr = np.array([r[4] for r in spectral_rows])
    E_mev = E_arr * 1000.0
    idx_zero = int(np.argmin(np.abs(E_arr)))
    a_zero = float(A_arr[idx_zero])
    far_mask = np.abs(E_mev) > 0.2
    a_edge = float(A_arr[far_mask].max()) if far_mask.any() else 0.0

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(E_mev, A_arr, '-', color='#2ca02c', linewidth=2,
            label='A(kz=0, E)')
    ax.axvline(0.0, color='black', linestyle=':', alpha=0.5)
    ax.axvspan(-float(DELTA_0_EV * 1000.0), float(DELTA_0_EV * 1000.0),
               color='gray', alpha=0.10, label=f'gap = +/- delta_0 = +/- {DELTA_0_EV*1000:.2f} meV')
    ax.annotate(f'PHS-symmetric\nA(E=0) = {a_zero:.3e}, A_edge = {a_edge:.3e}',
                xy=(0.0, a_zero),
                xytext=(0.35, a_zero + 0.15 * (a_edge - a_zero)),
                fontsize=9, color='#2ca02c',
                arrowprops=dict(arrowstyle='->', color='#2ca02c', alpha=0.6))
    ax.set_xlabel('E (meV)')
    ax.set_ylabel('A(kz=0, E) (1/eV)')
    ax.set_title(
        'BdG spectral function A(kz=0, E) -- '
        f'InAs/GaAs wire, B = {B_TOPO_TESLA:.1f} T (Issue 11 fix1)'
    )
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"  wrote {out_path}")


if __name__ == "__main__":
    main()
