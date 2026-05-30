"""Shared utilities for standard-star benchmark scripts.

Provides common functions for running executables, parsing output files,
comparing values, and formatting benchmark table rows.

Usage by standard-star scripts:
    from star_helpers import (
        run_executable, parse_eigenvalues, parse_gfactor,
        parse_absorption, parse_topology_result, compare_value,
        format_benchmark_row, print_benchmark_header,
    )
"""

import os
import shutil
import subprocess
import tempfile

try:
    import numpy as np
except ImportError:
    raise ImportError(
        "numpy is required for standard-star benchmarks. "
        "Install with: pip install numpy"
    )


# ---------------------------------------------------------------------------
# Tolerance tiers (KD6)
# ---------------------------------------------------------------------------
TOL_EXACT = 1e-12
TOL_ANALYTICAL = 0.01   # 1% default for analytical comparisons
TOL_NUMERICAL = 0.05    # 5% default for numerical comparisons

# ---------------------------------------------------------------------------
# Physical constants (CODATA 2018)
# ---------------------------------------------------------------------------
HBAR_EV_S = 6.582119569e-16   # hbar in eV*s
HBAR_J_S  = 1.054571817e-34   # hbar in J*s
E_CHARGE  = 1.602176634e-19   # elementary charge in C
M0_KG     = 9.1093837015e-31  # free electron mass in kg

# numpy 2.0 compatibility: np.trapezoid was introduced in 2.0, replacing np.trapz
# np.trapz was removed in numpy 2.0
trapz_fn = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)
if trapz_fn is None:
    raise ImportError("numpy integration function not found (need numpy >= 1.6)")


# ---------------------------------------------------------------------------
# Executable runner
# ---------------------------------------------------------------------------

def run_executable(exe_path, config_path, work_dir, timeout=300):
    """Run a Fortran executable in work_dir with the given config.

    Copies config to <work_dir>/input.toml, runs the executable,
    returns the output directory path.

    Returns:
        (returncode, output_dir) — returncode is -1 on timeout.

    Raises:
        FileNotFoundError: if exe_path does not exist
    """
    dst_cfg = os.path.join(work_dir, "input.toml")
    shutil.copy2(config_path, dst_cfg)

    output_dir = os.path.join(work_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    try:
        result = subprocess.run(
            [exe_path],
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        return result.returncode, output_dir
    except subprocess.TimeoutExpired:
        return -1, output_dir


def run_exe(build_dir, name, config_path, work_dir, timeout=300):
    """Resolve and run a Fortran executable by name.

    Resolves <build_dir>/src/<name>, checks it exists, then calls
    run_executable.

    Args:
        build_dir: path to build/ directory
        name: executable name (e.g. 'bandStructure')
        config_path: path to config file
        work_dir: temporary working directory
        timeout: execution timeout in seconds

    Returns:
        (returncode, output_dir)

    Raises:
        FileNotFoundError: if executable not found
    """
    exe_path = os.path.abspath(os.path.join(build_dir, "src", name))
    if not os.path.isfile(exe_path):
        raise FileNotFoundError(f"Executable not found: {exe_path}")
    return run_executable(exe_path, config_path, work_dir, timeout)


def roth_gfactor(ep, eg, delta_so):
    """Roth g-factor formula (Winkler 2003, Eq. 6.42).

    g = 2 - 2*EP*DeltaSO / (3*Eg*(Eg + DeltaSO))

    Args:
        ep: Kane interband matrix element EP in eV
        eg: band gap in eV
        delta_so: spin-orbit splitting in eV

    Returns:
        g-factor (dimensionless)
    """
    return 2.0 - 2.0 * ep * delta_so / (3.0 * eg * (eg + delta_so))


# ---------------------------------------------------------------------------
# Output parsers
# ---------------------------------------------------------------------------

def parse_eigenvalues(filepath):
    """Parse eigenvalues file, returning list of (k, [eigenvalues]).

    Each line: |k| eigenval1 eigenval2 ...
    Comment lines start with '#'.

    Returns:
        list of (float, list[float])
    """
    results = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals = [float(x) for x in line.split()]
            if len(vals) >= 2:
                k = vals[0]
                evals = vals[1:]
                results.append((k, evals))
    return results


def parse_gfactor(filepath):
    """Parse gfactor.dat, returning (gx, gy, gz).

    File format: single line with 3 Fortran-formatted floats.

    Returns:
        tuple of (float, float, float)
    """
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals = [float(x) for x in line.split()]
            if len(vals) >= 3:
                return (vals[0], vals[1], vals[2])
    raise RuntimeError(f"No g-factor data found in {filepath}")


def parse_absorption(filepath):
    """Parse absorption spectrum file.

    File format: comment-prefixed header + data lines with E(eV) alpha(cm^-1).

    Returns:
        list of (energy, absorption) tuples
    """
    data = np.loadtxt(filepath, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return [(row[0], row[1]) for row in data]


def parse_topology_result(filepath):
    """Parse topology_result.dat, returning dict of key-value pairs.

    File format: '# key: value' lines.

    Returns:
        dict mapping key strings to value strings
    """
    results = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line.startswith('#'):
                continue
            # Remove leading '# '
            content = line[1:].strip()
            if ':' in content:
                key, value = content.split(':', 1)
                results[key.strip()] = value.strip()
    return results


# ---------------------------------------------------------------------------
# Effective mass extraction (adaptive parabolic fit)
# ---------------------------------------------------------------------------

HBAR2_OVER_2M0 = 3.80998  # eV * Angstrom^2 (hbar^2 / (2*m0))


def extract_effective_mass(eig_path, cb_index=-1, r2_threshold=0.9999):
    """Extract CB effective mass via adaptive parabolic fit near k=0.

    Progressively narrows the k-range until R^2 > threshold, giving
    the parabolic-limit mass even for strongly non-parabolic bands.
    Falls back to the best fit if threshold is never met.

    Args:
        eig_path: path to eigenvalues.dat
        cb_index: eigenvalue index for CB (default -1 = highest)
        r2_threshold: R^2 quality threshold

    Returns:
        (m_star, E0, r2, k_max) or None on failure
    """
    data = parse_eigenvalues(eig_path)
    if len(data) < 3:
        return None

    ks = np.array([d[0] for d in data])
    es = np.array([d[1][cb_index] for d in data])

    # Use k >= 0 only, sorted
    mask = ks >= 0
    ks = ks[mask]
    es = es[mask]
    idx = np.argsort(ks)
    ks = ks[idx]
    es = es[idx]

    if len(ks) < 3:
        return None

    best = None
    for npts in range(len(ks), 2, -1):
        kf = ks[:npts]
        ef = es[:npts]
        A = np.column_stack([kf**2, np.ones_like(kf)])
        result = np.linalg.lstsq(A, ef, rcond=None)
        c2, c0 = result[0]

        e_pred = A @ np.array([c2, c0])
        ss_res = np.sum((ef - e_pred)**2)
        ss_tot = np.sum((ef - np.mean(ef))**2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

        if c2 > 0:
            m_star = HBAR2_OVER_2M0 / c2
        else:
            continue

        if best is None or r2 > best[2]:
            best = (m_star, c0, r2, kf[-1])

        if r2 >= r2_threshold:
            return m_star, c0, r2, kf[-1]

    return best


# ---------------------------------------------------------------------------
# Bir-Pikus strain formulas (biaxial [001])
# ---------------------------------------------------------------------------

def bir_pikus_biaxial_001(a0, a_sub, C11, C12, ac, av, b_dp,
                          delta_so=0.341, Eg=1.519):
    """Bir-Pikus band-edge shifts for biaxial [001] strain with LH-SO coupling.

    Sign convention matches compute_bp_scalar in strain_solver.f90:
      P_eps = -av * Tr(eps)
      Q_eps = -(b_dp/2) * (eps_zz - 0.5*(eps_xx + eps_yy))
      delta_Ec  = ac * Tr(eps)
      delta_EHH = -P_eps + Q_eps
      delta_ELH = -P_eps - Q_eps
      delta_ESO = -P_eps

    The LH-SO off-diagonal coupling QT2/sqrt(2) produces a 2x2 eigenvalue
    problem. Under compressive strain (InAs on GaAs): LHSO_low < LHSO_high < HH < CB.
    Under tensile strain (GaAs on larger substrate): LHSO_low < HH < LHSO_high < CB.

    Args:
        a0: film lattice constant (Angstrom)
        a_sub: substrate lattice constant (Angstrom)
        C11, C12: elastic constants (kbar)
        ac: CB hydrostatic deformation potential (eV)
        av: VB hydrostatic deformation potential (eV)
        b_dp: shear deformation potential (eV)
        delta_so: spin-orbit splitting (eV), default GaAs
        Eg: band gap (eV), default GaAs

    Returns:
        dict with strain tensor, Bir-Pikus shifts, LH-SO coupled eigenvalues,
        and strained band edges (all in eV unless noted).
    """
    eps_xx = (a_sub - a0) / a0
    eps_yy = eps_xx
    eps_zz = -2.0 * C12 / C11 * eps_xx

    Tr_eps = eps_xx + eps_yy + eps_zz
    P_eps = -av * Tr_eps
    Q_eps = -b_dp / 2.0 * (eps_zz - 0.5 * (eps_xx + eps_yy))

    delta_Ec  = ac * Tr_eps
    delta_EHH = -P_eps + Q_eps
    delta_ELH = -P_eps - Q_eps
    delta_ESO = -P_eps

    # LH-SO off-diagonal coupling
    QT2 = 2.0 * Q_eps
    coupling = QT2 / np.sqrt(2.0)

    E_SO_unstrained = -delta_so
    E_CB_unstrained = Eg

    E_HH_strained = delta_EHH
    E_CB_strained = E_CB_unstrained + delta_Ec

    # LH-SO coupled 2x2 system
    a = delta_ELH           # E_LH unstrained (= 0) + shift
    b = E_SO_unstrained + delta_ESO
    c = abs(coupling)
    mid = (a + b) / 2.0
    half_diff = (a - b) / 2.0
    root = np.sqrt(half_diff**2 + c**2)
    E_LHSO_low  = mid - root
    E_LHSO_high = mid + root

    HH_LH_splitting = E_HH_strained - E_LHSO_high

    return {
        "eps_xx": eps_xx, "eps_yy": eps_yy, "eps_zz": eps_zz,
        "Tr_eps": Tr_eps, "P_eps": P_eps, "Q_eps": Q_eps,
        "delta_Ec": delta_Ec, "delta_EHH": delta_EHH,
        "delta_ELH": delta_ELH, "delta_ESO": delta_ESO,
        "E_CB": E_CB_strained, "E_HH": E_HH_strained,
        "E_LHSO_low": E_LHSO_low, "E_LHSO_high": E_LHSO_high,
        "HH_LH_splitting": HH_LH_splitting,
        "QT2": QT2, "coupling": coupling,
    }


# ---------------------------------------------------------------------------
# Value comparison and benchmark formatting
# ---------------------------------------------------------------------------

def compare_value(actual, expected, tolerance, name, unit=""):
    """Compare actual vs expected value within tolerance.

    Uses relative tolerance when |expected| > 1e-14 eV (typical energy scale),
    absolute tolerance otherwise (for near-zero values).

    Args:
        actual: computed value
        expected: reference value
        tolerance: relative tolerance (fraction, e.g. 0.01 for 1%)
        name: observable name for reporting
        unit: physical unit string (e.g. "eV", "m0")

    Returns:
        (passed, delta, row_dict) where delta is relative error,
        and row_dict has keys for benchmark table formatting.
    """
    if abs(expected) > 1e-14:
        delta = abs(actual - expected) / abs(expected)
    else:
        delta = abs(actual - expected)

    passed = delta <= tolerance

    row = {
        'computed': actual,
        'expected': expected,
        'tolerance': tolerance,
        'delta': delta,
        'name': name,
        'unit': unit,
        'status': 'PASS' if passed else 'FAIL',
    }

    return passed, delta, row


def format_benchmark_row(material, observable, computed, expected,
                         reference, tolerance, delta, status):
    """Format a single benchmark table row as markdown.

    Args:
        material: material system name
        observable: observable name
        computed: computed value
        expected: expected/reference value
        reference: literature citation
        tolerance: tolerance tier description
        delta: relative error
        status: 'PASS' or 'FAIL'

    Returns:
        markdown table row string
    """
    def fmt(v):
        if abs(v) < 1e-3 or abs(v) > 1e4:
            return f"{v:.4e}"
        return f"{v:.6f}"

    return (
        f"| {material} | {observable} | {fmt(computed)} | {fmt(expected)} "
        f"| {reference} | {tolerance} | {delta:.2e} | {status} |"
    )


def print_benchmark_header():
    """Print standard benchmark table header."""
    print("| Material | Observable | Computed | Expected | Reference "
          "| Tolerance | Delta | Status |")
    print("|----------|-----------|----------|----------|----------"
          "|-----------|-------|--------|")
