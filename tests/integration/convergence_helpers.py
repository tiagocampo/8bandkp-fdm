"""Shared utilities for Richardson extrapolation convergence tests.

Provides core functions for Richardson extrapolation, Grid Convergence Index
(GCI) computation, convergence rate extraction, and JSON diagnostics output.

Used by convergence test scripts (test_qw_grid_convergence.py, etc.).
Complements star_helpers.py (which provides Fortran I/O and parsers).
"""

import json
import os

import numpy as np


# ---------------------------------------------------------------------------
# Richardson extrapolation
# ---------------------------------------------------------------------------

def richardson_extrapolate(h_vals, obs_vals, order=2):
    """Richardson extrapolation to h=0 using two-finest-grid formula.

    For a method of order p:
        E_exact = E(h1) + (E(h1) - E(h2)) / ((h2/h1)^p - 1)

    Args:
        h_vals: grid spacing values (any order)
        obs_vals: observable values at each grid spacing
        order: expected convergence order p (e.g. 2 for FD order 2)

    Returns:
        Extrapolated value at h=0

    Raises:
        ValueError: if fewer than 2 data points provided
    """
    h = np.asarray(h_vals, dtype=float)
    E = np.asarray(obs_vals, dtype=float)

    if len(h) < 2:
        raise ValueError("Need at least 2 data points for Richardson extrapolation")

    # Sort by h ascending (finest grid first)
    idx = np.argsort(h)
    h = h[idx]
    E = E[idx]

    # Use the two finest grids
    h1, h2 = h[0], h[1]
    E1, E2 = E[0], E[1]

    p = order
    ratio = (h2 / h1) ** p - 1.0
    if abs(ratio) < 1e-12:
        return E1

    return E1 + (E1 - E2) / ratio


# ---------------------------------------------------------------------------
# Three-grid Richardson extrapolation (higher accuracy)
# ---------------------------------------------------------------------------

def richardson_extrapolate_3grid(h_vals, obs_vals, order=2):
    """Three-grid Richardson extrapolation using Romberg-like refinement.

    Applies two rounds of Richardson: first eliminates O(h^p), then O(h^{2p}).
    Requires at least 3 data points sorted by h ascending.

    Args:
        h_vals: grid spacing values
        obs_vals: observable values
        order: convergence order p

    Returns:
        Extrapolated value using three finest grids
    """
    h = np.asarray(h_vals, dtype=float)
    E = np.asarray(obs_vals, dtype=float)

    idx = np.argsort(h)
    h = h[idx]
    E = E[idx]

    if len(h) < 3:
        return richardson_extrapolate(h, E, order)

    # First Richardson: eliminate O(h^p) from pairs
    r1 = richardson_extrapolate(h[:2], E[:2], order)
    r2 = richardson_extrapolate(h[1:3], E[1:3], order)

    # Second Richardson: eliminate O(h^{2p}) using h[0] and h[2]
    # Ratio of effective spacings
    hr = (h[1] / h[0]) / (h[2] / h[1])
    denom = hr ** order - 1.0
    if abs(denom) < 1e-12:
        return r1
    return r1 + (r1 - r2) / denom


# ---------------------------------------------------------------------------
# Grid Convergence Index (GCI)
# ---------------------------------------------------------------------------

def compute_gci(h_vals, obs_vals, order=2, safety_factor=1.25):
    """Compute Grid Convergence Index (Roache 1998).

    GCI = Fs * |E1 - E2| / (E1 * (r^p - 1))

    where r = h_coarse/h_fine, p = convergence order, Fs = safety factor.
    For two-grid Richardson: Fs=3. For three-grid: Fs=1.25.

    Returns GCI as a fraction (e.g. 0.01 = 1%).

    Args:
        h_vals: grid spacing values
        obs_vals: observable values
        order: convergence order p
        safety_factor: Roache safety factor (1.25 for 3-grid, 3 for 2-grid)

    Returns:
        GCI as a fraction, or None if computation fails
    """
    h = np.asarray(h_vals, dtype=float)
    E = np.asarray(obs_vals, dtype=float)

    if len(h) < 2:
        return None

    idx = np.argsort(h)
    h = h[idx]
    E = E[idx]

    h1, h2 = h[0], h[1]
    E1, E2 = E[0], E[1]

    r = h2 / h1
    if r <= 1.0:
        return None

    p = order
    denominator = E1 * (r ** p - 1.0)
    if abs(denominator) < 1e-30:
        return None

    return abs(safety_factor * (E1 - E2) / denominator)


# ---------------------------------------------------------------------------
# Convergence rate extraction
# ---------------------------------------------------------------------------

def extract_convergence_rates(h_vals, obs_vals, order=2):
    """Extract pairwise convergence rates between consecutive grid levels.

    For consecutive grids i, i+1 (h_i < h_{i+1}):
        p_obs = log(|E_i - E_{i+1}| / |E_{i+1} - E_{i+2}|) / log(h_{i+1}/h_i)

    Uses three consecutive levels to estimate rate.

    Args:
        h_vals: grid spacing values
        obs_vals: observable values
        order: expected order (unused, for API consistency)

    Returns:
        list of (h_mid, rate) tuples for each consecutive triple,
        or empty list if fewer than 3 points
    """
    h = np.asarray(h_vals, dtype=float)
    E = np.asarray(obs_vals, dtype=float)

    idx = np.argsort(h)
    h = h[idx]
    E = E[idx]

    if len(h) < 3:
        return []

    rates = []
    for i in range(len(h) - 2):
        e21 = abs(E[i+1] - E[i])
        e32 = abs(E[i+2] - E[i+1])
        r21 = h[i+1] / h[i]

        if e21 < 1e-30 or e32 < 1e-30:
            continue

        rate = np.log(e32 / e21) / np.log(r21)
        rates.append((h[i+1], rate))

    return rates


def max_convergence_rate(h_vals, obs_vals, order=2):
    """Maximum observed convergence rate across all consecutive grid pairs.

    Uses the max-rate strategy: round-off causes rates to drop, never increase,
    so the maximum is the best estimator of the true asymptotic rate.

    Args:
        h_vals: grid spacing values
        obs_vals: observable values
        order: expected order

    Returns:
        (max_rate, all_rates) tuple, or (0.0, []) if insufficient data
    """
    rates = extract_convergence_rates(h_vals, obs_vals, order)
    if not rates:
        return 0.0, []
    rate_values = [r[1] for r in rates]
    return max(rate_values), rate_values


# ---------------------------------------------------------------------------
# Monotonicity check
# ---------------------------------------------------------------------------

def check_monotonic(obs_vals, direction='both'):
    """Check if observable values converge monotonically.

    Args:
        obs_vals: values ordered from finest to coarsest grid
        direction: 'increasing', 'decreasing', or 'both' (either direction)

    Returns:
        (is_monotonic, direction_found) tuple
    """
    E = np.asarray(obs_vals, dtype=float)
    if len(E) < 2:
        return True, 'flat'

    diffs = np.diff(E)
    if np.all(diffs >= -1e-14):
        return True, 'increasing'
    if np.all(diffs <= 1e-14):
        return True, 'decreasing'
    return False, 'non-monotonic'


# ---------------------------------------------------------------------------
# Absorption edge extraction
# ---------------------------------------------------------------------------

def extract_absorption_edge(spectrum, fraction=0.1):
    """Extract absorption edge from spectrum data.

    The edge is defined as the energy where absorption crosses a fixed
    fraction of the peak value, using linear interpolation.

    Args:
        spectrum: list of (energy, absorption) tuples
        fraction: fraction of peak value defining the edge (default 0.1)

    Returns:
        edge energy in same units as spectrum, or None on failure
    """
    if not spectrum:
        return None

    energies = np.array([s[0] for s in spectrum])
    absorption = np.array([s[1] for s in spectrum])

    # Only consider positive absorption
    mask = absorption > 0
    if not np.any(mask):
        return None

    peak = np.max(absorption[mask])
    threshold = fraction * peak

    # Find first crossing above threshold
    for i in range(len(absorption) - 1):
        if absorption[i] < threshold <= absorption[i + 1]:
            # Linear interpolation
            e0, e1 = energies[i], energies[i + 1]
            a0, a1 = absorption[i], absorption[i + 1]
            if abs(a1 - a0) < 1e-30:
                return e0
            return e0 + (threshold - a0) * (e1 - e0) / (a1 - a0)

    return None


# ---------------------------------------------------------------------------
# Exciton output parsing
# ---------------------------------------------------------------------------

def parse_exciton_file(filepath):
    """Parse output/exciton.dat file.

    Returns dict with lambda_opt_AA, E_binding_meV, mu_over_m0, eps_r.
    """
    try:
        data = np.loadtxt(filepath, comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] >= 2:
            return {
                'lambda_opt_AA': float(data[0, 0]),
                'E_binding_meV': float(data[0, 1]),
                'mu_over_m0': float(data[0, 2]) if data.shape[1] > 2 else None,
                'eps_r': float(data[0, 3]) if data.shape[1] > 3 else None,
            }
    except (ValueError, OSError, IndexError):
        pass
    return None


# ---------------------------------------------------------------------------
# JSON diagnostics output
# ---------------------------------------------------------------------------

def write_convergence_json(results, output_path):
    """Write convergence results to JSON file.

    Args:
        results: dict with convergence data
        output_path: path to output JSON file
    """
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)


def make_convergence_report(system, observable, h_vals, obs_vals, order=2,
                            theoretical_rate=None, rate_tolerance=None,
                            gci_threshold=0.1):
    """Build a complete convergence report dict for one observable.

    Args:
        system: system name (e.g. 'S4_GaAs_AlGaAs_QW')
        observable: observable name (e.g. 'CB1_energy')
        h_vals: grid spacings
        obs_vals: observed values at each spacing
        order: FD order used
        theoretical_rate: expected theoretical rate (e.g. 2.0)
        rate_tolerance: allowed deviation from theoretical rate

    Returns:
        dict with all convergence diagnostics
    """
    h = np.asarray(h_vals, dtype=float)
    E = np.asarray(obs_vals, dtype=float)

    report = {
        'system': system,
        'observable': observable,
        'fd_order': order,
        'grid_data': [
            {'h': float(h[i]), 'value': float(E[i])}
            for i in range(len(h))
        ],
    }

    # Richardson extrapolation
    try:
        rich_val = richardson_extrapolate(h, E, order)
        report['richardson_extrapolated'] = float(rich_val)
    except (ValueError, IndexError) as e:
        report['richardson_extrapolated'] = None
        report['richardson_error'] = str(e)

    # Three-grid Richardson if enough data
    if len(h) >= 3:
        try:
            rich3 = richardson_extrapolate_3grid(h, E, order)
            report['richardson_3grid'] = float(rich3)
        except (ValueError, IndexError):
            pass

    # GCI
    gci = compute_gci(h, E, order)
    report['gci_fraction'] = float(gci) if gci is not None else None
    report['gci_percent'] = float(gci * 100) if gci is not None else None

    # Convergence rates
    max_rate, all_rates = max_convergence_rate(h, E, order)
    report['observed_rates'] = [float(r) for r in all_rates]
    report['max_observed_rate'] = float(max_rate)

    # Rate assertion
    if theoretical_rate is not None and max_rate > 0:
        if rate_tolerance is not None:
            rate_pass = bool(abs(max_rate - theoretical_rate) <= rate_tolerance)
            report['rate_assertion'] = {
                'theoretical': float(theoretical_rate),
                'observed_max': float(max_rate),
                'tolerance': float(rate_tolerance),
                'passed': rate_pass,
            }

    # Monotonicity
    mono, direction = check_monotonic(E)
    report['monotonic'] = bool(mono)
    report['monotonic_direction'] = direction

    # Overall pass/fail
    failures = []
    if not mono:
        failures.append('non-monotonic convergence')
    if theoretical_rate is not None and rate_tolerance is not None and max_rate > 0:
        if abs(max_rate - theoretical_rate) > rate_tolerance:
            failures.append(f'rate {max_rate:.2f} outside tolerance of {theoretical_rate} +/- {rate_tolerance}')
    if gci is not None and gci > gci_threshold:
        failures.append(f'GCI {gci*100:.1f}% exceeds {gci_threshold*100:.0f}%')
    report['passed'] = bool(len(failures) == 0)
    report['failures'] = failures

    return report
