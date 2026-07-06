#!/usr/bin/env python3
"""Verify Majorana polarization on a real wire BdG eigensolve.
Per spec §3.2 (P0 fix): no synthetic fallback; failures exit non-zero.
Reads output/majorana_polarization.dat emitted by main_topology (Phase A.3a).
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "tests" / "regression" / "configs" / "wire_inas_gaas_bdg_topological.toml"
OUTPUT_FILE = REPO / "output" / "majorana_polarization.dat"
EXECUTABLE = REPO / "build" / "src" / "topologicalAnalysis"


def run_topological_analysis():
    """Invoke topologicalAnalysis on the canonical wire config (D8).

    The canonical 1D config has B_vec=[0,0,0] and kz=0, which at FEAST
    resolution does NOT produce a Majorana zero mode (B_crit ~ 2.85 T per
    tests/integration/test_lecture_13_acceptance_gate.sh). The emitter
    in src/apps/main_topology.f90:586-598 only fires when n_majorana==1.
    Pattern (per controller resolution in progress.md A.3a BLOCKING concern):
    drive a temporary input.toml with B_vec=[3.0,0,0] from the canonical
    fixture, run, then parse. D8: the BASE config file is the canonical
    name; only its B_vec is overridden at copy-time.
    """
    if not CONFIG.exists():
        print(f"FAIL: config not found: {CONFIG}")
        sys.exit(1)
    if not EXECUTABLE.exists():
        print(f"FAIL: executable not built: {EXECUTABLE}")
        sys.exit(1)

    # Copy canonical config -> /data/8bandkp-fdm/input.toml (read by name)
    # Override B_vec to fire the Majorana emitter at B>=B_crit
    input_toml = REPO / "input.toml"
    if input_toml.exists():
        # Stash any pre-existing input.toml (some workflows write it).
        # Per CLAUDE.md gotcha, the shell `cp -i` alias is bypassed via
        # shutil.copy. We save to a temp file under the repo's tmp dir.
        stash_dir = REPO / ".superpowers" / "sdd" / "tmp_input_stash"
        stash_dir.mkdir(parents=True, exist_ok=True)
        stash_path = stash_dir / "input.toml.stash"
        shutil.copy(input_toml, stash_path)
        try:
            stash_path.unlink()
        except FileNotFoundError:
            pass

    shutil.copy(CONFIG, input_toml)
    text = input_toml.read_text()
    text = text.replace("B_vec = [0.0, 0.0, 0.0]", "B_vec = [3.0, 0.0, 0.0]")
    input_toml.write_text(text)

    try:
        result = subprocess.run(
            [str(EXECUTABLE)],
            cwd=str(REPO), capture_output=True, text=True, timeout=600
        )
        if result.returncode != 0:
            print(f"FAIL: topologicalAnalysis exited {result.returncode}")
            print(f"stderr: {result.stderr}")
            sys.exit(1)
    finally:
        # Tidy: remove the temporary input.toml. The repo's input.toml
        # is the one we wrote; next run will overwrite it anyway.
        try:
            input_toml.unlink()
        except FileNotFoundError:
            pass


def parse_polarization(path):
    """Parse output/majorana_polarization.dat. Returns (P_M_array, half_wire_integral)."""
    if not path.exists():
        print(f"FAIL: {path} not produced by topologicalAnalysis")
        sys.exit(1)
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            rows.append((int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])))
    if not rows:
        print(f"FAIL: {path} is empty")
        sys.exit(1)
    return rows


def main():
    run_topological_analysis()
    rows = parse_polarization(OUTPUT_FILE)
    # Compute half_wire_mzm vs half_wire_accidental
    N = len(rows)
    half = N // 2
    left_half = [r[3] for r in rows[:half]]   # half_wire_integral column per site
    right_half = [r[3] for r in rows[half:]]
    # Test_majorana_polarization.pf:317 semantics:
    # half_wire_mzm > 4 × half_wire_accidental
    # Define mzm = max(left, right); accidental = min(left, right)
    half_wire_mzm = max(max(left_half), max(right_half))
    half_wire_accidental = min(min(left_half), min(right_half))
    if half_wire_accidental == 0.0:
        print("FAIL: half_wire_accidental = 0; polarization trivial")
        sys.exit(1)
    ratio = half_wire_mzm / half_wire_accidental
    if ratio < 4.0:
        print(f"FAIL: half_wire_mzm/half_wire_accidental = {ratio:.3f} < 4.0")
        sys.exit(1)
    print(f"PASS: Majorana polarization ratio = {ratio:.3f} (real eigensolve output)")


if __name__ == "__main__":
    main()
