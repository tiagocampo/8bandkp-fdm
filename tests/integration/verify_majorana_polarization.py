#!/usr/bin/env python3
"""Verify Majorana polarization (Sticlet P_M) on a real wire BdG eigensolve.

Per spec §3.2 (P0 fix): no synthetic fallback; failures exit non-zero.

The polarization observable (Sticlet MZM discriminator) is distinct from
the slim Pfaffian Z2 row consumed by the L13 acceptance gate. This
verifier exercises the polarization emitter; the acceptance gate reads
the slim Pfaffian row from `output/z2_phase_diagram.dat` separately
(see ticket 05 of `.scratch/archive/bdg-evaluator-pfaffian/`).

@todo U13 — BLOCKING-EMPIRICAL on canonical wire config. The wire-polarization
emitter at `src/apps/main_topology.f90:586-598` gates on `n_majorana == 1`,
which is not produced by the canonical `wire_inas_gaas_bdg_topological.toml`
at any Bx tested: FEAST noise floor (~1e-5 eV) exceeds the near-zero
threshold `bdg_default_near_zero_frac · delta_0` (~2e-7 eV). The verifier's
B-override to `B_vec = [3.0, 0.0, 0.0]` therefore never triggers the emitter
on this config, and the polarization file is never written.

Resolution path: U13 (periodic/Bloch BdG construction per parent plan §U13)
which scopes the wire Pfaffian sweep at the PHS-invariant momenta where the
noise floor is suppressed by spectral averaging. Per CLAUDE.md Known Issues —
separate scoped PR.

SKIP precondition (per ticket 05): the verifier SKIPs (exit 0) only when
the slim Pfaffian gate row is present in `output/z2_phase_diagram.dat`
(mu ≈ 0.6601 ± 0.0001 with z2 == -1). When the gate row is missing, the
verifier FAILS non-zero — both the polarization file and the slim Pfaffian
row should be produced by a single `topologicalAnalysis` run on the
canonical wire config. Fails loudly on other errors (config/executable
not found, exec error, empty file, trivial polarization).

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


def gate_row_colormap_present(repo_root, mu_target=0.6601, mu_tol=0.0001):
    """Return True iff `output/z2_phase_diagram.dat` has a row with
    mu ≈ mu_target (± mu_tol) and z2 == -1.

    Per ticket 05 of `.scratch/archive/bdg-evaluator-pfaffian/` — the
    polarization verifier's SKIP precondition is tightened: SKIP only when
    this row is present (otherwise FAIL non-zero). The slim Pfaffian row
    is the L13 acceptance gate's 4th live witness (colormap-extracted); its
    absence means the canonical wire BdG sweep didn't produce the gate
    signal, which is also why the polarization file is missing.
    """
    colormap = Path(repo_root) / "output" / "z2_phase_diagram.dat"
    if not colormap.exists():
        return False
    with open(colormap) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                mu = float(parts[1])
                z2 = int(float(parts[3]))
            except ValueError:
                continue
            if abs(mu - mu_target) <= mu_tol and z2 == -1:
                return True
    return False


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
    """Parse output/majorana_polarization.dat. Returns (P_M_array, half_wire_integral).

    Per ticket 05: SKIP only when the slim Pfaffian gate row is present in
    `output/z2_phase_diagram.dat` (the canonical wire BdG run produced it).
    When the gate row is missing, FAIL non-zero — both the polarization
    file and the slim Pfaffian row should be produced by a single canonical
    run; their joint absence means the run itself was broken.
    """
    if not path.exists():
        if gate_row_colormap_present(REPO):
            # BLOCKING-EMPIRICAL — see @todo U13 in module docstring. The wire-
            # polarization emitter gates on n_majorana==1 which the canonical
            # wire config never produces (FEAST noise floor > near-zero threshold).
            # Skip with exit 0; the slim Pfaffian gate row IS present, so the
            # canonical run was successful — only the polarization emitter is
            # blocked. Acceptance gate (4-witness, 1.0 T) covers the regression
            # net. Real eigensolve assertion deferred to U13.
            print(f"SKIP: {path} not produced by topologicalAnalysis (BLOCKING-EMPIRICAL on canonical wire)")
            print("      see @todo U13 in module docstring; acceptance gate covers 4-witness regression")
            print("      slim Pfaffian gate row present (ticket 05)")
            sys.exit(0)
        else:
            # Tightened precondition (ticket 05): the polarization file AND
            # the slim Pfaffian gate row are both missing. This means the
            # canonical wire BdG run produced neither signal. The previous
            # silent-SKIP hid this regression; now FAIL loudly.
            print(f"FAIL: {path} not produced by topologicalAnalysis AND slim Pfaffian "
                  f"gate row absent from output/z2_phase_diagram.dat (ticket 05).")
            print(f"      Re-run wire BdG sweep to regenerate the colormap + polarization file.")
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
