#!/usr/bin/env python3
# COVERAGE: observable=bdg_spectrum geometry=wire material=InAs/GaAs tier=numerical
# Regression test for Issue #04 — strain-aware wire setup type.
#
# Runs topologicalAnalysis (mode="bdg") on a strained InAs/GaAs core/shell
# wire TWICE: once with [strain] enabled (as the config ships) and once with
# the [strain] block removed (no strain field). The BdG eigenvalue spectra
# MUST differ — strain shifts the band edges (Bir-Pikus deformation
# potentials on a lattice-mismatched InAs/GaAs wire), which shifts the BdG
# spectrum.
#
# Before #04 the wire BdG driver (run_bdg_wire in main_topology.f90) called
# confinementInitialization_2d but NEVER compute_strain +
# compute_bir_pikus_blocks, so cfg%strain_blocks was left unallocated and
# ZB8bandGeneralized silently skipped the strain insertion. Both runs then
# produced identical spectra = the latent physics defect this test catches.
#
# Usage: verify_wire_bdg_strain_shift.py <topologicalAnalysis_exe> <config.toml>
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from star_helpers import run_executable  # noqa: E402


def _strip_strain_block(text):
    """Remove the [strain] section and its key lines from a TOML config."""
    lines = text.splitlines()
    out, in_block = [], False
    for line in lines:
        s = line.strip()
        if s.startswith("[") and s.endswith("]"):
            in_block = s == "[strain]"
        if in_block:
            continue
        out.append(line)
    while out and out[-1].strip() == "":
        out.pop()
    return "\n".join(out) + "\n"


def parse_bdg_eigenvalues(filepath):
    """Parse a bdg_eigenvalues.dat file (columns: index, energy)."""
    vals = []
    if not Path(filepath).is_file():
        return vals
    for line in Path(filepath).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) >= 2:
            try:
                vals.append(float(parts[1]))
            except ValueError:
                pass
    return vals


def run_one(exe, config_text, workdir, tag):
    """Write config, run topologicalAnalysis, return parsed BdG eigenvalues."""
    cfg_path = workdir / f"input_{tag}.toml"
    cfg_path.write_text(config_text)
    sub = workdir / tag
    sub.mkdir(parents=True, exist_ok=True)
    rc, out_dir = run_executable(exe, cfg_path, sub, timeout=600)
    if rc != 0:
        print(f"FAIL: {tag} run exited with rc={rc}")
        return None
    return parse_bdg_eigenvalues(Path(out_dir) / "bdg_eigenvalues.dat")


def main():
    EXE = sys.argv[1]
    CONFIG = Path(sys.argv[2])
    base_text = CONFIG.read_text()
    if "[strain]" not in base_text:
        print("FAIL: config has no [strain] block to toggle")
        return 1

    strained_text = base_text
    unstrained_text = _strip_strain_block(base_text)
    # The stripper must remove the [strain] table header line. (Comments that
    # merely mention the word "strain" are fine — the parser ignores them.)
    if any(ln.strip() == "[strain]" for ln in unstrained_text.splitlines()):
        print("FAIL: strain-block stripper did not remove the [strain] header")
        return 1

    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        eigs_strained = run_one(EXE, strained_text, workdir, "strained")
        if not eigs_strained:
            print("FAIL: strained run produced no BdG eigenvalues")
            return 1
        eigs_unstrained = run_one(EXE, unstrained_text, workdir, "unstrained")
        if not eigs_unstrained:
            print("FAIL: unstrained run produced no BdG eigenvalues")
            return 1

    n = min(len(eigs_strained), len(eigs_unstrained))
    if n == 0:
        print("FAIL: no overlapping eigenvalues to compare")
        return 1

    max_abs_diff = 0.0
    for i in range(n):
        d = abs(eigs_strained[i] - eigs_unstrained[i])
        if d > max_abs_diff:
            max_abs_diff = d

    print(f"  compared {n} BdG eigenvalues (strained vs unstrained)")
    print(f"  max |E_strained - E_unstrained| = {max_abs_diff:.6e} eV")
    print(f"  strained   first 3: {eigs_strained[:3]}")
    print(f"  unstrained first 3: {eigs_unstrained[:3]}")

    # The strain shift on InAs/GaAs is tens of meV; require a clearly
    # non-zero shift well above FEAST's 1e-10 tolerance.
    if max_abs_diff < 1.0e-3:
        print(f"FAIL: strain did NOT shift the BdG spectrum "
              f"(max diff {max_abs_diff:.3e} eV < 1e-3 eV)")
        print("  This is the latent #04 strain-omission bug: the wire BdG "
              "driver skipped compute_strain + compute_bir_pikus_blocks, "
              "so cfg%strain_blocks was unallocated and the strain "
              "insertion in ZB8bandGeneralized was silently skipped.")
        return 1

    print("PASS: strained-wire BdG spectrum shifts with strain applied "
          "(strain-omission bug #04 fixed)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
